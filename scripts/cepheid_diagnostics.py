"""Cepheid-level diagnostics from SH0ES R22.

Read the master Cepheid catalogue (R22_orig19_NIR.out: ~3000 Cepheids in
~37 anchor galaxies) and look at:

  - Period-Luminosity-Metallicity (PLZ) residuals per host
  - Slope of residuals vs metallicity ([O/H] - 8.69)
  - Slope of residuals vs period
  - Are any hosts outliers?
  - Effect of dropping each host: H0 sensitivity
"""
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

ROOT = Path(__file__).resolve().parents[1]
PATH = ROOT / "external/SH0ES_Data/R22_orig19_NIR.out"


def load():
    rows = []
    with open(PATH) as fh:
        for ln in fh:
            ln = ln.rstrip()
            if not ln or ln.startswith("Host") or ln.startswith("---"):
                continue
            parts = ln.split()
            if len(parts) < 11:
                continue
            host = parts[0]
            try:
                ra = float(parts[1]); dec = float(parts[2])
                period = float(parts[4]); VI = float(parts[5]); sigVI = float(parts[6])
                H = float(parts[7]); sigH = float(parts[8]); metal = float(parts[9])
            except ValueError:
                continue
            rows.append((host, ra, dec, period, VI, sigVI, H, sigH, metal))
    dt = [("host", "U16"), ("ra", "f8"), ("dec", "f8"), ("period", "f8"),
          ("VI", "f8"), ("sigVI", "f8"), ("H", "f8"), ("sigH", "f8"),
          ("metal", "f8")]
    return np.array(rows, dtype=dt)


def main():
    cep = load()
    print(f"Loaded {len(cep)} Cepheids in {len(np.unique(cep['host']))} hosts")
    hosts = np.unique(cep["host"])
    print(f"Hosts: {', '.join(hosts)}\n")

    # The Wesenheit magnitude is the key observable used by SH0ES:
    # m_W = m_H - R*(V-I), where R is the reddening coefficient
    # R = 0.386 (Riess+2022, Wesenheit NIR)
    R = 0.386
    mW = cep["H"] - R * cep["VI"]
    sigmW = np.sqrt(cep["sigH"] ** 2 + (R * cep["sigVI"]) ** 2)

    # Per-host best-fit Wesenheit magnitude (zero-point + period slope)
    # Standard PLZ: m_W = mu_host + b * (log10 P - 1) + Z * [O/H] + M_W
    # We use: log10(P_d / 10) for centred period
    print("=== Per-host PL fit:  m_W = mu_host + b*(logP - 1) + ε ===")
    print(f"  {'host':<10s} {'n':>4s} {'mu_host':>10s} {'b':>8s} {'rms':>8s}")
    mu_per_host = {}
    for h in hosts:
        m = cep["host"] == h
        if m.sum() < 6:
            continue
        x = np.log10(cep["period"][m]) - 1.0
        y = mW[m]
        w = 1.0 / sigmW[m] ** 2
        # Weighted linear fit y = a + b x
        Sw = w.sum()
        Swx = (w * x).sum()
        Swy = (w * y).sum()
        Swxx = (w * x * x).sum()
        Swxy = (w * x * y).sum()
        det = Sw * Swxx - Swx * Swx
        a = (Swxx * Swy - Swx * Swxy) / det
        b = (Sw * Swxy - Swx * Swy) / det
        rms = np.sqrt(((y - a - b * x) ** 2).mean())
        mu_per_host[h] = a
        print(f"  {h:<10s} {m.sum():4d}  {a:10.4f}  {b:+.4f}  {rms:.4f}")

    # ------------------------------------------------------------
    # Metallicity dependence: do Cepheids in metal-richer hosts have
    # systematically different m_W vs PL relation?
    # ------------------------------------------------------------
    print("\n=== Global fit: m_W = mu_host + b (logP - 1) + Z*[O/H-8.69] ===")
    # Fit globally with per-host intercept + universal slope b + metallicity Z
    valid = sigmW > 0
    cep_v = cep[valid]; mW_v = mW[valid]; sigmW_v = sigmW[valid]
    hosts_v = cep_v["host"]
    unique_hosts_v = np.unique(hosts_v)
    n_h = len(unique_hosts_v)
    n_obs = len(cep_v)
    # Build design matrix: columns = [host_indicator (n_h cols), logP, [O/H]]
    A = np.zeros((n_obs, n_h + 2))
    host_idx = {h: i for i, h in enumerate(unique_hosts_v)}
    for j, h in enumerate(hosts_v):
        A[j, host_idx[h]] = 1.0
    A[:, n_h] = np.log10(cep_v["period"]) - 1.0
    A[:, n_h + 1] = cep_v["metal"]  # already [O/H]-8.69
    W = 1.0 / sigmW_v ** 2
    AtWA = A.T @ (A * W[:, None])
    AtWy = A.T @ (W * mW_v)
    params = np.linalg.solve(AtWA, AtWy)
    cov = np.linalg.inv(AtWA)
    err = np.sqrt(np.diag(cov))
    pred = A @ params
    chi2 = ((mW_v - pred) ** 2 * W).sum()
    print(f"  global b (logP slope) = {params[n_h]:.4f} ± {err[n_h]:.4f}")
    print(f"  global Z (metal coef) = {params[n_h+1]:.4f} ± {err[n_h+1]:.4f}")
    print(f"  fit chi2 = {chi2:.1f}, n_obs = {n_obs}, n_par = {n_h+2}, "
          f"reduced = {chi2/(n_obs-n_h-2):.3f}")
    # SH0ES: b_NIR ~ -3.30, Z ~ -0.20 (significantly negative)
    print("  SH0ES R22 reports b ~ -3.30, Z ~ -0.20 (so metal-rich Cep are brighter).")

    # ------------------------------------------------------------
    # Sensitivity test: drop one host at a time, re-fit, see Δ(intercept)
    # ------------------------------------------------------------
    print("\n=== Drop-one-host sensitivity ===")
    print("(robust check: does any single host dominate the calibration?)")
    print(f"  {'dropped':<10s} {'b':>8s} {'Z':>8s}")
    for h_drop in unique_hosts_v:
        keep = hosts_v != h_drop
        if keep.sum() < 50:
            continue
        Ak = A[keep]; mk = mW_v[keep]; Wk = W[keep]
        # Drop the column corresponding to that host (it has no rows now)
        # Easier: just re-build
        cv = cep_v[keep]
        unique_h = np.unique(cv["host"])
        n_h2 = len(unique_h)
        A2 = np.zeros((keep.sum(), n_h2 + 2))
        idx2 = {h: i for i, h in enumerate(unique_h)}
        for j, h in enumerate(cv["host"]):
            A2[j, idx2[h]] = 1.0
        A2[:, n_h2] = np.log10(cv["period"]) - 1.0
        A2[:, n_h2 + 1] = cv["metal"]
        try:
            p = np.linalg.solve(A2.T @ (A2 * Wk[:, None]), A2.T @ (Wk * mk))
            print(f"  {h_drop:<10s}  {p[n_h2]:+.4f}  {p[n_h2+1]:+.4f}")
        except Exception:
            pass


if __name__ == "__main__":
    main()
