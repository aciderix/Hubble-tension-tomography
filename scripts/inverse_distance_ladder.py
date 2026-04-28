"""Inverse distance ladder: take H0=73 from SH0ES, fit Pantheon+ + DESI DR2,
extract the required sound horizon r_d. Compare to BBN+matter-density prediction.

If the inferred r_d differs significantly from r_d^Planck ~ 147 Mpc, that
quantifies the *minimal* shift that any pre-recombination physics must produce.
"""
import sys
from pathlib import Path
import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.data_loaders import load_pantheon_plus
from src.cosmo import distance_modulus, rd_aubourg
from src.desi_bao import load_desi_dr2_all, predict_obs


def chi2_pantheon_plus(H0, Om, M, dat, use_calib=True):
    mu_th = distance_modulus(dat["zHEL"], dat["zHD"], H0, Om)
    mu_model = np.where(dat["is_calibrator"], dat["ceph_dist"], mu_th)
    resid = dat["mb"] - (M + mu_model)
    if not use_calib:
        m = ~dat["is_calibrator"]
        resid = resid[m]
        cov = dat["cov"][np.ix_(m, m)]
    else:
        cov = dat["cov"]
    L = np.linalg.cholesky(cov)
    y = np.linalg.solve(L, resid)
    return float(y @ y)


def chi2_bao(H0, Om, rd, bao):
    pred = predict_obs(bao["z"], bao["obs"], H0, Om, rd)
    d = bao["val"] - pred
    return float(d @ bao["invcov"] @ d)


def main():
    pp = load_pantheon_plus()
    bao = load_desi_dr2_all()
    print(f"Pantheon+: {pp['n']} SNe, {pp['is_calibrator'].sum()} calibrators")
    print(f"DESI DR2: {bao['n']} BAO measurements\n")

    # ---------- Fit 1: Pantheon+SH0ES alone (baseline) ----------
    print("=== 1. Pantheon+SH0ES only (LCDM) ===")
    f1 = lambda x: chi2_pantheon_plus(x[0], x[1], x[2], pp)
    r1 = minimize(f1, [73, 0.33, -19.25], method="Nelder-Mead",
                  options=dict(xatol=1e-4, fatol=1e-3, maxiter=2000))
    H0_sh, Om_sh, M_sh = r1.x
    print(f"  H0={H0_sh:.3f}  Om={Om_sh:.4f}  M={M_sh:.4f}  chi2={r1.fun:.2f}")

    # ---------- Fit 2: DESI DR2 alone, free (H0, Om, rd) ----------
    # NOTE: DESI BAO alone is degenerate in H0*rd.
    # Fix H0=70 (mid) and fit Om, rd to show the degeneracy
    print("\n=== 2. DESI DR2 alone, free (H0, Om, rd) — degenerate ===")
    for H0_fix in [67.4, 70.0, 73.0]:
        f2 = lambda x: chi2_bao(H0_fix, x[0], x[1], bao)
        r2 = minimize(f2, [0.3, 147], method="Nelder-Mead",
                      options=dict(xatol=1e-4, fatol=1e-3, maxiter=1000))
        Om2, rd2 = r2.x
        print(f"  H0={H0_fix}: Om*={Om2:.4f}, rd*={rd2:.2f} Mpc, chi2={r2.fun:.2f}, "
              f"H0*rd={H0_fix*rd2:.1f}")

    # ---------- Fit 3: Pantheon+SH0ES + DESI DR2 jointly, free rd ----------
    print("\n=== 3. JOINT Pantheon+SH0ES + DESI DR2, free (H0, Om, M, rd) ===")
    def f3(x):
        H0, Om, M, rd = x
        if not (50 < H0 < 90 and 0.1 < Om < 0.6 and -22 < M < -18 and 100 < rd < 200):
            return 1e30
        return chi2_pantheon_plus(H0, Om, M, pp) + chi2_bao(H0, Om, rd, bao)

    r3 = minimize(f3, [73, 0.33, -19.25, 140.0], method="Nelder-Mead",
                  options=dict(xatol=1e-4, fatol=1e-3, maxiter=4000))
    H0j, Omj, Mj, rdj = r3.x
    chi2_pp = chi2_pantheon_plus(H0j, Omj, Mj, pp)
    chi2_bao_j = chi2_bao(H0j, Omj, rdj, bao)
    print(f"  H0={H0j:.3f}  Om={Omj:.4f}  M={Mj:.4f}  rd={rdj:.2f} Mpc")
    print(f"  chi2: Pantheon+SH0ES = {chi2_pp:.2f},  DESI = {chi2_bao_j:.2f},  "
          f"total = {r3.fun:.2f}")

    # ---------- Compare inferred rd to BBN+matter predicted rd ----------
    h = H0j / 100.0
    omh2 = Omj * h ** 2
    # Cooke+ 2018 BBN: Omega_b h^2 = 0.02218 +/- 0.00055
    obh2_bbn = 0.02218
    rd_pred = rd_aubourg(omh2, obh2_bbn)
    rd_planck = 147.05
    print(f"\n  inferred rd from joint fit: {rdj:.2f} Mpc")
    print(f"  predicted rd (omh2={omh2:.4f}, BBN obh2={obh2_bbn}): {rd_pred:.2f} Mpc")
    print(f"  Planck 2018 rd: {rd_planck:.2f} Mpc")
    print(f"  shift required: rd_joint / rd_predicted = {rdj/rd_pred:.4f} "
          f"({100*(rdj/rd_pred-1):+.2f}%)")
    print(f"  shift required: rd_joint / rd_Planck    = {rdj/rd_planck:.4f} "
          f"({100*(rdj/rd_planck-1):+.2f}%)")

    # ---------- Fit 4: Pantheon+ (no SH0ES) + DESI DR2 + Planck rd ----------
    print("\n=== 4. Pantheon+ HF only + DESI DR2 + rd=147 (no SH0ES anchor) ===")
    def f4(x):
        H0, Om, M = x
        if not (50 < H0 < 90 and 0.1 < Om < 0.6 and -22 < M < -18):
            return 1e30
        return (chi2_pantheon_plus(H0, Om, M, pp, use_calib=False)
                + chi2_bao(H0, Om, 147.05, bao))

    r4 = minimize(f4, [70, 0.3, -19.4], method="Nelder-Mead",
                  options=dict(xatol=1e-4, fatol=1e-3, maxiter=4000))
    print(f"  H0={r4.x[0]:.3f}  Om={r4.x[1]:.4f}  M={r4.x[2]:.4f}  "
          f"chi2={r4.fun:.2f}")
    print("  -> this is the 'CMB-side' inverse-ladder H0 with rd=Planck")

    # ---------- Fit 5: Pantheon+SH0ES + DESI DR2 + rd=Planck (forced consistency) ----------
    print("\n=== 5. Pantheon+SH0ES + DESI DR2 + rd FIXED to Planck 147.05 ===")
    def f5(x):
        H0, Om, M = x
        if not (50 < H0 < 90 and 0.1 < Om < 0.6 and -22 < M < -18):
            return 1e30
        return (chi2_pantheon_plus(H0, Om, M, pp)
                + chi2_bao(H0, Om, 147.05, bao))

    r5 = minimize(f5, [70, 0.32, -19.3], method="Nelder-Mead",
                  options=dict(xatol=1e-4, fatol=1e-3, maxiter=4000))
    chi2_pp5 = chi2_pantheon_plus(r5.x[0], r5.x[1], r5.x[2], pp)
    chi2_bao_5 = chi2_bao(r5.x[0], r5.x[1], 147.05, bao)
    print(f"  H0={r5.x[0]:.3f}  Om={r5.x[1]:.4f}  M={r5.x[2]:.4f}")
    print(f"  chi2: Pantheon+SH0ES = {chi2_pp5:.2f},  DESI = {chi2_bao_5:.2f},  "
          f"total = {r5.fun:.2f}")
    print(f"  Δχ² vs free-rd fit: {r5.fun - r3.fun:+.2f}  "
          f"(if >9 => >3σ tension)")


if __name__ == "__main__":
    main()
