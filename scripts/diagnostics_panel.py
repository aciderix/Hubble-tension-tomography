"""Multi-angle diagnostics on Pantheon+SH0ES.

A machine going through the data and looking for structure that LCDM does
not predict:

  1. H0(z) running: bin SNe in z, fit (H0, M) per bin, see if H0 trends with z.
  2. Anisotropic H0: bin SNe by sky direction (NGP/SGP, hemispheres along
     RA/DEC axes, dipole projection along CMB-dipole direction).
  3. Residual correlations: Pearson r between Hubble residual and every
     other column (host mass, x1, c, MWEBV, peculiar velocity, ...).
  4. Survey-specific H0: split by IDSURVEY.
  5. Per-anchor drop-one: drop each calibrator galaxy, see how H0 moves.

Aim: find any column that drives the tension or any subsample that
disagrees with the rest.
"""
import sys
from pathlib import Path
import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.data_loaders import load_pantheon_plus
from src.cosmo import distance_modulus_grid


def fit_HM(zHEL, zCMB, mb, dmb, M_init=-19.3, H0_init=70.0):
    """Fit (H0, M) for a fixed Om=0.334; uses diagonal errors only (fast).
    Returns (H0, M, chi2, n)."""
    n = len(zHEL)

    def chi2(theta):
        H0, M = theta
        if H0 < 30 or H0 > 200:
            return 1e30
        mu = distance_modulus_grid(zHEL, zCMB, H0, 0.334)
        return float(((mb - M - mu) ** 2 / dmb ** 2).sum())

    r = minimize(chi2, [H0_init, M_init], method="Nelder-Mead",
                 options=dict(xatol=1e-4, fatol=1e-3, maxiter=2000))
    return r.x[0], r.x[1], r.fun, n


def main():
    pp = load_pantheon_plus()
    arr = pp["data"]
    zHEL = pp["zHEL"]; zCMB = pp["zCMB"]
    mb = pp["mb"]; is_calib = pp["is_calibrator"]
    is_HF = arr["USED_IN_SH0ES_HF"].astype(bool)
    dmb = arr["m_b_corr_err_DIAG"]
    RA = arr["RA"]; DEC = arr["DEC"]

    # Reference: full HF fit
    H0_ref, M_ref, chi2_ref, n_ref = fit_HM(
        zHEL[is_HF], zCMB[is_HF], mb[is_HF], dmb[is_HF])
    print(f"Reference (HF subset, n={n_ref}, Om=0.334 fixed): "
          f"H0={H0_ref:.3f}  M={M_ref:.4f}  chi2={chi2_ref:.2f}")
    # Don't expect this to match SH0ES exactly because we use diag errors and
    # don't include the Cepheid anchor - this is just a relative reference.
    # NOTE: M and H0 are degenerate without anchor; minimize will pick a
    # combination. We subtract off this M when comparing subsamples.

    # ============================================================
    # 1. H0(z) running (Krishnan 2020 style, but on HF SNe with calibrators
    #    fixed to global M)
    # ============================================================
    print("\n=== 1. H0(z) running ===")
    print("Fixing M from full Pantheon+SH0ES anchor; fit H0 in z bins.")
    # Anchor M from calibrators: M_anchor = mean(mb - CEPH_DIST) over calibrators
    M_anchor = np.mean(mb[is_calib] - arr["CEPH_DIST"][is_calib])
    print(f"  Anchored M = {M_anchor:.4f} mag")

    # bins
    z_edges = [0.01, 0.03, 0.06, 0.1, 0.2, 0.4, 0.7, 1.5]
    print(f"  {'z range':<14s} {'n':>5s} {'H0(z)':>8s} {'+/-':>6s} {'chi2/dof':>10s}")
    for i in range(len(z_edges) - 1):
        lo, hi = z_edges[i], z_edges[i + 1]
        m = (zCMB >= lo) & (zCMB < hi) & (~is_calib)
        if m.sum() < 5:
            continue
        # at fixed M_anchor, fit only H0
        def chi2_H0(H0):
            mu = distance_modulus_grid(zHEL[m], zCMB[m], H0[0], 0.334)
            return float(((mb[m] - M_anchor - mu) ** 2 / dmb[m] ** 2).sum())
        r = minimize(chi2_H0, [73.0], method="Nelder-Mead",
                     options=dict(xatol=1e-4, fatol=1e-3))
        # 1-sigma
        H0_b = r.x[0]
        # bracket
        chi_lo = chi2_H0(np.array([H0_b - 5]))
        chi_hi = chi2_H0(np.array([H0_b + 5]))
        # Numeric Δχ²=1 width
        from scipy.optimize import brentq
        def df(dH):
            return chi2_H0(np.array([H0_b + dH])) - r.fun - 1.0
        try:
            sig = abs(brentq(df, 0.01, 10))
        except Exception:
            sig = float('nan')
        print(f"  [{lo:5.3f},{hi:5.3f})  {m.sum():5d}  {H0_b:8.3f}  ±{sig:5.3f}  "
              f"{r.fun:8.2f}/{m.sum()-1}")

    # ============================================================
    # 2. Anisotropic H0 — RA/DEC hemispheres + CMB-dipole projection
    # ============================================================
    print("\n=== 2. Anisotropic H0 (RA/DEC/CMB-dipole hemispheres) ===")
    # Use only Hubble-flow SNe so the calibrator weights don't dominate
    m_HF = is_HF & (~is_calib)
    # Convert to galactic-like cartesian for dipole projection
    ra_rad = np.deg2rad(RA[m_HF])
    dec_rad = np.deg2rad(DEC[m_HF])
    # CMB dipole direction (galactic l,b ≈ 264,48 → ICRS RA,DEC ≈ 167, -7)
    cmb_RA_deg, cmb_DEC_deg = 167.94, -6.94
    cmb_ra = np.deg2rad(cmb_RA_deg)
    cmb_dec = np.deg2rad(cmb_DEC_deg)
    # cosine of angle between SN and CMB dipole
    cos_th = (np.sin(dec_rad) * np.sin(cmb_dec)
              + np.cos(dec_rad) * np.cos(cmb_dec) * np.cos(ra_rad - cmb_ra))

    splits = {
        "RA<180 (E hemisphere)":  RA[m_HF] < 180,
        "RA>=180 (W hemisphere)": RA[m_HF] >= 180,
        "DEC>0 (N gal)":          DEC[m_HF] > 0,
        "DEC<=0 (S gal)":         DEC[m_HF] <= 0,
        "Toward CMB dipole":      cos_th > 0,
        "Away from CMB dipole":   cos_th < 0,
    }
    print(f"  {'subset':<26s} {'n':>5s} {'H0':>8s} {'M':>8s} {'chi2':>8s}")
    for label, mask in splits.items():
        idx = np.where(m_HF)[0][mask]
        H0_s, M_s, chi2_s, n_s = fit_HM(
            zHEL[idx], zCMB[idx], mb[idx], dmb[idx])
        print(f"  {label:<26s} {n_s:5d}  {H0_s:8.3f}  {M_s:8.4f}  {chi2_s:8.2f}")

    # ============================================================
    # 3. Residual correlations — look at every numerical column
    # ============================================================
    print("\n=== 3. Pearson r between residual and other columns ===")
    H0_full = 73.0  # use SH0ES value as reference
    mu_pred = distance_modulus_grid(zHEL, zCMB, H0_full, 0.334)
    resid = mb - M_anchor - mu_pred
    # Only Hubble-flow + non-calibrators for clean test
    m = m_HF
    print(f"  Using {m.sum()} HF-non-calibrator SNe, residual = m_b - M_anchor"
          f" - mu(H0=73, Om=0.334).  RMS = {resid[m].std():.3f} mag")
    cols_to_test = ["c", "x1", "HOST_LOGMASS", "MWEBV", "VPEC", "zHD",
                    "biasCor_m_b", "biasCorErr_m_b", "FITPROB",
                    "m_b_corr_err_RAW", "m_b_corr_err_VPEC", "PKMJDERR"]
    for col in cols_to_test:
        if col not in arr.dtype.names:
            continue
        x = arr[col][m]
        good = np.isfinite(x) & np.isfinite(resid[m]) & (x > -900)
        if good.sum() < 30:
            continue
        x = x[good]
        r_pearson = np.corrcoef(x, resid[m][good])[0, 1]
        print(f"  {col:<24s}  n={good.sum():4d}  r = {r_pearson:+.4f}")

    # ============================================================
    # 4. Survey-by-survey H0
    # ============================================================
    print("\n=== 4. Per-survey H0 (HF only, non-calibrator) ===")
    surveys = {
        4:  "SDSS",  5:  "SNLS", 51: "SNLS-redo", 53: "PS1MD",
        61: "Foundation", 100: "SuperNova Legacy",
        # Pantheon+ uses many ID codes; let's just enumerate existing ones
    }
    ids = arr["IDSURVEY"]
    unique_ids = np.unique(ids[m_HF])
    print(f"  {'IDSURVEY':<10s} {'n':>5s} {'H0':>8s} {'M':>8s}")
    for s_id in unique_ids:
        mm = m_HF & (ids == s_id)
        if mm.sum() < 10:
            continue
        H0_s, M_s, _, n_s = fit_HM(zHEL[mm], zCMB[mm], mb[mm], dmb[mm])
        print(f"  {s_id:<10d}  {n_s:5d}  {H0_s:8.3f}  {M_s:8.4f}")


if __name__ == "__main__":
    main()
