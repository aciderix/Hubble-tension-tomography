"""H0 estimate using ZERO disputed anchor.

Combines:
  - BBN: Cooke+ 2018 omega_b h^2 = 0.02218 +/- 0.00055 (D/H from primordial
    deuterium abundance)
  - DESI DR2 BAO: 13 measurements (rd inferred from BBN+CDM model, not CMB)
  - Pantheon+ HF SNe (NO Cepheid calibrators)
  - Moresco cosmic chronometers (full covariance)
  - NO Planck CMB
  - NO SH0ES Cepheid anchor

If the tension is between SH0ES and Planck, this fit deliberately uses
NEITHER. Should give a clean intermediate value.

The full Pantheon+ likelihood with calibrators MASKED OUT requires marginalizing
over the (M, H0) degeneracy in the SN-only term. We do this analytically
following Conley+ 2011 / Brout+ 2022: chi2 = d^T Cinv d - (1^T Cinv d)^2
/ (1^T Cinv 1), which marginalizes over the global magnitude offset.
"""
import sys
from pathlib import Path
import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.fast_lh import HubbleLH
from src.cosmic_chronometers import chi2_cc_full
from src.cosmo import (rd_planck_anchored, distance_modulus_grid, E_LCDM,
                       comoving_distance_grid, C_KMS)


def main():
    lh = HubbleLH()

    # Prepare Pantheon+ HF non-calibrator subset and analytically marginalize M
    is_calib = lh.is_calib
    is_HF = lh.bao_z * 0  # placeholder; we need USED_IN_SH0ES_HF from arr
    # Get USED_IN_SH0ES_HF directly
    from src.data_loaders import load_pantheon_plus
    pp = load_pantheon_plus()
    is_HF = pp["data"]["USED_IN_SH0ES_HF"].astype(bool)
    mask = is_HF & ~is_calib
    sub_idx = np.where(mask)[0]
    n_sub = mask.sum()
    print(f"Pantheon+ HF subset (no calibrators): n = {n_sub}")
    cov_sub = lh.cov_pp[np.ix_(sub_idx, sub_idx)]
    Cinv_sub = np.linalg.inv(cov_sub)
    zHEL_sub = lh.zHEL[sub_idx]
    zHD_sub = lh.zHD[sub_idx]
    mb_sub = lh.mb[sub_idx]

    def chi2_sn_marg(H0, Om):
        """Pantheon+ HF chi2 analytically marginalized over global M."""
        mu = distance_modulus_grid(zHEL_sub, zHD_sub, H0, Om)
        d = mb_sub - mu  # = M_eff(z)
        # Marginalize: chi2_M = d^T Cinv d - (1^T Cinv d)^2 / (1^T Cinv 1)
        Cinv_d = Cinv_sub @ d
        ones = np.ones(n_sub)
        Cinv_one = Cinv_sub @ ones
        a = d @ Cinv_d
        b = ones @ Cinv_d
        c = ones @ Cinv_one
        return a - b * b / c

    # BBN prior (Cooke+ 2018)
    BBN_obh2 = 0.02218
    BBN_sig = 0.00055

    def chi2_bbn(obh2):
        return ((obh2 - BBN_obh2) / BBN_sig) ** 2

    # Joint anchor-free chi2
    def chi2_joint(theta):
        H0, Om, obh2 = theta
        if not (40 < H0 < 100 and 0.05 < Om < 0.7
                and 0.005 < obh2 < 0.05):
            return 1e30
        rd = rd_planck_anchored(Om * (H0 / 100) ** 2, obh2)
        return (chi2_sn_marg(H0, Om)
                + lh.chi2_bao(H0, Om, rd)
                + chi2_bbn(obh2)
                + chi2_cc_full(H0, Om))

    # Best fit
    print("\n=== Joint anchor-free fit (BBN + DESI BAO + SN-HF marginalised + CC-full) ===")
    best = None
    for H0_start in [65, 68, 70, 72, 74]:
        r = minimize(chi2_joint, [H0_start, 0.30, BBN_obh2], method="Nelder-Mead",
                     options=dict(xatol=1e-4, fatol=1e-3, maxiter=5000))
        if best is None or r.fun < best.fun:
            best = r
    H0_b, Om_b, obh2_b = best.x
    rd_b = rd_planck_anchored(Om_b * (H0_b / 100) ** 2, obh2_b)
    print(f"  H0       = {H0_b:.3f}")
    print(f"  Om       = {Om_b:.4f}")
    print(f"  ωb       = {obh2_b:.5f}")
    print(f"  rd (Brieden) = {rd_b:.2f} Mpc")
    print(f"  chi2_total = {best.fun:.2f}")
    # Per-dataset breakdown
    chi2_sn = chi2_sn_marg(H0_b, Om_b)
    chi2_bao_v = lh.chi2_bao(H0_b, Om_b, rd_b)
    chi2_bbn_v = chi2_bbn(obh2_b)
    chi2_cc_v = chi2_cc_full(H0_b, Om_b)
    print(f"    SN(HF)  = {chi2_sn:.2f}  (n = {n_sub})")
    print(f"    BAO     = {chi2_bao_v:.2f}  (n = 13)")
    print(f"    BBN     = {chi2_bbn_v:.2f}  (n = 1)")
    print(f"    CC      = {chi2_cc_v:.2f}  (n = 15)")

    # 1D profile chi2(H0) -> 1-sigma error
    print("\nProfile chi2(H0):")
    chi2_min = best.fun

    def prof(H0_t):
        r = minimize(lambda x: chi2_joint([H0_t, x[0], x[1]]),
                     [Om_b, obh2_b], method="Nelder-Mead",
                     options=dict(xatol=1e-4, fatol=1e-3))
        return r.fun

    H0_grid = np.arange(H0_b - 5, H0_b + 5.01, 0.5)
    print(f"  {'H0':>6s} {'Δχ²':>8s}")
    for h in H0_grid:
        print(f"  {h:6.2f} {prof(h)-chi2_min:>8.3f}")

    # Find 1σ width  Δχ² = 1
    from scipy.optimize import brentq
    def fH(dH):
        return prof(H0_b + dH) - chi2_min - 1.0
    try:
        sig_p = brentq(fH, 0.01, 5.0)
    except Exception:
        sig_p = float('nan')
    try:
        sig_m = -brentq(fH, -5.0, -0.01)
    except Exception:
        sig_m = float('nan')
    print(f"\n  H0 = {H0_b:.3f} +{sig_p:.3f} -{sig_m:.3f}  (1σ profile, no CMB, no SH0ES)")
    print()
    print("Comparison:")
    print("  Planck (CMB+lensing)         : 67.36 +/- 0.54")
    print("  ACT+DESI                     : 68.22 +/- 0.30")
    print("  DESI BAO + BBN (Aubourg-like): ~68.5  +/- 0.7")
    print("  This work (BBN+BAO+SN+CC)    : ", end="")
    print(f"{H0_b:.2f} +{sig_p:.2f} -{sig_m:.2f}")
    print("  SH0ES (Cepheid anchor)       : 73.04 +/- 1.04")


if __name__ == "__main__":
    main()
