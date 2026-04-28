"""Systematic test of candidate Hubble-tension resolutions.

Each test fits a different *minimal* extension and reports:
  - best-fit parameters
  - chi2 (Pantheon+SH0ES, DESI DR2, joint)
  - Δχ² vs the LCDM joint fit with rd FIXED to Planck
  - whether the fit is internally consistent (i.e. each dataset's chi2 is OK)

Falsification: if the reduction in joint chi2 < 9.49 for 4 DoF (95 %CL), the
modification is *not* a detected resolution.
"""
import sys
from pathlib import Path
import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.fast_lh import HubbleLH
from src.cosmo import (rd_aubourg, comoving_distance_w0wa, E_w0wa, C_KMS,
                       comoving_distance_grid, E_LCDM, distance_modulus_grid)


def opt(f, x0, bounds=None, maxiter=10000):
    return minimize(f, x0, method="Nelder-Mead",
                    options=dict(xatol=1e-5, fatol=1e-4, maxiter=maxiter))


def main():
    lh = HubbleLH()

    # ============================================================
    # A) BASELINE: LCDM with rd FIXED to Planck (147.05 Mpc)
    # ============================================================
    print("\n" + "=" * 70)
    print("A) LCDM joint fit — rd FIXED to Planck (147.05 Mpc).")
    print("   This is the 'no new physics' answer. Both ladders forced together.")
    print("=" * 70)
    rd_p = 147.05
    f = lambda x: lh.chi2_joint(x[0], x[1], x[2], rd_p)
    rA = opt(f, [70, 0.32, -19.3])
    H0_A, Om_A, M_A = rA.x
    chi2_pp_A = lh.chi2_pp(*rA.x)
    chi2_bao_A = lh.chi2_bao(H0_A, Om_A, rd_p)
    print(f"  H0={H0_A:.3f}  Om={Om_A:.4f}  M={M_A:.4f}  rd={rd_p}")
    print(f"  chi2: PP+SH0ES={chi2_pp_A:.2f}  DESI={chi2_bao_A:.2f}  "
          f"total={rA.fun:.2f}")
    chi2_baseline = rA.fun

    # Decompose: what is Pantheon+SH0ES alone preferring vs DESI alone?
    fa1 = lambda x: lh.chi2_pp(x[0], x[1], x[2])
    rA1 = opt(fa1, [73, 0.33, -19.25])
    print(f"  [PP+SH0ES alone -> H0={rA1.x[0]:.3f}, chi2={rA1.fun:.2f}]")

    # ============================================================
    # B) FREE rd (early-time sound-horizon shift, EDE-effective)
    # ============================================================
    print("\n" + "=" * 70)
    print("B) LCDM with FREE r_d  (effective early-time modification).")
    print("=" * 70)
    f = lambda x: lh.chi2_joint(x[0], x[1], x[2], x[3])
    rB = opt(f, [73, 0.31, -19.25, 140.0])
    H0_B, Om_B, M_B, rd_B = rB.x
    chi2_pp_B = lh.chi2_pp(H0_B, Om_B, M_B)
    chi2_bao_B = lh.chi2_bao(H0_B, Om_B, rd_B)
    h_B = H0_B / 100
    rd_pred_BBN_B = rd_aubourg(Om_B * h_B**2, 0.02218)
    print(f"  H0={H0_B:.3f}  Om={Om_B:.4f}  M={M_B:.4f}  rd={rd_B:.2f}")
    print(f"  chi2: PP+SH0ES={chi2_pp_B:.2f}  DESI={chi2_bao_B:.2f}  total={rB.fun:.2f}")
    print(f"  Δχ² vs (A): {rB.fun - chi2_baseline:+.2f}  (1 extra param: 95% if <-3.84)")
    print(f"  rd_inferred / rd_BBN = {rd_B/rd_pred_BBN_B:.4f}  ({100*(rd_B/rd_pred_BBN_B-1):+.2f}%)")
    print(f"  rd_inferred / rd_Planck = {rd_B/rd_p:.4f}  ({100*(rd_B/rd_p-1):+.2f}%)")
    print("  Verdict: this measures the *size* of pre-recombination physics needed.")

    # ============================================================
    # C) w0waCDM (late-time evolving DE) — DESI's own headline model
    # ============================================================
    print("\n" + "=" * 70)
    print("C) w0waCDM with rd FIXED to Planck (DESI 2024 headline).")
    print("=" * 70)
    # We need a w0wa joint chi2
    from src.cosmo import comoving_distance_w0wa_grid
    def chi2_pp_w0wa(H0, Om, M, w0, wa):
        dC = comoving_distance_w0wa_grid(lh.zHD, H0, Om, w0, wa)
        dL = (1 + lh.zHEL) * dC
        mu_th = 5 * np.log10(dL) + 25
        mu_model = np.where(lh.is_calib, lh.ceph_dist, mu_th)
        resid = lh.mb - (M + mu_model)
        return float(resid @ lh.Cinv_pp @ resid)

    def chi2_bao_w0wa(H0, Om, w0, wa, rd):
        unique_z = lh.bao_unique_z
        dM = comoving_distance_w0wa_grid(unique_z, H0, Om, w0, wa)
        Hz = H0 * E_w0wa(unique_z, Om, w0, wa)
        dM_map = dict(zip(unique_z, dM))
        Hz_map = dict(zip(unique_z, Hz))
        pred = np.empty(lh.n_bao)
        for i, (z, ob) in enumerate(zip(lh.bao_z, lh.bao_obs)):
            if ob == "DM_over_rs":
                pred[i] = dM_map[z] / rd
            elif ob == "DH_over_rs":
                pred[i] = (C_KMS / Hz_map[z]) / rd
            elif ob == "DV_over_rs":
                dv = (z * dM_map[z]**2 * C_KMS / Hz_map[z])**(1/3)
                pred[i] = dv / rd
        d = lh.bao_val - pred
        return float(d @ lh.bao_invcov @ d)

    def f_C(x):
        H0, Om, M, w0, wa = x
        if not (50 < H0 < 90 and 0.1 < Om < 0.6 and -22 < M < -18
                and -3 < w0 < 1 and -5 < wa < 5):
            return 1e30
        return chi2_pp_w0wa(H0, Om, M, w0, wa) + chi2_bao_w0wa(H0, Om, w0, wa, rd_p)
    rC = opt(f_C, [70, 0.32, -19.3, -1, 0])
    H0_C, Om_C, M_C, w0_C, wa_C = rC.x
    chi2_pp_C = chi2_pp_w0wa(H0_C, Om_C, M_C, w0_C, wa_C)
    chi2_bao_C = chi2_bao_w0wa(H0_C, Om_C, w0_C, wa_C, rd_p)
    print(f"  H0={H0_C:.3f}  Om={Om_C:.4f}  M={M_C:.4f}  w0={w0_C:.3f}  wa={wa_C:.3f}")
    print(f"  chi2: PP+SH0ES={chi2_pp_C:.2f}  DESI={chi2_bao_C:.2f}  total={rC.fun:.2f}")
    print(f"  Δχ² vs (A): {rC.fun - chi2_baseline:+.2f}  (2 extra params: 95% if <-5.99)")

    # ============================================================
    # D) Free rd + free w0wa (most flexible model)
    # ============================================================
    print("\n" + "=" * 70)
    print("D) w0waCDM + FREE rd (combined late+early).")
    print("=" * 70)
    def f_D(x):
        H0, Om, M, w0, wa, rd = x
        if not (50 < H0 < 90 and 0.1 < Om < 0.6 and -22 < M < -18
                and -3 < w0 < 1 and -5 < wa < 5 and 100 < rd < 200):
            return 1e30
        return chi2_pp_w0wa(H0, Om, M, w0, wa) + chi2_bao_w0wa(H0, Om, w0, wa, rd)
    rD = opt(f_D, [73, 0.31, -19.25, -1, 0, 140])
    print(f"  H0={rD.x[0]:.3f}  Om={rD.x[1]:.4f}  M={rD.x[2]:.4f}  "
          f"w0={rD.x[3]:.3f}  wa={rD.x[4]:.3f}  rd={rD.x[5]:.2f}")
    print(f"  total chi2={rD.fun:.2f}  Δχ² vs (A): {rD.fun - chi2_baseline:+.2f}")

    # ============================================================
    # E) Cepheid systematic: separate M for calibrators vs Hubble-flow
    # ============================================================
    print("\n" + "=" * 70)
    print("E) LCDM + r_d=147 + 'Cepheid offset' ΔM_calib (free M_calib).")
    print("=" * 70)
    def f_E(x):
        H0, Om, M_hf, M_cal = x
        if not (50 < H0 < 90 and 0.1 < Om < 0.6 and -22 < M_hf < -18
                and -22 < M_cal < -18):
            return 1e30
        return (lh.chi2_pp(H0, Om, M_hf, M_calib=M_cal)
                + lh.chi2_bao(H0, Om, rd_p))
    rE = opt(f_E, [70, 0.31, -19.4, -19.25])
    H0_E, Om_E, M_hf_E, M_cal_E = rE.x
    chi2_pp_E = lh.chi2_pp(H0_E, Om_E, M_hf_E, M_calib=M_cal_E)
    chi2_bao_E = lh.chi2_bao(H0_E, Om_E, rd_p)
    print(f"  H0={H0_E:.3f}  Om={Om_E:.4f}  M_HF={M_hf_E:.4f}  M_cal={M_cal_E:.4f}")
    print(f"  ΔM = M_cal - M_HF = {M_cal_E - M_hf_E:+.4f} mag")
    print(f"  chi2: PP+SH0ES={chi2_pp_E:.2f}  DESI={chi2_bao_E:.2f}  total={rE.fun:.2f}")
    print(f"  Δχ² vs (A): {rE.fun - chi2_baseline:+.2f}")
    print("  Interpretation: ΔM>0 means Cepheid distances overestimated.")

    # ============================================================
    # F) M(z) drift: M = M0 + alpha * z, free alpha (SN evolution)
    # ============================================================
    print("\n" + "=" * 70)
    print("F) LCDM + r_d=147 + linear SN-mag drift M(z) = M + alpha*z.")
    print("=" * 70)
    def chi2_pp_drift(H0, Om, M, alpha):
        mu_th = distance_modulus_grid(lh.zHEL, lh.zHD, H0, Om)
        mu_model = np.where(lh.is_calib, lh.ceph_dist, mu_th)
        # alpha multiplies HD-redshift; calibrators have z<0.011 mostly so small impact
        M_eff = M + alpha * lh.zHD
        resid = lh.mb - (M_eff + mu_model)
        return float(resid @ lh.Cinv_pp @ resid)

    def f_F(x):
        H0, Om, M, alpha = x
        if not (50 < H0 < 90 and 0.1 < Om < 0.6 and -22 < M < -18
                and -1 < alpha < 1):
            return 1e30
        return chi2_pp_drift(H0, Om, M, alpha) + lh.chi2_bao(H0, Om, rd_p)
    rF = opt(f_F, [70, 0.31, -19.4, 0.0])
    print(f"  H0={rF.x[0]:.3f}  Om={rF.x[1]:.4f}  M={rF.x[2]:.4f}  α={rF.x[3]:.4f} mag/z")
    print(f"  total chi2={rF.fun:.2f}  Δχ² vs (A): {rF.fun - chi2_baseline:+.2f}")

    # ============================================================
    # SUMMARY
    # ============================================================
    print("\n" + "#" * 70)
    print("# SUMMARY (lower Δχ² is better; need < -k for k extra params at 1σ)")
    print("#" * 70)
    print(f"  A) LCDM, rd=Planck            chi2={chi2_baseline:.2f}  (baseline)")
    print(f"  B) LCDM + free rd             chi2={rB.fun:.2f}  Δχ²={rB.fun-chi2_baseline:+.2f}  "
          f"(1 param)  rd*={rd_B:.1f}")
    print(f"  C) w0waCDM, rd=Planck          chi2={rC.fun:.2f}  Δχ²={rC.fun-chi2_baseline:+.2f}  "
          f"(2 params)")
    print(f"  D) w0waCDM + free rd           chi2={rD.fun:.2f}  Δχ²={rD.fun-chi2_baseline:+.2f}  "
          f"(3 params)")
    print(f"  E) Cepheid offset M_calib free chi2={rE.fun:.2f}  Δχ²={rE.fun-chi2_baseline:+.2f}  "
          f"(1 param)")
    print(f"  F) SN drift α(z)               chi2={rF.fun:.2f}  Δχ²={rF.fun-chi2_baseline:+.2f}  "
          f"(1 param)")
    print(f"\n  H0_SH0ES_alone = 73.42  vs  H0_inverse_ladder_rd147 = 68.66")
    print("  Resolution requires moving these into agreement.")


if __name__ == "__main__":
    main()
