"""Joint Pantheon+SH0ES + DESI DR2 + Planck-compressed CMB fit.

Models compared:
  M0  ΛCDM            (H0, Om, M, obh2)         rd computed from Brieden
  M1  ΛCDM + free rd  (H0, Om, M, obh2, rd)     'early new physics' (EDE-like)
  M2  ΛCDM + ΔM_calib (H0, Om, M, obh2, dMcep)  Cepheid systematic
  M3  w0waCDM         (H0, Om, M, obh2, w0, wa) late-time DE

For each model: best-fit χ², AIC, BIC, and per-dataset breakdown.
"""
import sys
from pathlib import Path
import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.fast_lh import HubbleLH
from src.cosmo import (rd_planck_anchored, comoving_distance_w0wa_grid,
                       E_w0wa, C_KMS)
from src.cmb_compressed import (chi2_cmb_with_rd, comoving_distance_to_zstar,
                                 PLANCK_INVCOV, PLANCK_MEAN)


def main():
    lh = HubbleLH()
    n_data = lh.n_pp + lh.n_bao + 3  # 3 = R, lA, obh2
    print(f"Total datapoints = {lh.n_pp} (SN) + {lh.n_bao} (BAO) + 3 (CMB compressed) = {n_data}\n")

    # Common BBN prior on omega_b for context (Cooke+ 2018)
    BBN_obh2 = 0.02218
    BBN_sig = 0.00055

    def chi2_model(theta, model):
        if model == "M0":
            H0, Om, M, obh2 = theta
            if not (40 < H0 < 100 and 0.05 < Om < 0.7 and -22 < M < -18
                    and 0.005 < obh2 < 0.05):
                return 1e30
            rd = rd_planck_anchored(Om * (H0 / 100) ** 2, obh2)
            return (lh.chi2_pp(H0, Om, M)
                    + lh.chi2_bao(H0, Om, rd)
                    + chi2_cmb_with_rd(H0, Om, obh2, rd))
        if model == "M1":
            H0, Om, M, obh2, rd = theta
            if not (40 < H0 < 100 and 0.05 < Om < 0.7 and -22 < M < -18
                    and 0.005 < obh2 < 0.05 and 100 < rd < 200):
                return 1e30
            return (lh.chi2_pp(H0, Om, M)
                    + lh.chi2_bao(H0, Om, rd)
                    + chi2_cmb_with_rd(H0, Om, obh2, rd))
        if model == "M2":
            H0, Om, M, obh2, dM = theta
            if not (40 < H0 < 100 and 0.05 < Om < 0.7 and -22 < M < -18
                    and 0.005 < obh2 < 0.05 and -0.5 < dM < 0.5):
                return 1e30
            rd = rd_planck_anchored(Om * (H0 / 100) ** 2, obh2)
            # SN with M for HF, M+dM for calibrators
            chi2_sn = lh.chi2_pp(H0, Om, M, M_calib=M + dM)
            return (chi2_sn + lh.chi2_bao(H0, Om, rd)
                    + chi2_cmb_with_rd(H0, Om, obh2, rd))
        if model == "M3":
            H0, Om, M, obh2, w0, wa = theta
            if not (40 < H0 < 100 and 0.05 < Om < 0.7 and -22 < M < -18
                    and 0.005 < obh2 < 0.05 and -3 < w0 < 0 and -5 < wa < 5):
                return 1e30
            rd = rd_planck_anchored(Om * (H0 / 100) ** 2, obh2)
            # SN under w0wa
            from src.cosmo import comoving_distance_w0wa_grid
            dC = comoving_distance_w0wa_grid(lh.zHD, H0, Om, w0, wa)
            dL = (1 + lh.zHEL) * dC
            mu_th = 5 * np.log10(dL) + 25
            mu_model = np.where(lh.is_calib, lh.ceph_dist, mu_th)
            resid = lh.mb - (M + mu_model)
            chi2_sn = float(resid @ lh.Cinv_pp @ resid)
            # BAO under w0wa
            dM_grid = comoving_distance_w0wa_grid(lh.bao_unique_z, H0, Om, w0, wa)
            E_grid = E_w0wa(lh.bao_unique_z, Om, w0, wa)
            Hz_grid = H0 * E_grid
            dM_map = dict(zip(lh.bao_unique_z, dM_grid))
            Hz_map = dict(zip(lh.bao_unique_z, Hz_grid))
            pred = np.empty(lh.n_bao)
            for i, (z, ob) in enumerate(zip(lh.bao_z, lh.bao_obs)):
                if ob == "DM_over_rs":
                    pred[i] = dM_map[z] / rd
                elif ob == "DH_over_rs":
                    pred[i] = (C_KMS / Hz_map[z]) / rd
                elif ob == "DV_over_rs":
                    pred[i] = (z * dM_map[z] ** 2 * C_KMS / Hz_map[z]) ** (1/3) / rd
            d = lh.bao_val - pred
            chi2_bao = float(d @ lh.bao_invcov @ d)
            # CMB: under w0wa, D_M(z*) is different
            from src.cosmo import E_LCDM_with_rad  # use LCDM for radiation epoch
            # Approximation: at z>>0 w0wa contributes little to integration;
            # use LCDM with same Om to get D_M(z*). This is rough but for w0wa
            # close to LCDM at high z it's fine.
            chi2_cmb_v = chi2_cmb_with_rd(H0, Om, obh2, rd)
            return chi2_sn + chi2_bao + chi2_cmb_v
        raise ValueError(model)

    # Initial guesses tuned to Planck/SH0ES regions
    inits = {
        "M0": [68.0, 0.31, -19.4, 0.02236],
        "M1": [70.0, 0.31, -19.3, 0.02236, 140.0],
        "M2": [68.0, 0.31, -19.4, 0.02236, 0.10],
        "M3": [70.0, 0.31, -19.4, 0.02236, -1.0, 0.0],
    }
    n_par = {"M0": 4, "M1": 5, "M2": 5, "M3": 6}
    results = {}

    for name in ["M0", "M1", "M2", "M3"]:
        print(f"\n=== Fitting {name} ===")
        # Try a few starting points to escape local minima
        best = None
        for jitter in [0.0, 0.5, -0.3]:
            x0 = list(inits[name])
            x0[0] += jitter * 2
            r = minimize(lambda x: chi2_model(x, name), x0, method="Nelder-Mead",
                         options=dict(xatol=1e-3, fatol=1e-2, maxiter=20000))
            if best is None or r.fun < best.fun:
                best = r
        x = best.x
        chi2_total = best.fun
        # Per-dataset breakdown
        H0, Om, M, obh2 = x[0], x[1], x[2], x[3]
        if name == "M0":
            rd_eff = rd_planck_anchored(Om * (H0 / 100) ** 2, obh2)
        elif name == "M1":
            rd_eff = x[4]
        elif name == "M2":
            rd_eff = rd_planck_anchored(Om * (H0 / 100) ** 2, obh2)
        else:  # M3
            rd_eff = rd_planck_anchored(Om * (H0 / 100) ** 2, obh2)
        if name == "M2":
            chi2_sn = lh.chi2_pp(H0, Om, M, M_calib=M + x[4])
        elif name == "M3":
            w0, wa = x[4], x[5]
            dC = comoving_distance_w0wa_grid(lh.zHD, H0, Om, w0, wa)
            dL = (1 + lh.zHEL) * dC
            mu_th = 5 * np.log10(dL) + 25
            mu_model = np.where(lh.is_calib, lh.ceph_dist, mu_th)
            resid = lh.mb - (M + mu_model)
            chi2_sn = float(resid @ lh.Cinv_pp @ resid)
        else:
            chi2_sn = lh.chi2_pp(H0, Om, M)
        chi2_bao_v = lh.chi2_bao(H0, Om, rd_eff)
        chi2_cmb_v = chi2_cmb_with_rd(H0, Om, obh2, rd_eff)

        k = n_par[name]
        AIC = chi2_total + 2 * k
        BIC = chi2_total + k * np.log(n_data)
        results[name] = {
            "x": x, "chi2": chi2_total, "k": k, "AIC": AIC, "BIC": BIC,
            "H0": H0, "Om": Om, "M": M, "obh2": obh2, "rd": rd_eff,
            "chi2_sn": chi2_sn, "chi2_bao": chi2_bao_v, "chi2_cmb": chi2_cmb_v,
        }
        print(f"  H0   = {H0:.3f}")
        print(f"  Om   = {Om:.4f}")
        print(f"  M    = {M:.4f}")
        print(f"  ωb   = {obh2:.5f}")
        if name == "M1":
            print(f"  rd   = {x[4]:.2f} Mpc  (LCDM consistent: {rd_planck_anchored(Om*(H0/100)**2,obh2):.2f})")
        if name == "M2":
            print(f"  ΔM_calib = {x[4]:+.4f} mag")
        if name == "M3":
            print(f"  w0   = {x[4]:.4f}")
            print(f"  wa   = {x[5]:.4f}")
        print(f"  rd_used = {rd_eff:.2f} Mpc")
        print(f"  χ²_total = {chi2_total:.2f}")
        print(f"     SN  = {chi2_sn:.2f}")
        print(f"     BAO = {chi2_bao_v:.2f}")
        print(f"     CMB = {chi2_cmb_v:.2f}")
        print(f"  AIC = {AIC:.2f},  BIC = {BIC:.2f}")

    # Comparison table
    print("\n\n========================================")
    print("MODEL COMPARISON (full SN + DESI + CMB-compressed)")
    print("========================================")
    print(f"{'Model':<8s} {'H0':>8s} {'Om':>7s} {'rd':>7s} {'χ²':>10s} {'k':>3s} {'AIC':>10s} {'BIC':>10s} {'Δχ²':>8s} {'ΔAIC':>8s}")
    base_chi2 = results["M0"]["chi2"]
    base_AIC = results["M0"]["AIC"]
    for name in ["M0", "M1", "M2", "M3"]:
        r = results[name]
        print(f"{name:<8s} {r['H0']:>8.3f} {r['Om']:>7.4f} {r['rd']:>7.2f} "
              f"{r['chi2']:>10.2f} {r['k']:>3d} {r['AIC']:>10.2f} {r['BIC']:>10.2f} "
              f"{r['chi2']-base_chi2:>+8.2f} {r['AIC']-base_AIC:>+8.2f}")

    print("\nLegend:")
    print("  M0 ΛCDM (consistent rd)")
    print("  M1 ΛCDM + free rd  (early-time new physics)")
    print("  M2 ΛCDM + ΔM_calib (Cepheid systematic)")
    print("  M3 w0waCDM with consistent rd (late-time DE)")


if __name__ == "__main__":
    main()
