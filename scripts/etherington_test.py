"""Distance-duality (Etherington reciprocity) test.

Standard GR + photon conservation gives D_L(z) = (1+z)^2 D_A(z), so
  η(z) ≡ D_L(z) / [(1+z) D_M(z)] = 1  for all z.

A breaking would point to:
  - photon-axion mixing  (η > 1, SNe dim)
  - cosmic opacity (dust between sources and observer)
  - modified gravity (non-metric)
  - photon number non-conservation

Method:
  - From SN: bin Pantheon+ HF SNe in z, compute mean (mb - M_anchor) -> mu_obs
    -> D_L^SN = 10^((mu_obs - 25)/5) Mpc.
  - From DESI BAO at the same z: D_M^BAO = (D_M/r_d)_obs × r_d_assumed.
  - Form η_z = D_L^SN / [(1+z) D_M^BAO].
  - Test η = 1.  Fit η(z) = η0 + η1 (z) to look for trend.
"""
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.data_loaders import load_pantheon_plus
from src.desi_bao import load_desi_dr2_all


def main():
    pp = load_pantheon_plus()
    arr = pp["data"]
    is_calib = pp["is_calibrator"]
    is_HF = arr["USED_IN_SH0ES_HF"].astype(bool)
    M_anchor = np.mean(pp["mb"][is_calib] - arr["CEPH_DIST"][is_calib])
    print(f"Calibrator-anchored M = {M_anchor:.4f} mag")

    bao = load_desi_dr2_all()

    # SN: use all non-calibrator SNe (z up to 2.26) since DESI z goes up to 2.33
    m_use = (~is_calib)
    z_SN = pp["zCMB"][m_use]
    mu_SN = pp["mb"][m_use] - M_anchor
    err_SN = arr["m_b_corr_err_DIAG"][m_use]

    # For each DESI z, find SNe within +/- delta z and compute weighted mean mu
    print("\n=== Etherington test ===\n")
    print("For each DESI BAO z, derive D_L from nearby SNe, D_M from BAO,")
    print("compute η = D_L / [(1+z) D_M]. With r_d assumed.\n")

    # Work in two scenarios: r_d = 147 (Planck) and r_d = 137 (free-rd best fit)
    for rd_assumed, label in [(147.05, "rd = 147 Mpc (Planck)"),
                              (137.10, "rd = 137 Mpc (joint best-fit)")]:
        print(f"--- {label} ---")
        unique_z = np.unique(bao["z"])
        rows = []
        for z in unique_z:
            # Find DM from BAO at this z
            mask = (bao["z"] == z) & (bao["obs"] == "DM_over_rs")
            if not mask.any():
                # Try DV
                mask = (bao["z"] == z) & (bao["obs"] == "DV_over_rs")
                if not mask.any():
                    continue
                DV = bao["val"][mask][0] * rd_assumed
                # We can't separate DM here; compute (1+z)D_M from D_V indirectly
                # Skip for now
                continue
            DM = bao["val"][mask][0] * rd_assumed
            # SN window
            dz = max(0.05, 0.1 * z)
            sn_mask = np.abs(z_SN - z) < dz
            if sn_mask.sum() < 5:
                continue
            # Inverse-variance weighted mean of mu = D_L SN
            mu_local = mu_SN[sn_mask]
            err_local = err_SN[sn_mask]
            w = 1.0 / err_local ** 2
            mu_mean = np.sum(mu_local * w) / w.sum()
            mu_err = 1.0 / np.sqrt(w.sum())
            DL_SN = 10 ** ((mu_mean - 25) / 5)  # Mpc
            DL_SN_err = 0.4 * np.log(10) * DL_SN * mu_err
            DL_BAO = (1 + z) * DM
            eta = DL_SN / DL_BAO
            eta_err = DL_SN_err / DL_BAO
            rows.append((z, DM, DL_SN, eta, eta_err, sn_mask.sum()))

        print(f"  {'z':>6s} {'D_M(BAO)':>10s} {'D_L(SN)':>10s} {'η':>8s} {'±':>6s} {'n_SN':>6s}")
        eta_arr = []
        eta_err_arr = []
        z_arr = []
        for z, DM, DL_SN, eta, eta_err, n in rows:
            print(f"  {z:6.3f}  {DM:10.2f}  {DL_SN:10.2f}  {eta:8.4f} ±{eta_err:6.4f}  {n:6d}")
            eta_arr.append(eta); eta_err_arr.append(eta_err); z_arr.append(z)
        eta_arr = np.array(eta_arr); eta_err_arr = np.array(eta_err_arr); z_arr = np.array(z_arr)
        # Test η = 1 + η1 z
        # Weighted least-squares fit
        # η(z) = a + b z
        A = np.vstack([np.ones_like(z_arr), z_arr]).T
        W = np.diag(1.0 / eta_err_arr ** 2)
        cov = np.linalg.inv(A.T @ W @ A)
        params = cov @ A.T @ W @ eta_arr
        err = np.sqrt(np.diag(cov))
        chi2 = np.sum(((eta_arr - A @ params) / eta_err_arr) ** 2)
        print(f"  fit η(z) = a + b z:")
        print(f"    a = {params[0]:.4f} ± {err[0]:.4f}  ({(params[0]-1)/err[0]:+.2f}σ from 1)")
        print(f"    b = {params[1]:.4f} ± {err[1]:.4f}  ({params[1]/err[1]:+.2f}σ from 0)")
        print(f"    chi2 = {chi2:.2f}, dof = {len(z_arr)-2}\n")


if __name__ == "__main__":
    main()
