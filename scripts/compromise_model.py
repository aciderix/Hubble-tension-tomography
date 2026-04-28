"""Contrarian compromise model: two small concurrent systematics.

Hypothesis (contre-pied of orthodoxy EDE / Calibration / Late-time):
the H0 tension reflects TWO sub-publication-grade systematics that the
SH0ES camp and Planck camp respectively dismiss in isolation:

  C1: small Cepheid photometric crowding bias dM_cep ~ +0.04 mag
      (Mortsell-Goobar-Friedman 2022, rejected by Riess+ 2022 but the
      rejection only covers nearest hosts; high-z host crowding is
      observationally hard to constrain at this level)
  C2: small varying-me-equivalent r_d shift, dr_d/r_d ~ -0.015
      (Hart-Chluba 2020, Sekiguchi-Takahashi 2021; under-tested because
      eclipsed by EDE in the literature; predicts NO sigma_8 enhancement
      unlike pure EDE)

Combined target: H0 ~ 70 km/s/Mpc -- the TRGB / Freedman CCHP value that
both 67-camp and 73-camp avoid, and which sits at the median of every
truly independent measurement.

This script tests the model with PHYSICALLY-MOTIVATED narrow priors:
  Gaussian dM_cep = 0.00 +/- 0.05 mag  (centered on no-bias, allow physical
                                         envelope per Mortsell+22)
  Gaussian dr_d/r_d = 0.00 +/- 0.02    (centered on no-shift, allow Hart-
                                         Chluba envelope)
This penalizes large excursions, so M_compromise must IMPROVE the fit
over M0 LCDM through MODEST values of both parameters, not by maxing
out a single one (which would mimic M1 or M2).

Compare with M0, M1, M2 already-run nautilus evidence.
"""
import sys
from pathlib import Path
import numpy as np
from nautilus import Prior, Sampler

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.fast_lh import HubbleLH
from src.cosmo import rd_planck_anchored
from src.cmb_compressed import chi2_cmb_with_rd
from src.cosmic_chronometers import chi2_cc_full

LH = HubbleLH()


def loglike_compromise(theta_dict):
    H0, Om, M, obh2 = (theta_dict["H0"], theta_dict["Om"],
                       theta_dict["M"], theta_dict["obh2"])
    dM = theta_dict["dMcep"]
    frd = theta_dict["frd"]  # rd = rd_LCDM * frd

    rd_lcdm = rd_planck_anchored(Om * (H0 / 100) ** 2, obh2)
    rd_eff = rd_lcdm * frd

    # SN with calibrators offset by dM
    chi2 = LH.chi2_pp(H0, Om, M, M_calib=M + dM)
    chi2 += LH.chi2_bao(H0, Om, rd_eff)
    chi2 += chi2_cmb_with_rd(H0, Om, obh2, rd_eff)
    chi2 += chi2_cc_full(H0, Om)

    # Physical Gaussian priors on the two systematic parameters
    chi2 += (dM / 0.05) ** 2          # crowding envelope: 1-sigma = 0.05 mag
    chi2 += ((frd - 1.0) / 0.02) ** 2  # me-variation envelope: 1-sigma = 2%

    return -0.5 * chi2


def make_prior():
    p = Prior()
    p.add_parameter("H0", dist=(50.0, 90.0))
    p.add_parameter("Om", dist=(0.10, 0.50))
    p.add_parameter("M", dist=(-20.0, -18.5))
    p.add_parameter("obh2", dist=(0.018, 0.027))
    # Wide top-hat prior; the Gaussian penalty in loglike does the physical work
    p.add_parameter("dMcep", dist=(-0.20, 0.20))
    p.add_parameter("frd",   dist=(0.94, 1.06))
    return p


def main():
    prior = make_prior()
    print("=== Sampling M_compromise (Cepheid crowding + me-variation, narrow priors) ===")
    print("Priors: dMcep N(0, 0.05),  frd = rd/rd_LCDM N(1, 0.02)")
    s = Sampler(prior, loglike_compromise, n_live=600, pool=4, n_networks=4, seed=42)
    s.run(verbose=False, discard_exploration=True)
    points, log_w, log_l = s.posterior()
    log_Z = s.log_z
    idx = np.argmax(log_l)
    map_p = points[idx]
    map_chi2 = -2 * log_l[idx]

    keys = list(prior.keys)
    print(f"\n  log Z = {log_Z:.3f}")
    print(f"  MAP χ² (incl. priors) = {map_chi2:.2f}")
    print(f"  MAP point:")
    for k, v in zip(keys, map_p):
        print(f"     {k:8s} = {v:.4f}")

    # Effective rd at MAP
    Hm, Om, M, obh2, dM, frd = map_p
    rd_eff = rd_planck_anchored(Om * (Hm / 100) ** 2, obh2) * frd
    print(f"     rd_eff   = {rd_eff:.2f} Mpc")
    print(f"     dM_cep   = {dM:+.4f} mag  ({dM/0.05:+.2f}σ in physical prior)")
    print(f"     dr_d/r_d = {frd-1:+.4f}     ({(frd-1)/0.02:+.2f}σ in physical prior)")

    # Marginal H0
    weights = np.exp(log_w - log_w.max())
    weights /= weights.sum()
    h0 = points[:, keys.index("H0")]
    h0_mean = (h0 * weights).sum()
    h0_std = np.sqrt(((h0 - h0_mean) ** 2 * weights).sum())
    print(f"     H0 marginal: {h0_mean:.3f} ± {h0_std:.3f}")

    # Marginal dM and frd
    dM_arr = points[:, keys.index("dMcep")]
    frd_arr = points[:, keys.index("frd")]
    dM_m = (dM_arr * weights).sum(); dM_s = np.sqrt(((dM_arr - dM_m) ** 2 * weights).sum())
    frd_m = (frd_arr * weights).sum(); frd_s = np.sqrt(((frd_arr - frd_m) ** 2 * weights).sum())
    print(f"     dM_cep marg: {dM_m:+.4f} ± {dM_s:.4f}  ({dM_m/dM_s:+.2f}σ from 0)")
    print(f"     frd-1 marg : {frd_m-1:+.4f} ± {frd_s:.4f}  ({(frd_m-1)/frd_s:+.2f}σ from 0)")

    print("\n=== Comparison vs reference models (from previous nautilus runs) ===")
    refs = {
        "M0 ΛCDM (rd consistent)":      (-804.78, "68.97 ± 0.30"),
        "M1 free rd (early new physics)": (-795.16, "73.45 ± 0.97"),
        "M2 ΔM_cep (Cepheid systematic)": (-794.37, "68.51 ± 0.31"),
        "M_compromise (this work)":      (log_Z,   f"{h0_mean:.2f} ± {h0_std:.2f}"),
    }
    print(f"  {'Model':<35s}  {'log Z':>10s}  {'ln B vs M0':>12s}  {'H0':>16s}")
    z0 = refs["M0 ΛCDM (rd consistent)"][0]
    for name, (z, h0_str) in refs.items():
        print(f"  {name:<35s}  {z:>10.3f}  {z-z0:>+12.3f}  {h0_str:>16s}")

    # Comparison vs M2
    z2 = refs["M2 ΔM_cep (Cepheid systematic)"][0]
    print(f"\n  ln B(M_compromise / M2) = {log_Z - z2:+.3f}")
    print(f"  ln B(M_compromise / M1) = {log_Z - refs['M1 free rd (early new physics)'][0]:+.3f}")

    # Save points for later cross-checks
    out_path = Path(__file__).resolve().parents[1] / "out_compromise.npz"
    np.savez(out_path, points=points, log_w=log_w, log_l=log_l,
             keys=np.array(keys, dtype="U10"), log_Z=log_Z)
    print(f"\n  Saved posterior to {out_path}")


if __name__ == "__main__":
    main()
