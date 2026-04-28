"""Bayesian evidence for the four key models via nautilus nested sampling.

Models share priors on (H0, Om, M, omega_b); each adds different params.
Bayes factor B_ij = exp(log Z_i - log Z_j) tells us posterior odds for
model i vs j (assuming equal prior probabilities).

Reference scale: |ln B| > 5 = strong, |ln B| > 2.5 = moderate, < 1 = inconclusive.

Datasets:
  - Pantheon+SH0ES (1701 SNe + Cepheid anchors)
  - DESI DR2 BAO (13 measurements)
  - Planck 2018 distance priors (R, lA, omega_b)
  - Moresco cosmic chronometers (15 H(z) measurements) [NEW]
"""
import sys
from pathlib import Path
import numpy as np
from nautilus import Prior, Sampler

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.fast_lh import HubbleLH
from src.cosmo import (rd_planck_anchored, comoving_distance_w0wa_grid,
                       E_w0wa, C_KMS)
from src.cmb_compressed import chi2_cmb_with_rd
from src.cosmic_chronometers import chi2_cc_full as chi2_cc

LH = HubbleLH()


def loglike_M0(theta_dict):
    H0, Om, M, obh2 = (theta_dict["H0"], theta_dict["Om"],
                       theta_dict["M"], theta_dict["obh2"])
    rd = rd_planck_anchored(Om * (H0 / 100) ** 2, obh2)
    chi2 = (LH.chi2_pp(H0, Om, M)
            + LH.chi2_bao(H0, Om, rd)
            + chi2_cmb_with_rd(H0, Om, obh2, rd)
            + chi2_cc(H0, Om))
    return -0.5 * chi2


def loglike_M1(theta_dict):
    H0, Om, M, obh2, rd = (theta_dict["H0"], theta_dict["Om"],
                            theta_dict["M"], theta_dict["obh2"],
                            theta_dict["rd"])
    chi2 = (LH.chi2_pp(H0, Om, M)
            + LH.chi2_bao(H0, Om, rd)
            + chi2_cmb_with_rd(H0, Om, obh2, rd)
            + chi2_cc(H0, Om))
    return -0.5 * chi2


def loglike_M2(theta_dict):
    H0, Om, M, obh2, dM = (theta_dict["H0"], theta_dict["Om"],
                            theta_dict["M"], theta_dict["obh2"],
                            theta_dict["dMcep"])
    rd = rd_planck_anchored(Om * (H0 / 100) ** 2, obh2)
    chi2 = (LH.chi2_pp(H0, Om, M, M_calib=M + dM)
            + LH.chi2_bao(H0, Om, rd)
            + chi2_cmb_with_rd(H0, Om, obh2, rd)
            + chi2_cc(H0, Om))
    return -0.5 * chi2


def loglike_M5(theta_dict):
    """Combined: free rd AND free ΔM_calib."""
    H0, Om, M, obh2, rd, dM = (theta_dict["H0"], theta_dict["Om"],
                                theta_dict["M"], theta_dict["obh2"],
                                theta_dict["rd"], theta_dict["dMcep"])
    chi2 = (LH.chi2_pp(H0, Om, M, M_calib=M + dM)
            + LH.chi2_bao(H0, Om, rd)
            + chi2_cmb_with_rd(H0, Om, obh2, rd)
            + chi2_cc(H0, Om))
    return -0.5 * chi2


def make_prior(extra=None):
    p = Prior()
    p.add_parameter("H0", dist=(50.0, 90.0))
    p.add_parameter("Om", dist=(0.10, 0.50))
    p.add_parameter("M", dist=(-20.0, -18.5))
    p.add_parameter("obh2", dist=(0.018, 0.027))
    if extra:
        for k, rng in extra.items():
            p.add_parameter(k, dist=rng)
    return p


def run(name, ll, prior, n_live=400):
    print(f"\n=== Sampling {name} (n_live={n_live}) ===")
    s = Sampler(prior, ll, n_live=n_live, pool=4, n_networks=4, seed=42)
    s.run(verbose=False, discard_exploration=True)
    points, log_w, log_l = s.posterior()
    log_Z = s.log_z
    # MAP
    idx = np.argmax(log_l)
    map_p = points[idx]
    print(f"  log Z = {log_Z:.3f}")
    print(f"  MAP χ² = {-2*log_l[idx]:.2f}")
    print(f"  MAP point: {dict(zip(prior.keys, map_p))}")
    # Marginal H0
    weights = np.exp(log_w - log_w.max())
    weights /= weights.sum()
    h0_idx = prior.keys.index("H0")
    h0 = points[:, h0_idx]
    h0_mean = (h0 * weights).sum()
    h0_std = np.sqrt(((h0 - h0_mean) ** 2 * weights).sum())
    print(f"  H0 marginal: {h0_mean:.3f} ± {h0_std:.3f}")
    return {"log_Z": log_Z, "map": map_p, "keys": prior.keys,
            "H0": (h0_mean, h0_std), "points": points, "weights": weights,
            "log_l": log_l, "log_w": log_w}


def main():
    priors = {
        "M0": make_prior(),
        "M1": make_prior(extra={"rd": (110.0, 170.0)}),
        "M2": make_prior(extra={"dMcep": (-0.5, 0.5)}),
    }
    lls = {"M0": loglike_M0, "M1": loglike_M1, "M2": loglike_M2}

    res = {}
    for name in ["M0", "M1", "M2"]:
        res[name] = run(name, lls[name], priors[name], n_live=400)

    print("\n\n=== Bayesian comparison (with cosmic chronometers added) ===")
    print(f"{'Model':<8s} {'log Z':>10s}  {'ln B vs M0':>12s}  {'ln B vs M1':>12s}  {'H0':>12s}")
    z0 = res["M0"]["log_Z"]
    z1 = res["M1"]["log_Z"]
    for name in ["M0", "M1", "M2"]:
        r = res[name]
        h0m, h0s = r["H0"]
        print(f"{name:<8s} {r['log_Z']:>10.3f}  {r['log_Z']-z0:>+12.3f}  "
              f"{r['log_Z']-z1:>+12.3f}  {h0m:>6.2f}±{h0s:.2f}")

    print("\nInterpretation:")
    print("  ln B > 5  : strong preference")
    print("  ln B > 2.5: moderate preference")
    print("  ln B < 1  : inconclusive")


if __name__ == "__main__":
    main()
