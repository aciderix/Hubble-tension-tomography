"""Final test: Bayesian evidence comparison of M1 under TWO physical
interpretations of the rd shift, with KiDS-1000 S8 added to the likelihood.

  M1_EDE: rd shift achieved by Early Dark Energy. Standard literature
          claim (Hill+ 2020): EDE boosts sigma_8 by ~6% at f_EDE = 10%.
          Therefore S8_predicted = S8_geom_only x 1.06.

  M1_me:  rd shift achieved by varying electron mass at recombination
          (Hart-Chluba 2020). Recombination time shifts but no energy is
          injected, no boost to growth. S8_predicted = S8_geom_only.

  M_compromise: small Cepheid bias + small rd shift, with the rd shift
          purely from me-variation (no boost).

We add the KiDS-1000 weak-lensing S8 likelihood:
  S8_KiDS = 0.766 +/- 0.020

The Bayesian evidence comparison is now:
  - Does M1_EDE survive S8?  (literature says barely, ~2sigma penalty)
  - Does M1_me survive S8?   (CONTRARIAN PREDICTION: yes, no penalty)
  - Does M_compromise survive S8?
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
from src.growth import sigma8 as sigma8_calc

LH = HubbleLH()
SIG_DM = 0.05
SIG_FRD = 0.02

# Planck 2018 tight constraint on omega_m h^2 (matter physical density).
# This is the actual CMB measurement; the compressed (R, lA, omega_b)
# distance priors are too loose on omega_m, so we impose it explicitly.
OMH2_PLANCK = 0.1430
OMH2_PLANCK_SIG = 0.0011

# KiDS-1000 weak lensing
S8_KIDS = 0.766
S8_KIDS_SIG = 0.020


def S8_pred(H0, Om, obh2, ns, growth_boost=0.0):
    """Predict S8 = sigma8 * sqrt(Om/0.3).

    growth_boost = fractional sigma8 boost from the rd-shift physics.
    0.0 for me-variation, 0.06 for EDE-like.
    """
    h = H0 / 100
    s8 = sigma8_calc(Om, obh2, h, ns)
    s8 *= (1.0 + growth_boost)
    return s8 * np.sqrt(Om / 0.3)


def chi2_S8(H0, Om, obh2, ns=0.9649, growth_boost=0.0):
    S = S8_pred(H0, Om, obh2, ns, growth_boost=growth_boost)
    return ((S - S8_KIDS) / S8_KIDS_SIG) ** 2


def chi2_omh2(H0, Om):
    """Tight Planck 2018 prior on omega_m h^2 = 0.1430 +/- 0.0011."""
    omh2 = Om * (H0 / 100) ** 2
    return ((omh2 - OMH2_PLANCK) / OMH2_PLANCK_SIG) ** 2


def chi2_data(H0, Om, M, obh2, rd_eff, dM=None):
    M_calib = None if dM is None else M + dM
    return (LH.chi2_pp(H0, Om, M, M_calib=M_calib)
            + LH.chi2_bao(H0, Om, rd_eff)
            + chi2_cmb_with_rd(H0, Om, obh2, rd_eff)
            + chi2_cc_full(H0, Om)
            + chi2_omh2(H0, Om))


def loglike_M0(td):
    H0, Om, M, obh2 = td["H0"], td["Om"], td["M"], td["obh2"]
    rd = rd_planck_anchored(Om * (H0 / 100) ** 2, obh2)
    chi2 = chi2_data(H0, Om, M, obh2, rd)
    chi2 += chi2_S8(H0, Om, obh2)  # LCDM = no boost
    return -0.5 * chi2


def loglike_M1_EDE(td):
    H0, Om, M, obh2 = td["H0"], td["Om"], td["M"], td["obh2"]
    frd = td["frd"]
    rd_lcdm = rd_planck_anchored(Om * (H0 / 100) ** 2, obh2)
    rd = rd_lcdm * frd
    # EDE-style: sigma8 boost ∝ |1-frd| / 0.07 * 0.06 (approx Hill+2020)
    boost = abs(1 - frd) / 0.07 * 0.06
    chi2 = chi2_data(H0, Om, M, obh2, rd)
    chi2 += chi2_S8(H0, Om, obh2, growth_boost=boost)
    chi2 += ((frd - 1.0) / SIG_FRD) ** 2  # physical prior on rd shift
    return -0.5 * chi2


def loglike_M1_me(td):
    H0, Om, M, obh2 = td["H0"], td["Om"], td["M"], td["obh2"]
    frd = td["frd"]
    rd_lcdm = rd_planck_anchored(Om * (H0 / 100) ** 2, obh2)
    rd = rd_lcdm * frd
    chi2 = chi2_data(H0, Om, M, obh2, rd)
    chi2 += chi2_S8(H0, Om, obh2, growth_boost=0.0)  # me-var: no boost
    chi2 += ((frd - 1.0) / SIG_FRD) ** 2
    return -0.5 * chi2


def loglike_M2(td):
    H0, Om, M, obh2 = td["H0"], td["Om"], td["M"], td["obh2"]
    dM = td["dMcep"]
    rd = rd_planck_anchored(Om * (H0 / 100) ** 2, obh2)
    chi2 = chi2_data(H0, Om, M, obh2, rd, dM=dM)
    chi2 += chi2_S8(H0, Om, obh2)  # no rd shift, no boost
    chi2 += (dM / SIG_DM) ** 2
    return -0.5 * chi2


def loglike_Mco_me(td):
    """Compromise with rd shift via me-variation (no growth boost)."""
    H0, Om, M, obh2 = td["H0"], td["Om"], td["M"], td["obh2"]
    dM = td["dMcep"]
    frd = td["frd"]
    rd_lcdm = rd_planck_anchored(Om * (H0 / 100) ** 2, obh2)
    rd = rd_lcdm * frd
    chi2 = chi2_data(H0, Om, M, obh2, rd, dM=dM)
    chi2 += chi2_S8(H0, Om, obh2, growth_boost=0.0)
    chi2 += (dM / SIG_DM) ** 2
    chi2 += ((frd - 1.0) / SIG_FRD) ** 2
    return -0.5 * chi2


def loglike_Mco_EDE(td):
    """Compromise with rd shift via EDE (with growth boost)."""
    H0, Om, M, obh2 = td["H0"], td["Om"], td["M"], td["obh2"]
    dM = td["dMcep"]
    frd = td["frd"]
    rd_lcdm = rd_planck_anchored(Om * (H0 / 100) ** 2, obh2)
    rd = rd_lcdm * frd
    boost = abs(1 - frd) / 0.07 * 0.06
    chi2 = chi2_data(H0, Om, M, obh2, rd, dM=dM)
    chi2 += chi2_S8(H0, Om, obh2, growth_boost=boost)
    chi2 += (dM / SIG_DM) ** 2
    chi2 += ((frd - 1.0) / SIG_FRD) ** 2
    return -0.5 * chi2


def make_prior(extra_keys=None):
    p = Prior()
    p.add_parameter("H0", dist=(50.0, 90.0))
    p.add_parameter("Om", dist=(0.10, 0.50))
    p.add_parameter("M", dist=(-20.0, -18.5))
    p.add_parameter("obh2", dist=(0.018, 0.027))
    if extra_keys is None:
        extra_keys = []
    if "frd" in extra_keys:
        p.add_parameter("frd", dist=(0.92, 1.08))
    if "dMcep" in extra_keys:
        p.add_parameter("dMcep", dist=(-0.20, 0.20))
    return p


def run(name, ll, prior, n_live=400):
    print(f"\n=== {name} ===")
    s = Sampler(prior, ll, n_live=n_live, pool=4, n_networks=4, seed=42)
    s.run(verbose=False, discard_exploration=True)
    points, log_w, log_l = s.posterior()
    log_Z = s.log_z
    keys = list(prior.keys)
    weights = np.exp(log_w - log_w.max()); weights /= weights.sum()
    h0 = points[:, keys.index("H0")]
    h0m = (h0 * weights).sum()
    h0s = np.sqrt(((h0 - h0m) ** 2 * weights).sum())
    om = points[:, keys.index("Om")]
    obh2 = points[:, keys.index("obh2")]
    om_m = (om * weights).sum()
    obh2_m = (obh2 * weights).sum()
    # Compute S8 at posterior mean
    s8_mean = sigma8_calc(om_m, obh2_m, h0m / 100, 0.9649)
    S8_mean = s8_mean * np.sqrt(om_m / 0.3)
    print(f"  log Z = {log_Z:.3f}")
    print(f"  H0 = {h0m:.3f} ± {h0s:.3f}")
    print(f"  Om = {om_m:.4f}  obh2 = {obh2_m:.5f}")
    print(f"  S8 at posterior mean (no rd boost) = {S8_mean:.4f}")
    return {"log_Z": log_Z, "H0": (h0m, h0s), "Om": om_m, "obh2": obh2_m,
            "S8_mean": S8_mean}


def main():
    print("Datasets: SN+DESI BAO+CMB-compressed+CC+KiDS-1000 S8")
    print(f"S8_KiDS = {S8_KIDS} ± {S8_KIDS_SIG}")
    print(f"Priors: dMcep ~ N(0, {SIG_DM}), frd ~ N(1, {SIG_FRD})")

    res = {}
    res["M0"]      = run("M0  LCDM",                                   loglike_M0,     make_prior())
    res["M1_EDE"]  = run("M1 rd-shift via EDE physics (sigma8 boost)", loglike_M1_EDE, make_prior(["frd"]))
    res["M1_me"]   = run("M1 rd-shift via me-variation (no boost)",    loglike_M1_me,  make_prior(["frd"]))
    res["M2"]      = run("M2 Cepheid systematic (no rd, no boost)",    loglike_M2,     make_prior(["dMcep"]))
    res["Mco_me"]  = run("Mco compromise via me-variation",            loglike_Mco_me, make_prior(["frd", "dMcep"]))
    res["Mco_EDE"] = run("Mco compromise via EDE",                     loglike_Mco_EDE, make_prior(["frd", "dMcep"]))

    z0 = res["M0"]["log_Z"]
    print("\n\n=== Bayesian comparison (with KiDS S8 included) ===")
    print(f"{'Model':<40s} {'log Z':>10s} {'lnB vs M0':>12s} {'H0':>10s} {'S8':>8s}")
    for name, label in [
        ("M0",      "M0  LCDM"),
        ("M1_EDE",  "M1 rd via EDE physics"),
        ("M1_me",   "M1 rd via me-variation"),
        ("M2",      "M2 Cepheid bias"),
        ("Mco_me",  "Mco compromise (me-variation)"),
        ("Mco_EDE", "Mco compromise (EDE)"),
    ]:
        r = res[name]
        h0m, h0s = r["H0"]
        print(f"{label:<40s} {r['log_Z']:>10.3f} {r['log_Z']-z0:>+12.3f} {h0m:>5.2f}±{h0s:.2f} {r['S8_mean']:>8.4f}")

    print("\nKey Bayes factors:")
    print(f"  ln B(M1_me / M1_EDE)   = {res['M1_me']['log_Z']  - res['M1_EDE']['log_Z']:+.3f}  (CONTRARIAN PREDICTION)")
    print(f"  ln B(Mco_me / Mco_EDE) = {res['Mco_me']['log_Z'] - res['Mco_EDE']['log_Z']:+.3f}")
    print(f"  ln B(M1_me / M2)       = {res['M1_me']['log_Z']  - res['M2']['log_Z']:+.3f}")
    print(f"  ln B(Mco_me / M2)      = {res['Mco_me']['log_Z'] - res['M2']['log_Z']:+.3f}")


if __name__ == "__main__":
    main()
