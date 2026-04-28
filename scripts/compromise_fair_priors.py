"""Fair Bayesian model comparison with PHYSICALLY MOTIVATED priors on
ALL extension parameters (not just M_compromise).

Without informative priors on extension parameters, the original nautilus
run gave artificially high evidence to M1 and M2 because their parameter
spaces are big and the data picks out a specific point with high
likelihood. Adding informative physical priors makes the comparison fair:
each extension parameter must come from a physically-motivated envelope.

Physical envelopes:
  Δr_d / r_d   ~ N(0, 0.02):
      EDE-class models in literature give 5-10% rd shifts at f_EDE = 10%;
      a 2σ envelope on rd at the 4% level is realistic given current
      CMB damping-tail constraints. Wider would be unphysical.

  ΔM_cep        ~ N(0, 0.05) mag:
      Mortsell-Goobar-Friedman 2022 estimate worst-case crowding bias
      0.05 mag in distant hosts; this is the 1σ envelope.

We re-fit M0, M1, M2, M_compromise with these priors:

  M0    : LCDM consistent, no extensions
  M1    : free rd with N(rd_LCDM, 2%) Gaussian prior
  M2    : free ΔM_cep with N(0, 0.05) Gaussian prior
  M_co  : both, both N priors

Each datum N(...) prior here is the SAME physical envelope; the test
is which structure (one big effect, two modest effects, or no effect)
the data prefer.
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
SIG_DM = 0.05      # Cepheid crowding envelope (mag)
SIG_FRD = 0.02     # rd shift envelope (fractional)


def chi2_data(H0, Om, M, obh2, rd_eff, dM=None):
    M_calib = None if dM is None else M + dM
    return (LH.chi2_pp(H0, Om, M, M_calib=M_calib)
            + LH.chi2_bao(H0, Om, rd_eff)
            + chi2_cmb_with_rd(H0, Om, obh2, rd_eff)
            + chi2_cc_full(H0, Om))


def loglike_M0(td):
    H0, Om, M, obh2 = td["H0"], td["Om"], td["M"], td["obh2"]
    rd = rd_planck_anchored(Om * (H0 / 100) ** 2, obh2)
    return -0.5 * chi2_data(H0, Om, M, obh2, rd)


def loglike_M1(td):
    H0, Om, M, obh2 = td["H0"], td["Om"], td["M"], td["obh2"]
    frd = td["frd"]
    rd_lcdm = rd_planck_anchored(Om * (H0 / 100) ** 2, obh2)
    rd = rd_lcdm * frd
    chi2 = chi2_data(H0, Om, M, obh2, rd)
    chi2 += ((frd - 1.0) / SIG_FRD) ** 2  # physical prior
    return -0.5 * chi2


def loglike_M2(td):
    H0, Om, M, obh2 = td["H0"], td["Om"], td["M"], td["obh2"]
    dM = td["dMcep"]
    rd = rd_planck_anchored(Om * (H0 / 100) ** 2, obh2)
    chi2 = chi2_data(H0, Om, M, obh2, rd, dM=dM)
    chi2 += (dM / SIG_DM) ** 2
    return -0.5 * chi2


def loglike_Mco(td):
    H0, Om, M, obh2 = td["H0"], td["Om"], td["M"], td["obh2"]
    dM = td["dMcep"]
    frd = td["frd"]
    rd_lcdm = rd_planck_anchored(Om * (H0 / 100) ** 2, obh2)
    rd = rd_lcdm * frd
    chi2 = chi2_data(H0, Om, M, obh2, rd, dM=dM)
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
        p.add_parameter("frd", dist=(0.92, 1.08))   # 4σ wide
    if "dMcep" in extra_keys:
        p.add_parameter("dMcep", dist=(-0.20, 0.20))  # 4σ wide
    return p


def run(name, ll, prior, n_live=600, seed=42):
    print(f"\n=== {name} ===")
    s = Sampler(prior, ll, n_live=n_live, pool=4, n_networks=4, seed=seed)
    s.run(verbose=False, discard_exploration=True)
    points, log_w, log_l = s.posterior()
    log_Z = s.log_z
    keys = list(prior.keys)
    weights = np.exp(log_w - log_w.max()); weights /= weights.sum()
    h0 = points[:, keys.index("H0")]
    h0m = (h0 * weights).sum()
    h0s = np.sqrt(((h0 - h0m) ** 2 * weights).sum())
    idx = np.argmax(log_l)
    map_p = dict(zip(keys, points[idx]))
    map_chi2 = -2 * log_l[idx]
    print(f"  log Z = {log_Z:.3f}  MAP χ²(incl prior) = {map_chi2:.2f}")
    print(f"  H0 = {h0m:.3f} ± {h0s:.3f}")
    extras = {}
    for k in ("frd", "dMcep"):
        if k in keys:
            arr = points[:, keys.index(k)]
            mu = (arr * weights).sum()
            sg = np.sqrt(((arr - mu) ** 2 * weights).sum())
            extras[k] = (mu, sg)
            print(f"  {k}: {mu:+.4f} ± {sg:.4f}  ({mu/sg if sg > 0 else 0:+.2f}σ from 0)")
    print(f"  MAP point: {map_p}")
    return {"log_Z": log_Z, "H0": (h0m, h0s), "extras": extras, "MAP": map_p,
            "MAP_chi2": map_chi2}


def main():
    print(f"Physical priors:  Δr_d/r_d ~ N(0, {SIG_FRD}),  ΔM_cep ~ N(0, {SIG_DM}) mag")
    print(f"All four models share the same physical envelopes when their")
    print(f"extension parameters are free; otherwise extensions are zero.\n")

    res = {}
    res["M0"]  = run("M0  LCDM (no extensions)",         loglike_M0,  make_prior())
    res["M1"]  = run("M1  Δr_d/r_d ~ N(0, 2%)",          loglike_M1,  make_prior(["frd"]))
    res["M2"]  = run("M2  ΔM_cep ~ N(0, 0.05 mag)",      loglike_M2,  make_prior(["dMcep"]))
    res["Mco"] = run("Mco both, narrow physical priors", loglike_Mco, make_prior(["frd", "dMcep"]))

    print("\n\n=== Bayesian model comparison (consistent physical priors) ===")
    z0 = res["M0"]["log_Z"]
    print(f"  {'Model':<35s} {'log Z':>10s} {'ln B vs M0':>12s} {'H0':>16s}")
    for name, label in [
        ("M0",  "M0  LCDM"),
        ("M1",  "M1  rd shift only (early phys)"),
        ("M2",  "M2  ΔM_cep only (Cepheid bias)"),
        ("Mco", "Mco both, joint compromise"),
    ]:
        r = res[name]
        h0m, h0s = r["H0"]
        print(f"  {label:<35s} {r['log_Z']:>10.3f} {r['log_Z']-z0:>+12.3f} {h0m:>5.2f} ± {h0s:.2f}")

    print(f"\n  ln B(Mco/M1) = {res['Mco']['log_Z']-res['M1']['log_Z']:+.3f}")
    print(f"  ln B(Mco/M2) = {res['Mco']['log_Z']-res['M2']['log_Z']:+.3f}")
    print(f"  ln B(M2/M1)  = {res['M2']['log_Z']-res['M1']['log_Z']:+.3f}")

    print("\nInterpretation guide:")
    print("  ln B > 5  : strong  | > 2.5: moderate | < 1: inconclusive")


if __name__ == "__main__":
    main()
