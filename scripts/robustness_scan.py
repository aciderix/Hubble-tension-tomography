"""Strict robustness analysis: how does the conclusion depend on the
omega_m prior width?

Per ChatGPT-recommended methodology:
  - 3 omega_m priors: tight (Planck ΛCDM σ=0.001), intermediate (σ=0.005),
    loose / weakly informative (σ=0.02 or no prior)
  - 4 representative models: M0 LCDM, M1 me-var (rd-shift), M2 Cepheid bias,
    Mco compromise (me-var + Cepheid)
  - Decomposition tests:
    a) M2 with Cepheid prior OFF (only data drives the dM_cep)
    b) M2 + M1 with omega_m prior OFF (only compressed CMB constrains)
  - Convergence: repeat each best-vs-worst case with a different seed,
    report log Z stability.

For each (model, prior-width) combination report:
  - log Z (with statistical error from nautilus)
  - H0 marginal posterior mean & std
  - omega_m marginal posterior mean & std
  - S8 at MAP point (CAMB-validated externally)
  - max|frd-1| and max|dM_cep| posterior

Output: a structured table that lets us cleanly identify which conclusions
are robust vs prior-dependent.
"""
import sys
from pathlib import Path
import numpy as np
import json
from nautilus import Prior, Sampler

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.fast_lh import HubbleLH
from src.cosmo import rd_planck_anchored
from src.cmb_compressed import chi2_cmb_with_rd
from src.cosmic_chronometers import chi2_cc_full
from src.growth import sigma8 as sigma8_calc

LH = HubbleLH()

# Constants
S8_KIDS = 0.766
S8_KIDS_SIG = 0.020
OMH2_PLANCK = 0.1430


def chi2_omh2(H0, Om, sigma_omh2):
    if sigma_omh2 is None or sigma_omh2 <= 0:
        return 0.0
    omh2 = Om * (H0 / 100) ** 2
    return ((omh2 - OMH2_PLANCK) / sigma_omh2) ** 2


def chi2_S8(H0, Om, obh2, ns=0.9649, growth_boost=0.0, include=True):
    if not include:
        return 0.0
    h = H0 / 100
    s = sigma8_calc(Om, obh2, h, ns) * (1.0 + growth_boost)
    S = s * np.sqrt(Om / 0.3)
    return ((S - S8_KIDS) / S8_KIDS_SIG) ** 2


def make_loglike(model, *, sigma_omh2, sig_dM, sig_frd, include_S8, growth_boost):
    """Factory for likelihood functions with configurable nuisance priors."""
    def chi2_data(H0, Om, M, obh2, dM=None, rd_shift_factor=1.0):
        rd_lcdm = rd_planck_anchored(Om * (H0 / 100) ** 2, obh2)
        rd_eff = rd_lcdm * rd_shift_factor
        M_calib = None if dM is None else M + dM
        chi2 = (LH.chi2_pp(H0, Om, M, M_calib=M_calib)
                + LH.chi2_bao(H0, Om, rd_eff)
                + chi2_cmb_with_rd(H0, Om, obh2, rd_eff)
                + chi2_cc_full(H0, Om))
        chi2 += chi2_omh2(H0, Om, sigma_omh2)
        chi2 += chi2_S8(H0, Om, obh2, growth_boost=growth_boost, include=include_S8)
        return chi2

    if model == "M0":
        def ll(td):
            chi2 = chi2_data(td["H0"], td["Om"], td["M"], td["obh2"])
            return -0.5 * chi2
    elif model == "M1":
        def ll(td):
            chi2 = chi2_data(td["H0"], td["Om"], td["M"], td["obh2"],
                              rd_shift_factor=td["frd"])
            if sig_frd > 0:
                chi2 += ((td["frd"] - 1.0) / sig_frd) ** 2
            return -0.5 * chi2
    elif model == "M2":
        def ll(td):
            chi2 = chi2_data(td["H0"], td["Om"], td["M"], td["obh2"],
                              dM=td["dMcep"])
            if sig_dM > 0:
                chi2 += (td["dMcep"] / sig_dM) ** 2
            return -0.5 * chi2
    elif model == "Mco":
        def ll(td):
            chi2 = chi2_data(td["H0"], td["Om"], td["M"], td["obh2"],
                              dM=td["dMcep"], rd_shift_factor=td["frd"])
            if sig_dM > 0:
                chi2 += (td["dMcep"] / sig_dM) ** 2
            if sig_frd > 0:
                chi2 += ((td["frd"] - 1.0) / sig_frd) ** 2
            return -0.5 * chi2
    return ll


def make_prior(model):
    p = Prior()
    p.add_parameter("H0", dist=(50.0, 90.0))
    p.add_parameter("Om", dist=(0.10, 0.50))
    p.add_parameter("M", dist=(-20.0, -18.5))
    p.add_parameter("obh2", dist=(0.018, 0.027))
    if model in ("M1", "Mco"):
        p.add_parameter("frd", dist=(0.92, 1.08))
    if model in ("M2", "Mco"):
        p.add_parameter("dMcep", dist=(-0.20, 0.20))
    return p


def run_one(name, model, *, sigma_omh2, sig_dM, sig_frd, include_S8,
            growth_boost, n_live=400, seed=42):
    ll = make_loglike(model, sigma_omh2=sigma_omh2, sig_dM=sig_dM,
                      sig_frd=sig_frd, include_S8=include_S8,
                      growth_boost=growth_boost)
    prior = make_prior(model)
    s = Sampler(prior, ll, n_live=n_live, pool=4, n_networks=4, seed=seed)
    s.run(verbose=False, discard_exploration=True)
    pts, log_w, log_l = s.posterior()
    keys = list(prior.keys)
    log_Z = s.log_z
    weights = np.exp(log_w - log_w.max()); weights /= weights.sum()
    H0 = pts[:, keys.index("H0")]
    Om = pts[:, keys.index("Om")]
    h0m = (H0 * weights).sum()
    h0s = np.sqrt(((H0 - h0m) ** 2 * weights).sum())
    omm = (Om * weights).sum()
    omhs = (Om * (H0 / 100) ** 2 * weights).sum()
    omh2_std = np.sqrt(((Om * (H0/100)**2 - omhs) ** 2 * weights).sum())
    out = {"name": name, "model": model, "log_Z": float(log_Z),
           "H0": float(h0m), "H0_std": float(h0s),
           "Om": float(omm),
           "omh2": float(omhs), "omh2_std": float(omh2_std),
           "n_par": len(keys), "n_eff": int(s.effective_sample_size())}
    if "frd" in keys:
        frd = pts[:, keys.index("frd")]
        m = (frd * weights).sum(); s_ = np.sqrt(((frd-m)**2*weights).sum())
        out["frd"] = float(m); out["frd_std"] = float(s_)
    if "dMcep" in keys:
        dM = pts[:, keys.index("dMcep")]
        m = (dM * weights).sum(); s_ = np.sqrt(((dM-m)**2*weights).sum())
        out["dMcep"] = float(m); out["dMcep_std"] = float(s_)
    return out


def main():
    print("=== Robustness scan: 4 models x 3 omega_m priors + decomposition ===\n")

    OMH2_PRIORS = {
        "tight":   ("σ_omh2 = 0.001 (Planck-LCDM)", 0.001),
        "interm":  ("σ_omh2 = 0.005 (5x looser)",   0.005),
        "loose":   ("σ_omh2 = inf (compressed CMB only)", None),
    }

    SIG_DM_DEFAULT = 0.05
    SIG_FRD_DEFAULT = 0.02
    INCLUDE_S8 = True
    GROWTH_BOOST = 0.0  # me-variation interpretation

    runs = []

    print("--- Phase 1: 4 models x 3 omega_m priors (with S8, with Cep prior 0.05) ---")
    for prior_name, (prior_label, sig_om) in OMH2_PRIORS.items():
        print(f"\n=== Prior {prior_name} ({prior_label}) ===")
        for model in ["M0", "M1", "M2", "Mco"]:
            tag = f"{model}_{prior_name}"
            r = run_one(tag, model, sigma_omh2=sig_om,
                        sig_dM=SIG_DM_DEFAULT, sig_frd=SIG_FRD_DEFAULT,
                        include_S8=INCLUDE_S8, growth_boost=GROWTH_BOOST)
            runs.append(r)
            print(f"  {tag:<12s}: logZ={r['log_Z']:8.2f}  H0={r['H0']:6.2f}±{r['H0_std']:.2f}"
                  f"  omh2={r['omh2']:.4f}±{r['omh2_std']:.4f}  ESS={r['n_eff']}")

    # Decomposition tests at INTERMEDIATE prior
    print("\n--- Phase 2: Decomposition (intermediate omega_m prior) ---")
    sig_om_decomp = 0.005
    print("\nTest A: M2 with Cepheid prior OFF")
    r = run_one("M2_noCepPrior", "M2", sigma_omh2=sig_om_decomp,
                sig_dM=0.0, sig_frd=SIG_FRD_DEFAULT, include_S8=INCLUDE_S8,
                growth_boost=0.0)
    runs.append(r)
    print(f"  logZ={r['log_Z']:.2f}  H0={r['H0']:.2f}  dMcep={r.get('dMcep', 0):+.4f}±{r.get('dMcep_std', 0):.4f}")

    print("\nTest B: M1 me-var with frd prior OFF")
    r = run_one("M1_noFrdPrior", "M1", sigma_omh2=sig_om_decomp,
                sig_dM=SIG_DM_DEFAULT, sig_frd=0.0, include_S8=INCLUDE_S8,
                growth_boost=0.0)
    runs.append(r)
    print(f"  logZ={r['log_Z']:.2f}  H0={r['H0']:.2f}  frd={r.get('frd', 0):.4f}±{r.get('frd_std', 0):.4f}")

    print("\nTest C: ALL priors OFF on extensions (data only)")
    r = run_one("Mco_dataOnly", "Mco", sigma_omh2=sig_om_decomp,
                sig_dM=0.0, sig_frd=0.0, include_S8=INCLUDE_S8,
                growth_boost=0.0)
    runs.append(r)
    print(f"  logZ={r['log_Z']:.2f}  H0={r['H0']:.2f}  dMcep={r.get('dMcep', 0):+.4f}  frd={r.get('frd', 0):.4f}")

    # Convergence stability: re-run M2_tight with seed=99
    print("\n--- Phase 3: Convergence stability (M2_tight repeat seed=99) ---")
    r = run_one("M2_tight_repeat", "M2", sigma_omh2=0.001,
                sig_dM=SIG_DM_DEFAULT, sig_frd=SIG_FRD_DEFAULT,
                include_S8=INCLUDE_S8, growth_boost=0.0, seed=99)
    runs.append(r)
    print(f"  logZ={r['log_Z']:.2f}  H0={r['H0']:.2f}  vs original "
          f"logZ={[x['log_Z'] for x in runs if x['name']=='M2_tight'][0]:.2f}")

    # Save all
    out_path = Path(__file__).resolve().parents[1] / "out_robustness.json"
    with open(out_path, "w") as f:
        json.dump(runs, f, indent=2)
    print(f"\nSaved {len(runs)} runs to {out_path}")

    # Summary table
    print("\n\n=========================================================")
    print("SUMMARY TABLE  (S8 included, growth_boost = 0 = me-variation)")
    print("=========================================================")
    z_ref = {p: next(r['log_Z'] for r in runs if r['name']==f"M0_{p}") for p in ["tight","interm","loose"]}
    for prior_name in ["tight", "interm", "loose"]:
        print(f"\n*** Prior {prior_name} ***  (M0 reference logZ = {z_ref[prior_name]:.2f})")
        print(f"{'Model':<10s} {'logZ':>8s} {'lnB':>8s} {'H0':>10s} {'omh2':>15s} {'frd or dMcep':>18s}")
        for model in ["M0", "M1", "M2", "Mco"]:
            r = next(rr for rr in runs if rr['name']==f"{model}_{prior_name}")
            extra = ""
            if 'frd' in r: extra += f"frd={r['frd']:.4f}"
            if 'dMcep' in r: extra += f"  dM={r['dMcep']:+.4f}" if extra else f"dM={r['dMcep']:+.4f}"
            print(f"{model:<10s} {r['log_Z']:>8.2f} {r['log_Z']-z_ref[prior_name]:>+8.2f} "
                  f"{r['H0']:>5.2f}±{r['H0_std']:.2f}  {r['omh2']:.5f}±{r['omh2_std']:.5f}  {extra}")


if __name__ == "__main__":
    main()
