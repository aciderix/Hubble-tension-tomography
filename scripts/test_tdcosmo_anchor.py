"""Quick test: add TDCOSMO_NFW (Wong+ 2020) H0 = 73.3 ± 1.8 as Gaussian
prior to our two best models, see how Bayes factors shift.

Analytical prediction:
  Delta chi2 on M2 (H0=68.81 +/- 0.28) = (68.81 - 73.3)^2 / (0.28^2 + 1.8^2) = 6.07
  Delta chi2 on M1_me (H0=70.0 +/- 0.60) = (70 - 73.3)^2 / (0.6^2 + 1.8^2) = 3.03

Expected logZ shifts:
  M2 logZ -> -798.43 - (6.07 - 3.03)/2 + correction = approx -801.0
  M1_me logZ -> -806.78 stays + small correction from posterior shift

Expected ln B(M2/M1_me) shift: from +8.35 to approximately +5

If TDCOSMO_NFW IS RIGHT (mass-sheet not relaxed), our M2 conclusion would
weaken from "decisive" to "moderate".
If TDCOSMO_relaxed (Birrer+ 2020) is right (H0 = 67.4 +/- 4), conclusion
unchanged.
"""
import sys
from pathlib import Path
import numpy as np
from nautilus import Sampler

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts.robustness_scan import run_one, make_loglike, make_prior
import scripts.robustness_scan as RS
from src.external_anchors import chi2_anchor


def make_loglike_with_anchor(model, anchor_name, **kwargs):
    base = make_loglike(model, **kwargs)
    def ll(td):
        return base(td) - 0.5 * chi2_anchor(td["H0"], anchor_name)
    return ll


def run_with_anchor(name, model, anchor, **kwargs):
    n_live = kwargs.pop("n_live", 400)
    seed = kwargs.pop("seed", 42)
    ll = make_loglike_with_anchor(model, anchor, **kwargs)
    prior = make_prior(model)
    s = Sampler(prior, ll, n_live=n_live, pool=4, n_networks=4, seed=seed)
    s.run(verbose=False, discard_exploration=True)
    pts, log_w, log_l = s.posterior()
    keys = list(prior.keys)
    log_Z = s.log_z
    weights = np.exp(log_w - log_w.max()); weights /= weights.sum()
    H0 = pts[:, keys.index("H0")]
    h0m = (H0 * weights).sum()
    h0s = np.sqrt(((H0 - h0m) ** 2 * weights).sum())
    return {"name": name, "log_Z": float(log_Z),
            "H0": float(h0m), "H0_std": float(h0s)}


def main():
    print("=== TDCOSMO_NFW (Wong+ 2020 H0=73.3+/-1.8) test on M2, M1_me ===")
    print("Baseline (no TDCOSMO):")
    print("  M2_KIDS    : logZ = -798.43, H0 = 68.81 +/- 0.28")
    print("  M1_me_KIDS : logZ = -806.78, H0 = 70.00 +/- 0.60")
    print("  ln B(M2/M1_me) = +8.35")
    print()

    common = dict(sigma_omh2=0.005, sig_dM=0.05, sig_frd=0.02,
                  include_S8=True, growth_boost=0.0)

    print("Adding TDCOSMO_NFW Gaussian prior on H0...")
    rM2_W = run_with_anchor("M2_TDCOSMO_NFW", "M2", "TDCOSMO_W", **common)
    print(f"  M2 + TDCOSMO_W : logZ = {rM2_W['log_Z']:.2f}, H0 = {rM2_W['H0']:.2f} +/- {rM2_W['H0_std']:.2f}")
    rM1_W = run_with_anchor("M1_me_TDCOSMO_NFW", "M1", "TDCOSMO_W", **common)
    print(f"  M1_me + TDCOSMO_W : logZ = {rM1_W['log_Z']:.2f}, H0 = {rM1_W['H0']:.2f} +/- {rM1_W['H0_std']:.2f}")
    print()
    print(f"  ln B(M2/M1_me) with TDCOSMO_NFW = {rM2_W['log_Z'] - rM1_W['log_Z']:+.2f}")
    print()

    print("Adding TDCOSMO_relaxed Gaussian prior on H0...")
    rM2_B = run_with_anchor("M2_TDCOSMO_B", "M2", "TDCOSMO_B", **common)
    print(f"  M2 + TDCOSMO_B : logZ = {rM2_B['log_Z']:.2f}, H0 = {rM2_B['H0']:.2f} +/- {rM2_B['H0_std']:.2f}")
    rM1_B = run_with_anchor("M1_me_TDCOSMO_B", "M1", "TDCOSMO_B", **common)
    print(f"  M1_me + TDCOSMO_B : logZ = {rM1_B['log_Z']:.2f}, H0 = {rM1_B['H0']:.2f} +/- {rM1_B['H0_std']:.2f}")
    print()
    print(f"  ln B(M2/M1_me) with TDCOSMO_relaxed = {rM2_B['log_Z'] - rM1_B['log_Z']:+.2f}")


if __name__ == "__main__":
    main()
