"""Targeted me-var vs EDE comparison at intermediate omega_m prior.

To be run AFTER robustness_scan.py finishes. Adds the EDE-physics variant
(growth_boost = +6%) for M1 and Mco at sigma_omh2 = 0.005, with KiDS S8.

This isolates the rd-shift PHYSICS (me-var vs EDE) from everything else.
"""
import sys
from pathlib import Path
import numpy as np
import json

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts.robustness_scan import run_one


def main():
    runs = []
    print("=== me-var vs EDE comparison (sigma_omh2=0.005, S8 included) ===\n")

    for model in ["M1", "Mco"]:
        for label, boost in [("me", 0.0), ("EDE", 0.06)]:
            tag = f"{model}_{label}_interm"
            r = run_one(tag, model, sigma_omh2=0.005,
                        sig_dM=0.05, sig_frd=0.02,
                        include_S8=True, growth_boost=boost)
            runs.append(r)
            print(f"  {tag:<18s}: logZ={r['log_Z']:8.2f}  H0={r['H0']:.2f}±{r['H0_std']:.2f}  "
                  f"omh2={r['omh2']:.4f}")

    # Bayes factors me-var vs EDE
    print("\nKey Bayes factors (me-var vs EDE):")
    for model in ["M1", "Mco"]:
        z_me = next(r['log_Z'] for r in runs if r['name']==f"{model}_me_interm")
        z_ede = next(r['log_Z'] for r in runs if r['name']==f"{model}_EDE_interm")
        print(f"  ln B({model}_me / {model}_EDE) = {z_me - z_ede:+.3f}")

    out = Path(__file__).resolve().parents[1] / "out_me_vs_ede.json"
    with open(out, "w") as f:
        json.dump(runs, f, indent=2)
    print(f"\nSaved to {out}")


if __name__ == "__main__":
    main()
