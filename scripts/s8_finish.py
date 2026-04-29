"""Finish the S8 sensitivity scan: remaining 6 runs.

Already done:  KIDS x4, DES x4, KIDSDES x4, PLENS x (M0, M1_me) = 14
Missing: PLENS x (M1_EDE, M2), NONE x4
"""
import sys
from pathlib import Path
import json

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts.robustness_scan import run_one
import scripts.robustness_scan as RS


def main():
    runs = []

    # Planck lensing - finish
    RS.S8_KIDS = 0.832
    RS.S8_KIDS_SIG = 0.013
    print(f"\n=== S8 source: Planck CMB lensing (S8=0.832 ± 0.013), finishing ===")
    for tag, model, boost in [("M1_EDE_PLENS", "M1", 0.06),
                                ("M2_PLENS", "M2", 0.0)]:
        r = run_one(tag, model, sigma_omh2=0.005,
                    sig_dM=0.05, sig_frd=0.02,
                    include_S8=True, growth_boost=boost)
        runs.append(r)
        print(f"  {tag:<22s}: logZ={r['log_Z']:8.2f}  H0={r['H0']:.2f}±{r['H0_std']:.2f}")

    # No S8 (skip S8 likelihood)
    print(f"\n=== S8 source: No S8 likelihood ===")
    for tag, model, boost in [("M0_NONE", "M0", 0.0),
                                ("M1_me_NONE", "M1", 0.0),
                                ("M1_EDE_NONE", "M1", 0.06),
                                ("M2_NONE", "M2", 0.0)]:
        r = run_one(tag, model, sigma_omh2=0.005,
                    sig_dM=0.05, sig_frd=0.02,
                    include_S8=False, growth_boost=boost)
        runs.append(r)
        print(f"  {tag:<22s}: logZ={r['log_Z']:8.2f}  H0={r['H0']:.2f}±{r['H0_std']:.2f}")

    out = Path(__file__).resolve().parents[1] / "out_s8_finish.json"
    with open(out, "w") as f:
        json.dump(runs, f, indent=2)
    print(f"\nSaved to {out}")


if __name__ == "__main__":
    main()
