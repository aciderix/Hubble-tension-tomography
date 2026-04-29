"""S8 dataset sensitivity scan.

Tests how the M2/M1 ranking depends on which S8 dataset is used:

  KiDS-1000:           S8 = 0.766 ± 0.020   (Asgari+ 2021)
  DES Y3 cosmic shear: S8 = 0.776 ± 0.017   (Amon+ 2022 / Secco+ 2022)
  KiDS+DES combined:   S8 = 0.770 ± 0.013   (Amon, Heymans 2022)
  Planck lensing:      S8 = 0.832 ± 0.013   (Planck 2018 collab + Carron+ 2022)
                                             *consistent with Planck-LCDM*
  No S8 (skip)

We re-run only the 4 most informative models (M0, M2, M1_me, M1_EDE)
under intermediate omega_m prior, varying the S8 likelihood.

The contrarian thesis predicts:
  - KiDS / DES / combined → M2 wins, M1_EDE worst
  - Planck lensing       → M0 LCDM wins, M2 disfavored, M1_EDE viable
  - No S8                → M2 still wins (from anchor-free arguments) but by less
"""
import sys
from pathlib import Path
import numpy as np
import json

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts.robustness_scan import run_one
from src.growth import sigma8 as sigma8_calc

# Patch run_one to allow alternative S8 dataset by monkey-patching the
# chi2_S8 function. Cleaner: override S8_KIDS, S8_KIDS_SIG in the module.
import scripts.robustness_scan as RS


def main():
    runs = []
    print("=== S8 dataset sensitivity scan ===\n")

    S8_DATASETS = {
        "KIDS":     ("KiDS-1000",          0.766, 0.020),
        "DES":      ("DES Y3 cosmic shear", 0.776, 0.017),
        "KIDSDES":  ("KiDS+DES combined",   0.770, 0.013),
        "PLENS":    ("Planck CMB lensing",  0.832, 0.013),
        "NONE":     ("No S8 likelihood",    None, None),
    }

    for ds_key, (label, s8_val, s8_sig) in S8_DATASETS.items():
        # Override module constants
        if s8_val is not None:
            RS.S8_KIDS = s8_val
            RS.S8_KIDS_SIG = s8_sig
            include_S8 = True
        else:
            include_S8 = False

        print(f"\n=== S8 source: {label} (S8={s8_val} ± {s8_sig}) ===")
        for model in ["M0", "M1", "M2"]:
            for boost_label, boost in ([("me", 0.0), ("EDE", 0.06)]
                                        if model == "M1"
                                        else [("", 0.0)]):
                tag = f"{model}{'_'+boost_label if boost_label else ''}_{ds_key}"
                r = run_one(tag, model, sigma_omh2=0.005,
                            sig_dM=0.05, sig_frd=0.02,
                            include_S8=include_S8, growth_boost=boost)
                runs.append(r)
                print(f"  {tag:<22s}: logZ={r['log_Z']:8.2f}  H0={r['H0']:.2f}±{r['H0_std']:.2f}")

    out = Path(__file__).resolve().parents[1] / "out_s8_sensitivity.json"
    with open(out, "w") as f:
        json.dump(runs, f, indent=2)

    # Summary table
    print("\n\n=== S8 sensitivity SUMMARY: ln B vs M0_<dataset> ===")
    print(f"  {'Model':<8s}  {'KiDS':>10s}  {'DES':>10s}  {'KiDS+DES':>10s}  {'Planck-lens':>11s}  {'No-S8':>10s}")
    for model_disp, model_tag in [("M2", "M2"), ("M1_me", "M1_me"), ("M1_EDE", "M1_EDE")]:
        row = []
        for ds in ["KIDS", "DES", "KIDSDES", "PLENS", "NONE"]:
            r0 = next((rr for rr in runs if rr['name']==f"M0_{ds}"), None)
            rm = next((rr for rr in runs if rr['name']==f"{model_tag}_{ds}"), None)
            if r0 and rm:
                row.append(f"{rm['log_Z']-r0['log_Z']:+.2f}")
            else:
                row.append("--")
        print(f"  {model_disp:<8s}  {row[0]:>10s}  {row[1]:>10s}  {row[2]:>10s}  {row[3]:>11s}  {row[4]:>10s}")

    print(f"\nSaved to {out}")


if __name__ == "__main__":
    main()
