"""Post-hoc analysis of out_robustness.json from robustness_scan.py.

Produces:
  - Summary table per omega_m prior
  - Bayes factors with appropriate references
  - Stratified conclusion: robust vs prior-dependent vs inconclusive
"""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main():
    runs = json.load(open(ROOT / "out_robustness.json"))
    by_name = {r["name"]: r for r in runs}

    print("=" * 78)
    print("ROBUSTNESS ANALYSIS — STRATIFIED CONCLUSION")
    print("=" * 78)
    print()

    # 1. Primary table: 4 models x 3 priors
    print("Primary scan (with KiDS S8 + me-variation growth_boost = 0):")
    print()
    for prior_name in ["tight", "interm", "loose"]:
        if f"M0_{prior_name}" not in by_name:
            continue
        z_M0 = by_name[f"M0_{prior_name}"]["log_Z"]
        print(f"--- omega_m prior {prior_name:<6s} (M0 logZ = {z_M0:.2f}) ---")
        print(f"  {'Model':<8s} {'logZ':>8s} {'lnB/M0':>8s} {'H0':>14s} {'omh2':>15s} {'frd/dM':>16s}")
        for model in ["M0", "M1", "M2", "Mco"]:
            r = by_name.get(f"{model}_{prior_name}")
            if r is None: continue
            extra = ""
            if 'frd' in r:
                extra += f"frd={r['frd']:.4f}±{r['frd_std']:.4f}"
            if 'dMcep' in r:
                e2 = f"dM={r['dMcep']:+.4f}±{r['dMcep_std']:.4f}"
                extra = e2 if not extra else f"{extra} {e2}"
            print(f"  {model:<8s} {r['log_Z']:>8.2f} {r['log_Z']-z_M0:>+8.2f} "
                  f"{r['H0']:>5.2f}±{r['H0_std']:.2f}  "
                  f"{r['omh2']:.5f}±{r['omh2_std']:.5f}  {extra}")
        print()

    # 2. Bayes factor evolution
    print("=" * 78)
    print("How Bayes factors EVOLVE with omega_m prior tightness:")
    print()
    print(f"  {'Comparison':<25s} {'tight':>10s} {'interm':>10s} {'loose':>10s}")
    for tag, m1, m2 in [
        ("M2 vs M0",   "M2", "M0"),
        ("M1 vs M0",   "M1", "M0"),
        ("Mco vs M0",  "Mco", "M0"),
        ("M2 vs M1",   "M2", "M1"),
        ("M2 vs Mco",  "M2", "Mco"),
    ]:
        row = []
        for p in ["tight", "interm", "loose"]:
            r1 = by_name.get(f"{m1}_{p}"); r2 = by_name.get(f"{m2}_{p}")
            if r1 and r2:
                row.append(f"{r1['log_Z']-r2['log_Z']:+.2f}")
            else:
                row.append("--")
        print(f"  {tag:<25s} {row[0]:>10s} {row[1]:>10s} {row[2]:>10s}")
    print()

    # 3. H0 marginal evolution
    print("=" * 78)
    print("H0 marginal under each prior:")
    print()
    print(f"  {'Model':<8s} {'tight':>14s} {'interm':>14s} {'loose':>14s}")
    for model in ["M0", "M1", "M2", "Mco"]:
        row = []
        for p in ["tight", "interm", "loose"]:
            r = by_name.get(f"{model}_{p}")
            if r:
                row.append(f"{r['H0']:>5.2f}±{r['H0_std']:.2f}")
            else:
                row.append("--")
        print(f"  {model:<8s}  {row[0]}    {row[1]}    {row[2]}")
    print()

    # 4. Decomposition tests
    print("=" * 78)
    print("Decomposition tests (intermediate omega_m prior):")
    for tag in ["M2_noCepPrior", "M1_noFrdPrior", "Mco_dataOnly"]:
        if tag not in by_name: continue
        r = by_name[tag]
        print(f"\n  {tag}:")
        print(f"    logZ = {r['log_Z']:.2f}  H0 = {r['H0']:.2f}±{r['H0_std']:.2f}  "
              f"omh2 = {r['omh2']:.4f}±{r['omh2_std']:.4f}")
        if 'dMcep' in r:
            print(f"    dMcep posterior = {r['dMcep']:+.4f}±{r['dMcep_std']:.4f}")
        if 'frd' in r:
            print(f"    frd posterior = {r['frd']:.4f}±{r['frd_std']:.4f}")
    print()

    # 5. Convergence
    print("=" * 78)
    if "M2_tight_repeat" in by_name and "M2_tight" in by_name:
        z1 = by_name["M2_tight"]["log_Z"]
        z2 = by_name["M2_tight_repeat"]["log_Z"]
        print(f"Convergence stability (M2_tight): logZ seed=42 = {z1:.3f}, seed=99 = {z2:.3f}, "
              f"|Δ| = {abs(z2-z1):.3f}")
        if abs(z2-z1) < 0.5:
            print("  → log Z stable to seed change (Δ < 0.5 nats)")
        else:
            print(f"  → log Z FLUCTUATING by {abs(z2-z1):.2f} nats — increase n_live or check sampler")
    print()

    # 6. Stratified conclusion
    print("=" * 78)
    print("STRATIFIED CONCLUSION:")
    print("=" * 78)

    # Robust: M2 vs M0
    bs = []
    for p in ["tight", "interm", "loose"]:
        r1 = by_name.get(f"M2_{p}"); r0 = by_name.get(f"M0_{p}")
        if r1 and r0:
            bs.append(r1["log_Z"] - r0["log_Z"])
    print(f"\n[ROBUST] M2 (Cepheid bias) vs M0 (LCDM):")
    print(f"  ln B = {bs[0]:.2f} (tight), {bs[1]:.2f} (interm), {bs[2]:.2f} (loose)")
    print(f"  → {'consistently strong' if min(bs) > 5 else 'prior-dependent'}")

    bs = []
    for p in ["tight", "interm", "loose"]:
        r1 = by_name.get(f"M1_{p}"); r0 = by_name.get(f"M0_{p}")
        if r1 and r0:
            bs.append(r1["log_Z"] - r0["log_Z"])
    print(f"\n[??] M1 (rd shift) vs M0 (LCDM):")
    print(f"  ln B = {bs[0]:.2f} (tight), {bs[1]:.2f} (interm), {bs[2]:.2f} (loose)")
    if min(bs) > 5: print("  → ROBUST: M1 always preferred")
    elif max(bs) > 5 and min(bs) < 1: print("  → PRIOR-DEPENDENT: outcome depends on omega_m prior")
    else: print("  → INCONCLUSIVE")

    bs = []
    for p in ["tight", "interm", "loose"]:
        r1 = by_name.get(f"M2_{p}"); rm1 = by_name.get(f"M1_{p}")
        if r1 and rm1:
            bs.append(r1["log_Z"] - rm1["log_Z"])
    print(f"\n[??] M2 vs M1:")
    print(f"  ln B = {bs[0]:.2f} (tight), {bs[1]:.2f} (interm), {bs[2]:.2f} (loose)")
    if min(bs) > 2.5: print("  → ROBUST: M2 always preferred")
    elif max(bs) > 2.5 and min(bs) < -2.5: print("  → SIGN-FLIPS: omega_m prior decides")
    elif abs(max(bs)) < 1 and abs(min(bs)) < 1: print("  → INCONCLUSIVE")
    else: print(f"  → PRIOR-DEPENDENT (not sign-flip)")


if __name__ == "__main__":
    main()
