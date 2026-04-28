# ⚠️ ERRATUM — Level 4 Step D σ8 calculation was wrong

**Date:** April 28, 2026 (same day as Level 4 commit)
**Status:** Critical correction. Read this BEFORE citing LEVEL4_RESULTS.md.

---

## The bug

In Level 4 Step D, I computed σ8 for each H0 solution using a **naive scaling**:

```
σ8(model) ≈ σ8(Planck) · √(Ωm/0.315)   ← WRONG
```

This is **not how σ8 scales with Ωm** in LCDM. The actual dependence comes from the linear growth equation, which depends nonlinearly on Ωm via the matter-Λ transition redshift.

## The correct calculation

Level 5 Step E computed σ8 properly:
1. Linear growth ODE: D'' + (2 + H'/H)·D' − (3/2)·Ωm(a)·μ(k,a)·D = 0
2. Eisenstein-Hu transfer function T(k)
3. Top-hat window W(kR) at R = 8 Mpc/h
4. Calibrated against Planck (Ωm=0.315, h=0.6736 → σ8=0.811)

**Result:** σ8 at Ωm=0.260 (the f(R) Step C best-fit Ωm) is **0.629**, not 0.798.

## Impact on the rankings

| Model | Ωm | S8 (Step D, claimed) | S8 (Step E, real) | Verdict change |
|---|---|---|---|---|
| LCDM-Planck | 0.315 | 0.831 | 0.811 (calibrated) | small |
| LCDM-SH0ES | 0.296 | 0.802 | ~0.78 | minor |
| **f(R) Ωm-free** | **0.260** | **0.780 ✓** | **0.585 ✗** | **INVERTED** |
| EDE r_d-shrink | 0.304 | 0.888 | ~0.70 (with proper σ8) | also overestimated |

**The "f(R) wins" conclusion of Step D was an artifact of the wrong σ8 scaling.**

## What Step F+G found

- **f(R) with Ωm=0.315 FIXED:** can match S8 (≈0.76, −0.2σ KiDS) but H0 collapses to 60. Cannot resolve H0 tension.
- **f(R) with Ωm=0.260 (Step C fit):** good H0 (73.0) but S8≈0.59 (−9σ KiDS). Dead.

## The real conclusion

When the σ8 calculation is done properly:

1. **No simple late-time modification (f(R), DGP, wCDM) simultaneously satisfies H0=73 AND S8≈0.76 with Ωm≈0.31.**
2. The reason late-time fits achieve H0=73 in our scaffold is by lowering Ωm to ~0.26, which then crashes S8.
3. **Pre-recombination physics (EDE) remains the leading candidate** for the H0 problem because it changes r_d without requiring low Ωm.
4. EDE itself does enhance σ8 by ~5–7% (Hill+2020), worsening S8 modestly — but not by 6σ as Step D claimed.

## The thesis update

Lead with: **"Levels 1-4 confirmed EDE-like behavior (r_d shrinks to 140 Mpc agnostically, Step B). Late-time MG ansatzes can mimic this in background but fail S8 once growth is computed properly. The simultaneous H0+S8 problem requires either: (a) full MGCAMB f(R) computation we couldn't run in the sandbox, (b) interacting DE / decaying DM models, or (c) acceptance that EDE solves H0 at the cost of a moderate (~2σ) S8 pull."**

This is the honest framing. The "f(R) wins" headline was wrong.

## Lessons learned (for future runs)

- **Always compute σ8 from first principles** (growth ODE + transfer + top-hat). Never use scaling shortcuts.
- **Cross-check at every level** with an independent code (CAMB or CLASS) before claiming a result.
- **Sandbox limitations are real:** without MGCAMB, full f(R) growth with k-dependence cannot be properly tested. Move this work to a local machine with gcc/make.

## Files affected

- `LEVEL4_RESULTS.md` — Step D table is wrong, this erratum supersedes it
- `level4_data.json` — `D_S8_crosscheck.predictions` and `combined_ranking` are wrong
- `CLAUDE_HANDOFF.md` — "f(R) wins" line should be replaced with this erratum's conclusion
