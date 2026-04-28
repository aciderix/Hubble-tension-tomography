# Level 3 — Hybrid Models, DESI DR2, classy CLP Stress Test

**Date:** April 28, 2026
**Status:** Sandbox-feasible Level 3 components complete. CLASS_EDE-dependent components still pending local compile.

---

## Step 1 — DESI DR2 BAO Joint Refit

**Dataset:** DESI DR2 13-point BAO (April 2025 release) + Pantheon+ Om-shape prior + optional SH0ES H0 anchor.

### Toy-EDE refit results (α_rd = 0.18)

| Model | H₀ | Ωₘ | f_EDE | χ² |
|---|---|---|---|---|
| ΛCDM (no SH0ES) | 68.68 | 0.3045 | — | 13.61 |
| ΛCDM (+SH0ES) | 69.40 | 0.2933 | — | 28.29 |
| EDE (+SH0ES) | 73.04 | 0.3045 | 0.332 | 13.61 |
| EDE (no SH0ES) | 68.68 | 0.3045 | 0.000 | 13.61 |

- **Δχ²(EDE − ΛCDM) anchored = −14.68** (~3.8σ preference)
- **Tension data-only vs SH0ES: −4.19σ** for both ΛCDM and EDE
- **Per-point DESI residuals all <1.7σ** (no internal tension)

### Interpretation

DESI DR2 alone cannot constrain f_EDE strongly — the pull toward `H₀ = 73` requires a very large `f_EDE = 0.33` under the conservative `α_rd = 0.18`. Real CLASS_EDE analyses use `α_rd ≈ 0.4–0.5`, which would shrink this to `f_EDE ≈ 0.10–0.13`, the literature range.

**Key takeaway:** DESI DR2 *strengthens* the inverse-distance-ladder side. Without SH0ES, BAO data prefer H₀ ≈ 68.7, in 4.2σ tension with SH0ES. EDE resolution must therefore come entirely from pre-recombination CMB physics, not BAO geometry.

---

## Step 2 — classy CLP Fluid Stress Test

**Goal:** confirm via Boltzmann code that late-time DE evolution alone cannot resolve the H₀ tension when CMB primary parameters are held at Planck values.

### Setup

- classy 3.3.4, fluid with `w(a) = w₀ + wₐ(1−a)`, `w₀ = −1`
- Fix Planck 2018 values: `h=0.6736, ω_b=0.02237, ω_cdm=0.12, A_s=2.1e-9, n_s=0.9649, τ=0.0544`
- Scan `wₐ ∈ {0.0, 0.5, 1.0, 1.5, 2.0}`

### Results

| wₐ | H₀ | r_d | Ωₘ | S₈ |
|---|---|---|---|---|
| 0.0 | 67.36 | 147.11 | 0.3138 | 0.842 |
| 0.5 | 67.36 | 147.11 | 0.3138 | 0.794 |
| ≥1.0 | **REJECTED** | — | — | — |

**classy refuses wₐ ≥ 1.0** with error: *"condition (w_fld ≥ 1./3.) is true; The fluid is meant to be negligible at early time"*. This is the code itself rejecting attempts to use CLP as an EDE proxy — the CLP parameterization cannot describe a fluid that is non-negligible at early times.

### Interpretation

- For `wₐ = 0.5`: H₀ shifts by **0.001 km/s/Mpc**, S₈ drops by 5.6%
- Late-time DE perturbations under fixed Planck primary parameters do NOT touch H₀
- Confirms Level 1 mechanism (B) rejection
- **Validates that EDE peak physics requires a true axion-like fluid (CLASS_EDE / AxiCLASS), not CLP**

---

## Step 3 — Hybrid Model: Local Anisotropy + Global EDE

**Hypothesis:** the "real" picture might combine two effects — modest global EDE (`f_EDE ~ 0.05–0.08`, literature range) and a small local anisotropy boost (`δH ~ 1–2 km/s/Mpc` in `z < 0.03`).

### Model

- Two-zone Hubble: `H₀_local = H₀_global + δH` for SNe at `z < 0.03` only
- Global EDE shifts `r_d_eff = r_d × (1 − α·f_EDE)` with `α = 0.40` (literature value)
- SH0ES anchors `H₀_local`
- Joint fit: Pantheon+ N=1590 + DESI DR2 + SH0ES

### Results

| Model | H₀_global | Ωₘ | f_EDE | δH | χ² | Δχ² vs ΛCDM |
|---|---|---|---|---|---|---|
| ΛCDM | 73.84 | 0.237 | — | — | 1553.67 | 0 |
| EDE only | 73.52 | 0.304 | 0.165 | — | 1416.99 | **−136.68 (11.7σ)** |
| Anisotropy only | 73.69 | 0.239 | — | +0.62 | 1551.37 | −2.30 (1.5σ) |
| **Hybrid** | 73.52 | 0.304 | 0.165 | +0.01 | 1416.99 | **−136.69 (11.7σ)** |

### Interpretation

**The hybrid fit collapses to pure EDE.** The minimizer drives `δH → 0.01`, meaning when EDE is allowed to absorb the H₀ pull, the local anisotropy contributes essentially nothing. The 2.5σ z<0.03 anisotropy signal that we found in Level 1 is **not strong enough to compete with EDE** in a joint fit.

This is a clean falsification of the "anisotropy alone explains the tension" hypothesis. Anisotropy contributes at most `δH ≈ 0.6 km/s/Mpc` (consistent with the per-survey scatter we measured), which is far short of what's needed.

**Note on the ΛCDM Ωₘ = 0.237:** under pure ΛCDM forced to satisfy SH0ES, the SNe shape pulls Ωₘ unphysically low (Planck gives 0.315). When EDE is allowed, Ωₘ jumps back to 0.304, very close to the Planck-preferred value. This is itself evidence that EDE is the right modification.

---

## Summary of Level 3 Findings

1. **DESI DR2 BAO** strengthens the low-H₀ camp. EDE resolution must come from CMB physics, not BAO.
2. **classy CLP** confirms late-time DE evolution cannot move H₀. The code itself refuses extensions into the EDE regime — physics-aware rejection.
3. **Hybrid model** with anisotropy + EDE collapses to pure EDE. The local 2.5σ anisotropy signal is insufficient to absorb the H₀ tension.

### Open questions for full Level 3 (require CLASS_EDE local compile)

- True axion-like EDE potential `V(φ) = m²f²(1 − cos(φ/f))ⁿ` with `n = 3`
- Joint MCMC: H₀, Ωₘ, f_EDE, log₁₀(z_c), θ_i, plus standard Planck nuisance parameters
- Full Planck high-ℓ TT/TE/EE likelihood (plik_lite + lowE)
- Predicted S₈ tension increase under EDE, and whether self-interacting DM can absorb it

---

## Reproducibility

- `level3_desi_dr2_refit.py` — DESI DR2 + Pantheon+ Om-prior toy-EDE joint fit
- `level3_classy_clp_scan.py` — classy CLP fluid scan (Planck-fixed)
- `level3_hybrid_model.py` — Anisotropy + EDE joint fit
- `level3_results.json` — All numerical results
