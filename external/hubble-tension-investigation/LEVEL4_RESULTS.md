# Level 4 Results — Robustness, Anisotropy, H(z) Reconstruction, MG vs EDE, S₈ Cross-check

**Date:** April 28, 2026
**Repo:** [aciderix/hubble-tension-investigation](https://github.com/aciderix/hubble-tension-investigation)
**Status:** Level 4 complete. Level 5 (MGCAMB / growth-equation validation) queued.

---

## TL;DR

Four independent diagnostics (A, B, C, D) were executed against Pantheon+ (1380 SNe), DESI DR2 BAO (6 points), Cosmic Chronometers (31 H(z) measurements), and the SH0ES H₀ anchor.

**Result:** Standard ΛCDM with the SH0ES anchor is rejected at Δχ² ≈ +24 vs the best alternative. When we add the S₈ cross-check (KiDS-1000, DES Y3), the picture is **inverted from the EDE-favored consensus**:

- **f(R) late-time modification** is the best simultaneous fit to H₀ and S₈ (Δχ²_total = −18 vs ΛCDM-Planck).
- **EDE r_d-shrink** resolves H₀ (Δχ² = −10.8) but worsens S₈ tension to **+6σ** vs KiDS/DES.
- **DGP brane** achieves H₀=73 but requires Ωₘ=0.21 (excluded by Planck CMB) and crashes S₈ at −9.9σ.

---

## Step A — Hemispheric Anisotropy & Dipole

**Method:** Pantheon+ z<0.5 split into 100 random hemispheres + best dipole fit.

**Result:**
- N(z<0.5) = 1380 SNe
- Max hemisphere asymmetry: ΔH₀ = 0.58 (1.24σ; look-elsewhere corrected ~0σ)
- Best dipole: H₀_iso = 73.71±1.36, amplitude = 104±234 (0.45σ), direction RA=269.5°, DEC=−0.4°
- Misalignment from CMB dipole: 79°

**Conclusion:** Pantheon+ SH0ES anchor is isotropic. The H₀ tension is **cosmological, not local/directional**.

---

## Step B — Non-parametric H(z) Reconstruction

**Method:** Piecewise-linear H(z) on 8 z-bins. Three fits.

| Fit | χ² | r_d (Mpc) | H(z=0) |
|---|---|---|---|
| 1 (SH0ES + r_d=147) | 724.08 | 147.05 (fixed) | 69.43 |
| 2 (no anchor) | 704.67 | 147.05 (fixed) | 67.24 |
| 3 (SH0ES + r_d free) | **712.77** | **🎯 139.78** | 71.52 |

When SH0ES is forced AND r_d is free, the data prefer **r_d = 139.78 Mpc** — a 4.94% shrinkage from Planck. Tension localized at **z ≈ 0.3–1.0**.

---

## Step C — Modified Gravity vs EDE (Background H(z) only)

| Model | k | χ² | H₀ |
|---|---|---|---|
| ΛCDM-Planck (no anchor) | 2 | 715.87 | 68.52 |
| ΛCDM (SH0ES, r_d=147) | 2 | **731.71** | 69.25 |
| DGP brane | 3 | 716.42 | 73.04 |
| **f(R) late-time** | 4 | **707.76** | **73.01** |
| wCDM | 3 | 729.88 | 68.83 |
| EDE proxy (r_d free) | 3 | 720.95 | 71.77 |

**Caveat:** f(R) here is `H_LCDM(z) · [1 − δ·exp(−z/z*)]` — phenomenological background-only. Level 5 needs MGCAMB.

---

## Step D — S₈ Cross-Check

| Solution | S₈ predicted | σ vs KiDS/DES avg |
|---|---|---|
| ΛCDM-Planck | 0.831 | +3.23σ |
| ΛCDM-SH0ES | 0.802 | +1.82σ |
| DGP brane | 0.567 | **−9.95σ** (dead) |
| **f(R) late-time** | **0.780** | **+0.71σ** |
| EDE r_d-shrink | 0.888 | **+6.11σ** |

**Combined ranking (χ²_H₀ + χ²_S₈):**

| Rank | Solution | Total χ² |
|---|---|---|
| ★ 1 | **f(R) late-time** | **708.26** |
| 2 | ΛCDM-Planck | 726.28 |
| 3 | ΛCDM-SH0ES | 735.01 |
| 4 | EDE r_d-shrink | 758.26 |
| 5 | DGP brane | 815.38 |

⚠️ σ₈ for f(R) used a phenomenological 5% boost factor. Must be replaced with MGCAMB.

---

## Roadmap — Level 5

1. **MGCAMB f(R) full run.** Hu-Sawicki, n=1, vary |f_R0|. Proper σ₈, μ(k,z), Σ(k,z).
2. **RSD test.** f·σ₈(z) vs eBOSS/DESI.
3. **CMB ISW.** f(R) vs Planck low-ℓ TT × lensing.
4. **Full Pantheon+ covariance** (1701×1701).
5. **Joint MCMC (cobaya).** f(R) + Pantheon+ + DESI + Planck + KiDS.
6. **If it survives:** preprint draft.

---

## Honest disclaimer

Exploratory work in a Python sandbox by an AI assistant + the researcher. Treat the f(R)-wins-over-EDE conclusion as a **hypothesis worth testing rigorously**, not as a proof. Publishable claim requires independent re-implementation in CLASS/MGCAMB, full systematic budget, and peer review.
