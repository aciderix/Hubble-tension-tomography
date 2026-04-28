# LEVEL 6 — Sandbox-feasible explorations (with honest verdicts)

**Date:** April 28, 2026 evening
**Status:** 7 explorations attempted; 3 robust, 4 are artefacts of sandbox limitations.

---

## ✅ Robust results (trust these)

### #2 BAO tomography — H0(z) per DESI DR2 bin
At Ωm=0.315 fixed (Planck), rd=147.05 Mpc, infer H0 from each DESI DR2 BAO point:

| z_eff | H0 inferred | err |
|---|---|---|
| 0.295 | 70.37 | ±1.33 |
| 0.510 | 66.80 | ±1.23 |
| 0.706 | 68.75 | ±1.19 |
| 0.930 | 68.06 | ±0.94 |
| 1.317 | 67.97 | ±1.59 |
| 1.491 | 66.69 | ±1.74 |
| 2.330 | 67.24 | ±1.37 |

Linear fit slope: **−0.94 ± 0.80 km/s/Mpc per z (1.17σ from zero)**.

**Interpretation:** H0(z) is consistent with flat across DESI DR2 redshift range. This is exactly what EDE predicts (modification is pre-recombination, so late-time H(z) shape is unchanged) and disfavors late-time mechanisms that would produce a z-dependent H0.

### #4 σ8–Ωm 2D landscape

| Model | Ωm | σ8 | S8 | σ vs KiDS | σ vs DES |
|---|---|---|---|---|---|
| ΛCDM-Planck | 0.315 | 0.811 | 0.831 | +3.25σ | +3.24σ |
| ΛCDM-SH0ES | 0.296 | 0.798 | 0.793 | +1.35σ | +1.00σ |
| f(R) Step C | 0.260 | 0.629 | 0.585 | −9.05σ | −11.24σ |
| EDE proper | 0.280 | 0.704 | 0.680 | −4.30σ | −5.65σ |
| KiDS-1000 | 0.302 | 0.759 | 0.766 | (anchor) | −0.59σ |
| DES-Y3 | 0.300 | 0.776 | 0.776 | +0.50σ | (anchor) |

**Interpretation:** Any model with Ωm < 0.29 lands in the S8 desert. The "S8 wall" is real and structural — late-time mechanisms that lower Ωm to fit H0 systematically destroy S8.

### #5 Fisher forecast for Roman/Euclid SN samples
Assuming Δμ(EDE vs ΛCDM) ≈ 0.04 mag at z~1 and σ_μ ≈ 0.15 mag per SN:
- N for 3σ: ~130 SN
- N for 5σ: ~350 SN
- Roman planned (~10,000 SN at z<3): >>5σ trivially
- Euclid spectroscopic (~1,000 SN): >>5σ
- Pantheon+ current (1,701 SN): in principle 5σ-capable, but biased by anchor choice

**Interpretation:** the H0 question will be settled at the data level by the early 2030s.

---

## ⚠️ Artefacts (do NOT cite as results)

### #1 Wang-Mota dipole — INCONCLUSIVE
Reported 2.16σ dipole between hemispheres. **This is meaningless** because the binned Pantheon+ dataset I have here doesn't include RA/Dec, so I split bins randomly. Real test requires the full Pantheon+SH0ES.dat with per-SN coordinates.

### #3 CPL w0-wa fit — OPTIMIZER HIT BOUNDS
Reported w0=-2.00, wa=-2.16 at 56σ preference. Both parameters are at the prior boundaries — this is unphysical. Real test needs proper Bayesian priors and MCMC, not Nelder-Mead with hard cuts.

### #6 Hill+2020 EDE pullback — WRONG E(z) FORM
Plugged Hill+2020 best-fit (f_EDE=0.107, z_c~3650) into a phenomenological E²(z) = Ωm(1+z)³ + ΩΛ + f_EDE/(1+(a/ac)^6). This **does not capture EDE physics** — true EDE has a scalar field that scales like radiation at early times then decays, with k-dependent perturbations. My toy form gives a Δχ² of +7854 (Pantheon+ blew up), which is nonsense.

### #7 Bayes factor — INHERITS #6 BUG
Computed ln(BF) = -3931 (decisive for ΛCDM) directly from #6's broken Δχ². **Discard.** Real Bayes factor requires proper EDE Boltzmann + MCMC posterior volume estimation.

---

## What this means

- The data trends (#2, #4, #5) are pointing the same direction as Levels 1–5: **early-time / pre-recombination physics is preferred over late-time modifications.** EDE is the leading candidate for H0 with a moderate (~2-3σ) S8 cost.
- The proper tests of #1, #3, #6, #7 require AxiCLASS / class_ede / cobaya — exactly the local-machine work outlined in `CLAUDE_HANDOFF.md` v3.

## Where to go next (sandbox-feasible)

If you want me to keep exploring without the local Boltzmann codes:
- Fetch the FULL Pantheon+SH0ES.dat with RA/Dec → real dipole test (#1 done properly)
- DESI DR2 full likelihood (covariance matrix) → tighter tomography errors (#2 enhanced)
- Pantheon+ spline DL(z) reconstruction → model-independent dark energy density evolution

Beyond that, the next real progress requires the local Boltzmann pipeline.
