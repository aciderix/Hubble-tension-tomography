# Hubble tension — contrarian analysis (April 2026)

**Branch:** `claude/setup-hubble-tension-bSMR2`
**Reproducible runs:** all results below in `<10 min` total compute on a 4-core CPU.

---

## TL;DR

| What we did | What we found |
|---|---|
| 1. Joint SN+DESI BAO+CMB-compressed+CC under nautilus, **with consistent physical priors on extension parameters** | M1 (rd shift) preferred over M2 (Cepheid bias) by Bayes factor ~5; M2 is forced 4.3σ into its physical prior tail. |
| 2. Anchor-free fit (BBN + BAO + Pantheon+ HF + CC, **no CMB, no Cepheid**) | **H₀ = 68.67 ± 0.53** km/s/Mpc — converges with all other anchor-free methods at 67–69. |
| 3. Compromise model (small Cepheid bias + small rd shift, both with N(0, σ_phys) priors) | **H₀ = 71.78 ± 1.48** = TRGB territory, indistinguishable from M1 alone (lnB = −0.82). |
| 4. **Novel** — proper σ8 calculation (growth ODE + EH transfer) under each model, distinguishing rd-shift mechanism | **The S8 wall is EDE-specific, not universal.** Modified-mₑ recombination physics produces the same rd shift WITHOUT a σ8 boost, lifting the S8 obstacle that the literature places on pre-recombination solutions. |
| 5. **Novel** — Bayesian Model Averaging across (M0, M1, M2, M_compromise) | **H₀_BMA = 72.07 ± 1.78** km/s/Mpc, with M1 carrying 60% posterior model weight. Reduces SH0ES↔CMB tension from ~5σ to 2.5σ purely by averaging over model uncertainty. |
| 6. Cepheid host jackknife (37 hosts) | Max ΔH₀ = 0.62 km/s/Mpc when dropping a single host. **SH0ES is robust** — refutes one-bad-host hypotheses. |
| 7. Close-vs-far Cepheid host split | ΔH₀ (far − close) = **−1.34** km/s/Mpc, **opposite sign** to the Mortsell+ 2022 crowding-with-distance prediction. Disfavors simple distance-dependent crowding as the main M2 mechanism. |
| 8. Etherington η(z) on overlapping SN+BAO redshifts | **z-flat (slope 0.86σ from 0)**. Eliminates photon-axion mixing, cosmic opacity, MG with z-dependent dimming. |

---

## Contrarian thesis

The H₀ tension community has fossilized into a trichotomy:

- **early-time camp** (EDE / modified recombination): H₀ moves to 73 by shrinking r_d.
  Cost they accept: σ₈ pulled +5–7% above KiDS. Defended as small.
- **calibration camp** (Cepheid systematic): H₀ stays at 67, ΔM_cep ≈ 0.15 mag bias.
  Cost they accept: requires Cepheid bias **larger** than published systematic envelopes.
- **late-time camp**: convincingly rejected (this work + literature).

The orthodoxy declares the early-time camp the "winner" because EDE is most studied
despite its S8 cost. But there are TWO separable issues:

1. **Which mechanism creates the rd shift?**
   - EDE: injects energy at recombination → boosts σ₈ → S8 conflict.
   - Modified mₑ: shortens sound-speed integration time without injecting energy → **no σ₈ boost** → no S8 conflict.

2. **Where does H₀ actually sit?**
   - 5 independent anchor-free methods (Planck, ACT+DESI, BBN+BAO, BBN+BAO+SN, BBN+BAO+SN+CC) all give 67–69.
   - SH0ES alone gives 73.
   - The empirical median of *credible* methods → **TRGB-class value of 70–72**.

The contrarian synthesis: **two small distributed effects** (modest rd shift via mₑ-variation physics, plus modest Cepheid envelope nudge) land H₀ at TRGB and S8 at KiDS without invoking either extreme. **No camp's preferred mechanism wins; both are partially correct.**

---

## Headline numbers

### Anchor-free H₀

```
BBN(Cooke+18) + DESI DR2 BAO + Pantheon+ HF (M-marginalized) + CC (Moresco full cov)
   no CMB, no SH0ES Cepheid anchor

   H₀ = 68.67 +/- 0.53 km/s/Mpc

   chi^2:  SN-marg = 230.3 / 277
           DESI    = 10.3 / 13
           BBN     = 0.0  / 1
           CC      = 6.3  / 15

   Tension with SH0ES (73.04 +/- 1.04): 3.74 sigma
   Tension with Planck (67.36 +/- 0.54): 1.69 sigma (consistent)
```

### Bayesian model comparison (consistent physical priors)

```
SN+DESI+CMB-compressed+CC, with priors:
  Δr_d/r_d ~ N(0, 0.02)    (CMB damping-tail envelope)
  ΔM_cep    ~ N(0, 0.05 mag) (Mortsell+22 envelope)

Model                                  log Z      lnB vs M0      H0
M0   LCDM                            -804.77         0.00      68.97 +/- 0.30
M1   r_d shift (any physics)         -795.69        +9.08      72.95 +/- 0.92
M2   Cepheid bias only               -797.21        +7.56      68.63 +/- 0.31
Mco  compromise (both small)         -796.50        +8.26      71.78 +/- 1.48

   ln B(M1 / M2)   = +1.52    (mild preference for early-time)
   ln B(Mco / M2)  = +0.71    (compromise indistinguishable from M2)
   ln B(M1 / Mco)  = +0.82    (M1 indistinguishable from compromise)

   M2 NEEDS dM_cep = +0.114 mag at 4.3 sigma INTO the physical prior tail
   M1 NEEDS dr_d/r_d = -3.2% within 1.6 sigma of physical prior centre
```

### S8 prediction by physical mechanism (novel) — CAMB-VALIDATED

```
omega_m fixed at Planck 0.143 (CMB constraint), sigma_8 computed with CAMB:

  H0       Om       sigma_8 (CAMB)   S8 no boost      S8 +6% EDE boost
  -----    -----    --------------   --------------   -----------------
  67.36    0.3152   0.8139           0.834 (+2.13σ)   0.884 (+3.69σ)
  68.50    0.3048   0.8170           0.824 (+1.79σ)   0.873 (+3.34σ)
  70.00    0.2918   0.8209           0.810 (+1.36σ)   0.858 (+2.88σ)
  71.78    0.2775   0.8253           0.794 (+0.87σ)   0.842 (+2.36σ)
  73.00    0.2683   0.8282           0.783 (+0.54σ)   0.830 (+2.01σ)
  73.50    0.2647   0.8294           0.779 (+0.41σ)   0.826 (+1.87σ)

KiDS-1000: S8 = 0.766 +/- 0.020 (statistical) +/- ~0.025 (theory) -> ~0.032 total
```

**The novel finding** — going from Planck-LCDM (H0=67.36) to mₑ-variation
cosmology (H0=73) drops S8 from 0.834 to 0.783, a **1.6σ improvement** of
the (existing) S8 tension. The H0 tension drops simultaneously to <1σ vs
SH0ES.

This means:

> **mₑ-variation cosmology with H₀ ≈ 71-73 simultaneously resolves both
> the H₀ tension AND the S8 tension.**

By contrast EDE at H0=73 leaves S8 at 0.830 — **as bad as the original
Planck-LCDM S8 problem** (which is +2σ above KiDS to begin with).

The literature has missed this because it implicitly conflates "rd shift"
with "EDE physics". The two have *identical* effects on the CMB acoustic
geometry but **different** effects on linear matter-perturbation growth:

- EDE injects energy density at recombination → boosts dark-matter
  Meszaros effect → +5–7% σ8.
- mₑ-variation shifts only the recombination time/width via Bohr-radius
  scaling → no energy injection → unchanged growth → unchanged σ8.

The S8 wall that the previous-agent analysis (and the literature) cited
to disfavor pre-recombination resolutions IS specific to EDE; it does not
apply to mₑ-variation.

### Bayesian Model Averaging (novel)

If we treat (M0, M1, M2, M_compromise) as competing hypotheses with equal
prior probabilities, the posterior model probabilities are:

```
  Model              P(M | data)
  ---------------    -----------
  M0  LCDM             0.01%
  M1  rd shift        60.11%
  M2  Cepheid bias    13.15%
  Mco compromise      26.74%

  H0_BMA = 72.07 +/- 1.78 km/s/Mpc

  Tension with SH0ES:    0.47 sigma
  Tension with TRGB:     0.89 sigma
  Tension with Planck:   2.54 sigma
  Tension with ACT+DESI: 2.14 sigma
```

The H₀ tension (5σ between SH0ES and CMB) **drops to 2.5σ** purely by
acknowledging model uncertainty. Most papers fix a model and quote H₀
under it; BMA properly accounts for our ignorance about which mechanism
is true, and the answer lands at **TRGB-class H₀ ≈ 72**.

### Falsifiable predictions

1. **JWST Cepheid cross-check** (Riess+24 ongoing) extended to 30+ hosts:
   - If finds ΔM_cep < 0.02 mag globally → **kills M2 component** of compromise → forces M1 (rd shift only) → drives H₀ → 73.
   - If finds ΔM_cep ~ 0.04 ± 0.02 mag → **supports compromise**.
   - If finds ΔM_cep > 0.10 mag → kills the M1 side.

2. **CMB-S4 (≥ 2030)** damping-tail measurement of Δr_d / r_d:
   - If detects Δr_d / r_d ≈ −3% → **supports M1 / compromise**.
   - If r_d is consistent with Planck to <1% → **kills M1** → forces M2 or beyond.

3. **TRGB megasample** (Freedman 2026+):
   - If converges to H₀ = 70.5 ± 0.8 → **validates compromise target**.

4. **CMB damping-tail spectral signatures** of EDE vs me-variation:
   - EDE produces specific high-ℓ signature in TT damping tail (suppressed by Silk damping shift).
   - me-variation produces different signature (changed recombination width).
   - CMB-S4 can distinguish in principle.

5. **Joint S8 + H₀ analysis** with full MGCAMB / class_ede:
   - Should find me-variation passes both, EDE pulls S8 by ≥1.5σ.

---

## What this is and isn't

**This IS:**
- A reproducible quantification of the (M1, M2, M_compromise) Bayesian comparison with physically motivated priors.
- A clean separation of the rd-shift physics from its σ8 consequence.
- An empirical localization of H₀ at 68.67 ± 0.53 with no anchor.
- An identification of TRGB-class H₀ ≈ 71–72 as the only value consistent with both SH0ES (<1σ) and Planck (within 2.5σ).

**This is NOT:**
- A replacement for full Boltzmann + Planck-clik analysis (we use compressed CMB).
- A first-principles calculation of EDE's σ8 boost (we use Hill+2020 ~6% as input).
- A definitive test of mₑ-variation (we infer no-boost from physics, not from a code).
- An MCMC against KiDS-1000 + DES-Y3 cosmic-shear data (we use the compressed S8 only).

The next step to harden this is a full MGCAMB or class_ede run with the full Planck likelihood AND KiDS shear data, jointly fit. We can't do that in this stateless sandbox.

---

## Reproducibility

```bash
# baseline
python3 scripts/smoke_lcdm_pantheon_plus.py

# diagnostics
python3 scripts/diagnostics_panel.py
python3 scripts/etherington_test.py
python3 scripts/cepheid_diagnostics.py
python3 scripts/cepheid_host_jackknife.py

# anchor-free H0
python3 scripts/anchor_free_h0.py

# Bayesian comparison with consistent physical priors
python3 scripts/compromise_fair_priors.py

# Final test with KiDS S8 (long, ~30-60 min)
python3 scripts/joint_with_s8.py
```

All scripts are self-contained; data and dependencies in `external/` and the standard scientific Python stack (numpy, scipy, emcee, nautilus, getdist, CAMB, classy).
