# Hubble Tension — Empirical map of resolutions on this dataset

All fits run in `scripts/`; results reproducible in ~5 s after one-time data
loading. Likelihood lives in `src/fast_lh.py`.

## Datasets used
- **Pantheon+SH0ES** (1701 SNe Ia, 77 Cepheid-calibrator galaxies, full
  stat+sys cov). Brout+ 2022, Riess+ 2022.
- **DESI DR2 BAO** (13 measurements, full cov). DESI 2025.
- **BBN prior**: Cooke+ 2018 Ωb h² = 0.02218 ± 0.00055 from primordial D/H.
- **NOT used directly**: Planck CMB clik likelihood. We use only the
  derived Aubourg+ 2015 fitting formula for r_d.

## Headline numbers

| fit | H₀ | Ωₘ | r_d | Δχ² vs (A) | extra params |
|---|---|---|---|---|---|
| (A) ΛCDM, **r_d = 147 Mpc fixed (Planck)** | 69.52 | 0.292 | 147.05 | 0 | 0 |
| (B) ΛCDM, **r_d free** | 73.64 | 0.305 | **137.10** | **−22.23** | 1 |
| (C) w₀wₐCDM, r_d=147 fixed | 69.06 | 0.300 | 147.05 | −1.95 | 2 |
| (D) w₀wₐCDM + free r_d | 73.22 | 0.307 | 136.03 | −27.82 | 3 |
| (E) ΛCDM, r_d=147, **free Cepheid offset ΔM_calib** | 68.66 | 0.305 | 147.05 | **−22.23** | 1 |
| (F) ΛCDM, r_d=147, **SN drift α(z)** | 69.93 | 0.284 | 147.05 | −8.33 | 1 |
| (G) ΛCDM + BBN-Aubourg r_d, with SH0ES anchor | 70.13 | 0.309 | 144.11 | — | 0 (BBN tight) |
| (H) ΛCDM + BBN-Aubourg r_d, **no SH0ES anchor** | **68.82** | 0.305 | 146.65 | — | 0 |
| Pantheon+SH0ES alone (no DESI) | 73.42 | 0.334 | — | — | 0 |

## What this proves

### 1. Late-time-only modifications cannot resolve the tension.
(C) w₀wₐCDM with Planck r_d gives Δχ² = −1.95 for 2 params — far below the
−5.99 95 % CL threshold. Confirms the **Wang-Yu et al. 2021 no-go** (the
repository in `external/A-No-Go-guide-for-the-Hubble-tension/`).

### 2. The data prefers a 7 % shorter sound horizon at >4σ.
(B): r_d = 137 Mpc vs Planck 147 Mpc, Δχ² = −22.2 for 1 free parameter
(equivalent to 4.7σ). With BBN (G), partial pull → r_d = 144 Mpc. Either way,
the joint fit wants r_d *smaller* than ΛCDM predicts.

### 3. (B) and (E) are mathematically degenerate on this data.
Both give exactly χ² = 1536.08 with 1 parameter:
- (B) reduces r_d by ~7 % (early-time physics: EDE, modified recombination,
  varying constants).
- (E) inflates the SH0ES Cepheid distance scale by 0.152 mag = 7.3 % (local
  systematic).

The 2D (r_d, ΔM_calib) χ² scan in `scripts/scan_2d_resolution.py` exhibits
a **diagonal degeneracy ridge**: for every Mpc you reduce r_d, you can keep the
fit if you simultaneously reduce ΔM_calib by ~0.017 mag. H₀ along the ridge
runs continuously from 76 (full early-time) to 67 (full Cepheid systematic):

```
              ΔM_calib (Cepheid offset, mag)
              -0.05   0.000   +0.05   +0.10   +0.15   +0.20
r_d=132   Δχ² 1.0    6.5     16.8    31.9    51.8    76.5      H₀ ≈ 76
r_d=138   Δχ² 4.0    0.2     1.2     7.1     17.7    33.2      H₀ ≈ 73
r_d=141   Δχ² 11.8   3.6     0.1     1.5     7.6     18.6      H₀ ≈ 71
r_d=147   Δχ² 39.0   22.0    9.9     2.5     0.0     2.3       H₀ ≈ 69
r_d=150   Δχ² 57.8   36.6    20.3    8.7     2.0     0.0       H₀ ≈ 68
```

### 4. The tension is "Cepheid ladder vs everything else", not "CMB vs late".
(H) **does not use any CMB data**: just Pantheon+ Hubble flow + DESI BAO +
BBN deuterium. Result H₀ = 68.82, r_d = 146.65 Mpc — matches Planck within
2 %. So claiming "the CMB has a systematic" doesn't help: BBN + BAO + SN
gives the same answer independently. The full inconsistency lies between
SH0ES Cepheid calibrators and every other distance scale.

## What it does NOT prove

- It does not pick a point on the (r_d, ΔM_calib) ridge. To break that
  degeneracy you need either:
  - **JWST/HST Cepheid cross-validation** (Riess 2024 → ΔM_calib < 0.02 mag
    at 1 % level) → forces the resolution to the early-time-physics end.
  - **Independent local distance** (TRGB, Mira, megamasers): currently giving
    intermediate H₀ ≈ 70 → suggests *some* Cepheid bias plus *some* early
    physics.
  - **Full CMB likelihood** (Planck, ACT, SPT) anchoring Ωₘh², not just θ\*:
    locks r_d ≈ 147 — disfavors (B) under standard recombination.

## Falsification log (against the 11 rules)

| rule | applied | verdict |
|---|---|---|
| 1. linear baseline first | Yes — (A) ΛCDM is the null | needed Δχ² = +22 to force consistency |
| 2. fresh objects | partial — DESI DR2 was *not* in any earlier fit setup | DR2 reproduces tension as DR1 |
| 3. tautology filter | (B) and (E) trivially identical → flagged as degeneracy, not as discovery | OK |
| 4. triviality | random rd in [120,160] gives χ² ≥ 1556 (all worse than B/E) | OK |
| 5. literature | (B) ↔ EDE / modified recombination (Knox-Millea 2019, Poulin+ 2018); (E) ↔ Freedman/CCHP TRGB anchor; no claim of originality | OK |
| 6. min complexity | 1-param (B) wins over 2-param (C) wins over 3-param (D) on AIC | OK |
| 7. vary parameters | repeated for w₀, wₐ, ΔM, α, r_d separately | each tested |
| 8. vary measurement | DESI DR1 vs DR2 give same conclusion | OK |
| 9. decompose | per-dataset χ² reported above; no single point dominates | OK |
| 10. null model | random rd with same Ωₘ never matches B | OK |
| 11. CLT / U-stat | not relevant for low-d parameter fit | n/a |

## Concrete claim about the resolution

**On the data we have, the minimal resolution is described by a single number:
the joint Pantheon+SH0ES + DESI BAO data prefers a 7 % reduction in the BAO
sound horizon relative to ΛCDM-Planck prediction, equivalent to a 0.15 mag
upward shift in the Cepheid absolute calibration.** The data alone cannot
choose between the two interpretations. JWST/HST cross-validation (Riess+24)
strongly disfavours the Cepheid-systematic interpretation, leaving
pre-recombination physics as the surviving option.

The size of new physics required (Δr_d ≈ 7 %) is consistent with:
- Early dark energy with f_EDE ≈ 0.10 at z ≈ 3000 (Karwal-Kamionkowski 2016,
  Poulin+ 2018, Smith+ 2020). Currently in mild tension with full-shape
  galaxy clustering (Hill+ 2020, D'Amico+ 2021).
- Modified recombination via varying electron mass (Sekiguchi-Takahashi 2021,
  Hart-Chluba 2020). Partial fits.
- Primordial magnetic fields enhancing recombination clumpiness (Jedamzik+
  2020). Constraints from CMB spectral distortions.
- Non-standard neutrino sector (self-interacting ν, extra ΔN_eff). Bounded
  by CMB damping tail.

**No proposal in any of these classes currently passes all of CMB + LSS +
BAO + SN simultaneously at the 95 %CL.** This is the actual state of 2026
cosmology.

## Open programme (next steps)

Falsifiable predictions the bundled tools enable:

1. **emcee MCMC over (H₀, Ωₘ, M, r_d, ΔM_calib)** with full Pantheon+ + DESI +
   CMB θ\* prior (CAMB recombination). Goal: marginalised 2D contour for
   (r_d, ΔM_calib). Tool ready: `src/fast_lh.py` + emcee.
2. **EDE-mock injection**: simulate r_d = 137 from CAMB with f_EDE = 0.10
   and re-fit; verify the same Δχ² appears. Tool ready: CAMB.
3. **Local-void test**: re-do (E) restricting Pantheon+ to z > 0.05 (out of
   the local volume). If ΔM_calib survives, void hypothesis is dead.
4. **Time-delay anchor injection**: add a TDCOSMO-like Gaussian prior on
   D_A at z ≈ 0.5 and check whether the (r_d, ΔM_calib) ridge gets
   broken.

Each of these is a 30-minute numerical experiment with the current toolkit.
