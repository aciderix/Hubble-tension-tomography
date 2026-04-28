# Hubble Tension Investigation

A pedagogical computational study of the **H₀ tension** — the >5σ
discrepancy between local (SH0ES: 73.04 ± 1.04 km/s/Mpc) and CMB
(Planck: 67.4 ± 0.5 km/s/Mpc) measurements of the Hubble constant.

This repository documents an end-to-end original analysis using public
data (Pantheon+, DESI BAO, Planck CMB, KiDS-1000 / DES Y3 S₈) and tests
four candidate mechanisms:
- **(A) Early Dark Energy** — energy injection before recombination
- **(B) Late-time w(z) evolution** — modified dark-energy equation of state
- **(C) Calibration bias** — systematic offset in distance-ladder rungs
- **(D) Local void** — anisotropic local underdensity (KBC-like)

## Headline results

| Mechanism | Verdict |
|-----------|---------|
| (A) EDE   | **Viable**. Toy parameterization yields f_EDE ≈ 9.8%, Δχ²=−33 vs ΛCDM (>5σ preference) when SH0ES anchor is included |
| (B) w(z)  | **Rejected**. No viable region in (w₀, wₐ) space resolves tension while respecting BAO+CMB |
| (C) Calibration | **Rejected**. Per-survey and drop-one tests show signal is distributed |
| (D) Local void | **Rejected**. Inconsistent with Pantheon+ residual structure at z>0.05 |

## Repository structure

```
docs/
  level0_pedagogical_walkthrough.md   ← Read this first
  level1_anisotropy_diagnosis.md      ← Per-survey + drop-one robustness
  level2_boltzmann_validation.md      ← classy/CAMB validation + ΔN_eff scan
scripts/
  level1_*.py                          ← Mechanism A–D quantitative tests
  anisotropy_*.py                      ← Anisotropy diagnostics
  level2_neff_scan.py                  ← Boltzmann ΔN_eff stress test
results/
  *.json                               ← Raw analysis outputs
```

## Honest scientific status

What this repository **does**:
- Validate the toy-EDE Δχ²=−33 result with a simplified `r_d → r_d(1−α·f_EDE)` parameterization
- Diagnose that residual low-z anisotropy is real but sub-publication-grade (2.5σ)
- Demonstrate that pure ΔN_eff dark radiation (the simplest pre-recombination
  energy injection) is excluded by Planck high-ℓ at >5σ — motivating axion-like EDE

What this repository **does not** do:
- Run a real MCMC against the Planck likelihood (requires CLASS_EDE/AxiCLASS forks)
- Compute model evidence for axion EDE vs ΛCDM
- Fit DESI DR2 BAO data

These are the natural next steps for a full publication-grade analysis.

## Reproducibility

All scripts are self-contained; data sources documented in each script header.
Required: Python 3.9+, numpy, scipy, pandas, emcee, classy 3.3.4+, camb 1.6+.

## License
MIT. See `LICENSE`.
