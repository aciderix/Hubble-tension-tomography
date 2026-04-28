# Hubble Tension Tomography — Working Notes for Claude

## Mission
Resolve (or sharply localise the failure of) the Hubble tension between
- **CMB / early-universe**: Planck 2018 ΛCDM → **H₀ = 67.4 ± 0.5** km/s/Mpc, ACT+DESI ≈ 68.2, DESI+CMB joint ≈ 67.97.
- **Late-universe ladder**: SH0ES (Riess+22) → **H₀ = 73.04 ± 1.04** km/s/Mpc; TDCOSMO ≈ 71.6.
- **Sound-horizon-free**: Freedman/CCHP varies by indicator (~70). DES Y5, Union3 SNe shift Ωₘ but not H₀.

The discrepancy is now ≥ 5σ. Either *something* in one ladder is wrong at the 1% level, or new physics shifts the inferred H₀ between epochs.

## Behavioural Rules (never drop these across compaction)

### Philosophy
Behave like an **exploratory system mapping unknown space**, not a mathematician
chasing a proof.
- Explore unknown, unexpected, counter-intuitive paths.
- Think like a machine, not like a human.
- Results in chat, not in report files (except very long content).
- Do **not** cheerleader — attack every positive result *immediately*.
- Falsify systematically and immediately.
- Document failures as clearly as successes.
- Do not blindly follow suggestions from other LLMs.
- Copy heavy data to `/tmp/` before processing (FUSE mount `/agent/` is slow).
- Pure numpy by default (sklearn does not compile on Alpine; numpy 1.26 is fine).

### 11 Mandatory Falsification Rules
1. **Linear baseline first** — before any sophisticated method.
2. **Fresh objects** for any validation (no peeking at the test set).
3. **Tautology filter** — drop trivial algebraic identities.
4. **Triviality test** — compare against random objects.
5. **Literature check** before any originality claim.
6. **Minimum complexity** — if a simple explanation fits, prefer it.
7. **Vary parameters** for any universality claim.
8. **Vary measurement method** for any algorithm-based metric.
9. **Decompose** the result into components before concluding.
10. **Null model** — if a null model gives the same answer, it is artefactual.
11. **CLT / U-statistics** check for any weak scaling law.

### Sweet-spot loop
generate → measure → patterns → falsify

## Repo layout

```
.
├── data/                         # Pantheon (2018) + Pantheon+SH0ES likelihood inputs
├── external/                     # bundled third-party data (light slices)
│   ├── A-No-Go-guide-for-the-Hubble-tension/   # Cai/Wang+ 2024, arXiv:2405.07039
│   ├── PantheonPlus/             # Brout+ 2022 distance/cov + cosmology chains (light)
│   ├── Pantheon/                 # Scolnic+ 2018 lcparam + sys
│   ├── SH0ES_Data/               # Riess+ 2022 Cepheid+TRGB anchor data
│   ├── TomoPost-BAO/             # Wang+ 2017 tomographic BAO
│   ├── cobaya_bao_data/          # all SDSS BOSS/eBOSS + DESI DR1/DR2 BAO
│   └── cobaya_sn_data/           # DES-Dovekie, DES Y5, Union3 likelihood files
├── src/
│   ├── cosmo.py                  # pure-numpy distances (E, dC, dL, μ, DV)
│   └── data_loaders.py           # all dataset loaders
├── scripts/
│   └── smoke_lcdm_pantheon_plus.py   # baseline ΛCDM fit reproducing SH0ES
├── fig1.nb / fig2-3.nb / fig4.nb # original Mathematica notebooks (Cai+ 2024)
└── CLAUDE.md                     # this file
```

## Toolchain installed
- `numpy 1.26.4`, `scipy 1.12`, `matplotlib`, `pandas`, `astropy`, `h5py`
- MCMC / sampling: `emcee 3.1.6`, `dynesty 3.0`, `nautilus-sampler 1.0.6`
- Bayesian framework: `cobaya 3.6.2`
- Boltzmann solvers: **CAMB 1.6.6** *and* **CLASS (`classy 3.3.4`)** — both compiled and importable
- Posterior tools: `getdist 1.7.6`, `corner 2.2.3`
- Gaussian processes: `GPy 1.13.2`

GAP (the discrete-algebra system) is **not** installed — almost surely irrelevant
to cosmology. If it is needed for a specific algebraic check, install via
`apt-get install gap`.

## Baseline anchor (verified 2026-04-28)
Running `python3 scripts/smoke_lcdm_pantheon_plus.py`:

| dataset | H₀ | Ωₘ | M | χ² | dof |
|---|---|---|---|---|---|
| Pantheon+SH0ES (with Cepheid calibrators) | **73.42** | 0.334 | -19.25 | 1522.4 | 1698 |
| Pantheon+ (no calibrators) | *unconstrained* | 0.334 | M-5log₁₀H₀ degenerate | 1456.4 | 1621 |

Δχ² profile gives 1σ ≈ ±0.55 — consistent with Brout+22 / Riess+22 within rounding.
This anchors any modification: **a successful resolution must move H₀ down to ~67
without breaking this fit**, or move CMB-inferred H₀ up to ~73 without breaking
DESI BAO + CMB acoustic-scale fits.

## Working principles for this project
- Every claim about a model must be falsified against the table above first.
- Linear / single-parameter modifications come before exotic physics.
- A fit improvement of Δχ² < 6.18 over ΛCDM with 2 extra params is **not** a
  detection — that is the 95 %-CL Δχ² for two free parameters.
- Mixing local + CMB likelihoods requires a Boltzmann solver (use CAMB), not just
  the pure-numpy distance code in `src/cosmo.py`.
