# External datasets bundled in this repo

All material below is included under each project's original licence; cite
the upstream paper, not this fork.

| folder | upstream | what / why |
|---|---|---|
| `A-No-Go-guide-for-the-Hubble-tension/` | [Wang-weiYu/A-No-Go-guide-for-the-Hubble-tension](https://github.com/Wang-weiYu/A-No-Go-guide-for-the-Hubble-tension) — Cai, Guo, Wang, Yu, Zhou, [arXiv:2107.13286](https://arxiv.org/abs/2107.13286) | Their `code/page.py` runs an emcee fit of a 2-param late-time deformation against SN+BAO+H(z); useful as a known-bad-fit baseline. Heavy MCMC chains were dropped — regenerate with `python3 external/A-No-Go-guide-for-the-Hubble-tension/code/page.py`. |
| `PantheonPlus/` | [PantheonPlusSH0ES/DataRelease](https://github.com/PantheonPlusSH0ES/DataRelease) — Brout+ 2022, Riess+ 2022 | `Pantheon+SH0ES.dat` distance vector + Cobaya `5_COSMOLOGY/` posterior chains and configs. The 442 MB systematic-grouping covariances and 365 MB SALT2 trainsets were *not* copied — fetch from upstream if needed. |
| `Pantheon/` | [dscolnic/Pantheon](https://github.com/dscolnic/Pantheon) — Scolnic+ 2018 | `lcparam_full_long*.txt` and `sys_full_long.txt` — the Wang et al. 2024 (arXiv:2405.07039) likelihood reads these. |
| `SH0ES_Data/` | [PantheonPlusSH0ES/DataRelease](https://github.com/PantheonPlusSH0ES/DataRelease) | Cepheid + TRGB anchor data behind H₀ = 73.04. `run_mcmc.py` reproduces the SH0ES posterior. |
| `TomoPost-BAO/` | [ytcosmo/TomoPost-BAO](https://github.com/ytcosmo/TomoPost-BAO) — Wang+ 2017, [arXiv:1607.03154](https://arxiv.org/abs/1607.03154) | Tomographic BAO covariance referenced by Cai+ 2024 (`tomographic_BAO_invcov.txt`). |
| `cobaya_bao_data/` | [CobayaSampler/bao_data](https://github.com/CobayaSampler/bao_data) | All BAO likelihood data files used by the standard Cobaya `bao.*` modules — SDSS DR12, DR16, eBOSS, **DESI DR1 (2024)**, **DESI DR2 (2025)**. |
| `cobaya_sn_data/` | [CobayaSampler/sn_data](https://github.com/CobayaSampler/sn_data) | DES Y5, DES-Dovekie, Union3 SNe likelihood files (Pantheon+ already at top-level `data/`). |

## What was deliberately *not* bundled

- Planck 2018 baseline likelihood (heavy, lives in the `clik` framework — install
  separately if you want a CMB likelihood; otherwise use the Planck-derived
  Gaussian compressed likelihoods in `cobaya_bao_data/sdss_DR16_*.dat`).
- DESI DR2 raw galaxy catalogues — the BAO-only compressed likelihood is enough
  for cosmological inference and is in `cobaya_bao_data/desi_bao_dr2/`.
- Pantheon SALT2 trainset and Pantheon+ systematic-grouping covariances — only
  needed if you reprocess lightcurves, not for cosmology.
