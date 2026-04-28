# CLAUDE_HANDOFF — Updated April 28, 2026 (Level 5 erratum + Steps 3-4 roadmap)

## ⚠️ Read LEVEL5_ERRATUM.md FIRST

The "f(R) wins on H0+S8" conclusion from Level 4 (Step D) was based on a naive σ8 scaling. When σ8 is computed properly via the growth ODE + Eisenstein-Hu transfer (Level 5 Step E), the result inverts: f(R) at Ωm=0.260 gives S8≈0.59 (not 0.78), which is −9σ vs KiDS/DES.

## Honest current status

- **Level 1–3:** Background fits valid (Pantheon+ + DESI + CC + SH0ES anchor reject ΛCDM at Δχ²~24).
- **Level 4 Step B:** Valid. r_d shrinks to 139.78 Mpc agnostically — strong EDE-like signature.
- **Level 4 Steps C, D:** Background fits valid; σ8/S8 column is WRONG.
- **Level 5 Steps E, F, G:** σ8 corrected. Late-time MG cannot solve H0+S8 simultaneously with Ωm≈0.31.

## Real thesis to defend

EDE-like pre-recombination physics is preferred for H0. Tradeoff: ~2-3σ pull on S8 (modest, not catastrophic). This matches the published consensus (Hill+2020, Murgia+2021, Smith+2022).

---

## Steps 3-4 roadmap — what we still want to test (NOT done in sandbox)

### Step 3 — Exotic late-time models

**Why not done:** These require modified Boltzmann solvers (CLASS variants) that need C/Fortran compilation. Stateless sandbox can't compile/link CLASS.

**Models to test:**
1. **Interacting Dark Energy (IDE)** — coupled DM-DE Q = ξ·H·ρ_DM
   - Code: `class_idecoupling` (Di Valentino+2020 fork)
   - Free params: ξ (coupling), w_DE
   - Goal: see if DM→DE energy transfer can ease H0+S8 jointly
2. **Decaying DM (DDM)** — DM decays to dark radiation
   - Code: `class_DDM` (Audren+2014, Pandey+2020)
   - Free params: Γ (decay rate), f (decaying fraction)
   - Goal: late-time DM depletion → less structure growth → lowers S8 while preserving early Universe
3. **Negative DE density (sign-flip DE)** — Akarsu+2023
   - Code: standard CLASS with custom w(z)
   - Free param: z_flip
   - Goal: 5σ preference reported, worth replicating

**Joint fit:** cobaya MCMC with Pantheon+1701 + DESI DR2 BAO + Planck high-ℓ TT,TE,EE + KiDS-1000 + DES Y3.

### Step 4 — CMB anisotropy / dipole

**Why not done:** Needs full Planck SMICA maps (~50 GB) + healpy stack with custom masking, not feasible in sandbox.

**Tests to run:**
1. **Wang-Mota dipole** — directional H0 in Pantheon+ subsamples (Migkas+2021 reported ~5° anisotropy at 4σ in galaxy clusters)
   - Reanalyze with Pantheon+ split by hemispheres
2. **CMB-frame vs heliocentric** dipole consistency check (Secrest+2021 reported tension with quasars)
3. **Planck high-ℓ residuals in supercluster directions** (looking for void-induced ISW)

### Local machine setup script

Save as `setup_local.sh`:

```bash
#!/usr/bin/env bash
set -e

# 1. Dependencies (Ubuntu/Debian; for macOS use brew)
sudo apt-get install -y build-essential gfortran git cython3 python3-dev \
                        libfftw3-dev libgsl-dev libcfitsio-dev healpy

# 2. Python env
python3 -m venv hubble_env
source hubble_env/bin/activate
pip install numpy scipy matplotlib emcee getdist cobaya healpy \
            astropy camb classy

# 3. AxiCLASS (CLASS with EDE)
git clone https://github.com/PoulinV/AxiCLASS.git
cd AxiCLASS && make -j 4 && cd python && pip install -e . && cd ../..

# 4. MGCAMB (CAMB with f(R), DGP, etc)
git clone https://github.com/sfu-cosmo/MGCAMB.git
cd MGCAMB && python setup.py install && cd ..

# 5. class_ede (alternative EDE implementation)
git clone https://github.com/mwt5345/class_ede.git
cd class_ede && make -j 4 && cd python && pip install -e . && cd ../..

# 6. Cobaya likelihoods
cobaya-install cosmo -p ./packages

# 7. Data
mkdir -p data
cd data
wget https://github.com/PantheonPlusSH0ES/DataRelease/raw/main/Pantheon%2B_Data/4_DISTANCES_AND_COVAR/Pantheon%2BSH0ES.dat
wget https://github.com/PantheonPlusSH0ES/DataRelease/raw/main/Pantheon%2B_Data/4_DISTANCES_AND_COVAR/Pantheon%2BSH0ES_STAT%2BSYS.cov
wget https://data.desi.lbl.gov/public/dr2/vac/dr2/cosmology/v2.1/desi_dr2_bao.dat
cd ..

echo "Setup complete. Activate env: source hubble_env/bin/activate"
```

### Cobaya YAML for joint MCMC (`hubble_joint_mcmc.yaml`)

```yaml
likelihood:
  pantheonplusshoes.PantheonPlusSH0ES:
    path: ./data/Pantheon+SH0ES.dat
    covmat: ./data/Pantheon+SH0ES_STAT+SYS.cov
  bao.desi_dr2:
    path: ./data/desi_dr2_bao.dat
  planck_2018_lowl.TT_native: null
  planck_2018_lowl.EE_native: null
  planck_2018_highl_plik.TTTEEE_lite_native: null
  planck_2018_lensing.native: null
  des_y3_clustering_and_galaxy_galaxy_lensing: null
  kids1000_3x2pt: null

theory:
  classy:
    extra_args:
      Omega_Lambda: 0  # use fluid for DE
      Omega_fld: 0.6889
      Omega_scf: -1
      attractor_ic_scf: 'no'
      scf_potential: 'axion'
      n_axion: 3.0
      log10_axion_ac: -3.55
      scf_parameters: '2.83, 0'
      adptative_stepsize: 100
      scf_tuning_index: 0
      do_shooting: 'yes'
      do_shooting_scf: 'yes'

params:
  H0:
    prior: {min: 60, max: 80}
    proposal: 0.5
    latex: H_0
  omega_b:
    prior: {min: 0.02, max: 0.025}
    proposal: 0.0001
    latex: \omega_b
  omega_cdm:
    prior: {min: 0.10, max: 0.14}
    proposal: 0.001
    latex: \omega_{cdm}
  log10_axion_ac:
    prior: {min: -4, max: -3}
    proposal: 0.05
    latex: \log_{10}(a_c)
  fEDE:
    prior: {min: 0.001, max: 0.30}
    proposal: 0.01
    latex: f_{EDE}

sampler:
  mcmc:
    Rminus1_stop: 0.05
    max_samples: 200000

output: chains/hubble_joint_EDE
```

### Run command

```bash
cobaya-run hubble_joint_mcmc.yaml -r  # resumable run
# When done:
getdist-gui chains/hubble_joint_EDE  # plot triangle plots
```

---

## What survival looks like → preprint

If the joint MCMC gives **f_EDE > 0.05 at >3σ** AND **S8 within 2σ of KiDS/DES** AND **Δχ² < −10 vs ΛCDM**, that's a publishable result. Target: arxiv astro-ph.CO, ~10-15 pages.

Skeleton sections:
1. Intro: H0 tension status quo (cite Riess 2024, DESI 2024, ACT DR6 2024)
2. Methods: AxiCLASS EDE + cobaya joint MCMC
3. Results: corner plots, Δχ² table, σ8/S8 prediction vs KiDS/DES
4. Discussion: comparison to Hill+2020, Murgia+2021, recent EDE literature
5. Conclusions: EDE preferred at Xσ, residual S8 tension at Yσ

---

## Repo state (Apr 28, 2026)

- LEVEL4_RESULTS.md (background fits valid, σ8 column WRONG - see erratum)
- LEVEL5_ERRATUM.md (critical correction)
- level4_data.json
- level5_data.json (Step E numbers)
- CLAUDE_HANDOFF.md (this file)

Public GitHub: https://github.com/aciderix/hubble-tension-investigation
