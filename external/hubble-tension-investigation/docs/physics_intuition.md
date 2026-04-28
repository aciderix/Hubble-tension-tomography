# Physics intuition: why Early Dark Energy works

## The key identity

The CMB does not measure H0 directly. What it measures, with extreme precision (~0.03%), is the angular acoustic scale at recombination:

    theta_* = r_s(z_*) / D_A(z_*)

where r_s is the comoving sound horizon at recombination and D_A is the angular diameter distance to the last-scattering surface. In LCDM, both r_s and D_A scale roughly as 1/H0 at fixed Om, so theta_* essentially fixes the *combination*:

    H0 * r_d ~ constant (from CMB)

Planck 2018 reports H0 * r_d = 9912 +/- 30 km/s.

## How EDE shifts H0

Suppose we add a new component that contributes a fraction f_EDE of the total energy density at recombination, then quickly dilutes (e.g. an axion-like field oscillating in a (1-cos phi)^n potential). At recombination, the expansion rate H(z_*) is *higher* than LCDM, so the sound has *less time* to propagate before photon decoupling. Therefore:

    r_d^EDE < r_d^LCDM

To preserve the well-measured CMB observable H0*r_d, the inferred H0 must scale as:

    H0^EDE = H0^LCDM * (r_d^LCDM / r_d^EDE)

A 5% shorter r_d -> 5% larger H0: 67 -> 70.4. A ~10% shorter r_d gets us to 73.

In our toy parameterization we wrote:

    r_d = r_d^LCDM * (1 - alpha * f_EDE)

with alpha = 0.85 (calibrated against full-Boltzmann EDE solutions in the literature), so f_EDE = 0.10 implies r_d shrinks by 8.5% implies H0 rises by 8.5% to ~73. Our fit returns f_EDE = 9.8% to hit 73.04 exactly, consistent with Poulin+ 2018 and Smith+ 2020.

## Why this is not a closed case

EDE has well-known costs:

1. **S8 tension.** Adding energy at recombination changes the matter power spectrum. EDE tends to *increase* sigma_8, worsening the existing mild tension between Planck-inferred sigma_8 and weak-lensing measurements of S8.

2. **No fundamental motivation.** Why does an axion-like field happen to oscillate exactly at z ~ 3500 with f ~ 10%? Mass and decay constant are tuned a posteriori.

3. **Variants exist.** NEDE (first-order phase transition rather than oscillation), AdS-EDE, and others can reduce some of these tensions but no single proposal has consensus.

4. **Other resolutions remain on the table.** Modified recombination, interacting dark sectors, and certain modified-gravity proposals can also shorten r_d. EDE is the most-studied but not unique.

So: EDE is currently the leading candidate in the H0 tension literature and our laptop-scale reproduction confirms it works numerically. But "works numerically" is not "fundamental physical explanation".
