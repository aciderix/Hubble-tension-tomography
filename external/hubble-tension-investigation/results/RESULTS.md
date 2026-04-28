# Numerical results

All numbers below are from one self-contained run on April 28, 2026, on the public Pantheon+SH0ES dataset (1701 SNe Ia, 1578 non-calibrators after cuts) and the 7 DESI BAO 2024 D_V/r_d data points.

## Test A — Early Dark Energy (the one that works)

```
LCDM (f_EDE=0)  : H0 = 66.96 +/- 0.50, Om = 0.288, chi^2 = 357.14
LCDM + EDE      : H0 = 71.09 +/- 0.20, Om = 0.288, f_EDE = 0.068 (6.8%)
                  r_d effective = 138.51 Mpc (vs 147.05 LCDM)

Tension vs SH0ES (73.04 +/- 1.04): 1.84 sigma
f_EDE required to reach H0 = 73.04 exactly: 9.8%
```

Comparison with literature (Poulin+ 2018, Smith+ 2020, Murgia+ 2021): they find f_EDE ~ 10-12% required to fully bridge the tension. Our 9.8% is consistent.

## Test B — Evolving dark energy w0waCDM

```
LCDM      : H0 = 63.51, Om = 0.309, chi^2 = 314.09
w0waCDM   : H0 = 62.27, Om = 0.309, w0 = -0.388, wa = -2.753, chi^2 = 288.77
            Delta(chi^2) = -25.31  (preference for w(z) =/= -1)
```

DESI 2024 reports w0 = -0.45 +/- 0.21, wa = -1.79 +/- 0.48. Our central values are in the same direction; wa is more negative because we use Pantheon+ rather than DESY5 SNe.

But H0 stays near 62 km/s/Mpc. Tension with SH0ES grows from ~5σ to ~9.5σ. Late-time evolving dark energy does not solve the Hubble tension.

## Test C — Distance-ladder internal consistency

```
H0 by redshift bin:
  z in [0.005, 0.020] : H0 = 71.27 +/- 0.67  (N=190)
  z in [0.020, 0.040] : H0 = 72.40 +/- 0.39  (N=318)
  z in [0.040, 0.060] : H0 = 73.32 +/- 0.67  (N=85)
  z in [0.060, 0.080] : H0 = 73.42 +/- 0.91  (N=49)
  z in [0.100, 0.150] : H0 = 73.30 +/- 0.66  (N=85)

Trend: dH0/dz = +14.5 +/- 7.2 km/s/Mpc per unit z (2.0 sigma)

Mass-step test:
  M* < 10  : H0 = 72.55 +/- 0.41 (N=342)
  M* >= 10 : H0 = 72.63 +/- 0.33 (N=399)
  Delta(H0) = +0.08 km/s/Mpc (0.16 sigma -- null result)
```

A genuine local underdensity ("KBC void") would produce a high H0 locally that DECREASES with redshift toward the cosmological value at the void boundary. The observed sign is the opposite, more consistent with residual peculiar-velocity contamination at the lowest redshifts than with a void.

## Summary table

| Test | Mechanism | H0 | Tension | Verdict |
|------|---|---|---|---|
| A | EDE (f_EDE = 6.8%) | 71.1 | 1.8σ | resolves |
| A | EDE target (f_EDE = 9.8%) | 73.04 | 0σ | resolves |
| B | w0waCDM | 62.3 | 9.5σ | worsens |
| C | Calibration consistency | n/a | n/a | no inconsistency |
