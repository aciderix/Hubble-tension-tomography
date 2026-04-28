"""Cosmic chronometers H(z) likelihood (Moresco compilation).

H(z) measurements from differential ages of passively evolving early-type
galaxies — model-independent: doesn't assume any cosmology, doesn't use a
distance ladder, doesn't need rd or Cepheids.

Data: 15 points from Moresco+ 2012, 2015, 2016 (file MM_BC03 in
mmoresco/CCcovariance gitlab). z range 0.18 to 1.97.

Two likelihoods provided:
  - chi2_cc(H0, Om):    DIAGONAL errors only.
  - chi2_cc_full(H0, Om): FULL covariance per Moresco+ 2020 prescription
                           (diagonal + IMF systematic + SPS one-of-one).
"""
import numpy as np
from .cosmo import E_LCDM, E_LCDM_with_rad, E_w0wa


# (z, H, sigma_H)  -- Moresco BC03 stellar-library version
CC_DATA = np.array([
    (0.1791,  74.91,  3.8069262),
    (0.1993,  74.96,  4.9001352),
    (0.3519,  82.78, 13.94843),
    (0.3802,  83.0,  13.54),
    (0.4004,  76.97, 10.18),
    (0.4247,  87.08, 11.24),
    (0.4497,  92.78, 12.9),
    (0.4783,  80.91,  9.044),
    (0.5929, 103.8,  12.49752),
    (0.6797,  91.6,   7.961872),
    (0.7812, 104.5,  12.19515),
    (0.8754, 125.1,  16.70085),
    (1.037,  153.7,  19.6736),
    (1.363,  160.0,  32.63),
    (1.965,  186.5,  49.58),
], dtype=[("z", "f8"), ("H", "f8"), ("sig", "f8")])

# Systematic components from Moresco+ 2020 (data_MM20.dat). Each value is
# the percent contribution to the H(z) error from that systematic, vs a
# fiducial H_fid(z). We interpolate to our 15 z values.
_SYS_TABLE = np.array([
    # z,    imf%, slib%, sps%, spsooo%
    (0.075, 0.47, 7.40, 15.86, 9.91),
    (0.125, 0.47, 7.40, 14.23, 6.98),
    (0.175, 0.47, 7.40, 13.34, 5.40),
    (0.225, 0.47, 7.40, 13.21, 5.40),
    (0.275, 0.47, 7.40, 13.29, 5.40),
    (0.325, 0.47, 7.40, 12.20, 5.40),
    (0.375, 0.47, 7.40, 12.99, 5.40),
    (0.425, 0.47, 7.40, 10.29, 6.20),
    (0.475, 0.46, 7.39,  8.91, 5.86),
    (0.525, 0.23, 7.40,  9.99, 6.51),
    (0.575, 0.28, 6.87, 10.09, 6.12),
    (0.625, 0.47, 6.65, 11.17, 6.21),
    (0.675, 0.47, 6.57, 11.12, 5.71),
    (0.725, 0.47, 5.90, 10.81, 5.16),
    (0.775, 0.45, 6.03, 10.75, 5.05),
    (0.825, 0.47, 6.10, 10.75, 5.05),
    (0.875, 0.47, 5.89,  9.08, 2.79),
    (0.925, 0.44, 5.80,  8.62, 3.70),
    (0.975, 0.40, 5.94,  7.32, 3.65),
    (1.025, 0.27, 6.07,  5.84, 3.37),
    (1.075, 0.20, 6.08,  6.02, 3.49),
    (1.125, 0.20, 6.07,  4.72, 2.33),
    (1.175, 0.19, 6.09,  4.31, 2.33),
    (1.225, 0.19, 6.09,  3.90, 2.33),
    (1.275, 0.19, 6.09,  3.90, 2.33),
    (1.325, 0.20, 6.09,  3.91, 2.34),
    (1.375, 0.19, 6.09,  3.90, 2.34),
    (1.425, 0.19, 6.09,  3.90, 2.33),
    (1.475, 0.20, 6.09,  3.91, 2.34),
])


def _build_full_cov():
    z = CC_DATA["z"]
    H = CC_DATA["H"]
    sig = CC_DATA["sig"]
    n = len(z)
    # Interpolate systematic fractions at our z's (cap to nearest at high z)
    z_sys = _SYS_TABLE[:, 0]
    imf_pct = np.interp(z, z_sys, _SYS_TABLE[:, 1]) / 100.0
    spsooo_pct = np.interp(z, z_sys, _SYS_TABLE[:, 4]) / 100.0
    cov = np.diag(sig ** 2)  # diagonal stat (already includes met)
    # Add fully-correlated systematic outer products
    delta_imf = H * imf_pct
    delta_sps = H * spsooo_pct
    cov += np.outer(delta_imf, delta_imf)
    cov += np.outer(delta_sps, delta_sps)
    return cov


_FULL_COV = _build_full_cov()
_FULL_INV = np.linalg.inv(_FULL_COV)


def chi2_cc(H0, Om):
    """Diagonal chi2 against Moresco H(z) under flat LCDM."""
    z = CC_DATA["z"]
    H_pred = H0 * E_LCDM(z, Om)
    d = CC_DATA["H"] - H_pred
    return float(((d / CC_DATA["sig"]) ** 2).sum())


def chi2_cc_full(H0, Om):
    """Chi2 with full Moresco+2020 covariance (diagonal + IMF + SPS_ooo
    correlated systematic). Inflated relative to diagonal."""
    z = CC_DATA["z"]
    H_pred = H0 * E_LCDM(z, Om)
    d = CC_DATA["H"] - H_pred
    return float(d @ _FULL_INV @ d)


def chi2_cc_w0wa(H0, Om, w0, wa):
    z = CC_DATA["z"]
    H_pred = H0 * E_w0wa(z, Om, w0, wa)
    d = CC_DATA["H"] - H_pred
    return float(((d / CC_DATA["sig"]) ** 2).sum())
