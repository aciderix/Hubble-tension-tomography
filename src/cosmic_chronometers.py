"""Cosmic chronometers H(z) likelihood (Moresco compilation).

H(z) measurements from differential ages of passively evolving early-type
galaxies — model-independent: doesn't assume any cosmology, doesn't use a
distance ladder, doesn't need rd or Cepheids.

Data: 15 points from Moresco+ 2012, 2015, 2016 (file MM_BC03 in
mmoresco/CCcovariance gitlab). z range 0.18 to 1.97.

Errors are split into stat_contr (statistical) and met_contr (metallicity
systematic). For diagonal-only treatment we add them in quadrature, which
gives the same total variance as `errHz`. Full covariance treatment per
Moresco+ 2020 prescription is left for future work; the off-diagonal
correlations *inflate* the effective error so the diagonal fit gives a
slightly tighter constraint.
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


def chi2_cc(H0, Om):
    """Diagonal chi2 against Moresco H(z) under flat LCDM."""
    z = CC_DATA["z"]
    H_pred = H0 * E_LCDM(z, Om)
    d = CC_DATA["H"] - H_pred
    return float(((d / CC_DATA["sig"]) ** 2).sum())


def chi2_cc_w0wa(H0, Om, w0, wa):
    z = CC_DATA["z"]
    H_pred = H0 * E_w0wa(z, Om, w0, wa)
    d = CC_DATA["H"] - H_pred
    return float(((d / CC_DATA["sig"]) ** 2).sum())
