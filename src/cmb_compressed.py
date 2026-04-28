"""Planck 2018 compressed CMB likelihood (distance priors).

Reference: Chen, Huang, Wang 2019 (arXiv:1808.05724), Planck 2018 TTTEEE+lowE.

Constrains (R, l_A, omega_b) at z* = 1089.92:

  R     = sqrt(Omega_m H0^2) * D_M(z*) / c    (shift parameter)
  l_A   = pi * D_M(z*) / r_s(z*)              (acoustic angular scale)
  omega_b = Omega_b h^2

This is independent of low-z H0; it ties together early-universe physics
with the integrated geometry to z*=1090. Adding it to BAO+SN breaks the
(r_d, calibration-offset) degeneracy seen in the SN+BAO-only analysis.
"""
import numpy as np
from scipy.integrate import quad

from .cosmo import E_LCDM_with_rad, C_KMS, rd_planck_anchored

# Planck 2018 distance priors (TTTEEE+lowE+lensing+BAO consensus, Chen+ 2019)
PLANCK_R = 1.7502
PLANCK_LA = 301.471
PLANCK_OBH2 = 0.02236
PLANCK_SIG_R = 0.0046
PLANCK_SIG_LA = 0.090
PLANCK_SIG_OBH2 = 0.00015
PLANCK_CORR = np.array([
    [1.00,  0.46, -0.66],
    [0.46,  1.00, -0.34],
    [-0.66, -0.34, 1.00],
])
PLANCK_SIG = np.array([PLANCK_SIG_R, PLANCK_SIG_LA, PLANCK_SIG_OBH2])
PLANCK_COV = PLANCK_CORR * (PLANCK_SIG[:, None] * PLANCK_SIG[None, :])
PLANCK_INVCOV = np.linalg.inv(PLANCK_COV)
PLANCK_MEAN = np.array([PLANCK_R, PLANCK_LA, PLANCK_OBH2])

Z_STAR = 1089.92  # Planck 2018 best-fit redshift of last scattering


def comoving_distance_to_zstar(H0, Om, h=None, Neff=3.046):
    """D_C(z*) including radiation (necessary for accuracy at z*~1100)."""
    if h is None:
        h = H0 / 100.0
    val, _ = quad(lambda z: 1.0 / E_LCDM_with_rad(z, Om, h, Neff=Neff),
                  0.0, Z_STAR, limit=300, epsabs=1e-9, epsrel=1e-9)
    return (C_KMS / H0) * val


def compute_R_lA(H0, Om, obh2, Neff=3.046):
    """Compute (R, l_A, omega_b) for given LCDM parameters.

    r_s(z*) is approximated as r_d * 1.0184  (Planck 2018 ratio, Aubourg+15).
    """
    h = H0 / 100.0
    omh2 = Om * h ** 2
    DM_zstar = comoving_distance_to_zstar(H0, Om, h=h, Neff=Neff)
    rd = rd_planck_anchored(omh2, obh2, Neff=Neff)
    rs_zstar = rd / 1.0184
    R = np.sqrt(omh2) * 100 * DM_zstar / C_KMS  # = sqrt(omh2) * H0/c * DM, but c in km/s, H0 in km/s/Mpc
    lA = np.pi * DM_zstar / rs_zstar
    return R, lA, obh2


def chi2_cmb(H0, Om, obh2, Neff=3.046):
    R, lA, _ = compute_R_lA(H0, Om, obh2, Neff=Neff)
    d = np.array([R, lA, obh2]) - PLANCK_MEAN
    return float(d @ PLANCK_INVCOV @ d)


def chi2_cmb_with_rd(H0, Om, obh2, rd, Neff=3.046):
    """Variant where r_d is supplied as a free parameter (not derived from
    obh2/omh2). Useful to test 'what if pre-recombination physics shifts r_d
    independently of standard LCDM'."""
    h = H0 / 100.0
    omh2 = Om * h ** 2
    DM_zstar = comoving_distance_to_zstar(H0, Om, h=h, Neff=Neff)
    rs_zstar = rd / 1.0184  # ratio still assumed
    R = np.sqrt(omh2) * 100 * DM_zstar / C_KMS
    lA = np.pi * DM_zstar / rs_zstar
    d = np.array([R, lA, obh2]) - PLANCK_MEAN
    return float(d @ PLANCK_INVCOV @ d)
