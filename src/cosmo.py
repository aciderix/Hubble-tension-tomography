"""Pure-numpy cosmological distance utilities.

No dependency on CAMB/CLASS so it stays trivially auditable. Use this for
baselines and quick falsification tests, then cross-check with CAMB in
src/cosmo_camb.py.
"""
import numpy as np
from scipy import integrate

C_KMS = 299792.458  # km/s


def E_LCDM(z, Om):
    return np.sqrt(Om * (1.0 + z) ** 3 + (1.0 - Om))


def H_of_z(z, H0, Om):
    return H0 * E_LCDM(z, Om)


def comoving_distance(z, H0, Om):
    z = np.atleast_1d(z)
    out = np.empty_like(z, dtype=float)
    for i, zi in enumerate(z):
        v, _ = integrate.quad(lambda zp: 1.0 / E_LCDM(zp, Om), 0.0, zi)
        out[i] = v
    return (C_KMS / H0) * out


def luminosity_distance(zhel, zcmb, H0, Om):
    dC = comoving_distance(zcmb, H0, Om)
    return (1.0 + np.atleast_1d(zhel)) * dC


def distance_modulus(zhel, zcmb, H0, Om):
    dL = luminosity_distance(zhel, zcmb, H0, Om)
    return 5.0 * np.log10(dL) + 25.0


def angular_diameter_distance(z, H0, Om):
    return comoving_distance(z, H0, Om) / (1.0 + np.atleast_1d(z))


def DV(z, H0, Om):
    """BAO spherically-averaged distance."""
    z = np.atleast_1d(z)
    dM = comoving_distance(z, H0, Om)
    Hz = H_of_z(z, H0, Om)
    return (z * dM ** 2 * C_KMS / Hz) ** (1.0 / 3.0)
