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


def DH(z, H0, Om):
    """BAO Hubble distance c / H(z)."""
    return C_KMS / H_of_z(z, H0, Om)


def rd_aubourg(omh2, obh2, omnuh2=0.0):
    """Aubourg+ 2015 fitting formula for the sound horizon at the drag epoch.

    Accurate to ~0.1% across the Planck parameter region.
    """
    return (55.154 * np.exp(-72.3 * (omnuh2 + 0.0006) ** 2)
            / (omh2 ** 0.25351 * obh2 ** 0.12807))


def E_w0wa(z, Om, w0, wa):
    """w0waCDM expansion: w(z) = w0 + wa * z/(1+z) (CPL)."""
    z = np.atleast_1d(z)
    a = 1.0 / (1.0 + z)
    de = (1.0 - Om) * a ** (-3.0 * (1.0 + w0 + wa)) * np.exp(-3.0 * wa * (1.0 - a))
    return np.sqrt(Om * (1.0 + z) ** 3 + de)


def comoving_distance_w0wa(z, H0, Om, w0, wa):
    z = np.atleast_1d(z)
    out = np.empty_like(z, dtype=float)
    for i, zi in enumerate(z):
        v, _ = integrate.quad(
            lambda zp: 1.0 / E_w0wa(np.array([zp]), Om, w0, wa)[0], 0.0, zi)
        out[i] = v
    return (C_KMS / H0) * out


def comoving_distance_w0wa_grid(z_eval, H0, Om, w0, wa, n_grid=4000):
    """Vectorized comoving distance under CPL w0waCDM."""
    z_eval = np.atleast_1d(z_eval).astype(float)
    z_max = max(float(z_eval.max()), 0.001) * 1.01
    zg = np.linspace(0.0, z_max, n_grid)
    invE = 1.0 / E_w0wa(zg, Om, w0, wa)
    dz = zg[1] - zg[0]
    cum = np.concatenate([[0.0], np.cumsum(0.5 * (invE[:-1] + invE[1:]) * dz)])
    return (C_KMS / H0) * np.interp(z_eval, zg, cum)


def comoving_distance_grid(z_eval, H0, Om, n_grid=4000):
    """Vectorized D_C using a precomputed cumtrapz table on a log-z grid.

    Evaluates D_C at arbitrary z_eval in O(N) instead of N quad calls.
    Accuracy: ~1e-6 over z in [1e-4, 3] for n_grid=4000.
    """
    z_eval = np.atleast_1d(z_eval).astype(float)
    z_max = max(float(z_eval.max()), 0.001) * 1.01
    # Use a fine linear grid; could go log but linear is fine here
    zg = np.linspace(0.0, z_max, n_grid)
    invE = 1.0 / E_LCDM(zg, Om)
    # cumulative integral via trapezoid
    dz = zg[1] - zg[0]
    cum = np.concatenate([[0.0], np.cumsum(0.5 * (invE[:-1] + invE[1:]) * dz)])
    return (C_KMS / H0) * np.interp(z_eval, zg, cum)


def distance_modulus_grid(zhel, zcmb, H0, Om, n_grid=4000):
    """Vectorized mu = 5 log10(dL) + 25 with dL = (1+zhel) * D_C(zcmb)."""
    dC = comoving_distance_grid(zcmb, H0, Om, n_grid=n_grid)
    dL = (1.0 + np.atleast_1d(zhel)) * dC
    return 5.0 * np.log10(dL) + 25.0
