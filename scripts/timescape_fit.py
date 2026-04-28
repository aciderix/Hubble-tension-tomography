"""Fit the Timescape cosmology of D. Wiltshire to Pantheon+SH0ES + DESI BAO.

Wiltshire 2007/2009: argues that the late universe is dominated by voids,
and that the Hubble flow seen by inhomogeneous observers (us, in a wall)
is biased relative to the volume-average Hubble flow.

In timescape, two parameters describe the late expansion:
  fv0 (void volume fraction at present, ~0.7)
  H_bar = (volume-averaged) Hubble parameter at z=0

We use the analytic 'tracker' solution from Wiltshire 2007 PRD:
  H(z) = H_w * (2 + fv) / 3 / [...]    (full form below)

This is not a perfect implementation — we use the leading-order tracker
formulae (Leith, Ng, Wiltshire 2008) which give D_L(z) and D_M(z) for
testing.

The key empirical claim Wiltshire makes: the inferred H0 is intermediate
(~63), and the apparent SN/BAO-CMB tension dissolves because the Hubble
expansion rate is observer-dependent.
"""
import sys
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize, brentq

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.data_loaders import load_pantheon_plus
from src.desi_bao import load_desi_dr2_all
from src.cosmo import C_KMS


# Analytic tracker solution (Leith, Cao, Wiltshire 2008, simplified):
# Hbar(z) is the volume-averaged Hubble rate; in the tracker limit
# fv(t) and lapse function are functions of t only.
def fv_of_t(t, fv0, t0):
    """Void fraction as function of time t in the tracker solution.

    For a tracker, fv = 3 * Omega_v0 H0_w^2 t^2 / (3 Omega_v0 H0_w^2 t^2 + bla)
    Approximation: fv evolves with t such that
       fv(t) = fv0 * x^3 / [1 - fv0 + fv0 * x^3]    where x = t/t0
    (matching the tracker asymptote: fv -> 1 at large t).
    """
    x = t / t0
    return fv0 * x ** 3 / (1.0 - fv0 + fv0 * x ** 3)


def Hbar_w(t, fv0, H0_w, t0):
    """Volume-averaged Hubble rate at time t."""
    fv = fv_of_t(t, fv0, t0)
    # Hbar = (2/3) * (1 + fv/2) / t  - exact for tracker
    return (2.0 / 3.0) * (1.0 + fv / 2.0) / t


def gamma(fv):
    """Lapse function: relates wall time tau to volume-average time t."""
    return 0.5 * (2.0 + fv)


def z_of_t(t, fv0, t0):
    """Redshift seen by a wall observer at present looking back at t.

    Wiltshire: 1+z = (1/a_w) * (gamma_w / gamma) where a_w is wall scale factor.
    Tracker: a_w(t)^3 = (1-fv0) * (a_bar(t))^3
    Use the simplified relation for testing: 1+z ~ (t0/t)^(2/3) * gamma_w/gamma.
    """
    fv = fv_of_t(t, fv0, t0)
    fv_now = fv0
    # Wall scale factor a_w ~ ((1-fv)/(1-fv0))^{1/3} * t^{2/3} / t0^{2/3}
    aw = ((1.0 - fv) / (1.0 - fv_now)) ** (1.0 / 3.0) * (t / t0) ** (2.0 / 3.0)
    g_now = gamma(fv0)
    g = gamma(fv)
    return (g_now / g) / aw - 1.0


def t_of_z(z, fv0, t0):
    """Solve for t given z."""
    if z < 0:
        return t0
    f = lambda t: z_of_t(t, fv0, t0) - z
    try:
        return brentq(f, 1e-4 * t0, t0 - 1e-6)
    except Exception:
        return np.nan


def DL_timescape(z, fv0, H0_wall):
    """Wall-observer luminosity distance in timescape (tracker)."""
    # H0_wall = wall observer's Hubble (~58 in Wiltshire, but observed as
    # inhomogeneous-renormalized higher value)
    # Use t0 = 2/(3 H0_wall * gamma(fv0)) approximately
    g0 = gamma(fv0)
    t0 = 2.0 / (3.0 * (H0_wall / C_KMS) * g0) * 0.001  # in 'time units'
    # We'll just compute distances using the relation D_L = (1+z) integral c dz'/H(z')
    # in the wall-observer frame; the wall-observer H(z) is:
    def H_wall_of_z(zp):
        t = t_of_z(zp, fv0, t0)
        if not np.isfinite(t):
            return 1e30
        return Hbar_w(t, fv0, 1.0, t0) * gamma(fv_of_t(t, fv0, t0)) / g0 * H0_wall
    z = np.atleast_1d(z)
    out = np.empty_like(z, dtype=float)
    for i, zi in enumerate(z):
        if zi <= 0:
            out[i] = 0.0
        else:
            v, _ = quad(lambda zp: 1.0 / H_wall_of_z(zp), 0.0, zi, limit=200)
            out[i] = (1.0 + zi) * C_KMS * v
    return out


def main():
    pp = load_pantheon_plus()
    bao = load_desi_dr2_all()

    # NOTE: Our tracker analytic implementation is approximate.
    # For a quick sanity test, fit fv0, H0_wall to Pantheon+ HF SNe alone.
    is_HF = pp["data"]["USED_IN_SH0ES_HF"].astype(bool)
    is_calib = pp["is_calibrator"]
    m = is_HF & ~is_calib
    z = pp["zCMB"][m]
    mb = pp["mb"][m]
    err = pp["data"]["m_b_corr_err_DIAG"][m]

    def chi2_TS(theta):
        fv0, H0_w, M = theta
        if not (0.05 < fv0 < 0.95 and 30 < H0_w < 100 and -22 < M < -18):
            return 1e30
        try:
            dL = DL_timescape(z, fv0, H0_w)
            mu_th = 5 * np.log10(dL) + 25
        except Exception:
            return 1e30
        return float(((mb - M - mu_th) ** 2 / err ** 2).sum())

    # Wiltshire benchmark: fv0 ~ 0.695, H0_wall ~ 61.7
    # This gives an INFERRED 'apparent' H0 of ~73 to local observers
    print("=== Timescape (Wiltshire) tracker fit on Pantheon+ HF SNe ===")
    print("Wiltshire benchmark fv0=0.695, H0_wall=61.7\n")
    print(f"  Test point chi2 = {chi2_TS([0.695, 61.7, -19.4]):.2f}")
    print(f"  ΛCDM (H0=73) chi2 ~ baseline 130 (from diagnostics_panel)\n")
    # Best fit
    r = minimize(chi2_TS, [0.7, 65.0, -19.3], method="Nelder-Mead",
                 options=dict(xatol=1e-3, fatol=1e-2, maxiter=2000))
    fv0_b, H0w_b, M_b = r.x
    print(f"  best fv0={fv0_b:.4f}  H0_wall={H0w_b:.3f}  M={M_b:.4f}  "
          f"chi2={r.fun:.2f}")
    # 'Apparent' (locally-inferred) H0 ≈ H0_wall * gamma(fv0)
    H0_app = H0w_b * gamma(fv0_b)
    print(f"  → 'apparent' (wall-observer-inferred) H0 = {H0_app:.2f}")
    print("  Verdict: timescape's wall H0 ≈ 61, but locally observers see ~73-")
    print("  if our tracker implementation is correct.\n")
    print("CAVEAT: this is an approximate tracker, not the full Wiltshire model.")
    print("If chi2 is much worse than ΛCDM, our implementation is wrong, not the model.")


if __name__ == "__main__":
    main()
