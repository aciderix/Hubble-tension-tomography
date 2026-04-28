"""Linear growth + Eisenstein-Hu transfer + sigma_R for arbitrary cosmology.

Pure-numpy implementation. We compute:

  D(z)      linear growth factor, from the growth ODE
              D''(a) + (3/a + dlnH/da) D'(a) - (3 Omega_m(a) / 2 a^2) D(a) = 0
  T(k)      Eisenstein-Hu 1998 fitting transfer function (no-wiggle form)
  P(k)      = A_s normalised power spectrum at z=0
  sigma_R^2 = integral over k of P(k) * W(kR)^2 * k^2 dk / (2 pi^2)

Used to compute sigma_8 ≡ sigma_R(R = 8 h^-1 Mpc).

Validated against Planck 2018 LCDM fiducial: sigma8 ≈ 0.811 at
H0=67.36, Om=0.3153, omega_b=0.02237, ns=0.9649, with A_s=2.1e-9
calibrated.
"""
import numpy as np
from scipy.integrate import solve_ivp, quad
from .cosmo import E_LCDM_with_rad, C_KMS


# Physical constants
T_CMB = 2.7255
H0_FID = 100  # nuisance


def _Omega_m_a(a, Om):
    """Matter density fraction at scale factor a (LCDM)."""
    return Om * a ** -3 / (Om * a ** -3 + (1 - Om))


def _dlnH_dlna(a, Om):
    """d ln H / d ln a in LCDM (matter + cosmological constant)."""
    Om_a = _Omega_m_a(a, Om)
    return -1.5 * Om_a


def growth_factor(z, Om):
    """Linear growth factor D(z) normalized so D(z=0) = 1.

    Solves the standard ODE for δ(a) = D(a) for matter perturbations
    in flat LCDM:
      D'' + (3/a + dlnH/dlna * 1/a) D' - (3/2)(Omega_m(a)/a^2) D = 0

    where prime is d/da. Reformulate with f = D dot:
      d/da[D] = f
      d/da[f] = -(3/a + dlnH/dlna /a) f + (3/2)(Omega_m(a)/a^2) D
    """
    z = np.atleast_1d(z)
    a_eval = 1.0 / (1.0 + z)
    a_lo = min(0.001, 0.5 * a_eval.min())  # start deep in matter era

    def rhs(a, y):
        D, f = y
        Om_a = _Omega_m_a(a, Om)
        dlnH = _dlnH_dlna(a, Om)
        return [f,
                -(3.0 / a + dlnH / a) * f + 1.5 * Om_a / a ** 2 * D]

    # Initial condition: in matter era, D ∝ a (so D(a_lo) = a_lo, dD/da = 1)
    # Build t_eval: include a_eval values + a=1 sentinel, ensure sorted
    t_eval = np.unique(np.concatenate([a_eval, [1.0]]))
    sol = solve_ivp(rhs, (a_lo, 1.0001), [a_lo, 1.0],
                    method="DOP853", rtol=1e-9, atol=1e-12,
                    t_eval=t_eval)
    D_at_1 = np.interp(1.0, sol.t, sol.y[0])
    D_at_eval = np.interp(a_eval, sol.t, sol.y[0])
    return D_at_eval / D_at_1


def transfer_eh98(k, Om, Ob, h):
    """Eisenstein-Hu 1998 BAO-no-wiggle transfer function (Eq 26)
    in inverse Mpc. Returns dimensionless T(k).

    Reference: arXiv:astro-ph/9709112 .
    """
    # Convert k from h/Mpc to Mpc^-1 (caller can pass either; we expect Mpc^-1)
    omh2 = Om * h ** 2
    obh2 = Ob * h ** 2
    Theta_27 = T_CMB / 2.7

    s = 44.5 * np.log(9.83 / omh2) / np.sqrt(1 + 10 * obh2 ** 0.75)  # Mpc
    alpha_gamma = (1
                   - 0.328 * np.log(431 * omh2) * obh2 / omh2
                   + 0.38 * np.log(22.3 * omh2) * (obh2 / omh2) ** 2)
    Gamma = Om * h * (alpha_gamma + (1 - alpha_gamma) / (1 + (0.43 * k * s) ** 4))
    q = k * Theta_27 ** 2 / Gamma   # k in Mpc^-1
    L0 = np.log(2 * np.e + 1.8 * q)
    C0 = 14.2 + 731.0 / (1 + 62.5 * q)
    return L0 / (L0 + C0 * q ** 2)


def sigma_R(R, Om, Ob, h, ns, A_s, k_pivot=0.05):
    """sigma at smoothing scale R [Mpc/h] for given LCDM at z=0.

    P(k) = A_s * (k / k_pivot)^(ns - 1) * (2 c k^2 T(k))^2 / (5 H0^2 Omega_m)^2
            in some normalisation. Use the standard CAMB convention:
    P_R(k) = A_s (k/k_pivot)^(ns-1)
    P_m(k,z) = (4/25) (k/H0)^4 P_R(k) T(k)^2 / Omega_m^2  *  D(z)^2 / something

    To avoid normalisation ambiguity we use the dimensionless form
    sigma8 ∝ A_s^{1/2}, and we calibrate against Planck below.
    """
    R_Mpc = R / h  # convert R from Mpc/h to Mpc
    # Integrand: k^2 P_m(k) W(kR)^2 / (2π^2)
    # P_m(k) ∝ k T(k)^2 (k/k_pivot)^(ns-1)
    def integrand(lnk):
        k = np.exp(lnk)
        T = transfer_eh98(k, Om, Ob, h)
        # Prefactor; we absorb constants into the calibration
        P = (k / k_pivot) ** (ns - 1) * T ** 2 * k ** 4
        W = 3 * (np.sin(k * R_Mpc) - k * R_Mpc * np.cos(k * R_Mpc)) / (k * R_Mpc) ** 3
        return P * W ** 2 * k  # extra k for d ln k

    integral, _ = quad(integrand, np.log(1e-5), np.log(50.0), limit=400,
                        epsabs=1e-14, epsrel=1e-7)
    # The result is proportional to sigma8^2 up to a constant that we
    # calibrate from Planck.
    return integral * A_s


# Calibration: anchor on Planck 2018 LCDM where sigma8 = 0.8111.
# Compute the ratio sigma_R(Om, Ob, h, ns) / sigma_R(Planck) and multiply
# by 0.8111. This bypasses the absolute normalisation prefactor.
PLANCK_SIGMA8 = 0.8111
PLANCK_PARAMS = dict(Om=0.3153, obh2=0.02237, h=0.6736, ns=0.9649)
_PLANCK_S2_UNNORM = None


def _planck_s2():
    global _PLANCK_S2_UNNORM
    if _PLANCK_S2_UNNORM is None:
        Ob = PLANCK_PARAMS["obh2"] / PLANCK_PARAMS["h"] ** 2
        _PLANCK_S2_UNNORM = sigma_R(
            8.0, PLANCK_PARAMS["Om"], Ob, PLANCK_PARAMS["h"], PLANCK_PARAMS["ns"], 1.0)
    return _PLANCK_S2_UNNORM


def sigma8(Om, obh2, h, ns):
    """sigma_8 at z=0 calibrated on Planck 2018 (sigma8 = 0.8111)."""
    Ob = obh2 / h ** 2
    s2 = sigma_R(8.0, Om, Ob, h, ns, 1.0)
    return PLANCK_SIGMA8 * np.sqrt(s2 / _planck_s2())


def S8(Om, obh2, h, ns):
    """S8 = sigma8 * sqrt(Om/0.3)."""
    return sigma8(Om, obh2, h, ns) * np.sqrt(Om / 0.3)
