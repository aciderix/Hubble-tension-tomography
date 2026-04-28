"""Smoke-test ΛCDM fit to Pantheon+SH0ES.

Goal: anchor the baseline. Reproduce H0 ≈ 73 km/s/Mpc when including
Cepheid calibrators (SH0ES anchor), and verify that without the calibrators
the SN magnitudes are degenerate with H0 (i.e. only Om is constrained).

Following Brout et al. 2022 / Riess et al. 2022 likelihood structure.
"""
import sys
from pathlib import Path
import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.data_loaders import load_pantheon_plus
from src.cosmo import distance_modulus


def chi2(theta, dat, use_calibrators=True):
    H0, Om, M = theta
    if H0 <= 20 or H0 >= 200 or Om <= 0 or Om >= 1.0:
        return 1e30
    mu_th = distance_modulus(dat["zHEL"], dat["zHD"], H0, Om)
    # For calibrators, the model distance is replaced by CEPH_DIST
    mu_model = np.where(dat["is_calibrator"], dat["ceph_dist"], mu_th)
    resid = dat["mb"] - (M + mu_model)
    if not use_calibrators:
        mask = ~dat["is_calibrator"]
        resid = resid[mask]
        cov = dat["cov"][np.ix_(mask, mask)]
    else:
        cov = dat["cov"]
    L = np.linalg.cholesky(cov)
    y = np.linalg.solve(L, resid)
    return float(y @ y)


def main():
    dat = load_pantheon_plus()
    print(f"loaded {dat['n']} SNe, {dat['is_calibrator'].sum()} calibrators")

    print("\n=== Fit 1: SN-only (no calibrators) — H0 should be unconstrained ===")
    # Fix M=-19.25 by hand (degenerate); fit Om at H0=70
    for H0 in [60, 70, 80]:
        f = lambda x: chi2((H0, x[0], x[1]), dat, use_calibrators=False)
        r = minimize(f, x0=[0.3, -19.25], method="Nelder-Mead",
                     options=dict(xatol=1e-4, fatol=1e-3))
        print(f"  H0={H0}: Om*={r.x[0]:.4f}, M*={r.x[1]:.4f}, chi2={r.fun:.2f}")

    print("\n=== Fit 2: SN + SH0ES Cepheid calibrators ===")
    f = lambda x: chi2(x, dat, use_calibrators=True)
    r = minimize(f, x0=[73.0, 0.3, -19.25], method="Nelder-Mead",
                 options=dict(xatol=1e-4, fatol=1e-3, maxiter=2000))
    H0_best, Om_best, M_best = r.x
    print(f"  H0  = {H0_best:.3f}")
    print(f"  Om  = {Om_best:.4f}")
    print(f"  M   = {M_best:.4f}")
    print(f"  chi2 = {r.fun:.2f}, dof ≈ {dat['n'] - 3}")

    # Quick uncertainty: profile-likelihood on H0 with Δχ² = 1
    print("\n=== Profile chi2(H0) ===")
    for H0 in np.arange(70.0, 76.1, 0.5):
        g = lambda x: chi2((H0, x[0], x[1]), dat)
        rg = minimize(g, x0=[Om_best, M_best], method="Nelder-Mead",
                      options=dict(xatol=1e-4, fatol=1e-3))
        print(f"  H0={H0:.1f}: Δχ² = {rg.fun - r.fun:+.3f}")


if __name__ == "__main__":
    main()
