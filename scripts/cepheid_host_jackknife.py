"""Cepheid host drop-one jackknife on Pantheon+SH0ES H0.

Tests robustness: if dropping a single Cepheid host shifts H0 by more
than the formal sigma, that host dominates the constraint - a red flag
for systematic dependence.

We group the 77 SH0ES calibrator SNe by their CEPH_DIST value (shared
within a host). Each unique CEPH_DIST is one host.

For each host, we (a) drop all calibrators in that host, (b) refit
(H0, Om, M) jointly against the full Pantheon+ likelihood, (c) record
H0_drop - H0_full.
"""
import sys
from pathlib import Path
import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.fast_lh import HubbleLH
from src.data_loaders import load_pantheon_plus
from src.cosmo import distance_modulus_grid


def main():
    lh = HubbleLH()
    pp = load_pantheon_plus()
    arr = pp["data"]
    calib_mask = pp["is_calibrator"]
    ceph_dist = arr["CEPH_DIST"]

    # Cluster by ceph_dist value (same host = same value)
    calib_dists = ceph_dist[calib_mask]
    unique_dists, counts = np.unique(calib_dists, return_counts=True)
    print(f"Found {len(unique_dists)} distinct Cepheid hosts in 77 calibrators")
    print(f"  ceph_dist range: [{unique_dists.min():.3f}, {unique_dists.max():.3f}]")

    # Establish the 'host index' mapping for each calibrator
    # Two calibrators are in the same host iff they share ceph_dist
    host_idx_full = np.full(len(arr), -1, dtype=int)
    for h, d in enumerate(unique_dists):
        host_idx_full[calib_mask & (ceph_dist == d)] = h
    print(f"  Calibrators per host (count distribution): "
          f"min={counts.min()}, max={counts.max()}, median={int(np.median(counts))}")

    # Pre-compute Cinv for the full sample and for each drop-host scenario
    print("\nPre-computing inverse covariance for full + 37 drop-host masks...")
    Cinv_full = lh.Cinv_pp
    cache = {-1: (np.ones(len(arr), dtype=bool), Cinv_full,
                  lh.zHEL, lh.zHD, lh.mb, lh.is_calib, lh.ceph_dist)}
    for h in range(len(unique_dists)):
        keep = ~(calib_mask & (host_idx_full == h))
        Cinv_sub = np.linalg.inv(lh.cov_pp[np.ix_(keep, keep)])
        cache[h] = (keep, Cinv_sub,
                    lh.zHEL[keep], lh.zHD[keep], lh.mb[keep],
                    lh.is_calib[keep], lh.ceph_dist[keep])
    print("  done.\n")

    def chi2(theta, drop_host=-1):
        H0, Om, M = theta
        if not (40 < H0 < 100 and 0.05 < Om < 0.7 and -22 < M < -18):
            return 1e30
        keep, Cinv_sub, zHEL, zHD, mb, ic, cd = cache[drop_host]
        mu_th = distance_modulus_grid(zHEL, zHD, H0, Om)
        mu_model = np.where(ic, cd, mu_th)
        resid = mb - (M + mu_model)
        return float(resid @ Cinv_sub @ resid)

    print("\n=== Reference (all calibrators) ===")
    r = minimize(lambda t: chi2(t, -1), [73, 0.30, -19.25], method="Nelder-Mead",
                 options=dict(xatol=1e-4, fatol=1e-3, maxiter=2000))
    H0_ref, Om_ref, M_ref = r.x
    print(f"  H0 = {H0_ref:.3f}, Om = {Om_ref:.4f}, M = {M_ref:.4f}, chi2 = {r.fun:.2f}")

    # 1-sigma error estimate via quick profile
    # Δχ²=1 width from this same data
    def prof(H0_t):
        rp = minimize(lambda x: chi2([H0_t, x[0], x[1]], -1),
                      [Om_ref, M_ref], method="Nelder-Mead",
                      options=dict(xatol=1e-4, fatol=1e-3))
        return rp.fun
    from scipy.optimize import brentq
    try:
        sig = brentq(lambda dH: prof(H0_ref + dH) - r.fun - 1.0, 0.01, 5.0)
    except Exception:
        sig = float('nan')
    print(f"  Profile 1σ ≈ ±{sig:.3f}")

    # Drop-one
    print("\n=== Drop-one Cepheid host ===")
    print(f"  {'host_idx':>8s} {'ceph_d':>8s} {'n_drop':>7s} {'H0_drop':>9s} {'dH0':>9s} {'(sig)':>6s}")
    deltas = []
    for h, d in enumerate(unique_dists):
        n_in_host = (calib_mask & (host_idx_full == h)).sum()
        rd_h = minimize(lambda t: chi2(t, h), [H0_ref, Om_ref, M_ref],
                        method="Nelder-Mead",
                        options=dict(xatol=1e-4, fatol=1e-3, maxiter=2000))
        H0_d = rd_h.x[0]
        dH = H0_d - H0_ref
        in_sig = dH / sig if sig > 0 else 0
        deltas.append((d, h, n_in_host, H0_d, dH, in_sig))
        print(f"  {h:>8d} {d:>8.3f} {n_in_host:>7d} {H0_d:>9.3f} {dH:>+8.3f} {in_sig:>+6.2f}")

    # Sort by largest |ΔH0|
    deltas.sort(key=lambda x: -abs(x[4]))
    print("\n=== Top 5 hosts by |ΔH0| (largest impact) ===")
    print(f"  {'rank':>4s} {'ceph_d':>8s} {'n':>4s} {'dH0':>9s} {'(sig_ref)':>10s}")
    for rank, (d, h, n, H0d, dH, ins) in enumerate(deltas[:5], 1):
        print(f"  {rank:>4d} {d:>8.3f} {n:>4d} {dH:>+8.3f} {ins:>+8.2f}")

    print("\n=== Summary statistics ===")
    arr_dH = np.array([t[4] for t in deltas])
    print(f"  mean ΔH0      = {arr_dH.mean():+.3f}")
    print(f"  std  ΔH0      = {arr_dH.std():.3f}")
    print(f"  max |ΔH0|     = {np.abs(arr_dH).max():.3f}  ({np.abs(arr_dH).max()/sig:+.2f} σ_ref)")
    if np.abs(arr_dH).max() > 2 * sig:
        print(f"  WARNING: at least one host shifts H0 by > 2σ — calibration NOT robust")
    elif np.abs(arr_dH).max() > sig:
        print(f"  Notable: at least one host shifts H0 by > 1σ — borderline robust")
    else:
        print(f"  GOOD: no single host shifts H0 by more than 1σ — calibration is robust")


if __name__ == "__main__":
    main()
