"""2D scan: fix rd, scan ΔM_calib (Cepheid systematic offset).

Tests whether the fit at (rd, ΔM_calib) is degenerate or whether the data
prefers a particular split between early-time (rd) and local (ΔM) shifts.
"""
import sys
from pathlib import Path
import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.fast_lh import HubbleLH


def main():
    lh = HubbleLH()

    rd_grid = np.linspace(132.0, 150.0, 13)
    dM_grid = np.linspace(-0.05, 0.20, 11)
    chi2_grid = np.zeros((len(rd_grid), len(dM_grid)))
    H0_grid = np.zeros_like(chi2_grid)
    Om_grid = np.zeros_like(chi2_grid)

    for i, rd in enumerate(rd_grid):
        for j, dM in enumerate(dM_grid):
            # M is the Hubble-flow absolute mag; M_calib = M + dM
            def f(x):
                H0, Om, M = x
                if not (50 < H0 < 90 and 0.1 < Om < 0.6 and -22 < M < -18):
                    return 1e30
                return (lh.chi2_pp(H0, Om, M, M_calib=M + dM)
                        + lh.chi2_bao(H0, Om, rd))
            r = minimize(f, [70, 0.31, -19.4], method="Nelder-Mead",
                         options=dict(xatol=1e-4, fatol=1e-3, maxiter=2000))
            chi2_grid[i, j] = r.fun
            H0_grid[i, j] = r.x[0]
            Om_grid[i, j] = r.x[1]

    chi2_min = chi2_grid.min()
    print(f"Joint min chi2 = {chi2_min:.2f}")
    i0, j0 = np.unravel_index(chi2_grid.argmin(), chi2_grid.shape)
    print(f"  at rd={rd_grid[i0]:.1f}, ΔM={dM_grid[j0]:+.3f}, "
          f"H0={H0_grid[i0,j0]:.2f}, Om={Om_grid[i0,j0]:.4f}")

    print("\nΔχ² grid (rows=rd, cols=ΔM_calib):")
    print("ΔM:   " + "  ".join(f"{m:+.3f}" for m in dM_grid))
    for i, rd in enumerate(rd_grid):
        row = "  ".join(f"{chi2_grid[i,j]-chi2_min:5.1f}" for j in range(len(dM_grid)))
        print(f"rd={rd:5.1f}  {row}")

    print("\nH0 grid (rows=rd, cols=ΔM):")
    print("ΔM:   " + "  ".join(f"{m:+.3f}" for m in dM_grid))
    for i, rd in enumerate(rd_grid):
        row = "  ".join(f"{H0_grid[i,j]:5.2f}" for j in range(len(dM_grid)))
        print(f"rd={rd:5.1f}  {row}")

    print("\n# rd=147 (Planck), ΔM=0:  forced consistency baseline")
    print(f"  Δχ² there = "
          f"{chi2_grid[np.argmin(np.abs(rd_grid-147)),np.argmin(np.abs(dM_grid-0))]-chi2_min:.2f}")


if __name__ == "__main__":
    main()
