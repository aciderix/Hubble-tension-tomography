"""ΛCDM joint fit with rd derived from BBN-priored Ωb h² + matter density.

This is the 'inverse distance ladder using only BBN + BAO + SN' result.
Independent of CMB. The inferred H0 should match the published value (~68.5).

Then test whether SH0ES anchor is statistically consistent.
"""
import sys
from pathlib import Path
import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.fast_lh import HubbleLH
from src.cosmo import rd_aubourg


def main():
    lh = HubbleLH()

    # Cooke+ 2018 BBN: Omega_b h^2 = 0.02218 +/- 0.00055
    obh2_BBN = 0.02218
    sig_obh2 = 0.00055

    # ============================================================
    # G) LCDM + BBN-prior on Ωb h², free Ωb h² with Gaussian prior
    # ============================================================
    print("\n" + "=" * 70)
    print("G) LCDM joint fit, rd = rd_Aubourg(Ωmh², Ωbh²) with BBN prior on Ωbh².")
    print("=" * 70)
    def f(x):
        H0, Om, M, obh2 = x
        if not (50 < H0 < 90 and 0.1 < Om < 0.6 and -22 < M < -18
                and 0.018 < obh2 < 0.026):
            return 1e30
        h = H0 / 100
        rd = rd_aubourg(Om * h**2, obh2)
        chi2 = lh.chi2_pp(H0, Om, M) + lh.chi2_bao(H0, Om, rd)
        chi2 += ((obh2 - obh2_BBN) / sig_obh2) ** 2
        return chi2

    rG = minimize(f, [73, 0.31, -19.25, obh2_BBN], method="Nelder-Mead",
                  options=dict(xatol=1e-5, fatol=1e-4, maxiter=5000))
    H0_G, Om_G, M_G, obh2_G = rG.x
    h_G = H0_G / 100
    rd_G = rd_aubourg(Om_G * h_G**2, obh2_G)
    chi2_pp_G = lh.chi2_pp(H0_G, Om_G, M_G)
    chi2_bao_G = lh.chi2_bao(H0_G, Om_G, rd_G)
    chi2_bbn = ((obh2_G - obh2_BBN) / sig_obh2) ** 2
    print(f"  H0={H0_G:.3f}  Om={Om_G:.4f}  M={M_G:.4f}  Ωb h²={obh2_G:.5f}")
    print(f"  rd_inferred = {rd_G:.2f} Mpc")
    print(f"  chi2: PP+SH0ES={chi2_pp_G:.2f}  DESI={chi2_bao_G:.2f}  "
          f"BBN={chi2_bbn:.3f}  total={rG.fun:.2f}")

    # ============================================================
    # H) Same as G but DROP the SH0ES Cepheid anchor (Pantheon+ HF only)
    # ============================================================
    print("\n" + "=" * 70)
    print("H) Same with NO SH0ES anchor (Pantheon+ Hubble flow only) + DESI + BBN.")
    print("=" * 70)
    def f_H(x):
        H0, Om, M, obh2 = x
        if not (50 < H0 < 90 and 0.1 < Om < 0.6 and -22 < M < -18
                and 0.018 < obh2 < 0.026):
            return 1e30
        h = H0 / 100
        rd = rd_aubourg(Om * h**2, obh2)
        chi2 = lh.chi2_pp(H0, Om, M, use_calib=False) + lh.chi2_bao(H0, Om, rd)
        chi2 += ((obh2 - obh2_BBN) / sig_obh2) ** 2
        return chi2

    rH = minimize(f_H, [70, 0.31, -19.4, obh2_BBN], method="Nelder-Mead",
                  options=dict(xatol=1e-5, fatol=1e-4, maxiter=5000))
    H0_H, Om_H, M_H, obh2_H = rH.x
    h_H = H0_H / 100
    rd_H = rd_aubourg(Om_H * h_H**2, obh2_H)
    print(f"  H0={H0_H:.3f}  Om={Om_H:.4f}  M={M_H:.4f}  Ωb h²={obh2_H:.5f}")
    print(f"  rd = {rd_H:.2f} Mpc")
    print(f"  total chi2 = {rH.fun:.2f}")
    print("  -> THE 'CMB-INDEPENDENT' early-universe-side answer")

    # Tension between (G) and (H) is the actual data tension
    print(f"\n  Δχ² between (G) including SH0ES vs (H) excluding SH0ES "
          f"= {rG.fun - rH.fun:+.2f}")
    print(f"  H0 shift: {H0_G - H0_H:+.3f}  ({100*(H0_G - H0_H)/H0_H:+.2f}%)")

    # ============================================================
    # I) Internal consistency — reduced chi2 per dataset at (G)
    # ============================================================
    print("\n" + "=" * 70)
    print("I) Internal consistency at G")
    print("=" * 70)
    n_pp = lh.n_pp; n_bao = lh.n_bao
    print(f"  PP+SH0ES: χ²/dof = {chi2_pp_G:.1f}/{n_pp-3} = "
          f"{chi2_pp_G/(n_pp-3):.3f}")
    print(f"  DESI:     χ²/dof = {chi2_bao_G:.1f}/{n_bao-3} = "
          f"{chi2_bao_G/(n_bao-3):.3f}  -> p ~"
          f" {100*(1-_chi2_cdf(chi2_bao_G, n_bao-3)):.1f}% (low if reduced<<1)")


def _chi2_cdf(x, k):
    from scipy.stats import chi2
    return chi2.cdf(x, k)


if __name__ == "__main__":
    main()
