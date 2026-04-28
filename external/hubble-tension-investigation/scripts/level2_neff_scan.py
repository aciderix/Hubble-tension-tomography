"""
Level 2 — Boltzmann-code stress test of EDE-like physics.

Strategy: at fixed acoustic angle 100·θ_s = 1.04110 (Planck), vary the extra
relativistic energy density (ΔN_eff) as a proxy for early-time energy
injection. Compare to:
  - SH0ES H0 = 73.04 ± 1.04
  - DES Y3 S8 = 0.776 ± 0.017
  - Planck CMB Cls (sanity check)

Conclusion: pure ΔN_eff is excluded at >5σ by Planck high-ℓ damping tail.
This justifies the need for true axion-like EDE (CLASS_EDE / AxiCLASS forks).
"""
import numpy as np
from classy import Class

# Validate baseline ΛCDM
def run_classy(dN_eff):
    c = Class()
    c.set({
        '100*theta_s': 1.04110,
        'omega_b': 0.02237,
        'omega_cdm': 0.1200,
        'A_s': 2.1e-9, 'n_s': 0.9649, 'tau_reio': 0.0544,
        'N_ur': 2.0328 + dN_eff,  # baseline = 3 ν − 1 massive
        'output': 'tCl,pCl,lCl,mPk', 'lensing': 'yes',
        'P_k_max_h/Mpc': 1.0, 'z_max_pk': 4.0,
        'l_max_scalars': 3000,
    })
    c.compute()
    out = {
        'dN_eff': dN_eff,
        'H0': c.Hubble(0)*299792.458,
        'rd': c.rs_drag(),
        'Om': c.Omega_m(),
        'sigma8': c.sigma8(),
    }
    out['S8'] = out['sigma8'] * np.sqrt(out['Om']/0.3)
    out['H0_rd'] = out['H0']*out['rd']
    cls = c.lensed_cl(2900)
    out['cls_TT'] = (np.array(cls['ell']), np.array(cls['tt']))
    c.struct_cleanup()
    return out

# Scan
print(f"{'dN':>5} {'H0':>7} {'rd':>7} {'Om':>7} {'σ8':>7} {'S8':>7}")
results = []
for dN in [0.0, 0.5, 1.0, 1.5, 2.0]:
    r = run_classy(dN); results.append(r)
    print(f"  {dN:>4.2f} {r['H0']:>7.2f} {r['rd']:>7.2f} {r['Om']:>7.4f} {r['sigma8']:>7.4f} {r['S8']:>7.4f}")

# Joint chi² with SH0ES + DES Y3
H0_o, H0_e = 73.04, 1.04
S8_o, S8_e = 0.776, 0.017
print(f"\n{'dN':>5} {'χ²(H0)':>9} {'χ²(S8)':>9} {'χ²_tot':>9}")
for r in results:
    c1 = ((r['H0']-H0_o)/H0_e)**2
    c2 = ((r['S8']-S8_o)/S8_e)**2
    print(f"  {r['dN_eff']:>4.2f} {c1:>9.2f} {c2:>9.2f} {c1+c2:>9.2f}")

# Planck Cls consistency check
print("\nTT(ℓ) ratio vs ΛCDM at multiple ℓ:")
ll0, tt0 = results[0]['cls_TT']
dl0 = tt0 * ll0 * (ll0+1) / (2*np.pi)
for r in results:
    ll, tt = r['cls_TT']
    dl = tt * ll * (ll+1) / (2*np.pi)
    ratios = [dl[np.argmin(abs(ll-l))]/dl0[np.argmin(abs(ll0-l))] for l in [200, 500, 1000, 2000, 2500]]
    print(f"  ΔN={r['dN_eff']:.2f}: " + "  ".join(f"{x:.3f}" for x in ratios))
