"""
Level 1 — Per-survey H0 in z<0.03 Pantheon+ subset.

Tests whether the residual ~2.5σ H0 anisotropy (post peculiar-velocity
correction) is driven by a single photometric survey calibration error,
or is distributed across multiple independent samples.
"""
import numpy as np
import pandas as pd

# Load Pantheon+ (zHD already PV-corrected by 2M++ flow model)
data = pd.read_csv("Pantheon+SH0ES.dat", sep=r"\s+")
mask = (data['zHD'] >= 0.01) & (data['zHD'] < 0.03) & (data['IS_CALIBRATOR'] == 0)
sub = data[mask].copy()
print(f"N (z<0.03, non-calibrators): {len(sub)}")

# Per-survey H0 fit (assuming flat ΛCDM with Om=0.315 fixed)
def H0_from_subset(z, mu, mu_err):
    """Fit H0 from distance-modulus relation:
        mu = 5*log10(d_L/Mpc) + 25
        d_L = (c/H0)*(1+z)*∫dz'/E(z')
    For Om fixed, this becomes mu = M(Om,z) + 5*log10(H0_ref/H0)
    """
    from scipy.optimize import minimize_scalar
    Om = 0.315; c = 299792.458
    def E(z): return np.sqrt(Om*(1+z)**3 + (1-Om))
    def d_L(z, H0):
        from scipy.integrate import quad
        I = np.array([quad(lambda zp: 1/E(zp), 0, zi)[0] for zi in z])
        return c/H0*(1+z)*I
    def chi2(H0):
        mu_th = 5*np.log10(d_L(z, H0)) + 25
        return np.sum(((mu - mu_th)/mu_err)**2)
    res = minimize_scalar(chi2, bounds=(50, 90), method='bounded')
    return res.x

# Group by survey
results = []
for sid, grp in sub.groupby('IDSURVEY'):
    if len(grp) < 20: continue
    H0 = H0_from_subset(grp['zHD'].values, grp['m_b_corr'].values - grp['M_B'].mean(), grp['MU_SH0ES_ERR_DIAG'].values)
    results.append((int(sid), len(grp), H0))
    print(f"  Survey {sid:>4}: N={len(grp):>3}, H0={H0:.2f}")

H0s = [r[2] for r in results]
print(f"\nSpan: {min(H0s):.2f} — {max(H0s):.2f} ({max(H0s)-min(H0s):.2f} km/s/Mpc)")
print(f"Std:  {np.std(H0s):.2f}")
