"""
Level 1 — Drop-one-survey robustness test for low-z H0 anisotropy.

For each photometric survey contributing to the z<0.03 sample, remove it
and recompute the full-sky H0 hemisphere span. If anisotropy is driven by
a single bad-calibration survey, dropping it should collapse the span.
"""
import numpy as np
import pandas as pd

def hemisphere_H0_span(z, mu, mu_err, ra, dec, n_directions=20):
    """Return (max H0 − min H0) over n random sky hemispheres."""
    from scipy.optimize import minimize_scalar
    from scipy.integrate import quad
    Om, c = 0.315, 299792.458
    def E(zp): return np.sqrt(Om*(1+zp)**3 + (1-Om))
    def d_L(z, H0):
        I = np.array([quad(lambda zp: 1/E(zp), 0, zi)[0] for zi in z])
        return c/H0*(1+z)*I
    def fit_H0(zs, mus, errs):
        def chi2(H0):
            mu_th = 5*np.log10(d_L(zs, H0)) + 25
            return np.sum(((mus - mu_th)/errs)**2)
        return minimize_scalar(chi2, bounds=(50,90), method='bounded').x

    # Convert RA/Dec to unit vectors
    cra, cdec = np.cos(np.radians(ra)), np.cos(np.radians(dec))
    sra, sdec = np.sin(np.radians(ra)), np.sin(np.radians(dec))
    pts = np.stack([cdec*cra, cdec*sra, sdec], axis=1)

    H0_max, H0_min = -np.inf, np.inf
    rng = np.random.default_rng(42)
    for _ in range(n_directions):
        n = rng.normal(size=3); n /= np.linalg.norm(n)
        m_pos = pts @ n > 0
        if m_pos.sum() < 30 or (~m_pos).sum() < 30: continue
        h_pos = fit_H0(z[m_pos], mu[m_pos], mu_err[m_pos])
        h_neg = fit_H0(z[~m_pos], mu[~m_pos], mu_err[~m_pos])
        H0_max = max(H0_max, h_pos, h_neg)
        H0_min = min(H0_min, h_pos, h_neg)
    return H0_max - H0_min

data = pd.read_csv("Pantheon+SH0ES.dat", sep=r"\s+")
mask = (data['zHD'] >= 0.01) & (data['zHD'] < 0.03) & (data['IS_CALIBRATOR'] == 0)
sub = data[mask].copy()

# Baseline span
span0 = hemisphere_H0_span(sub['zHD'].values, sub['m_b_corr'].values - sub['M_B'].mean(),
                           sub['MU_SH0ES_ERR_DIAG'].values, sub['RA'].values, sub['DEC'].values)
print(f"Baseline span (all surveys): {span0:.2f} km/s/Mpc\n")

# Drop-one
print("Drop-one-survey results:")
for sid in sorted(sub['IDSURVEY'].unique()):
    grp = sub[sub['IDSURVEY'] != sid]
    if len(grp) < 100: continue
    span = hemisphere_H0_span(grp['zHD'].values, grp['m_b_corr'].values - grp['M_B'].mean(),
                              grp['MU_SH0ES_ERR_DIAG'].values, grp['RA'].values, grp['DEC'].values)
    print(f"  Drop survey {sid:>4}: span={span:.2f}, Δ={span-span0:+.2f}")
