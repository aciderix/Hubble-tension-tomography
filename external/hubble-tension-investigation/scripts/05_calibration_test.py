"""
05_calibration_test.py — TEST C: distance-ladder internal consistency.
Bin H0 by redshift; mass-step test. Falsifies the local-void interpretation
of the anisotropy signal: H0(z) RISES locally instead of falling.
"""
import urllib.request, numpy as np
from scipy.optimize import minimize_scalar

URL = "https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/Pantheon%2B_Data/4_DISTANCES_AND_COVAR/Pantheon%2BSH0ES.dat"
c_kms = 299792.458

def load():
    req = urllib.request.Request(URL, headers={"User-Agent":"Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=30) as r: data = r.read().decode("utf-8")
    lines = data.strip().split("\n"); h = lines[0].split()
    rows = [l.split() for l in lines[1:]]; col = h.index
    z   = np.array([float(r[col("zCMB")])              for r in rows])
    mu  = np.array([float(r[col("MU_SH0ES")])          for r in rows])
    err = np.array([float(r[col("MU_SH0ES_ERR_DIAG")]) for r in rows])
    cal = np.array([int(r[col("IS_CALIBRATOR")])       for r in rows])
    hm  = np.array([float(r[col("HOST_LOGMASS")])      for r in rows])
    return z, mu, err, cal, hm

def fit_H0(z, mu, err, q0=-0.55):
    def chi2(H0):
        dL = (c_kms/H0) * z * (1 + 0.5*(1-q0)*z)
        return np.sum(((mu - 5*np.log10(dL) - 25)/err)**2)
    res = minimize_scalar(chi2, bounds=(50,90), method="bounded")
    h_step = 0.1
    sig = np.sqrt(2/((chi2(res.x+h_step) - 2*chi2(res.x) + chi2(res.x-h_step))/h_step**2))
    return res.x, sig, len(z)

if __name__ == "__main__":
    z, mu, err, cal, hm = load()
    m = (cal==0) & (z>0.005) & (z<0.15)
    z, mu, err, hm = z[m], mu[m], err[m], hm[m]

    print("H0 by redshift bin:")
    for zmin, zmax in [(0.005,0.02),(0.02,0.04),(0.04,0.06),(0.06,0.08),(0.08,0.10),(0.10,0.15)]:
        sel = (z>=zmin) & (z<zmax)
        if sel.sum() > 20:
            H0, s, n = fit_H0(z[sel], mu[sel], err[sel])
            print(f"  z in [{zmin:.3f},{zmax:.3f}] : H0 = {H0:.2f} +/- {s:.2f} (N={n})")

    print("\nMass-step:")
    for label, sel in [("M*<10", hm<10), ("M*>=10", hm>=10)]:
        H0, s, n = fit_H0(z[sel], mu[sel], err[sel])
        print(f"  {label}: H0 = {H0:.2f} +/- {s:.2f} (N={n})")
