"""
02_anisotropy_test.py — Directional H0 map from Pantheon+ low-z SNe.
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
    ra  = np.array([float(r[col("RA")])                for r in rows])
    dec = np.array([float(r[col("DEC")])               for r in rows])
    cal = np.array([int(r[col("IS_CALIBRATOR")])       for r in rows])
    return z, mu, err, ra, dec, cal

def fit_H0(z, mu, err, q0=-0.55):
    def chi2(H0):
        dL = (c_kms/H0) * z * (1 + 0.5*(1-q0)*z)
        return np.sum(((mu - 5*np.log10(dL) - 25)/err)**2)
    return minimize_scalar(chi2, bounds=(50,90), method="bounded").x

if __name__ == "__main__":
    z, mu, err, ra, dec, cal = load()
    mask = (cal==0) & (z>0.01) & (z<0.15)
    z, mu, err, ra, dec = [a[mask] for a in (z,mu,err,ra,dec)]
    print(f"Global H0 (z<0.15): {fit_H0(z, mu, err):.2f}")

    rng = np.random.default_rng(0)
    H0_map = []
    for _ in range(3000):
        u = rng.uniform(-1,1); ra_c = rng.uniform(0,360); dec_c = np.degrees(np.arcsin(u))
        cs = (np.sin(np.radians(dec_c))*np.sin(np.radians(dec)) +
              np.cos(np.radians(dec_c))*np.cos(np.radians(dec))*np.cos(np.radians(ra-ra_c)))
        sel = cs > np.cos(np.radians(30))
        if sel.sum() > 30:
            H0_map.append(fit_H0(z[sel], mu[sel], err[sel]))
    arr = np.array(H0_map)
    print(f"Sky map: min={arr.min():.2f}, max={arr.max():.2f}, span={arr.max()-arr.min():.2f}, std={arr.std():.2f}")
