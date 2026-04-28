"""
04_w0wa_DESI.py — TEST B: w0waCDM joint fit on Pantheon+ + DESI BAO + Planck Om prior.
Confirms DESI 2024 w(z) signal but does NOT help the Hubble tension.
"""
import urllib.request, numpy as np
from scipy.optimize import minimize
from scipy.integrate import quad

URL = "https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/Pantheon%2B_Data/4_DISTANCES_AND_COVAR/Pantheon%2BSH0ES.dat"
c_kms = 299792.458; rd = 147.05
DESI_BAO = [(0.295,7.93,0.15),(0.510,13.62,0.25),(0.706,16.85,0.32),
            (0.930,21.71,0.28),(1.317,27.79,0.69),(1.491,26.07,0.67),(2.330,39.71,0.94)]
OM_P, OM_S = 0.315, 0.007

def load():
    req = urllib.request.Request(URL, headers={"User-Agent":"Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=30) as r: data = r.read().decode("utf-8")
    lines = data.strip().split("\n"); h = lines[0].split()
    rows = [l.split() for l in lines[1:]]; col = h.index
    z   = np.array([float(r[col("zCMB")])              for r in rows])
    mu  = np.array([float(r[col("MU_SH0ES")])          for r in rows])
    err = np.array([float(r[col("MU_SH0ES_ERR_DIAG")]) for r in rows])
    cal = np.array([int(r[col("IS_CALIBRATOR")])       for r in rows])
    return z, mu, err, cal

def E(z, Om, w0, wa):
    return np.sqrt(Om*(1+z)**3 + (1-Om)*((1+z)**(3*(1+w0+wa)))*np.exp(-3*wa*z/(1+z)))
def DM(z, H0, Om, w0, wa):
    return (c_kms/H0) * quad(lambda zp: 1/E(zp,Om,w0,wa), 0, z, limit=80)[0]
def DV_rd(z, H0, Om, w0, wa):
    dm = DM(z,H0,Om,w0,wa); dh = c_kms/(H0*E(z,Om,w0,wa))
    return ((z*dm**2*dh)**(1/3))/rd
def mu_th(zarr, H0, Om, w0, wa):
    return np.array([5*np.log10((1+zi)*DM(zi,H0,Om,w0,wa))+25 for zi in zarr])

def chi2(p, zsn, msn, esn):
    H0, Om, w0, wa = p
    if not (50<H0<85 and 0.15<Om<0.50 and -3<w0<0.5 and -3<wa<2): return 1e10
    try:
        m = mu_th(zsn,H0,Om,w0,wa); w_arr = 1/esn**2
        d = msn-m; M = (w_arr*d).sum()/w_arr.sum()
        chi2_sn = (w_arr*(d-M)**2).sum()
        chi2_bao = sum(((DV_rd(ze,H0,Om,w0,wa)-dv)/s)**2 for ze,dv,s in DESI_BAO)
        chi2_om = ((Om-OM_P)/OM_S)**2
        return chi2_sn + chi2_bao + chi2_om
    except: return 1e10

if __name__ == "__main__":
    z, mu, err, cal = load()
    mask = (cal==0) & (z>0.01) & (z<2.3)
    z, mu, err = z[mask], mu[mask], err[mask]
    rng = np.random.default_rng(42); idx = rng.choice(len(z), 500, replace=False)
    zsn, msn, esn = z[idx], mu[idx], err[idx]

    res_l = minimize(lambda p: chi2([p[0],p[1],-1,0],zsn,msn,esn), [70,0.315],
                     method="Nelder-Mead", options={"xatol":0.005,"fatol":0.1,"maxiter":1000})
    print(f"LCDM:    H0 = {res_l.x[0]:.2f}, Om = {res_l.x[1]:.3f}, chi2 = {res_l.fun:.2f}")

    best = (1e10, None)
    for s in [[70,0.315,-1,0],[68,0.31,-0.5,-1.5],[73,0.32,-0.8,-0.5],[69,0.315,-0.45,-1.79]]:
        r = minimize(chi2, s, args=(zsn,msn,esn), method="Nelder-Mead",
                     options={"xatol":0.005,"fatol":0.1,"maxiter":1500})
        if r.fun < best[0]: best = (r.fun, r.x)
    H0,Om,w0,wa = best[1]
    print(f"w0waCDM: H0 = {H0:.2f}, Om = {Om:.3f}, w0 = {w0:.3f}, wa = {wa:.3f}")
    print(f"         chi2 = {best[0]:.2f}, dchi2 vs LCDM = {best[0]-res_l.fun:+.2f}")
