"""
06_early_dark_energy.py — TEST A: Early Dark Energy effective fit.
Models EDE only via its dominant effect on r_d at recombination.
Result: f_EDE ~ 7% gives H0 ~ 71; ~10% gives H0 ~ 73 (consistent with Poulin+ 2018).

CAVEAT: This is an effective parameterization. A full analysis would solve the
scalar-field dynamics and fit the CMB power spectrum.
"""
import urllib.request, numpy as np
from scipy.optimize import minimize
from scipy.integrate import quad

URL = "https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/Pantheon%2B_Data/4_DISTANCES_AND_COVAR/Pantheon%2BSH0ES.dat"
c_kms = 299792.458; rd_LCDM = 147.05; ALPHA = 0.85
H0rd_P, H0rd_S = 67.4*147.05, 30.0
OM_P, OM_S = 0.315, 0.007
DESI_BAO = [(0.295,7.93,0.15),(0.510,13.62,0.25),(0.706,16.85,0.32),
            (0.930,21.71,0.28),(1.317,27.79,0.69),(1.491,26.07,0.67),(2.330,39.71,0.94)]

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

def E(z, Om): return np.sqrt(Om*(1+z)**3 + (1-Om))
def DM(z, H0, Om): return (c_kms/H0)*quad(lambda zp: 1/E(zp,Om), 0, z, limit=80)[0]
def DV_rd(z, H0, Om, fEDE):
    rd = rd_LCDM*(1 - ALPHA*fEDE)
    dm = DM(z,H0,Om); dh = c_kms/(H0*E(z,Om))
    return ((z*dm**2*dh)**(1/3))/rd
def mu_th(zarr, H0, Om):
    return np.array([5*np.log10((1+zi)*DM(zi,H0,Om))+25 for zi in zarr])

def chi2(p, zsn, msn, esn):
    H0, Om, fEDE = p
    if not (60<H0<80 and 0.25<Om<0.40 and 0<=fEDE<0.20): return 1e10
    try:
        m = mu_th(zsn,H0,Om); w_arr = 1/esn**2
        d = msn-m; M = (w_arr*d).sum()/w_arr.sum()
        chi2_sn  = (w_arr*(d-M)**2).sum()
        chi2_bao = sum(((DV_rd(ze,H0,Om,fEDE)-dv)/s)**2 for ze,dv,s in DESI_BAO)
        chi2_om  = ((Om-OM_P)/OM_S)**2
        rd = rd_LCDM*(1 - ALPHA*fEDE)
        chi2_cmb = ((H0*rd - H0rd_P)/H0rd_S)**2
        return chi2_sn + chi2_bao + chi2_om + chi2_cmb
    except: return 1e10

if __name__ == "__main__":
    z, mu, err, cal = load()
    mask = (cal==0) & (z>0.01) & (z<2.3)
    z, mu, err = z[mask], mu[mask], err[mask]
    rng = np.random.default_rng(42); idx = rng.choice(len(z), 500, replace=False)
    zsn, msn, esn = z[idx], mu[idx], err[idx]

    res_l = minimize(lambda p: chi2([p[0],p[1],0],zsn,msn,esn), [67,0.315],
                     method="Nelder-Mead", options={"xatol":0.005,"fatol":0.1,"maxiter":1000})
    print(f"LCDM (fEDE=0): H0 = {res_l.x[0]:.2f}, Om = {res_l.x[1]:.3f}, chi2 = {res_l.fun:.2f}")

    best = (1e10, None)
    for s in [[67,0.315,0],[70,0.315,0.05],[73,0.315,0.10],[69,0.32,0.03],[72,0.31,0.08]]:
        r = minimize(chi2, s, args=(zsn,msn,esn), method="Nelder-Mead",
                     options={"xatol":0.005,"fatol":0.1,"maxiter":1500})
        if r.fun < best[0]: best = (r.fun, r.x)
    H0, Om, fEDE = best[1]
    rd_eff = rd_LCDM*(1 - ALPHA*fEDE)
    print(f"EDE: H0 = {H0:.2f}, Om = {Om:.3f}, f_EDE = {fEDE:.3f} ({100*fEDE:.1f}%)")
    print(f"     r_d effective = {rd_eff:.2f} Mpc (vs {rd_LCDM:.2f} LCDM)")
    print(f"     Tension vs SH0ES (73.04+/-1.04): {abs(H0-73.04)/np.sqrt(0.2**2+1.04**2):.2f} sigma")

    res_t = minimize(lambda p: chi2([73.04,p[0],p[1]],zsn,msn,esn), [0.315,0.05],
                     method="Nelder-Mead", options={"xatol":0.005,"fatol":0.1,"maxiter":1500})
    print(f"To reach H0=73.04 exactly: f_EDE = {100*res_t.x[1]:.1f}%  (literature: 10-12%)")
