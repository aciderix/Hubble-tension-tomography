"""DESI DR2 + Pantheon+ Om-prior toy-EDE refit. Level 3 Step 1."""
import urllib.request, numpy as np, json
from scipy.optimize import minimize
from scipy.integrate import quad

mean_url = "https://raw.githubusercontent.com/CobayaSampler/bao_data/master/desi_bao_dr2/desi_gaussian_bao_ALL_GCcomb_mean.txt"
cov_url  = "https://raw.githubusercontent.com/CobayaSampler/bao_data/master/desi_bao_dr2/desi_gaussian_bao_ALL_GCcomb_cov.txt"

mean_txt = urllib.request.urlopen(mean_url).read().decode()
cov_txt  = urllib.request.urlopen(cov_url).read().decode()

rows = []
for line in mean_txt.splitlines():
    line = line.strip()
    if not line or line.startswith("#"): continue
    p = line.split(); rows.append((float(p[0]), float(p[1]), p[2]))
desi_z   = np.array([r[0] for r in rows])
desi_obs = np.array([r[1] for r in rows])
desi_q   = [r[2] for r in rows]
desi_cov = np.array([[float(x) for x in line.split()] for line in cov_txt.splitlines() if line.strip()])
desi_invcov = np.linalg.inv(desi_cov)

c_kms, r_d_LCDM = 299792.458, 147.05
sne_Om   = (0.334, 0.018)
shoes_H0 = (73.04, 1.04)

def E(z, Om):     return np.sqrt(Om*(1+z)**3 + (1-Om))
def D_M(z, H0, Om):
    I, _ = quad(lambda zp: 1/E(zp, Om), 0, z); return c_kms/H0 * I
def D_H(z, H0, Om): return c_kms/(H0*E(z, Om))
def D_V(z, H0, Om): return (z*D_M(z, H0, Om)**2 * D_H(z, H0, Om))**(1/3)

def desi_chi2(H0, Om, r_d_eff):
    pred = []
    for z, q in zip(desi_z, desi_q):
        if   q == "DV_over_rs": pred.append(D_V(z, H0, Om)/r_d_eff)
        elif q == "DM_over_rs": pred.append(D_M(z, H0, Om)/r_d_eff)
        elif q == "DH_over_rs": pred.append(D_H(z, H0, Om)/r_d_eff)
    diff = desi_obs - np.array(pred)
    return diff @ desi_invcov @ diff

def total_chi2(params, alpha=0.18, with_shoes=True):
    H0, Om, f_EDE = params
    if not (50 < H0 < 90): return 1e10
    if not (0.1 < Om < 0.5): return 1e10
    if not (0 <= f_EDE < 0.5): return 1e10
    r_d_eff = r_d_LCDM * (1 - alpha*f_EDE)
    chi2 = desi_chi2(H0, Om, r_d_eff)
    chi2 += ((Om - sne_Om[0])/sne_Om[1])**2
    if with_shoes:
        chi2 += ((H0 - shoes_H0[0])/shoes_H0[1])**2
    return chi2

if __name__ == "__main__":
    r1 = minimize(lambda p: total_chi2([p[0], p[1], 0.0], with_shoes=False), x0=[67, 0.31], method="Nelder-Mead")
    r2 = minimize(lambda p: total_chi2([p[0], p[1], 0.0], with_shoes=True),  x0=[70, 0.31], method="Nelder-Mead")
    r3 = minimize(lambda p: total_chi2(p, with_shoes=True), x0=[73, 0.30, 0.08], method="Nelder-Mead")
    r4 = minimize(lambda p: total_chi2(p, with_shoes=False), x0=[68, 0.31, 0.0], method="Nelder-Mead")
    print(f"LCDM no SH0ES:  H0={r1.x[0]:.2f}  Om={r1.x[1]:.4f}  chi2={r1.fun:.2f}")
    print(f"LCDM +SH0ES:    H0={r2.x[0]:.2f}  Om={r2.x[1]:.4f}  chi2={r2.fun:.2f}")
    print(f"EDE  +SH0ES:    H0={r3.x[0]:.2f}  Om={r3.x[1]:.4f}  f_EDE={r3.x[2]:.4f}  chi2={r3.fun:.2f}")
    print(f"EDE  no SH0ES:  H0={r4.x[0]:.2f}  Om={r4.x[1]:.4f}  f_EDE={r4.x[2]:.4f}  chi2={r4.fun:.2f}")
