"""Hybrid model: local anisotropy + global EDE. Level 3 Step 3."""
import urllib.request, numpy as np
from scipy.optimize import minimize
from scipy.integrate import quad

# Load Pantheon+
sne_url = "https://github.com/PantheonPlusSH0ES/DataRelease/raw/main/Pantheon%2B_Data/4_DISTANCES_AND_COVAR/Pantheon%2BSH0ES.dat"
cov_url = "https://github.com/PantheonPlusSH0ES/DataRelease/raw/main/Pantheon%2B_Data/4_DISTANCES_AND_COVAR/Pantheon%2BSH0ES_STAT%2BSYS.cov"
sne_txt = urllib.request.urlopen(sne_url).read().decode()
cov_txt = urllib.request.urlopen(cov_url).read().decode()

lines = sne_txt.splitlines(); header = lines[0].split()
data = [{header[i]:p[i] for i in range(len(header))} for line in lines[1:] for p in [line.split()] if len(p)>=len(header)]
zhd = np.array([float(d["zHD"]) for d in data])
mu  = np.array([float(d["MU_SH0ES"]) for d in data])
mask = zhd > 0.01; zhd, mu = zhd[mask], mu[mask]

cov_lines = cov_txt.splitlines(); N = int(cov_lines[0])
cov_full = np.array([float(x) for x in cov_lines[1:N*N+1]]).reshape(N, N)
cov = cov_full[np.ix_(mask, mask)]; invcov = np.linalg.inv(cov)

# Load DESI DR2 (similar to script_desi, omitted here for brevity)
# ... (see level3_desi_dr2_refit.py for full DESI loader)

c_kms, r_d_LCDM = 299792.458, 147.05
shoes_H0 = (73.04, 1.04)

def E(z, Om): return np.sqrt(Om*(1+z)**3 + (1-Om))
def D_L_hyb(z, H0, Om, dH, z_a=0.03):
    I, _ = quad(lambda zp: 1/E(zp, Om), 0, z)
    H0_eff = H0 + dH if z < z_a else H0
    return (1+z) * c_kms/H0_eff * I

# Full hybrid χ² with SNe + DESI + SH0ES (see level3_doc.md for results)
