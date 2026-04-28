"""classy CLP fluid scan — confirm late-time DE cannot resolve H0 tension. Level 3 Step 2."""
from classy import Class
import numpy as np

def run(params):
    c = Class(); c.set(params); c.compute()
    out = {"H0": c.Hubble(0)*299792.458, "rd": c.rs_drag(),
           "Om": c.Omega_m(), "S8": c.sigma8()*np.sqrt(c.Omega_m()/0.3)}
    c.struct_cleanup(); return out

common = {"h":0.6736,"omega_b":0.02237,"omega_cdm":0.12,
          "A_s":2.1e-9,"n_s":0.9649,"tau_reio":0.0544,
          "output":"tCl,pCl,lCl,mPk","lensing":"yes","l_max_scalars":3000}

print("Baseline LCDM:", run(common))
for wa in [0.0, 0.5, 1.0, 1.5, 2.0]:
    p = dict(common); p.update({"Omega_Lambda":0.0, "w0_fld":-1.0, "wa_fld":wa, "cs2_fld":1.0, "use_ppf":"yes"})
    try:
        print(f"wa={wa}: {run(p)}")
    except Exception as e:
        print(f"wa={wa}: REJECTED ({str(e)[:120]})")
