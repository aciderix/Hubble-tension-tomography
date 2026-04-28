"""DESI DR2 BAO Gaussian likelihood.

Each measurement is a tracer-bin row (z, value, observable) with a
covariance. Observables encode which combination of distances:
  DM_over_rs : D_M(z) / r_d
  DH_over_rs : D_H(z) / r_d  with D_H = c/H(z)
  DV_over_rs : D_V(z) / r_d  with D_V = (z * D_M^2 * c/H)^(1/3)
"""
from pathlib import Path
import numpy as np

from .cosmo import comoving_distance, DH, DV, comoving_distance_w0wa, C_KMS

ROOT = Path(__file__).resolve().parents[1]
DR2 = ROOT / "external" / "cobaya_bao_data" / "desi_bao_dr2"
DR1 = ROOT / "external" / "cobaya_bao_data"


def _load_block(mean_file):
    cov = np.loadtxt(str(mean_file).replace("_mean.txt", "_cov.txt"))
    rows = np.atleast_1d(np.genfromtxt(mean_file, dtype=None, encoding="utf-8"))
    z = np.array([r[0] for r in rows])
    val = np.array([r[1] for r in rows])
    obs = np.array([r[2] for r in rows])
    if cov.ndim == 0:
        cov = np.array([[float(cov)]])
    return z, val, obs, cov


def load_desi_dr2_all():
    """Returns the combined DESI DR2 vector (13 measurements)."""
    z, val, obs, cov = _load_block(DR2 / "desi_gaussian_bao_ALL_GCcomb_mean.txt")
    return {"z": z, "val": val, "obs": obs, "cov": cov,
            "invcov": np.linalg.inv(cov), "n": len(z)}


def load_desi_dr1_all():
    z, val, obs, cov = _load_block(
        DR1 / "desi_2024_gaussian_bao_ALL_GCcomb_mean.txt")
    return {"z": z, "val": val, "obs": obs, "cov": cov,
            "invcov": np.linalg.inv(cov), "n": len(z)}


def predict_obs(z_arr, obs_arr, H0, Om, rd):
    """Predict the BAO observables under flat LCDM with given (H0, Om, rd)."""
    out = np.empty_like(z_arr, dtype=float)
    # Group by z to avoid recomputing comoving_distance
    unique_z = np.unique(z_arr)
    dM = {z: comoving_distance(np.array([z]), H0, Om)[0] for z in unique_z}
    Hz = {z: H0 * np.sqrt(Om * (1 + z) ** 3 + (1 - Om)) for z in unique_z}
    for i, (z, ob) in enumerate(zip(z_arr, obs_arr)):
        if ob == "DM_over_rs":
            out[i] = dM[z] / rd
        elif ob == "DH_over_rs":
            out[i] = (C_KMS / Hz[z]) / rd
        elif ob == "DV_over_rs":
            dv = (z * dM[z] ** 2 * C_KMS / Hz[z]) ** (1.0 / 3.0)
            out[i] = dv / rd
        else:
            raise ValueError(f"unknown observable {ob}")
    return out


def chi2_desi(H0, Om, rd, dr2=True):
    bao = load_desi_dr2_all() if dr2 else load_desi_dr1_all()
    pred = predict_obs(bao["z"], bao["obs"], H0, Om, rd)
    d = bao["val"] - pred
    return float(d @ bao["invcov"] @ d)


def predict_obs_w0wa(z_arr, obs_arr, H0, Om, w0, wa, rd):
    out = np.empty_like(z_arr, dtype=float)
    unique_z = np.unique(z_arr)
    from .cosmo import E_w0wa
    dM = {z: comoving_distance_w0wa(np.array([z]), H0, Om, w0, wa)[0]
          for z in unique_z}
    Hz = {z: H0 * E_w0wa(np.array([z]), Om, w0, wa)[0] for z in unique_z}
    for i, (z, ob) in enumerate(zip(z_arr, obs_arr)):
        if ob == "DM_over_rs":
            out[i] = dM[z] / rd
        elif ob == "DH_over_rs":
            out[i] = (C_KMS / Hz[z]) / rd
        elif ob == "DV_over_rs":
            dv = (z * dM[z] ** 2 * C_KMS / Hz[z]) ** (1.0 / 3.0)
            out[i] = dv / rd
    return out
