"""Loaders for the data products bundled in this repo.

Sources:
- data/Pantheon+SH0ES.txt        Pantheon+ SNe with SH0ES-calibrator flags
- data/Pantheon+SH0ES_STAT+SYS.txt   full covariance (1701x1701)
- data/Pantheon+SH0ES_STATONLY.txt   diagonal-stat covariance (1701x1701)
- data/lcparam_full_long_zhel.txt    Pantheon (Scolnic 2018) lightcurve params
- data/sys_full_long.txt             Pantheon systematic covariance (1048x1048)

External (in external/):
- external/cobaya_bao_data/desi_bao_dr2/...
- external/SH0ES_Data/                R22 cepheid + HST data
- external/TomoPost-BAO/              tomographic BAO Wang+ 2017
"""
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
EXT = ROOT / "external"


def load_pantheon_plus(stat_only=False):
    """Returns dict with arrays for Pantheon+SH0ES."""
    path = DATA / "Pantheon+SH0ES.txt"
    # genfromtxt with names=True parses the header
    arr = np.genfromtxt(path, names=True, dtype=None, encoding="utf-8")
    cov_path = DATA / ("Pantheon+SH0ES_STATONLY.txt" if stat_only
                       else "Pantheon+SH0ES_STAT+SYS.txt")
    # First line is the matrix size; remaining are flattened entries
    with open(cov_path) as fh:
        n = int(fh.readline().strip())
    flat = np.loadtxt(cov_path, skiprows=1)
    cov = flat.reshape(n, n)
    assert n == arr.shape[0], f"cov dim {n} != data rows {arr.shape[0]}"
    return {
        "data": arr,
        "cov": cov,
        "n": n,
        "zHD": arr["zHD"],
        "zCMB": arr["zCMB"],
        "zHEL": arr["zHEL"],
        "mb": arr["m_b_corr"],
        "mu_sh0es": arr["MU_SH0ES"],
        "ceph_dist": arr["CEPH_DIST"],
        "is_calibrator": arr["IS_CALIBRATOR"].astype(bool),
        "used_in_sh0es_hf": arr["USED_IN_SH0ES_HF"].astype(bool),
    }


def load_pantheon_old():
    """Pantheon (Scolnic 2018) — used by Wang et al. 2024."""
    zcmb, zhel, mb, dmb = np.genfromtxt(
        DATA / "lcparam_full_long_zhel.txt",
        skip_header=1, delimiter=" ", usecols=(1, 2, 4, 5), unpack=True,
    )
    sys_path = DATA / "sys_full_long.txt"
    with open(sys_path) as fh:
        n_sys = int(fh.readline().strip())
    cov = np.loadtxt(sys_path, skiprows=1).reshape(n_sys, n_sys)
    cov_full = cov + np.diag(dmb ** 2)
    return {"zcmb": zcmb, "zhel": zhel, "mb": mb, "dmb": dmb,
            "cov": cov_full, "n": len(zcmb)}


def load_desi_dr2_bao():
    """DESI DR2 BAO measurements (Gaussian likelihood, Cobaya format).

    Returns a dict mapping tracer-bin label -> (mean array, cov array).
    The mean files have rows: value tracer_z observable
    """
    base = EXT / "cobaya_bao_data" / "desi_bao_dr2"
    out = {}
    for f in sorted(base.glob("*_mean.txt")):
        label = f.stem.replace("desi_gaussian_bao_", "").replace("_mean", "")
        # Use object dtype to keep the observable string column
        rows = np.genfromtxt(f, dtype=None, encoding="utf-8")
        # rows can be 1-D structured if single row
        rows = np.atleast_1d(rows)
        cov = np.loadtxt(str(f).replace("_mean.txt", "_cov.txt"))
        if cov.ndim == 0:
            cov = np.array([[float(cov)]])
        out[label] = {"rows": rows, "cov": cov}
    return out


def load_tomographic_bao():
    """Tomographic BAO from Wang et al. 2017 (TomoPost-BAO)."""
    inv = np.loadtxt(EXT / "TomoPost-BAO" / "tomographic_BAO_invcov.txt")
    # tomographic_BAO.txt: z, observable, value, error
    obs = np.genfromtxt(EXT / "TomoPost-BAO" / "tomographic_BAO.txt",
                        dtype=None, encoding="utf-8")
    return {"obs": obs, "invcov": inv}


def load_cepheid_anchors():
    """SH0ES R22 — return path to FITS file with the H0 chain."""
    return EXT / "SH0ES_Data"
