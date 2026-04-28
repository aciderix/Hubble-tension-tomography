"""Fast joint likelihood for Pantheon+SH0ES + DESI DR2 BAO.

- Pre-compute inv-cov of Pantheon+ once, store residual quadratic form
  as resid @ Cinv @ resid.
- Use grid-interpolated comoving distance (1e-5 accuracy).
- Each chi2 evaluation is sub-millisecond.
"""
from pathlib import Path
import numpy as np

from .data_loaders import load_pantheon_plus
from .desi_bao import load_desi_dr2_all
from .cosmo import (distance_modulus_grid, comoving_distance_grid,
                    rd_aubourg, C_KMS, E_LCDM)


class HubbleLH:
    def __init__(self):
        pp = load_pantheon_plus()
        self.zHEL = pp["zHEL"]
        self.zHD = pp["zHD"]
        self.mb = pp["mb"]
        self.is_calib = pp["is_calibrator"]
        self.ceph_dist = pp["ceph_dist"]
        self.cov_pp = pp["cov"]
        self.Cinv_pp = np.linalg.inv(self.cov_pp)
        self.n_pp = pp["n"]

        bao = load_desi_dr2_all()
        self.bao_z = bao["z"]
        self.bao_val = bao["val"]
        self.bao_obs = bao["obs"]
        self.bao_invcov = bao["invcov"]
        self.n_bao = bao["n"]
        self.bao_unique_z = np.unique(self.bao_z)

    # ---------- Pantheon+SH0ES ----------
    def chi2_pp(self, H0, Om, M, use_calib=True, M_calib=None):
        """If M_calib is None, the calibrators use the same M; else a separate
        calibrator-only nuisance (tests a calibrator-vs-Hubble-flow offset)."""
        mu_th = distance_modulus_grid(self.zHEL, self.zHD, H0, Om)
        mu_model = np.where(self.is_calib, self.ceph_dist, mu_th)
        if M_calib is not None:
            M_eff = np.where(self.is_calib, M_calib, M)
        else:
            M_eff = M
        resid = self.mb - (M_eff + mu_model)
        if not use_calib:
            m = ~self.is_calib
            return float(resid[m] @ self.Cinv_pp[np.ix_(m, m)] @ resid[m])
        return float(resid @ self.Cinv_pp @ resid)

    # ---------- DESI BAO ----------
    def _bao_pred_lcdm(self, H0, Om, rd):
        # Get D_M at unique z
        dM = comoving_distance_grid(self.bao_unique_z, H0, Om)
        Hz = H0 * E_LCDM(self.bao_unique_z, Om)
        dM_map = dict(zip(self.bao_unique_z, dM))
        Hz_map = dict(zip(self.bao_unique_z, Hz))
        out = np.empty(self.n_bao)
        for i, (z, ob) in enumerate(zip(self.bao_z, self.bao_obs)):
            if ob == "DM_over_rs":
                out[i] = dM_map[z] / rd
            elif ob == "DH_over_rs":
                out[i] = (C_KMS / Hz_map[z]) / rd
            elif ob == "DV_over_rs":
                dv = (z * dM_map[z] ** 2 * C_KMS / Hz_map[z]) ** (1.0 / 3.0)
                out[i] = dv / rd
        return out

    def chi2_bao(self, H0, Om, rd):
        d = self.bao_val - self._bao_pred_lcdm(H0, Om, rd)
        return float(d @ self.bao_invcov @ d)

    # ---------- CMB compressed (Planck 2018 distance priors) ----------
    def chi2_cmb(self, H0, Om, obh2, rd=None, Neff=3.046):
        """Compressed Planck CMB (R, lA, omega_b). If rd is None, uses
        Brieden anchored fitting form; else uses supplied rd directly.
        """
        from .cmb_compressed import chi2_cmb as _cmb, chi2_cmb_with_rd
        if rd is None:
            return _cmb(H0, Om, obh2, Neff=Neff)
        return chi2_cmb_with_rd(H0, Om, obh2, rd, Neff=Neff)

    # ---------- Joint ----------
    def chi2_joint(self, H0, Om, M, rd, use_calib=True):
        return self.chi2_pp(H0, Om, M, use_calib=use_calib) + self.chi2_bao(H0, Om, rd)

    def chi2_full(self, H0, Om, M, rd, obh2, use_calib=True, Neff=3.046):
        """Full SN + BAO + CMB-compressed joint chi2.

        rd is supplied externally; if you want LCDM-consistent rd, pass
        rd_planck_anchored(Om*h^2, obh2). The CMB term then "sees" the same rd
        used for BAO, which is the self-consistent LCDM case.
        """
        return (self.chi2_pp(H0, Om, M, use_calib=use_calib)
                + self.chi2_bao(H0, Om, rd)
                + self.chi2_cmb(H0, Om, obh2, rd=rd, Neff=Neff))
