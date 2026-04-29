"""External H0 anchor likelihoods.

Independent (anchor-free, no rd, no Cepheid) measurements that CAN break
the M2/M1 degeneracy. We implement them as Gaussian priors on the marginal
H0 with appropriate uncertainties.

Datasets (April 2026 best-published values):

  TRGB_F24   Freedman+ 2024 CCHP TRGB-only H0    = 69.85 +/- 1.75
             arXiv:2408.06153
             - Anchor: NGC 4258 megamaser + LMC. NO Cepheids.
             - The "neutral" measurement that should agree with both
               SH0ES and CMB if either is correct.

  SBF_B21    Blakeslee+ 2021 surface-brightness-fluctuation H0 = 73.3 +/- 2.4
             arXiv:2101.02221
             - Calibrated on TRGB. NO Cepheids in the absolute scale.

  TDCOSMO_W  Wong+ 2020 H0LiCOW 6-lens time-delay H0 = 73.3 +1.7/-1.8
             arXiv:1907.04869
             - Time-delay cosmography. Independent of rd AND Cepheid.
             - Assumes NFW + stars mass profile (mass-sheet not relaxed).

  TDCOSMO_B  Birrer+ 2020 (relaxed mass-sheet) H0 = 67.4 +4.1/-3.2
             arXiv:2007.02941
             - Same lenses, but mass-sheet degeneracy is fully accounted for
             - The MOST conservative time-delay cosmography number.

  LVC_O3     Abbott+ 2023 LIGO/Virgo/KAGRA O3 dark+bright sirens
             H0 = 68.0 +12/-8 (asymmetric; modelled as 68.0 +/- 9 here)
             arXiv:2111.03604
             - Pure GR distance from gravitational-wave amplitude.
             - Anchor-free in the full sense.

  ALL_DESC   Combination of TRGB + Mira + Sirens + TDCOSMO Birrer-relaxed:
             ~70 +/- 1   (a "consensus anchor-free" value)
"""
import numpy as np

# (H0, sigma)  ; asymmetric upper/lower errors symmetrised for Gaussian prior
ANCHORS = {
    "TRGB_F24":    (69.85, 1.75),
    "SBF_B21":     (73.30, 2.40),
    "TDCOSMO_W":   (73.30, 1.80),
    "TDCOSMO_B":   (67.40, 3.65),
    "LVC_O3":      (68.00, 9.00),
    "TRGB+Mira+Sirens": (70.00, 1.00),  # synthetic consensus from Freedman 2024 review
}


def chi2_anchor(H0, anchor_name):
    H0_obs, sig = ANCHORS[anchor_name]
    return ((H0 - H0_obs) / sig) ** 2


def chi2_anchors(H0, anchor_names):
    """Sum chi2 from a list of anchors."""
    return sum(chi2_anchor(H0, a) for a in anchor_names)


def all_anchors():
    return list(ANCHORS.keys())
