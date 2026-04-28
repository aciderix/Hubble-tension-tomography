# Level 1 — Anisotropy Diagnostic at z<0.03

## Question
Pantheon+ shows a residual H₀ anisotropy of ~2.5σ at z<0.03 even **after**
peculiar-velocity correction (zHD vs zCMB redshifts). Is this:
- (A) A real physical signal (low-z anisotropic expansion)?
- (B) A calibration systematic from a specific photometric survey?

## Method
Two complementary tests on the N=386 Pantheon+ SNe at z<0.03 (non-calibrators):

### Test 1 — Per-survey H₀
Compute H₀ separately within each photometric survey contributing >20 SNe.
If anisotropy is calibration-driven, one survey should show an outlier H₀.

### Test 2 — Drop-one-survey
Recompute the full-sky H₀ span (max − min over hemispheres) after removing
each individual survey. If anisotropy is dominated by one survey, dropping
it should collapse the span.

## Results

### Per-survey H₀ (z<0.03, N>20)

| Survey ID | N    | H₀ [km/s/Mpc] |
|-----------|------|---------------|
| 5         | 47   | 73.42         |
| 18        | 31   | 74.83         |
| 50        | 24   | 73.05         |
| 51        | 38   | 73.91         |
| 56        | 28   | 72.95         |
| 57        | 35   | 73.66         |
| 61        | 22   | 74.12         |
| 62        | 26   | 73.58         |
| 63        | 33   | 74.02         |
| 64        | 21   | 73.77         |
| 65        | 25   | 73.31         |
| 66        | 29   | 74.18         |
| 150       | 27   | 73.49         |

**Range: 72.95 — 74.83** (window of 1.88 km/s/Mpc, all surveys mutually
compatible within their ~1.5 km/s/Mpc statistical uncertainties).

### Drop-one-survey H₀ span

| Survey dropped | H₀ span | Δ vs baseline |
|----------------|---------|---------------|
| (none)         | 4.88    | —             |
| 5              | 5.12    | +0.24         |
| 18             | 4.49    | −0.39         |
| 50             | 5.31    | +0.43         |
| 51             | 5.93    | +1.05         |
| 56             | 4.78    | −0.10         |
| 57             | 5.02    | +0.14         |
| 61             | 5.41    | +0.53         |
| 62             | 4.61    | −0.27         |
| 63             | 4.94    | +0.06         |
| 64             | 5.18    | +0.30         |
| 65             | 4.83    | −0.05         |
| 66             | 5.07    | +0.19         |
| 150            | 4.71    | −0.17         |

Span remains in **[4.49, 5.93]** regardless of which survey is removed.
No single survey dominates the signal.

## Interpretation

**The residual 2.5σ anisotropy is NOT a single-survey calibration artifact.**

If it were, dropping that survey would shrink the span to <2 km/s/Mpc.
Instead, the signal is **distributed** across all 13 contributing surveys.

Two physical interpretations remain consistent with this:
1. **Sub-dominant real anisotropy** — small but genuine departure from FLRW
   isotropy at the local scale (z<0.03 ≈ 130 Mpc).
2. **Common low-level systematic** — affecting all surveys similarly (e.g.
   the Pantheon+ PV correction model itself has residual bias).

With only 386 SNe at z<0.03, the 2.5σ signal cannot be promoted to a
detection. **A future test would require independent SNe samples (DES,
ZTF, LSST) with different selection functions.**

## Files
- `scripts/anisotropy_per_survey.py` — Test 1 implementation
- `scripts/anisotropy_drop_one.py` — Test 2 implementation
- `results/anisotropy_z_lt_003.json` — raw outputs

## Conclusion for Level 2
Anisotropy is real-but-marginal. Not a calibration artifact, not a
publication-grade detection. **Level 2 proceeds without modeling
anisotropy** (focus on early-time physics: EDE).
