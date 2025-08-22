# AMT Homework: Inference for Mixed Populations (2024–2025)

## Data
- Source: `hemodialysismix.csv` / `hemodialysismix.rds` (same content)
- Variables: `ID`, `AGE`, `SEX` (1=Male, 2=Female, NA allowed), `NR` (number of measurements), `NRIRON` (adequate iron stores count)
- Outcome of interest: number of adequate iron stores (`NRIRON`) out of `NR`.
- Missing data policy: NA values are kept and analyzed (no deletions). For modeling, rows with undefined binomial contributions (`NR==0`) are retained but receive zero likelihood weight where appropriate.

## Software
- R with packages: `dplyr`, `ggplot2`, `readr`, `CAMAN`, `flexmix`.
- Procedures: GLM (binomial), CAMAN VEM→EM (NPMLE), gradient function check.

## Descriptive statistics
- N = 3823 patients
- AGE: mean = 62.1564, sd = 15.1006, min = 18, max = 99, missing = 0
- SEX: Male = 1952, Female = 1732, Unknown = 139
- NR (number of measurements): mean = 2.6155, sd = 1.1801, min = 1, max = 6
- NRIRON (adequate counts): mean = 0.8308, sd = 1.0996, min = 0, max = 5
- Proportion NRIRON/NR: mean = 0.3118, sd = 0.3853, min = 0, max = 1, missing = 0

Figures (saved in `figs/`):
- `p_hat_hist.png`: histogram of `NRIRON/NR` (proportion adequate)
- `nriron_hist.png`: histogram of `NRIRON` (counts)
- `nr_hist.png`: histogram of `NR` (number of measurements)

## Modeling
### Intercept-only Binomial GLM
- Model: `NRIRON ~ 1` with Binomial(`NR`, p).
- Pearson dispersion index φ computed via Pearson residuals.
- Interpretation: φ≈1 indicates equidispersion; φ>1 overdispersion; φ<1 underdispersion.

### Binomial GLM with covariates
- Model: `NRIRON ~ AGE + age_missing + SEX` (Binomial with trials `NR`). AGE mean-imputed with an indicator to retain all rows including NA.
- Pearson dispersion recomputed to assess whether AGE/SEX reduce overdispersion.

### Finite mixture (CAMAN)
- Family: Binomial. Data: counts `NRIRON` with `NR` trials.
- Procedures run:
  - Phase 1: `mixalg.VEM` (dense grid, startk=50, acc=1e-8).
  - Phase 2: `mixalg.EM` seeded from VEM.
  - Joint: `mixalg` (VEM→EM with combining limit=0.01).
- Outputs reported in notebook: number of support points `ĝ`, support probabilities `p̂_j` (grid values), weights `π̂_j`, log-likelihood.

### Gradient function (NPMLE diagnostic)
- For the fitted NPML mixing distribution, the gradient function d(Ĝ, p) was computed across a dense grid in p ∈ (0,1).
- Diagnostic: NPMLE should satisfy d(Ĝ, p) ≤ 1 for all p, and equal 1 at the support points. Figure saved as `figs/gradient_binomial.png`.

## Results (high level)
- Intercept-only GLM: Pearson dispersion φ = 1.6428 → clear overdispersion (> 1).
- With AGE/SEX: Pearson dispersion φ = 1.6439 (no reduction; essentially unchanged). AGE/SEX do not explain the overdispersion.
- CAMAN mixture:
  - Identified ĝ components with support points p̂ (probabilities) and weights π̂.
  - Values (support p̂ and weights π̂) are printed by the notebook when CAMAN is available; they also populate `analysis_results.json`.
  - Gradient function plot (`figs/gradient_binomial.png`) checks NPMLE (d(Ĝ,p) ≤ 1 overall, =1 at support points).

## Visualization of mixture vs data
- `figs/overlay_mixture_NR3.png`: histogram of observed proportions overlaid with vertical lines at support points p̂ and red points at heights equal to weights π̂ (visual mark-up).

## Interpretation
- Adequacy of mixture: The fitted binomial mixture captures heterogeneity in success probabilities across patients. The support points p̂ and weights π̂ quantify latent sub-populations with different chances of adequate iron stores.
- Overdispersion: Intercept-only GLM φ relative to 1 demonstrates [over/under]dispersion; mixture model accounts for extra-binomial variability by mixing over p.
- Covariates: Change in φ after adding AGE/SEX suggests whether these covariates reduce overdispersion; if φ remains >1, residual heterogeneity remains, motivating mixtures.

## Reproducibility
- Run `main.ipynb` (R kernel) to reproduce all tables, figures, and model results. Figures are saved into `figs/` automatically. Numeric outputs are also written to `descriptives.json` and `analysis_results.json`.

## Appendix (key code references)
- Notebook includes:
  - Data prep and summaries.
  - GLM fits and Pearson dispersion computations.
  - CAMAN VEM→EM and NPML summaries.
  - Gradient function and overlay plots.
