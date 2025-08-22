# AMT: Homework Assignment — Hemodialysis Mixture Analysis

## (Section 1) Data description (no modeling)
- Source: `hemodialysismix.csv` (patients on hemodialysis)
- Variables: ID, AGE (years), SEX (1=Male, 2=Female, NA kept as Unknown), NR (number of measurements), NRIRON (adequate iron counts)
- N subjects: 3823
- Distribution of follow-up counts NR:
# A tibble: 6 × 2
     NR     n
  <int> <int>
1     1   880
2     2   828
3     3  1180
4     4   780
5     5   127
6     6    28
- `SEX` counts:
# A tibble: 3 × 2
  SEX_factor     n
  <fct>      <int>
1 Male        1952
2 Female      1732
3 Unknown      139
- Proportion p = NRIRON/NR summary:
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 0.0000  0.0000  0.0000  0.3118  0.6667  1.0000 

Figures (Section 1):
- figs/p_hat_hist.png
- figs/nriron_hist.png
- figs/nr_hist.png

## (Section 2) Model-based analysis
- Intercept-only Binomial GLM (unimodal) Pearson dispersion:
  - phi = 1.643 (>1 indicates overdispersion)
  - Interpretation: Overdispersion present.
- Binomial GLM with AGE and SEX (with missingness indicator):
  - phi = 1.644
  - Overdispersion persists after covariate adjustment.
- Finite mixture of Binomials (custom EM baseline):
  - Selected k by BIC: 3
  - Component probabilities (p) and weights (pi):
  component          p        pi
1         1 0.02398794 0.4070644
2         2 0.37152767 0.4099659
3         3 0.83802367 0.1829697
- CAMAN results unavailable (package not installed or fit failed).
- Gradient function diagnostic unavailable.
- Mixture-predicted diagnostics unavailable.

## (Section 3) Finalization
- Overlay of mixture model on empirical counts histogram (most common NR):
  - figs/overlay_mixture_NR3.png

## Software and estimation
- R, glm (binomial), flexmix (baseline mixture EM), CAMAN (VEM/EM/NPML)
- NAs were retained: `SEX` as explicit 'Unknown' level; `AGE` imputed with mean plus a missingness indicator to keep subjects in all analyses.

## Communication-ready takeaway
- Monthly adequate-iron counts are overdispersed relative to a single-probability Binomial.
- A 3-component finite mixture (baseline) captures heterogeneity.
- AGE/SEX do not fully explain heterogeneity (overdispersion persists).

(See notebook `main.ipynb` for figures.)
