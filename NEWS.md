# simcow 0.2.0

## Breaking changes

* `create_cow()` now computes the weights of the paper (Behrens et al., under revision at
  PLOS ONE). For each external cohort g, a logistic membership model (target cohort vs
  cohort g) uses the covariates and the outcome. External observation j gets the weight
  `J_g * p_jg / mean(p_g)`: its fitted probability of target membership, divided by the mean
  within its cohort, times the subgroup weight `J_g = 1 - |2 AUC_g - 1|`. `AUC_g` comes from
  5-fold cross-validated membership predictions. Target observations get weight 1.
* New arguments: `outcome` (outcome column for the membership models), `K_folds` (folds of the
  cross-validated AUC), and `seed` (seed for the fold assignment; leaves the global random
  number stream unchanged).
* Removed arguments: `trunc_weights`, `AUC`, and `overlap_threshold`.
* The returned data keep the input row order (0.1.0 returned the external cohorts first and
  the target cohort last). `ps_weight` is NA for target observations. The attribute
  `"subgroups"` gives `n`, `auc`, and `J` per external cohort.

## Old behavior (0.1.0)

`create_cow()` 0.1.0 weighted external observation j by `min(1, p_jg / AUC_g)` (without the
cap if `trunc_weights = FALSE`), with `p_jg` from a membership model on the covariates in
`baseline_covs` and `AUC_g` from the fitted (not cross-validated) membership predictions
(`AUC = FALSE` dropped the AUC). `overlap_threshold` set the weights below a quantile to 0.
To use it, install version 0.1.0:
`remotes::install_github("maxjonasbehrens/simcow@c1182c9")`.

## Other changes

* `paper/` holds the scripts and the aggregated simulation results that reproduce the
  paper's simulation study, Tables 1 and 2, and Figs 1-3, and the code of the clinical
  application (the clinical data are not public).
* Unit tests compare `create_cow()` with the reference implementation of the paper.

# simcow 0.1.0

* First version.
