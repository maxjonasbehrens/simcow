# Reproducing the paper

Scripts and aggregated results for Behrens et al., *Improving prediction models by incorporating
external data with weights based on similarity* (under revision at PLOS ONE). The weighted
sample uses `create_cow()` from this package (`R/create_cow.R`, sourced directly, so no
installation is needed).

## Requirements

R 4.1 or later and the packages dplyr, lme4, tidyr, ggplot2, gridExtra, and ragg (readr for the
clinical application). The results in `results/` were computed with R 4.4.0 on macOS; see the
end of this file for the package versions.

## Simulation study

Run from the repository root (or any directory):

```sh
Rscript paper/run_simulation.R   # Table 1, bias-variance decomposition
Rscript paper/run_table2.R       # Table 2 (needs the data sets cached by run_simulation.R)
Rscript paper/make_fig1.R        # Fig 1
Rscript paper/make_figures.R     # Figs 2 and 3 (Fig 2B needs the cached data sets)
```

The scripts write to `paper/output/` (not tracked by git). `run_simulation.R` replays the
random number stream of the paper (`set.seed(123)`), caches the 32,400 simulated data sets in
`paper/output/cache/`, and stops with `CHECK FAILED` unless its results equal those in
`results/` (absolute difference below 1e-9). `run_table2.R` checks the same. The cache takes
about 600 MB. Run times on an Apple Silicon Mac with 6 cores, with both scripts running at the
same time: `run_simulation.R` 16 min (3.5 min of it replay the data sets), `run_table2.R`
10 min; the figure scripts take under 1 min. On Windows the scripts run on one core.

| File | Script | Content |
|---|---|---|
| `../R/create_cow.R` (package) | all | Weights of the paper (method D) |
| `R/methods.R` | all | Weight variants A–D and the covariate-shift comparator |
| `R/simulate_scenario.R` | simulation | Data-generating process (Fig 1) |
| `results/sim_runs.csv` | `run_simulation.R` | One row per run (32,400): RMSE and ESS ratio per method |
| `results/sim_subgroups.csv` | `run_simulation.R` | One row per external subgroup and run: size, J_g, shift |
| `results/bias_variance.csv` | `run_simulation.R` | Bias-variance decomposition per scenario and method |
| `results/table2_runs.csv` | `run_table2.R` | One row per run: RMSE per weight-handling method |
| `results/table2_summary.csv` | `run_table2.R` | Table 2 |

Codes in the result files:

- `scenario`: `AX` covariate shift, `AY` outcome shift, `AH` covariate + outcome shift.
- `to`: largest shift of an external subgroup; 3 dissimilar, 2 medium, 1 similar.
- `n_cat`: number of subgroups (target included); `n_train`: training observations per
  subgroup; `n_target`: target observations; `n0`: target observations in the training data.
- Weight variants (`rmse_*`, `ess_*`, `J_*`): `D` the weights of the paper (`create_cow()`,
  "weighted" in the paper); `A` without J_g; `B` J_g from the in-sample AUC with ties
  counted as in the first submission; `C` J_g from the in-sample AUC; `IW` the covariate-shift
  importance weights; `global`, `local`, and `mixed` the unweighted comparators.
- `J_cf`, `J_std`, `J_code`: J_g from the cross-validated AUC (method D), the in-sample AUC,
  and the in-sample AUC of the first submission.
- Table 2: `none`, `truncate`, `stabilize`, `calibrate` are the weight-handling methods;
  suffix `_cf` uses the cross-validated J_g (`stabilize_cf` = method D), `_code` the
  in-sample J_g of the first submission.

## Clinical application (SENIOR)

The SENIOR data are not public (see the paper's data availability statement). With access to
the data file, run

```sh
Rscript paper/run_senior.R <data.csv> <outdir>   # median errors, ESS, J_g per center
Rscript paper/make_fig4.R <data.csv> <outdir>    # Fig 4
```

Both scripts write every output, including per-patient errors and Fig 4, to `<outdir>`,
and stop if `<outdir>` lies inside this repository.

## Software versions

The scripts reproduce `results/` with R 4.4.0 (aarch64-apple-darwin20) and dplyr 1.1.4,
lme4 1.1-35.3 (Matrix 1.7-0), tidyr 1.3.1, ggplot2 3.5.1, gridExtra 2.3, ragg 1.3.2, and
readr 2.1.5; simcow imports ggplot2 and Hmisc 5.1-3. The parallel package ships with R.
