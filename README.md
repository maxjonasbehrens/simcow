<img src='logo/simcow_logo.png' width='250'>

# Sim<sub>ilarity-based</sub> Co<sub>hort</sub> W<sub>eighting</sub> (simcow)

This R package provides functions to generate cohort weights for datasets with multiple cohorts. The weights let a prediction model for a target cohort borrow observations from external cohorts in proportion to their similarity to the target cohort. Version 0.2.0 implements the method of the paper below; `paper/` holds the code that reproduces its results.

## Citation

Behrens M, et al. Improving prediction models by incorporating external data with weights based on similarity. Under revision at PLOS ONE.

## Installation

You can install this package from GitHub using the `remotes` package:

```R
# Install the remotes package if not already installed
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

# Install the package (version 0.2.0)
remotes::install_github("maxjonasbehrens/simcow@v0.2.0")
```

## Functions

### 1. create_cow()

Generate cohort weights for datasets with multiple cohorts.

Target observations get weight 1. For each external cohort g, a logistic membership model discriminates the target cohort from cohort g, using the covariates and the outcome. External observation j gets the weight

w<sub>j</sub> = J<sub>g</sub> &middot; p&#770;<sub>jg</sub> / p&#772;<sub>g</sub>, &nbsp; J<sub>g</sub> = 1 &minus; |2 AUC<sub>g</sub> &minus; 1|,

where p&#770;<sub>jg</sub> is the fitted probability that j belongs to the target cohort, p&#772;<sub>g</sub> its mean within cohort g, and AUC<sub>g</sub> the AUC of 5-fold cross-validated membership predictions. J<sub>g</sub> is 1 if the membership model cannot tell the two cohorts apart and 0 if it separates them completely.

#### Arguments

- `data`: A dataframe containing the multi-cohort data (no missing values in the columns used).
- `cohort_id`: A character string indicating the column name which specifies cohort labels in the `data`.
- `target_cohort`: A character or numeric value specifying which cohort is the target.
- `baseline_covs`: A character vector containing the column names of the covariates used in the membership models.
- `outcome`: The column name of the outcome; it enters the membership models together with `baseline_covs`, as in the paper. Default is NULL (covariates only).
- `K_folds`: Number of cross-validation folds for AUC<sub>g</sub>. Default is 5.
- `seed`: If set, seed for the fold assignment; the global random number stream stays unchanged. Default is NULL.

#### Returns

The input `data` in its original row order with additional columns:
- `c_flag`: 1 for the target cohort, 0 for external cohorts.
- `ps_weight`: The fitted probability of target membership of an external observation (NA for target observations).
- `weight`: The computed cohort weight for each observation.

The attribute `"subgroups"` gives the size `n`, the cross-validated `auc`, and the subgroup weight `J` of each external cohort.

Version 0.1.0 computed different weights; see [NEWS.md](NEWS.md).

#### Example

```R
# Generate a mock dataset
set.seed(123)
data <- data.frame(cohort_id = sample(1:3, 100, replace = TRUE),
                   covariate1 = runif(100),
                   covariate2 = rnorm(100))
data$outcome <- data$covariate1 + data$covariate2 + rnorm(100)

# Use the create_cow function
weighted_data <- create_cow(data = data,
                           cohort_id = "cohort_id",
                           target_cohort = 2,
                           baseline_covs = c("covariate1", "covariate2"),
                           outcome = "outcome",
                           seed = 1)
attr(weighted_data, "subgroups")

# Fit the prediction model for the target cohort on the weighted data
fit <- lm(outcome ~ covariate1 + covariate2, data = weighted_data, weights = weight)
```

### 2. plot_cdfs()

Plot Cumulative Distribution Functions (CDFs) for a given variable in a dataset.

#### Arguments

- `data`: A dataframe containing the data to be visualized.
- `variable`: A character string indicating the column name of the continuous variable to visualize.
- `cohort_id`: A character string indicating the column name that differentiates cohorts.
- `target_cohort`: The value in the `cohort_id` column that defines the target cohort.
- `weight_var`: A character string indicating the column name for the weights.

#### Returns

A `ggplot` object representing the CDFs.

#### Details

The function produces three lines on the CDF plot:
- Target cohort CDF, colored in "#FF0155".
- Overall cohort CDF, colored in "#617989".
- Weighted global cohort CDF, colored in "#9447FB".

#### Example

```R
# Generate a mock dataset
set.seed(123)
data <- data.frame(variable = rnorm(100),
                   target_var = sample(1:3, 100, replace = TRUE),
                   weight_var = runif(100))

# Use the plot_cdfs function
p <- plot_cdfs(data = data,
               variable = "variable",
               cohort_id = "target_var",
               target_cohort = 2,
               weight_var = "weight_var")
print(p)
```

## Reproducing the paper

See [paper/README.md](paper/README.md). The folder `simulation_study/` holds an earlier version of the simulation (2024) and does not reproduce the paper's results.
