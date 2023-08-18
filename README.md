# Similarity-based Cohort Weighting (simcow)

This R package provides functions to generate cohort weights for datasets with multiple cohorts. These weights are useful for adjusting datasets to make them more comparable across different cohorts.

## Installation

You can install this package from GitHub using the `devtools` package:

```R
# Install the devtools package if not already installed
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# Install the package
devtools::install_github("maxjonasbehrens/simcow")
```

## Functions

### 1. create_cow()

Generate cohort weights for datasets with multiple cohorts.

#### Arguments

- `data`: A dataframe containing the multi-cohort data.
- `cohort_id`: A character string indicating the column name which specifies cohort labels in the `data`.
- `target_cohort`: A character or numeric value specifying which cohort is the target.
- `baseline_covs`: A character vector containing the column names of baseline covariates used in computing propensity scores.
- `trunc_weights`: A logical value. If TRUE, the computed weights are truncated to 1. Default is TRUE.

#### Returns

A dataframe similar to the input `data` but with additional columns:
- `c_flag`: A binary flag indicating the current cohort versus the target cohort.
- `ps_weight`: The computed propensity score for each observation.
- `weight`: The computed cohort weight for each observation.

#### Example

```R
# Generate a mock dataset
set.seed(123)
data <- data.frame(cohort_id = sample(1:3, 100, replace = TRUE),
                   covariate1 = runif(100),
                   covariate2 = rnorm(100))

# Use the create_cow function
weighted_data <- create_cow(data = data,
                           cohort_id = "cohort_id",
                           target_cohort = 2,
                           baseline_covs = c("covariate1", "covariate2"))
```

### 2. plot_cdfs()

Plot Cumulative Distribution Functions (CDFs) for a given variable in a dataset.

#### Arguments

- `data`: A dataframe containing the data to be visualized.
- `variable`: A character string indicating the column name of the continuous variable to visualize.
- `target_var`: A character string indicating the column name that differentiates cohorts.
- `target`: The value in the `target_var` column that defines the target cohort.
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
               target_var = "target_var",
               target = 2,
               weight_var = "weight_var")
print(p)
```

---

You may want to modify the installation instructions to include the actual GitHub URL where you host your package. This README provides a good starting point for potential users of your R package to understand its purpose, installation process, and usage.
