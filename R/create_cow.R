#' Create Cohort Weights (COW) for Multi-Cohort Data
#'
#' This function generates cohort weights for datasets with multiple cohorts. It utilizes propensity scores and the C-index to compute the weights.
#'
#' @param data A dataframe containing the multi-cohort data.
#' @param cohort_id A character string indicating the column name which specifies cohort labels in the `data`.
#' @param target_cohort A character or numeric value specifying which cohort is the target.
#' @param baseline_covs A character vector containing the column names of baseline covariates used in computing propensity scores.
#' @param trunc_weights A logical value. If TRUE, the computed weights are truncated to 1. Default is TRUE.
#'
#' @return A dataframe similar to the input `data` but with additional columns:
#' \itemize{
#'   \item \strong{c_flag}: A binary flag indicating the current cohort versus the target cohort.
#'   \item \strong{ps_weight}: The computed propensity score for each observation.
#'   \item \strong{weight}: The computed cohort weight for each observation.
#' }
#'
#' @details The function begins by identifying the cohorts in the `data` and renaming the cohort column to "cohort".
#' It then computes the propensity scores for each cohort relative to the target cohort.
#' Following that, the function computes the inverse C-index using the propensity scores and cohort flags.
#' Finally, the cohort weights are calculated and either truncated or not based on the `trunc_weights` parameter.
#'
#' @examples
#' # Generate a mock dataset
#' set.seed(123)
#' data <- data.frame(cohort_id = sample(1:3, 100, replace = TRUE),
#'                    covariate1 = runif(100),
#'                    covariate2 = rnorm(100))
#'
#' # Use the create_cow function
#' weighted_data <- create_cow(data = data,
#'                            cohort_id = "cohort_id",
#'                            target_cohort = 2,
#'                            baseline_covs = c("covariate1", "covariate2"))
#'
#' @importFrom stats glm predict
#' @export
create_cow <- function(data,
                       cohort_id,
                       target_cohort,
                       baseline_covs,
                       trunc_weights = TRUE){

  # C Index functions
  # Function to calculate the number of concordant and discordant pairs
  calculate_concordance_pairs <- function(predicted, actual) {
    n <- length(predicted)
    concordant <- 0
    discordant <- 0
    ties <- 0

    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        if (actual[i] < actual[j]) {
          if (predicted[i] < predicted[j]) {
            concordant <- concordant + 1
          } else if (predicted[i] > predicted[j]) {
            discordant <- discordant + 1
          } else {
            ties <- ties + 1
          }
        } else if (actual[i] > actual[j]) {
          if (predicted[i] > predicted[j]) {
            concordant <- concordant + 1
          } else if (predicted[i] < predicted[j]) {
            discordant <- discordant + 1
          } else {
            ties <- ties + 1
          }
        } else {
          # Ties in actual values
          if (predicted[i] == predicted[j]) {
            ties <- ties + 1
          }
        }
      }
    }

    return(list(concordant = concordant, discordant = discordant, ties = ties))
  }

  # Function to calculate C-index
  calculate_c_index <- function(predicted, actual) {
    pairs <- calculate_concordance_pairs(predicted, actual)
    concordant_pairs <- pairs$concordant
    discordant_pairs <- pairs$discordant
    ties <- pairs$ties

    total_pairs <- concordant_pairs + discordant_pairs + ties

    c_index <- (concordant_pairs + 0.5 * ties) / total_pairs
    return(c_index)
  }


  # Cohorts present in the data
  cohorts <- unique(data[,cohort_id])
  cohorts <- cohorts[cohorts!=target_cohort]

  # Empty list to store data combinations
  data_list <- list()

  # Formula for glm
  ps_formula <- as.formula(paste("c_flag ~ ",paste(baseline_covs,collapse = "+")))

  # Compute PS and IPW for each cohort
  for (c in cohorts){

    data_list[[c]] = data[data[, cohort_id] %in% c(target_cohort, c), ]
    data_list[[c]]$c_flag = ifelse(data_list[[c]][, cohort_id] == c, 0, 1)

    ps.fit.temp = glm(ps_formula, data = data_list[[c]], family = binomial(link = "logit"))
    ps.temp  = predict(ps.fit.temp, type = "response")

    # Propensity scores
    data_list[[c]]$ps_weight = ps.temp

    # Data with only the target site included
    target_dat = data_list[[c]][data_list[[c]][, cohort_id] == target_cohort, ]
    target_dat$weight = 1

    # Compute inverse c index
    inverse_c <- 1/calculate_c_index(predicted = data_list[[c]]$ps_weight,
                                     actual = data_list[[c]]$c_flag)

    # Adjust weights
    if (trunc_weights){
      data_list[[c]]$weight <- ifelse(data_list[[c]]$ps_weight * inverse_c > 1, 1, data_list[[c]]$ps_weight * inverse_c)
    } else {
      data_list[[c]]$weight <- data_list[[c]]$ps_weight * inverse_c
    }
    data_list[[c]] <- data_list[[c]][data_list[[c]][, cohort_id] != target_cohort,]


  }

  # Put everything together
  final_dat <- do.call(rbind.data.frame, data_list)
  final_dat <- rbind(final_dat, target_dat)

  return(final_dat)
}
