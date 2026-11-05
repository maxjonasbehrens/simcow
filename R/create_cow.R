#' Create Cohort Weights (COW) for Multi-Cohort Data
#'
#' This function computes similarity-based weights for borrowing observations from
#' external subgroups (cohorts) when a prediction model is fitted for a target subgroup.
#' Target observations get weight 1. Each external observation gets a weight that combines
#' an individual similarity (how much the observation resembles the target subgroup) and a
#' subgroup similarity (how much its subgroup overlaps with the target subgroup).
#'
#' @param data A dataframe containing the multi-cohort data. The columns in `cohort_id`,
#'   `baseline_covs`, and `outcome` must not contain missing values.
#' @param cohort_id A character string indicating the column name which specifies cohort labels in the `data`.
#' @param target_cohort A character or numeric value specifying which cohort is the target.
#' @param baseline_covs A character vector containing the column names of the covariates used in the membership models.
#' @param outcome A character string with the column name of the outcome. The outcome enters the
#'   membership models together with `baseline_covs`, so that the weights reflect differences in
#'   the outcome given the covariates. This gives the weights of Behrens et al. If NULL (default),
#'   the membership models use `baseline_covs` only.
#' @param K_folds Number of cross-validation folds for the subgroup AUC. Default is 5.
#' @param seed An integer or NULL. If set, the random fold assignment uses this seed and the
#'   global random number stream is left unchanged. If NULL (default), the folds are drawn from
#'   the current random number stream.
#'
#' @return The input `data` in its original row order with three additional columns:
#' \itemize{
#'   \item \strong{c_flag}: 1 for observations of the target cohort, 0 for external observations.
#'   \item \strong{ps_weight}: The fitted probability of target membership of an external observation
#'     (from the membership model of its cohort); NA for target observations.
#'   \item \strong{weight}: The cohort weight; 1 for target observations.
#' }
#' The attribute \code{"subgroups"} holds one row per external cohort with its size \code{n},
#' the cross-validated \code{auc}, and the subgroup weight \code{J}.
#'
#' @details For each external cohort g, a logistic membership model discriminates the target
#' cohort from cohort g, using the covariates in `baseline_covs` and the outcome in `outcome`.
#' The model is fitted on all observations of the two cohorts; \eqn{\hat p_{jg}} is its fitted
#' probability that external observation j belongs to the target cohort. The weight of j is
#' \deqn{w_j = J_g \, \hat p_{jg} / \bar p_g, \qquad J_g = 1 - |2\,\mathrm{AUC}_g - 1|,}
#' where \eqn{\bar p_g} is the mean of \eqn{\hat p_{jg}} within cohort g, so the individual
#' weights in cohort g average to \eqn{J_g}. \eqn{\mathrm{AUC}_g} is the two-class AUC of the
#' membership predictions from `K_folds`-fold cross-validation (ties count one half).
#' \eqn{J_g} is 1 if the membership model cannot tell the two cohorts apart (AUC 0.5) and 0 if it
#' separates them completely (AUC 0 or 1).
#'
#' Version 0.1.0 computed \eqn{\min(1, \hat p_{jg} / \mathrm{AUC}_g)} instead, with the AUC from the
#' fitted (not cross-validated) predictions; see NEWS.md.
#'
#' @references Behrens M, et al. Improving prediction models by incorporating external data
#'   with weights based on similarity. Under revision at PLOS ONE.
#'
#' @examples
#' # Generate a mock dataset
#' set.seed(123)
#' data <- data.frame(cohort_id = sample(1:3, 100, replace = TRUE),
#'                    covariate1 = runif(100),
#'                    covariate2 = rnorm(100))
#' data$outcome <- data$covariate1 + data$covariate2 + rnorm(100)
#'
#' # Use the create_cow function
#' weighted_data <- create_cow(data = data,
#'                            cohort_id = "cohort_id",
#'                            target_cohort = 2,
#'                            baseline_covs = c("covariate1", "covariate2"),
#'                            outcome = "outcome",
#'                            seed = 1)
#' attr(weighted_data, "subgroups")
#'
#' @importFrom stats glm predict binomial as.formula
#' @export
create_cow <- function(data,
                       cohort_id,
                       target_cohort,
                       baseline_covs,
                       outcome = NULL,
                       K_folds = 5,
                       seed = NULL){

  covs <- c(baseline_covs, outcome)
  missing_cols <- setdiff(c(cohort_id, covs), names(data))
  if (length(missing_cols) > 0) {
    stop("Columns not found in data: ", paste(missing_cols, collapse = ", "))
  }
  if (anyNA(data[, c(cohort_id, covs)])) {
    stop("Missing values in the cohort, covariate, or outcome columns; remove or impute them first.")
  }
  groups <- as.character(data[[cohort_id]])
  target <- as.character(target_cohort)
  if (!target %in% groups) stop("target_cohort not found in column ", cohort_id)

  # Two-class AUC: probability that a target observation scores higher than an external one
  auc_two_class <- function(p, y) {
    r <- rank(p); n1 <- sum(y == 1); n0 <- sum(y == 0)
    (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  }

  # With a seed, draw the folds from their own stream and restore the global one afterwards
  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
    on.exit(if (is.null(old_seed)) rm(".Random.seed", envir = globalenv())
            else assign(".Random.seed", old_seed, envir = globalenv()))
    set.seed(seed)
  }

  ps_formula <- as.formula(paste("c_flag ~", paste(covs, collapse = " + ")))
  weight <- rep(1, nrow(data))
  ps_weight <- rep(NA_real_, nrow(data))
  cohorts <- setdiff(unique(groups), target)
  subgroups <- data.frame(cohort = cohorts, n = NA_integer_, auc = NA_real_, J = NA_real_)

  for (i in seq_along(cohorts)) {
    g <- cohorts[i]
    pair <- groups %in% c(target, g)
    pair_data <- as.data.frame(data[pair, covs, drop = FALSE])
    pair_data$c_flag <- ifelse(groups[pair] == g, 0, 1)

    # Membership probabilities from the fit on all observations of the pair
    fit <- suppressWarnings(glm(ps_formula, data = pair_data, family = binomial))
    p <- predict(fit, type = "response")

    # AUC from cross-validated membership predictions
    folds <- sample(rep(seq_len(K_folds), length.out = nrow(pair_data)))
    p_cv <- numeric(nrow(pair_data))
    for (f in seq_len(K_folds)) {
      if (any(folds == f)) {
        fit_f <- suppressWarnings(glm(ps_formula, data = pair_data[folds != f, , drop = FALSE], family = binomial))
        p_cv[folds == f] <- predict(fit_f, newdata = pair_data[folds == f, , drop = FALSE], type = "response")
      }
    }
    auc <- auc_two_class(p_cv, pair_data$c_flag)
    J <- 1 - abs(2 * auc - 1)

    p_ext <- p[pair_data$c_flag == 0]
    idx <- which(groups == g)
    ps_weight[idx] <- p_ext
    weight[idx] <- J * p_ext / mean(p_ext)
    subgroups[i, c("n", "auc", "J")] <- list(length(idx), auc, J)
  }

  data$c_flag <- ifelse(groups == target, 1, 0)
  data$ps_weight <- ps_weight
  data$weight <- weight
  attr(data, "subgroups") <- subgroups
  return(data)
}
