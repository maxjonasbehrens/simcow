#' Plot Cumulative Distribution Functions (CDFs)
#'
#' This function plots the CDFs for a given variable in a dataset, and showcases the distributions
#' of the target cohort, the overall cohort, and a weighted global cohort.
#'
#' @param data A dataframe containing the data to be visualized.
#' @param variable A character string indicating the column name of the continuous variable to visualize.
#' @param target_var A character string indicating the column name that differentiates cohorts.
#' @param target The value in the `target_var` column that defines the target cohort.
#' @param weight_var A character string indicating the column name for the weights.
#'
#' @return A `ggplot` object representing the CDFs.
#'
#' @details
#' The function produces three lines on the CDF plot:
#' \itemize{
#'   \item Target cohort CDF, colored in "#FF0155".
#'   \item Overall cohort CDF, colored in "#617989".
#'   \item Weighted global cohort CDF, colored in "#9447FB".
#' }
#' The function uses weighted empirical CDF calculations (`wtd.Ecdf`) and visualizes the results using `ggplot2`.
#'
#' @examples
#' # Generate a mock dataset
#' set.seed(123)
#' data <- data.frame(variable = rnorm(100),
#'                    target_var = sample(1:3, 100, replace = TRUE),
#'                    weight_var = runif(100))
#'
#' # Use the plot_cdfs function
#' p <- plot_cdfs(data = data,
#'                variable = "variable",
#'                target_var = "target_var",
#'                target = 2,
#'                weight_var = "weight_var")
#' print(p)
#'
#' @import ggplot2
#' @import Hmisc
#' @export
plot_cdfs <- function(data, variable, target_var, target, weight_var) {
  # Ensure the inputs are correct variable names in the dataframe
  if (!(variable %in% names(data)) | !(target_var %in% names(data)) | !(weight_var %in% names(data))) {
    stop("One or more of the provided variable names do not exist in the data.")
  }

  # Use the provided variable names to subset and weight the data
  subset_data <- data
  subset_data[[weight_var]] <- ifelse(data[[target_var]] == target,1,0)
  overall_data <- data[[variable]]
  weighted_data <- data[[variable]]

  # Calculate sample sizes
  subset_size <- nrow(subset(subset_data, subset = subset_data[[target_var]] == target))
  overall_size <- length(overall_data)
  weighted_size <- sum(data[[weight_var]])

  t = unique(wtd.Ecdf(overall_data)$x)
  t = c(min(t)-1,t)

  cdf_target = wtd.Ecdf(subset_data[[variable]], weights = subset_data[[weight_var]])$ecdf
  if (length(t)!=length(cdf_target)){
    cdf_target = c(0.0,cdf_target)
  }
  cdf_global = wtd.Ecdf(overall_data)$ecdf
  cdf_weighted = wtd.Ecdf(weighted_data, weights = data[[weight_var]])$ecdf

  # create ggplot object
  p <- ggplot() +
    geom_line(aes(x = t, y = cdf_target), color = "#FF0155") +
    geom_line(aes(x = t, y = cdf_global), color = "#617989") +
    geom_line(aes(x = t, y = cdf_weighted), color = "#9447FB") +
    labs(y = "(weighted) CDF") +
    theme_minimal() +
    theme(legend.position = "topleft") +
    scale_color_manual(values = c("#FF0155", "#617989", "#9447FB")) +
    guides(color = guide_legend(override.aes = list(linetype = c("solid", "solid", "solid")))) +
    labs(color = "Legend",
         title = "",
         subtitle = "",
         caption = paste("Target (n =", subset_size, ")",
                         "Global (n =", overall_size, ")",
                         "Weighted Global (n =", round(weighted_size, 2), ")"))

  # return the ggplot object
  return(p)
}
