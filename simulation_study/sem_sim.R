simulate_scenario <- function(type, target_environment, n_target, n = 30, n_cat = 4, num_global_predictors = 1, num_env_predictors = 1, from = -2, to = 2) {

  # Categorical distribution with n_cat categories, range controlled by from and to
  A <- sample(seq(from = from, to = to, length.out = n_cat), size = n, replace = TRUE, prob = rep(1/n_cat, n_cat))

  # Generate global predictors which don't depend on A
  epsilon_X_global <- matrix(rnorm(n * num_global_predictors), ncol = num_global_predictors)

  # Generate environment-specific predictors
  epsilon_X_env <- matrix(rnorm(n * num_env_predictors), ncol = num_env_predictors)

  # Normal noise terms
  epsilon_H <- matrix(rnorm(n, 0, 1), ncol = 1)
  epsilon_Y <- rnorm(n, 0, 1)

  # Equations based on the type of distribution shift
  switch(type,
         "AX" = {
           H <- epsilon_H

           X_env <- epsilon_X_env
           for (i in 1:num_env_predictors) {
             X_env[,i] <- X_env[,i] + A + H
           }

           X_global <- epsilon_X_global
           for (i in 1:num_global_predictors) {
             X_global[,i] <- X_global[,i] + epsilon_H
           }

           Y <- rowSums(X_global) + rowSums(X_env) + 2 * H + epsilon_Y
         },
         "AY" = {
           H <- epsilon_H

           X_env <- epsilon_X_env
           for (i in 1:num_env_predictors) {
             X_env[,i] <- X_env[,i] + H
           }

           X_global <- epsilon_X_global
           for (i in 1:num_global_predictors) {
             X_global[,i] <- X_global[,i] + epsilon_H
           }

           Y <- rowSums(X_global) + rowSums(X_env) + 2 * H + epsilon_Y + A
         },
         "AH" = {
           H <- epsilon_H + A

           X_env <- epsilon_X_env
           for (i in 1:num_env_predictors) {
             X_env[,i] <- X_env[,i] + H
           }

           X_global <- epsilon_X_global
           for (i in 1:num_global_predictors) {
             X_global[,i] <- X_global[,i] + epsilon_H
           }

           Y <- rowSums(X_global) + rowSums(X_env) + 2 * H + epsilon_Y
         }
  )

  # Combine into data frame
  A <- as.factor(A)
  levels(A) <- c(1:n_cat)
  colnames(X_env) <- paste0("X_env_", 1:num_env_predictors)
  colnames(X_global) <- paste0("X_global_", 1:num_global_predictors)
  data_df <- data.frame(A, X_global, X_env, H = H, Y)

  t_dat <- data_df |>
    filter(A == target_environment) |>
    head(n_target)

  o_dat <- data_df |>
    filter(A != target_environment)

  data_df <- rbind(t_dat, o_dat)

  return(data_df)
}

# test <- simulate_scenario(type = "AX",
#                           target_environment = 1,
#                           n_target = 100,
#                           n = 1000,
#                           n_cat = 4,
#                           num_global_predictors = 1,
#                           num_env_predictors = 4,
#                           to = 3,from = 0)
#
# test |>
#   ggplot(aes(X_env_1,Y,color=A))+
#   geom_point()+
#   geom_smooth(method = "lm")
#
# test |>
#   ggplot(aes(X_env_1,Y))+
#   geom_point()+
#   geom_smooth(method = "lm")
#
# test |>
#   ggplot(aes(Y, fill = A)) +
#   geom_density()
#
# summary(lm(Y~X_global_1, test, subset = A==1))
# summary(lm(Y~X_global_1, test))

# Produce results from simulation
evaluate_methods <- function(scenario_type, target_environment, n_target, num_simulations = 1000, n_train = 30, n_test = 30,
                             n_cat = 4, num_global_predictors = 1, num_env_predictors = 1, from = -2, to = 2, AUC = TRUE,
                             overlap_threshold = NULL) {

  rmse <- function(y_true, y_pred) {
    return(sqrt(mean((y_true - y_pred)^2)))
  }

  # Prediction performance
  rmse_global <- numeric(num_simulations)
  rmse_env_specific <- numeric(num_simulations)
  rmse_weighted_global <- numeric(num_simulations)
  ess_ratio <- numeric(num_simulations)

  # Weights as similarity measure
  subgroup_weights <- vector("list", num_simulations)
  subgroup_similarities <- vector("list", num_simulations)

  global_predictors <- paste("X_global_", 1:num_global_predictors, sep = "", collapse = " + ")
  env_predictors <- paste("X_env_", 1:num_env_predictors, sep = "", collapse = " + ")
  formula_str <- paste("Y ~", global_predictors, "+", env_predictors)

  for (i in 1:num_simulations) {
    # Generate training and test data
    train_data <- simulate_scenario(type = scenario_type,
                                    target_environment = target_environment,
                                    n_target = n_target, n = n_train, n_cat = n_cat,
                                    num_global_predictors = num_global_predictors,
                                    num_env_predictors = num_env_predictors,
                                    from = from, to = to)

    test_data <- simulate_scenario(type = scenario_type,
                                   target_environment = target_environment,
                                   n_target = n_test/n_cat, n = n_test, n_cat = n_cat,
                                   num_global_predictors = num_global_predictors,
                                   num_env_predictors = num_env_predictors,
                                   from = from, to = to)

    test_data_env <- subset(test_data, A == target_environment)

    # Create weighted train data
    train_w <- create_cow(data = train_data, cohort_id = "A",
                          target_cohort = target_environment,
                          baseline_covs = c(global_predictors, env_predictors,"Y"),
                          trunc_weights = FALSE,
                          AUC = AUC,
                          overlap_threshold = overlap_threshold)

    # Store weights
    subgroup_weights[[i]] <- train_w$weight

    # Compute and store similarity
    subgroup_similarities[[i]] <- (as.numeric(train_w$A) - min(as.numeric(train_w$A))) /
      (max(as.numeric(train_w$A)) - min(as.numeric(train_w$A))) * to

    # Create weighted test data
    test_data_w <- rbind(subset(test_data, A == target_environment),
                         subset(train_data, A != target_environment))

    test_data_w <- create_cow(data = test_data_w, cohort_id = "A",
                              target_cohort = target_environment,
                              baseline_covs = c(global_predictors, env_predictors,"Y"),
                              trunc_weights = FALSE,
                              AUC = AUC,
                              overlap_threshold = overlap_threshold)

    test_data_env_w <- subset(test_data_w, A == target_environment)

    # Sample size
    ess_weighted_global <- sum(train_w$weight)
    ess_env_specific <- nrow(subset(train_data, A == target_environment))
    ess_global <- nrow(train_data)

    # ESS ratio
    ess_ratio[i] <- ess_weighted_global/ess_env_specific

    # Global model
    model_global <- lm(formula_str, data = train_w)
    predictions_global <- predict(model_global, newdata = test_data_env)
    rmse_global[i] <- rmse(test_data_env$Y, predictions_global)

    # Weighted global model
    model_weighted_global <- lm(formula_str, data = train_w, weights = weight, x = TRUE)
    predictions_weighted_global <- predict(model_weighted_global, newdata = test_data_env_w)
    rmse_weighted_global[i] <- rmse(test_data_env_w$Y, predictions_weighted_global)

    # Environment-specific model
    train_data_env <- subset(train_data, A == target_environment)

    if(nrow(train_data_env) > 0) {
      model_env <- lm(formula_str, data = train_data_env)
      predictions_env <- predict(model_env, newdata = test_data_env)
      rmse_env_specific[i] <- rmse(test_data_env$Y, predictions_env)
    } else {
      rmse_env_specific[i] <- NA # NA for iterations where target environment data is missing
    }


  }

  close(pb)
  return(list(rmse_global = rmse_global,
              rmse_env_specific = rmse_env_specific,
              rmse_weighted_global = rmse_weighted_global,
              ess_ratio = ess_ratio,
              weights = subgroup_weights,
              similarities = subgroup_similarities))
}


### Simulation run
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)

# Parameter Grid for prediction performance
param_grid <- expand.grid(scenario = c("AX","AY","AH"),
                          n_cat = c(2,4,6,8),
                          n_train = round(seq(10,50,length.out = 10)),
                          n_target = c(10,12,15,20),
                          to = c(3,2,1),
                          num_env_predictors = c(3)
)

# Initialize a list to hold RMSE values
param_grid$results <- vector("list", nrow(param_grid))

# Progress
pb = txtProgressBar(min = 0, max = nrow(param_grid), initial = 0, style = 3)

# Loop through each row of the data frame
set.seed(123)
for (i in 1:nrow(param_grid)) {
  params <- param_grid[i, ]

  eval_result <- evaluate_methods(scenario_type = params$scenario,
                                  target_environment = 1,
                                  num_simulations = 100,
                                  num_env_predictors = params$num_env_predictors,
                                  num_global_predictors = 4-params$num_env_predictors,
                                  n_train = params$n_train*params$n_cat,
                                  n_test = 50*params$n_cat,
                                  n_target = params$n_target,
                                  n_cat = params$n_cat,
                                  from = 0,
                                  to = params$to,
                                  AUC = TRUE, # Set to false for result for p hat
                                  overlap_threshold = NULL)

  # Store RMSE values in a data frame and then store this data frame in the list
  rmse_df <- data.frame(Global = eval_result$rmse_global,
                        Target = eval_result$rmse_env_specific,
                        Weighted = eval_result$rmse_weighted_global)

  # Create an ID for each element in the nested lists
  similarity_ids <- rep(seq_along(eval_result$similarities),
                        sapply(eval_result$similarities, length))

  similarity_df <- data.frame(ID = similarity_ids,
                              scenario = params$scenario,
                              similarity = unlist(eval_result$similarities),
                              weights = unlist(eval_result$weights))

  param_grid$plotting[[i]] <- similarity_df
  param_grid$results[[i]] <- rmse_df
  param_grid$ess_ratio[i] <- mean(eval_result$ess_ratio)

  setTxtProgressBar(pb,i)
}

param_grid$SID <- rownames(param_grid)

# Unnest the results for easier plotting
param_grid_unnested <- param_grid %>%
  unnest(results) %>%
  gather(key = "Method", value = "RMSE", Global, Target, Weighted)

levels(param_grid_unnested$scenario) <- c("covariate","outcome","out. + cov.")
param_grid_unnested$scenario <- factor(param_grid_unnested$scenario, levels = c("outcome","covariate","out. + cov."))

param_grid_unnested$Method <- factor(param_grid_unnested$Method, levels = c("Target","Global","Weighted"))

library(ggsci)

# Boxplots
param_grid_unnested |>
  filter(ess_ratio < 4.5) |>
  mutate(to_fct = ifelse(to == 1,"similar",ifelse(to == 2,"medium","dissimilar"))) |>
  mutate(ess_ratio = as.factor(round(ess_ratio))) |>
  ggplot(aes(x = ess_ratio, y = RMSE, fill = as.factor(Method))) +
  geom_boxplot(outlier.shape = NA) +  # Hide outliers
  scale_fill_manual(values = c("#DDAA33","#BB5566","#004488"), labels = c("local","global","weighted")) +
  facet_grid(scenario~to_fct, scales = "free", space = "free") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.position = "bottom"
  ) +
  labs(
    # title = expression(paste(n["target"])),
    title = "Subgroup Similarity",
    x = expression(paste(ESS["weighted"] / ESS["local"])),
    y = "RMSE",
    fill = "Sample",
    alpha = NULL
  ) + coord_cartesian(ylim = c(0.7,3))

# Save plot
# ggsave("PLOTS/sim_ps/sem_sim_results_123.pdf", width = 12, height = 8, dpi= 700, device = "pdf")

# Results for Table 1
mean(unlist(param_grid_unnested[param_grid_unnested$Method=="Weighted","ess_ratio"]))

mean(unlist(param_grid_unnested[param_grid_unnested$Method=="Global","RMSE"]))
mean(unlist(param_grid_unnested[param_grid_unnested$Method=="Target","RMSE"]))
mean(unlist(param_grid_unnested[param_grid_unnested$Method=="Weighted","RMSE"]))

