# Data-generating process of the simulation study (three types of distribution shift).
# Copied verbatim from the analysis code of the paper, so that the random number stream
# and with it every simulated data set match the paper. Uses dplyr::filter().
suppressPackageStartupMessages(library(dplyr))

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
           coeff_env <- runif(ncol(X_env),min = -1, max = 1)
           # Y <- rowSums(X_global) + as.vector(X_env %*% coeff_env)  + 2 * H + epsilon_Y
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
           coeff_env <- runif(ncol(X_env),min = -1, max = 1)
           # Y <- rowSums(X_global) + as.vector(X_env %*% coeff_env)  + 2 * H + epsilon_Y + A
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
           coeff_env <- runif(ncol(X_env),min = -1, max = 1)
           # Y <- rowSums(X_global) + as.vector(X_env %*% coeff_env)  + 2 * H + epsilon_Y
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
