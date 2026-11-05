# Clinical application of the paper (SENIOR cohort; Fig 4): prediction of the radiotherapy dose
# in five target centers from the local, global, and weighted samples.
# Resampling: for each target center, 1000 repetitions; 63.2 % of the target patients drawn
# with replacement for training, the out-of-bag target patients as test set (seed 42 per
# center). The weights (create_cow(), method D) are re-estimated in each repetition from the
# training target patients (each patient once) and all external patients.
# The data are not public. All outputs, including per-patient errors, go to <outdir>, which must
# lie outside this repository.
# Usage: Rscript paper/run_senior.R <data.csv> <outdir>
suppressPackageStartupMessages({ library(dplyr); library(readr); library(parallel) })
here <- normalizePath(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(here, "R", "methods.R"))         # kish_ess()
source(file.path(here, "R", "read_clinical.R"))
source(file.path(here, "..", "R", "create_cow.R"))
args <- commandArgs(TRUE)
if (length(args) != 2) stop("Usage: Rscript paper/run_senior.R <data.csv> <outdir>")
outdir <- outdir_outside_repo(args[2], file.path(here, ".."))
cores <- if (.Platform$OS.type == "windows") 1L else 5L

senior <- read_clinical(args[1])
covs <- c("BeginRadioTherapie", "CCI", "Hb", "Leukozyten"); outcome <- "DosisEQDW"
centers <- c(2, 3, 6, 7, 11)
tmp <- senior |> filter(Zentrum %in% centers)
f_y <- as.formula(paste0(outcome, "~", paste(covs, collapse = "+")))
sq_rmse <- function(y, p) sqrt(mean((y - p)^2))
gy <- function(x) sqrt(pmax(x, 0))   # back-transform from the squared dose
M <- c("local", "global", "weighted")

run_target <- function(target) {
  # Row order of the paper's run: external patients by center, then the target patients
  ext_rows <- which(tmp$Zentrum != target)
  data <- tmp[c(ext_rows[order(tmp$Zentrum[ext_rows])], which(tmp$Zentrum == target)), ]
  names(data)[names(data) == "Zentrum"] <- "cohort"
  data$ID <- seq_len(nrow(data))
  target_data <- data %>% filter(cohort == target)
  not_target <- data %>% filter(cohort != target)
  n0 <- nrow(target_data)
  set.seed(42)
  per_rep <- vector("list", 1000); per_obs <- vector("list", 1000); jg <- vector("list", 1000)
  for (rep in 1:1000) {
    train_t <- target_data %>% ungroup() %>% slice_sample(prop = 0.632, replace = TRUE)
    train_set <- rbind(train_t, not_target)
    test_set <- target_data %>% filter(!(ID %in% train_set$ID))
    est <- rbind(distinct(train_t, ID, .keep_all = TRUE), not_target)    # each training patient once
    cw <- create_cow(est, "cohort", target, covs, outcome = outcome, seed = 3e6 + 1000 * target + rep)
    ext <- train_set$cohort != target
    w <- rep(1, nrow(train_set)); w[ext] <- cw$weight[est$cohort != target]
    fits <- list(local = lm(f_y, data = train_set[train_set$cohort == target, ]),
                 global = lm(f_y, data = train_set),
                 weighted = { d <- train_set; d$.w <- w; lm(f_y, data = d, weights = .w) })
    pred <- sapply(fits, predict, newdata = test_set)
    if (is.null(dim(pred))) pred <- matrix(pred, nrow = 1, dimnames = list(NULL, M))
    y <- test_set[[outcome]]
    per_rep[[rep]] <- data.frame(rep = rep, model = M,
      err_sq = apply(pred, 2, function(p) sq_rmse(y, p)),          # squared-dose scale
      rmse_gy = apply(pred, 2, function(p) sq_rmse(gy(y), gy(p))),
      ess = c(NA, NA, n0 + kish_ess(w[ext])))
    per_obs[[rep]] <- data.frame(ID = rep(test_set$ID, length(M)), model = rep(M, each = nrow(test_set)),
      abs_sq = as.vector(abs(y - pred)), abs_gy = as.vector(abs(gy(y) - gy(pred))))
    jg[[rep]] <- data.frame(rep = rep, attr(cw, "subgroups"))
  }
  per_rep <- do.call(rbind, per_rep); per_obs <- do.call(rbind, per_obs); jg <- do.call(rbind, jg)
  obs <- per_obs |> group_by(model, ID) |> summarise(AVG = mean(abs_sq), mae_gy = mean(abs_gy), .groups = "drop")
  write.csv(cbind(center = target, obs), file.path(outdir, paste0("senior_obs_errors_center", target, ".csv")), row.names = FALSE)
  metrics <- per_rep |> group_by(model) |>
    summarise(RMSE_code = median(sqrt(err_sq)), RMSE_gy_median = median(rmse_gy),
              RMSE_gy_q25 = quantile(rmse_gy, .25), RMSE_gy_q75 = quantile(rmse_gy, .75),
              ESS_median = median(ess), .groups = "drop") |>
    left_join(obs |> group_by(model) |>
                summarise(MAE_code = median(sqrt(AVG)), MAE_gy_median = median(mae_gy),
                          MAE_gy_q25 = quantile(mae_gy, .25), MAE_gy_q75 = quantile(mae_gy, .75), .groups = "drop"),
              by = "model")
  list(metrics = cbind(center = target, n0 = n0, metrics),
       jg = jg |> group_by(cohort) |> summarise(J_median = median(J), J_q25 = quantile(J, .25),
                                                J_q75 = quantile(J, .75), J_zero = mean(J == 0), .groups = "drop") |>
         mutate(target = target, .before = 1))
}

res <- mclapply(centers, run_target, mc.cores = cores)
errs <- Filter(function(r) inherits(r, "try-error"), res)
if (length(errs)) { print(errs[[1]]); quit(status = 1) }
metrics <- do.call(rbind, lapply(res, `[[`, "metrics"))
jg <- do.call(rbind, lapply(res, `[[`, "jg"))
write.csv(metrics, file.path(outdir, "senior_metrics.csv"), row.names = FALSE)
write.csv(jg, file.path(outdir, "senior_jg.csv"), row.names = FALSE)

cat("\nMedian MAE in Gy (median over patients of the mean absolute error across resamples):\n")
print(tidyr::pivot_wider(metrics[, c("center", "n0", "model", "MAE_gy_median")], names_from = model, values_from = MAE_gy_median), digits = 3)
cat("\nMedian RMSE in Gy across resamples:\n")
print(tidyr::pivot_wider(metrics[, c("center", "model", "RMSE_gy_median")], names_from = model, values_from = RMSE_gy_median), digits = 3)
cat("\nMedian ESS (n0 + Kish ESS of the external weights) and ratio to n0:\n")
print(metrics |> filter(model == "weighted") |> transmute(center, ESS_median, ratio = ESS_median / n0), digits = 3)
cat("\nCross-validated J_g per pair (median over resamples):\n")
print(as.data.frame(jg), digits = 3)
cat("\nWritten to", outdir, "\n")
