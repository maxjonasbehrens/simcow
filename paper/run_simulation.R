# Simulation study of the paper: main grid (Table 1, Figs 2 and 3) and bias-variance
# decomposition. The weighted sample uses simcow's create_cow() (method D); variants A-C and
# the covariate-shift comparator come from R/methods.R.
# Main grid: 3 scenarios x 4 numbers of subgroups x 3 n_train x 3 n_target x 3 similarity
# levels, 100 data sets each (32,400 runs). The data sets replay the random number stream of
# the paper's run (set.seed(123); the grid also lists the four weight-handling methods of
# Table 2, whose data sets are drawn in the same stream).
# Usage (from any directory): Rscript paper/run_simulation.R
# Writes paper/output/sim_runs.csv, sim_subgroups.csv, bias_variance.csv, and the replayed data
# sets to paper/output/cache/ (about 600 MB; needed by run_table2.R and make_figures.R). Run time
# under 20 min on 6 cores; stops if the results differ from paper/results/.
suppressPackageStartupMessages({ library(dplyr); library(lme4); library(parallel) })
here <- normalizePath(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(here, "R", "methods.R"))
source(file.path(here, "R", "simulate_scenario.R"))
source(file.path(here, "..", "R", "create_cow.R"))
outdir <- file.path(here, "output"); cachedir <- file.path(outdir, "cache")
dir.create(cachedir, showWarnings = FALSE, recursive = TRUE)
cores <- if (.Platform$OS.type == "windows") 1L else min(6L, detectCores(), na.rm = TRUE)
t0 <- Sys.time()

rmse <- function(y, yhat) sqrt(mean((y - yhat)^2))
fit_mixed <- function(d, f_y, f_mixed) {   # lm fallback as in the paper; NULL after an error
  gc <- table(d$A)
  tryCatch({
    if (min(gc) >= 5 && length(gc) >= 3) {
      mm <- suppressMessages(suppressWarnings(lmer(f_mixed, data = d, control = lmerControl(
        optimizer = "bobyqa", optCtrl = list(maxfun = 2e5), check.conv.singular = "ignore"))))
      if (isSingular(mm) || !mm@optinfo$conv$opt == 0) list(m = lm(f_y, data = d), fb = 1) else list(m = mm, fb = 0)
    } else list(m = lm(f_y, data = d), fb = 1)
  }, error = function(e) NULL)
}
pred_mixed <- function(fm, te) if (inherits(fm$m, "merMod")) predict(fm$m, newdata = te, allow.new.levels = TRUE) else predict(fm$m, newdata = te)
fit_weighted <- function(d, f_y, w) { d$.w <- w; lm(f_y, data = d, weights = .w) }

# ---------------- Main grid ----------------
param_grid <- expand.grid(scenario = c("AX","AY","AH"), n_cat = c(2,4,8,10), n_train = c(10,20,30),
                          n_target = c(10,20,30), to = c(3,2,1), num_env_predictors = c(3),
                          stabilization_method = c("none","truncate","stabilize","calibrate"),
                          stringsAsFactors = FALSE)
stab_rows <- which(param_grid$stabilization_method == "stabilize")
none_rows <- which(param_grid$stabilization_method == "none")
cache <- file.path(cachedir, "replay_kept.rds")             # data sets of the main grid
cache_none <- file.path(cachedir, "replay_kept_none.rds")   # data sets of the "none" row of Table 2
if (file.exists(cache) && file.exists(cache_none)) kept <- readRDS(cache) else {
  sim <- function(p, n, n_target) simulate_scenario(type = p$scenario, target_environment = 1,
    n_target = n_target, n = n, n_cat = p$n_cat, num_global_predictors = 1, num_env_predictors = 3,
    from = 0, to = p$to)
  set.seed(123); kept <- list(); kept_none <- list()
  for (i in seq_len(max(stab_rows))) { p <- param_grid[i, ]
    for (s in 1:100) { tr <- sim(p, p$n_train * p$n_cat, p$n_target); te <- sim(p, 50 * p$n_cat, 50)
      x <- list(i = i, s = s, tr = tr, te = subset(te, A == 1))
      if (i %in% stab_rows) kept[[length(kept) + 1]] <- x
      if (i %in% none_rows) kept_none[[length(kept_none) + 1]] <- x } }
  saveRDS(kept, cache, compress = FALSE); saveRDS(kept_none, cache_none, compress = FALSE)
  cat("Data sets replayed:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
}

covs <- c("X_global_1", paste0("X_env_", 1:3), "Y")
f_y <- Y ~ X_global_1 + X_env_1 + X_env_2 + X_env_3
f_mixed <- Y ~ X_global_1 + X_env_1 + X_env_2 + X_env_3 + (1 | A)
V <- c("A", "B", "C", "D")

one <- function(k) {
  x <- kept[[k]]; d <- x$tr; te <- x$te; p <- param_grid[x$i, ]
  mw <- membership_weights(d, "A", 1, covs, seed = 1e6 + k)
  cw <- create_cow(d, "A", 1, covs[1:4], outcome = "Y", seed = 1e6 + k)
  if (max(abs(cw$weight - mw$w[, "D"])) > 1e-12) stop("create_cow() differs from variant D in run ", k)
  W <- mw$w; W[, "D"] <- cw$weight
  w_iw <- iw_weights(d, "A", 1, covs[1:4])
  n0 <- sum(d$A == "1")
  r <- sapply(V, function(v) rmse(te$Y, predict(fit_weighted(d, f_y, W[, v]), newdata = te)))
  fm <- fit_mixed(d, f_y, f_mixed)
  info <- mw$info
  info$shift <- (as.numeric(info$group) - 1) / (p$n_cat - 1) * p$to
  list(run = data.frame(row = x$i, sim = x$s, p[, c("scenario","n_cat","n_train","n_target","to")], n0 = n0,
         t(setNames(r, paste0("rmse_", V))),
         rmse_IW = rmse(te$Y, predict(fit_weighted(d, f_y, w_iw), newdata = te)),
         rmse_global = rmse(te$Y, predict(lm(f_y, data = d), newdata = te)),
         rmse_local = if (n0 > 0) rmse(te$Y, predict(lm(f_y, data = d[d$A == "1", ]), newdata = te)) else NA,
         rmse_mixed = if (is.null(fm)) 0 else rmse(te$Y, pred_mixed(fm, te)),   # the paper's code leaves 0 after an error
         mixed_fallback = if (is.null(fm)) NA else fm$fb,
         t(setNames(sapply(V, function(v) kish_ess(W[, v]) / n0), paste0("ess_", V))),
         ess_IW = kish_ess(w_iw) / n0),
       info = cbind(row = x$i, sim = x$s, scenario = p$scenario, n_target = p$n_target, info))
}
res <- mclapply(seq_along(kept), one, mc.cores = cores)
errs <- Filter(function(r) inherits(r, "try-error"), res)
if (length(errs)) { cat("Errors in", length(errs), "jobs; first:\n"); print(errs[[1]]); quit(status = 1) }
runs <- do.call(rbind, lapply(res, `[[`, "run"))
info <- do.call(rbind, lapply(res, `[[`, "info"))
write.csv(runs, file.path(outdir, "sim_runs.csv"), row.names = FALSE)
write.csv(info, file.path(outdir, "sim_subgroups.csv"), row.names = FALSE)
cat("Main grid:", nrow(runs), "runs,", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")

# Checks: mean RMSE of the methods unchanged since the first submission, and every value
# equal to the results shipped in paper/results/
max_diff <- function(new, shipped, by) {
  m <- merge(shipped, new, by = by, suffixes = c(".shipped", ""))
  if (nrow(m) != nrow(new) || nrow(m) != nrow(shipped)) return(Inf)
  num <- setdiff(names(new)[sapply(new, is.numeric)], by)
  max(sapply(num, function(v) { a <- m[[paste0(v, ".shipped")]]; b <- m[[v]]
    if (!identical(is.na(a), is.na(b))) Inf else max(c(0, abs(a - b)), na.rm = TRUE) }))
}
chk <- c(A = mean(runs$rmse_A) - 1.570686332376494,
         IW = mean(runs$rmse_IW) - 1.8113088131950021,
         global = mean(runs$rmse_global) - 1.6500943606489176,
         local = mean(runs$rmse_local, na.rm = TRUE) - 2.11557515410958,
         sim_runs = max_diff(runs, read.csv(file.path(here, "results", "sim_runs.csv")), c("row", "sim")),
         sim_subgroups = max_diff(info, read.csv(file.path(here, "results", "sim_subgroups.csv")), c("row", "sim", "group")))
cat("Checks (difference to the first submission / to paper/results):\n"); print(signif(chk, 3))
if (any(abs(chk) > 1e-9)) { cat("CHECK FAILED\n"); quit(status = 1) }

# ---------------- Bias-variance decomposition ----------------
bv_cov <- c("X_global_1", "X_global_2", "X_env_1", "X_env_2", "Y")
f_y_bv <- Y ~ X_global_1 + X_global_2 + X_env_1 + X_env_2
f_mixed_bv <- Y ~ X_global_1 + X_global_2 + X_env_1 + X_env_2 + (1 | A)
bv_of <- function(P, y) { P <- P[, colSums(is.na(P)) < nrow(P), drop = FALSE]
  c(bias_sq = mean((rowMeans(P, na.rm = TRUE) - y)^2), variance = mean(apply(P, 1, var, na.rm = TRUE))) }
set.seed(123)
bv <- list()
for (si in 1:3) {
  scn <- c("AX", "AY", "AH")[si]
  sim_bv <- function(n_target, n) simulate_scenario(type = scn, target_environment = 1, n_target = n_target,
    n = n, n_cat = 4, num_global_predictors = 2, num_env_predictors = 2, from = 0, to = 2)
  te <- subset(sim_bv(50 / 4, 50), A == 1)
  P <- setNames(replicate(8, matrix(NA, nrow(te), 50), simplify = FALSE),
                c("global","env_specific","weighted_A","weighted_D","weighted_B","weighted_C","covariate_shift","mixed_model"))
  for (k in 1:50) {
    d <- sim_bv(10, 100)
    mw <- membership_weights(d, "A", 1, bv_cov, seed = 2e6 + 100 * si + k)
    mw$w[, "D"] <- create_cow(d, "A", 1, bv_cov[1:4], outcome = "Y", seed = 2e6 + 100 * si + k)$weight
    P$global[, k] <- predict(lm(f_y_bv, data = d), newdata = te)
    if (any(d$A == "1")) P$env_specific[, k] <- predict(lm(f_y_bv, data = d[d$A == "1", ]), newdata = te)
    for (v in V) P[[paste0("weighted_", v)]][, k] <- predict(fit_weighted(d, f_y_bv, mw$w[, v]), newdata = te)
    P$covariate_shift[, k] <- predict(fit_weighted(d, f_y_bv, iw_weights(d, "A", 1, bv_cov[1:4])), newdata = te)
    fm <- fit_mixed(d, f_y_bv, f_mixed_bv)
    if (!is.null(fm)) P$mixed_model[, k] <- pred_mixed(fm, te)
  }
  bv[[scn]] <- data.frame(scenario = scn, model = names(P), t(sapply(P, bv_of, y = te$Y)))
}
bv <- do.call(rbind, bv); bv$mse <- bv$bias_sq + bv$variance
write.csv(bv, file.path(outdir, "bias_variance.csv"), row.names = FALSE)
bv_chk <- max_diff(bv, read.csv(file.path(here, "results", "bias_variance.csv")), c("scenario", "model"))
cat("\nBias-variance check vs paper/results (max abs diff over", nrow(bv), "rows):", signif(bv_chk, 3), "\n")
if (bv_chk > 1e-9) { cat("BV CHECK FAILED\n"); quit(status = 1) }

# ---------------- Summaries ----------------
lab_s <- c(AX = "covariate", AY = "outcome", AH = "covariate + outcome")
lab_t <- c(`3` = "dissimilar", `2` = "medium", `1` = "similar")
cols <- c(weighted = "rmse_D", without_Jg = "rmse_A", Jg_in_sample_code_auc = "rmse_B",
          Jg_in_sample_auc = "rmse_C", covariate_shift_IW = "rmse_IW", global = "rmse_global",
          local = "rmse_local", mixed = "rmse_mixed")
cat("\nOverall RMSE (Table 1; 32,400 runs):\n")
print(round(sapply(cols, function(v) c(mean = mean(runs[[v]], na.rm = TRUE), sd = sd(runs[[v]], na.rm = TRUE),
                                       median = median(runs[[v]], na.rm = TRUE))), 4))
cat("\nMean ESS ratio: weighted", round(mean(runs$ess_D), 3), "| without J_g", round(mean(runs$ess_A), 3),
    "| IW", round(mean(runs$ess_IW), 3), "\n")
cat("Mixed model lm fallbacks:", sum(runs$mixed_fallback, na.rm = TRUE), "| errors:", sum(is.na(runs$mixed_fallback)), "\n")
cat("J_g = 0 in", sum(info$J_cf == 0), "of", nrow(info), "external subgroups;",
    "runs with all external weight 0:", sum(tapply(info$J_cf, paste(info$row, info$sim), max) == 0), "\n")
cat("\nMean RMSE per cell:\n")
print(as.data.frame(runs |> mutate(scenario = lab_s[scenario], similarity = lab_t[as.character(to)]) |>
  group_by(scenario, similarity) |>
  summarise(weighted = mean(rmse_D), without_Jg = mean(rmse_A), IW = mean(rmse_IW), global = mean(rmse_global),
            mixed = mean(rmse_mixed), .groups = "drop")), digits = 3)
cat("\nBias-variance:\n")
print(bv[bv$model %in% c("weighted_D", "weighted_A", "global", "env_specific", "mixed_model"), ], digits = 3, row.names = FALSE)
cat("\nRuntime (min):", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "\n")
