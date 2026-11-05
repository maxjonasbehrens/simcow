# Table 2 (sensitivity to the handling of the weights). All four handling methods run on the
# 32,400 data sets of the main grid and share the cross-validated J_g, so "stabilize" equals
# the weighted sample of Table 1 (method D). Each method is also computed with J_g from the
# in-sample AUC (as in the first submission) on the same data, which separates the effect of
# the J_g estimate from the effect of the data sets.
# Usage (from any directory): Rscript paper/run_table2.R
# Needs the replayed data sets in paper/output/cache/ (written by run_simulation.R).
# Writes paper/output/table2_runs.csv and table2_summary.csv. Run time under 15 min on 6 cores;
# stops if the results differ from paper/results/.
suppressPackageStartupMessages({ library(dplyr); library(parallel) })
here <- normalizePath(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(here, "R", "methods.R"))
source(file.path(here, "R", "simulate_scenario.R"))
source(file.path(here, "..", "R", "create_cow.R"))
outdir <- file.path(here, "output"); cachedir <- file.path(outdir, "cache")
cores <- if (.Platform$OS.type == "windows") 1L else min(6L, detectCores(), na.rm = TRUE)
t0 <- Sys.time()
ok <- function(cond, msg) { cat(if (cond) "PASS" else "FAIL", msg, "\n"); if (!cond) quit(status = 1) }

H <- c("none", "truncate", "stabilize", "calibrate")

# Weights of the four handling methods, as in the first submission (weight_type = "ps"):
# r_j = J_g (n_g / n_0) p_jg, then
#   none: as is; truncate: capped at the 99th percentile within the subgroup;
#   calibrate: p_jg recalibrated by Platt scaling before weighting;
#   stabilize: J_g p_jg / mean(p_g), with J_g applied after the division
#              (cf = method D, code = variant B).
# J_g^code uses the in-sample AUC (of the calibrated p_jg for calibrate, as submitted);
# J_g^cf is the cross-validated J_g of membership_weights(), with the same folds.
# cf = FALSE skips the cross-validation and returns the J_g^code columns only.
handling_weights <- function(data, cohort_id, target, covs, K = 5, seed = 1, cf = TRUE) {
  f_ps <- as.formula(paste("c_flag ~", paste(covs, collapse = " + ")))
  groups <- as.character(data[[cohort_id]])
  js <- if (cf) c("cf", "code") else "code"
  w <- matrix(1, nrow(data), 4 * length(js), dimnames = list(NULL, as.vector(outer(H, js, paste, sep = "_"))))
  with_seed(seed, {
    for (g in setdiff(unique(groups), as.character(target))) {
      pair <- groups %in% c(as.character(target), g)
      pd <- data[pair, , drop = FALSE]
      pd$c_flag <- ifelse(groups[pair] == g, 0, 1)
      p <- predict(suppressWarnings(glm(f_ps, data = pd, family = binomial)), type = "response")
      J <- c(code = overlap_index(code_auc(p, pd$c_flag)))
      if (cf) {
        folds <- sample(rep(1:K, length.out = nrow(pd)))
        oof <- numeric(nrow(pd))
        for (f in 1:K) if (any(folds == f)) {
          fit <- suppressWarnings(glm(f_ps, data = pd[folds != f, , drop = FALSE], family = binomial))
          oof[folds == f] <- predict(fit, newdata = pd[folds == f, , drop = FALSE], type = "response")
        }
        J[["cf"]] <- overlap_index(std_auc(oof, pd$c_flag))
      }
      p_clip <- pmin(pmax(p, 1e-6), 1 - 1e-6)
      cal <- tryCatch(suppressWarnings(glm(y ~ s, data = data.frame(y = pd$c_flag, s = qlogis(p_clip)),
                                           family = binomial())), error = function(e) NULL)
      p_cal <- if (is.null(cal)) p else predict(cal, newdata = data.frame(s = qlogis(p_clip)), type = "response")
      J_cal <- c(code = overlap_index(code_auc(p_cal, pd$c_flag)), cf = if (cf) J[["cf"]] else NA)
      ext <- pd$c_flag == 0
      ratio <- sum(ext) / sum(!ext)
      idx <- which(groups == g)
      for (j in js) {
        r <- J[[j]] * ratio * p[ext]
        w[idx, paste0("none_", j)] <- r
        w[idx, paste0("truncate_", j)] <- pmin(r, quantile(r, 0.99, na.rm = TRUE))
        w[idx, paste0("stabilize_", j)] <- J[[j]] * p[ext] / mean(p[ext])
        w[idx, paste0("calibrate_", j)] <- J_cal[[j]] * ratio * p_cal[ext]
      }
    }
  })
  w
}

# ---------------- Tests: "stabilize" equals method D (create_cow) and variant B ----------------
set.seed(99)
covs <- c("X_global_1", paste0("X_env_", 1:3), "Y")
max_diff <- c(D = 0, B = 0)
for (r in 1:40) {
  scn <- sample(c("AX", "AY", "AH"), 1); n_cat <- sample(c(2, 4, 8), 1); to <- sample(1:3, 1)
  d <- simulate_scenario(type = scn, target_environment = 1, n_target = sample(c(10, 20), 1),
                         n = 20 * n_cat, n_cat = n_cat, num_global_predictors = 1,
                         num_env_predictors = 3, from = 0, to = to)
  hw <- handling_weights(d, "A", 1, covs, seed = r)
  cw <- create_cow(d, "A", 1, covs[1:4], outcome = "Y", seed = r)
  mw <- membership_weights(d, "A", 1, covs, seed = r)
  max_diff[["D"]] <- max(max_diff[["D"]], abs(cw$weight - hw[, "stabilize_cf"]))
  max_diff[["B"]] <- max(max_diff[["B"]], abs(mw$w[, "B"] - hw[, "stabilize_code"]))
}
ok(max_diff[["D"]] < 1e-8, sprintf("stabilize_cf == create_cow() (max diff %.1e)", max_diff[["D"]]))
ok(max_diff[["B"]] < 1e-8, sprintf("stabilize_code == variant B (max diff %.1e)", max_diff[["B"]]))

# ---------------- Data sets (replayed by run_simulation.R) ----------------
cache <- file.path(cachedir, "replay_kept.rds"); cache_none <- file.path(cachedir, "replay_kept_none.rds")
if (!file.exists(cache) || !file.exists(cache_none)) stop("Replayed data sets missing: run paper/run_simulation.R first")
kept <- readRDS(cache); kept_none <- readRDS(cache_none)
param_grid <- expand.grid(scenario = c("AX","AY","AH"), n_cat = c(2,4,8,10), n_train = c(10,20,30),
                          n_target = c(10,20,30), to = c(3,2,1), num_env_predictors = c(3),
                          stabilization_method = c("none","truncate","stabilize","calibrate"),
                          stringsAsFactors = FALSE)
cat("Data sets:", length(kept), "(main grid),", length(kept_none), "(none row of the first submission)\n")

# ---------------- Runs ----------------
rmse <- function(y, yhat) sqrt(mean((y - yhat)^2))
fit_weighted <- function(d, f_y, w) { d$.w <- w; lm(f_y, data = d, weights = .w) }
f_y <- Y ~ X_global_1 + X_env_1 + X_env_2 + X_env_3

one <- function(k) {
  x <- kept[[k]]; d <- x$tr; te <- x$te; p <- param_grid[x$i, ]
  hw <- handling_weights(d, "A", 1, covs, seed = 1e6 + k)                # same folds as method D
  r <- apply(hw, 2, function(w) rmse(te$Y, predict(fit_weighted(d, f_y, w), newdata = te)))
  data.frame(row = x$i, sim = x$s, p[, c("scenario","n_cat","n_train","n_target","to")],
             t(setNames(r, paste0("rmse_", colnames(hw)))),
             all_zero_cf = all(hw[d$A != "1", "stabilize_cf"] == 0))
}
one_none <- function(k) {   # "none" row of the first submission on its own data sets
  x <- kept_none[[k]]; d <- x$tr; te <- x$te
  hw <- handling_weights(d, "A", 1, covs, cf = FALSE)
  rmse(te$Y, predict(fit_weighted(d, f_y, hw[, "none_code"]), newdata = te))
}
res <- mclapply(seq_along(kept), one, mc.cores = cores)
errs <- Filter(function(r) inherits(r, "try-error"), res)
if (length(errs)) { cat("Errors in", length(errs), "jobs; first:\n"); print(errs[[1]]); quit(status = 1) }
runs <- do.call(rbind, res)
res_none <- mclapply(seq_along(kept_none), one_none, mc.cores = cores)
if (any(sapply(res_none, inherits, "try-error"))) { cat("Errors in none block\n"); quit(status = 1) }
rmse_none_submitted <- unlist(res_none)
write.csv(runs, file.path(outdir, "table2_runs.csv"), row.names = FALSE)
cat("Runs:", nrow(runs), "+", length(rmse_none_submitted), ",",
    round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")

# ---------------- Checks ----------------
main <- read.csv(file.path(here, "results", "sim_runs.csv"))
m <- merge(runs, main[, c("row", "sim", "rmse_D", "rmse_B")], by = c("row", "sim"))
ok(nrow(m) == 32400, "all 32,400 runs matched to paper/results/sim_runs.csv")
ok(max(abs(m$rmse_stabilize_cf - m$rmse_D)) < 1e-9, "stabilize_cf == weighted sample of Table 1")
ok(max(abs(m$rmse_stabilize_code - m$rmse_B)) < 1e-9, "stabilize_code == variant B of the main run")
# "none" row of Table 2 in the first submission: mean, SD, median, quartiles
q <- quantile(rmse_none_submitted, c(.25, .5, .75))
none_chk <- c(mean(rmse_none_submitted) - 1.5937304545503801, sd(rmse_none_submitted) - 1.0492923067356801,
              q[[2]] - 1.49886348536314, q[[1]] - 1.3664107464708, q[[3]] - 1.67305669961919)
ok(max(abs(none_chk)) < 1e-9, sprintf("'none' row of the first submission reproduced (max diff %.1e)", max(abs(none_chk))))
shipped <- read.csv(file.path(here, "results", "table2_runs.csv"))
s <- merge(shipped, runs, by = c("row", "sim"), suffixes = c(".shipped", ""))
num <- grep("^rmse_", names(runs), value = TRUE)
ok(nrow(s) == nrow(runs) && identical(s$all_zero_cf.shipped, s$all_zero_cf) &&
   max(sapply(num, function(v) max(abs(s[[paste0(v, ".shipped")]] - s[[v]])))) < 1e-9,
   "every run equals paper/results/table2_runs.csv")

# ---------------- Summary ----------------
summ <- function(x) { q <- quantile(x, c(.25, .5, .75)); c(mean = mean(x), sd = sd(x), median = q[[2]], q25 = q[[1]], q75 = q[[3]]) }
tab <- do.call(rbind, lapply(paste0("rmse_", c(paste0(H, "_cf"), paste0(H, "_code"))), function(v) {
  x <- runs[[v]]; dd <- x - runs$rmse_stabilize_cf; ci <- t.test(dd)$conf.int
  data.frame(method = sub("rmse_", "", v), t(summ(x)), sd_without_all_zero = sd(x[!runs$all_zero_cf]),
             diff_vs_D = mean(dd), diff_lo = ci[1], diff_hi = ci[2], share_D_better = mean(dd > 0))
}))
write.csv(tab, file.path(outdir, "table2_summary.csv"), row.names = FALSE)
cat("\nRuns with all external weights zero under the cross-validated J_g:", sum(runs$all_zero_cf), "\n")
cat("\nTable 2: RMSE over 32,400 runs; diff = method minus stabilize_cf (method D), paired:\n")
print(cbind(tab[1], round(tab[-1], 3)), row.names = FALSE)
