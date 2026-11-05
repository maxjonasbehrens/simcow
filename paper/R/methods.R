# Weight variants of the paper. Variant D is the method of the paper and equals
# simcow::create_cow() (tested in tests/testthat/test-create_cow.R); variants A-C and the
# covariate-shift comparator are the ablations and comparators in Tables 1 and 2.

# Two-level weights for borrowing from external subgroups (revised method "D").
#
# For external subgroup g, a logistic membership model discriminates the target
# subgroup from g using covariates and outcome. Weight of external observation j:
#
#   w_j = J_g * p_jg / mean_{i in g}(p_ig),   J_g = 1 - |2 AUC_g - 1|,
#
# with AUC_g the standard two-class AUC of 5-fold cross-validated membership
# predictions. Target observations get weight 1.
#
# membership_weights() also returns the ablation and sensitivity variants:
#   A  p/mean(p), zero if J_g^code = 0      (submitted version; J_g cancels)
#   B  J_g^code * p/mean(p)                 (AUC as in the submitted code)
#   C  J_g^std  * p/mean(p)                 (standard AUC, in-sample)
#   D  J_g^cf   * p/mean(p)                 (standard AUC, cross-validated; revised method)

# Standard two-class AUC: P(score of target > score of external), ties count 1/2
std_auc <- function(p, y) {
  r <- rank(p); n1 <- sum(y == 1); n0 <- sum(y == 0)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# AUC as computed by calculate_c_index() in the submitted code: additionally counts
# pairs from the same group with equal scores as ties
code_auc <- function(p, y) {
  diffs <- outer(p[y == 1], p[y == 0], "-")
  conc <- sum(diffs > 0); disc <- sum(diffs < 0); ties <- sum(diffs == 0)
  same_ties <- sum(sapply(split(p, y), function(v) sum(choose(rle(sort(v))$lengths, 2))))
  total <- conc + disc + ties + same_ties
  if (total == 0) return(0.5)
  (conc + 0.5 * (ties + same_ties)) / total
}

overlap_index <- function(auc) 1 - abs(2 * auc - 1)

# Evaluate expr with its own seed and leave the global RNG stream untouched
with_seed <- function(seed, expr) {
  old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else assign(".Random.seed", old, envir = globalenv()))
  set.seed(seed)
  expr
}

membership_weights <- function(data, cohort_id, target, covs, K = 5, seed = 1) {
  f_ps <- as.formula(paste("c_flag ~", paste(covs, collapse = " + ")))
  groups <- as.character(data[[cohort_id]])
  w <- matrix(1, nrow(data), 4, dimnames = list(NULL, c("A", "B", "C", "D")))
  info <- list()
  with_seed(seed, {
    for (g in setdiff(unique(groups), as.character(target))) {
      pair <- groups %in% c(as.character(target), g)
      pd <- data[pair, , drop = FALSE]
      pd$c_flag <- ifelse(groups[pair] == g, 0, 1)
      p <- predict(suppressWarnings(glm(f_ps, data = pd, family = binomial)), type = "response")
      folds <- sample(rep(1:K, length.out = nrow(pd)))
      oof <- numeric(nrow(pd))
      for (f in 1:K) if (any(folds == f)) {
        fit <- suppressWarnings(glm(f_ps, data = pd[folds != f, , drop = FALSE], family = binomial))
        oof[folds == f] <- predict(fit, newdata = pd[folds == f, , drop = FALSE], type = "response")
      }
      J <- c(code = overlap_index(code_auc(p, pd$c_flag)),
             std = overlap_index(std_auc(p, pd$c_flag)),
             cf = overlap_index(std_auc(oof, pd$c_flag)))
      pe <- p[pd$c_flag == 0]
      base <- pe / mean(pe)
      idx <- which(groups == g)
      w[idx, "A"] <- if (J[["code"]] > 0) base else 0
      w[idx, "B"] <- J[["code"]] * base
      w[idx, "C"] <- J[["std"]] * base
      w[idx, "D"] <- J[["cf"]] * base
      info[[g]] <- data.frame(group = g, n_g = length(idx), J_code = J[["code"]],
                              J_std = J[["std"]], J_cf = J[["cf"]])
    }
  })
  list(w = w, info = do.call(rbind, info))
}

# Covariate-shift comparator as in the submitted code: odds of an X-only membership
# model, standardized within each external subgroup, zero if J_g^code = 0
iw_weights <- function(data, cohort_id, target, covs_x) {
  f_ps <- as.formula(paste("c_flag ~", paste(covs_x, collapse = " + ")))
  groups <- as.character(data[[cohort_id]])
  w <- rep(1, nrow(data))
  for (g in setdiff(unique(groups), as.character(target))) {
    pair <- groups %in% c(as.character(target), g)
    pd <- data[pair, , drop = FALSE]
    pd$c_flag <- ifelse(groups[pair] == g, 0, 1)
    p <- predict(suppressWarnings(glm(f_ps, data = pd, family = binomial)), type = "response")
    odds <- p / (1 - p + 1e-8)
    oe <- odds[pd$c_flag == 0]
    w[groups == g] <- if (overlap_index(code_auc(p, pd$c_flag)) > 0) oe / mean(oe) else 0
  }
  w
}

kish_ess <- function(w) sum(w)^2 / sum(w^2)
