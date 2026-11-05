# Small simulated data: target cohort 1 and three external cohorts with increasing shift in
# the covariates and in the outcome given the covariates.
sim_data <- function(n = c(15, 20, 25, 30), shift = c(0, 0.5, 1.5, 4)) {
  cohort <- rep(seq_along(n), n)
  s <- shift[cohort]
  x1 <- rnorm(sum(n), mean = s)
  x2 <- rnorm(sum(n))
  y <- x1 + x2 + s + rnorm(sum(n))
  data.frame(cohort = cohort, x1 = x1, x2 = x2, y = y)
}

test_that("create_cow() equals the reference implementation of method D", {
  set.seed(2026)
  for (r in 1:20) {
    d <- sim_data()
    ref <- membership_weights(d, "cohort", 1, c("x1", "x2", "y"), seed = r)
    w <- create_cow(d, "cohort", 1, c("x1", "x2"), outcome = "y", seed = r)
    expect_equal(w$weight, unname(ref$w[, "D"]), tolerance = 1e-12)
    expect_equal(attr(w, "subgroups")$J, ref$info$J_cf, tolerance = 1e-12)
  }
})

test_that("the outcome can also enter through baseline_covs", {
  set.seed(1)
  d <- sim_data()
  w1 <- create_cow(d, "cohort", 1, c("x1", "x2"), outcome = "y", seed = 5)
  w2 <- create_cow(d, "cohort", 1, c("x1", "x2", "y"), seed = 5)
  expect_identical(w1$weight, w2$weight)
})

test_that("target weights are 1, external weights average to J_g, row order is kept", {
  set.seed(3)
  d <- sim_data()
  d <- d[sample(nrow(d)), ]
  w <- create_cow(d, "cohort", 1, c("x1", "x2"), outcome = "y", seed = 1)
  expect_identical(w[, names(d)], d)
  expect_true(all(w$weight[d$cohort == 1] == 1))
  expect_true(all(is.na(w$ps_weight[d$cohort == 1])))
  expect_identical(w$c_flag, ifelse(d$cohort == 1, 1, 0))
  sg <- attr(w, "subgroups")
  for (i in seq_len(nrow(sg))) {
    expect_equal(mean(w$weight[d$cohort == sg$cohort[i]]), sg$J[i], tolerance = 1e-12)
  }
  expect_true(all(sg$J >= 0 & sg$J <= 1))
})

test_that("a seed leaves the global random number stream unchanged", {
  set.seed(4)
  d <- sim_data()
  before <- .Random.seed
  create_cow(d, "cohort", 1, c("x1", "x2"), outcome = "y", seed = 10)
  expect_identical(.Random.seed, before)
})

test_that("missing columns and missing values give errors", {
  set.seed(5)
  d <- sim_data()
  expect_error(create_cow(d, "cohort", 1, c("x1", "x3")), "Columns not found")
  expect_error(create_cow(d, "cohort", 9, c("x1", "x2")), "target_cohort not found")
  d$x1[3] <- NA
  expect_error(create_cow(d, "cohort", 1, c("x1", "x2")), "Missing values")
})
