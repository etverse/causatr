# Bootstrap variance tests for gcomp and IPW transport (chunk 17d).
#
# DGP (simulate_transport, helper-dgp.R):
#   L ~ N(0,1)
#   P(S=1|L) = expit(-0.5 + L)    [S=1 over-represents high-L]
#   A | L, S=1 ~ Bern(expit(0.2 + 0.3L))
#   Y | A, L ~ N(2 + 3A + 1.5L + A*L, 1)   [S=0 rows: Y = NA]
#
# Target ATE for transportability (target = S=0):
#   3 + E[L | S=0]  (< 3 because S=0 over-represents low-L)
# Target ATE for generalizability (target = all rows):
#   3 + E[L] ≈ 3 (L ~ N(0,1))
#
# In all bootstrap tests the point estimate equals the sandwich estimate
# (ci_method only affects the variance, not the plug-in); truth-based
# assertions check both the point estimate and that the bootstrap CI
# brackets the true value.

test_that("gcomp transport bootstrap (transportability): point estimate near truth", {
  d <- simulate_transport(n = 4000, seed = 42)
  truth_ate <- 3 + mean(d$L[d$S == 0])

  fit <- causat(
    d,
    "Y",
    "A",
    ~ L + A:L,
    estimator = "gcomp",
    target = "S",
    target_subset = "target"
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "bootstrap",
    n_boot = 200
  )

  est <- coef(res)["a1"] - coef(res)["a0"]
  expect_lt(abs(est - truth_ate), 0.15)
  expect_false(anyNA(coef(res)))
})

test_that("gcomp transport bootstrap (generalizability): point estimate near truth", {
  d <- simulate_transport(n = 4000, seed = 99)
  truth_ate <- 3 + mean(d$L)

  fit <- causat(
    d,
    "Y",
    "A",
    ~ L + A:L,
    estimator = "gcomp",
    target = "S",
    target_subset = "all"
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "bootstrap",
    n_boot = 200
  )

  est <- coef(res)["a1"] - coef(res)["a0"]
  expect_lt(abs(est - truth_ate), 0.15)
})

test_that("gcomp transport bootstrap: CI brackets truth", {
  d <- simulate_transport(n = 4000, seed = 13)
  truth_ate <- 3 + mean(d$L[d$S == 0])

  fit <- causat(
    d,
    "Y",
    "A",
    ~ L + A:L,
    estimator = "gcomp",
    target = "S",
    target_subset = "target"
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "bootstrap",
    n_boot = 200
  )

  v <- vcov(res)
  se <- sqrt(v["a1", "a1"] + v["a0", "a0"] - 2 * v["a1", "a0"])
  est <- coef(res)["a1"] - coef(res)["a0"]
  expect_true(est - 1.96 * se < truth_ate && truth_ate < est + 1.96 * se)
})

test_that("gcomp transport bootstrap: point estimate equals sandwich estimate", {
  # ci_method controls only the variance; the plug-in mean is identical.
  d <- simulate_transport(n = 2000, seed = 7)
  fit <- causat(
    d,
    "Y",
    "A",
    ~ L + A:L,
    estimator = "gcomp",
    target = "S",
    target_subset = "target"
  )
  ivs <- list(a1 = static(1), a0 = static(0))

  res_sw <- contrast(
    fit,
    interventions = ivs,
    type = "difference",
    ci_method = "sandwich"
  )
  res_bt <- contrast(
    fit,
    interventions = ivs,
    type = "difference",
    ci_method = "bootstrap",
    n_boot = 100
  )

  expect_equal(coef(res_sw), coef(res_bt))
})

test_that("IPW transport bootstrap (transportability): point estimate near truth", {
  # IPW is noisier than gcomp; larger n and wider tolerance.
  d <- simulate_transport(n = 10000, seed = 42)
  truth_ate <- 3 + mean(d$L[d$S == 0])

  fit <- causat(
    d,
    "Y",
    "A",
    ~L,
    estimator = "ipw",
    target = "S",
    target_subset = "target"
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "bootstrap",
    n_boot = 200
  )

  est <- coef(res)["a1"] - coef(res)["a0"]
  expect_lt(abs(est - truth_ate), 0.2)
  expect_false(anyNA(coef(res)))
})

test_that("IPW transport bootstrap (generalizability): point estimate near truth", {
  d <- simulate_transport(n = 10000, seed = 99)
  truth_ate <- 3 + mean(d$L)

  fit <- causat(
    d,
    "Y",
    "A",
    ~L,
    estimator = "ipw",
    target = "S",
    target_subset = "all"
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "bootstrap",
    n_boot = 200
  )

  est <- coef(res)["a1"] - coef(res)["a0"]
  expect_lt(abs(est - truth_ate), 0.2)
})

test_that("IPW transport bootstrap: CI brackets truth", {
  d <- simulate_transport(n = 8000, seed = 17)
  truth_ate <- 3 + mean(d$L[d$S == 0])

  fit <- causat(
    d,
    "Y",
    "A",
    ~L,
    estimator = "ipw",
    target = "S",
    target_subset = "target"
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "bootstrap",
    n_boot = 200
  )

  v <- vcov(res)
  se <- sqrt(v["a1", "a1"] + v["a0", "a0"] - 2 * v["a1", "a0"])
  est <- coef(res)["a1"] - coef(res)["a0"]
  expect_true(est - 1.96 * se < truth_ate && truth_ate < est + 1.96 * se)
})

test_that("IPW transport bootstrap: sampling model refitted per replicate", {
  # Directly exercise refit_ipw() to verify that each bootstrap replicate
  # gets a fresh sampling model fit (different gamma_hat from the original
  # because the bootstrap sample has a different empirical S distribution).
  d <- simulate_transport(n = 500, seed = 42)
  fit <- causat(d, "Y", "A", ~L, estimator = "ipw", target = "S")

  set.seed(1)
  n <- nrow(fit$data)
  idx <- sample(n, n, replace = TRUE)
  d_b <- fit$data[idx]
  fit_b <- refit_ipw(fit, d_b)

  expect_false(is.null(fit_b$details$sampling_model))
  # gamma_hat is the vector of estimated log-odds coefficients; they differ
  # between the original data and the bootstrap resample.
  expect_false(identical(
    fit_b$details$sampling_model$gamma_hat,
    fit$details$sampling_model$gamma_hat
  ))
})

test_that("gcomp vs IPW transport bootstrap: point estimates agree", {
  d <- simulate_transport(n = 10000, seed = 55)
  ivs <- list(a1 = static(1), a0 = static(0))

  # gcomp requires the A:L interaction to be specified; IPW needs only L for
  # the propensity model (treatment model logistic(A ~ L)).
  fit_gc <- causat(
    d,
    "Y",
    "A",
    ~ L + A:L,
    estimator = "gcomp",
    target = "S",
    target_subset = "target"
  )
  fit_ipw <- causat(
    d,
    "Y",
    "A",
    ~L,
    estimator = "ipw",
    target = "S",
    target_subset = "target"
  )

  res_gc <- contrast(
    fit_gc,
    interventions = ivs,
    type = "difference",
    ci_method = "bootstrap",
    n_boot = 200
  )
  res_ipw <- contrast(
    fit_ipw,
    interventions = ivs,
    type = "difference",
    ci_method = "bootstrap",
    n_boot = 200
  )

  est_gc <- coef(res_gc)["a1"] - coef(res_gc)["a0"]
  est_ipw <- coef(res_ipw)["a1"] - coef(res_ipw)["a0"]
  expect_lt(abs(est_gc - est_ipw), 0.3)
})
