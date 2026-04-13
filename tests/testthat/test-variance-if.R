# Phase A unit tests for the IF variance engine primitives.
#
# These test bread_inv(), iv_design_matrix(), correct_model(),
# vcov_from_if() (with and without cluster), and variance_if_numeric()
# in isolation. The dispatcher variance_if() and correct_propensity() are
# tested in later phases through the public causat() / contrast() API.

# ── correct_model() — non-canonical link agreement with sandwich::estfun ──

test_that("correct_model() matches sandwich::estfun on a probit GLM", {
  set.seed(101)
  n <- 400
  L <- rnorm(n)
  A <- rbinom(n, 1, plogis(0.4 * L))
  Y <- rbinom(n, 1, pnorm(-0.3 + 0.7 * A + 0.5 * L))
  dat <- data.frame(Y = Y, A = A, L = L)

  m <- stats::glm(Y ~ A + L, data = dat, family = binomial(link = "probit"))

  # An arbitrary gradient g (the role of J in correct_model).
  set.seed(7)
  g <- runif(length(stats::coef(m)), -1, 1)

  res <- correct_model(m, g, fit_idx = seq_len(n), n_total = n)

  # Expected: n * (estfun %*% bread/nobs %*% g), per row.
  ef <- sandwich::estfun(m)
  br <- sandwich::bread(m) / stats::nobs(m)
  expected <- n * as.numeric(ef %*% br %*% g)

  expect_equal(res$correction, expected, tolerance = 1e-10)
})


# ── correct_model() — Gamma-log link agreement with sandwich::estfun ──

test_that("correct_model() matches sandwich::estfun on a Gamma-log GLM", {
  set.seed(202)
  n <- 500
  L <- rnorm(n)
  A <- rbinom(n, 1, 0.5)
  mu <- exp(0.4 + 0.3 * A + 0.2 * L)
  Y <- rgamma(n, shape = 5, rate = 5 / mu)
  dat <- data.frame(Y = Y, A = A, L = L)

  m <- stats::glm(Y ~ A + L, data = dat, family = Gamma(link = "log"))

  set.seed(8)
  g <- runif(length(stats::coef(m)), -1, 1)

  res <- correct_model(m, g, fit_idx = seq_len(n), n_total = n)

  ef <- sandwich::estfun(m)
  br <- sandwich::bread(m) / stats::nobs(m)
  expected <- n * as.numeric(ef %*% br %*% g)

  expect_equal(res$correction, expected, tolerance = 1e-10)
})


# ── Regression guard: response-residual formula must disagree on probit ──

test_that("response-residual formula disagrees materially on probit (regression guard)", {
  set.seed(303)
  n <- 400
  L <- rnorm(n)
  A <- rbinom(n, 1, 0.5)
  Y <- rbinom(n, 1, pnorm(-0.2 + 0.6 * A + 0.4 * L))
  dat <- data.frame(Y = Y, A = A, L = L)

  m <- stats::glm(Y ~ A + L, data = dat, family = binomial(link = "probit"))
  g <- c(0.1, 0.2, 0.3)

  res <- correct_model(m, g, fit_idx = seq_len(n), n_total = n)

  # Old formula: n * d_i * r_response (no mu_eta/V scale factor).
  X <- stats::model.matrix(m)
  B_inv <- bread_inv(m, X)
  h <- as.numeric(B_inv %*% g)
  d <- as.numeric(X %*% h)
  r_response <- stats::residuals(m, type = "response")
  old_correction <- n * d * r_response

  # The two must differ materially. Mean absolute difference normalised
  # by the magnitude of the new (correct) correction.
  rel_diff <- mean(abs(res$correction - old_correction)) /
    mean(abs(res$correction))
  expect_gt(rel_diff, 0.05)
})


# ── correct_model() — canonical-link sanity (logistic, identity, Poisson-log) ──

test_that("correct_model() canonical links match sandwich::estfun", {
  set.seed(404)
  n <- 300
  L <- rnorm(n)
  A <- rbinom(n, 1, 0.5)

  for (cfg in list(
    list(
      family = binomial(),
      Y = rbinom(n, 1, plogis(0.2 + 0.5 * A + 0.3 * L))
    ),
    list(family = gaussian(), Y = 1 + 0.5 * A + 0.3 * L + rnorm(n)),
    list(family = poisson(), Y = rpois(n, exp(0.2 + 0.4 * A + 0.2 * L)))
  )) {
    dat <- data.frame(Y = cfg$Y, A = A, L = L)
    m <- stats::glm(Y ~ A + L, data = dat, family = cfg$family)
    g <- c(0.1, -0.2, 0.3)
    res <- correct_model(m, g, fit_idx = seq_len(n), n_total = n)
    ef <- sandwich::estfun(m)
    br <- sandwich::bread(m) / stats::nobs(m)
    expect_equal(
      res$correction,
      n * as.numeric(ef %*% br %*% g),
      tolerance = 1e-10
    )
  }
})


# ── vcov_from_if() — independent aggregation ──

test_that("vcov_from_if() reproduces the (1/n^2) sum-of-products formula", {
  set.seed(505)
  n <- 50
  IFa <- rnorm(n)
  IFb <- rnorm(n)

  V <- vcov_from_if(list(a = IFa, b = IFb), n = n, int_names = c("a", "b"))

  expect_equal(V["a", "a"], sum(IFa * IFa) / n^2)
  expect_equal(V["b", "b"], sum(IFb * IFb) / n^2)
  expect_equal(V["a", "b"], sum(IFa * IFb) / n^2)
  expect_equal(V["a", "b"], V["b", "a"])
})


# ── vcov_from_if() — cluster aggregation sums within cluster ──

test_that("vcov_from_if() cluster path sums within cluster before squaring", {
  IF <- c(1, 2, 3, 4, 5, 6)
  cl <- c("A", "A", "A", "B", "B", "B")
  n <- length(IF)

  V <- vcov_from_if(list(x = IF), n = n, int_names = "x", cluster = cl)

  # Sum within clusters: A -> 6, B -> 15. Sum of squares -> 36 + 225 = 261.
  expect_equal(V["x", "x"], 261 / n^2)

  # Singleton clusters reduce to the standard formula.
  cl2 <- as.character(seq_along(IF))
  V2 <- vcov_from_if(list(x = IF), n = n, int_names = "x", cluster = cl2)
  V_indep <- vcov_from_if(list(x = IF), n = n, int_names = "x")
  expect_equal(V2, V_indep)
})


# ── Matching: IF vs vcovCL + J V_beta J^T at Channel-1 = 0 ──

test_that("variance_if matching path matches vcovCL+J V_beta J^T at saturated MSM", {
  set.seed(808)
  n <- 800
  L <- rnorm(n)
  A <- rbinom(n, 1, plogis(-0.7 + 0.4 * L))
  Y <- 1 + 0.5 * A + 0.3 * L + rnorm(n)
  df <- data.frame(Y = Y, A = A, L = L)

  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    method = "matching",
    estimand = "ATT"
  )

  res_if <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich"
  )
  V_if <- res_if$vcov

  # Reference: vcovCL on the matched glm propagated through the Jacobian
  # of mu_a = g^{-1}(beta_0 + beta_1 a). For a saturated linear MSM with
  # static interventions, channel 1 is identically zero and J V_beta J^T
  # is the entire variance.
  matched <- as.data.frame(fit$details$matched_data)
  # cadjust = FALSE disables sandwich::vcovCL's G/(G-1) finite-sample
  # correction, matching the raw 1/n^2 cluster-sum aggregation in
  # vcov_from_if(cluster = ...).
  V_beta <- sandwich::vcovCL(
    fit$model,
    cluster = matched$subclass,
    cadjust = FALSE
  )
  J <- matrix(c(1, 1, 1, 0), nrow = 2, byrow = TRUE) # rows: a1, a0
  V_ref <- J %*% V_beta %*% t(J)

  expect_equal(unname(V_if), unname(V_ref), tolerance = 1e-10)
})


# ── variance_if_numeric() Tier 1 — recovers analytic IF on a glm ──

test_that("variance_if_numeric() Tier 1 matches the analytic IF on a logistic glm", {
  set.seed(606)
  n <- 300
  L <- rnorm(n)
  A <- rbinom(n, 1, plogis(0.3 * L))
  Y <- rbinom(n, 1, plogis(-0.2 + 0.6 * A + 0.4 * L))
  dat <- data.frame(Y = Y, A = A, L = L)
  m <- stats::glm(Y ~ A + L, data = dat, family = binomial)

  fit <- list(
    model = m,
    details = list(fit_rows = rep(TRUE, n))
  )

  data_a_list <- list(
    a1 = transform(dat, A = 1),
    a0 = transform(dat, A = 0)
  )
  preds_list <- lapply(data_a_list, function(df) {
    stats::predict(m, newdata = df, type = "response")
  })
  mu_hat <- vapply(preds_list, mean, numeric(1))
  target_idx <- rep(TRUE, n)

  V_num <- variance_if_numeric(
    fit,
    data_a_list,
    preds_list,
    mu_hat,
    target_idx
  )

  # Analytic reference: build IF by hand using correct_model().
  beta <- stats::coef(m)
  IF_list <- lapply(seq_along(data_a_list), function(j) {
    df <- as.data.frame(data_a_list[[j]])
    X_star <- iv_design_matrix(m, df)
    eta_star <- as.numeric(X_star %*% beta)
    mu_eta_star <- m$family$mu.eta(eta_star)
    g <- as.numeric(crossprod(X_star, mu_eta_star)) / n
    Ch1 <- preds_list[[j]] - mu_hat[j]
    Ch2 <- correct_model(m, g, fit_idx = seq_len(n), n_total = n)$correction
    Ch1 + Ch2
  })
  V_ref <- vcov_from_if(IF_list, n = n, int_names = names(data_a_list))

  expect_equal(unname(V_num), unname(V_ref), tolerance = 1e-6)
})


# ── correct_propensity() Branch A — IPW WeightIt path ──

test_that("variance_if IPW path matches WeightIt vcov + J V_beta J^T at saturated MSM", {
  skip_if_not_installed("WeightIt")
  set.seed(909)
  n <- 1000
  L <- rnorm(n)
  A <- rbinom(n, 1, plogis(0.4 * L))
  Y <- 1 + 0.5 * A + 0.3 * L + rnorm(n)
  df <- data.frame(Y = Y, A = A, L = L)

  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    method = "ipw"
  )

  res_if <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich"
  )
  V_if <- res_if$vcov

  # Reference: Y ~ A is saturated, A is static, so Channel 1 = 0 and the
  # IF variance equals J V_beta J^T with V_beta = vcov(glm_weightit) (the
  # M-estimation sandwich that already accounts for weight uncertainty).
  V_beta <- stats::vcov(fit$model)
  J <- matrix(c(1, 1, 1, 0), nrow = 2, byrow = TRUE)
  V_ref <- J %*% V_beta %*% t(J)

  expect_equal(unname(V_if), unname(V_ref), tolerance = 1e-10)
})


# ── correct_propensity() Branch B — Phase 4 scaffold aborts informatively ──

test_that("correct_propensity() Branch B aborts with a Phase 4 message", {
  fake_fit <- list(
    model = structure(list(), class = "lm"),
    method = "ipw"
  )
  expect_error(
    correct_propensity(fake_fit, J = c(1, 0), fit_idx = 1L, n_total = 1L),
    "Branch B"
  )
})


# ── Mparts guardrail — warning at fit time for non-Mparts WeightIt method ──

test_that("fit_ipw warns when WeightIt method lacks Mparts", {
  skip_if_not_installed("WeightIt")

  set.seed(1010)
  n <- 100
  L <- rnorm(n)
  A <- rbinom(n, 1, plogis(0.3 * L))
  Y <- 1 + 0.5 * A + rnorm(n)
  df <- data.frame(Y = Y, A = A, L = L)

  fake_w <- WeightIt::weightit(A ~ L, data = df, method = "glm")
  attr(fake_w, "Mparts") <- NULL

  testthat::local_mocked_bindings(
    weightit = function(...) fake_w,
    .package = "WeightIt"
  )

  expect_warning(
    causat(
      df,
      outcome = "Y",
      treatment = "A",
      confounders = ~L,
      method = "ipw"
    ),
    "Mparts"
  )
})


# ── variance_if_numeric() Tier 2 — warns and returns V1 + J V_beta J^T ──

test_that("variance_if_numeric() Tier 2 warns and returns V1 + J V_beta J^T", {
  set.seed(707)
  n <- 200
  L <- rnorm(n)
  A <- rbinom(n, 1, 0.5)
  Y <- 1 + 0.4 * A + 0.3 * L + rnorm(n)
  dat <- data.frame(Y = Y, A = A, L = L)
  m <- stats::glm(Y ~ A + L, data = dat, family = gaussian())

  fit <- list(
    model = m,
    details = list(fit_rows = rep(TRUE, n))
  )
  data_a_list <- list(
    a1 = transform(dat, A = 1),
    a0 = transform(dat, A = 0)
  )
  preds_list <- lapply(data_a_list, function(df) {
    stats::predict(m, newdata = df, type = "response")
  })
  mu_hat <- vapply(preds_list, mean, numeric(1))
  target_idx <- rep(TRUE, n)

  # Force Tier 2 by stubbing sandwich::estfun to error for any input.
  testthat::local_mocked_bindings(
    estfun = function(x, ...) stop("no estfun for this class"),
    .package = "sandwich"
  )

  expect_warning(
    V <- variance_if_numeric(
      fit,
      data_a_list,
      preds_list,
      mu_hat,
      target_idx
    ),
    "drops the IF cross-term"
  )

  # Hand-computed reference: V1 (Channel 1 sum of squares) + J V_beta J^T.
  # Channel 1 enters as t_i * (p_i - mu_hat) with the 1/n^2 aggregation.
  V1 <- vapply(
    seq_along(preds_list),
    function(j) {
      sum((preds_list[[j]] - mu_hat[j])^2) / n^2
    },
    numeric(1)
  )
  beta <- stats::coef(m)
  pred_fun <- function(beta) {
    mt <- m
    mt$coefficients <- beta
    vapply(
      data_a_list,
      function(df) {
        mean(stats::predict(mt, newdata = df, type = "response"))
      },
      numeric(1)
    )
  }
  J <- numDeriv::jacobian(pred_fun, x = beta)
  V_beta <- stats::vcov(m)
  V_ref <- J %*% V_beta %*% t(J) + diag(V1, nrow = 2)

  expect_equal(unname(V), unname(V_ref), tolerance = 1e-8)
})
