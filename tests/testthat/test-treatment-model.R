# Unit tests for R/treatment_model.R (Phase 4 foundation layer).
#
# The treatment density model + density evaluator underpin the
# self-contained IPW engine. These tests pin down:
#   - family classification (binary vs continuous vs categorical)
#   - the fit pipeline end-to-end against `stats::glm` internals
#   - Bernoulli / Gaussian density evaluation against hand formulas
#   - the "categorical not yet" abort
#   - sigma extraction across GLM / LM fit shapes
#
# The categorical branch is explicitly deferred in this first Phase 4
# chunk, and the test locks in the abort message so a future PR that
# lands the multinomial branch has to remove this test alongside the
# abort.

test_that("detect_treatment_family() classifies each input type", {
  expect_identical(detect_treatment_family(c(0L, 1L, 0L, 1L)), "bernoulli")
  expect_identical(detect_treatment_family(c(0, 1, 0, 1)), "bernoulli")
  expect_identical(detect_treatment_family(c(TRUE, FALSE, TRUE)), "bernoulli")
  expect_identical(
    detect_treatment_family(rnorm(20)),
    "gaussian"
  )
  expect_identical(
    detect_treatment_family(factor(c("a", "b", "c", "a"))),
    "categorical"
  )
  expect_identical(
    detect_treatment_family(c("a", "b", "c", "a")),
    "categorical"
  )
})

test_that("detect_treatment_family() rejects unsupported types", {
  expect_error(
    detect_treatment_family(as.Date("2020-01-01") + 0:3),
    "is not supported"
  )
})

test_that("fit_treatment_model() on binary treatment returns a Bernoulli fit", {
  d <- simulate_binary_continuous(n = 500, seed = 1)
  dt <- data.table::as.data.table(d)

  tm <- fit_treatment_model(
    dt,
    treatment = "A",
    confounders = ~L,
    model_fn = stats::glm
  )

  expect_s3_class(tm, "causatr_treatment_model")
  expect_identical(tm$family, "bernoulli")
  expect_identical(tm$treatment, "A")
  expect_null(tm$sigma)
  expect_equal(sum(tm$fit_rows), nrow(d))
  # alpha_hat length matches propensity design matrix width
  expect_equal(length(tm$alpha_hat), ncol(tm$X_prop))
  # X_prop should be the fitted model's design matrix verbatim
  ref_fit <- stats::glm(A ~ L, data = d, family = stats::binomial())
  expect_equal(
    unname(tm$alpha_hat),
    unname(stats::coef(ref_fit)),
    tolerance = 1e-10
  )
})

test_that("fit_treatment_model() on continuous treatment recovers sigma", {
  d <- simulate_continuous_continuous(n = 500, seed = 2)
  dt <- data.table::as.data.table(d)

  tm <- fit_treatment_model(dt, treatment = "A", confounders = ~L)

  expect_identical(tm$family, "gaussian")
  expect_true(is.numeric(tm$sigma) && tm$sigma > 0)

  # Compare sigma and alpha against a plain glm() fit on the same data.
  # Note: `summary(glm)$sigma` is NULL for a Gaussian GLM — the
  # residual-variance estimator lives on `$dispersion`. Take sqrt to
  # get the residual SD.
  ref_fit <- stats::glm(A ~ L, data = d, family = stats::gaussian())
  expect_equal(
    tm$sigma,
    sqrt(summary(ref_fit)$dispersion),
    tolerance = 1e-10
  )
  expect_equal(
    unname(tm$alpha_hat),
    unname(stats::coef(ref_fit)),
    tolerance = 1e-10
  )
})

test_that("fit_treatment_model() aborts on categorical (Phase 4 pending)", {
  dt <- data.table::data.table(
    A = factor(sample(c("a", "b", "c"), 50, replace = TRUE)),
    L = rnorm(50)
  )
  expect_error(
    fit_treatment_model(dt, treatment = "A", confounders = ~L),
    class = "causatr_phase4_categorical_pending"
  )
})

test_that("fit_treatment_model() drops rows with NA treatment / NA confounders", {
  d <- simulate_binary_continuous(n = 200, seed = 3)
  d$A[1:5] <- NA
  d$L[10:12] <- NA
  dt <- data.table::as.data.table(d)

  tm <- fit_treatment_model(dt, treatment = "A", confounders = ~L)
  # 5 NA in A + 3 NA in L (non-overlapping) = 8 dropped rows
  expect_equal(sum(tm$fit_rows), nrow(d) - 8L)
})

test_that("evaluate_density() binary returns Bernoulli pmf per row", {
  d <- simulate_binary_continuous(n = 300, seed = 4)
  dt <- data.table::as.data.table(d)
  tm <- fit_treatment_model(dt, treatment = "A", confounders = ~L)

  fit_data <- dt[tm$fit_rows]
  a_obs <- fit_data$A
  p <- stats::predict(tm$model, newdata = fit_data, type = "response")

  # Evaluate at the observed treatment: ifelse(a == 1, p, 1 - p)
  f_obs <- evaluate_density(tm, a_obs, fit_data)
  expect_equal(f_obs, ifelse(a_obs == 1, p, 1 - p), tolerance = 1e-12)

  # Evaluate at a constant 1: f_int is p for every row
  f_int_1 <- evaluate_density(tm, rep(1, nrow(fit_data)), fit_data)
  expect_equal(f_int_1, unname(p), tolerance = 1e-12)

  # Evaluate at a constant 0: f_int is 1 - p
  f_int_0 <- evaluate_density(tm, rep(0, nrow(fit_data)), fit_data)
  expect_equal(f_int_0, unname(1 - p), tolerance = 1e-12)
})

test_that("evaluate_density() continuous returns Gaussian pdf per row", {
  d <- simulate_continuous_continuous(n = 300, seed = 5)
  dt <- data.table::as.data.table(d)
  tm <- fit_treatment_model(dt, treatment = "A", confounders = ~L)

  fit_data <- dt[tm$fit_rows]
  a_obs <- fit_data$A
  mu <- stats::predict(tm$model, newdata = fit_data, type = "response")

  f_obs <- evaluate_density(tm, a_obs, fit_data)
  expect_equal(
    f_obs,
    stats::dnorm(a_obs, mean = mu, sd = tm$sigma),
    tolerance = 1e-12
  )

  # Density at a shifted treatment matches dnorm at the shifted point
  f_shift <- evaluate_density(tm, a_obs - 0.5, fit_data)
  expect_equal(
    f_shift,
    stats::dnorm(a_obs - 0.5, mean = mu, sd = tm$sigma),
    tolerance = 1e-12
  )
})

test_that("evaluate_density() rejects length mismatch", {
  d <- simulate_binary_continuous(n = 100, seed = 6)
  dt <- data.table::as.data.table(d)
  tm <- fit_treatment_model(dt, treatment = "A", confounders = ~L)

  fit_data <- dt[tm$fit_rows]
  expect_error(
    evaluate_density(tm, rep(1, 3), fit_data),
    "length"
  )
})

test_that("extract_sigma() recovers sigma from a Gaussian GLM via dispersion", {
  set.seed(7)
  df <- data.frame(A = rnorm(100), L = rnorm(100))
  glm_fit <- stats::glm(A ~ L, data = df, family = stats::gaussian())
  # `summary(glm)$sigma` is NULL on a Gaussian GLM (common gotcha);
  # `summary(glm)$dispersion` is the estimator of sigma^2.
  expect_equal(
    extract_sigma(glm_fit),
    sqrt(summary(glm_fit)$dispersion),
    tolerance = 1e-12
  )
})

test_that("extract_sigma() recovers sigma from an lm fit via $sigma", {
  set.seed(8)
  df <- data.frame(A = rnorm(100), L = rnorm(100))
  lm_fit <- stats::lm(A ~ L, data = df)
  expect_equal(extract_sigma(lm_fit), summary(lm_fit)$sigma, tolerance = 1e-12)
})
