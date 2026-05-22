# Bootstrap variance tests for SNM estimator (chunk 18i).
# Validates that bootstrap SEs agree with sandwich SEs across
# point / longitudinal / treatment-free / treatment_values paths.

test_that("SNM point bootstrap: SE consistent with sandwich", {
  set.seed(42)
  n <- 1000
  L <- rnorm(n)
  M <- as.integer(L > 0)
  A <- 0.5 * L + rnorm(n)
  Y <- 2 + 3 * A + 1.5 * L + 2 * A * M + rnorm(n)
  d <- data.table::data.table(L = L, M = M, A = A, Y = Y)

  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L + A:M,
    estimator = "snm",
    propensity_model_fn = stats::glm
  )

  set.seed(123)
  res_boot <- contrast(fit, ci_method = "bootstrap", n_boot = 500)
  res_sand <- contrast(fit, ci_method = "sandwich")

  expect_equal(res_boot$estimates$estimate, res_sand$estimates$estimate)

  se_ratio <- res_boot$estimates$se / res_sand$estimates$se
  expect_true(all(se_ratio > 0.5 & se_ratio < 2.0))

  expect_equal(res_boot$estimates$estimate[1], 3, tolerance = 0.3)
  expect_equal(res_boot$estimates$estimate[2], 2, tolerance = 0.5)
  expect_equal(res_boot$ci_method, "bootstrap")

  eig <- eigen(res_boot$vcov, only.values = TRUE)$values
  expect_true(all(eig >= -1e-10))
})


test_that("SNM point bootstrap: no EM (single psi)", {
  set.seed(42)
  n <- 800
  L <- rnorm(n)
  A <- 0.5 * L + rnorm(n)
  Y <- 2 + 3 * A + 1.5 * L + rnorm(n)
  d <- data.table::data.table(L = L, A = A, Y = Y)

  expect_message(
    fit <- causat(
      d,
      outcome = "Y",
      treatment = "A",
      confounders = ~L,
      estimator = "snm",
      propensity_model_fn = stats::glm
    ),
    class = "causatr_snm_no_em"
  )

  set.seed(123)
  res_boot <- contrast(fit, ci_method = "bootstrap", n_boot = 300)
  res_sand <- contrast(fit, ci_method = "sandwich")

  expect_equal(res_boot$estimates$estimate, res_sand$estimates$estimate)
  se_ratio <- res_boot$estimates$se / res_sand$estimates$se
  expect_true(se_ratio > 0.5 && se_ratio < 2.0)
  expect_equal(res_boot$estimates$estimate, 3, tolerance = 0.2)
})


test_that("SNM point bootstrap + treatment-free: SE consistent", {
  set.seed(42)
  n <- 1000
  L <- rnorm(n)
  M <- as.integer(L > 0)
  A <- 0.5 * L + rnorm(n)
  Y <- 2 + 3 * A + 1.5 * L + 2 * A * M + rnorm(n)
  d <- data.table::data.table(L = L, M = M, A = A, Y = Y)

  fit_tf <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L + A:M,
    estimator = "snm",
    propensity_model_fn = stats::glm,
    treatment_free = ~L
  )

  set.seed(123)
  res_boot <- contrast(fit_tf, ci_method = "bootstrap", n_boot = 500)
  res_sand <- contrast(fit_tf, ci_method = "sandwich")

  expect_equal(res_boot$estimates$estimate, res_sand$estimates$estimate)
  se_ratio <- res_boot$estimates$se / res_sand$estimates$se
  expect_true(all(se_ratio > 0.5 & se_ratio < 2.0))

  fit_no_tf <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L + A:M,
    estimator = "snm",
    propensity_model_fn = stats::glm
  )
  res_no_tf <- contrast(fit_no_tf, ci_method = "sandwich")
  expect_true(all(res_sand$estimates$se < res_no_tf$estimates$se))
})


test_that("SNM point bootstrap with treatment_values", {
  set.seed(42)
  n <- 1000
  L <- rnorm(n)
  M <- as.integer(L > 0)
  A <- 0.5 * L + rnorm(n)
  Y <- 2 + 3 * A + 1.5 * L + 2 * A * M + rnorm(n)
  d <- data.table::data.table(L = L, M = M, A = A, Y = Y)

  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L + A:M,
    estimator = "snm",
    propensity_model_fn = stats::glm
  )

  set.seed(123)
  res_boot <- contrast(
    fit,
    treatment_values = c(0, 1),
    ci_method = "bootstrap",
    n_boot = 500
  )
  res_sand <- contrast(fit, treatment_values = c(0, 1), ci_method = "sandwich")

  expect_equal(res_boot$estimates$estimate, res_sand$estimates$estimate)

  se_ratio <- res_boot$estimates$se / res_sand$estimates$se
  expect_true(se_ratio > 0.5 && se_ratio < 2.0)

  # CI brackets truth: avg blip = 3 + 2 * P(M=1) ~ 4
  expect_true(
    res_boot$estimates$ci_lower < 4 &&
      res_boot$estimates$ci_upper > 4
  )
})


test_that("SNM point bootstrap: binary treatment + EM", {
  set.seed(42)
  n <- 1000
  L <- rnorm(n)
  M <- as.integer(L > 0)
  A <- rbinom(n, 1, plogis(0.5 * L))
  Y <- 2 + 1 * A + 1.5 * L + 0.8 * A * M + rnorm(n)
  d <- data.table::data.table(L = L, M = M, A = A, Y = Y)

  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L + A:M,
    estimator = "snm",
    propensity_model_fn = stats::glm
  )

  set.seed(123)
  res_boot <- contrast(fit, ci_method = "bootstrap", n_boot = 400)
  res_sand <- contrast(fit, ci_method = "sandwich")

  se_ratio <- res_boot$estimates$se / res_sand$estimates$se
  expect_true(all(se_ratio > 0.5 & se_ratio < 2.0))

  expect_equal(res_boot$estimates$estimate[1], 1, tolerance = 0.4)
  expect_equal(res_boot$estimates$estimate[2], 0.8, tolerance = 0.6)
})


test_that("SNM longitudinal bootstrap: SE consistent with sandwich", {
  dgp <- simulate_snm_longitudinal(n = 800, seed = 42)

  expect_message(
    expect_warning(
      fit <- causat(
        dgp$data,
        outcome = "Y",
        treatment = "A",
        confounders = ~L,
        confounders_tv = ~L,
        estimator = "snm",
        type = "longitudinal",
        id = "id",
        time = "time",
        propensity_model_fn = stats::glm
      )
    )
  )

  set.seed(123)
  res_boot <- contrast(fit, ci_method = "bootstrap", n_boot = 400)
  res_sand <- contrast(fit, ci_method = "sandwich")

  expect_equal(res_boot$estimates$estimate, res_sand$estimates$estimate)

  se_ratio <- res_boot$estimates$se / res_sand$estimates$se
  expect_true(all(se_ratio > 0.4 & se_ratio < 2.5))

  expect_equal(res_boot$estimates$estimate[1], 3.15, tolerance = 0.4)
  expect_equal(res_boot$estimates$estimate[2], 3.0, tolerance = 0.3)
  expect_equal(res_boot$ci_method, "bootstrap")

  expect_equal(nrow(res_boot$vcov), 2)
  eig <- eigen(res_boot$vcov, only.values = TRUE)$values
  expect_true(all(eig >= -1e-10))
})


test_that("SNM longitudinal bootstrap + TF: SE consistent", {
  dgp <- simulate_snm_longitudinal(n = 800, seed = 42)

  expect_message(
    expect_warning(
      fit_tf <- causat(
        dgp$data,
        outcome = "Y",
        treatment = "A",
        confounders = ~L,
        confounders_tv = ~L,
        estimator = "snm",
        type = "longitudinal",
        id = "id",
        time = "time",
        propensity_model_fn = stats::glm,
        treatment_free = ~L
      )
    )
  )

  set.seed(123)
  res_boot <- contrast(fit_tf, ci_method = "bootstrap", n_boot = 400)
  res_sand <- contrast(fit_tf, ci_method = "sandwich")

  expect_equal(res_boot$estimates$estimate, res_sand$estimates$estimate)
  se_ratio <- res_boot$estimates$se / res_sand$estimates$se
  expect_true(all(se_ratio > 0.4 & se_ratio < 2.5))
})


test_that("SNM longitudinal bootstrap with TV-EM", {
  dgp <- simulate_snm_longitudinal_tv_em(n = 800, seed = 42)

  fit <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders_outcome = ~ A:M,
    confounders_tv = ~ L + M,
    estimator = "snm",
    type = "longitudinal",
    id = "id",
    time = "time",
    propensity_model_fn = stats::glm,
    history = 0
  )

  set.seed(123)
  res_boot <- contrast(fit, ci_method = "bootstrap", n_boot = 400)
  res_sand <- contrast(fit, ci_method = "sandwich")

  expect_equal(res_boot$estimates$estimate, res_sand$estimates$estimate)

  expect_equal(nrow(res_boot$estimates), 4)
  se_ratio <- res_boot$estimates$se / res_sand$estimates$se
  expect_true(all(se_ratio > 0.4 & se_ratio < 2.5))
})


test_that("SNM bootstrap no longer rejected", {
  set.seed(42)
  n <- 100
  L <- rnorm(n)
  A <- 0.5 * L + rnorm(n)
  Y <- 2 + 3 * A + L + rnorm(n)
  d <- data.table::data.table(L = L, A = A, Y = Y)

  expect_message(
    fit <- causat(
      d,
      outcome = "Y",
      treatment = "A",
      confounders = ~L,
      estimator = "snm",
      propensity_model_fn = stats::glm
    ),
    class = "causatr_snm_no_em"
  )

  set.seed(1)
  res <- contrast(fit, ci_method = "bootstrap", n_boot = 50)
  expect_s3_class(res, "causatr_result")
  expect_equal(res$ci_method, "bootstrap")
  expect_true(all(is.finite(res$estimates$se)))
})
