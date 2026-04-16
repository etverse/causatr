test_that("causat(estimator = 'ipw') fits on simple data", {
  df <- data.frame(
    Y = c(1, 2, 3, 4, 5, 6, 7, 8),
    A = c(0, 0, 0, 0, 1, 1, 1, 1),
    L = c(1, 2, 3, 4, 1, 2, 3, 4)
  )
  fit <- suppressWarnings(causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw"
  ))
  expect_s3_class(fit, "causatr_fit")
  expect_equal(fit$estimator, "ipw")
  # Self-contained IPW: the treatment density model lives in `$details`;
  # `$weights_obj` is NULL and `$model` is only a display placeholder.
  expect_true(!is.null(fit$details$treatment_model))
  expect_true(!is.null(fit$details$propensity_model))
  expect_null(fit$weights_obj)
})

test_that("IPW with external weights stores in details", {
  d <- simulate_binary_continuous(n = 500, seed = 42)
  w <- runif(nrow(d), 0.5, 2)
  # `suppressWarnings()` here silences GLM's "non-integer #successes
  # in a binomial glm" warning, which fires whenever external
  # continuous-valued weights enter a logistic propensity fit. It is
  # harmless (the GLM fits exactly the same quasibinomial-style
  # weighted score equations) and appears for every weighted IPW
  # fit, so we keep the signal in the test clean.
  fit <- suppressWarnings(causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw",
    weights = w
  ))
  expect_true(!is.null(fit$details$weights))
  expect_equal(length(fit$details$weights), nrow(d))
  expect_false(".causatr_w" %in% names(fit$data))
})

test_that("IPW bootstrap with external weights gives finite SE", {
  d <- simulate_binary_continuous(n = 500, seed = 42)
  w <- runif(nrow(d), 0.5, 2)
  fit <- suppressWarnings(causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw",
    weights = w
  ))
  res <- suppressWarnings(contrast(
    fit,
    interventions = list(a0 = static(0), a1 = static(1)),
    ci_method = "bootstrap",
    n_boot = 50
  ))
  expect_true(all(is.finite(res$contrasts$se)))
  expect_true(res$contrasts$se > 0)
})

test_that("IPW continuous treatment fits and contrasts", {
  d <- simulate_continuous_continuous(n = 2000, seed = 42)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw"
  )
  expect_s3_class(fit, "causatr_fit")
  expect_equal(fit$estimator, "ipw")
  expect_equal(fit$details$treatment_model$family, "gaussian")
})

test_that("causat(estimator = 'ipw') rejects longitudinal data", {
  df <- data.frame(
    Y = c(0, 1, 0, 1),
    A = c(0, 1, 0, 1),
    L = c(1, 1, 2, 2),
    id = c(1, 1, 2, 2),
    time = c(0, 1, 0, 1)
  )
  expect_error(
    causat(
      df,
      outcome = "Y",
      treatment = "A",
      confounders = ~L,
      estimator = "ipw",
      id = "id",
      time = "time"
    ),
    "longitudinal"
  )
})

test_that("IPW with propensity_model_fn = mgcv::gam fits", {
  skip_if_not_installed("mgcv")
  d <- simulate_binary_continuous(n = 400, seed = 11)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~ s(L),
    estimator = "ipw",
    propensity_model_fn = mgcv::gam
  )
  expect_s3_class(fit, "causatr_fit")
  expect_true(inherits(fit$details$propensity_model, "gam"))
})

test_that("IPW rejects static(1) on a continuous treatment", {
  d <- simulate_continuous_continuous(n = 400, seed = 1)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw"
  )
  expect_snapshot(
    error = TRUE,
    contrast(
      fit,
      interventions = list(a = static(1), b = static(0)),
      ci_method = "sandwich"
    )
  )
})

test_that("IPW rejects A:modifier interaction terms in confounders", {
  # EM terms are detected but not yet supported under IPW. The error
  # class will change from `causatr_em_unsupported` to supported
  # behavior once chunk 6b lands.
  d <- simulate_binary_continuous(n = 200, seed = 1)
  d$sex <- rbinom(nrow(d), 1, 0.5)
  expect_error(
    causat(
      d,
      outcome = "Y",
      treatment = "A",
      confounders = ~ L + sex + A:sex,
      estimator = "ipw"
    ),
    class = "causatr_em_unsupported"
  )
})

test_that("IPW rejects bare treatment in confounders", {
  d <- simulate_binary_continuous(n = 200, seed = 2)
  expect_error(
    causat(
      d,
      outcome = "Y",
      treatment = "A",
      confounders = ~ L + A,
      estimator = "ipw"
    ),
    class = "causatr_bare_treatment_in_confounders"
  )
})
