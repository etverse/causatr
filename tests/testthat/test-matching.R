test_that("causat(method = 'matching') fits on simple data", {
  set.seed(1)
  df <- data.frame(
    Y = rnorm(100),
    A = rep(c(0, 1), 50),
    L = rnorm(100)
  )
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    method = "matching",
    estimand = "ATT"
  )
  expect_s3_class(fit, "causatr_fit")
  expect_equal(fit$method, "matching")
})

test_that("matching recovers ATT in the right direction", {
  d <- simulate_binary_continuous(n = 2000, seed = 42)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    method = "matching",
    estimand = "ATT"
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    type = "difference"
  )
  # 1:1 NN matching has O(n^{-1/p}) bias for continuous confounders
  # (Abadie & Imbens, 2006). True ATT = 3; matching is biased upward
  # but should be in the right ballpark.
  expect_equal(res$contrasts$estimate, 3, tolerance = 0.3)
  expect_gt(res$contrasts$se, 0)
})

test_that("matching with external weights stores in details", {
  d <- simulate_binary_continuous(n = 500, seed = 42)
  w <- runif(nrow(d), 0.5, 2)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    method = "matching",
    estimand = "ATT",
    weights = w
  )
  expect_true(!is.null(fit$details$weights))
  expect_false(".causatr_w" %in% names(fit$data))
})

test_that("matching bootstrap with external weights gives finite SE", {
  d <- simulate_binary_continuous(n = 500, seed = 42)
  w <- runif(nrow(d), 0.5, 2)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    method = "matching",
    estimand = "ATT",
    weights = w
  )
  res <- suppressWarnings(contrast(
    fit,
    interventions = list(a0 = static(0), a1 = static(1)),
    ci_method = "bootstrap",
    n_boot = 50
  ))
  expect_true(all(is.finite(res$contrasts$se)))
  expect_true(res$contrasts$se > 0)
})

test_that("causat(method = 'matching') rejects longitudinal data", {
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
      method = "matching",
      id = "id",
      time = "time"
    ),
    "longitudinal"
  )
})


test_that("matching aborts on continuous treatment", {
  # Locks in the rejection path (MatchIt itself errors). causatr
  # does not pre-check this; we rely on MatchIt's own validation,
  # but the FEATURE_COVERAGE_MATRIX requires this case to be
  # tested so we can detect a regression if a future MatchIt
  # version changes the behavior or if causatr ever adds
  # continuous-matching support without the gating check.
  set.seed(1)
  df <- data.frame(
    Y = stats::rnorm(200),
    A = stats::rnorm(200),
    L = stats::rnorm(200)
  )
  expect_error(
    causat(
      df,
      outcome = "Y",
      treatment = "A",
      confounders = ~L,
      method = "matching"
    )
  )
})
