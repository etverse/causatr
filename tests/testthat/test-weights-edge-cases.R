# Edge-case coverage for the external-weights path across all methods.
# Targets the check_weights() gate, the gcomp / IPW / matching / ICE
# weight plumbing, and the downstream Channel-1 / Mparts IF correction.
# Added after the 2026-04-15 critical-review sweep to cover every
# branch the review flagged as sensitive to weights.

test_that("check_weights rejects NA, Inf, negative, non-numeric, wrong length", {
  df <- simulate_binary_continuous(n = 50, seed = 1L)
  expect_error(
    causat(
      df,
      outcome = "Y",
      treatment = "A",
      confounders = ~L,
      estimator = "gcomp",
      weights = c(NA_real_, rep(1, nrow(df) - 1))
    ),
    regexp = "missing"
  )
  expect_error(
    causat(
      df,
      outcome = "Y",
      treatment = "A",
      confounders = ~L,
      estimator = "gcomp",
      weights = c(Inf, rep(1, nrow(df) - 1))
    ),
    regexp = "non-finite|Inf"
  )
  expect_error(
    causat(
      df,
      outcome = "Y",
      treatment = "A",
      confounders = ~L,
      estimator = "gcomp",
      weights = c(-1, rep(1, nrow(df) - 1))
    ),
    regexp = "non-negative"
  )
  expect_error(
    causat(
      df,
      outcome = "Y",
      treatment = "A",
      confounders = ~L,
      estimator = "gcomp",
      weights = rep("1", nrow(df))
    ),
    regexp = "numeric"
  )
  expect_error(
    causat(
      df,
      outcome = "Y",
      treatment = "A",
      confounders = ~L,
      estimator = "gcomp",
      weights = rep(1, nrow(df) - 1)
    ),
    regexp = "length"
  )
})

test_that("gcomp with zero weights treats those rows as excluded", {
  df <- simulate_binary_continuous(n = 300, seed = 2L)
  w <- rep(1, nrow(df))
  w[1:50] <- 0
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "gcomp",
    weights = w
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference"
  )
  # Re-fit dropping the zero-weight rows explicitly; point estimate
  # must agree. This is a truth test on the pass-through semantics:
  # zero weights must be equivalent to row exclusion.
  df2 <- df[-(1:50), ]
  fit2 <- causat(
    df2,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "gcomp"
  )
  res2 <- contrast(
    fit2,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference"
  )
  # Point estimates close; SEs differ slightly because zero-weight
  # rows still enter the target population for the marginal mean
  # under the full-data fit (different `n` in the IF).
  expect_equal(
    res$contrasts$estimate[1],
    res2$contrasts$estimate[1],
    tolerance = 0.15
  )
})

test_that("gcomp uniform weights give exactly the same SE as no weights", {
  df <- simulate_binary_continuous(n = 400, seed = 3L)
  fit_u <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "gcomp"
  )
  fit_w <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "gcomp",
    weights = rep(1, nrow(df))
  )
  r_u <- contrast(
    fit_u,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference"
  )
  r_w <- contrast(
    fit_w,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference"
  )
  expect_equal(
    r_u$contrasts$estimate[1],
    r_w$contrasts$estimate[1],
    tolerance = 1e-10
  )
  expect_equal(r_u$contrasts$se[1], r_w$contrasts$se[1], tolerance = 1e-8)
})

test_that("gcomp heterogeneous weights change SE (not a no-op)", {
  df <- simulate_binary_continuous(n = 400, seed = 4L)
  w <- runif(nrow(df), 0.1, 4)
  fit_u <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "gcomp"
  )
  fit_w <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "gcomp",
    weights = w
  )
  r_u <- contrast(
    fit_u,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference"
  )
  r_w <- contrast(
    fit_w,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference"
  )
  # Weighted SE should differ from unweighted SE by a non-trivial
  # amount — otherwise the weights aren't reaching Channel 1 or the
  # sandwich meat.
  expect_gt(abs(r_u$contrasts$se[1] - r_w$contrasts$se[1]), 1e-3)
})

test_that("IPW external weights reach the treatment density fit", {
  df <- simulate_binary_continuous(n = 400, seed = 5L)
  w <- runif(nrow(df), 0.5, 2)
  fit <- suppressWarnings(causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw",
    weights = w
  ))
  # Self-contained IPW: the external weights enter the propensity
  # GLM as prior.weights, so they are observable via the GLM's
  # own `$prior.weights` slot. They are also stashed on
  # `fit$details$weights` for Channel 1 inside `contrast()`.
  pw <- fit$details$propensity_model$prior.weights
  expect_equal(length(pw), sum(fit$details$fit_rows))
  expect_equal(as.numeric(pw), w[fit$details$fit_rows])
  expect_equal(fit$details$weights, w)
})

test_that("IPW survey-weighted sandwich SE is finite and positive", {
  # Regression: with external weights, the self-contained IPW
  # sandwich must propagate propensity uncertainty through the
  # cross-derivative. We check the sandwich SE is finite and
  # non-zero on a realistic weighted design.
  df <- simulate_binary_continuous(n = 500, seed = 6L)
  w <- runif(nrow(df), 0.5, 2)
  fit <- suppressWarnings(causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw",
    weights = w
  ))
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference"
  )
  expect_true(is.finite(res$contrasts$se[1]))
  expect_gt(res$contrasts$se[1], 0)
})

test_that("matching combines match weights with external weights correctly", {
  df <- simulate_binary_continuous(n = 400, seed = 7L)
  w <- runif(nrow(df), 0.5, 2)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "matching",
    weights = w
  )
  # The outcome model on the matched sample must carry non-trivial
  # weights (match weights * external weights).
  mw <- stats::weights(fit$model)
  expect_true(all(mw > 0))
  expect_false(all(mw == 1))
})

test_that("ICE with uniform weights matches ICE without weights (both SE and est)", {
  set.seed(8L)
  n <- 150
  id <- rep(seq_len(n), each = 2)
  t <- rep(1:2, n)
  L0 <- rnorm(n)[id]
  L <- rnorm(2 * n) + L0
  A <- rbinom(2 * n, 1, plogis(0.5 * L))
  Y <- 1 + 0.8 * A + 0.5 * L + rnorm(2 * n)
  d <- data.frame(id = id, t = t, A = A, L = L, L0 = L0, Y = Y)
  args <- list(
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    id = "id",
    time = "t",
    estimator = "gcomp",
    type = "longitudinal"
  )
  fit_u <- do.call(causat, c(list(data = d), args))
  fit_w <- do.call(causat, c(list(data = d, weights = rep(1, nrow(d))), args))
  r_u <- contrast(
    fit_u,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference"
  )
  r_w <- contrast(
    fit_w,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference"
  )
  expect_equal(
    r_u$contrasts$estimate[1],
    r_w$contrasts$estimate[1],
    tolerance = 1e-10
  )
  expect_equal(r_u$contrasts$se[1], r_w$contrasts$se[1], tolerance = 1e-10)
})

test_that("ICE heterogeneous weights produce finite, non-zero SE", {
  set.seed(9L)
  n <- 200
  id <- rep(seq_len(n), each = 2)
  t <- rep(1:2, n)
  L0 <- rnorm(n)[id]
  L <- rnorm(2 * n) + L0
  A <- rbinom(2 * n, 1, plogis(0.5 * L))
  Y <- 1 + 0.8 * A + 0.5 * L + rnorm(2 * n)
  d <- data.frame(id = id, t = t, A = A, L = L, L0 = L0, Y = Y)
  # Weight by id so both rows of each subject get the same weight.
  w_id <- runif(n, 0.5, 2)
  w <- w_id[id]
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    id = "id",
    time = "t",
    estimator = "gcomp",
    type = "longitudinal",
    weights = w
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference"
  )
  expect_true(is.finite(res$contrasts$se[1]))
  expect_gt(res$contrasts$se[1], 0)
  expect_true(is.finite(res$contrasts$estimate[1]))
})
