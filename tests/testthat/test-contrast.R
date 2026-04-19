test_that("contrast() rejects non-causatr_fit input", {
  expect_snapshot(error = TRUE, contrast(list(), list(a = static(1))))
})

test_that("contrast() rejects unnamed intervention list", {
  fit <- structure(list(estimator = "gcomp"), class = "causatr_fit")
  expect_snapshot(
    error = TRUE,
    contrast(fit, list(static(1), static(0)))
  )
})

test_that("contrast() rejects duplicated intervention names", {
  fit <- structure(list(estimator = "gcomp"), class = "causatr_fit")
  expect_snapshot(
    error = TRUE,
    contrast(fit, list(a = static(1), a = static(0)))
  )
})

test_that("contrast() rejects an empty intervention list", {
  fit <- structure(list(estimator = "gcomp"), class = "causatr_fit")
  expect_snapshot(
    error = TRUE,
    contrast(fit, interventions = list())
  )
})

test_that("contrast() rejects non-static interventions for matching", {
  fit <- structure(list(estimator = "matching"), class = "causatr_fit")
  expect_snapshot(
    error = TRUE,
    contrast(fit, list(a1 = dynamic(\(d, a) 1), a0 = static(0)))
  )
})

test_that("contrast() rejects estimand and subset together", {
  fit <- structure(
    list(estimator = "gcomp", estimand = "ATE"),
    class = "causatr_fit"
  )
  expect_snapshot(
    error = TRUE,
    contrast(
      fit,
      list(a1 = static(1), a0 = static(0)),
      estimand = "ATT",
      subset = quote(A == 1)
    )
  )
})

test_that("contrast() rejects non-language subset", {
  fit <- structure(
    list(estimator = "gcomp", estimand = "ATE"),
    class = "causatr_fit"
  )
  expect_snapshot(
    error = TRUE,
    contrast(
      fit,
      list(a1 = static(1), a0 = static(0)),
      subset = "age > 50"
    )
  )
})

test_that("contrast() aborts when IPW estimand is changed", {
  fit <- structure(
    list(estimator = "ipw", estimand = "ATE", treatment = "A", type = "point"),
    class = "causatr_fit"
  )
  expect_snapshot(
    error = TRUE,
    contrast(fit, list(a1 = static(1), a0 = static(0)), estimand = "ATT")
  )
})

test_that("contrast() aborts when matching estimand is changed", {
  fit <- structure(
    list(
      estimator = "matching",
      estimand = "ATT",
      treatment = "A",
      type = "point"
    ),
    class = "causatr_fit"
  )
  expect_snapshot(
    error = TRUE,
    contrast(fit, list(a1 = static(1), a0 = static(0)), estimand = "ATE")
  )
})

test_that("contrast() rejects ATT for longitudinal fit", {
  fit <- structure(
    list(
      estimator = "gcomp",
      estimand = "ATE",
      treatment = "A",
      type = "longitudinal"
    ),
    class = "causatr_fit"
  )
  expect_snapshot(
    error = TRUE,
    contrast(fit, list(a1 = static(1), a0 = static(0)), estimand = "ATT")
  )
})

test_that("contrast() rejects reference not in interventions", {
  fit <- structure(
    list(
      estimator = "gcomp",
      estimand = "ATE",
      treatment = "A",
      type = "point"
    ),
    class = "causatr_fit"
  )
  expect_snapshot(
    error = TRUE,
    contrast(fit, list(a1 = static(1), a0 = static(0)), reference = "a2")
  )
})

test_that("contrast() handles NULL intervention (natural course)", {
  df <- data.frame(
    Y = rnorm(50),
    A = rnorm(50),
    L = rnorm(50)
  )
  fit <- causat(df, outcome = "Y", treatment = "A", confounders = ~L)
  result <- contrast(
    fit,
    list(shifted = shift(-1), observed = NULL),
    ci_method = "sandwich"
  )
  expect_s3_class(result, "causatr_result")
  expect_equal(nrow(result$estimates), 2L)
})

test_that("contrast() rejects multivariate intervention, missing trt var", {
  fit <- structure(
    list(
      estimator = "gcomp",
      estimand = "ATE",
      treatment = "A",
      type = "point"
    ),
    class = "causatr_fit"
  )
  expect_snapshot(
    error = TRUE,
    contrast(fit, list(a1 = list(static(1), static(0))))
  )
})

# ---- compute_contrast() validation aborts for ratio / OR ----------------
#
# When the user requests `type = "ratio"` or `"or"` on data that
# violates the log-scale delta method's positivity assumptions, the
# result would be NaN or -Inf CIs. The validator aborts up front
# instead of returning silently broken numbers. Each branch here pins
# one of those abort messages.

test_that("contrast() aborts on ratio when an alternative mu_hat is non-positive", {
  # Gaussian outcome with negative counterfactual mean: mu_hat <= 0 for
  # the static(0) arm of a centred-at-zero DGP. log(R) = log(neg) = NaN,
  # so ratio is undefined.
  set.seed(301)
  d <- data.frame(
    Y = rnorm(200, mean = -2),
    A = rbinom(200, 1, 0.5),
    L = rnorm(200)
  )
  fit <- causat(d, outcome = "Y", treatment = "A", confounders = ~L)
  expect_error(
    contrast(
      fit,
      interventions = list(a1 = static(1), a0 = static(0)),
      type = "ratio"
    ),
    "non-positive estimate"
  )
})

test_that("contrast() aborts on OR when a marginal mean is at the (0,1) boundary", {
  # Bernoulli outcome where every prediction lands at exactly 1.0 --
  # the OR validator catches it via either the `mu_ref` boundary check
  # or the `mu_hat >= 1` check (whichever arm trips first). Both are
  # legitimate paths that protect against NaN OR CIs.
  set.seed(302)
  n <- 100
  d <- data.frame(
    Y = rep(1, n),
    A = rbinom(n, 1, 0.5),
    L = rnorm(n)
  )
  fit <- suppressWarnings(causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    family = "binomial"
  ))
  expect_error(
    contrast(
      fit,
      interventions = list(a1 = static(1), a0 = static(0)),
      type = "or"
    ),
    "odds ratio is undefined|lie strictly in \\(0, 1\\)"
  )
})

test_that("contrast() aborts on empty target population", {
  # Force an empty target via a `subset` expression that never matches.
  set.seed(303)
  d <- simulate_binary_continuous(n = 100, seed = 303)
  fit <- causat(d, outcome = "Y", treatment = "A", confounders = ~L)
  expect_error(
    contrast(
      fit,
      interventions = list(a1 = static(1), a0 = static(0)),
      subset = quote(L > 100) # no rows satisfy this
    ),
    class = "causatr_empty_target"
  )
})
