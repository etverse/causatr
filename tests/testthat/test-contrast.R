test_that("contrast() rejects non-causatr_fit input", {
  expect_snapshot(error = TRUE, contrast(list(), list(a = static(1))))
})

test_that("contrast() rejects unnamed intervention list", {
  fit <- structure(list(method = "gcomp"), class = "causatr_fit")
  expect_snapshot(
    error = TRUE,
    contrast(fit, list(static(1), static(0)))
  )
})

test_that("contrast() rejects non-static interventions for IPW", {
  fit <- structure(list(method = "ipw"), class = "causatr_fit")
  expect_snapshot(
    error = TRUE,
    contrast(fit, list(a1 = shift(1), a0 = static(0)))
  )
})

test_that("contrast() rejects non-static interventions for matching", {
  fit <- structure(list(method = "matching"), class = "causatr_fit")
  expect_snapshot(
    error = TRUE,
    contrast(fit, list(a1 = dynamic(\(d, a) 1), a0 = static(0)))
  )
})

test_that("contrast() rejects estimand and subset together", {
  fit <- structure(
    list(method = "gcomp", estimand = "ATE"),
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
    list(method = "gcomp", estimand = "ATE"),
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
    list(method = "ipw", estimand = "ATE", treatment = "A", type = "point"),
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
      method = "matching",
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
      method = "gcomp",
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
    list(method = "gcomp", estimand = "ATE", treatment = "A", type = "point"),
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
    list(method = "gcomp", estimand = "ATE", treatment = "A", type = "point"),
    class = "causatr_fit"
  )
  expect_snapshot(
    error = TRUE,
    contrast(fit, list(a1 = list(static(1), static(0))))
  )
})
