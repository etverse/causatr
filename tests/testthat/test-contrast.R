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

test_that("contrast() allows NULL intervention (natural course)", {
  fit <- structure(
    list(method = "gcomp", estimand = "ATE", treatment = "A", type = "point"),
    class = "causatr_fit"
  )
  expect_snapshot(
    error = TRUE,
    contrast(fit, list(shifted = shift(-5), observed = NULL))
  )
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
