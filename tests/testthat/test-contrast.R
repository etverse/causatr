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
