test_that("causat_mice() aborts when mice is not installed gracefully", {
  skip_if(requireNamespace("mice", quietly = TRUE), "mice is installed")
  expect_snapshot(
    error = TRUE,
    causat_mice(list(), outcome = "Y", treatment = "A",
                confounders = ~L, interventions = list(a1 = static(1)))
  )
})

test_that("causat_mice() rejects non-mids input", {
  skip_if_not_installed("mice")
  expect_snapshot(
    error = TRUE,
    causat_mice(list(), outcome = "Y", treatment = "A",
                confounders = ~L, interventions = list(a1 = static(1)))
  )
})

test_that("causat_mice() errors on unimplemented stub", {
  skip_if_not_installed("mice")
  imp <- structure(list(), class = "mids")
  expect_snapshot(
    error = TRUE,
    causat_mice(imp, outcome = "Y", treatment = "A",
                confounders = ~L, interventions = list(a1 = static(1)))
  )
})
