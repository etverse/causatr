test_that("g-computation for point treatments is not yet implemented", {
  df <- data.frame(Y = c(0, 1), A = c(0, 1), L = c(1, 2))
  expect_snapshot(
    error = TRUE,
    causat(
      df,
      outcome = "Y",
      treatment = "A",
      confounders = ~L,
      method = "gcomp"
    )
  )
})
