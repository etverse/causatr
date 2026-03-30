test_that("ICE g-computation is not yet implemented", {
  df <- data.frame(
    id = c(1, 1, 2, 2),
    time = c(0, 1, 0, 1),
    Y = c(0, 1, 0, 0),
    A = c(0, 1, 1, 1),
    L = c(1, 2, 3, 4)
  )
  expect_snapshot(
    error = TRUE,
    causat(
      df,
      outcome = "Y",
      treatment = "A",
      confounders = ~L,
      id = "id",
      time = "time"
    )
  )
})
