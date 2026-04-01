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
