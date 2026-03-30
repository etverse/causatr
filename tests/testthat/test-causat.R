test_that("causat() rejects missing outcome column", {
  df <- data.frame(A = c(0, 1), L = c(1, 2))
  expect_snapshot(
    error = TRUE,
    causat(df, outcome = "Y", treatment = "A", confounders = ~L)
  )
})

test_that("causat() rejects missing treatment column", {
  df <- data.frame(Y = c(0, 1), L = c(1, 2))
  expect_snapshot(
    error = TRUE,
    causat(df, outcome = "Y", treatment = "A", confounders = ~L)
  )
})

test_that("causat() rejects missing confounder column", {
  df <- data.frame(Y = c(0, 1), A = c(0, 1))
  expect_snapshot(
    error = TRUE,
    causat(df, outcome = "Y", treatment = "A", confounders = ~L)
  )
})

test_that("causat() rejects id without time", {
  df <- data.frame(Y = c(0, 1), A = c(0, 1), L = c(1, 2), id = c(1, 2))
  expect_snapshot(
    error = TRUE,
    causat(df, outcome = "Y", treatment = "A", confounders = ~L, id = "id")
  )
})

test_that("causat() rejects time without id", {
  df <- data.frame(Y = c(0, 1), A = c(0, 1), L = c(1, 2), t = c(0, 1))
  expect_snapshot(
    error = TRUE,
    causat(df, outcome = "Y", treatment = "A", confounders = ~L, time = "t")
  )
})

test_that("causat() errors on unimplemented gcomp", {
  df <- data.frame(Y = c(0, 1), A = c(0, 1), L = c(1, 2))
  expect_snapshot(
    error = TRUE,
    causat(df, outcome = "Y", treatment = "A", confounders = ~L)
  )
})
