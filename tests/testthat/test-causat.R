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

test_that("causat() rejects ATT for continuous treatment", {
  df <- data.frame(Y = rnorm(10), A = rnorm(10), L = rnorm(10))
  expect_snapshot(
    error = TRUE,
    causat(df, outcome = "Y", treatment = "A", confounders = ~L,
           estimand = "ATT")
  )
})

test_that("causat() rejects ATT for multivariate treatment", {
  df <- data.frame(Y = c(0, 1), A1 = c(0, 1), A2 = c(1, 0), L = c(1, 2))
  expect_snapshot(
    error = TRUE,
    causat(df, outcome = "Y", treatment = c("A1", "A2"), confounders = ~L,
           estimand = "ATT")
  )
})

test_that("causat() aborts when treatment has NAs and no censoring", {
  df <- data.frame(Y = c(0, 1, 0), A = c(0, NA, 1), L = c(1, 2, 3))
  expect_snapshot(
    error = TRUE,
    causat(df, outcome = "Y", treatment = "A", confounders = ~L)
  )
})

test_that("causat() rejects missing confounders_tv column", {
  df <- data.frame(
    Y = c(0, 1, 0, 1), A = c(0, 1, 0, 1), L = c(1, 1, 2, 2),
    id = c(1, 1, 2, 2), time = c(0, 1, 0, 1)
  )
  expect_snapshot(
    error = TRUE,
    causat(df, outcome = "Y", treatment = "A", confounders = ~L,
           confounders_tv = ~CD4, id = "id", time = "time")
  )
})

test_that("causat() rejects invalid history value", {
  df <- data.frame(
    Y = c(0, 1, 0, 1), A = c(0, 1, 0, 1), L = c(1, 1, 2, 2),
    id = c(1, 1, 2, 2), time = c(0, 1, 0, 1)
  )
  expect_snapshot(
    error = TRUE,
    causat(df, outcome = "Y", treatment = "A", confounders = ~L,
           id = "id", time = "time", history = 0)
  )
})
