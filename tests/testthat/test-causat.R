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

test_that("causat() rejects NA / non-finite / negative / mis-sized weights", {
  df <- data.frame(Y = stats::rnorm(50), A = stats::rbinom(50, 1, 0.5), L = stats::rnorm(50))

  expect_snapshot(
    error = TRUE,
    causat(df, outcome = "Y", treatment = "A", confounders = ~L,
           weights = c(NA, rep(1, 49)))
  )
  expect_snapshot(
    error = TRUE,
    causat(df, outcome = "Y", treatment = "A", confounders = ~L,
           weights = c(Inf, rep(1, 49)))
  )
  expect_snapshot(
    error = TRUE,
    causat(df, outcome = "Y", treatment = "A", confounders = ~L,
           weights = c(-1, rep(1, 49)))
  )
  expect_snapshot(
    error = TRUE,
    causat(df, outcome = "Y", treatment = "A", confounders = ~L,
           weights = rep(1, 40))
  )
  expect_snapshot(
    error = TRUE,
    causat(df, outcome = "Y", treatment = "A", confounders = ~L,
           weights = as.character(rep(1, 50)))
  )
})

test_that("causat() successfully fits a gcomp model", {
  df <- data.frame(
    Y = c(1, 2, 3, 4, 5, 6, 7, 8),
    A = c(0, 0, 0, 0, 1, 1, 1, 1),
    L = c(1, 2, 3, 4, 1, 2, 3, 4)
  )
  fit <- causat(df, outcome = "Y", treatment = "A", confounders = ~L)
  expect_s3_class(fit, "causatr_fit")
  expect_equal(fit$method, "gcomp")
})

test_that("causat() rejects ATT for continuous treatment", {
  df <- data.frame(Y = rnorm(10), A = rnorm(10), L = rnorm(10))
  expect_snapshot(
    error = TRUE,
    causat(
      df,
      outcome = "Y",
      treatment = "A",
      confounders = ~L,
      estimand = "ATT"
    )
  )
})

test_that("causat() rejects ATT for multivariate treatment", {
  df <- data.frame(Y = c(0, 1), A1 = c(0, 1), A2 = c(1, 0), L = c(1, 2))
  expect_snapshot(
    error = TRUE,
    causat(
      df,
      outcome = "Y",
      treatment = c("A1", "A2"),
      confounders = ~L,
      estimand = "ATT"
    )
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
    Y = c(0, 1, 0, 1),
    A = c(0, 1, 0, 1),
    L = c(1, 1, 2, 2),
    id = c(1, 1, 2, 2),
    time = c(0, 1, 0, 1)
  )
  expect_snapshot(
    error = TRUE,
    causat(
      df,
      outcome = "Y",
      treatment = "A",
      confounders = ~L,
      confounders_tv = ~CD4,
      id = "id",
      time = "time"
    )
  )
})

test_that("causat() rejects invalid history value", {
  df <- data.frame(
    Y = c(0, 1, 0, 1),
    A = c(0, 1, 0, 1),
    L = c(1, 1, 2, 2),
    id = c(1, 1, 2, 2),
    time = c(0, 1, 0, 1)
  )
  expect_snapshot(
    error = TRUE,
    causat(
      df,
      outcome = "Y",
      treatment = "A",
      confounders = ~L,
      id = "id",
      time = "time",
      history = 0
    )
  )
})


test_that("causat_survival() aborts on non-NULL competing= until Phase 6", {
  # Regression guard: previously `competing` was stored in
  # fit$details$competing but never consumed by the pooled-logistic
  # fit, silently fitting a cause-deleted hazard model. That gives
  # biased cumulative-incidence estimates whenever a competing event
  # is present — abort loudly until Phase 6 implements proper
  # sub-distribution / Aalen-Johansen handling.
  set.seed(1)
  long <- data.table::data.table(
    id = rep(1:5, each = 3),
    t = rep(0:2, times = 5),
    A = rep(c(0, 1, 0, 1, 0), each = 3),
    L = stats::rnorm(15),
    Y = c(0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1),
    cmp = 0
  )
  expect_snapshot(
    error = TRUE,
    causat_survival(
      long,
      outcome = "Y",
      treatment = "A",
      confounders = ~L,
      id = "id",
      time = "t",
      competing = "cmp"
    )
  )
})
