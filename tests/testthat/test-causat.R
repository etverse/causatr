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
  df <- data.frame(
    Y = stats::rnorm(50),
    A = stats::rbinom(50, 1, 0.5),
    L = stats::rnorm(50)
  )

  expect_snapshot(
    error = TRUE,
    causat(
      df,
      outcome = "Y",
      treatment = "A",
      confounders = ~L,
      weights = c(NA, rep(1, 49))
    )
  )
  expect_snapshot(
    error = TRUE,
    causat(
      df,
      outcome = "Y",
      treatment = "A",
      confounders = ~L,
      weights = c(Inf, rep(1, 49))
    )
  )
  expect_snapshot(
    error = TRUE,
    causat(
      df,
      outcome = "Y",
      treatment = "A",
      confounders = ~L,
      weights = c(-1, rep(1, 49))
    )
  )
  expect_snapshot(
    error = TRUE,
    causat(
      df,
      outcome = "Y",
      treatment = "A",
      confounders = ~L,
      weights = rep(1, 40)
    )
  )
  expect_snapshot(
    error = TRUE,
    causat(
      df,
      outcome = "Y",
      treatment = "A",
      confounders = ~L,
      weights = as.character(rep(1, 50))
    )
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
  expect_equal(fit$estimator, "gcomp")
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


test_that("causat() rejects na.action = na.exclude forwarded through ...", {
  # Regression guard for the 2026-04-15 critical review Issue #7.
  # Under na.exclude, `residuals(model, 'working')` is padded with NAs
  # to the original data length, while `model.matrix(model)` drops NA
  # rows. `prepare_model_if()`'s `r_score <- residuals * weights`
  # then has a different length than `X_fit`, and
  # `apply_model_correction()`'s `d_fit * r_score` silently recycles
  # — R emits only a "longer object length is not a multiple" warning
  # and returns a mathematically wrong correction vector. Reproduced
  # numerically in /tmp/causatr_repro_issue7.R: baseline na.omit gave
  # SE = 0.0590, na.exclude triggered the recycling warning. Fix is
  # a hard abort at `check_dots_na_action()` in `causat()` /
  # `causat_survival()` — refuse `na.exclude` at the causatr boundary
  # rather than hardening every residuals() call site.
  set.seed(1)
  n <- 200
  L <- stats::rnorm(n)
  A <- stats::rbinom(n, 1, stats::plogis(0.3 * L))
  Y <- 1 + 0.5 * L + 0.8 * A + stats::rnorm(n, sd = 0.5)
  L[sample(n, 20)] <- NA
  d <- data.frame(L = L, A = A, Y = Y)

  # Function form of na.exclude.
  expect_error(
    causat(
      d,
      outcome = "Y",
      treatment = "A",
      confounders = ~L,
      estimator = "gcomp",
      na.action = stats::na.exclude
    ),
    class = "causatr_bad_na_action"
  )

  # String form — same rejection path.
  expect_error(
    causat(
      d,
      outcome = "Y",
      treatment = "A",
      confounders = ~L,
      estimator = "gcomp",
      na.action = "na.exclude"
    ),
    class = "causatr_bad_na_action"
  )

  # na.omit (the default) and na.fail must still be accepted.
  expect_no_error(
    causat(
      d,
      outcome = "Y",
      treatment = "A",
      confounders = ~L,
      estimator = "gcomp",
      na.action = stats::na.omit
    )
  )
})
