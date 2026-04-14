# S3 method tests for causatr_fit, causatr_result, causatr_diag.
# DGP functions are defined in helper-dgp.R (auto-loaded by testthat).

tol_bf <- 0.5


test_that("confint() respects level argument", {
  df <- simulate_effect_mod()
  fit <- causat(df, outcome = "Y", treatment = "A", confounders = ~ L + sex)
  result <- contrast(fit, list(a1 = static(1), a0 = static(0)))

  ci_95 <- confint(result, level = 0.95)
  ci_99 <- confint(result, level = 0.99)

  expect_true(ci_99[1, "lower"] < ci_95[1, "lower"])
  expect_true(ci_99[1, "upper"] > ci_95[1, "upper"])

  z95 <- qnorm(0.975)
  z99 <- qnorm(0.995)
  se <- result$estimates$se[1]
  est <- result$estimates$estimate[1]
  expect_equal(ci_95[1, "lower"], est - z95 * se)
  expect_equal(ci_99[1, "lower"], est - z99 * se)
})

test_that("confint() returns correct structure", {
  df <- simulate_effect_mod(n = 100)
  fit <- causat(df, outcome = "Y", treatment = "A", confounders = ~L)
  result <- contrast(fit, list(a1 = static(1), a0 = static(0)))

  ci <- confint(result)
  expect_equal(nrow(ci), 2L)
  expect_equal(colnames(ci), c("lower", "upper"))
  expect_equal(rownames(ci), c("a1", "a0"))
})

test_that("by parameter recovers subgroup-specific ATEs", {
  df <- simulate_effect_mod()
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L + sex + A:sex
  )

  result <- contrast(
    fit,
    list(a1 = static(1), a0 = static(0)),
    ci_method = "sandwich",
    by = "sex"
  )

  expect_s3_class(result, "causatr_result")
  expect_true("by" %in% names(result$estimates))
  expect_equal(sort(unique(result$estimates$by)), c("0", "1"))
  expect_equal(nrow(result$estimates), 4L)
  expect_equal(nrow(result$contrasts), 2L)

  ate_sex0 <- abs(result$contrasts$estimate[result$contrasts$by == "0"])
  ate_sex1 <- abs(result$contrasts$estimate[result$contrasts$by == "1"])
  expect_true(abs(ate_sex0 - 3.0) < tol_bf)
  expect_true(abs(ate_sex1 - 4.2) < tol_bf)

  expect_true("n_by" %in% names(result$estimates))
  expect_true("n_by" %in% names(result$contrasts))
  expect_true(all(result$estimates$n_by > 0))
  expect_equal(
    sum(unique(result$estimates$n_by)),
    result$n
  )
})

test_that("by parameter rejects missing variable", {
  df <- simulate_effect_mod(n = 50)
  fit <- causat(df, outcome = "Y", treatment = "A", confounders = ~L)
  expect_snapshot(
    error = TRUE,
    contrast(
      fit,
      list(a1 = static(1), a0 = static(0)),
      by = "nonexistent"
    )
  )
})

test_that("causat_survival() does not mutate input data.table", {
  dt <- data.table::data.table(
    id = rep(1:20, each = 3),
    time = rep(1:3, 20),
    event = rbinom(60, 1, 0.1),
    A = rep(rbinom(20, 1, 0.5), each = 3),
    L = rnorm(60)
  )
  original_names <- copy(names(dt))
  original_ncol <- ncol(dt)

  # Use `~ factor(time)` instead of the default `splines::ns(time, 4)`:
  # the test data has only 3 unique time points, which is below the
  # minimum needed for ns(df=4) and triggers a "shoving interior knots"
  # warning. factor() is the saturated baseline-hazard alternative
  # documented in causat_survival()'s `time_formula` parameter.
  fit <- causat_survival(
    dt,
    outcome = "event",
    treatment = "A",
    confounders = ~L,
    id = "id",
    time = "time",
    time_formula = ~ factor(time)
  )

  expect_equal(names(dt), original_names)
  expect_equal(ncol(dt), original_ncol)
  expect_false("prev_event" %in% names(dt))
})

test_that("survival type aborts in contrast()", {
  dt <- data.table::data.table(
    id = rep(1:20, each = 3),
    time = rep(1:3, 20),
    event = rbinom(60, 1, 0.1),
    A = rep(rbinom(20, 1, 0.5), each = 3),
    L = rnorm(60)
  )

  # Use `~ factor(time)` instead of the default `splines::ns(time, 4)`:
  # the test data has only 3 unique time points, which is below the
  # minimum needed for ns(df=4) and triggers a "shoving interior knots"
  # warning. factor() is the saturated baseline-hazard alternative
  # documented in causat_survival()'s `time_formula` parameter.
  fit <- causat_survival(
    dt,
    outcome = "event",
    treatment = "A",
    confounders = ~L,
    id = "id",
    time = "time",
    time_formula = ~ factor(time)
  )

  expect_snapshot(
    error = TRUE,
    contrast(fit, list(a1 = static(1), a0 = static(0)))
  )
})

test_that("multivariate treatment blocked for IPW", {
  expect_snapshot(
    error = TRUE,
    causat(
      data.frame(Y = 1:10, A1 = 1:10, A2 = 1:10, L = 1:10),
      outcome = "Y",
      treatment = c("A1", "A2"),
      confounders = ~L,
      method = "ipw"
    )
  )
})

test_that("multivariate treatment blocked for matching", {
  expect_snapshot(
    error = TRUE,
    causat(
      data.frame(Y = 1:10, A1 = 1:10, A2 = 1:10, L = 1:10),
      outcome = "Y",
      treatment = c("A1", "A2"),
      confounders = ~L,
      method = "matching"
    )
  )
})

test_that("vcov.causatr_result returns correct matrix", {
  df <- data.frame(
    Y = rnorm(100),
    A = rbinom(100, 1, 0.5),
    L = rnorm(100)
  )
  fit <- causat(df, outcome = "Y", treatment = "A", confounders = ~L)
  result <- contrast(
    fit,
    list(a1 = static(1), a0 = static(0)),
    ci_method = "sandwich"
  )

  v <- vcov(result)
  expect_true(is.matrix(v))
  expect_equal(nrow(v), 2L)
  expect_equal(ncol(v), 2L)
  expect_equal(rownames(v), c("a1", "a0"))
})

test_that("vcov.causatr_result returns per-stratum list when by is used", {
  df <- simulate_effect_mod()
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L + sex + A:sex
  )
  result <- contrast(
    fit,
    list(a1 = static(1), a0 = static(0)),
    ci_method = "sandwich",
    by = "sex"
  )

  v <- vcov(result)
  expect_true(is.list(v))
  expect_equal(length(v), 2L)
  expect_true(all(vapply(v, is.matrix, logical(1))))
  expect_equal(nrow(v[[1]]), 2L)
  expect_equal(ncol(v[[1]]), 2L)
})

test_that("tidy.causatr_result returns data frame", {
  df <- data.frame(
    Y = rnorm(100),
    A = rbinom(100, 1, 0.5),
    L = rnorm(100)
  )
  fit <- causat(df, outcome = "Y", treatment = "A", confounders = ~L)
  result <- contrast(
    fit,
    list(a1 = static(1), a0 = static(0)),
    ci_method = "sandwich"
  )

  td <- tidy(result)
  expect_s3_class(td, "data.frame")
  expect_true("term" %in% names(td))
  expect_true("estimate" %in% names(td))
  expect_true("std.error" %in% names(td))
  expect_true("conf.low" %in% names(td))
  expect_true("conf.high" %in% names(td))

  td_means <- tidy(result, which = "means")
  expect_equal(nrow(td_means), 2L)
  expect_true(all(td_means$type == "mean"))

  td_all <- tidy(result, which = "all")
  expect_equal(nrow(td_all), 3L)
})

test_that("glance.causatr_result returns one-row data frame", {
  df <- data.frame(
    Y = rnorm(100),
    A = rbinom(100, 1, 0.5),
    L = rnorm(100)
  )
  fit <- causat(df, outcome = "Y", treatment = "A", confounders = ~L)
  result <- contrast(
    fit,
    list(a1 = static(1), a0 = static(0)),
    ci_method = "sandwich"
  )

  gl <- glance(result)
  expect_s3_class(gl, "data.frame")
  expect_equal(nrow(gl), 1L)
  expect_equal(gl$method, "gcomp")
  expect_equal(gl$n_interventions, 2L)
})

test_that("print.causatr_fit shows enhanced output", {
  df <- data.frame(
    Y = rnorm(50),
    A = rbinom(50, 1, 0.5),
    L = rnorm(50)
  )
  fit <- causat(df, outcome = "Y", treatment = "A", confounders = ~L)

  out <- capture.output(print(fit))
  expect_true(any(grepl("G-computation", out)))
  expect_true(any(grepl("Estimand", out)))
  expect_true(any(grepl("gaussian", out)))
  expect_true(any(grepl("Confounders", out)))
})

test_that("print.causatr_result shows estimand", {
  df <- data.frame(
    Y = rnorm(50),
    A = rbinom(50, 1, 0.5),
    L = rnorm(50)
  )
  fit <- causat(df, outcome = "Y", treatment = "A", confounders = ~L)
  result <- contrast(
    fit,
    list(a1 = static(1), a0 = static(0)),
    ci_method = "sandwich"
  )

  out <- capture.output(print(result))
  expect_true(any(grepl("Estimand", out)))
  expect_true(any(grepl("ATE", out)))
})

test_that("weight summary handles continuous treatment IPW", {
  set.seed(42)
  df <- data.frame(
    Y = rnorm(100),
    A = rnorm(100),
    L = rnorm(100)
  )
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    method = "ipw"
  )
  diag <- diagnose(fit)

  expect_true(!is.null(diag$weights))
  expect_true(any(diag$weights$group == "overall"))
})

test_that("coef.causatr_result returns named numeric vector", {
  df <- simulate_effect_mod(n = 200)
  fit <- causat(df, outcome = "Y", treatment = "A", confounders = ~L)
  result <- contrast(
    fit,
    list(a1 = static(1), a0 = static(0)),
    ci_method = "sandwich"
  )

  cf <- coef(result)
  expect_true(is.numeric(cf))
  expect_equal(length(cf), 2L)
  expect_equal(names(cf), c("a1", "a0"))
  expect_equal(unname(cf), result$estimates$estimate)
})

test_that("plot.causatr_result runs without error", {
  skip_if_not_installed("forrest")
  df <- simulate_effect_mod(n = 200)
  fit <- causat(df, outcome = "Y", treatment = "A", confounders = ~L)
  result <- contrast(
    fit,
    list(a1 = static(1), a0 = static(0)),
    ci_method = "sandwich"
  )

  expect_no_error(plot(result))
  expect_no_error(plot(result, which = "means"))
})

test_that("summary.causatr_result includes vcov output", {
  df <- simulate_effect_mod(n = 200)
  fit <- causat(df, outcome = "Y", treatment = "A", confounders = ~L)
  result <- contrast(
    fit,
    list(a1 = static(1), a0 = static(0)),
    ci_method = "sandwich"
  )

  out <- capture.output(summary(result))
  expect_true(any(grepl("Variance-covariance", out)))
})

test_that("summary.causatr_fit includes model summary", {
  df <- simulate_effect_mod(n = 200)
  fit <- causat(df, outcome = "Y", treatment = "A", confounders = ~L)

  out <- capture.output(summary(fit))
  expect_true(any(
    grepl("Coefficients", out) | grepl("coef", out, ignore.case = TRUE)
  ))
})
