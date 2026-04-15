# Regression tests for the 2026-04-15 critical review (B1, B2, B5, B6, B7,
# B8, R3, R6, R7, R8, R12). Each block asserts the specific failure mode
# documented in the review and is expected to FAIL on the pre-fix code.

test_that("B1: subset expression resolves session-scoped variables", {
  # Pre-fix: `parent.frame()` inside compute_contrast() was the internal
  # dispatch frame, not the caller's; `cutoff` would fail to resolve.
  df <- simulate_binary_continuous(n = 500, seed = 1L)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "gcomp"
  )
  cutoff <- 0
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    subset = quote(L > cutoff)
  )
  expect_s3_class(res, "causatr_result")
  # Point estimate on the `L > 0` subgroup is still ~3 (constant effect).
  ate <- res$contrasts$estimate[1]
  expect_true(abs(ate - 3) < 0.5)
})

test_that("B1: subset works inside bootstrap workers too", {
  df <- simulate_binary_continuous(n = 400, seed = 2L)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "gcomp"
  )
  cutoff <- -0.5
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "bootstrap",
    n_boot = 50L,
    subset = quote(L > cutoff)
  )
  expect_s3_class(res, "causatr_result")
  expect_false(is.na(res$contrasts$se[1]))
})

test_that("B1: subset length mismatch aborts with a clear message", {
  df <- simulate_binary_continuous(n = 200, seed = 3L)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "gcomp"
  )
  # A scalar `TRUE` would historically recycle silently to nrow(data).
  expect_error(
    contrast(
      fit,
      interventions = list(a1 = static(1), a0 = static(0)),
      type = "difference",
      subset = quote(TRUE)
    ),
    regexp = "length",
    fixed = FALSE
  )
})

test_that("B2: bootstrap refit replays user's ... (gcomp quasipoisson)", {
  # Pre-fix: refit_gcomp dropped user `...`, so a gcomp fit with a
  # non-default family silently bootstrapped with the default family.
  df <- simulate_binary_continuous(n = 400, seed = 4L)
  df$Y <- pmax(0, round(df$Y + 5))  # positive count-like outcome
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "gcomp",
    family = stats::quasipoisson()
  )
  # Bootstrap should run without aborting on family/response mismatch.
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "ratio",
    ci_method = "bootstrap",
    n_boot = 30L
  )
  expect_s3_class(res, "causatr_result")
  expect_true(all(res$estimates$estimate > 0))
})

test_that("B2: IPW bootstrap replays stashed WeightIt dots", {
  df <- simulate_binary_continuous(n = 500, seed = 5L)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw",
    method = "glm",
    stabilize = TRUE
  )
  expect_true(!is.null(fit$details$dots$stabilize))
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "bootstrap",
    n_boot = 25L
  )
  expect_s3_class(res, "causatr_result")
})

test_that("B5: causat_survival drops all rows at/after first censor", {
  # Build a tiny long panel where id=1 is censored at t=2.
  d <- data.table::data.table(
    id = rep(1:3, each = 3),
    t = rep(1:3, 3),
    Y = c(0, 0, 0,  0, 0, 1,  0, 0, 0),
    A = c(1, 1, 1,  0, 0, 0,  1, 1, 1),
    L = rnorm(9),
    C = c(0, 1, 0,  0, 0, 0,  0, 0, 0)
  )
  fit <- suppressMessages(causat_survival(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    id = "id",
    time = "t",
    censoring = "C",
    time_formula = ~ factor(t)
  ))
  # The row at (id=1, t=3) must NOT be in the fit set: id=1 was censored
  # at t=2 so subsequent rows must be dropped. Pre-fix only the t=2 row
  # was excluded.
  fit_rows <- fit$details$fit_rows
  row_at_1_3 <- which(d$id == 1 & d$t == 3)
  expect_false(fit_rows[row_at_1_3])
})

test_that("B6: external weights enter WeightIt via s.weights", {
  df <- simulate_binary_continuous(n = 400, seed = 6L)
  w <- runif(nrow(df), 0.5, 2)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw",
    method = "glm",
    weights = w
  )
  # Survey weights must now be attached to the weightit object as
  # s.weights rather than silently post-multiplied into w$weights.
  expect_true(!is.null(fit$weights_obj$s.weights))
  expect_equal(
    length(fit$weights_obj$s.weights),
    sum(fit$details$fit_rows)
  )
})

test_that("B7: ICE Ch1 IF uniform-weighted matches unweighted (unification)", {
  # With uniform weights, the weighted and unweighted Ch1 formulas must
  # agree exactly. Pre-fix the two branches used (n/n_target) and
  # n*(w/sum_w) which drift unless sum(w) == n_target.
  set.seed(42)
  n <- 150
  id <- rep(seq_len(n), each = 2)
  t <- rep(1:2, n)
  L0 <- rnorm(n)[id]
  L <- rnorm(2 * n) + L0
  A <- rbinom(2 * n, 1, plogis(0.5 * L))
  Y <- 1 + 0.8 * A + 0.5 * L + rnorm(2 * n)
  d <- data.frame(id = id, t = t, A = A, L = L, L0 = L0, Y = Y)
  # Only the last-time Y is the outcome; earlier Ys are intermediate.
  args <- list(
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    id = "id",
    time = "t",
    estimator = "gcomp",
    type = "longitudinal"
  )
  fit_u <- do.call(causat, c(list(data = d), args))
  fit_w <- do.call(causat, c(list(data = d, weights = rep(1, nrow(d))), args))
  r_u <- contrast(
    fit_u,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference"
  )
  r_w <- contrast(
    fit_w,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference"
  )
  expect_equal(r_u$contrasts$se[1], r_w$contrasts$se[1], tolerance = 1e-6)
})

test_that("B8: `by` skips empty strata instead of aborting", {
  df <- simulate_binary_continuous(n = 400, seed = 8L)
  df$g <- sample(c("a", "b", "c"), nrow(df), replace = TRUE)
  # Knock out `c` from the outcome model fit rows by setting NA outcomes.
  df$Y[df$g == "c"] <- NA
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "gcomp"
  )
  # Expect a warning about skipping level `c`, not a hard abort.
  expect_warning(
    res <- contrast(
      fit,
      interventions = list(a1 = static(1), a0 = static(0)),
      type = "difference",
      by = "g"
    ),
    regexp = "Skipped",
    fixed = FALSE
  )
  expect_s3_class(res, "causatr_result")
  expect_true(all(c("a", "b") %in% res$estimates$by))
  expect_false("c" %in% res$estimates$by)
})

test_that("R6: OR validation does not abort on NA mu_hat", {
  df <- simulate_binary_continuous(n = 200, seed = 9L)
  df$Y <- as.integer(df$Y > 3)  # binary
  df$Y[1:10] <- NA
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "gcomp",
    family = stats::binomial()
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "or"
  )
  expect_s3_class(res, "causatr_result")
})

test_that("R12: reserved column name is rejected up front", {
  df <- simulate_binary_continuous(n = 100, seed = 10L)
  df$.pseudo_y <- 0
  expect_error(
    causat(
      df,
      outcome = "Y",
      treatment = "A",
      confounders = ~L,
      estimator = "gcomp"
    ),
    regexp = "reserved",
    fixed = FALSE
  )
})
