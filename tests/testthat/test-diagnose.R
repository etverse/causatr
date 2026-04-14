simulate_binary_continuous <- function(n = 2000, seed = 42) {
  set.seed(seed)
  L <- rnorm(n)
  ps <- plogis(0.5 * L)
  A <- rbinom(n, 1, ps)
  Y <- 2 + 3 * A + 1.5 * L + rnorm(n)
  data.frame(Y = Y, A = A, L = L)
}

# ============================================================
# BASIC STRUCTURE
# ============================================================

test_that("diagnose() rejects non-causatr_fit input", {
  expect_error(diagnose("not_a_fit"), "causatr_fit")
})

test_that("diagnose() rejects longitudinal fits with a clear error", {
  # Coverage gap: running diagnose() on a longitudinal fit used to
  # produce a degenerate single-row positivity table (logistic on all
  # time points pooled) and then crash inside cobalt's print method.
  # Reject explicitly until per-period diagnostics land in a future
  # phase. Snapshot locks the message so a future refactor can't
  # silently downgrade the error or accidentally re-enable the broken
  # code path.
  set.seed(1)
  n <- 60
  L0 <- stats::rnorm(n)
  A0 <- stats::rbinom(n, 1, plogis(0.4 * L0))
  L1 <- L0 + 0.3 * A0 + stats::rnorm(n)
  A1 <- stats::rbinom(n, 1, plogis(0.4 * L1 + 0.5 * A0))
  Y <- 1 + A0 + A1 + 0.5 * L1 + stats::rnorm(n)
  long <- data.table::data.table(
    id = rep(seq_len(n), each = 2),
    time = rep(0:1, times = n),
    A = as.numeric(rbind(A0, A1)),
    L = as.numeric(rbind(L0, L1)),
    Y = rep(Y, each = 2)
  )
  fit <- causat(
    long,
    outcome = "Y",
    treatment = "A",
    confounders = ~1,
    confounders_tv = ~L,
    id = "id",
    time = "time"
  )
  expect_snapshot(error = TRUE, diagnose(fit))
})

test_that("diagnose() returns causatr_diag for gcomp", {
  df <- simulate_binary_continuous(n = 500)
  fit <- causat(df, outcome = "Y", treatment = "A", confounders = ~L)
  diag <- diagnose(fit)

  expect_s3_class(diag, "causatr_diag")
  expect_equal(diag$method, "gcomp")
  expect_false(is.null(diag$positivity))
  expect_false(is.null(diag$balance))
  expect_null(diag$weights)
  expect_null(diag$match_quality)
})

test_that("diagnose() returns causatr_diag for IPW", {
  df <- simulate_binary_continuous(n = 500)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    method = "ipw"
  )
  diag <- diagnose(fit)

  expect_s3_class(diag, "causatr_diag")
  expect_equal(diag$method, "ipw")
  expect_false(is.null(diag$positivity))
  expect_false(is.null(diag$balance))
  expect_false(is.null(diag$weights))
  expect_null(diag$match_quality)
})

test_that("diagnose() returns causatr_diag for matching", {
  df <- simulate_binary_continuous(n = 500)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    method = "matching",
    estimand = "ATT"
  )
  diag <- diagnose(fit)

  expect_s3_class(diag, "causatr_diag")
  expect_equal(diag$method, "matching")
  expect_false(is.null(diag$positivity))
  expect_false(is.null(diag$balance))
  expect_null(diag$weights)
  expect_false(is.null(diag$match_quality))
})

# ============================================================
# POSITIVITY
# ============================================================

test_that("positivity table has expected rows for binary treatment", {
  df <- simulate_binary_continuous(n = 500)
  fit <- causat(df, outcome = "Y", treatment = "A", confounders = ~L)
  diag <- diagnose(fit)

  pos <- diag$positivity
  expect_s3_class(pos, "data.table")
  expect_true("statistic" %in% names(pos))
  expect_true("value" %in% names(pos))
  expect_true("min" %in% pos$statistic)
  expect_true("max" %in% pos$statistic)
  expect_true("n_violations" %in% pos$statistic)
})

test_that("positivity uses propensity scores from weightit for IPW", {
  df <- simulate_binary_continuous(n = 500)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    method = "ipw"
  )
  diag <- diagnose(fit)

  pos <- diag$positivity
  ps_min <- pos$value[pos$statistic == "min"]
  ps_max <- pos$value[pos$statistic == "max"]
  expect_gt(ps_min, 0)
  expect_lt(ps_max, 1)
})

# ============================================================
# BALANCE (cobalt)
# ============================================================

test_that("balance uses cobalt::bal.tab for IPW when available", {
  skip_if_not_installed("cobalt")
  df <- simulate_binary_continuous(n = 500)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    method = "ipw"
  )
  diag <- diagnose(fit)

  expect_s3_class(diag$balance, "bal.tab")
})

test_that("balance uses cobalt::bal.tab for matching when available", {
  skip_if_not_installed("cobalt")
  df <- simulate_binary_continuous(n = 500)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    method = "matching",
    estimand = "ATT"
  )
  diag <- diagnose(fit)

  expect_s3_class(diag$balance, "bal.tab")
})

test_that("balance uses cobalt::bal.tab for gcomp (unadjusted)", {
  skip_if_not_installed("cobalt")
  df <- simulate_binary_continuous(n = 500)
  fit <- causat(df, outcome = "Y", treatment = "A", confounders = ~L)
  diag <- diagnose(fit)

  expect_s3_class(diag$balance, "bal.tab")
})

# ============================================================
# WEIGHT DISTRIBUTION (IPW)
# ============================================================

test_that("weight summary has treated/control/overall for IPW", {
  df <- simulate_binary_continuous(n = 500)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    method = "ipw"
  )
  diag <- diagnose(fit)

  w <- diag$weights
  expect_s3_class(w, "data.table")
  expect_equal(nrow(w), 3L)
  expect_true(all(c("treated", "control", "overall") %in% w$group))
  expect_true("ess" %in% names(w))
  expect_true(all(w$ess > 0))
})

# ============================================================
# MATCH QUALITY (matching)
# ============================================================

test_that("match_quality has expected fields for matching", {
  df <- simulate_binary_continuous(n = 500)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    method = "matching",
    estimand = "ATT"
  )
  diag <- diagnose(fit)

  mq <- diag$match_quality
  expect_s3_class(mq, "data.table")
  expect_true("n_total" %in% mq$statistic)
  expect_true("n_matched" %in% mq$statistic)
  expect_true("n_discarded" %in% mq$statistic)
  expect_true("pct_retained" %in% mq$statistic)
})

# ============================================================
# PRINT AND SUMMARY
# ============================================================

test_that("print.causatr_diag outputs method info", {
  df <- simulate_binary_continuous(n = 500)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    method = "ipw"
  )
  diag <- diagnose(fit)

  expect_output(print(diag), "causatr_diag")
  expect_output(print(diag), "ipw")
})

test_that("summary.causatr_diag works", {
  df <- simulate_binary_continuous(n = 500)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    method = "ipw"
  )
  diag <- diagnose(fit)

  expect_output(summary(diag), "causatr_diag")
})

# ============================================================
# PLOT (cobalt love.plot)
# ============================================================

test_that("plot.causatr_diag produces a love plot for IPW", {
  skip_if_not_installed("cobalt")
  skip_if_not_installed("ggplot2")
  df <- simulate_binary_continuous(n = 500)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    method = "ipw"
  )
  diag <- diagnose(fit)

  p <- plot(diag)
  expect_true(inherits(p, "gg") || inherits(p, "gtable"))
})

test_that("plot.causatr_diag produces a love plot for matching", {
  skip_if_not_installed("cobalt")
  skip_if_not_installed("ggplot2")
  df <- simulate_binary_continuous(n = 500)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    method = "matching",
    estimand = "ATT"
  )
  diag <- diagnose(fit)

  p <- plot(diag)
  expect_true(inherits(p, "gg") || inherits(p, "gtable"))
})

test_that("plot.causatr_diag errors for gcomp (no love plot available)", {
  df <- simulate_binary_continuous(n = 500)
  fit <- causat(df, outcome = "Y", treatment = "A", confounders = ~L)
  diag <- diagnose(fit)

  expect_error(plot(diag), "IPW or matching")
})

# ============================================================
# CUSTOM PARAMETERS
# ============================================================

test_that("diagnose() respects custom thresholds", {
  skip_if_not_installed("cobalt")
  df <- simulate_binary_continuous(n = 500)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    method = "ipw"
  )
  diag <- diagnose(fit, thresholds = c(m = 0.05, v = 1.5))

  expect_s3_class(diag$balance, "bal.tab")
})

test_that("diagnose() respects custom ps_bounds", {
  df <- simulate_binary_continuous(n = 500)
  fit <- causat(df, outcome = "Y", treatment = "A", confounders = ~L)

  diag_wide <- diagnose(fit, ps_bounds = c(0.01, 0.99))
  diag_narrow <- diagnose(fit, ps_bounds = c(0.1, 0.9))

  n_viol_wide <- diag_wide$positivity$value[
    diag_wide$positivity$statistic == "n_violations"
  ]
  n_viol_narrow <- diag_narrow$positivity$value[
    diag_narrow$positivity$statistic == "n_violations"
  ]
  expect_lte(n_viol_wide, n_viol_narrow)
})
