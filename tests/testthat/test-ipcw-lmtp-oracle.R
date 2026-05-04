# Cross-checks: causatr IPCW vs lmtp with censoring.
#
# lmtp::lmtp_sdr() handles censoring via its `cens =` argument. These
# tests verify that causatr's built-in IPCW produces estimates
# consistent with lmtp's doubly-robust estimator under the same DGP.
#
# Uses DGP-M2b (non-linear outcome with interaction, differential
# censoring) so that IPCW genuinely matters — a misspecified linear
# outcome model is biased under complete-case analysis but corrected
# by IPCW.
#
# These cross-checks apply to all IPCW chunks (14b, 14c) and should
# NOT be removed from the plan when individual chunks are completed.


# ── DGP-M2b: non-linear outcome with differential censoring ──────

test_that("IPCW eliminates censoring bias with misspecified model", {
  d <- simulate_mar_outcome_complex(n = 5000, seed = 400)
  dt <- data.table::as.data.table(d)

  # Misspecified model (linear, no interaction/quadratic)
  # WITH IPCW: should be near truth despite misspecification
  fit_ipcw <- causat(
    dt,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L1 + L2,
    estimator = "gcomp",
    censoring = "C",
    ipcw = TRUE
  )
  r_ipcw <- contrast(
    fit_ipcw,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "sandwich"
  )

  # WITHOUT IPCW (naive complete-case): biased
  fit_naive <- causat(
    dt,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L1 + L2,
    estimator = "gcomp",
    censoring = "C",
    ipcw = FALSE
  )
  r_naive <- contrast(
    fit_naive,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "sandwich"
  )

  ate_ipcw <- -r_ipcw$contrasts$estimate
  ate_naive <- -r_naive$contrasts$estimate

  # IPCW should be closer to truth (ATE = 3) than naive
  expect_lt(abs(ate_ipcw - 3), abs(ate_naive - 3))
  expect_equal(ate_ipcw, 3, tolerance = 0.15)
})


test_that("IPW without IPCW is biased; with IPCW is corrected", {
  d <- simulate_mar_outcome_complex(n = 5000, seed = 401)
  dt <- data.table::as.data.table(d)

  fit_ipcw <- causat(
    dt,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L1 + L2,
    estimator = "ipw",
    censoring = "C",
    ipcw = TRUE
  )
  r_ipcw <- contrast(
    fit_ipcw,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "sandwich"
  )

  fit_naive <- causat(
    dt,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L1 + L2,
    estimator = "ipw",
    censoring = "C",
    ipcw = FALSE
  )
  r_naive <- contrast(
    fit_naive,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "sandwich"
  )

  ate_ipcw <- -r_ipcw$contrasts$estimate
  ate_naive <- -r_naive$contrasts$estimate

  expect_lt(abs(ate_ipcw - 3), abs(ate_naive - 3))
  expect_equal(ate_ipcw, 3, tolerance = 0.20)
})


# ── lmtp cross-checks on DGP-M2b ─────────────────────────────────

test_that("gcomp+IPCW agrees with lmtp_sdr on DGP-M2b", {
  skip_if_not_installed("lmtp")

  d <- simulate_mar_outcome_complex(n = 5000, seed = 402)
  dt <- data.table::as.data.table(d)

  fit <- causat(
    dt,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L1 + L2,
    estimator = "gcomp",
    censoring = "C",
    ipcw = TRUE
  )
  r <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "sandwich"
  )
  est_causatr <- -r$contrasts$estimate

  # lmtp with censoring (C_lmtp = 1 - C: lmtp uses 1 = observed)
  d_lmtp <- d
  d_lmtp$C_lmtp <- 1L - d_lmtp$C

  lmtp_1 <- tryCatch(
    suppressWarnings(suppressMessages(lmtp::lmtp_sdr(
      data = d_lmtp,
      trt = "A",
      outcome = "Y",
      baseline = c("L1", "L2"),
      cens = "C_lmtp",
      shift = function(data, trt) 1,
      outcome_type = "continuous",
      mtp = FALSE,
      learners_outcome = "SL.glm",
      learners_trt = "SL.glm",
      folds = 1
    ))),
    error = function(e) NULL
  )
  lmtp_0 <- tryCatch(
    suppressWarnings(suppressMessages(lmtp::lmtp_sdr(
      data = d_lmtp,
      trt = "A",
      outcome = "Y",
      baseline = c("L1", "L2"),
      cens = "C_lmtp",
      shift = function(data, trt) 0,
      outcome_type = "continuous",
      mtp = FALSE,
      learners_outcome = "SL.glm",
      learners_trt = "SL.glm",
      folds = 1
    ))),
    error = function(e) NULL
  )

  skip_if(
    is.null(lmtp_1) || is.null(lmtp_0),
    "lmtp::lmtp_sdr() failed"
  )

  e1 <- tryCatch(lmtp_1$estimate@x, error = function(e) lmtp_1$theta)
  e0 <- tryCatch(lmtp_0$estimate@x, error = function(e) lmtp_0$theta)
  est_lmtp <- e1 - e0

  # Both near truth and near each other
  expect_equal(est_causatr, 3, tolerance = 0.15)
  expect_equal(est_lmtp, 3, tolerance = 0.15)
  expect_lt(abs(est_causatr - est_lmtp), 0.15)
})


test_that("IPW+IPCW agrees with lmtp_sdr on DGP-M2b", {
  skip_if_not_installed("lmtp")

  d <- simulate_mar_outcome_complex(n = 5000, seed = 403)
  dt <- data.table::as.data.table(d)

  fit <- causat(
    dt,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L1 + L2,
    estimator = "ipw",
    censoring = "C",
    ipcw = TRUE
  )
  r <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "sandwich"
  )
  est_causatr <- -r$contrasts$estimate

  d_lmtp <- d
  d_lmtp$C_lmtp <- 1L - d_lmtp$C

  lmtp_1 <- tryCatch(
    suppressWarnings(suppressMessages(lmtp::lmtp_sdr(
      data = d_lmtp,
      trt = "A",
      outcome = "Y",
      baseline = c("L1", "L2"),
      cens = "C_lmtp",
      shift = function(data, trt) 1,
      outcome_type = "continuous",
      mtp = FALSE,
      learners_outcome = "SL.glm",
      learners_trt = "SL.glm",
      folds = 1
    ))),
    error = function(e) NULL
  )
  lmtp_0 <- tryCatch(
    suppressWarnings(suppressMessages(lmtp::lmtp_sdr(
      data = d_lmtp,
      trt = "A",
      outcome = "Y",
      baseline = c("L1", "L2"),
      cens = "C_lmtp",
      shift = function(data, trt) 0,
      outcome_type = "continuous",
      mtp = FALSE,
      learners_outcome = "SL.glm",
      learners_trt = "SL.glm",
      folds = 1
    ))),
    error = function(e) NULL
  )

  skip_if(
    is.null(lmtp_1) || is.null(lmtp_0),
    "lmtp::lmtp_sdr() failed"
  )

  e1 <- tryCatch(lmtp_1$estimate@x, error = function(e) lmtp_1$theta)
  e0 <- tryCatch(lmtp_0$estimate@x, error = function(e) lmtp_0$theta)
  est_lmtp <- e1 - e0

  expect_equal(est_causatr, 3, tolerance = 0.20)
  expect_equal(est_lmtp, 3, tolerance = 0.20)
  expect_lt(abs(est_causatr - est_lmtp), 0.25)
})


# ── Also keep simple DGP-M2 truth checks ─────────────────────────

test_that("point MAR (DGP-M2): all estimators near truth", {
  d <- simulate_mar_outcome(n = 5000, seed = 404)
  dt <- data.table::as.data.table(d)

  fit_g <- causat(
    dt,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "gcomp",
    censoring = "C",
    ipcw = TRUE
  )
  r_g <- contrast(
    fit_g,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "sandwich"
  )

  fit_i <- causat(
    dt,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw",
    censoring = "C",
    ipcw = TRUE
  )
  r_i <- contrast(
    fit_i,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "sandwich"
  )

  expect_equal(-r_g$contrasts$estimate, 3, tolerance = 0.10)
  expect_equal(-r_i$contrasts$estimate, 3, tolerance = 0.15)
})


# ── Longitudinal MAR (DGP-M5) ────────────────────────────────────
# Requires chunk 14c (longitudinal IPCW). Enabled once 14c lands.

test_that("longitudinal MAR (DGP-M5): ICE+IPCW vs lmtp_sdr", {
  skip("awaiting chunk 14c: longitudinal IPCW")
  skip_if_not_installed("lmtp")
})
