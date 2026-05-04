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

test_that("longitudinal MAR (DGP-M5): ICE+IPCW near truth", {
  d <- simulate_longitudinal_mar_outcome(n = 5000, seed = 500)
  dt <- data.table::as.data.table(d)

  fit <- causat(
    dt,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    estimator = "gcomp",
    type = "longitudinal",
    id = "id",
    time = "time",
    censoring = "C",
    ipcw = TRUE
  )

  r <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "sandwich"
  )

  expect_equal(-r$contrasts$estimate, 5, tolerance = 0.15)
})


test_that("longitudinal MAR (DGP-M5): IPW+IPCW near truth", {
  d <- simulate_longitudinal_mar_outcome(n = 5000, seed = 501)
  dt <- data.table::as.data.table(d)

  fit <- causat(
    dt,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    estimator = "ipw",
    type = "longitudinal",
    id = "id",
    time = "time",
    censoring = "C",
    ipcw = TRUE
  )

  r <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "sandwich"
  )

  expect_equal(-r$contrasts$estimate, 5, tolerance = 0.25)
})


# ── Longitudinal lmtp cross-checks (DGP-M5) ─────────────────────

# Helper: convert DGP-M5 long format to lmtp wide format.
# causatr: long (one row per person-period), C=1 censored.
# lmtp:    wide (one row per person), C=1 observed, NAs after censoring.
dgp_m5_to_lmtp_wide <- function(d) {
  dt <- data.table::as.data.table(d)
  t0 <- as.data.frame(dt[dt$time == 0, c("id", "L0", "A")])
  names(t0)[3] <- "A_0"
  t1 <- as.data.frame(dt[dt$time == 1, c("id", "L", "A", "C", "Y")])
  names(t1) <- c("id", "L_1", "A_1", "C_causatr", "Y")
  wide <- merge(t0, t1, by = "id")
  # lmtp convention: C=1 observed. No censoring at t=0; censoring
  # happens between A_1 and Y (outcome not observed).
  wide$C_0 <- 1L
  wide$C_1 <- 1L - wide$C_causatr
  wide$Y[wide$C_1 == 0L] <- NA_real_
  wide$C_causatr <- NULL
  wide[, c("L0", "A_0", "C_0", "L_1", "A_1", "C_1", "Y")]
}


test_that("longitudinal lmtp cross-check: ICE+IPCW agrees with lmtp_sdr", {
  skip_if_not_installed("lmtp")

  d <- simulate_longitudinal_mar_outcome(n = 5000, seed = 500)
  dt <- data.table::as.data.table(d)

  # causatr ICE + IPCW
  fit <- causat(
    dt,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    estimator = "gcomp",
    type = "longitudinal",
    id = "id",
    time = "time",
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

  # lmtp SDR with censoring
  df <- dgp_m5_to_lmtp_wide(d)
  lmtp_1 <- tryCatch(
    suppressWarnings(suppressMessages(lmtp::lmtp_sdr(
      data = df,
      trt = c("A_0", "A_1"),
      outcome = "Y",
      baseline = "L0",
      time_vary = list(NULL, "L_1"),
      cens = c("C_0", "C_1"),
      shift = function(data, trt) 1,
      outcome_type = "continuous",
      learners_outcome = "SL.glm",
      learners_trt = "SL.glm",
      folds = 1
    ))),
    error = function(e) NULL
  )
  lmtp_0 <- tryCatch(
    suppressWarnings(suppressMessages(lmtp::lmtp_sdr(
      data = df,
      trt = c("A_0", "A_1"),
      outcome = "Y",
      baseline = "L0",
      time_vary = list(NULL, "L_1"),
      cens = c("C_0", "C_1"),
      shift = function(data, trt) 0,
      outcome_type = "continuous",
      learners_outcome = "SL.glm",
      learners_trt = "SL.glm",
      folds = 1
    ))),
    error = function(e) NULL
  )

  skip_if(
    is.null(lmtp_1) || is.null(lmtp_0),
    "lmtp::lmtp_sdr() failed on longitudinal censored data"
  )

  est_lmtp <- lmtp_1$estimate@x - lmtp_0$estimate@x

  expect_equal(est_causatr, 5, tolerance = 0.15)
  expect_equal(est_lmtp, 5, tolerance = 0.15)
  expect_lt(abs(est_causatr - est_lmtp), 0.15)
})


test_that("longitudinal lmtp cross-check: IPW+IPCW agrees with lmtp_sdr", {
  skip_if_not_installed("lmtp")

  d <- simulate_longitudinal_mar_outcome(n = 5000, seed = 501)
  dt <- data.table::as.data.table(d)

  # causatr IPW + IPCW
  fit <- causat(
    dt,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    estimator = "ipw",
    type = "longitudinal",
    id = "id",
    time = "time",
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

  # lmtp SDR with censoring
  df <- dgp_m5_to_lmtp_wide(d)
  lmtp_1 <- tryCatch(
    suppressWarnings(suppressMessages(lmtp::lmtp_sdr(
      data = df,
      trt = c("A_0", "A_1"),
      outcome = "Y",
      baseline = "L0",
      time_vary = list(NULL, "L_1"),
      cens = c("C_0", "C_1"),
      shift = function(data, trt) 1,
      outcome_type = "continuous",
      learners_outcome = "SL.glm",
      learners_trt = "SL.glm",
      folds = 1
    ))),
    error = function(e) NULL
  )
  lmtp_0 <- tryCatch(
    suppressWarnings(suppressMessages(lmtp::lmtp_sdr(
      data = df,
      trt = c("A_0", "A_1"),
      outcome = "Y",
      baseline = "L0",
      time_vary = list(NULL, "L_1"),
      cens = c("C_0", "C_1"),
      shift = function(data, trt) 0,
      outcome_type = "continuous",
      learners_outcome = "SL.glm",
      learners_trt = "SL.glm",
      folds = 1
    ))),
    error = function(e) NULL
  )

  skip_if(
    is.null(lmtp_1) || is.null(lmtp_0),
    "lmtp::lmtp_sdr() failed on longitudinal censored data"
  )

  est_lmtp <- lmtp_1$estimate@x - lmtp_0$estimate@x

  # Wider tolerance for IPW (higher variance)
  expect_equal(est_causatr, 5, tolerance = 0.30)
  expect_equal(est_lmtp, 5, tolerance = 0.30)
  expect_lt(abs(est_causatr - est_lmtp), 0.30)
})
