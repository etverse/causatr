# Longitudinal IPW (Phase 10 chunk 10a)
#
# Per-period density chain + cumulative product weight + final-period
# weighted Hajek MSM + stacked sandwich + id-clustered bootstrap.
#
# Truth-based tests use `make_continuous_scm()` and `make_linear_scm()`
# from helper-dgp.R, which both implement the canonical 2-period
# linear SCM with treatment-confounder feedback. Cross-method
# agreement against ICE g-computation is the primary correctness
# check; `lmtp::lmtp_sdr()` provides an external point-estimate
# oracle when installed.
#
# Rejection snapshots cover the deferred combinations (multivariate,
# stabilize, EM, IPSI, threshold) that ship in follow-up chunks.

# ----------------------------------------------------------------------
# T-long-ipw1: 2-period continuous DGP, shift recovers ICE point
# ----------------------------------------------------------------------
test_that("T-long-ipw1: longitudinal IPW shift agrees with ICE g-comp on linear-Gaussian DGP", {
  # `make_continuous_scm()` is the canonical 2-period continuous SCM
  # with treatment-confounder feedback (E[A_k] = 1, ATE_shift_per_unit
  # = 4 from gcomp). On the same data both methods should converge to
  # the same point estimate up to MC error.
  set.seed(7001)
  d <- make_continuous_scm(n = 1500, seed = 7001)
  d <- data.table::as.data.table(d)

  fit_ipw_long <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    estimator = "ipw"
  )
  res_ipw <- contrast(
    fit_ipw_long,
    interventions = list(shifted = shift(0.5), nat = NULL),
    type = "difference",
    reference = "nat",
    ci_method = "sandwich"
  )

  # Suppress the structural "dropping all-NA column" inform from ICE
  # at t = 0 (L is NA at the first period in this DGP); it's
  # informational, not a correctness signal.
  fit_ice <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    estimator = "gcomp"
  )
  res_ice <- suppressMessages(contrast(
    fit_ice,
    interventions = list(shifted = shift(0.5), nat = NULL),
    type = "difference",
    reference = "nat"
  ))

  # Cross-method agreement on the contrast point estimate. IPW SE is
  # typically larger than ICE SE under correct outcome-model
  # specification; we don't compare SEs across methods.
  est_ipw <- res_ipw$contrasts$estimate[1]
  est_ice <- res_ice$contrasts$estimate[1]
  expect_lt(abs(est_ipw - est_ice), 0.5)

  # Marginal means recoverable: nat course should match the observed
  # E[Y] (= 15.9 in this seed).
  expect_lt(
    abs(
      res_ipw$estimates$estimate[res_ipw$estimates$intervention == "nat"] -
        mean(d$Y, na.rm = TRUE)
    ),
    0.1
  )

  # Sandwich SE is finite and positive.
  expect_true(all(is.finite(res_ipw$contrasts$se) & res_ipw$contrasts$se > 0))
})


# ----------------------------------------------------------------------
# T-long-ipw2: binary treatment, static intervention recovers ICE
# ----------------------------------------------------------------------
test_that("T-long-ipw2: longitudinal IPW static recovers ICE on binary DGP", {
  # `make_linear_scm()` is the canonical 2-period binary-treatment SCM
  # with treatment-confounder feedback. Under always vs never the
  # gcomp truth is 5 (Hernan & Robins Ch. 21). Cross-method agreement
  # IPW vs ICE pins the longitudinal IPW pipeline.
  set.seed(7002)
  d <- make_linear_scm(n = 2000, seed = 7002)
  d <- data.table::as.data.table(d)

  fit_ipw_long <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    estimator = "ipw"
  )
  res_ipw <- contrast(
    fit_ipw_long,
    interventions = list(always = static(1), never = static(0)),
    type = "difference",
    reference = "never",
    ci_method = "sandwich"
  )

  fit_ice <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    estimator = "gcomp"
  )
  res_ice <- suppressMessages(contrast(
    fit_ice,
    interventions = list(always = static(1), never = static(0)),
    type = "difference",
    reference = "never"
  ))

  est_ipw <- res_ipw$contrasts$estimate[1]
  est_ice <- res_ice$contrasts$estimate[1]

  # Cross-method agreement. IPW under static can have heavier
  # variance from rare-cell HT weights; widen the tolerance compared
  # to the continuous case but still tighter than 1.
  expect_lt(abs(est_ipw - est_ice), 0.8)
  # IPW estimate is in the right ballpark of the truth (5).
  expect_lt(abs(est_ipw - 5), 1.2)
  expect_true(all(is.finite(res_ipw$contrasts$se) & res_ipw$contrasts$se > 0))
})


# ----------------------------------------------------------------------
# T-long-ipw3: bootstrap variance ≈ sandwich variance
# ----------------------------------------------------------------------
test_that("T-long-ipw3: longitudinal IPW bootstrap SE within MC tolerance of sandwich SE", {
  # Bootstrap resamples ids (entire person-period trajectories
  # together). With B = 200 we expect ~30% MC tolerance on the SE.
  set.seed(7003)
  d <- make_continuous_scm(n = 800, seed = 7003)
  d <- data.table::as.data.table(d)

  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    estimator = "ipw"
  )

  res_sand <- contrast(
    fit,
    interventions = list(s = shift(0.5), n = NULL),
    reference = "n",
    ci_method = "sandwich"
  )
  res_boot <- suppressWarnings(contrast(
    fit,
    interventions = list(s = shift(0.5), n = NULL),
    reference = "n",
    ci_method = "bootstrap",
    n_boot = 200
  ))

  se_sand <- res_sand$contrasts$se[1]
  se_boot <- res_boot$contrasts$se[1]
  expect_true(is.finite(se_boot) && se_boot > 0)
  ratio <- se_boot / se_sand
  # 0.5 to 2.0 covers most legitimate sandwich-vs-bootstrap MC noise
  # at this sample size; if either SE is > 2x the other something is
  # structurally wrong with the variance pipeline.
  expect_gt(ratio, 0.5)
  expect_lt(ratio, 2.0)

  # Point estimates should agree across variance methods (they share
  # the fitting path; only the SE construction differs).
  expect_equal(
    res_sand$contrasts$estimate[1],
    res_boot$contrasts$estimate[1],
    tolerance = 1e-6
  )
})


# ----------------------------------------------------------------------
# T-long-ipw4: lmtp_sdr point-estimate cross-check (continuous shift)
# ----------------------------------------------------------------------
test_that("T-long-ipw4: longitudinal IPW shift point estimate agrees with lmtp::lmtp_sdr", {
  skip_if_not_installed("lmtp")
  skip_if_not_installed("SuperLearner")

  # lmtp does its own per-period density chain (TMLE/SDR with the
  # same parametric learners we use). On a correctly specified linear
  # DGP both should converge to the same target parameter. lmtp uses
  # EIF SE so we compare point estimates only.
  set.seed(7004)
  d_long <- make_continuous_scm(n = 1500, seed = 7004)

  # lmtp wants wide format with separate columns per period.
  d_wide <- data.frame(
    L0 = d_long$L0[d_long$time == 0],
    A_0 = d_long$A[d_long$time == 0],
    L_1 = d_long$L[d_long$time == 1],
    A_1 = d_long$A[d_long$time == 1],
    Y = d_long$Y[d_long$time == 1]
  )

  shift_fn <- function(data, trt) data[[trt]] + 0.5

  lmtp_res <- tryCatch(
    suppressWarnings(suppressMessages(lmtp::lmtp_sdr(
      data = d_wide,
      trt = c("A_0", "A_1"),
      outcome = "Y",
      baseline = "L0",
      time_vary = list("L_1"),
      shift = shift_fn,
      outcome_type = "continuous",
      learners_trt = "SL.glm",
      learners_outcome = "SL.glm",
      folds = 1
    ))),
    error = function(e) NULL
  )

  # If lmtp errors at fit time (e.g. SuperLearner unavailable on this
  # platform) skip rather than cascade a failure. The test exists to
  # cross-check an oracle, not to test lmtp itself.
  skip_if(is.null(lmtp_res), "lmtp::lmtp_sdr() unavailable in this environment")

  d_dt <- data.table::as.data.table(d_long)
  fit_ipw <- causat(
    d_dt,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    estimator = "ipw"
  )
  res_ipw <- contrast(
    fit_ipw,
    interventions = list(s = shift(0.5)),
    ci_method = "sandwich"
  )

  est_ipw <- res_ipw$estimates$estimate[1]
  # `lmtp >= 1.5` stores the shifted-treatment estimate in
  # `lmtp_res$estimate@x` (S4 influence_func_estimate object).
  est_lmtp <- tryCatch(lmtp_res$estimate@x, error = function(e) lmtp_res$theta)
  expect_lt(abs(est_ipw - est_lmtp), 0.5)
})


# ----------------------------------------------------------------------
# T-long-ipw5: sequential positivity warning fires on heavy-tail DGP
# ----------------------------------------------------------------------
test_that("T-long-ipw5: warn_seq_positivity fires above threshold and not below", {
  # Direct unit test of the helper. The contrast pipeline path would
  # also surface the upstream `warn_intervened_density_near_zero` from
  # `compute_density_ratio_weights()` first on a heavy-shift DGP, so
  # exercising the helper in isolation keeps the test focused on the
  # longitudinal-specific signal.
  set.seed(7005)

  # Heavy-tail per-period weight matrix: one extreme weight at each
  # period. Threshold 100 catches the period-1 extreme but not the
  # mostly-uniform period-2 column.
  W_bad <- matrix(rnorm(20, 1, 0.1), nrow = 10, ncol = 2)
  W_bad[1, 1] <- 250 # period 1 has a big tail entry

  withr::with_options(list(rlib_message_verbosity = "verbose"), {
    expect_warning(
      warn_seq_positivity(
        W_bad,
        time_points = c("0", "1"),
        threshold = 100
      ),
      class = "causatr_longitudinal_seq_positivity"
    )
  })

  # Below threshold: silent.
  W_ok <- matrix(rnorm(20, 1, 0.5), nrow = 10, ncol = 2)
  expect_no_warning(
    warn_seq_positivity(
      W_ok,
      time_points = c("0", "1"),
      threshold = 100
    ),
    class = "causatr_longitudinal_seq_positivity"
  )
})


# ----------------------------------------------------------------------
# T-long-ipw6: positivity warning does NOT fire on a clean DGP
# ----------------------------------------------------------------------
test_that("T-long-ipw6: positivity warning does NOT fire on a clean DGP under a small shift", {
  set.seed(7006)
  d <- make_continuous_scm(n = 1000, seed = 7006)
  d <- data.table::as.data.table(d)

  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    estimator = "ipw"
  )

  expect_no_warning(
    contrast(
      fit,
      interventions = list(small = shift(0.1), n = NULL),
      ci_method = "sandwich"
    ),
    class = "causatr_longitudinal_seq_positivity"
  )
})


# ----------------------------------------------------------------------
# T-long-ipw7: natural course recovers observed marginal mean
# ----------------------------------------------------------------------
test_that("T-long-ipw7: natural course (NULL intervention) returns observed marginal mean", {
  # Under natural course every per-period weight is identically 1, so
  # the cumulative product is 1 and the Hajek mean is the unweighted
  # mean of Y at the final period.
  set.seed(7007)
  d <- make_linear_scm(n = 800, seed = 7007)
  d <- data.table::as.data.table(d)

  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    estimator = "ipw"
  )
  res <- contrast(
    fit,
    interventions = list(nat = NULL),
    ci_method = "sandwich"
  )

  observed_mean <- mean(d$Y[d$time == 1], na.rm = TRUE)
  expect_equal(
    res$estimates$estimate[res$estimates$intervention == "nat"],
    observed_mean,
    tolerance = 1e-6
  )
  # Sandwich SE for natural course should still be finite
  # (Channel 1 covers the marginal-mean variability even though the
  # propensity correction is zero).
  expect_true(all(is.finite(res$estimates$se)))
})


# ----------------------------------------------------------------------
# Rejection snapshots
# ----------------------------------------------------------------------

test_that("R-long-ipw1: multivariate longitudinal IPW is rejected", {
  set.seed(7100)
  d <- make_linear_scm(n = 100, seed = 7100)
  d <- data.table::as.data.table(d)
  d$A2 <- d$A + stats::rnorm(nrow(d), 0, 0.1)

  expect_error(
    causat(
      d,
      outcome = "Y",
      treatment = c("A", "A2"),
      confounders = ~L0,
      confounders_tv = ~L,
      id = "id",
      time = "time",
      estimator = "ipw"
    ),
    class = "causatr_longitudinal_multivariate_pending"
  )
})


test_that("R-long-ipw2: stabilize = 'marginal' is rejected under longitudinal IPW", {
  set.seed(7101)
  d <- make_linear_scm(n = 100, seed = 7101)
  d <- data.table::as.data.table(d)

  expect_error(
    causat(
      d,
      outcome = "Y",
      treatment = "A",
      confounders = ~L0,
      confounders_tv = ~L,
      id = "id",
      time = "time",
      estimator = "ipw",
      stabilize = "marginal"
    ),
    class = "causatr_longitudinal_stabilize_pending"
  )
})


test_that("R-long-ipw3: numerator = formula is rejected under longitudinal IPW", {
  set.seed(7102)
  d <- make_linear_scm(n = 100, seed = 7102)
  d <- data.table::as.data.table(d)

  expect_error(
    causat(
      d,
      outcome = "Y",
      treatment = "A",
      confounders = ~L0,
      confounders_tv = ~L,
      id = "id",
      time = "time",
      estimator = "ipw",
      numerator = ~1
    ),
    class = "causatr_longitudinal_numerator_pending"
  )
})


test_that("R-long-ipw4: A:modifier (effect modification) is rejected under longitudinal IPW", {
  set.seed(7103)
  d <- make_em_ice_scm(n = 200, n_times = 2, seed = 7103)
  d <- data.table::as.data.table(d)

  expect_error(
    causat(
      d,
      outcome = "Y",
      treatment = "A",
      confounders = ~ L0 + sex + A:sex,
      confounders_tv = ~L,
      id = "id",
      time = "time",
      estimator = "ipw"
    ),
    class = "causatr_longitudinal_em_pending"
  )
})


test_that("R-long-ipw5: ipsi() is rejected under longitudinal IPW", {
  set.seed(7104)
  d <- make_linear_scm(n = 200, seed = 7104)
  d <- data.table::as.data.table(d)

  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    estimator = "ipw"
  )
  expect_error(
    contrast(fit, interventions = list(ip = ipsi(2.0))),
    class = "causatr_longitudinal_ipsi_pending"
  )
})


test_that("R-long-ipw6: ATT is rejected for longitudinal data (any estimator)", {
  # `check_estimand_trt_compat()` has always rejected ATT/ATC under
  # longitudinal data. Re-confirm under longitudinal IPW.
  set.seed(7105)
  d <- make_linear_scm(n = 100, seed = 7105)
  d <- data.table::as.data.table(d)

  expect_error(
    causat(
      d,
      outcome = "Y",
      treatment = "A",
      confounders = ~L0,
      confounders_tv = ~L,
      id = "id",
      time = "time",
      estimator = "ipw",
      estimand = "ATT"
    ),
    "ATT"
  )
})


test_that("R-long-ipw7: single-period data is rejected with a clear error", {
  d <- data.table::data.table(
    id = 1:10,
    time = 0L,
    A = c(0, 1, 0, 1, 0, 1, 0, 1, 0, 1),
    L = stats::rnorm(10),
    Y = stats::rnorm(10)
  )
  expect_error(
    causat(
      d,
      outcome = "Y",
      treatment = "A",
      confounders = ~L,
      id = "id",
      time = "time",
      estimator = "ipw",
      type = "longitudinal"
    ),
    class = "causatr_longitudinal_too_few_times"
  )
})
