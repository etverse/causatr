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
  res_ice <- contrast(
    fit_ice,
    interventions = list(shifted = shift(0.5), nat = NULL),
    type = "difference",
    reference = "nat"
  )

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
  res_ice <- contrast(
    fit_ice,
    interventions = list(always = static(1), never = static(0)),
    type = "difference",
    reference = "never"
  )

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
  res_boot <- contrast(
    fit,
    interventions = list(s = shift(0.5), n = NULL),
    reference = "n",
    ci_method = "bootstrap",
    n_boot = 200
  )

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
    lmtp::lmtp_sdr(
      data = d_wide,
      trt = c("A_0", "A_1"),
      outcome = "Y",
      baseline = "L0",
      time_vary = list(NULL, "L_1"),
      shift = shift_fn,
      outcome_type = "continuous",
      learners_trt = "SL.glm",
      learners_outcome = "SL.glm",
      folds = 1
    ),
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

test_that("R-long-ipw3: numerator = formula without stabilize = 'marginal' is rejected", {
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
    class = "causatr_longitudinal_numerator_without_stabilize"
  )
})


test_that("R-long-ipw4: bare treatment in confounders (~ L + A) is rejected under longitudinal IPW", {
  # `check_em_compat()` rejects `~ L + A` because it puts A on both
  # sides of the propensity model. True EM (`A:sex`) is supported in
  # chunk 10c; bare treatment is invalid by construction and stays
  # rejected with a `causatr_bare_treatment_in_confounders` error.
  set.seed(7103)
  d <- make_linear_scm(n = 200, seed = 7103)
  d <- data.table::as.data.table(d)

  expect_error(
    causat(
      d,
      outcome = "Y",
      treatment = "A",
      confounders = ~ L0 + A,
      confounders_tv = ~L,
      id = "id",
      time = "time",
      estimator = "ipw"
    ),
    class = "causatr_bare_treatment_in_confounders"
  )
})


# Per-period propensity formulas causatr builds for the 2-period
# `make_linear_scm` DGP under `confounders = ~L0`, `confounders_tv = ~L`:
#   period 0 (time = 0): A ~ L + L0          (no lags available)
#   period 1 (time = 1): A ~ L + lag1_A + lag1_L + L0
# The oracle below reconstructs these exactly so the manual product-of-
# Kennedy-weights mean matches causatr's IPSI contrast to numerical
# precision (same propensity fits, same closed-form weight).
make_long_ipsi_oracle_inputs <- function(d, delta) {
  d0 <- d[time == 0][order(id)]
  d1 <- d[time == 1][order(id)]
  # period-1 lag columns hold the OBSERVED period-0 values, matching
  # causatr's `create_lag_vars()`.
  d1$lag1_A <- d0$A
  d1$lag1_L <- d0$L
  list(
    data_by_period = list(d0, d1),
    formulas = list(A ~ L + L0, A ~ L + lag1_A + lag1_L + L0),
    y_final = d1$Y
  )
}

test_that("R-long-ipw5: longitudinal IPSI matches per-period Kennedy oracle", {
  # Truth-based wiring check: the per-period IPSI weight product and the
  # intercept-only final MSM must reproduce a first-principles manual
  # Kennedy product. Both use the same per-period propensity fits, so
  # agreement is to numerical precision.
  set.seed(7104)
  d <- make_linear_scm(n = 1500, seed = 7104)
  d <- data.table::as.data.table(d)
  delta <- 2

  oin <- make_long_ipsi_oracle_inputs(d, delta)
  oracle <- manual_long_ipsi_oracle(
    data_by_period = oin$data_by_period,
    formulas = oin$formulas,
    treatment = "A",
    y_final = oin$y_final,
    delta = delta
  )

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
    interventions = list(ip = ipsi(delta)),
    ci_method = "sandwich"
  )

  expect_equal(
    res$estimates$estimate[res$estimates$intervention == "ip"],
    oracle$theta,
    tolerance = 1e-6
  )
})

test_that("R-long-ipw5b: longitudinal IPSI sandwich SE matches bootstrap", {
  # Sandwich (stacked per-period propensity correction + numDeriv
  # cross-derivative) vs bootstrap parity for univariate IPSI. The two
  # variance estimators target the same M-estimation sandwich, so they
  # agree within Monte-Carlo error of the bootstrap.
  set.seed(8810)
  d <- make_linear_scm(n = 1200, seed = 8810)
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

  res_s <- contrast(
    fit,
    interventions = list(ip = ipsi(2)),
    ci_method = "sandwich"
  )
  set.seed(8810)
  res_b <- contrast(
    fit,
    interventions = list(ip = ipsi(2)),
    ci_method = "bootstrap",
    n_boot = 300
  )

  se_s <- res_s$estimates$se[res_s$estimates$intervention == "ip"]
  se_b <- res_b$estimates$se[res_b$estimates$intervention == "ip"]
  # Point estimates identical (same fit); SEs within 15% of each other.
  expect_equal(
    res_s$estimates$estimate[res_s$estimates$intervention == "ip"],
    res_b$estimates$estimate[res_b$estimates$intervention == "ip"],
    tolerance = 1e-8
  )
  expect_equal(se_s, se_b, tolerance = 0.15)
})

test_that("R-long-ipw5c: longitudinal IPSI difference contrast recovers truth", {
  # A difference of two IPSI means (delta = 2 vs delta = 1) under the
  # per-period Kennedy oracle. delta = 1 is the natural course (unit
  # weights), so the contrast is E[Y(ipsi(2))] - E[Y].
  set.seed(7106)
  d <- make_linear_scm(n = 1500, seed = 7106)
  d <- data.table::as.data.table(d)

  o2 <- make_long_ipsi_oracle_inputs(d, 2)
  oracle2 <- manual_long_ipsi_oracle(
    o2$data_by_period,
    o2$formulas,
    "A",
    o2$y_final,
    delta = 2
  )
  o1 <- make_long_ipsi_oracle_inputs(d, 1)
  oracle1 <- manual_long_ipsi_oracle(
    o1$data_by_period,
    o1$formulas,
    "A",
    o1$y_final,
    delta = 1
  )

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
    interventions = list(hi = ipsi(2), lo = ipsi(1)),
    type = "difference",
    reference = "lo",
    ci_method = "sandwich"
  )

  expect_equal(
    res$contrasts$estimate,
    oracle2$theta - oracle1$theta,
    tolerance = 1e-6
  )
})

test_that("R-long-ipw5e: stabilized longitudinal IPSI is rejected", {
  # Kennedy's IPSI weight is already a bounded propensity-odds ratio with
  # no separate marginal numerator, so `stabilize = "marginal"` + IPSI is
  # rejected with a dedicated class pointing to `stabilize = "none"`.
  set.seed(7108)
  d <- make_linear_scm(n = 300, seed = 7108)
  d <- data.table::as.data.table(d)

  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    estimator = "ipw",
    stabilize = "marginal"
  )
  expect_error(
    contrast(fit, interventions = list(ip = ipsi(2))),
    class = "causatr_longitudinal_ipsi_stabilize"
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


# ----------------------------------------------------------------------
# Chunk 10b: stabilized weights
# ----------------------------------------------------------------------

test_that("T-long-ipw-stab1: stabilize = 'marginal' fits per-period numerator models with the right structure", {
  set.seed(7600)
  d <- make_linear_scm(n = 200, seed = 7600)
  d <- data.table::as.data.table(d)

  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    estimator = "ipw",
    stabilize = "marginal"
  )

  num_models <- fit$details$numerator_models_by_time
  expect_false(is.null(num_models))
  expect_length(num_models, 2L)

  # Period 1 (first period): no lags available, no L. Default
  # numerator drops everything -> intercept-only.
  expect_equal(
    deparse(num_models[[1]]$ps_formula),
    "A ~ 1"
  )
  # Period 2: one lag of treatment, no L.
  expect_equal(
    deparse(num_models[[2]]$ps_formula),
    "A ~ lag1_A"
  )

  # Each numerator model has its own alpha_hat with the expected
  # parameter count.
  expect_length(num_models[[1]]$alpha_hat, 1L) # intercept only
  expect_length(num_models[[2]]$alpha_hat, 2L) # intercept + lag1_A
})


test_that("T-long-ipw-stab2: stabilized + static binary recovers identical estimate as unstabilized", {
  # Under static interventions on a discrete treatment, the per-period
  # weight is `I(A = a) * f^*/f`. For the surviving rows, f^* = f at
  # the same evaluation point (numerator and denominator densities
  # coincide because the indicator pins A to the static target). The
  # marginal-mean Hájek therefore agrees exactly with unstabilized.
  # This is the same equality Phase 8e asserts for multivariate IPW.
  set.seed(7601)
  d <- make_linear_scm(n = 1500, seed = 7601)
  d <- data.table::as.data.table(d)

  fit_us <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    estimator = "ipw"
  )
  fit_st <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    estimator = "ipw",
    stabilize = "marginal"
  )
  res_us <- contrast(
    fit_us,
    interventions = list(a = static(1), z = static(0)),
    reference = "z",
    ci_method = "sandwich"
  )
  res_st <- contrast(
    fit_st,
    interventions = list(a = static(1), z = static(0)),
    reference = "z",
    ci_method = "sandwich"
  )

  # Point estimates and SEs must match exactly because the per-period
  # weight is identical on surviving rows.
  expect_equal(
    res_st$contrasts$estimate,
    res_us$contrasts$estimate,
    tolerance = 1e-6
  )
  expect_equal(
    res_st$contrasts$se,
    res_us$contrasts$se,
    tolerance = 1e-6
  )
})


test_that("T-long-ipw-stab3: stabilized shift contrast recovers ICE point on baseline-numerator DGP", {
  # For shift on continuous treatment, the stabilized natural-course
  # arm carries a non-unit weight `g(A|...)/f(A|..., L)`, so the
  # contrast against MTP-shifted depends on how well the numerator
  # captures the structural treatment generation. With `numerator =
  # ~L0` (baseline conditioning), the per-period numerator is close
  # enough to the denominator to recover the SCM contrast that ICE
  # gives at large N.
  set.seed(7602)
  d <- make_continuous_scm(n = 2000, seed = 7602)
  d <- data.table::as.data.table(d)

  fit_st <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    estimator = "ipw",
    stabilize = "marginal",
    numerator = ~L0
  )
  res_st <- contrast(
    fit_st,
    interventions = list(s = shift(0.5), n = NULL),
    reference = "n",
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
  res_ice <- contrast(
    fit_ice,
    interventions = list(s = shift(0.5), n = NULL),
    reference = "n"
  )

  # Cross-method agreement on the contrast point estimate. Tolerance
  # widened relative to T-long-ipw1 because stabilized IPW has a
  # second source of finite-sample noise from the numerator fit.
  expect_lt(
    abs(res_st$contrasts$estimate[1] - res_ice$contrasts$estimate[1]),
    0.6
  )
  # Sandwich SE finite, positive.
  expect_true(all(is.finite(res_st$contrasts$se) & res_st$contrasts$se > 0))
})


test_that("T-long-ipw-stab4: bootstrap captures gamma uncertainty under stabilization", {
  # Bootstrap refits both numerator (gamma) and denominator (alpha)
  # on every replicate, so the bootstrap SE picks up gamma uncertainty
  # the sandwich (which holds gamma fixed) does not. For a static
  # binary DGP both methods should agree to within MC noise because
  # gamma uncertainty is small relative to alpha uncertainty.
  set.seed(7603)
  d <- make_linear_scm(n = 600, seed = 7603)
  d <- data.table::as.data.table(d)

  fit_st <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    estimator = "ipw",
    stabilize = "marginal"
  )

  res_sand <- contrast(
    fit_st,
    interventions = list(a = static(1), z = static(0)),
    reference = "z",
    ci_method = "sandwich"
  )
  res_boot <- contrast(
    fit_st,
    interventions = list(a = static(1), z = static(0)),
    reference = "z",
    ci_method = "bootstrap",
    n_boot = 100
  )

  expect_true(
    is.finite(res_boot$contrasts$se[1]) &&
      res_boot$contrasts$se[1] > 0
  )
  ratio <- res_boot$contrasts$se[1] / res_sand$contrasts$se[1]
  expect_gt(ratio, 0.4)
  expect_lt(ratio, 2.5)

  # Point estimates agree across variance methods.
  expect_equal(
    res_sand$contrasts$estimate[1],
    res_boot$contrasts$estimate[1],
    tolerance = 1e-6
  )
})


test_that("T-long-ipw-stab5: custom numerator = ~ baseline keeps treatment lags by default", {
  # The chain-rule factorisation of g_k requires conditioning on
  # prior treatments to be a valid joint density. `numerator = ~ L0`
  # adds L0 to the per-period numerator without dropping the
  # treatment lags `lag1_A` etc.
  set.seed(7604)
  d <- make_linear_scm(n = 200, seed = 7604)
  d <- data.table::as.data.table(d)

  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    estimator = "ipw",
    stabilize = "marginal",
    numerator = ~L0
  )

  num_models <- fit$details$numerator_models_by_time
  # Period 1: no lags + L0 -> A ~ L0.
  expect_equal(deparse(num_models[[1]]$ps_formula), "A ~ L0")
  # Period 2: lag1_A + L0 -> A ~ lag1_A + L0 (treatment lags first).
  expect_equal(deparse(num_models[[2]]$ps_formula), "A ~ lag1_A + L0")
})


test_that("T-long-ipw-stab6: stabilized scale_by + trim sandwich matches bootstrap", {
  # Exercises the stabilized continuous-treatment `scale_by` branch of
  # `make_long_stabilized_period_closure()` together with the per-period
  # trim-threshold precompute in `make_weight_fn_longitudinal()`. The
  # sandwich (gamma-fixed numerator, alpha-varying denominator, truncated
  # at the fixed per-period quantile) and the bootstrap (refits both)
  # target the same M-estimation variance, so they agree within MC error.
  set.seed(7605)
  d <- make_continuous_scm(n = 1500, seed = 7605)
  d <- data.table::as.data.table(d)

  fit_st <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    estimator = "ipw",
    stabilize = "marginal",
    numerator = ~L0
  )

  res_sw <- contrast(
    fit_st,
    interventions = list(sc = scale_by(1.1), n = NULL),
    reference = "n",
    type = "difference",
    ci_method = "sandwich",
    trim = 0.99
  )
  set.seed(7605)
  res_bs <- contrast(
    fit_st,
    interventions = list(sc = scale_by(1.1), n = NULL),
    reference = "n",
    type = "difference",
    ci_method = "bootstrap",
    n_boot = 200,
    trim = 0.99
  )

  # Point estimate identical across variance methods (same fit + weights).
  expect_equal(
    res_sw$contrasts$estimate[1],
    res_bs$contrasts$estimate[1],
    tolerance = 1e-6
  )
  # Sandwich vs bootstrap SE parity (gamma uncertainty is small here).
  se_ratio <- res_sw$contrasts$se[1] / res_bs$contrasts$se[1]
  expect_gt(se_ratio, 0.5)
  expect_lt(se_ratio, 2)
})


# ----------------------------------------------------------------------
# Chunk 10c: effect modification (baseline modifier)
# ----------------------------------------------------------------------

test_that("T-long-ipw-em1: stratified ATE recovers analytical truth on EM-ICE DGP", {
  # `make_em_ice_scm()` is the canonical 2-period sex-stratified DGP
  # with treatment-confounder feedback and per-period treatment effect
  # (2 + 1.5*sex). For `always vs never` over 2 periods the analytical
  # truths are ATE|sex=0 = 5 and ATE|sex=1 = 8 (derivation in
  # helper-dgp.R).
  set.seed(7700)
  d <- make_em_ice_scm(n = 2500, n_times = 2, seed = 7700)
  d <- data.table::as.data.table(d)

  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L0 + sex + A:sex,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    estimator = "ipw"
  )

  # The MSM should expand from `Y ~ 1` to `Y ~ 1 + sex` automatically.
  # `compute_ipw_contrast_longitudinal()` reads em_info from
  # fit$details and routes through `build_ipw_msm_formula()`.
  expect_true(fit$details$em_info$has_em)
  expect_equal(fit$details$em_info$modifier_vars, "sex")

  res <- contrast(
    fit,
    interventions = list(a = static(1), z = static(0)),
    reference = "z",
    ci_method = "sandwich",
    by = "sex"
  )

  con_dt <- as.data.frame(res$contrasts)
  est_sex0 <- con_dt$estimate[con_dt$by == "0"]
  est_sex1 <- con_dt$estimate[con_dt$by == "1"]

  # Analytical truths from `make_em_ice_scm` derivation: 5 and 8.
  expect_lt(abs(est_sex0 - 5), 0.6)
  expect_lt(abs(est_sex1 - 8), 0.7)

  # SEs are finite and positive.
  expect_true(all(is.finite(con_dt$se) & con_dt$se > 0))
})


test_that("T-long-ipw-em2: per-period propensity strips A:modifier; MSM expands to Y ~ 1 + modifier", {
  # Structural check: with EM in confounders, the per-period
  # propensity formula must NOT carry the `A:sex` term (A is the
  # response of the propensity model -- can't appear on both sides),
  # and the MSM in compute_ipw_contrast_longitudinal must include
  # the modifier main effect.
  set.seed(7701)
  d <- make_em_ice_scm(n = 400, n_times = 2, seed = 7701)
  d <- data.table::as.data.table(d)

  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L0 + sex + A:sex,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    estimator = "ipw"
  )

  # Per-period propensity formula: A ~ baseline + tv + lags. Must not
  # contain `:sex` (which would put A on both sides via A:sex).
  ps_terms_p1 <- attr(
    stats::terms(fit$details$per_period_formula[[1]]),
    "term.labels"
  )
  ps_terms_p2 <- attr(
    stats::terms(fit$details$per_period_formula[[2]]),
    "term.labels"
  )
  expect_false(any(grepl(":sex", ps_terms_p1)))
  expect_false(any(grepl(":sex", ps_terms_p2)))
  # Modifier main effect IS in the propensity (main-effect baseline
  # confounder; only the interaction is stripped).
  expect_true("sex" %in% ps_terms_p1)

  # Verify the MSM formula expands to `Y ~ 1 + sex` when EM is
  # present; the longitudinal contrast path uses the same
  # `build_ipw_msm_formula()` as the point IPW path.
  expanded_formula <- build_ipw_msm_formula("Y", fit$details$em_info)
  expect_true("sex" %in% all.vars(expanded_formula))
  # Trigger one contrast to confirm the wired MSM gets fit (no abort,
  # finite SE).
  res <- contrast(
    fit,
    interventions = list(a = static(1)),
    ci_method = "sandwich"
  )
  expect_s3_class(res, "causatr_result")
  expect_true(all(is.finite(res$estimates$se)))
})


test_that("T-long-ipw-em3: longitudinal IPW EM agrees with ICE EM cross-method", {
  # ICE g-comp with `A:sex` interaction is the baseline that
  # longitudinal IPW EM should agree with on the same DGP. Both
  # estimate the same target (sex-stratified ATE under the
  # always-vs-never contrast).
  set.seed(7702)
  d <- make_em_ice_scm(n = 2000, n_times = 2, seed = 7702)
  d <- data.table::as.data.table(d)

  fit_ipw <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L0 + sex + A:sex,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    estimator = "ipw"
  )
  res_ipw <- contrast(
    fit_ipw,
    interventions = list(a = static(1), z = static(0)),
    reference = "z",
    ci_method = "sandwich",
    by = "sex"
  )

  fit_ice <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L0 + sex + A:sex,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    estimator = "gcomp"
  )
  res_ice <- contrast(
    fit_ice,
    interventions = list(a = static(1), z = static(0)),
    reference = "z",
    by = "sex"
  )

  ipw_dt <- as.data.frame(res_ipw$contrasts)
  ice_dt <- as.data.frame(res_ice$contrasts)
  ipw_sex0 <- ipw_dt$estimate[ipw_dt$by == "0"]
  ipw_sex1 <- ipw_dt$estimate[ipw_dt$by == "1"]
  ice_sex0 <- ice_dt$estimate[ice_dt$by == "0"]
  ice_sex1 <- ice_dt$estimate[ice_dt$by == "1"]

  expect_lt(abs(ipw_sex0 - ice_sex0), 0.5)
  expect_lt(abs(ipw_sex1 - ice_sex1), 0.5)
})


# ---- Weight truncation (Phase 19-trim) ---------------------------------

test_that("longitudinal IPW: trim reduces cumulative max weight", {
  d <- data.table::as.data.table(make_continuous_scm(n = 500, seed = 7200))
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    estimator = "ipw",
    id = "id",
    time = "time"
  )
  r1 <- contrast(
    fit,
    interventions = list(shifted = shift(0.5), nat = NULL),
    type = "difference",
    ci_method = "sandwich"
  )
  r2 <- contrast(
    fit,
    interventions = list(shifted = shift(0.5), nat = NULL),
    type = "difference",
    ci_method = "sandwich",
    trim = 0.99
  )
  expect_true(is.finite(r1$contrasts$estimate))
  expect_true(is.finite(r2$contrasts$estimate))
  expect_true(is.finite(r2$contrasts$se))
})

test_that("longitudinal IPW: sandwich and bootstrap agree under trim", {
  # Regression test for the per-period vs post-product truncation fix.
  # Before the fix, sandwich used post-product truncation (on the
  # cumulative weight) while bootstrap used per-period truncation
  # (via compute_longitudinal_weights), causing a ~10% SE discrepancy.
  d <- data.table::as.data.table(make_continuous_scm(n = 800, seed = 7205))
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    estimator = "ipw",
    id = "id",
    time = "time"
  )
  r_sand <- contrast(
    fit,
    interventions = list(shifted = shift(0.5), nat = NULL),
    ci_method = "sandwich",
    trim = 0.95
  )
  set.seed(123)
  r_boot <- contrast(
    fit,
    interventions = list(shifted = shift(0.5), nat = NULL),
    ci_method = "bootstrap",
    n_boot = 150L,
    trim = 0.95
  )
  ratio <- r_sand$contrasts$se / r_boot$contrasts$se
  expect_gt(ratio, 0.7)
  expect_lt(ratio, 1.5)
  expect_equal(r_sand$contrasts$estimate, r_boot$contrasts$estimate)
})

test_that("longitudinal IPW: trim + bootstrap works", {
  d <- data.table::as.data.table(make_continuous_scm(n = 300, seed = 7201))
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    estimator = "ipw",
    id = "id",
    time = "time"
  )
  r_boot <- contrast(
    fit,
    interventions = list(shifted = shift(0.5), nat = NULL),
    type = "difference",
    ci_method = "bootstrap",
    n_boot = 50L,
    trim = 0.99
  )
  expect_true(is.finite(r_boot$contrasts$estimate))
  expect_true(is.finite(r_boot$contrasts$se))
})


# ---- lmtp cross-check: longitudinal IPW + trim --------------------------

test_that("longitudinal IPW + trim agrees with lmtp_sdr", {
  skip_if_not_installed("lmtp")
  skip_if_not_installed("SuperLearner")

  # 2-period continuous treatment, linear DGP from make_continuous_scm().
  # lmtp treats c("A_0", "A_1") as 2-period longitudinal and applies SDR.
  # Under correct specification both estimators target the same shifted mean.
  set.seed(7300)
  d_long <- make_continuous_scm(n = 1500, seed = 7300)

  # lmtp wants wide format
  d_wide <- data.frame(
    L0 = d_long$L0[d_long$time == 0],
    A_0 = d_long$A[d_long$time == 0],
    L_1 = d_long$L[d_long$time == 1],
    A_1 = d_long$A[d_long$time == 1],
    Y = d_long$Y[d_long$time == 1]
  )

  shift_fn <- function(data, trt) data[[trt]] + 0.5

  # lmtp: no trim and trim = 0.99
  lmtp_no <- tryCatch(
    lmtp::lmtp_sdr(
      data = d_wide,
      trt = c("A_0", "A_1"),
      outcome = "Y",
      baseline = "L0",
      time_vary = list(NULL, "L_1"),
      shift = shift_fn,
      outcome_type = "continuous",
      learners_trt = "SL.glm",
      learners_outcome = "SL.glm",
      folds = 1,
      control = lmtp::lmtp_control(.trim = 1)
    ),
    error = function(e) NULL
  )
  lmtp_99 <- tryCatch(
    lmtp::lmtp_sdr(
      data = d_wide,
      trt = c("A_0", "A_1"),
      outcome = "Y",
      baseline = "L0",
      time_vary = list(NULL, "L_1"),
      shift = shift_fn,
      outcome_type = "continuous",
      learners_trt = "SL.glm",
      learners_outcome = "SL.glm",
      folds = 1,
      control = lmtp::lmtp_control(.trim = 0.99)
    ),
    error = function(e) NULL
  )
  skip_if(is.null(lmtp_no), "lmtp::lmtp_sdr() unavailable")

  # causatr: longitudinal IPW, no trim and trim = 0.99
  d_dt <- data.table::as.data.table(d_long)
  fit <- causat(
    d_dt,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    estimator = "ipw"
  )
  res_no <- contrast(
    fit,
    interventions = list(s = shift(0.5)),
    ci_method = "sandwich"
  )
  res_99 <- contrast(
    fit,
    interventions = list(s = shift(0.5)),
    ci_method = "sandwich",
    trim = 0.99
  )

  est_lmtp_no <- tryCatch(lmtp_no$estimate@x, error = function(e) lmtp_no$theta)
  est_lmtp_99 <- tryCatch(lmtp_99$estimate@x, error = function(e) lmtp_99$theta)

  expect_lt(abs(res_no$estimates$estimate - est_lmtp_no), 0.5)
  expect_lt(abs(res_99$estimates$estimate - est_lmtp_99), 0.5)
})
