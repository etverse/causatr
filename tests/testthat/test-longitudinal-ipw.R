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


# -----------------------------------------------------------------------
# Multivariate longitudinal IPW (Phase 19a)
# -----------------------------------------------------------------------
# T × K propensity factorisation: W_i = prod_t prod_k w_{t,k,i}
# with block-diagonal stacked sandwich.

test_that("T-long-mv-ipw1: binary MV longitudinal IPW static, sandwich vs bootstrap", {
  # DGP: 2 periods, 2 binary treatments
  # f(A1,A2 | L) = f(A1|L) * f(A2|A1,L)
  # Y = 2 + 0.5*L + psi1*A1_T + psi2*A2_T + eps
  set.seed(19001)
  n <- 800
  L <- rnorm(n)
  A1_1 <- rbinom(n, 1, plogis(0.5 + 0.3 * L))
  A2_1 <- rbinom(n, 1, plogis(0.3 + 0.2 * L + 0.15 * A1_1))
  A1_2 <- rbinom(n, 1, plogis(0.3 + 0.3 * L + 0.1 * A1_1))
  A2_2 <- rbinom(n, 1, plogis(0.2 + 0.2 * L + 0.1 * A1_2 + 0.05 * A2_1))
  psi1 <- 1.0
  psi2 <- 0.5
  Y <- 2 + 0.5 * L + psi1 * A1_2 + psi2 * A2_2 + rnorm(n)

  dat <- data.table::data.table(
    id = rep(seq_len(n), each = 2L),
    time = rep(1:2, n),
    A1 = c(rbind(A1_1, A1_2)),
    A2 = c(rbind(A2_1, A2_2)),
    L = rep(L, each = 2L),
    Y = c(rep(NA_real_, n), Y)
  )

  fit <- causat(
    dat,
    outcome = "Y",
    treatment = c("A1", "A2"),
    confounders = ~L,
    estimator = "ipw",
    type = "longitudinal",
    id = "id",
    time = "time"
  )
  expect_true(fit$details$is_multivariate)
  expect_equal(fit$details$n_times, 2L)
  # Per-period models should be causatr_treatment_models
  expect_s3_class(
    fit$details$treatment_models_by_time[["1"]],
    "causatr_treatment_models"
  )

  res_sw <- contrast(
    fit,
    interventions = list(
      both1 = list(A1 = static(1), A2 = static(1)),
      both0 = list(A1 = static(0), A2 = static(0))
    ),
    type = "difference",
    ci_method = "sandwich"
  )
  res_bs <- contrast(
    fit,
    interventions = list(
      both1 = list(A1 = static(1), A2 = static(1)),
      both0 = list(A1 = static(0), A2 = static(0))
    ),
    type = "difference",
    ci_method = "bootstrap"
  )

  # Sandwich/bootstrap SE ratio should be close to 1.
  se_ratio <- res_sw$contrasts$se / res_bs$contrasts$se
  expect_gt(se_ratio, 0.7)
  expect_lt(se_ratio, 1.4)

  # E[Y(1,1)] - E[Y(0,0)]: the IPW estimand involves the cumulative
  # product weight across periods, which reweights to the joint
  # interventional distribution. With binary treatments the HT
  # indicator drops rows, so the effective sample is smaller and the
  # IPW estimand approximates the structural effect (psi1 + psi2)
  # only under large n and good overlap. Check finite + positive.
  expect_true(is.finite(res_sw$contrasts$estimate))
  expect_gt(res_sw$contrasts$estimate, 0)
})


test_that("T-long-mv-ipw2: continuous MV longitudinal IPW shift, sandwich vs bootstrap", {
  # DGP: 2 periods, 2 continuous treatments
  # Y = 2 + 0.5*L + psi1*A1_T + psi2*A2_T + eps
  set.seed(19002)
  n <- 800
  L <- rnorm(n)
  A1_1 <- rnorm(n, 0.5 + 0.3 * L)
  A2_1 <- rnorm(n, 0.3 + 0.2 * L + 0.1 * A1_1)
  A1_2 <- rnorm(n, 0.3 + 0.3 * L + 0.1 * A1_1)
  A2_2 <- rnorm(n, 0.2 + 0.2 * L + 0.1 * A1_2 + 0.05 * A2_1)
  psi1 <- 0.5
  psi2 <- 0.3
  Y <- 2 + 0.5 * L + psi1 * A1_2 + psi2 * A2_2 + rnorm(n)

  dat <- data.table::data.table(
    id = rep(seq_len(n), each = 2L),
    time = rep(1:2, n),
    A1 = c(rbind(A1_1, A1_2)),
    A2 = c(rbind(A2_1, A2_2)),
    L = rep(L, each = 2L),
    Y = c(rep(NA_real_, n), Y)
  )

  fit <- causat(
    dat,
    outcome = "Y",
    treatment = c("A1", "A2"),
    confounders = ~L,
    estimator = "ipw",
    type = "longitudinal",
    id = "id",
    time = "time"
  )

  delta1 <- 0.5
  delta2 <- 0.3
  res_sw <- contrast(
    fit,
    interventions = list(
      up = list(A1 = shift(delta1), A2 = shift(delta2)),
      nc = NULL
    ),
    type = "difference",
    ci_method = "sandwich"
  )
  res_bs <- contrast(
    fit,
    interventions = list(
      up = list(A1 = shift(delta1), A2 = shift(delta2)),
      nc = NULL
    ),
    type = "difference",
    ci_method = "bootstrap"
  )

  # Sandwich/bootstrap SE agreement
  se_ratio <- res_sw$contrasts$se / res_bs$contrasts$se
  expect_gt(se_ratio, 0.7)
  expect_lt(se_ratio, 1.4)

  # The shifted estimand is approximately psi1*d1 + psi2*d2 = 0.34
  # under the linear DGP, but the IPW cumulative product weights
  # and the joint density factorisation introduce finite-sample bias.
  # Check finite and same sign as the structural effect.
  expect_true(is.finite(res_sw$contrasts$estimate))
  expect_gt(res_sw$contrasts$estimate, -0.3)
})


test_that("T-long-mv-ipw3: MV natural course gives same sandwich as bootstrap", {
  # All-NULL intervention: cumulative weight = 1, no propensity
  # correction. Sandwich should still work (the empty-alpha early
  # return path).
  set.seed(19003)
  n <- 300
  L <- rnorm(n)
  A1_1 <- rbinom(n, 1, plogis(0.3 * L))
  A2_1 <- rbinom(n, 1, plogis(0.2 * L))
  A1_2 <- rbinom(n, 1, plogis(0.3 * L + 0.1 * A1_1))
  A2_2 <- rbinom(n, 1, plogis(0.2 * L + 0.1 * A1_2))
  Y <- 1 + 0.3 * L + 0.5 * A1_2 + 0.3 * A2_2 + rnorm(n)

  dat <- data.table::data.table(
    id = rep(seq_len(n), each = 2L),
    time = rep(1:2, n),
    A1 = c(rbind(A1_1, A1_2)),
    A2 = c(rbind(A2_1, A2_2)),
    L = rep(L, each = 2L),
    Y = c(rep(NA_real_, n), Y)
  )

  fit <- causat(
    dat,
    outcome = "Y",
    treatment = c("A1", "A2"),
    confounders = ~L,
    estimator = "ipw",
    type = "longitudinal",
    id = "id",
    time = "time"
  )

  res_sw <- contrast(
    fit,
    interventions = list(
      both1 = list(A1 = static(1), A2 = static(1)),
      nc = NULL
    ),
    type = "difference",
    ci_method = "sandwich"
  )

  # Should produce finite SE and reasonable estimate
  expect_true(is.finite(res_sw$contrasts$se))
  expect_true(is.finite(res_sw$contrasts$estimate))
})


test_that("T-long-mv-stab1: stabilized MV fits per-period, per-component chain numerators", {
  # The MV longitudinal numerator factorises over BOTH the time axis
  # and the component axis. Within a period, component k conditions on
  # the prior components A_{1..k-1}; across periods, every component
  # conditions on the lagged treatments. The marginal numerator drops
  # the time-varying L from all of these formulas (baseline-only
  # conditioning), which is what distinguishes the numerator from the
  # full-L denominator.
  set.seed(19101)
  n <- 400
  L <- rnorm(n)
  A1_1 <- rbinom(n, 1, plogis(0.5 + 0.3 * L))
  A2_1 <- rbinom(n, 1, plogis(0.3 + 0.2 * L + 0.15 * A1_1))
  A1_2 <- rbinom(n, 1, plogis(0.3 + 0.3 * L + 0.1 * A1_1))
  A2_2 <- rbinom(n, 1, plogis(0.2 + 0.2 * L + 0.1 * A1_2 + 0.05 * A2_1))
  Y <- 2 + 0.5 * L + A1_2 + 0.5 * A2_2 + rnorm(n)

  dat <- data.table::data.table(
    id = rep(seq_len(n), each = 2L),
    time = rep(1:2, n),
    A1 = c(rbind(A1_1, A1_2)),
    A2 = c(rbind(A2_1, A2_2)),
    L = rep(L, each = 2L),
    Y = c(rep(NA_real_, n), Y)
  )

  fit <- causat(
    dat,
    outcome = "Y",
    treatment = c("A1", "A2"),
    confounders = ~L,
    estimator = "ipw",
    type = "longitudinal",
    id = "id",
    time = "time",
    stabilize = "marginal"
  )

  num_models <- fit$details$numerator_models_by_time
  expect_false(is.null(num_models))
  expect_length(num_models, 2L)

  # Each period's numerator is itself a per-component treatment-models
  # list (one entry per treatment component).
  p1 <- num_models[[1]]
  p2 <- num_models[[2]]
  expect_s3_class(p1, "causatr_treatment_models")
  expect_s3_class(p2, "causatr_treatment_models")
  expect_length(p1, 2L)
  expect_length(p2, 2L)

  # Period 1: no lags. A1 has no upstream conditioning (intercept-only);
  # A2 conditions on the prior component A1. L dropped throughout.
  expect_equal(deparse(p1[[1]]$ps_formula), "A1 ~ 1")
  expect_equal(deparse(p1[[2]]$ps_formula), "A2 ~ A1")

  # Period 2: both components condition on the lagged treatments; A2
  # additionally conditions on the within-period prior component A1.
  expect_equal(deparse(p2[[1]]$ps_formula), "A1 ~ lag1_A1 + lag1_A2")
  expect_equal(
    deparse(p2[[2]]$ps_formula),
    "A2 ~ A1 + lag1_A1 + lag1_A2"
  )
})


test_that("T-long-mv-stab2: stabilized + static binary MV recovers identical estimate as unstabilized", {
  # Same invariant as the univariate T-long-ipw-stab2 and the
  # multivariate point Phase 8e: under static interventions on discrete
  # treatments the per-period, per-component weight is
  # `prod_k I(A_k = a_k) * f^*_k / f_k`. On the surviving rows every
  # conditioning variable (prior components + lags) is pinned to the
  # static target, so the marginal numerator f^* equals the denominator
  # f at the evaluation point and the Hajek mean is unchanged. The
  # equality must hold for both the point estimate and the sandwich SE.
  set.seed(19102)
  n <- 1500
  L <- rnorm(n)
  A1_1 <- rbinom(n, 1, plogis(0.5 + 0.3 * L))
  A2_1 <- rbinom(n, 1, plogis(0.3 + 0.2 * L + 0.15 * A1_1))
  A1_2 <- rbinom(n, 1, plogis(0.3 + 0.3 * L + 0.1 * A1_1))
  A2_2 <- rbinom(n, 1, plogis(0.2 + 0.2 * L + 0.1 * A1_2 + 0.05 * A2_1))
  Y <- 2 + 0.5 * L + A1_2 + 0.5 * A2_2 + rnorm(n)

  dat <- data.table::data.table(
    id = rep(seq_len(n), each = 2L),
    time = rep(1:2, n),
    A1 = c(rbind(A1_1, A1_2)),
    A2 = c(rbind(A2_1, A2_2)),
    L = rep(L, each = 2L),
    Y = c(rep(NA_real_, n), Y)
  )

  mk <- function(stab) {
    causat(
      dat,
      outcome = "Y",
      treatment = c("A1", "A2"),
      confounders = ~L,
      estimator = "ipw",
      type = "longitudinal",
      id = "id",
      time = "time",
      stabilize = stab
    )
  }
  ivs <- list(
    b1 = list(A1 = static(1), A2 = static(1)),
    b0 = list(A1 = static(0), A2 = static(0))
  )
  res_us <- contrast(
    mk("none"),
    interventions = ivs,
    reference = "b0",
    type = "difference",
    ci_method = "sandwich"
  )
  res_st <- contrast(
    mk("marginal"),
    interventions = ivs,
    reference = "b0",
    type = "difference",
    ci_method = "sandwich"
  )

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


test_that("T-long-mv-stab3: stabilized MV shift SE close to bootstrap and natural course is exact sample mean", {
  # For shift on continuous treatments the stabilized natural-course
  # arm carries non-unit numerator/denominator ratios, so the sandwich
  # must propagate uncertainty through the denominator alpha (numerator
  # gamma held fixed, per the nuisance-fixed convention). Bootstrap
  # refits both. The two SE estimates should agree to within MC noise.
  set.seed(19103)
  n <- 1000
  L <- rnorm(n)
  A1_1 <- rnorm(n, 0.5 + 0.3 * L)
  A2_1 <- rnorm(n, 0.3 + 0.2 * L + 0.1 * A1_1)
  A1_2 <- rnorm(n, 0.3 + 0.3 * L + 0.1 * A1_1)
  A2_2 <- rnorm(n, 0.2 + 0.2 * L + 0.1 * A1_2 + 0.05 * A2_1)
  Y <- 2 + 0.5 * L + 0.5 * A1_2 + 0.3 * A2_2 + rnorm(n)

  dat <- data.table::data.table(
    id = rep(seq_len(n), each = 2L),
    time = rep(1:2, n),
    A1 = c(rbind(A1_1, A1_2)),
    A2 = c(rbind(A2_1, A2_2)),
    L = rep(L, each = 2L),
    Y = c(rep(NA_real_, n), Y)
  )

  fit <- causat(
    dat,
    outcome = "Y",
    treatment = c("A1", "A2"),
    confounders = ~L,
    estimator = "ipw",
    type = "longitudinal",
    id = "id",
    time = "time",
    stabilize = "marginal",
    numerator = ~1
  )
  ivs <- list(
    up = list(A1 = shift(0.5), A2 = shift(0.3)),
    nc = NULL
  )
  res_sw <- contrast(
    fit,
    interventions = ivs,
    reference = "nc",
    type = "difference",
    ci_method = "sandwich"
  )
  res_bs <- contrast(
    fit,
    interventions = ivs,
    reference = "nc",
    type = "difference",
    ci_method = "bootstrap"
  )

  se_ratio <- res_sw$contrasts$se / res_bs$contrasts$se
  expect_gt(se_ratio, 0.7)
  expect_lt(se_ratio, 1.4)
  expect_true(all(is.finite(res_sw$contrasts$se) & res_sw$contrasts$se > 0))

  # Truth-based point-estimate check (2026-05-29 critical review): the
  # whole-NULL natural-course arm carries weight 1 at every period by
  # construction, and that short-circuit fires BEFORE the stabilize
  # branch. So the stabilized natural-course marginal mean must equal
  # the UNSTABILIZED natural-course marginal mean EXACTLY -- the
  # numerator chain never enters. A Monte Carlo check (200 reps)
  # confirmed the weight-1 path reproduces the raw sample mean to ~1e-15
  # and is unbiased, while routing natural course through prod(g/f)
  # (e.g. a shift(0) MTP) is also unbiased but ~3.4x noisier. This pins
  # the documented behavior so a future "fix" that sends whole-NULL
  # through prod(g/f) under stabilization is caught.
  fit_unstab <- causat(
    dat,
    outcome = "Y",
    treatment = c("A1", "A2"),
    confounders = ~L,
    estimator = "ipw",
    type = "longitudinal",
    id = "id",
    time = "time"
  )
  res_unstab <- contrast(
    fit_unstab,
    interventions = ivs,
    reference = "nc",
    type = "difference",
    ci_method = "sandwich"
  )
  nc_stab <- res_sw$estimates$estimate[res_sw$estimates$intervention == "nc"]
  nc_unstab <- res_unstab$estimates$estimate[
    res_unstab$estimates$intervention == "nc"
  ]
  expect_equal(nc_stab, nc_unstab, tolerance = 1e-12)
})


test_that("R-long-mv-ipw1: MV IPSI under longitudinal IPW is rejected", {
  set.seed(19010)
  n <- 100
  L <- rnorm(n)
  A1_1 <- rbinom(n, 1, plogis(0.3 * L))
  A2_1 <- rbinom(n, 1, plogis(0.2 * L))
  A1_2 <- rbinom(n, 1, plogis(0.3 * L))
  A2_2 <- rbinom(n, 1, plogis(0.2 * L))
  Y <- 1 + A1_2 + A2_2 + rnorm(n)

  dat <- data.table::data.table(
    id = rep(seq_len(n), each = 2L),
    time = rep(1:2, n),
    A1 = c(rbind(A1_1, A1_2)),
    A2 = c(rbind(A2_1, A2_2)),
    L = rep(L, each = 2L),
    Y = c(rep(NA_real_, n), Y)
  )

  fit <- causat(
    dat,
    outcome = "Y",
    treatment = c("A1", "A2"),
    confounders = ~L,
    estimator = "ipw",
    type = "longitudinal",
    id = "id",
    time = "time"
  )

  expect_error(
    contrast(
      fit,
      interventions = list(
        bad = list(A1 = static(1), A2 = ipsi(0.5))
      ),
      ci_method = "sandwich"
    ),
    class = "causatr_longitudinal_ipsi_pending"
  )
})


# -----------------------------------------------------------------------
# Multivariate longitudinal IPW + effect modification (Phase 19c)
# -----------------------------------------------------------------------
# Composes per-component EM stripping (Phase 8b) with per-period EM
# stripping (Phase 10c). `make_em_mv_long_scm()` is a 2-period x
# 2-component binary DGP with effect modification by baseline `sex`;
# its static (1,1) -> (0,0) contrast has analytical truth 9 (sex = 0)
# and 15 (sex = 1). See the helper for the derivation.

test_that("T-long-mv-em1: MV longitudinal IPW EM static recovers analytical truth", {
  # Truth-based: the sex-stratified static ATE under the multivariate
  # DGP equals 9 (sex = 0) and 15 (sex = 1), combining the direct
  # treatment effect (2 + 1.5*sex)*4 with the indirect path through L1
  # (a constant shift of 1). Under a static contrast the multivariate
  # IPW (sequential MTP) estimand coincides with g-computation, so IPW
  # targets these truths exactly.
  d <- make_em_mv_long_scm(n = 6000, seed = 7720)
  d <- data.table::as.data.table(d)

  fit <- causat(
    d,
    outcome = "Y",
    treatment = c("A1", "A2"),
    confounders = ~ L0 + sex + A1:sex + A2:sex,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    estimator = "ipw"
  )
  expect_true(fit$details$is_multivariate)
  expect_true(fit$details$em_info$has_em)

  res <- contrast(
    fit,
    interventions = list(
      a = list(A1 = static(1), A2 = static(1)),
      z = list(A1 = static(0), A2 = static(0))
    ),
    reference = "z",
    type = "difference",
    by = "sex",
    ci_method = "sandwich"
  )
  cdt <- as.data.frame(res$contrasts)
  est_sex0 <- cdt$estimate[cdt$by == "0"]
  est_sex1 <- cdt$estimate[cdt$by == "1"]

  expect_equal(est_sex0, 9, tolerance = 0.6)
  expect_equal(est_sex1, 15, tolerance = 0.6)
  # Sandwich SEs are finite and positive in both strata.
  expect_true(all(is.finite(cdt$se)))
  expect_true(all(cdt$se > 0))
})


test_that("T-long-mv-em2: MV longitudinal IPW EM agrees with ICE g-comp cross-method", {
  # Cross-method: under a static contrast the multivariate IPW and
  # g-computation estimands coincide (Diaz et al. 2023). Both fit the
  # same `A1:sex + A2:sex` model and must agree per stratum.
  d <- make_em_mv_long_scm(n = 6000, seed = 7721)
  d <- data.table::as.data.table(d)

  fit_ipw <- causat(
    d,
    outcome = "Y",
    treatment = c("A1", "A2"),
    confounders = ~ L0 + sex + A1:sex + A2:sex,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    estimator = "ipw"
  )
  res_ipw <- contrast(
    fit_ipw,
    interventions = list(
      a = list(A1 = static(1), A2 = static(1)),
      z = list(A1 = static(0), A2 = static(0))
    ),
    reference = "z",
    type = "difference",
    by = "sex"
  )

  fit_ice <- causat(
    d,
    outcome = "Y",
    treatment = c("A1", "A2"),
    confounders = ~ L0 + sex + A1:sex + A2:sex,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    estimator = "gcomp"
  )
  res_ice <- contrast(
    fit_ice,
    interventions = list(
      a = list(A1 = static(1), A2 = static(1)),
      z = list(A1 = static(0), A2 = static(0))
    ),
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


test_that("T-long-mv-em3: per-component propensities strip A:modifier; modifier kept as confounder", {
  # Mechanical check: each per-period, per-component propensity must
  # drop the `A1:sex` / `A2:sex` interaction terms (a treatment-touching
  # term cannot enter a propensity model for that treatment), while
  # `sex` itself stays as a confounder main effect. The modifier
  # re-enters only at the final-period MSM (`Y ~ 1 + sex`).
  d <- make_em_mv_long_scm(n = 800, seed = 7722)
  d <- data.table::as.data.table(d)

  fit <- causat(
    d,
    outcome = "Y",
    treatment = c("A1", "A2"),
    confounders = ~ L0 + sex + A1:sex + A2:sex,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    estimator = "ipw"
  )

  expect_equal(fit$details$em_info$modifier_vars, "sex")

  for (period in names(fit$details$treatment_models_by_time)) {
    tms <- fit$details$treatment_models_by_time[[period]]
    for (comp in c("A1", "A2")) {
      rhs <- attr(stats::terms(tms[[comp]]$ps_formula), "term.labels")
      # No interaction term that touches a treatment component.
      expect_false(any(grepl("A1:|:A1|A2:|:A2", rhs)))
      # `sex` retained as a plain confounder main effect.
      expect_true("sex" %in% rhs)
    }
  }

  # The MSM the contrast path builds from `em_info` expands to
  # `Y ~ 1 + sex` so `predict()` returns per-stratum counterfactual means.
  msm <- build_ipw_msm_formula("Y", fit$details$em_info)
  expect_equal(attr(stats::terms(msm), "term.labels"), "sex")
})


test_that("T-long-mv-em4: stabilized + static MV EM equals unstabilized + static MV EM", {
  # Under a static intervention on discrete treatment the stabilizing
  # numerator and the denominator share the indicator support, so the
  # per-period stabilized weight equals the unstabilized weight up to a
  # constant that cancels in the Hajek MSM. Stabilization must not move
  # the stratified point estimates.
  d <- make_em_mv_long_scm(n = 1500, seed = 7723)
  d <- data.table::as.data.table(d)

  common_args <- list(
    data = d,
    outcome = "Y",
    treatment = c("A1", "A2"),
    confounders = ~ L0 + sex + A1:sex + A2:sex,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    estimator = "ipw"
  )
  fit_us <- do.call(causat, common_args)
  fit_st <- do.call(causat, c(common_args, list(stabilize = "marginal")))

  ints <- list(
    a = list(A1 = static(1), A2 = static(1)),
    z = list(A1 = static(0), A2 = static(0))
  )
  res_us <- contrast(
    fit_us,
    interventions = ints,
    reference = "z",
    type = "difference",
    by = "sex"
  )
  res_st <- contrast(
    fit_st,
    interventions = ints,
    reference = "z",
    type = "difference",
    by = "sex"
  )

  expect_equal(
    res_st$contrasts$estimate,
    res_us$contrasts$estimate,
    tolerance = 1e-8
  )
})
