# Multivariate (vector-treatment) longitudinal IPW.
# Split from test-longitudinal-ipw.R for test-file-level parallelism.

# -----------------------------------------------------------------------
# Multivariate longitudinal IPW (Phase 19a)
# -----------------------------------------------------------------------
# T × K propensity factorisation: W_i = prod_t prod_k w_{t,k,i}
# with block-diagonal stacked sandwich.

test_that("T-long-mv-ipw1: binary MV longitudinal IPW static, sandwich vs bootstrap", {
  skip_if_fast()
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
    Y = c(rbind(rep(NA_real_, n), Y))
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

  # Cross-method truth: for a STATIC intervention the multivariate
  # sequential-MTP IPW estimand coincides with g-computation (the
  # deterministic joint transformation), so the IPW contrast must match
  # the ICE g-comp contrast on the same sample. This is a sharper check
  # than the structural psi1 + psi2 = 1.5, which the binary HT estimator
  # only approaches at large n (verified IPW = 1.20, ICE = 1.36 here).
  fit_ice <- causat(
    dat,
    outcome = "Y",
    treatment = c("A1", "A2"),
    confounders = ~L,
    estimator = "gcomp",
    id = "id",
    time = "time"
  )
  res_ice <- contrast(
    fit_ice,
    interventions = list(
      both1 = list(A1 = static(1), A2 = static(1)),
      both0 = list(A1 = static(0), A2 = static(0))
    ),
    type = "difference"
  )
  expect_equal(
    res_sw$contrasts$estimate,
    res_ice$contrasts$estimate,
    tolerance = 0.25
  )
})


test_that("T-long-mv-ipw2: continuous MV longitudinal IPW shift, sandwich vs bootstrap", {
  skip_if_fast()
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
    Y = c(rbind(rep(NA_real_, n), Y))
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
    reference = "nc",
    type = "difference",
    ci_method = "sandwich"
  )
  res_bs <- contrast(
    fit,
    interventions = list(
      up = list(A1 = shift(delta1), A2 = shift(delta2)),
      nc = NULL
    ),
    reference = "nc",
    type = "difference",
    ci_method = "bootstrap"
  )

  # Sandwich/bootstrap SE agreement
  se_ratio <- res_sw$contrasts$se / res_bs$contrasts$se
  expect_gt(se_ratio, 0.7)
  expect_lt(se_ratio, 1.4)

  # Structural truth: an additive shift on a linear-Gaussian outcome that
  # depends only on the final-period treatments adds psi1*delta1 +
  # psi2*delta2 = 0.5*0.5 + 0.3*0.3 = 0.34 (the shift is exogenous, so the
  # downstream feedback A_1 -> L,A_2 enters only through the natural
  # value and cancels in the contrast). The sequential-MTP IPW estimator
  # targets this; lmtp_sdr agrees to ~0.335 on this DGP. Verified
  # estimate = 0.359 at this seed/n.
  expect_equal(
    res_sw$contrasts$estimate,
    psi1 * delta1 + psi2 * delta2,
    tolerance = 0.1
  )
})


test_that("T-long-mv-ipw2b: continuous MV longitudinal IPW shift + trim, sandwich vs bootstrap", {
  skip_if_fast()
  # Same DGP as T-long-mv-ipw2 but with `trim < 1`, which exercises the
  # multivariate per-period trim-threshold precompute loop in
  # `make_weight_fn_longitudinal()` (the `is_mv` branch). Sandwich and
  # bootstrap both truncate per-period density ratios at the fixed
  # quantile before the cumulative product, so they agree within MC error.
  set.seed(19002)
  n <- 800
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
    Y = c(rbind(rep(NA_real_, n), Y))
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

  ivs <- list(up = list(A1 = shift(0.5), A2 = shift(0.3)), nc = NULL)
  res_sw <- contrast(
    fit,
    interventions = ivs,
    reference = "nc",
    type = "difference",
    ci_method = "sandwich",
    trim = 0.95
  )
  set.seed(19002)
  res_bs <- contrast(
    fit,
    interventions = ivs,
    reference = "nc",
    type = "difference",
    ci_method = "bootstrap",
    n_boot = 200,
    trim = 0.95
  )

  # Point estimate identical across variance methods (same fit + weights).
  expect_equal(
    res_sw$contrasts$estimate,
    res_bs$contrasts$estimate,
    tolerance = 1e-6
  )
  # Structural truth (psi1*d1 + psi2*d2 = 0.34) up to per-period trim bias.
  # Trim winsorizes each period's density ratio at its 95th percentile,
  # shrinking the shifted mean slightly toward natural course; verified
  # estimate = 0.316 at this seed/n, so a 0.12 tolerance is honest.
  expect_equal(
    res_sw$contrasts$estimate,
    0.5 * 0.5 + 0.3 * 0.3,
    tolerance = 0.12
  )
  se_ratio <- res_sw$contrasts$se / res_bs$contrasts$se
  expect_gt(se_ratio, 0.5)
  expect_lt(se_ratio, 2)
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
    Y = c(rbind(rep(NA_real_, n), Y))
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
    Y = c(rbind(rep(NA_real_, n), Y))
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
    Y = c(rbind(rep(NA_real_, n), Y))
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
  skip_if_fast()
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
    Y = c(rbind(rep(NA_real_, n), Y))
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
  # This stabilized MV shift deliberately produces extreme per-period
  # density-ratio weights (>100), so causatr's `causatr_longitudinal_seq_positivity`
  # warning fires here as incidental noise. Its correct firing-above-threshold /
  # silence-below-threshold is asserted directly in test-longitudinal-ipw.R
  # (T-long-ipw5/6); muffle that single class here -- a targeted class handler,
  # not blanket suppression, matching the inline pattern in test-mi-coverage.R --
  # so this SE-parity test stays warning-clean regardless of the 8-hour throttle
  # state. The warning is signalled once every 8 hours, so without this the test
  # leaks it on the first run of the day.
  res_sw <- withCallingHandlers(
    contrast(
      fit,
      interventions = ivs,
      reference = "nc",
      type = "difference",
      ci_method = "sandwich"
    ),
    causatr_longitudinal_seq_positivity = function(w) {
      invokeRestart("muffleWarning")
    }
  )
  res_bs <- withCallingHandlers(
    contrast(
      fit,
      interventions = ivs,
      reference = "nc",
      type = "difference",
      ci_method = "bootstrap"
    ),
    causatr_longitudinal_seq_positivity = function(w) {
      invokeRestart("muffleWarning")
    }
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
  # Univariate longitudinal IPSI is supported (R-long-ipw5); the
  # multivariate case stays rejected because Kennedy's closed form is
  # binary-univariate with no per-component density-chain shortcut.
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
    Y = c(rbind(rep(NA_real_, n), Y))
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
    class = "causatr_longitudinal_ipsi_multivariate"
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
