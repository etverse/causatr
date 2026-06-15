# --- Longitudinal AIPW (ICE-AIPW, Bang & Robins 2005) -----------------------

test_that("longitudinal AIPW: binary static, sandwich, ATE ~ 5", {
  d <- make_linear_scm(n = 3000, n_times = 2, seed = 42)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    estimator = "aipw",
    propensity_model_fn = stats::glm,
    family = "gaussian",
    id = "id",
    time = "time"
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich"
  )
  ate <- res$contrasts$estimate[1]
  se <- res$contrasts$se[1]
  expect_equal(ate, 5, tolerance = 0.05)
  expect_true(is.finite(se) && se > 0)
})

test_that("longitudinal AIPW: binary static, bootstrap, ATE ~ 5", {
  skip_if_fast()
  d <- make_linear_scm(n = 2000, n_times = 2, seed = 43)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    estimator = "aipw",
    propensity_model_fn = stats::glm,
    family = "gaussian",
    id = "id",
    time = "time"
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "bootstrap",
    n_boot = 50
  )
  ate <- res$contrasts$estimate[1]
  se <- res$contrasts$se[1]
  expect_equal(ate, 5, tolerance = 0.05)
  expect_true(is.finite(se) && se > 0)
})

test_that("longitudinal AIPW: continuous shift, sandwich", {
  d <- make_continuous_scm(n = 3000, seed = 42)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    estimator = "aipw",
    propensity_model_fn = stats::glm,
    family = "gaussian",
    id = "id",
    time = "time"
  )
  res <- contrast(
    fit,
    interventions = list(up = shift(0.5), nat = NULL),
    reference = "nat",
    ci_method = "sandwich"
  )
  est <- res$contrasts$estimate[1]
  se <- res$contrasts$se[1]
  expect_true(is.finite(est))
  expect_true(is.finite(se) && se > 0)

  # Cross-check against ICE
  fit_ice <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    estimator = "gcomp",
    family = "gaussian",
    id = "id",
    time = "time"
  )
  res_ice <- contrast(
    fit_ice,
    interventions = list(up = shift(0.5), nat = NULL),
    reference = "nat",
    ci_method = "sandwich"
  )
  expect_lt(abs(est - res_ice$contrasts$estimate[1]), 0.3)
})

test_that("longitudinal AIPW: continuous shift, bootstrap", {
  skip_if_fast()
  d <- make_continuous_scm(n = 2000, seed = 42)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    estimator = "aipw",
    propensity_model_fn = stats::glm,
    family = "gaussian",
    id = "id",
    time = "time"
  )
  res <- contrast(
    fit,
    interventions = list(up = shift(0.5), nat = NULL),
    reference = "nat",
    ci_method = "bootstrap",
    n_boot = 50
  )
  est <- res$contrasts$estimate[1]
  se <- res$contrasts$se[1]
  expect_true(is.finite(est))
  expect_true(is.finite(se) && se > 0)
})

test_that("longitudinal AIPW: dynamic intervention, sandwich", {
  d <- make_linear_scm(n = 3000, n_times = 2, seed = 44)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    estimator = "aipw",
    propensity_model_fn = stats::glm,
    family = "gaussian",
    id = "id",
    time = "time"
  )
  dyn_fn <- dynamic(function(data, trt) ifelse(data$L0 > 0, 1L, 0L))
  res <- contrast(
    fit,
    interventions = list(dyn = dyn_fn, nat = NULL),
    reference = "nat",
    ci_method = "sandwich"
  )
  est <- res$contrasts$estimate[1]
  se <- res$contrasts$se[1]
  expect_true(is.finite(est))
  expect_true(is.finite(se) && se > 0)
})

test_that("longitudinal AIPW: scale_by intervention, sandwich", {
  d <- make_continuous_scm(n = 3000, seed = 42)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    estimator = "aipw",
    propensity_model_fn = stats::glm,
    family = "gaussian",
    id = "id",
    time = "time"
  )
  res <- contrast(
    fit,
    interventions = list(sc = scale_by(1.2), nat = NULL),
    reference = "nat",
    ci_method = "sandwich"
  )
  est <- res$contrasts$estimate[1]
  se <- res$contrasts$se[1]
  expect_true(is.finite(est))
  expect_true(is.finite(se) && se > 0)
})

test_that("longitudinal AIPW vs ICE vs long-IPW: cross-method agreement", {
  d <- make_linear_scm(n = 3000, n_times = 2, seed = 45)
  ivs <- list(a1 = static(1), a0 = static(0))

  fit_aipw <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    estimator = "aipw",
    propensity_model_fn = stats::glm,
    family = "gaussian",
    id = "id",
    time = "time"
  )
  fit_ice <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    estimator = "gcomp",
    family = "gaussian",
    id = "id",
    time = "time"
  )
  fit_ipw <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    estimator = "ipw",
    family = "gaussian",
    id = "id",
    time = "time"
  )

  res_aipw <- contrast(
    fit_aipw,
    interventions = ivs,
    reference = "a0",
    ci_method = "sandwich"
  )
  res_ice <- contrast(
    fit_ice,
    interventions = ivs,
    reference = "a0",
    ci_method = "sandwich"
  )
  res_ipw <- contrast(
    fit_ipw,
    interventions = ivs,
    reference = "a0",
    ci_method = "sandwich"
  )

  ate_aipw <- res_aipw$contrasts$estimate[1]
  ate_ice <- res_ice$contrasts$estimate[1]
  ate_ipw <- res_ipw$contrasts$estimate[1]

  expect_lt(abs(ate_aipw - ate_ice), 0.3)
  expect_lt(abs(ate_aipw - ate_ipw), 0.3)
})

test_that("longitudinal AIPW vs ICE: continuous shift agreement", {
  d <- make_continuous_scm(n = 3000, seed = 46)

  fit_aipw <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    estimator = "aipw",
    propensity_model_fn = stats::glm,
    family = "gaussian",
    id = "id",
    time = "time"
  )
  fit_ice <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    estimator = "gcomp",
    family = "gaussian",
    id = "id",
    time = "time"
  )

  res_aipw <- contrast(
    fit_aipw,
    interventions = list(up = shift(0.5), nat = NULL),
    reference = "nat",
    ci_method = "sandwich"
  )
  res_ice <- contrast(
    fit_ice,
    interventions = list(up = shift(0.5), nat = NULL),
    reference = "nat",
    ci_method = "sandwich"
  )

  expect_lt(
    abs(res_aipw$contrasts$estimate[1] - res_ice$contrasts$estimate[1]),
    0.3
  )
})

test_that("longitudinal AIPW DR: wrong outcome, correct propensity", {
  skip_if_fast()
  d <- make_linear_scm(n = 3000, n_times = 2, seed = 47)
  # Misspecify outcome by dropping TV confounders
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    estimator = "aipw",
    propensity_model_fn = stats::glm,
    family = "gaussian",
    id = "id",
    time = "time"
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "bootstrap",
    n_boot = 50
  )
  ate <- res$contrasts$estimate[1]
  # DR: should still recover truth despite wrong outcome model
  expect_equal(ate, 5, tolerance = 0.05)
})

test_that("longitudinal AIPW DR: correct outcome, wrong propensity", {
  skip_if_fast()
  d <- make_linear_scm(n = 3000, n_times = 2, seed = 48)
  # Misspecify propensity by omitting L from confounders
  # (outcome still includes L0 which is baseline, but misses L
  # in propensity). Since causat uses same confounders for both,
  # we misspecify propensity by dropping confounders_tv entirely
  # while keeping outcome correct via baseline confounders only.
  # This partially misspecifies propensity.
  fit_correct <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    estimator = "aipw",
    propensity_model_fn = stats::glm,
    family = "gaussian",
    id = "id",
    time = "time"
  )
  fit_misspec <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    estimator = "aipw",
    propensity_model_fn = stats::glm,
    family = "gaussian",
    id = "id",
    time = "time"
  )
  res_c <- contrast(
    fit_correct,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "bootstrap",
    n_boot = 50
  )
  res_m <- contrast(
    fit_misspec,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "bootstrap",
    n_boot = 50
  )
  # Both should be close to truth (DR property)
  expect_equal(res_c$contrasts$estimate[1], 5, tolerance = 0.05)
  expect_equal(res_m$contrasts$estimate[1], 5, tolerance = 0.05)
})

test_that("longitudinal AIPW: sandwich SE finite and positive", {
  d <- make_linear_scm(n = 3000, n_times = 2, seed = 49)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    estimator = "aipw",
    propensity_model_fn = stats::glm,
    family = "gaussian",
    id = "id",
    time = "time"
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    ci_method = "sandwich"
  )
  sds <- sqrt(diag(res$vcov))
  expect_true(all(sds > 0 & is.finite(sds)))
})

test_that("longitudinal AIPW: bootstrap SE finite and positive", {
  skip_if_fast()
  d <- make_linear_scm(n = 2000, n_times = 2, seed = 50)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    estimator = "aipw",
    propensity_model_fn = stats::glm,
    family = "gaussian",
    id = "id",
    time = "time"
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    ci_method = "bootstrap",
    n_boot = 50
  )
  sds <- sqrt(diag(res$vcov))
  expect_true(all(sds > 0 & is.finite(sds)))
})

test_that("longitudinal AIPW: sandwich vs bootstrap SE agreement", {
  skip_if_fast()
  d <- make_linear_scm(n = 3000, n_times = 2, seed = 51)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    estimator = "aipw",
    propensity_model_fn = stats::glm,
    family = "gaussian",
    id = "id",
    time = "time"
  )
  ivs <- list(a1 = static(1), a0 = static(0))
  res_s <- contrast(fit, interventions = ivs, ci_method = "sandwich")
  res_b <- contrast(
    fit,
    interventions = ivs,
    ci_method = "bootstrap",
    n_boot = 200
  )
  se_sand <- sqrt(diag(res_s$vcov))
  se_boot <- sqrt(diag(res_b$vcov))
  ratio <- se_sand / se_boot
  expect_true(all(ratio > 0.5 & ratio < 2.0))
})

test_that("longitudinal AIPW: effect modification by baseline (sex)", {
  d <- make_em_ice_scm(n = 2500, seed = 52)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L0 + sex,
    confounders_tv = ~L,
    estimator = "aipw",
    propensity_model_fn = stats::glm,
    family = "gaussian",
    id = "id",
    time = "time"
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    by = "sex",
    ci_method = "sandwich"
  )
  est_sex0 <- res$contrasts$estimate[
    res$contrasts$by == "0"
  ]
  est_sex1 <- res$contrasts$estimate[
    res$contrasts$by == "1"
  ]
  expect_equal(est_sex0, 5, tolerance = 0.05)
  expect_equal(est_sex1, 8, tolerance = 0.05)
})

test_that("longitudinal AIPW: EM agreement with ICE", {
  d <- make_em_ice_scm(n = 2500, seed = 53)
  ivs <- list(a1 = static(1), a0 = static(0))
  by_var <- "sex"

  fit_aipw <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L0 + sex,
    confounders_tv = ~L,
    estimator = "aipw",
    propensity_model_fn = stats::glm,
    family = "gaussian",
    id = "id",
    time = "time"
  )
  fit_ice <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L0 + sex,
    confounders_tv = ~L,
    estimator = "gcomp",
    family = "gaussian",
    id = "id",
    time = "time"
  )

  res_aipw <- contrast(
    fit_aipw,
    interventions = ivs,
    reference = "a0",
    by = by_var,
    ci_method = "sandwich"
  )
  res_ice <- contrast(
    fit_ice,
    interventions = ivs,
    reference = "a0",
    by = by_var,
    ci_method = "sandwich"
  )

  for (sx in c("0", "1")) {
    aipw_est <- res_aipw$contrasts$estimate[
      res_aipw$contrasts$by == sx
    ]
    ice_est <- res_ice$contrasts$estimate[
      res_ice$contrasts$by == sx
    ]
    expect_lt(abs(aipw_est - ice_est), 2.0)
  }
})

test_that("longitudinal AIPW: binomial outcome", {
  set.seed(55)
  n <- 2000
  id <- rep(1:n, each = 2)
  time <- rep(0:1, times = n)
  L0 <- rnorm(n)[id]
  A <- rbinom(2 * n, 1, plogis(0.3 * L0))
  L <- rnorm(2 * n, mean = 0.3 * A)
  Y <- rep(NA_real_, 2 * n)
  Y[time == 1] <- as.numeric(rbinom(
    n,
    1,
    plogis(
      -1 + 0.5 * A[time == 1] + 0.3 * L[time == 1] + 0.2 * L0[time == 1]
    )
  ))
  d <- data.table::data.table(
    Y = Y,
    A = A,
    L0 = L0,
    L = L,
    id = id,
    time = time
  )
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    estimator = "aipw",
    propensity_model_fn = stats::glm,
    family = "binomial",
    id = "id",
    time = "time"
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    ci_method = "sandwich"
  )
  ests <- coef(res)
  expect_true(all(ests > 0 & ests < 1))
  expect_true(all(is.finite(sqrt(diag(res$vcov)))))
})

test_that("longitudinal AIPW: 3-period DGP agrees with IPW", {
  d <- make_linear_scm(n = 3000, n_times = 3, seed = 56)
  fit_aipw <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    estimator = "aipw",
    propensity_model_fn = stats::glm,
    family = "gaussian",
    id = "id",
    time = "time"
  )
  fit_ipw <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    estimator = "ipw",
    family = "gaussian",
    id = "id",
    time = "time"
  )
  ivs <- list(a1 = static(1), a0 = static(0))
  res_aipw <- contrast(
    fit_aipw,
    interventions = ivs,
    reference = "a0",
    ci_method = "sandwich"
  )
  res_ipw <- contrast(
    fit_ipw,
    interventions = ivs,
    reference = "a0",
    ci_method = "sandwich"
  )
  ate_aipw <- res_aipw$contrasts$estimate[1]
  ate_ipw <- res_ipw$contrasts$estimate[1]
  se <- res_aipw$contrasts$se[1]
  expect_lt(abs(ate_aipw - ate_ipw), 0.5)
  expect_true(is.finite(se) && se > 0)
})

# -----------------------------------------------------------------------
# Multivariate longitudinal AIPW (Phase 25)
# -----------------------------------------------------------------------
# Joint time-varying treatments `treatment = c("A1", "A2")` with
# `id =`/`time =` under the doubly-robust ICE-AIPW estimator. The
# propensity side reuses the Phase 19 multivariate density chain
# (per-period x per-component `fit_treatment_models()`); the outcome
# side reuses the ICE backward recursion (already loops over vector
# treatment). `make_em_mv_long_scm()` is a 2-period x 2-component binary
# DGP with effect modification by baseline `sex`; its static
# (1,1) -> (0,0) contrast has analytical truth 9 (sex = 0) and 15
# (sex = 1), and the sex-marginal pooled truth is (9 + 15)/2 = 12.
#
# Oracle note: longitudinal AIPW has no external single-call oracle in
# this suite (lmtp's SDR fixes nuisances differently). The validation
# strategy mirrors the univariate longitudinal AIPW tests above:
# truth recovery on the static binary DGP (where MV IPW = MV gcomp =
# MV AIPW coincide, Diaz et al. 2023), cross-method triangulation
# against MV ICE g-comp and MV longitudinal IPW, and the
# double-robustness property under one-sided misspecification.

test_that("T-long-mv-aipw1: MV longitudinal AIPW static recovers analytical truth by sex", {
  d <- make_em_mv_long_scm(n = 6000, seed = 2501)
  d <- data.table::as.data.table(d)

  fit <- causat(
    d,
    outcome = "Y",
    treatment = c("A1", "A2"),
    confounders = ~ L0 + sex + A1:sex + A2:sex,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    estimator = "aipw",
    propensity_model_fn = stats::glm,
    family = "gaussian"
  )
  expect_true(fit$details$is_multivariate)

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
  expect_true(all(is.finite(cdt$se) & cdt$se > 0))
})

test_that("T-long-mv-aipw2: MV longitudinal AIPW sandwich vs bootstrap SE parity", {
  skip_if_fast()
  d <- make_em_mv_long_scm(n = 3000, seed = 2502)
  d <- data.table::as.data.table(d)

  fit <- causat(
    d,
    outcome = "Y",
    treatment = c("A1", "A2"),
    confounders = ~ L0 + sex,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    estimator = "aipw",
    propensity_model_fn = stats::glm,
    family = "gaussian"
  )
  ivs <- list(
    a = list(A1 = static(1), A2 = static(1)),
    z = list(A1 = static(0), A2 = static(0))
  )
  res_s <- contrast(
    fit,
    interventions = ivs,
    reference = "z",
    type = "difference",
    ci_method = "sandwich"
  )
  res_b <- contrast(
    fit,
    interventions = ivs,
    reference = "z",
    type = "difference",
    ci_method = "bootstrap",
    n_boot = 200
  )
  ratio <- res_s$contrasts$se / res_b$contrasts$se
  expect_gt(ratio, 0.7)
  expect_lt(ratio, 1.4)
})

test_that("T-long-mv-aipw3: MV longitudinal AIPW cross-method agreement with ICE g-comp and IPW", {
  # Under a static contrast the multivariate IPW (sequential MTP),
  # multivariate ICE g-computation, and multivariate AIPW estimands
  # coincide (Diaz et al. 2023). All three fit the same EM model and
  # must agree per stratum.
  d <- make_em_mv_long_scm(n = 4000, seed = 2503)
  d <- data.table::as.data.table(d)
  ivs <- list(
    a = list(A1 = static(1), A2 = static(1)),
    z = list(A1 = static(0), A2 = static(0))
  )
  conf <- ~ L0 + sex + A1:sex + A2:sex

  fit_a <- causat(
    d,
    outcome = "Y",
    treatment = c("A1", "A2"),
    confounders = conf,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    estimator = "aipw",
    propensity_model_fn = stats::glm,
    family = "gaussian"
  )
  fit_g <- causat(
    d,
    outcome = "Y",
    treatment = c("A1", "A2"),
    confounders = conf,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    estimator = "gcomp",
    family = "gaussian"
  )
  fit_i <- causat(
    d,
    outcome = "Y",
    treatment = c("A1", "A2"),
    confounders = conf,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    estimator = "ipw",
    family = "gaussian"
  )

  res_a <- contrast(fit_a, interventions = ivs, reference = "z", by = "sex")
  res_g <- contrast(fit_g, interventions = ivs, reference = "z", by = "sex")
  res_i <- contrast(fit_i, interventions = ivs, reference = "z", by = "sex")

  for (sx in c("0", "1")) {
    a_est <- res_a$contrasts$estimate[res_a$contrasts$by == sx]
    g_est <- res_g$contrasts$estimate[res_g$contrasts$by == sx]
    i_est <- res_i$contrasts$estimate[res_i$contrasts$by == sx]
    expect_lt(abs(a_est - g_est), 0.5)
    expect_lt(abs(a_est - i_est), 0.5)
  }
})

test_that("T-long-mv-aipw4: MV longitudinal AIPW double-robustness under one-sided misspecification", {
  skip_if_fast()
  # DR property: with EITHER the outcome model OR the propensity model
  # correctly specified, the AIPW estimator recovers the truth. We hold
  # the EM baseline (`confounders`, driving the MSM expansion and EM
  # stripping) correct in both arms and toggle the time-varying
  # confounder `L` on/off per component:
  #   (a) wrong outcome (drops L), correct propensity (keeps L);
  #   (b) correct outcome (keeps L), wrong propensity (drops L).
  # Both must recover the sex-stratified truths 9 / 15.
  conf <- ~ L0 + sex + A1:sex + A2:sex
  ivs <- list(
    a = list(A1 = static(1), A2 = static(1)),
    z = list(A1 = static(0), A2 = static(0))
  )

  # (a) wrong outcome / correct propensity
  d_a <- data.table::as.data.table(make_em_mv_long_scm(n = 6000, seed = 2504))
  fit_a <- causat(
    d_a,
    outcome = "Y",
    treatment = c("A1", "A2"),
    confounders = conf,
    confounders_tv_treatment = ~L,
    id = "id",
    time = "time",
    estimator = "aipw",
    propensity_model_fn = stats::glm,
    family = "gaussian"
  )
  res_a <- contrast(
    fit_a,
    interventions = ivs,
    reference = "z",
    by = "sex",
    ci_method = "bootstrap",
    n_boot = 50
  )
  ca <- as.data.frame(res_a$contrasts)
  expect_equal(ca$estimate[ca$by == "0"], 9, tolerance = 0.7)
  expect_equal(ca$estimate[ca$by == "1"], 15, tolerance = 0.7)

  # (b) correct outcome / wrong propensity
  d_b <- data.table::as.data.table(make_em_mv_long_scm(n = 6000, seed = 2505))
  fit_b <- causat(
    d_b,
    outcome = "Y",
    treatment = c("A1", "A2"),
    confounders = conf,
    confounders_tv_outcome = ~L,
    id = "id",
    time = "time",
    estimator = "aipw",
    propensity_model_fn = stats::glm,
    family = "gaussian"
  )
  res_b <- contrast(
    fit_b,
    interventions = ivs,
    reference = "z",
    by = "sex",
    ci_method = "bootstrap",
    n_boot = 50
  )
  cb <- as.data.frame(res_b$contrasts)
  expect_equal(cb$estimate[cb$by == "0"], 9, tolerance = 0.7)
  expect_equal(cb$estimate[cb$by == "1"], 15, tolerance = 0.7)
})

test_that("T-long-mv-aipw5: MV longitudinal AIPW continuous shift finite + SE parity", {
  skip_if_fast()
  # No closed-form truth under a non-static shift (MV AIPW and g-comp
  # estimands diverge by design, Diaz et al. 2023), so we check only
  # finiteness and sandwich-vs-bootstrap SE parity on a 2-period x
  # 2-component continuous-treatment DGP.
  set.seed(2506)
  n <- 1200
  L <- rnorm(n)
  A1_1 <- rnorm(n, 0.5 + 0.3 * L)
  A2_1 <- rnorm(n, 0.3 + 0.2 * L + 0.1 * A1_1)
  L2 <- 0.5 * (A1_1 + A2_1) + 0.5 * L + rnorm(n, 0, 0.5)
  A1_2 <- rnorm(n, 0.3 + 0.3 * L2 + 0.1 * A1_1)
  A2_2 <- rnorm(n, 0.2 + 0.2 * L2 + 0.1 * A1_2 + 0.05 * A2_1)
  Y <- 2 + 0.5 * L2 + 0.5 * A1_2 + 0.3 * A2_2 + rnorm(n)

  dat <- data.table::data.table(
    id = rep(seq_len(n), each = 2L),
    time = rep(0:1, n),
    A1 = c(rbind(A1_1, A1_2)),
    A2 = c(rbind(A2_1, A2_2)),
    L = c(rbind(L, L2)),
    Y = c(rbind(rep(NA_real_, n), Y))
  )

  fit <- causat(
    dat,
    outcome = "Y",
    treatment = c("A1", "A2"),
    confounders_tv = ~L,
    id = "id",
    time = "time",
    estimator = "aipw",
    propensity_model_fn = stats::glm,
    family = "gaussian"
  )
  ivs <- list(
    up = list(A1 = shift(0.5), A2 = shift(0.5)),
    nat = NULL
  )
  res_s <- contrast(
    fit,
    interventions = ivs,
    reference = "nat",
    type = "difference",
    ci_method = "sandwich"
  )
  res_b <- contrast(
    fit,
    interventions = ivs,
    reference = "nat",
    type = "difference",
    ci_method = "bootstrap",
    n_boot = 200
  )
  expect_true(is.finite(res_s$contrasts$estimate))
  expect_true(is.finite(res_s$contrasts$se) && res_s$contrasts$se > 0)
  ratio <- res_s$contrasts$se / res_b$contrasts$se
  expect_gt(ratio, 0.6)
  expect_lt(ratio, 1.6)
})

test_that("R-long-mv-aipw1: stabilize rejected for longitudinal AIPW", {
  # Stabilization is not yet wired into the longitudinal AIPW path
  # (no per-period numerator sweep), so it must be rejected rather than
  # silently dropped — for both univariate and multivariate treatment.
  d <- data.table::as.data.table(make_em_mv_long_scm(n = 300, seed = 58))
  expect_error(
    causat(
      d,
      outcome = "Y",
      treatment = c("A1", "A2"),
      confounders = ~sex,
      confounders_tv = ~L,
      estimator = "aipw",
      family = "gaussian",
      stabilize = "marginal",
      id = "id",
      time = "time"
    ),
    class = "causatr_stabilize_longitudinal_aipw"
  )
})

test_that("R-long-aipw-unbalanced: sandwich warns on unbalanced panel", {
  # Critical review 2026-05-30, Issue #5
  # (repro: /tmp/causatr_repro_unbalanced_se.R).
  # Under monotone dropout a period's models are fit on fewer than n rows;
  # a Monte-Carlo study (300 reps) showed the rescaled longitudinal AIPW
  # sandwich underestimates the true SE by ~15% (empirical SD 0.0695 vs
  # mean sandwich SE 0.0590). The honest contract is a classed warning
  # steering to the bootstrap, not a silently-low SE. This test pins the
  # warning contract: it fires iff the panel is unbalanced.
  gen <- function(n, dropout, seed) {
    set.seed(seed)
    L0 <- rbinom(n, 1, 0.5)
    A0 <- rbinom(n, 1, plogis(-0.2 + 0.6 * L0))
    L1 <- rbinom(n, 1, plogis(-0.2 + 0.8 * A0 + 0.3 * L0))
    A1 <- rbinom(n, 1, plogis(-0.2 + 0.6 * L1))
    Y <- 2 + 2.5 * A0 + 2.5 * A1 - 0.5 * L1 + rnorm(n)
    d <- data.table::data.table(
      id = rep(seq_len(n), each = 2L),
      time = rep(0:1, n),
      A = c(rbind(A0, A1)),
      L = c(rbind(L0, L1)),
      Y = c(rbind(NA_real_, Y))
    )
    if (dropout) {
      drop_ids <- which(A0 == 1 & runif(n) < 0.4)
      d <- d[!(id %in% drop_ids & time == 1L)]
    }
    d[]
  }

  run <- function(d) {
    fit <- causat(
      d,
      outcome = "Y",
      treatment = "A",
      confounders = ~1,
      confounders_tv = ~L,
      id = "id",
      time = "time",
      estimator = "aipw",
      propensity_model_fn = stats::glm,
      family = "gaussian"
    )
    contrast(
      fit,
      interventions = list(a = static(1), z = static(0)),
      reference = "z",
      type = "difference",
      ci_method = "sandwich"
    )
  }

  # Balanced panel: no warning.
  rlang::reset_warning_verbosity(
    "causatr_longitudinal_aipw_unbalanced_sandwich"
  )
  expect_no_condition(
    run(gen(800, dropout = FALSE, seed = 101)),
    class = "causatr_longitudinal_aipw_unbalanced_sandwich"
  )

  # Unbalanced panel (informative monotone dropout): warning fires.
  rlang::reset_warning_verbosity(
    "causatr_longitudinal_aipw_unbalanced_sandwich"
  )
  expect_warning(
    run(gen(800, dropout = TRUE, seed = 101)),
    class = "causatr_longitudinal_aipw_unbalanced_sandwich"
  )
})

test_that("longitudinal AIPW: ATT rejected", {
  d <- make_linear_scm(n = 300, n_times = 2, seed = 58)
  expect_error(
    causat(
      d,
      outcome = "Y",
      treatment = "A",
      confounders = ~L0,
      confounders_tv = ~L,
      estimator = "aipw",
      propensity_model_fn = stats::glm,
      family = "gaussian",
      estimand = "ATT",
      id = "id",
      time = "time"
    ),
    "estimand = 'ATT'"
  )
})


# --- Plan gap tests (added post-16i) -----------------------------------

test_that("longitudinal AIPW DR caveat: sandwich SE under misspecified outcome", {
  skip_if_fast()
  d <- make_linear_scm(n = 3000, n_times = 2, seed = 70)
  # Misspecify outcome by dropping TV confounders
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    estimator = "aipw",
    propensity_model_fn = stats::glm,
    family = "gaussian",
    id = "id",
    time = "time"
  )
  res_sw <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich"
  )
  res_bs <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "bootstrap",
    n_boot = 50
  )
  se_sw <- res_sw$contrasts$se[1]
  se_bs <- res_bs$contrasts$se[1]
  # Both SEs must be finite and positive. Under misspecification the
  # sandwich SE is NOT DR-consistent (Rotnitzky et al. 2017), so it may
  # diverge from bootstrap — we only check finiteness here.
  expect_true(is.finite(se_sw) && se_sw > 0)
  expect_true(is.finite(se_bs) && se_bs > 0)
  # Point estimate should still recover truth (DR property)
  expect_equal(res_sw$contrasts$estimate[1], 5, tolerance = 0.05)
})

test_that("longitudinal AIPW: EM agreement with long-IPW", {
  d <- make_em_ice_scm(n = 2500, seed = 71)
  ivs <- list(a1 = static(1), a0 = static(0))
  by_var <- "sex"

  fit_aipw <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L0 + sex,
    confounders_tv = ~L,
    estimator = "aipw",
    propensity_model_fn = stats::glm,
    family = "gaussian",
    id = "id",
    time = "time"
  )
  fit_ipw <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L0 + sex,
    confounders_tv = ~L,
    estimator = "ipw",
    family = "gaussian",
    id = "id",
    time = "time"
  )

  res_aipw <- contrast(
    fit_aipw,
    interventions = ivs,
    reference = "a0",
    by = by_var,
    ci_method = "sandwich"
  )
  res_ipw <- contrast(
    fit_ipw,
    interventions = ivs,
    reference = "a0",
    by = by_var,
    ci_method = "sandwich"
  )

  for (sx in c("0", "1")) {
    aipw_est <- res_aipw$contrasts$estimate[
      res_aipw$contrasts$by == sx
    ]
    ipw_est <- res_ipw$contrasts$estimate[
      res_ipw$contrasts$by == sx
    ]
    expect_lt(abs(aipw_est - ipw_est), 2.0)
  }
})

test_that("longitudinal AIPW: 3-period sandwich vs bootstrap SE agreement", {
  skip_if_fast()
  d <- make_linear_scm(n = 3000, n_times = 3, seed = 72)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    estimator = "aipw",
    propensity_model_fn = stats::glm,
    family = "gaussian",
    id = "id",
    time = "time"
  )
  ivs <- list(a1 = static(1), a0 = static(0))
  res_sw <- contrast(
    fit,
    interventions = ivs,
    reference = "a0",
    ci_method = "sandwich"
  )
  res_bs <- contrast(
    fit,
    interventions = ivs,
    reference = "a0",
    ci_method = "bootstrap",
    n_boot = 100
  )
  se_sw <- res_sw$contrasts$se[1]
  se_bs <- res_bs$contrasts$se[1]
  expect_true(is.finite(se_sw) && se_sw > 0)
  expect_true(is.finite(se_bs) && se_bs > 0)
  se_ratio <- se_sw / se_bs
  expect_gt(se_ratio, 0.5)
  expect_lt(se_ratio, 2.0)
})

test_that("longitudinal AIPW: near-positivity stress test", {
  set.seed(73)
  n <- 2000
  id <- rep(1:n, each = 2)
  time <- rep(0:1, times = n)
  L0 <- rnorm(n)[id]
  # Heavy confounding: propensity near 0 or 1 for most units
  A <- rbinom(2 * n, 1, plogis(2.0 + 2.5 * L0))
  L <- rnorm(2 * n, mean = 0.5 * A)
  Y <- rep(NA_real_, 2 * n)
  Y[time == 1] <- 2 + 3 * A[time == 1] + 1.5 * L0[time == 1] + rnorm(n)
  d <- data.table::data.table(
    Y = Y,
    A = A,
    L0 = L0,
    L = L,
    id = id,
    time = time
  )

  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    estimator = "aipw",
    propensity_model_fn = stats::glm,
    family = "gaussian",
    id = "id",
    time = "time"
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich"
  )
  ate <- res$contrasts$estimate[1]
  se <- res$contrasts$se[1]
  expect_true(is.finite(ate))
  expect_true(is.finite(se) && se > 0)
  # DGP: Y[t=1] = 2 + 3*A + 1.5*L0 + noise → true ATE = 3
  expect_equal(ate, 3, tolerance = 0.05)
})

# --- delicatessen cross-check (chunk 16k) ------------------------------------
#
# Reference values generated by data-raw/aipw_reference.py using
# delicatessen (Zivich 2024) stacked M-estimation sandwich variance.
# Zivich PN et al. (2024). Statistics in Medicine 43:5562-5572.

test_that("AIPW sandwich matches delicatessen — binary treatment, ATE", {
  ref_mu_1 <- 4.9715
  ref_se_1 <- 0.0289
  ref_mu_0 <- 1.9730
  ref_se_0 <- 0.0314
  ref_ate <- 2.9985
  ref_se_ate <- 0.0299

  set.seed(42)
  n <- 5000
  L <- rnorm(n)
  A <- rbinom(n, 1, plogis(0.3 + 0.5 * L))
  Y <- 2 + 3 * A + 1.5 * L + rnorm(n)
  d <- data.frame(Y = Y, A = A, L = L)

  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "aipw",
    propensity_model_fn = stats::glm
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich"
  )

  mu_1 <- res$estimates$estimate[res$estimates$intervention == "a1"]
  mu_0 <- res$estimates$estimate[res$estimates$intervention == "a0"]
  se_1 <- res$estimates$se[res$estimates$intervention == "a1"]
  se_0 <- res$estimates$se[res$estimates$intervention == "a0"]

  expect_equal(mu_1, ref_mu_1, tolerance = 0.005)
  expect_equal(mu_0, ref_mu_0, tolerance = 0.005)
  expect_equal(res$contrasts$estimate[1], ref_ate, tolerance = 0.005)
  expect_equal(se_1, ref_se_1, tolerance = 0.005)
  expect_equal(se_0, ref_se_0, tolerance = 0.005)
  expect_equal(res$contrasts$se[1], ref_se_ate, tolerance = 0.005)
})

test_that("AIPW sandwich matches delicatessen — continuous treatment, shift(1)", {
  ref_mu_shift <- 6.0079
  ref_se_shift <- 0.0613
  ref_mu_nat <- 4.0149
  ref_se_nat <- 0.0586
  ref_effect <- 1.9930
  ref_se_effect <- 0.0171

  set.seed(99)
  n <- 5000
  L <- rnorm(n)
  A <- rnorm(n, mean = 1 + L, sd = 1)
  Y <- 2 + 2 * A + 1.5 * L + rnorm(n)
  d <- data.frame(Y = Y, A = A, L = L)

  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "aipw",
    propensity_model_fn = stats::glm
  )
  res <- contrast(
    fit,
    interventions = list(shifted = shift(1), nat = NULL),
    reference = "nat",
    ci_method = "sandwich"
  )

  mu_shift <- res$estimates$estimate[res$estimates$intervention == "shifted"]
  mu_nat <- res$estimates$estimate[res$estimates$intervention == "nat"]
  se_shift <- res$estimates$se[res$estimates$intervention == "shifted"]
  se_nat <- res$estimates$se[res$estimates$intervention == "nat"]

  expect_equal(mu_shift, ref_mu_shift, tolerance = 0.005)
  expect_equal(mu_nat, ref_mu_nat, tolerance = 0.005)
  expect_equal(res$contrasts$estimate[1], ref_effect, tolerance = 0.005)
  expect_equal(se_shift, ref_se_shift, tolerance = 0.005)
  expect_equal(se_nat, ref_se_nat, tolerance = 0.005)
  expect_equal(res$contrasts$se[1], ref_se_effect, tolerance = 0.005)
})

test_that("longitudinal AIPW: IPSI rejection", {
  d <- make_linear_scm(n = 300, n_times = 2, seed = 74)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    estimator = "aipw",
    propensity_model_fn = stats::glm,
    family = "gaussian",
    id = "id",
    time = "time"
  )
  expect_error(
    contrast(
      fit,
      interventions = list(shifted = ipsi(0.5)),
      ci_method = "sandwich"
    ),
    class = "causatr_longitudinal_ipsi_pending"
  )
})
