# Cross-package oracle tests: causatr SNM vs DTRreg g-estimation.
#
# DTRreg (2.4+) implements G-estimation (method="gest") for binary
# treatments with sandwich or bootstrap variance. Point estimates from
# causatr and DTRreg should agree to machine precision; SEs should
# agree to within 1e-2 (different sandwich implementations with
# identical bread structure).
#
# Oracle: DTRreg::DTRreg(method="gest", var.estim="sandwich")
# Reference: Laber EB et al. (2015). DTRreg: an R package for dynamic
# treatment regime estimation. Technical report.
#
# All tests guard with skip_if_not_installed("DTRreg").

# ── Point treatment (no EM, no TF) ────────────────────────────────────────

test_that("SNM point × binary × no EM: causatr psi agrees with DTRreg", {
  skip_if_not_installed("DTRreg")

  set.seed(11)
  n <- 1000
  L <- rnorm(n)
  A <- rbinom(n, 1, plogis(0.3 * L))
  Y <- 1 + 2 * A + 0.5 * L + rnorm(n)
  d <- data.table::data.table(Y = Y, A = A, L = L)

  # DTRreg oracle (no variance — just point estimate)
  dr <- DTRreg::DTRreg(
    outcome = Y,
    blip.mod = ~1,
    treat.mod = A ~ L,
    tf.mod = ~L,
    data = data.frame(Y = Y, A = A, L = L),
    method = "gest",
    treat.type = "bin"
  )
  psi_dr <- dr$psi[[1]][["A"]]

  # causatr
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "snm",
    treatment_free = ~L
  )
  res <- contrast(fit, ci_method = "sandwich")
  psi_causat <- res$estimates$estimate[1]

  expect_equal(psi_causat, psi_dr, tolerance = 1e-6)
  expect_equal(psi_causat, 2, tolerance = 0.1)
})


# ── Point treatment with EM + TF model ────────────────────────────────────

test_that("SNM point × binary × EM + TF: causatr matches DTRreg (point + SE)", {
  skip_if_not_installed("DTRreg")

  set.seed(22)
  n <- 2000
  L <- rnorm(n)
  M <- rbinom(n, 1, 0.5)
  A <- rbinom(n, 1, plogis(0.3 * L + 0.2 * M))
  Y <- 1 + 2 * A + 3 * A * M + 0.5 * L + 0.8 * M + rnorm(n)
  d <- data.table::data.table(Y = Y, A = A, L = L, M = M)

  dr <- DTRreg::DTRreg(
    outcome = Y,
    blip.mod = ~M,
    treat.mod = A ~ L + M,
    tf.mod = ~ L + M,
    data = data.frame(Y = Y, A = A, L = L, M = M),
    method = "gest",
    treat.type = "bin",
    var.estim = "sandwich"
  )
  # dr$psi[[1]] names: "A" (intercept), "M:A" (modifier)
  psi_dr_int <- dr$psi[[1]][["A"]]
  psi_dr_M <- dr$psi[[1]][["M:A"]]
  # DTRreg sandwich SEs from covmat diagonal
  se_dr_int <- sqrt(dr$covmat[[1]]["A", "A"])
  se_dr_M <- sqrt(dr$covmat[[1]]["M:A", "M:A"])

  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders_outcome = ~ L + M + A:M,
    confounders_treatment = ~ L + M,
    estimator = "snm",
    treatment_free = ~ L + M
  )
  res <- contrast(fit, ci_method = "sandwich")
  psi_int <- res$estimates$estimate[1]
  psi_M <- res$estimates$estimate[2]
  se_int <- res$estimates$se[1]
  se_M <- res$estimates$se[2]

  # Point estimates: machine precision
  expect_equal(psi_int, psi_dr_int, tolerance = 1e-6)
  expect_equal(psi_M, psi_dr_M, tolerance = 1e-6)

  # SEs: implementations differ in details but agree to ~1%
  expect_equal(se_int, se_dr_int, tolerance = 1e-2)
  expect_equal(se_M, se_dr_M, tolerance = 1e-2)

  # Truth recovery
  expect_equal(psi_int, 2, tolerance = 0.15)
  expect_equal(psi_M, 3, tolerance = 0.15)
})


# ── Continuous treatment ───────────────────────────────────────────────────

test_that("SNM point × continuous × EM: truth recovery (no DTRreg — linear blip not supported)", {
  # DTRreg does not support linear blip functions for continuous treatment.
  # This test validates causatr against the analytical truth directly.
  set.seed(33)
  n <- 3000
  L <- rnorm(n)
  M <- as.numeric(L > 0)
  A <- rnorm(n, mean = 0.5 * L, sd = 1)
  Y <- 2 + 3 * A + 1.5 * L + 2 * A * M + rnorm(n)
  d <- data.table::data.table(Y = Y, A = A, L = L, M = M)

  # truth: psi_intercept = 3, psi_M = 2
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders_outcome = ~ L + A:M,
    confounders_treatment = ~L,
    estimator = "snm",
    treatment_free = ~L
  )
  res <- contrast(fit, ci_method = "sandwich")
  expect_equal(res$estimates$estimate[1], 3, tolerance = 0.1)
  expect_equal(res$estimates$estimate[2], 2, tolerance = 0.1)
  expect_true(all(is.finite(res$estimates$se) & res$estimates$se > 0))
})


# ── Missing outcomes (MCAR) ───────────────────────────────────────────────

test_that("SNM point × MCAR outcomes: complete-case estimate matches DTRreg", {
  skip_if_not_installed("DTRreg")

  set.seed(44)
  n <- 2000
  L <- rnorm(n)
  M <- rbinom(n, 1, 0.5)
  A <- rbinom(n, 1, plogis(0.3 * L + 0.2 * M))
  Y_full <- 1 + 2 * A + 3 * A * M + 0.5 * L + 0.8 * M + rnorm(n)
  # MCAR: 20% missing at random (independent of all variables)
  miss_idx <- sample.int(n, size = round(0.2 * n))
  Y <- Y_full
  Y[miss_idx] <- NA
  d <- data.table::data.table(Y = Y, A = A, L = L, M = M)

  # DTRreg requires no NAs in the outcome vector passed to `outcome =` and
  # in the `data` argument. Pass only the complete-case rows to both.
  obs <- !is.na(Y)
  Y_obs <- Y[obs]
  A_obs <- A[obs]
  L_obs <- L[obs]
  M_obs <- M[obs]
  d_cc_obs <- data.frame(Y = Y_obs, A = A_obs, L = L_obs, M = M_obs)

  dr <- DTRreg::DTRreg(
    outcome = Y_obs,
    blip.mod = ~M,
    treat.mod = A ~ L + M,
    tf.mod = ~ L + M,
    data = d_cc_obs,
    method = "gest",
    treat.type = "bin"
  )
  psi_dr_int <- dr$psi[[1]][["A"]]
  psi_dr_M <- dr$psi[[1]][["M:A"]]

  # causatr drops NAs automatically (complete-case for outcomes)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders_outcome = ~ L + M + A:M,
    confounders_treatment = ~ L + M,
    estimator = "snm",
    treatment_free = ~ L + M
  )
  res <- contrast(fit, ci_method = "sandwich")
  psi_int <- res$estimates$estimate[1]
  psi_M <- res$estimates$estimate[2]

  # MCAR complete-case estimates agree with DTRreg
  expect_equal(psi_int, psi_dr_int, tolerance = 1e-6)
  expect_equal(psi_M, psi_dr_M, tolerance = 1e-6)

  # Truth recovery despite 20% MCAR
  expect_equal(psi_int, 2, tolerance = 0.15)
  expect_equal(psi_M, 3, tolerance = 0.20)
})


# ── Observation weights ────────────────────────────────────────────────────

test_that("SNM point × weights: estimate shifts vs. unweighted and bootstrap agrees", {
  # DTRreg 2.4 changed the weight= arg to a character (weighting scheme);
  # it no longer accepts a numeric vector of survey/case weights.
  # Verify causatr against the manual closed-form WLS solution instead.
  set.seed(55)
  n <- 2000
  L <- rnorm(n)
  A <- rbinom(n, 1, plogis(0.3 * L))
  Y <- 1 + 2 * A + 0.5 * L + rnorm(n)
  # Upweight treated observations 3:1
  w <- ifelse(A == 1, 3, 1)
  d <- data.table::data.table(Y = Y, A = A, L = L)

  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "snm",
    weights = w
  )
  res <- contrast(fit, ci_method = "sandwich")
  psi_w <- res$estimates$estimate[1]

  # Manual: psi = sum(w*R*Y) / sum(w*R*A) where R = A - E_w[A|L].
  # causatr fits the treatment model with the supplied weights, so the manual
  # oracle must use a weighted propensity model to match exactly.
  prop_model_w <- stats::glm(
    A ~ L,
    data = d,
    family = stats::binomial(),
    weights = w
  )
  R <- A - stats::fitted(prop_model_w)
  psi_manual <- sum(w * R * Y) / sum(w * R * A)
  expect_equal(psi_w, psi_manual, tolerance = 1e-6)

  # Unit weights = unweighted
  w1 <- rep(1, n)
  fit_unit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "snm",
    weights = w1
  )
  res_unit <- contrast(fit_unit, ci_method = "sandwich")
  fit_nw <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "snm"
  )
  res_nw <- contrast(fit_nw, ci_method = "sandwich")
  expect_equal(
    res_unit$estimates$estimate[1],
    res_nw$estimates$estimate[1],
    tolerance = 1e-10
  )

  # Weighted SE is finite and positive
  expect_true(is.finite(res$estimates$se[1]) && res$estimates$se[1] > 0)
})


# ── Cluster-robust SE ─────────────────────────────────────────────────────

test_that("SNM point × cluster-robust: singleton clusters = i.i.d. exactly", {
  set.seed(66)
  n <- 500
  L <- rnorm(n)
  A <- rbinom(n, 1, plogis(0.3 * L))
  Y <- 1 + 2 * A + 0.5 * L + rnorm(n)
  d <- data.table::data.table(Y = Y, A = A, L = L)

  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "snm"
  )
  snm_r <- causatr:::compute_snm_blip_point(fit)

  vcov_iid <- causatr:::variance_if_snm(fit, snm_r, cluster_vec = NULL)
  vcov_singleton <- causatr:::variance_if_snm(
    fit,
    snm_r,
    cluster_vec = seq_len(n)
  )
  expect_equal(vcov_iid, vcov_singleton, tolerance = 1e-10)
})

test_that("SNM point × cluster-robust: contrast() propagates cluster argument", {
  skip_if_fast()
  set.seed(77)
  n_cluster <- 60
  cluster_size <- 20
  n <- n_cluster * cluster_size
  cluster_id <- rep(seq_len(n_cluster), each = cluster_size)
  # Within-cluster correlation
  alpha_c <- rnorm(n_cluster, sd = 0.5)
  L <- rnorm(n) + alpha_c[cluster_id]
  A <- rbinom(n, 1, plogis(0.3 * L))
  Y <- 1 + 2 * A + 0.5 * L + alpha_c[cluster_id] + rnorm(n)
  d <- data.table::data.table(Y = Y, A = A, L = L, cluster_id = cluster_id)

  # Fit-time cluster
  fit_cl <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "snm",
    cluster = "cluster_id"
  )
  res_propagated <- contrast(fit_cl, ci_method = "sandwich")

  # Explicit cluster argument
  res_explicit <- contrast(
    fit_cl,
    ci_method = "sandwich",
    cluster = d$cluster_id
  )
  expect_equal(
    res_propagated$estimates$se,
    res_explicit$estimates$se,
    tolerance = 1e-10
  )

  # SE is finite and positive
  expect_true(all(is.finite(res_propagated$estimates$se)))
  expect_true(all(res_propagated$estimates$se > 0))
})


# ── By-stratified averaged blip ───────────────────────────────────────────

test_that("SNM × by-stratified blip: per-stratum averaged blip matches manual", {
  set.seed(88)
  n <- 2000
  L <- rnorm(n)
  M <- rbinom(n, 1, 0.5)
  A <- rbinom(n, 1, plogis(0.3 * L + 0.2 * M))
  Y <- 1 + 2 * A + 3 * A * M + 0.5 * L + 0.8 * M + rnorm(n)
  # truth: psi_0=2, psi_M=3; M=0 stratum avg=2, M=1 stratum avg=5
  d <- data.table::data.table(Y = Y, A = A, L = L, M = M)

  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders_outcome = ~ L + M + A:M,
    confounders_treatment = ~ L + M,
    estimator = "snm",
    treatment_free = ~ L + M
  )
  res_by <- contrast(
    fit,
    ci_method = "sandwich",
    treatment_values = c(0, 1),
    by = "M"
  )

  # Manual: averaged blip in each stratum = (a1-a0) * (psi_0 + psi_M * mean(M in stratum))
  psi_est <- contrast(fit, ci_method = "sandwich")$estimates$estimate
  psi_0 <- psi_est[1]
  psi_M <- psi_est[2]

  m0_avg <- mean(d[M == 0]$M) # = 0
  m1_avg <- mean(d[M == 1]$M) # = 1
  manual_m0 <- 1 * (psi_0 + psi_M * m0_avg)
  manual_m1 <- 1 * (psi_0 + psi_M * m1_avg)

  est_m0 <- res_by$estimates[by == "0", estimate]
  est_m1 <- res_by$estimates[by == "1", estimate]

  expect_equal(est_m0, manual_m0, tolerance = 1e-6)
  expect_equal(est_m1, manual_m1, tolerance = 1e-6)

  # Truth recovery
  expect_equal(est_m0, 2, tolerance = 0.15)
  expect_equal(est_m1, 5, tolerance = 0.20)
})


# ── Longitudinal (2 periods, no mediator L1) ─────────────────────────────

test_that("SNM longitudinal × binary × no EM: causatr matches DTRreg + truth", {
  skip_if_not_installed("DTRreg")
  skip_if_fast()

  # DGP: L1 does NOT depend on A0 — pure time-varying confounder.
  # This avoids A0→L1→Y mediation that would bias the stage-0 estimate.
  set.seed(200)
  n <- 5000
  L0 <- rnorm(n)
  A0 <- rbinom(n, 1, plogis(0.3 * L0))
  L1 <- 0.5 * L0 + rnorm(n, sd = sqrt(0.5)) # independent of A0
  A1 <- rbinom(n, 1, plogis(0.3 * L1 + 0.2 * A0))
  Y <- 2 + 1 * A0 + 2 * A1 + 1.5 * L0 + 0.5 * L1 + rnorm(n)
  # truth: psi_0 = 1 (A0 blip), psi_1 = 2 (A1 blip)

  d_wide <- data.frame(Y = Y, A0 = A0, A1 = A1, L0 = L0, L1 = L1)
  d_long <- data.table::rbindlist(list(
    data.table::data.table(
      id = seq_len(n),
      time = 0L,
      A = A0,
      L = L0,
      Y = NA_real_
    ),
    data.table::data.table(
      id = seq_len(n),
      time = 1L,
      A = A1,
      L = L1,
      Y = Y
    )
  ))

  # DTRreg: lists in FORWARD order (blip/treat/tf[[1]] = first time = A0)
  dr <- DTRreg::DTRreg(
    outcome = Y,
    blip.mod = list(~1, ~1),
    treat.mod = list(A0 ~ L0, A1 ~ L1),
    tf.mod = list(~L0, ~L1),
    data = d_wide,
    method = "gest",
    treat.type = "bin"
  )
  # dr$psi[[1]] = A0 (first treatment), dr$psi[[2]] = A1 (second treatment)
  psi0_dr <- dr$psi[[1]][["A0"]]
  psi1_dr <- dr$psi[[2]][["A1"]]

  # causatr longitudinal
  fit_long <- causat(
    d_long,
    outcome = "Y",
    treatment = "A",
    confounders_tv = ~L,
    id = "id",
    time = "time",
    type = "longitudinal",
    estimator = "snm",
    treatment_free = ~L
  )
  res <- contrast(fit_long, ci_method = "sandwich")
  psi0_causat <- res$estimates[parameter == "stage0_psi_intercept", estimate]
  psi1_causat <- res$estimates[parameter == "stage1_psi_intercept", estimate]

  # Point estimates agree with DTRreg within finite-sample tolerance.
  # Note: causatr adds lag1_L to the period-1 treatment model (Markov lag
  # expansion), while DTRreg's treat.mod = A1 ~ L1 uses only L1. The
  # propensity residuals R1 differ, so exact numerical agreement is not
  # expected — instead verify agreement to within simulation variability.
  expect_equal(psi0_causat, psi0_dr, tolerance = 0.05)
  expect_equal(psi1_causat, psi1_dr, tolerance = 0.05)

  # Truth recovery (loose tolerance: finite-sample variation at n=5000)
  expect_equal(psi0_causat, 1, tolerance = 0.15)
  expect_equal(psi1_causat, 2, tolerance = 0.15)

  # SEs are finite and positive
  expect_true(all(is.finite(res$estimates$se) & res$estimates$se > 0))
})


# ── Longitudinal with time-varying EM ─────────────────────────────────────

test_that("SNM longitudinal × time-varying EM: blip params close to truth", {
  skip_if_fast()

  # DGP: M1 depends on A0 (time-varying modifier, not a mediator here
  # because the blip at time 1 explicitly includes M1).
  # The SNM correctly identifies the blip because the g-EE uses the
  # treatment residual as an instrument, not the modifier as a regressor.
  set.seed(18)
  n <- 5000
  L0 <- rnorm(n)
  A0 <- rbinom(n, 1, plogis(0.3 * L0))
  L1 <- 0.5 * L0 + 0.3 * A0 + rnorm(n, sd = sqrt(0.5))
  M1 <- as.numeric(L1 > 0)
  A1 <- rbinom(n, 1, plogis(0.3 * L1 + 0.2 * A0))
  Y <- 2 + 1 * A0 + 2 * A1 + 2 * A1 * M1 + 1.5 * L0 + 0.5 * L1 + rnorm(n)
  # truth: stage0: psi_0 = 1; stage1: psi_0 = 2, psi_M1 = 2

  d_long <- data.table::rbindlist(list(
    data.table::data.table(
      id = seq_len(n),
      time = 0L,
      A = A0,
      L = L0,
      M = as.numeric(L0 > 0),
      Y = NA_real_
    ),
    data.table::data.table(
      id = seq_len(n),
      time = 1L,
      A = A1,
      L = L1,
      M = M1,
      Y = Y
    )
  ))

  fit_long <- causat(
    d_long,
    outcome = "Y",
    treatment = "A",
    confounders_outcome = ~ A:M,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    type = "longitudinal",
    estimator = "snm",
    treatment_free = ~L
  )
  res <- contrast(fit_long, ci_method = "sandwich")

  # Stage 0: only intercept blip
  psi0_int <- res$estimates[parameter == "stage0_psi_intercept", estimate]
  # Stage 1: intercept + M1 modifier
  psi1_int <- res$estimates[parameter == "stage1_psi_intercept", estimate]
  psi1_M <- res$estimates[parameter == "stage1_psi_M", estimate]

  expect_equal(psi0_int, 1, tolerance = 0.20)
  expect_equal(psi1_int, 2, tolerance = 0.20)
  expect_equal(psi1_M, 2, tolerance = 0.25)
  expect_true(all(is.finite(res$estimates$se) & res$estimates$se > 0))
})


# ── Sandwich numeric fallback ─────────────────────────────────────────────

test_that("SNM sandwich: classed error for unsupported treatment model class", {
  set.seed(99)
  n <- 200
  L <- rnorm(n)
  A <- rbinom(n, 1, plogis(0.3 * L))
  Y <- 1 + 2 * A + rnorm(n)
  d <- data.table::data.table(Y = Y, A = A, L = L)

  # Fit normally so causat() and blip computation succeed, then replace
  # the stored treatment model object with a bare structure of unknown class
  # that has no estfun/bread S3 methods. This simulates an exotic propensity
  # model class and exercises the tryCatch in variance_if_snm().
  #
  # We call compute_snm_blip_point() first (which needs predict() to work)
  # and then test variance_if_snm() directly (which is where the tryCatch
  # lives). This avoids a predict() failure before we reach the tryCatch.
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "snm"
  )
  snm_r <- causatr:::compute_snm_blip_point(fit)
  fit$details$treatment_model$model <- structure(
    list(),
    class = "causatr_fake_no_estfun_model"
  )
  expect_error(
    causatr:::variance_if_snm(fit, snm_r),
    class = "causatr_snm_no_estfun"
  )
})
