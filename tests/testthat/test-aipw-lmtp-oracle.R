# AIPW longitudinal lmtp cross-checks (external SDR oracle).
# Split from test-aipw.R so these slow blocks run on their own parallel worker.

test_that("longitudinal AIPW: lmtp cross-check (binary static)", {
  skip_if_not_installed("lmtp")
  skip_if_not_installed("SuperLearner")

  d <- make_linear_scm(n = 1500, n_times = 2, seed = 60)
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
  ate_aipw <- res$contrasts$estimate[1]

  d_wide <- data.frame(
    L0 = d$L0[d$time == 0],
    A_0 = d$A[d$time == 0],
    L_1 = d$L[d$time == 1],
    A_1 = d$A[d$time == 1],
    Y = d$Y[d$time == 1]
  )
  # lmtp_sdr: always-treated
  lmtp_a1 <- lmtp::lmtp_sdr(
    data = d_wide,
    trt = c("A_0", "A_1"),
    outcome = "Y",
    baseline = "L0",
    time_vary = list(NULL, "L_1"),
    shift = function(data, trt) rep(1, nrow(data)),
    outcome_type = "continuous",
    learners_trt = "SL.glm",
    learners_outcome = "SL.glm",
    folds = 1
  )
  # lmtp_sdr: never-treated
  lmtp_a0 <- lmtp::lmtp_sdr(
    data = d_wide,
    trt = c("A_0", "A_1"),
    outcome = "Y",
    baseline = "L0",
    time_vary = list(NULL, "L_1"),
    shift = function(data, trt) rep(0, nrow(data)),
    outcome_type = "continuous",
    learners_trt = "SL.glm",
    learners_outcome = "SL.glm",
    folds = 1
  )
  ate_lmtp <- lmtp_a1$estimate@x - lmtp_a0$estimate@x
  expect_lt(abs(ate_aipw - ate_lmtp), 0.5)

  # lmtp sums two marginal EIF SEs in quadrature (ignoring covariance),
  # causatr sandwich targets the contrast directly — ratio < 1 expected.
  se_aipw <- res$contrasts$se[1]
  expect_true(is.finite(se_aipw), label = "causatr SE is finite")
  expect_true(
    is.finite(lmtp_a1$estimate@std_error),
    label = "lmtp_a1 SE is finite"
  )
  expect_true(
    is.finite(lmtp_a0$estimate@std_error),
    label = "lmtp_a0 SE is finite"
  )
  se_lmtp <- sqrt(
    lmtp_a1$estimate@std_error^2 + lmtp_a0$estimate@std_error^2
  )
  se_ratio <- se_aipw / se_lmtp
  expect_gt(se_ratio, 0.15)
  expect_lt(se_ratio, 5.0)
})

test_that("longitudinal AIPW: lmtp cross-check (continuous shift)", {
  skip_if_not_installed("lmtp")
  skip_if_not_installed("SuperLearner")

  d <- make_continuous_scm(n = 1500, seed = 61)
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
  est_aipw <- res$contrasts$estimate[1]

  d_wide <- data.frame(
    L0 = d$L0[d$time == 0],
    A_0 = d$A[d$time == 0],
    L_1 = d$L[d$time == 1],
    A_1 = d$A[d$time == 1],
    Y = d$Y[d$time == 1]
  )
  shift_fn <- function(data, trt) data[[trt]] + 0.5
  lmtp_up <- lmtp::lmtp_sdr(
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
  )
  lmtp_nat <- lmtp::lmtp_sdr(
    data = d_wide,
    trt = c("A_0", "A_1"),
    outcome = "Y",
    baseline = "L0",
    time_vary = list(NULL, "L_1"),
    shift = NULL,
    outcome_type = "continuous",
    learners_trt = "SL.glm",
    learners_outcome = "SL.glm",
    folds = 1
  )
  est_lmtp <- lmtp_up$estimate@x - lmtp_nat$estimate@x
  expect_lt(abs(est_aipw - est_lmtp), 0.5)

  # lmtp SE is sum-in-quadrature of two marginal EIF SEs (ignores covariance),
  # while causatr sandwich targets the contrast directly — expect a ratio well
  # below 1.  Smoke-test bounds only.
  se_aipw <- res$contrasts$se[1]
  expect_true(is.finite(se_aipw), label = "causatr SE is finite")
  expect_true(
    is.finite(lmtp_up$estimate@std_error),
    label = "lmtp_up SE is finite"
  )
  expect_true(
    is.finite(lmtp_nat$estimate@std_error),
    label = "lmtp_nat SE is finite"
  )
  se_lmtp <- sqrt(
    lmtp_up$estimate@std_error^2 + lmtp_nat$estimate@std_error^2
  )
  se_ratio <- se_aipw / se_lmtp
  expect_gt(se_ratio, 0.15)
  expect_lt(se_ratio, 5.0)
})

test_that("longitudinal AIPW: lmtp cross-check (binary outcome)", {
  skip_if_not_installed("lmtp")
  skip_if_not_installed("SuperLearner")

  set.seed(62)
  n <- 1500
  id <- rep(1:n, each = 2)
  time <- rep(0:1, times = n)
  L0 <- rnorm(n)[id]
  A <- rbinom(2 * n, 1, plogis(0.3 * L0))
  L <- rnorm(2 * n, mean = 0.3 * A)
  Y <- rep(NA_real_, 2 * n)
  Y[time == 1] <- as.numeric(rbinom(
    n,
    1,
    plogis(-1 + 0.5 * A[time == 1] + 0.3 * L[time == 1] + 0.2 * L0[time == 1])
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
    reference = "a0",
    ci_method = "sandwich"
  )
  est_aipw <- res$contrasts$estimate[1]

  d_wide <- data.frame(
    L0 = d$L0[d$time == 0],
    A_0 = d$A[d$time == 0],
    L_1 = d$L[d$time == 1],
    A_1 = d$A[d$time == 1],
    Y = d$Y[d$time == 1]
  )
  lmtp_a1 <- lmtp::lmtp_sdr(
    data = d_wide,
    trt = c("A_0", "A_1"),
    outcome = "Y",
    baseline = "L0",
    time_vary = list(NULL, "L_1"),
    shift = function(data, trt) rep(1, nrow(data)),
    outcome_type = "binomial",
    learners_trt = "SL.glm",
    learners_outcome = "SL.glm",
    folds = 1
  )
  lmtp_a0 <- lmtp::lmtp_sdr(
    data = d_wide,
    trt = c("A_0", "A_1"),
    outcome = "Y",
    baseline = "L0",
    time_vary = list(NULL, "L_1"),
    shift = function(data, trt) rep(0, nrow(data)),
    outcome_type = "binomial",
    learners_trt = "SL.glm",
    learners_outcome = "SL.glm",
    folds = 1
  )
  est_lmtp <- lmtp_a1$estimate@x - lmtp_a0$estimate@x
  expect_lt(abs(est_aipw - est_lmtp), 0.3)

  # Same caveat as continuous shift: lmtp sums marginal SEs in quadrature
  # (ignoring covariance), causatr targets the contrast — ratio < 1 expected.
  se_aipw <- res$contrasts$se[1]
  expect_true(is.finite(se_aipw), label = "causatr SE is finite")
  expect_true(
    is.finite(lmtp_a1$estimate@std_error),
    label = "lmtp_a1 SE is finite"
  )
  expect_true(
    is.finite(lmtp_a0$estimate@std_error),
    label = "lmtp_a0 SE is finite"
  )
  se_lmtp <- sqrt(
    lmtp_a1$estimate@std_error^2 + lmtp_a0$estimate@std_error^2
  )
  se_ratio <- se_aipw / se_lmtp
  expect_gt(se_ratio, 0.15)
  expect_lt(se_ratio, 5.0)
})
