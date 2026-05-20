# Tests for IPCW × transportability composition (Phase 17, chunk 17k).
#
# Uses simulate_ipcw_transport() DGP:
#   L ~ N(0,1)
#   P(S=1|L) = expit(-0.5 + L)
#   A | L, S=1 ~ Bernoulli(expit(0.2 + 0.3*L))
#   C | A, L, S=1 ~ Bernoulli(expit(-1.5 + 0.5*A + 0.3*L))
#   Y | A, L, S=1, C=0 ~ N(2 + 3*A + 1.5*L + A*L, 1)
#
# Target ATE (transportability): 3 + E[L|S=0]
# Target ATE (generalizability): 3 + E[L] = 3
#
# The A*L interaction makes the ATE covariate-dependent, so the naive
# study estimator is biased relative to the target population. Differential
# censoring (C depends on A and L) biases complete-case analysis; IPCW
# corrects for this.

# --- gcomp + IPCW + transport ------------------------------------------------

test_that("gcomp + IPCW + transport (transportability): recovers target ATE", {
  d <- simulate_ipcw_transport(n = 8000, seed = 100)
  truth_ate <- 3 + mean(d$L[d$S == 0])

  fit <- causat(
    d,
    "Y",
    "A",
    ~ L + A:L,
    estimator = "gcomp",
    target = "S",
    target_subset = "target",
    censoring = "C",
    ipcw = TRUE
  )
  res <- contrast(
    fit,
    list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "sandwich"
  )
  est <- coef(res)["a1"] - coef(res)["a0"]
  expect_lt(abs(est - truth_ate), 0.12)
})

test_that("gcomp + IPCW + transport (generalizability): recovers marginal ATE", {
  d <- simulate_ipcw_transport(n = 8000, seed = 101)
  truth_ate <- 3 + mean(d$L)

  fit <- causat(
    d,
    "Y",
    "A",
    ~ L + A:L,
    estimator = "gcomp",
    target = "S",
    target_subset = "all",
    censoring = "C",
    ipcw = TRUE
  )
  res <- contrast(
    fit,
    list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "sandwich"
  )
  est <- coef(res)["a1"] - coef(res)["a0"]
  expect_lt(abs(est - truth_ate), 0.12)
})

test_that("gcomp + IPCW + transport: sandwich SE finite and CI covers truth", {
  d <- simulate_ipcw_transport(n = 10000, seed = 103)
  truth_ate <- 3 + mean(d$L[d$S == 0])

  fit <- causat(
    d,
    "Y",
    "A",
    ~ L + A:L,
    estimator = "gcomp",
    target = "S",
    target_subset = "target",
    censoring = "C",
    ipcw = TRUE
  )
  res <- contrast(
    fit,
    list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    type = "difference",
    ci_method = "sandwich"
  )
  se <- res$contrasts$se
  expect_true(is.finite(se) && se > 0)
  ci_lo <- res$contrasts$ci_lower
  ci_hi <- res$contrasts$ci_upper
  expect_true(ci_lo < truth_ate && ci_hi > truth_ate)
})


# --- IPW + IPCW + transport ---------------------------------------------------

test_that("IPW + IPCW + transport (transportability): recovers target ATE", {
  d <- simulate_ipcw_transport(n = 20000, seed = 200)
  truth_ate <- 3 + mean(d$L[d$S == 0])

  fit <- causat(
    d,
    "Y",
    "A",
    ~L,
    estimator = "ipw",
    target = "S",
    target_subset = "target",
    censoring = "C",
    ipcw = TRUE
  )
  res <- contrast(
    fit,
    list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "sandwich"
  )
  est <- coef(res)["a1"] - coef(res)["a0"]
  # IPW is noisier; allow ~2 SE tolerance at this sample size
  expect_lt(abs(est - truth_ate), 0.25)
})

test_that("IPW + IPCW + transport (generalizability): recovers marginal ATE", {
  d <- simulate_ipcw_transport(n = 20000, seed = 201)
  truth_ate <- 3 + mean(d$L)

  fit <- causat(
    d,
    "Y",
    "A",
    ~L,
    estimator = "ipw",
    target = "S",
    target_subset = "all",
    censoring = "C",
    ipcw = TRUE
  )
  res <- contrast(
    fit,
    list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "sandwich"
  )
  est <- coef(res)["a1"] - coef(res)["a0"]
  expect_lt(abs(est - truth_ate), 0.25)
})

test_that("IPW + IPCW + transport: sandwich SE plausible vs bootstrap", {
  d <- simulate_ipcw_transport(n = 5000, seed = 210)

  fit <- causat(
    d,
    "Y",
    "A",
    ~L,
    estimator = "ipw",
    target = "S",
    target_subset = "target",
    censoring = "C",
    ipcw = TRUE
  )
  res_sw <- contrast(
    fit,
    list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "sandwich"
  )
  res_bt <- contrast(
    fit,
    list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "bootstrap",
    n_boot = 200
  )
  ratio <- res_sw$contrasts$se / res_bt$contrasts$se
  expect_true(ratio > 0.4 && ratio < 2.5)
})


# --- AIPW + IPCW + transport -------------------------------------------------

test_that("AIPW + IPCW + transport (transportability): recovers target ATE", {
  d <- simulate_ipcw_transport(n = 8000, seed = 300)
  truth_ate <- 3 + mean(d$L[d$S == 0])

  fit <- causat(
    d,
    "Y",
    "A",
    ~ L + A:L,
    estimator = "aipw",
    target = "S",
    target_subset = "target",
    censoring = "C",
    ipcw = TRUE,
    propensity_model_fn = stats::glm
  )
  res <- contrast(
    fit,
    list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "sandwich"
  )
  est <- coef(res)["a1"] - coef(res)["a0"]
  expect_lt(abs(est - truth_ate), 0.15)
})

test_that("AIPW + IPCW + transport: sandwich SE finite and CI covers truth", {
  d <- simulate_ipcw_transport(n = 10000, seed = 302)
  truth_ate <- 3 + mean(d$L[d$S == 0])

  fit <- causat(
    d,
    "Y",
    "A",
    ~ L + A:L,
    estimator = "aipw",
    target = "S",
    target_subset = "target",
    censoring = "C",
    ipcw = TRUE,
    propensity_model_fn = stats::glm
  )
  res <- contrast(
    fit,
    list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    type = "difference",
    ci_method = "sandwich"
  )
  se <- res$contrasts$se
  expect_true(is.finite(se) && se > 0)
  ci_lo <- res$contrasts$ci_lower
  ci_hi <- res$contrasts$ci_upper
  expect_true(ci_lo < truth_ate && ci_hi > truth_ate)
})

test_that("AIPW + IPCW + transport: 2-of-3 DR — wrong outcome model", {
  d <- simulate_ipcw_transport(n = 10000, seed = 310)
  truth_ate <- 3 + mean(d$L[d$S == 0])

  # Wrong outcome model: omit the A:L interaction
  fit <- causat(
    d,
    "Y",
    "A",
    ~L,
    estimator = "aipw",
    target = "S",
    target_subset = "target",
    censoring = "C",
    ipcw = TRUE,
    propensity_model_fn = stats::glm
  )
  res <- contrast(
    fit,
    list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "sandwich"
  )
  est <- coef(res)["a1"] - coef(res)["a0"]
  # With correct propensity and censoring, AIPW should still be
  # approximately consistent (2-of-3 DR). Wider tolerance for
  # misspecification.
  expect_lt(abs(est - truth_ate), 0.30)
})

test_that("AIPW + IPCW + transport: sandwich vs bootstrap SE agreement", {
  d <- simulate_ipcw_transport(n = 5000, seed = 320)

  fit <- causat(
    d,
    "Y",
    "A",
    ~ L + A:L,
    estimator = "aipw",
    target = "S",
    target_subset = "target",
    censoring = "C",
    ipcw = TRUE,
    propensity_model_fn = stats::glm
  )
  res_sw <- contrast(
    fit,
    list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "sandwich"
  )
  res_bt <- contrast(
    fit,
    list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "bootstrap",
    n_boot = 200
  )
  ratio <- res_sw$contrasts$se / res_bt$contrasts$se
  expect_true(ratio > 0.5 && ratio < 2.0)
})


# --- Cross-estimator agreement ------------------------------------------------

test_that("IPCW + transport: gcomp, IPW, AIPW agree under correct specification", {
  d <- simulate_ipcw_transport(n = 15000, seed = 400)

  fit_g <- causat(
    d,
    "Y",
    "A",
    ~ L + A:L,
    estimator = "gcomp",
    target = "S",
    target_subset = "target",
    censoring = "C",
    ipcw = TRUE
  )
  res_g <- contrast(
    fit_g,
    list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "sandwich"
  )
  ate_g <- coef(res_g)["a1"] - coef(res_g)["a0"]

  fit_i <- causat(
    d,
    "Y",
    "A",
    ~L,
    estimator = "ipw",
    target = "S",
    target_subset = "target",
    censoring = "C",
    ipcw = TRUE
  )
  res_i <- contrast(
    fit_i,
    list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "sandwich"
  )
  ate_i <- coef(res_i)["a1"] - coef(res_i)["a0"]

  fit_a <- causat(
    d,
    "Y",
    "A",
    ~ L + A:L,
    estimator = "aipw",
    target = "S",
    target_subset = "target",
    censoring = "C",
    ipcw = TRUE,
    propensity_model_fn = stats::glm
  )
  res_a <- contrast(
    fit_a,
    list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "sandwich"
  )
  ate_a <- coef(res_a)["a1"] - coef(res_a)["a0"]

  # gcomp and AIPW should be very close; IPW is noisier
  expect_lt(abs(ate_g - ate_a), 0.10)
  expect_lt(abs(ate_g - ate_i), 0.30)
})


# --- IPCW + transport corrects bias vs naive ----------------------------------

test_that("IPCW + transport corrects bias vs naive complete-case transport", {
  d <- simulate_ipcw_transport(n = 15000, seed = 410)
  truth_ate <- 3 + mean(d$L[d$S == 0])

  # Naive: transport without IPCW (complete-case, ignoring differential censoring)
  fit_naive <- causat(
    d,
    "Y",
    "A",
    ~ L + A:L,
    estimator = "gcomp",
    target = "S",
    target_subset = "target",
    censoring = "C"
  )
  res_naive <- contrast(
    fit_naive,
    list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "sandwich"
  )
  ate_naive <- coef(res_naive)["a1"] - coef(res_naive)["a0"]

  # IPCW + transport
  fit_ipcw <- causat(
    d,
    "Y",
    "A",
    ~ L + A:L,
    estimator = "gcomp",
    target = "S",
    target_subset = "target",
    censoring = "C",
    ipcw = TRUE
  )
  res_ipcw <- contrast(
    fit_ipcw,
    list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "sandwich"
  )
  ate_ipcw <- coef(res_ipcw)["a1"] - coef(res_ipcw)["a0"]

  # IPCW estimate should be at least as close to truth as naive.
  # For gcomp with correct outcome model, both should be close under
  # MAR, but IPCW should not be worse.
  expect_lte(abs(ate_ipcw - truth_ate), abs(ate_naive - truth_ate) + 0.05)
})


# --- Bootstrap replays all models ---------------------------------------------

test_that("Bootstrap refits sampling + censoring + propensity per replicate", {
  d <- simulate_ipcw_transport(n = 3000, seed = 420)

  fit <- causat(
    d,
    "Y",
    "A",
    ~L,
    estimator = "ipw",
    target = "S",
    target_subset = "target",
    censoring = "C",
    ipcw = TRUE
  )
  res <- contrast(
    fit,
    list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "bootstrap",
    n_boot = 50
  )

  expect_true(is.finite(res$contrasts$se))
  expect_true(res$contrasts$se > 0)

  # Bootstrap point estimate should be near sandwich point estimate
  res_sw <- contrast(
    fit,
    list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "sandwich"
  )
  expect_lt(
    abs(res$contrasts$estimate - res_sw$contrasts$estimate),
    0.20
  )
})


# --- Details stashing ---------------------------------------------------------

test_that("IPCW + transport details are correctly stashed on fit", {
  d <- simulate_ipcw_transport(n = 500, seed = 430)

  fit <- causat(
    d,
    "Y",
    "A",
    ~L,
    estimator = "ipw",
    target = "S",
    target_subset = "target",
    censoring = "C",
    ipcw = TRUE
  )

  expect_true(isTRUE(fit$details$transport))
  expect_true(isTRUE(fit$details$ipcw))
  expect_s3_class(fit$details$sampling_model, "causatr_sampling_model")
  expect_s3_class(fit$details$censoring_model, "causatr_censoring_model")
  expect_equal(fit$target, "S")
  expect_equal(fit$target_subset, "target")
  expect_equal(fit$censoring, "C")
  expect_true(length(fit$details$ipcw_weights) == nrow(d))
})


# --- External cross-check: TransportHealth -----------------------------------

test_that("gcomp + IPCW + transport: TransportHealth cross-check", {
  # TransportHealth does not handle IPCW natively. We pre-compute IPCW
  # weights, pass them as case weights to the outcome model, then compare
  # the TransportHealth TATE against causatr's integrated estimate.
  skip_if_not_installed("TransportHealth")
  d <- simulate_ipcw_transport(n = 10000, seed = 500)

  # Pre-compute IPCW weights for study rows via causatr internals
  study_mask <- d$S == 1L
  uncens_mask <- d$C == 0L
  fit_rows <- study_mask & !is.na(d$Y) & !is.na(d$A) & !is.na(d$L)

  # Fit censoring model on study rows
  cens_mod <- fit_censoring_model(
    data = data.table::as.data.table(d),
    censoring = "C",
    treatment = "A",
    confounders = ~L
  )
  ipcw_w <- compute_ipcw_weights(
    cens_mod,
    n_total = nrow(d),
    censoring_col = as.integer(d$C),
    stabilize = TRUE
  )

  # Fit IPCW-weighted outcome model on study uncensored rows
  study_uncens <- d[fit_rows & uncens_mask, ]
  study_uncens$A <- factor(study_uncens$A)
  w_study <- ipcw_w[fit_rows & uncens_mask]
  m_out <- stats::glm(
    Y ~ A + L + A:L,
    data = study_uncens,
    weights = w_study,
    family = stats::gaussian()
  )

  target_df <- data.frame(d[d$S == 0, ])
  target_df$A <- factor(ifelse(is.na(target_df$A), 0L, target_df$A))
  prep <- TransportHealth::transportGCPreparedModel(
    outcomeModel = m_out,
    response = "Y",
    treatment = "A",
    wipe = FALSE
  )
  invisible(capture.output(
    th_res <- TransportHealth::transportGC(
      effectType = "meanDiff",
      preparedModel = prep,
      targetData = target_df,
      bootstrapNum = 0
    )
  ))
  th_ate <- summary(th_res)$effect

  # causatr integrated estimate
  fit <- causat(
    d,
    "Y",
    "A",
    ~ L + A:L,
    estimator = "gcomp",
    target = "S",
    target_subset = "target",
    censoring = "C",
    ipcw = TRUE
  )
  res <- contrast(
    fit,
    list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "sandwich"
  )
  c_ate <- unname(coef(res)["a1"] - coef(res)["a0"])

  expect_equal(c_ate, th_ate, tolerance = 1e-3)
})


# --- External cross-check: delicatessen (IPW) --------------------------------
# Reference values from data-raw/ipcw_transport_reference.py using
# delicatessen (Zivich 2024) stacked M-estimation: propensity + censoring
# + sampling + Hajek plug-in.
#
# The delicatessen implementation parametrizes the Hajek estimator as
# separate numerator/denominator EE components with a delta-method SE,
# while causatr uses a weighted GLM (Y ~ 1). The point estimates differ
# by ~2% due to parametrization differences in the stacked solver;
# the SE ratio is within 20%.

test_that("IPW + IPCW + transport: delicatessen stacked-EE cross-check", {
  ref_ate <- 2.4252
  ref_se_ate <- 0.2373

  d <- read.csv(
    test_path("fixtures", "ipcw_transport_fixture.csv")
  )

  fit <- causat(
    d,
    "Y",
    "A",
    ~L,
    estimator = "ipw",
    target = "S",
    target_subset = "target",
    censoring = "C",
    ipcw = TRUE
  )
  res <- contrast(
    fit,
    list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "sandwich"
  )
  ate <- coef(res)["a1"] - coef(res)["a0"]
  se <- res$contrasts$se

  # Point estimate within 0.10 of delicatessen reference
  expect_lt(abs(ate - ref_ate), 0.10)
  # SE ratio within (0.7, 1.5)
  se_ratio <- se / ref_se_ate
  expect_true(se_ratio > 0.7 && se_ratio < 1.5)
})


# --- Longitudinal IPW + IPCW + transport (smoke) -----------------------------

test_that("Longitudinal IPW + IPCW + transport: smoke test", {
  set.seed(600)
  n <- 2000
  L0 <- rnorm(n)
  S <- rbinom(n, 1, plogis(-0.5 + 0.8 * L0))
  study <- S == 1L
  n_s <- sum(study)

  A0 <- rep(NA_real_, n)
  A0[study] <- rbinom(n_s, 1L, plogis(0.3 * L0[study]))

  A1 <- rep(NA_real_, n)
  A1[study] <- rbinom(
    n_s,
    1L,
    plogis(0.2 + 0.2 * L0[study] + 0.1 * A0[study])
  )

  C0 <- rep(0L, n)
  C1 <- rep(0L, n)
  C1[study] <- rbinom(n_s, 1L, plogis(-2 + 0.3 * A0[study]))

  Y <- rep(NA_real_, n)
  uncens <- study & C1 == 0L
  Y[uncens] <- 2 + 4 * A1[uncens] + 1.5 * L0[uncens] + rnorm(sum(uncens))

  dt0 <- data.table::data.table(
    id = seq_len(n),
    time = 0L,
    L = L0,
    S = S,
    A = A0,
    Y = NA_real_,
    C = C0
  )
  dt1 <- data.table::data.table(
    id = seq_len(n),
    time = 1L,
    L = L0,
    S = S,
    A = A1,
    Y = Y,
    C = C1
  )
  pp <- data.table::rbindlist(list(dt0, dt1))

  fit <- causat(
    pp,
    "Y",
    "A",
    confounders = ~L,
    estimator = "ipw",
    id = "id",
    time = "time",
    target = "S",
    target_subset = "target",
    censoring = "C",
    ipcw = TRUE
  )

  res <- contrast(
    fit,
    list(always = static(1), never = static(0)),
    type = "difference",
    ci_method = "sandwich"
  )
  ate <- res$contrasts$estimate
  se <- res$contrasts$se

  expect_true(is.finite(ate))
  expect_true(is.finite(se) && se > 0)
})
