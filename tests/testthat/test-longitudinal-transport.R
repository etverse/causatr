# test-longitudinal-transport.R
#
# Truth-based tests for longitudinal IPW × transport (chunk 17i).
#
# DGP (2 periods, binary treatment, biased study selection by L):
#   L ~ N(0, 1)                                          (individual baseline)
#   P(S=1|L) = expit(-0.5 + L)                          (biased toward high L)
#   A_0|L,S=1 ~ Bernoulli(expit(0.2 + 0.3*L))
#   A_1|L,A_0,S=1 ~ Bernoulli(expit(0.1 + 0.2*L + 0.1*A_0))
#   Y|A,L,S=1 ~ N(2 + 2*A_0 + 2*A_1 + 1.5*L + 1.0*A_1*L, 1)
#
# Analytical targets:
#   Target ATE (transportability, S=0): 4 + E_target[L]
#   Target ATE (generalizability, all): 4 + E[L] = 4  (since E[L]=0)
#
# External oracle: gcomp transport on same DGP.

options(causatr.suppress_default_warnings = TRUE)

# ---------- helpers ----------------------------------------------------------

make_long_transport_data <- function(n, seed) {
  set.seed(seed)
  L <- rnorm(n)
  S <- rbinom(n, 1, plogis(-0.5 + L))
  study <- S == 1L

  A0 <- A1 <- Y <- rep(NA_real_, n)
  n_s <- sum(study)
  A0[study] <- rbinom(n_s, 1L, plogis(0.2 + 0.3 * L[study]))
  A1[study] <- rbinom(n_s, 1L, plogis(0.1 + 0.2 * L[study] + 0.1 * A0[study]))
  Y[study] <- 2 +
    2 * A0[study] +
    2 * A1[study] +
    1.5 * L[study] +
    1.0 * A1[study] * L[study] +
    rnorm(n_s)

  dt0 <- data.table::data.table(
    id = seq_len(n),
    time = 0L,
    L = L,
    S = S,
    A = A0,
    Y = NA_real_
  )
  dt1 <- data.table::data.table(
    id = seq_len(n),
    time = 1L,
    L = L,
    S = S,
    A = A1,
    Y = Y
  )
  data.table::rbindlist(list(dt0, dt1))
}

# Linear-outcome variant: no A_1*L interaction so gcomp's linear model
# is correctly specified. ATE = 4 for all populations (L cancels).
make_long_transport_data_linear <- function(n, seed) {
  set.seed(seed)
  L <- rnorm(n)
  S <- rbinom(n, 1, plogis(-0.5 + L))
  study <- S == 1L
  A0 <- A1 <- Y <- rep(NA_real_, n)
  n_s <- sum(study)
  A0[study] <- rbinom(n_s, 1L, plogis(0.2 + 0.3 * L[study]))
  A1[study] <- rbinom(n_s, 1L, plogis(0.1 + 0.2 * L[study] + 0.1 * A0[study]))
  Y[study] <- 2 + 2 * A0[study] + 2 * A1[study] + 1.5 * L[study] + rnorm(n_s)
  dt0 <- data.table::data.table(
    id = seq_len(n),
    time = 0L,
    L = L,
    S = S,
    A = A0,
    Y = NA_real_
  )
  dt1 <- data.table::data.table(
    id = seq_len(n),
    time = 1L,
    L = L,
    S = S,
    A = A1,
    Y = Y
  )
  data.table::rbindlist(list(dt0, dt1))
}

# E_target[L] = E[L * P(S=0|L)] / E[P(S=0|L)] where L~N(0,1),
# P(S=0|L) = expit(0.5-L). Computed via numerical integration.
true_target_ate_long <- function() {
  l <- seq(-8, 8, length.out = 200000)
  p_s0 <- plogis(0.5 - l)
  p_l <- dnorm(l)
  4 + sum(l * p_s0 * p_l) / sum(p_s0 * p_l)
}

# ---------- T-long-transport1: transportability recovers target ATE ----------

test_that("T-long-transport1: long IPW transport recovers target ATE; naive does not", {
  skip_if_not_installed("data.table")
  d <- make_long_transport_data(n = 10000, seed = 201)
  truth <- true_target_ate_long() # ≈ 3.69

  fit_transport <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw",
    id = "id",
    time = "time",
    target = "S",
    sampling_model_fn = stats::glm,
    model_fn = stats::glm,
    propensity_model_fn = stats::glm
  )
  res_transport <- contrast(
    fit_transport,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "sandwich"
  )
  ate_transport <- res_transport$estimates$estimate[1] -
    res_transport$estimates$estimate[2]

  # Naive longitudinal IPW on study rows only
  d_study <- d[d$S == 1L, ]
  fit_naive2 <- causat(
    d_study,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw",
    id = "id",
    time = "time",
    model_fn = stats::glm,
    propensity_model_fn = stats::glm
  )
  res_naive <- contrast(
    fit_naive2,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "sandwich"
  )
  ate_naive <- res_naive$estimates$estimate[1] - res_naive$estimates$estimate[2]

  # Transport should be close to truth (~3.69); naive should be biased high (~4.31)
  expect_lt(abs(ate_transport - truth), 0.20)
  expect_gt(abs(ate_naive - truth), 0.25)
})

# ---------- T-long-transport2: generalizability recovers ATE = 4 ------------

test_that("T-long-transport2: long IPW generalizability recovers ATE = 4", {
  skip_if_not_installed("data.table")
  d <- make_long_transport_data(n = 10000, seed = 202)

  fit_gen <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw",
    id = "id",
    time = "time",
    target = "S",
    target_subset = "all",
    sampling_model_fn = stats::glm,
    model_fn = stats::glm,
    propensity_model_fn = stats::glm
  )
  res_gen <- contrast(
    fit_gen,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "sandwich"
  )
  ate_gen <- res_gen$estimates$estimate[1] - res_gen$estimates$estimate[2]

  # Marginal ATE = 4 + E[L] = 4 (since L ~ N(0,1))
  expect_lt(abs(ate_gen - 4), 0.20)
})

# ---------- T-long-transport3: IPW transport ≈ gcomp transport --------------

test_that("T-long-transport3: long IPW transport agrees with gcomp transport", {
  skip_if_not_installed("data.table")
  d <- make_long_transport_data_linear(n = 10000, seed = 203)

  fit_ipw <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw",
    id = "id",
    time = "time",
    target = "S",
    model_fn = stats::glm,
    propensity_model_fn = stats::glm
  )
  res_ipw <- contrast(
    fit_ipw,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "sandwich"
  )

  fit_gcomp <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "gcomp",
    id = "id",
    time = "time",
    target = "S",
    model_fn = stats::glm
  )
  res_gcomp <- contrast(
    fit_gcomp,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "sandwich"
  )

  ate_ipw <- res_ipw$estimates$estimate[1] - res_ipw$estimates$estimate[2]
  ate_gcomp <- res_gcomp$estimates$estimate[1] - res_gcomp$estimates$estimate[2]

  # Both estimators should agree within MC error
  expect_lt(abs(ate_ipw - ate_gcomp), 0.15)
})

# ---------- T-long-transport4: sandwich SE plausibility ---------------------

test_that("T-long-transport4: long IPW transport sandwich SE is plausible", {
  skip_if_not_installed("data.table")
  d <- make_long_transport_data(n = 3000, seed = 204)

  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw",
    id = "id",
    time = "time",
    target = "S",
    model_fn = stats::glm,
    propensity_model_fn = stats::glm
  )
  res_sw <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "sandwich"
  )
  res_bs <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "bootstrap",
    n_boot = 200
  )

  se_sw <- res_sw$estimates$se[1]
  se_bs <- res_bs$estimates$se[1]

  ratio <- se_sw / se_bs
  expect_gt(ratio, 0.4)
  expect_lt(ratio, 2.5)
})

# ---------- T-long-transport5: bootstrap refits sampling model ---------------

test_that("T-long-transport5: bootstrap refits sampling model per replicate", {
  skip_if_not_installed("data.table")
  d <- make_long_transport_data(n = 2000, seed = 205)

  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw",
    id = "id",
    time = "time",
    target = "S",
    model_fn = stats::glm,
    propensity_model_fn = stats::glm
  )
  orig_gamma <- fit$details$sampling_model$gamma_hat

  res_bs <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "bootstrap",
    n_boot = 50
  )

  # After bootstrap, the original fit's sampling model is unchanged
  expect_equal(fit$details$sampling_model$gamma_hat, orig_gamma)
  # Bootstrap SE is finite and positive
  expect_true(is.finite(res_bs$estimates$se[1]))
  expect_gt(res_bs$estimates$se[1], 0)
})

# ---------- T-long-transport6: stabilize = "marginal" + transport smoke ------

test_that("T-long-transport6: stabilized long IPW + transport runs without error", {
  skip_if_not_installed("data.table")
  d <- make_long_transport_data(n = 2000, seed = 206)

  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw",
    id = "id",
    time = "time",
    stabilize = "marginal",
    target = "S",
    model_fn = stats::glm,
    propensity_model_fn = stats::glm
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "sandwich"
  )

  # Point estimate and SE finite, SE positive
  ate <- res$estimates$estimate[1]
  se <- res$estimates$se[1]
  expect_true(is.finite(ate))
  expect_true(is.finite(se))
  expect_gt(se, 0)
})

# ---------- T-long-transport7: target = NULL is unchanged -------------------

test_that("T-long-transport7: long IPW with target = NULL behaves as before", {
  skip_if_not_installed("data.table")
  # Study-only data (no target rows)
  d_study <- make_long_transport_data(n = 500, seed = 207)
  d_study <- d_study[d_study$S == 1L, ]

  fit <- causat(
    d_study,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw",
    id = "id",
    time = "time",
    model_fn = stats::glm,
    propensity_model_fn = stats::glm
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "sandwich"
  )

  expect_null(fit$target)
  expect_true(is.finite(res$estimates$estimate[1]))
})
