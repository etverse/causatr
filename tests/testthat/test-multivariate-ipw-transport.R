# Multivariate IPW × transport (chunk 17j).
#
# DGP: see simulate_mv_transport() in helper-dgp.R.
#   Y = 1 + 2*A1 + 1.5*A2 + L + A1*L + eps
#   Target ATE [(1,1) vs (0,0)] = 3.5 + E_target[L]
#   Generalizability (target=all): 3.5 (exact, E[L]=0)
#   Transportability (target=S=0): 3.5 + mean(d$L[d$S == 0])
#
# The study ATE = 3.5 + E[L|S=1] > 3.5, so the naive estimator
# (fit on S=1 only) is upward-biased relative to the target ATE.

ivs_11_00 <- list(
  a1 = list(A1 = static(1), A2 = static(1)),
  a0 = list(A1 = static(0), A2 = static(0))
)

# Delta-method SE for the a1-a0 difference from a 2×2 vcov.
se_diff <- function(r) {
  v <- vcov(r)
  sqrt(v["a1", "a1"] + v["a0", "a0"] - 2 * v["a1", "a0"])
}

test_that("mv IPW transport: fit_rows restricted to S=1 study rows", {
  d <- simulate_mv_transport(n = 500, seed = 1)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = c("A1", "A2"),
    confounders = ~L,
    estimator = "ipw",
    model_fn = stats::glm,
    propensity_model_fn = stats::glm,
    target = "S",
    target_subset = "target"
  )
  expected <- which(d$S == 1L & !is.na(d$Y))
  actual <- which(fit$details$fit_rows)
  expect_equal(sort(actual), sort(expected))
  expect_true(isTRUE(fit$details$transport))
  expect_true(isTRUE(fit$details$is_multivariate))
})

test_that("mv IPW transport (transportability): recovers target ATE", {
  # n=8000 gives ESS ~750 per HT arm; SE ~ 0.06. Tolerance 0.15 ~ 2.5 SE.
  d <- simulate_mv_transport(n = 8000, seed = 42)
  truth_transport <- 3.5 + mean(d$L[d$S == 0])

  fit <- causat(
    d,
    outcome = "Y",
    treatment = c("A1", "A2"),
    confounders = ~L,
    estimator = "ipw",
    model_fn = stats::glm,
    propensity_model_fn = stats::glm,
    target = "S",
    target_subset = "target"
  )
  res <- contrast(
    fit,
    interventions = ivs_11_00,
    type = "difference",
    ci_method = "sandwich"
  )
  ate <- coef(res)["a1"] - coef(res)["a0"]
  expect_lt(abs(ate - truth_transport), 0.15)
})

test_that("mv IPW transport: corrects study-population bias", {
  # Naive study estimator (fit on S=1 only) gives study ATE = 3.5 + E[L|S=1].
  # Target ATE = 3.5 + E[L|S=0]. Since E[L|S=1] > 0 > E[L|S=0], they differ
  # by roughly 1.4, well above the 0.1 threshold.
  d <- simulate_mv_transport(n = 8000, seed = 42)
  truth_transport <- 3.5 + mean(d$L[d$S == 0])

  # Restrict to study rows to avoid NA-treatment rejection from causat().
  d_study <- d[d$S == 1, ]
  fit_naive <- causat(
    d_study,
    outcome = "Y",
    treatment = c("A1", "A2"),
    confounders = ~L,
    estimator = "ipw",
    model_fn = stats::glm,
    propensity_model_fn = stats::glm
  )
  res_naive <- contrast(
    fit_naive,
    interventions = ivs_11_00,
    type = "difference",
    ci_method = "sandwich"
  )
  ate_naive <- coef(res_naive)["a1"] - coef(res_naive)["a0"]
  expect_gt(abs(ate_naive - truth_transport), 0.1)
})

test_that("mv IPW transport (generalizability): recovers marginal ATE", {
  # target_subset="all": target is the whole population, E[L]=0, ATE=3.5.
  # HT weights are noisy so tolerance is 0.20.
  d <- simulate_mv_transport(n = 8000, seed = 77)

  fit <- causat(
    d,
    outcome = "Y",
    treatment = c("A1", "A2"),
    confounders = ~L,
    estimator = "ipw",
    model_fn = stats::glm,
    propensity_model_fn = stats::glm,
    target = "S",
    target_subset = "all"
  )
  res <- contrast(
    fit,
    interventions = ivs_11_00,
    type = "difference",
    ci_method = "sandwich"
  )
  ate <- coef(res)["a1"] - coef(res)["a0"]
  expect_lt(abs(ate - 3.5), 0.20)
})

test_that("mv IPW transport: sandwich SE plausible vs bootstrap", {
  skip_if_fast()
  # Verify that the gamma-block correction lands the sandwich SE within
  # a reasonable range of the bootstrap SE. A large discrepancy would
  # indicate the sampling-model correction is wrong.
  d <- simulate_mv_transport(n = 3000, seed = 7)

  fit <- causat(
    d,
    outcome = "Y",
    treatment = c("A1", "A2"),
    confounders = ~L,
    estimator = "ipw",
    model_fn = stats::glm,
    propensity_model_fn = stats::glm,
    target = "S",
    target_subset = "target"
  )
  r_sw <- contrast(
    fit,
    interventions = ivs_11_00,
    type = "difference",
    ci_method = "sandwich"
  )
  r_bt <- contrast(
    fit,
    interventions = ivs_11_00,
    type = "difference",
    ci_method = "bootstrap",
    n_boot = 99
  )

  se_sw <- se_diff(r_sw)
  se_bt <- se_diff(r_bt)
  ratio <- se_sw / se_bt
  expect_gt(ratio, 0.4)
  expect_lt(ratio, 2.5)
})

test_that("mv IPW transport: bootstrap point estimate near truth", {
  skip_if_fast()
  d <- simulate_mv_transport(n = 4000, seed = 5)
  truth_transport <- 3.5 + mean(d$L[d$S == 0])

  fit <- causat(
    d,
    outcome = "Y",
    treatment = c("A1", "A2"),
    confounders = ~L,
    estimator = "ipw",
    model_fn = stats::glm,
    propensity_model_fn = stats::glm,
    target = "S",
    target_subset = "target"
  )
  r_bt <- contrast(
    fit,
    interventions = ivs_11_00,
    type = "difference",
    ci_method = "bootstrap",
    n_boot = 79
  )
  ate_bt <- coef(r_bt)["a1"] - coef(r_bt)["a0"]
  expect_lt(abs(ate_bt - truth_transport), 0.25)
})

test_that("mv IPW transport: target=NULL leaves behaviour unchanged", {
  # Use study-only data so NA treatments do not trigger the censoring guard.
  d_study <- simulate_mv_transport(n = 1000, seed = 3)
  d_study <- d_study[d_study$S == 1, ]

  fit_no_tgt <- causat(
    d_study,
    outcome = "Y",
    treatment = c("A1", "A2"),
    confounders = ~L,
    estimator = "ipw",
    model_fn = stats::glm,
    propensity_model_fn = stats::glm
  )
  expect_false(isTRUE(fit_no_tgt$details$transport))
  expect_null(fit_no_tgt$target)

  res <- contrast(
    fit_no_tgt,
    interventions = ivs_11_00,
    type = "difference",
    ci_method = "sandwich"
  )
  ate <- coef(res)["a1"] - coef(res)["a0"]
  expect_true(is.finite(ate))
  expect_true(is.finite(se_diff(res)))
})

test_that("mv IPW transport: stabilize='marginal' + transport recovers target ATE", {
  # Stabilized weights give a different (lower-variance) weight decomposition
  # but the same estimand. The transport correction (gamma block) applies
  # identically — only other_w = pw / w_S_fit differs structurally.
  d <- simulate_mv_transport(n = 6000, seed = 4)
  truth_transport <- 3.5 + mean(d$L[d$S == 0])

  fit <- causat(
    d,
    outcome = "Y",
    treatment = c("A1", "A2"),
    confounders = ~L,
    estimator = "ipw",
    model_fn = stats::glm,
    propensity_model_fn = stats::glm,
    stabilize = "marginal",
    target = "S",
    target_subset = "target"
  )
  res <- contrast(
    fit,
    interventions = ivs_11_00,
    type = "difference",
    ci_method = "sandwich"
  )
  ate <- coef(res)["a1"] - coef(res)["a0"]
  se <- se_diff(res)
  expect_lt(abs(ate - truth_transport), 0.20)
  expect_true(is.finite(se))
  expect_gt(se, 0)
})
