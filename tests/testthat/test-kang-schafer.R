# Adversarial double-robustness tests based on Kang & Schafer (2007).
#
# Four scenarios testing AIPW's DR property under model misspecification.
# The DGP has Y independent of A, so ATE = 0 and E[Y] = 210.
# The analyst has access to either the true Z covariates or nonlinearly
# transformed X covariates.

# -- Scenario 1: both models correctly specified -------------------------

test_that("KS S1: all estimators recover E[Y]=210 with correct models", {
  ks <- simulate_kang_schafer(n = 5000, seed = 42)
  d <- ks$data
  truth <- ks$truth

  fit_gc <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~ Z1 + Z2 + Z3 + Z4
  )
  res_gc <- contrast(
    fit_gc,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0"
  )

  fit_ipw <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~ Z1 + Z2 + Z3 + Z4,
    estimator = "ipw",
    propensity_model_fn = stats::glm
  )
  res_ipw <- contrast(
    fit_ipw,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0"
  )

  fit_aipw <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~ Z1 + Z2 + Z3 + Z4,
    estimator = "aipw",
    propensity_model_fn = stats::glm
  )
  res_aipw <- contrast(
    fit_aipw,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0"
  )

  # E[Y(1)] and E[Y(0)] should both be ~210 (ATE ~0).
  ey1_gc <- res_gc$estimates$estimate[res_gc$estimates$intervention == "a1"]
  ey0_gc <- res_gc$estimates$estimate[res_gc$estimates$intervention == "a0"]
  expect_equal(ey1_gc, truth, tolerance = 0.01)
  expect_equal(ey0_gc, truth, tolerance = 0.01)

  ey1_ipw <- res_ipw$estimates$estimate[res_ipw$estimates$intervention == "a1"]
  ey0_ipw <- res_ipw$estimates$estimate[res_ipw$estimates$intervention == "a0"]
  expect_equal(ey1_ipw, truth, tolerance = 0.01)
  expect_equal(ey0_ipw, truth, tolerance = 0.01)

  ey1_aipw <- res_aipw$estimates$estimate[
    res_aipw$estimates$intervention == "a1"
  ]
  ey0_aipw <- res_aipw$estimates$estimate[
    res_aipw$estimates$intervention == "a0"
  ]
  expect_equal(ey1_aipw, truth, tolerance = 0.01)
  expect_equal(ey0_aipw, truth, tolerance = 0.01)
})


# -- Scenario 2: outcome misspecified, PS correct -------------------------

test_that("KS S2: AIPW + IPW recover 210 with correct PS, wrong outcome", {
  ks <- simulate_kang_schafer(n = 5000, seed = 42)
  d <- ks$data
  truth <- ks$truth

  # G-comp uses misspecified X covariates for the outcome model.
  fit_gc <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~ X1 + X2 + X3 + X4
  )
  res_gc <- contrast(
    fit_gc,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0"
  )
  ey1_gc <- res_gc$estimates$estimate[res_gc$estimates$intervention == "a1"]

  # G-comp should be biased under outcome misspecification.
  expect_gt(abs(ey1_gc - truth), 3)

  # IPW with correct PS (Z covariates) should recover truth.
  fit_ipw <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~ Z1 + Z2 + Z3 + Z4,
    estimator = "ipw",
    propensity_model_fn = stats::glm
  )
  res_ipw <- contrast(
    fit_ipw,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0"
  )
  ey1_ipw <- res_ipw$estimates$estimate[res_ipw$estimates$intervention == "a1"]
  expect_equal(ey1_ipw, truth, tolerance = 0.02)

  # AIPW with correct PS (Z) + wrong outcome (X).
  # The propensity model uses Z confounders; the outcome model also uses the
  # same confounders by default. To get misspecified outcome + correct PS,
  # we use a two-step approach: fit AIPW with correct PS confounders (Z),
  # which means both outcome and PS use Z. This is actually S1.
  # For S2, we need outcome model on X, PS on Z. causatr uses the same
  # confounders formula for both. We work around this by fitting AIPW with Z
  # confounders (correct PS) and checking robustness.
  #
  # Since causatr uses the same confounders for outcome and PS, the canonical
  # S2/S3 separation requires a trick: we use IPW (PS-only) and g-comp
  # (outcome-only) directly, and check AIPW with Z (both correct = S1).
  # The key DR test is:
  #   - IPW with correct Z recovers 210 (PS-only is fine).
  #   - G-comp with wrong X is biased (outcome-only fails).
  ey0_ipw <- res_ipw$estimates$estimate[res_ipw$estimates$intervention == "a0"]
  expect_equal(ey0_ipw, truth, tolerance = 0.02)
})


# -- Scenario 3: PS misspecified, outcome correct -------------------------

test_that("KS S3: g-comp recovers 210 with correct outcome, wrong PS", {
  ks <- simulate_kang_schafer(n = 5000, seed = 42)
  d <- ks$data
  truth <- ks$truth

  # G-comp with correct outcome model (Z covariates) — recovers truth.
  fit_gc <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~ Z1 + Z2 + Z3 + Z4
  )
  res_gc <- contrast(
    fit_gc,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0"
  )
  ey1_gc <- res_gc$estimates$estimate[res_gc$estimates$intervention == "a1"]
  expect_equal(ey1_gc, truth, tolerance = 0.01)

  # IPW with misspecified PS (X covariates) — should be biased.
  fit_ipw <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~ X1 + X2 + X3 + X4,
    estimator = "ipw",
    propensity_model_fn = stats::glm
  )
  res_ipw <- contrast(
    fit_ipw,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0"
  )
  ey1_ipw <- res_ipw$estimates$estimate[res_ipw$estimates$intervention == "a1"]

  # IPW should be meaningfully biased with misspecified PS.
  expect_gt(abs(ey1_ipw - truth), 2)
})


# -- Scenario 4: both models misspecified (negative control) --------------

test_that("KS S4: all estimators biased with both models wrong", {
  ks <- simulate_kang_schafer(n = 5000, seed = 42)
  d <- ks$data
  truth <- ks$truth

  # G-comp with misspecified outcome (X).
  fit_gc <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~ X1 + X2 + X3 + X4
  )
  res_gc <- contrast(
    fit_gc,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0"
  )
  ey1_gc <- res_gc$estimates$estimate[res_gc$estimates$intervention == "a1"]
  expect_gt(abs(ey1_gc - truth), 3)

  # IPW with misspecified PS (X).
  fit_ipw <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~ X1 + X2 + X3 + X4,
    estimator = "ipw",
    propensity_model_fn = stats::glm
  )
  res_ipw <- contrast(
    fit_ipw,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0"
  )
  ey1_ipw <- res_ipw$estimates$estimate[res_ipw$estimates$intervention == "a1"]
  expect_gt(abs(ey1_ipw - truth), 2)

  # AIPW with both wrong (X for both outcome and PS) — also biased.
  fit_aipw <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~ X1 + X2 + X3 + X4,
    estimator = "aipw",
    propensity_model_fn = stats::glm
  )
  res_aipw <- contrast(
    fit_aipw,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0"
  )
  ey1_aipw <- res_aipw$estimates$estimate[
    res_aipw$estimates$intervention == "a1"
  ]
  # AIPW with both wrong should also be biased, though possibly less so.
  expect_gt(abs(ey1_aipw - truth), 2)
})


# -- AIPW DR property: correct confounders preserve consistency -----------

test_that("KS: AIPW with correct confounders (Z) has smaller bias than misspecified", {
  ks <- simulate_kang_schafer(n = 5000, seed = 42)
  d <- ks$data
  truth <- ks$truth

  # AIPW with correct Z confounders (S1 — both correct).
  fit_correct <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~ Z1 + Z2 + Z3 + Z4,
    estimator = "aipw",
    propensity_model_fn = stats::glm
  )
  res_correct <- contrast(
    fit_correct,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0"
  )

  # AIPW with misspecified X confounders (S4 — both wrong).
  fit_wrong <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~ X1 + X2 + X3 + X4,
    estimator = "aipw",
    propensity_model_fn = stats::glm
  )
  res_wrong <- contrast(
    fit_wrong,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0"
  )

  ey1_correct <- res_correct$estimates$estimate[
    res_correct$estimates$intervention == "a1"
  ]
  ey1_wrong <- res_wrong$estimates$estimate[
    res_wrong$estimates$intervention == "a1"
  ]

  # Correct AIPW should be close to truth.
  expect_equal(ey1_correct, truth, tolerance = 0.01)

  # Misspecified AIPW should be further from truth than correct AIPW.
  expect_lt(abs(ey1_correct - truth), abs(ey1_wrong - truth))
})


# -- Split-confounder DR tests (per-component formulas) --------------------

test_that("KS S2 split: AIPW recovers 210 with wrong outcome (X), correct PS (Z)", {
  ks <- simulate_kang_schafer(n = 5000, seed = 42)
  d <- ks$data
  truth <- ks$truth

  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders_outcome = ~ X1 + X2 + X3 + X4,
    confounders_treatment = ~ Z1 + Z2 + Z3 + Z4,
    estimator = "aipw",
    propensity_model_fn = stats::glm
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0"
  )

  ey1 <- res$estimates$estimate[res$estimates$intervention == "a1"]
  ey0 <- res$estimates$estimate[res$estimates$intervention == "a0"]
  expect_equal(ey1, truth, tolerance = 0.01)
  expect_equal(ey0, truth, tolerance = 0.01)

  # Confirm the outcome model actually used X (misspecified).
  outcome_vars <- all.vars(stats::formula(fit$model))
  expect_true("X1" %in% outcome_vars)
  expect_false("Z1" %in% outcome_vars)

  # Confirm the PS model actually used Z (correct).
  ps_vars <- all.vars(fit$details$treatment_model$ps_formula)
  expect_true("Z1" %in% ps_vars)
  expect_false("X1" %in% ps_vars)
})


test_that("KS S3 split: AIPW recovers 210 with correct outcome (Z), wrong PS (X)", {
  ks <- simulate_kang_schafer(n = 5000, seed = 42)
  d <- ks$data
  truth <- ks$truth

  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders_outcome = ~ Z1 + Z2 + Z3 + Z4,
    confounders_treatment = ~ X1 + X2 + X3 + X4,
    estimator = "aipw",
    propensity_model_fn = stats::glm
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0"
  )

  ey1 <- res$estimates$estimate[res$estimates$intervention == "a1"]
  ey0 <- res$estimates$estimate[res$estimates$intervention == "a0"]
  expect_equal(ey1, truth, tolerance = 0.01)
  expect_equal(ey0, truth, tolerance = 0.01)

  # Confirm routing.
  outcome_vars <- all.vars(stats::formula(fit$model))
  expect_true("Z1" %in% outcome_vars)
  expect_false("X1" %in% outcome_vars)
  ps_vars <- all.vars(fit$details$treatment_model$ps_formula)
  expect_true("X1" %in% ps_vars)
  expect_false("Z1" %in% ps_vars)
})


test_that("KS S4 split: AIPW biased when both components use X", {
  ks <- simulate_kang_schafer(n = 5000, seed = 42)
  d <- ks$data
  truth <- ks$truth

  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders_outcome = ~ X1 + X2 + X3 + X4,
    confounders_treatment = ~ X1 + X2 + X3 + X4,
    estimator = "aipw",
    propensity_model_fn = stats::glm
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0"
  )

  ey1 <- res$estimates$estimate[res$estimates$intervention == "a1"]
  expect_gt(abs(ey1 - truth), 2)
})


test_that("KS split DR: S2 and S3 both closer to truth than S4", {
  ks <- simulate_kang_schafer(n = 5000, seed = 42)
  d <- ks$data
  truth <- ks$truth

  ivs <- list(a1 = static(1), a0 = static(0))

  # S2: wrong outcome, correct PS.
  fit_s2 <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders_outcome = ~ X1 + X2 + X3 + X4,
    confounders_treatment = ~ Z1 + Z2 + Z3 + Z4,
    estimator = "aipw",
    propensity_model_fn = stats::glm
  )
  ey1_s2 <- contrast(fit_s2, ivs, reference = "a0")$estimates$estimate[
    contrast(fit_s2, ivs, reference = "a0")$estimates$intervention == "a1"
  ]

  # S3: correct outcome, wrong PS.
  fit_s3 <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders_outcome = ~ Z1 + Z2 + Z3 + Z4,
    confounders_treatment = ~ X1 + X2 + X3 + X4,
    estimator = "aipw",
    propensity_model_fn = stats::glm
  )
  ey1_s3 <- contrast(fit_s3, ivs, reference = "a0")$estimates$estimate[
    contrast(fit_s3, ivs, reference = "a0")$estimates$intervention == "a1"
  ]

  # S4: both wrong.
  fit_s4 <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders_outcome = ~ X1 + X2 + X3 + X4,
    confounders_treatment = ~ X1 + X2 + X3 + X4,
    estimator = "aipw",
    propensity_model_fn = stats::glm
  )
  ey1_s4 <- contrast(fit_s4, ivs, reference = "a0")$estimates$estimate[
    contrast(fit_s4, ivs, reference = "a0")$estimates$intervention == "a1"
  ]

  expect_lt(abs(ey1_s2 - truth), abs(ey1_s4 - truth))
  expect_lt(abs(ey1_s3 - truth), abs(ey1_s4 - truth))
})
