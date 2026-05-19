# Cross-estimator triangulation: multiple estimators on the same data
# must agree when all models are correctly specified.

# -- Heterogeneous effect DGP (DGP 5): ATE agreement ---------------------

test_that("triangulation: gcomp, IPW, AIPW, matching agree on ATE (het DGP)", {
  d <- simulate_het_effect(n = 5000, seed = 42)

  ivs <- list(a1 = static(1), a0 = static(0))
  conf <- ~ L + sex + A:L + A:sex

  fit_gc <- causat(d, outcome = "Y", treatment = "A", confounders = conf)
  fit_ipw <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L + sex,
    estimator = "ipw",
    propensity_model_fn = stats::glm
  )
  fit_aipw <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = conf,
    estimator = "aipw",
    propensity_model_fn = stats::glm
  )
  fit_match <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L + sex,
    estimator = "matching"
  )

  res_gc <- contrast(fit_gc, interventions = ivs, reference = "a0")
  res_ipw <- contrast(fit_ipw, interventions = ivs, reference = "a0")
  res_aipw <- contrast(fit_aipw, interventions = ivs, reference = "a0")
  res_match <- contrast(fit_match, interventions = ivs, reference = "a0")

  ates <- c(
    gc = res_gc$contrasts$estimate[1],
    ipw = res_ipw$contrasts$estimate[1],
    aipw = res_aipw$contrasts$estimate[1],
    match = res_match$contrasts$estimate[1]
  )

  # True marginal ATE ≈ 3.6. All should be close.
  for (a in ates) expect_equal(a, 3.6, tolerance = 0.5)

  # Pairwise agreement within 0.5.
  pairs <- combn(ates, 2)
  for (j in seq_len(ncol(pairs))) {
    expect_lt(abs(pairs[1, j] - pairs[2, j]), 0.5)
  }
})


# -- Heterogeneous effect: ATT agreement ----------------------------------

test_that("triangulation: gcomp, IPW, matching agree on ATT (het DGP)", {
  d <- simulate_het_effect(n = 5000, seed = 42)

  ivs <- list(a1 = static(1), a0 = static(0))

  fit_gc <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L + sex + A:L + A:sex
  )
  fit_ipw <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L + sex,
    estimator = "ipw",
    propensity_model_fn = stats::glm,
    estimand = "ATT"
  )
  fit_match <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L + sex,
    estimator = "matching",
    estimand = "ATT"
  )

  res_gc <- contrast(
    fit_gc,
    interventions = ivs,
    reference = "a0",
    estimand = "ATT"
  )
  res_ipw <- contrast(fit_ipw, interventions = ivs, reference = "a0")
  res_match <- contrast(fit_match, interventions = ivs, reference = "a0")

  atts <- c(
    gc = res_gc$contrasts$estimate[1],
    ipw = res_ipw$contrasts$estimate[1],
    match = res_match$contrasts$estimate[1]
  )

  # True ATT ≈ 4.58 (mean over sex). Pairwise agreement within 1.0.
  # Matching ATT is noisier due to finite-sample bias from caliper effects.
  pairs <- combn(atts, 2)
  for (j in seq_len(ncol(pairs))) {
    expect_lt(abs(pairs[1, j] - pairs[2, j]), 1.0)
  }
})


# -- Binary outcome: RD agreement ----------------------------------------

test_that("triangulation: gcomp, IPW, AIPW agree on RD (binary outcome)", {
  d <- simulate_het_binary(n = 5000, seed = 42)

  ivs <- list(a1 = static(1), a0 = static(0))
  conf <- ~ L + sex + A:L + A:sex

  fit_gc <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = conf,
    family = "binomial"
  )
  fit_ipw <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L + sex,
    estimator = "ipw",
    propensity_model_fn = stats::glm
  )
  fit_aipw <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = conf,
    estimator = "aipw",
    propensity_model_fn = stats::glm,
    family = "binomial"
  )

  res_gc <- contrast(fit_gc, interventions = ivs, reference = "a0")
  res_ipw <- contrast(fit_ipw, interventions = ivs, reference = "a0")
  res_aipw <- contrast(fit_aipw, interventions = ivs, reference = "a0")

  rds <- c(
    gc = res_gc$contrasts$estimate[1],
    ipw = res_ipw$contrasts$estimate[1],
    aipw = res_aipw$contrasts$estimate[1]
  )

  # True marginal RD ≈ 0.33. All should be close.
  for (rd in rds) expect_equal(rd, 0.33, tolerance = 0.06)

  # Pairwise agreement within 0.05.
  pairs <- combn(rds, 2)
  for (j in seq_len(ncol(pairs))) {
    expect_lt(abs(pairs[1, j] - pairs[2, j]), 0.05)
  }
})


# -- Continuous treatment: shift agreement --------------------------------

test_that("triangulation: gcomp, IPW, AIPW agree on shift(-1) (cont trt)", {
  d <- simulate_continuous_continuous(n = 3000, seed = 42)

  ivs <- list(
    shifted = shift(-1),
    natural = shift(0)
  )

  fit_gc <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L
  )
  fit_ipw <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw",
    propensity_model_fn = stats::glm
  )
  fit_aipw <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "aipw",
    propensity_model_fn = stats::glm
  )

  res_gc <- contrast(fit_gc, interventions = ivs, reference = "natural")
  res_ipw <- contrast(fit_ipw, interventions = ivs, reference = "natural")
  res_aipw <- contrast(fit_aipw, interventions = ivs, reference = "natural")

  ests <- c(
    gc = res_gc$contrasts$estimate[1],
    ipw = res_ipw$contrasts$estimate[1],
    aipw = res_aipw$contrasts$estimate[1]
  )

  # True shift(-1) effect = -2. All should be close.
  for (e in ests) expect_equal(e, -2, tolerance = 0.3)

  # Pairwise agreement within 0.3.
  pairs <- combn(ests, 2)
  for (j in seq_len(ncol(pairs))) {
    expect_lt(abs(pairs[1, j] - pairs[2, j]), 0.3)
  }
})


# -- Categorical treatment: gcomp vs IPW agreement -----------------------

test_that("triangulation: gcomp and IPW agree on categorical static", {
  d <- simulate_categorical_continuous(n = 5000, seed = 42)

  ivs <- list(
    set_b = static("b"),
    set_a = static("a")
  )

  fit_gc <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L
  )
  fit_ipw <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw"
  )

  res_gc <- contrast(fit_gc, interventions = ivs, reference = "set_a")
  res_ipw <- contrast(fit_ipw, interventions = ivs, reference = "set_a")

  ate_gc <- res_gc$contrasts$estimate[1]
  ate_ipw <- res_ipw$contrasts$estimate[1]

  # True ATE("b" vs "a") = 3.
  expect_equal(ate_gc, 3, tolerance = 0.3)
  expect_equal(ate_ipw, 3, tolerance = 0.3)
  expect_lt(abs(ate_gc - ate_ipw), 0.3)
})
