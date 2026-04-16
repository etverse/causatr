# lmtp and manual IPSI contrast-level oracle tests for the self-contained
# IPW engine's non-static interventions.
#
# T-oracle5 compares causatr's shift() point estimate against
# lmtp::lmtp_sdr() on the same data with the same parametric learners.
# T-oracle6 compares causatr's ipsi() point estimate against a manual
# Kennedy (2019) closed-form weight computation.
#
# Only point estimates are compared — lmtp's SE uses the EIF, not the
# M-estimation sandwich, so SEs are structurally different objects and
# only asymptotically comparable.
#
# lmtp >= 1.5.0 removed `lmtp_ipw()` (deprecated; IPW requires correctly
# pre-specified parametric models). `lmtp_sdr()` (sequentially doubly
# robust) with `SL.glm` learners and `folds = 1` is consistent for the
# same target parameter under correct specification of both treatment
# and outcome models (guaranteed by the linear DGPs used here).
#
# See PHASE_4_INTERVENTIONS_SELF_IPW.md §9 for the oracle design rationale.

test_that("T-oracle5: causatr shift(1) ≈ lmtp::lmtp_sdr shift (continuous trt)", {
  skip_if_not_installed("lmtp")

  # DGP 3: L ~ N(0,1), A = 1 + 0.5*L + N(0,1), Y = 1 + 2*A + L + N(0,1).
  # Both causatr and lmtp fit parametric GLMs for treatment and outcome.
  # Under correct specification, both estimators target the same
  # causal parameter E[Y(A + delta)].
  d <- simulate_continuous_continuous(n = 3000, seed = 77)
  delta <- 1

  # causatr path: fit treatment density model, compute shift(1) weights,
  # refit Y ~ 1 Hájek MSM, predict-then-average.
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw"
  )
  res <- contrast(
    fit,
    interventions = list(shifted = shift(delta), nat = NULL),
    reference = "nat",
    ci_method = "sandwich"
  )
  causatr_theta <- res$estimates$estimate[
    res$estimates$intervention == "shifted"
  ]

  # lmtp path: lmtp_sdr with SL.glm (standard GLM, no ensemble) and
  # folds = 1 (no cross-fitting) to match the parametric setup.
  oracle <- lmtp_oracle_shift(
    data = d,
    outcome = "Y",
    treatment = "A",
    baseline = "L",
    delta = delta
  )

  # Both target E[Y(A+1)] under a correctly specified linear DGP.
  # SDR and density-ratio IPW are different estimators, so finite-
  # sample agreement is approximate. Tolerance 0.15 accommodates
  # the estimator difference while still catching gross miscalibration
  # (the structural truth is ~5.0, so 0.15 is a ~3% relative bound).
  expect_equal(causatr_theta, oracle$theta, tolerance = 0.15)
})


test_that("T-oracle6: causatr ipsi(2) ≈ manual Kennedy Hájek mean (binary trt)", {
  # DGP 1: L ~ N(0,1), A ~ Bern(expit(0.5*L)), Y = 2 + 3*A + 1.5*L + N(0,1).
  # Both causatr and the manual oracle fit a logistic PS model A ~ L.
  # Under the same propensity model the IPSI Hájek means must agree
  # algebraically — this is a pure arithmetic identity check.
  d <- simulate_binary_continuous(n = 3000, seed = 88)
  delta <- 2

  # causatr path: fit propensity model, compute Kennedy closed-form
  # IPSI weights, refit Y ~ 1 Hájek MSM, predict-then-average.
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw"
  )
  res <- contrast(
    fit,
    interventions = list(boosted = ipsi(delta), nat = NULL),
    reference = "nat",
    ci_method = "sandwich"
  )
  causatr_theta <- res$estimates$estimate[
    res$estimates$intervention == "boosted"
  ]

  # Manual oracle: same logistic PS model, same Kennedy weights, same
  # Hájek weighted mean. No causatr code involved — pure first-principles.
  oracle <- manual_ipsi_oracle(
    data = d,
    outcome = "Y",
    treatment = "A",
    ps_formula = ~L,
    delta = delta
  )

  # Sanity: IPSI(2) shifts the propensity odds upward, so the
  # counterfactual mean should differ from the naive sample mean.
  expect_false(isTRUE(all.equal(oracle$theta, mean(d$Y))))

  # Both compute the identical Kennedy closed-form weights under the
  # same logistic propensity model and the same Hájek weighted mean.
  # Agreement is algebraic — tolerance ~1e-6 reflects floating-point
  # rounding only.
  expect_equal(causatr_theta, oracle$theta, tolerance = 1e-6)
})
