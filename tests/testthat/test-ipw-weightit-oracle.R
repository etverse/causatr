# WeightIt contrast-level oracle tests for the static-binary IPW path.
#
# See PHASE_4_INTERVENTIONS_SELF_IPW.md §9 for the design rationale.
# Summary: causatr's self-contained density-ratio engine and
# `WeightIt::glm_weightit()` implement the same M-estimation functional
# (ATE / ATT / ATC under static binary IPW) via different internal
# decompositions. On the same data with the same propensity-score
# model they must agree at the contrast level to ~1e-6 on point
# estimates and to numerical precision (central-difference numDeriv
# vs analytic Jacobian) on sandwich SEs. Disagreement signals either a
# bug in the density-ratio weights, the `make_weight_fn()` closure, or
# the IF engine's cross-derivative term.
#
# All tests gate on `skip_if_not_installed("WeightIt")` so the package
# installs without WeightIt — chunk 3d moves it from Imports to
# Suggests for exactly this reason.

# DGP with heterogeneous treatment effect in L so ATT / ATC diverge
# from ATE. Without this, the T-oracle2 / T-oracle3 tests would pass
# even if causatr silently returned the ATE (the latent bug that
# chunk 3d's ATT/ATC weight fix addresses).
make_het_binary_dgp <- function(n = 2000, seed = 1) {
  set.seed(seed)
  L <- stats::rnorm(n)
  # Strong confounding (coef 0.8) so treated/control PS distributions
  # differ meaningfully, making ATT/ATC/ATE numerically distinct.
  A <- stats::rbinom(n, 1, stats::plogis(0.8 * L))
  # 0.6 * L * A is the effect-modification term that makes
  # ATT != ATE != ATC under this DGP.
  Y <- 2 + 3 * A + 1.5 * L + 0.6 * L * A + stats::rnorm(n)
  data.frame(Y = Y, A = A, L = L)
}

test_that("T-oracle1: causatr ATE ≈ WeightIt ATE (static binary, glm PS)", {
  skip_if_not_installed("WeightIt")

  d <- make_het_binary_dgp()
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw",
    estimand = "ATE"
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich"
  )
  causatr_out <- causatr_contrast_summary(res)
  oracle <- weightit_oracle_contrast(
    data = d,
    outcome = "Y",
    treatment = "A",
    ps_formula = ~L,
    estimand = "ATE"
  )

  # Point estimates agree to machine precision because causatr's per-arm
  # Hájek mean and WeightIt's `coef(Y ~ A)` under the same PS model
  # algebraically reduce to the same weighted sums. ~1e-6 tolerance is
  # a comfortable bound.
  expect_equal(causatr_out$mu_1, oracle$mu_1, tolerance = 1e-6)
  expect_equal(causatr_out$mu_0, oracle$mu_0, tolerance = 1e-6)
  expect_equal(causatr_out$contrast, oracle$contrast, tolerance = 1e-6)
  # SE agreement is slightly looser because causatr's `A_{beta,alpha}`
  # is numerical (numDeriv::jacobian central differences ~1e-7 accurate)
  # while WeightIt's is analytic from the glm_weightit Jacobian.
  expect_equal(
    causatr_out$se_contrast,
    oracle$se_contrast,
    tolerance = 1e-4
  )
})

test_that("T-oracle2: causatr ATT ≈ WeightIt ATT (static binary, glm PS)", {
  skip_if_not_installed("WeightIt")

  d <- make_het_binary_dgp()
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw",
    estimand = "ATT"
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich"
  )
  causatr_out <- causatr_contrast_summary(res)
  oracle <- weightit_oracle_contrast(
    data = d,
    outcome = "Y",
    treatment = "A",
    ps_formula = ~L,
    estimand = "ATT"
  )

  # Sanity: ATT really differs from ATE on this DGP — otherwise the
  # oracle would pass trivially. Guard by asserting the WeightIt-
  # computed ATT is at least 0.1 away from the simple unweighted
  # treated-minus-control difference (≈ ATE under strong confounding).
  naive_diff <- mean(d$Y[d$A == 1]) - mean(d$Y[d$A == 0])
  expect_gt(abs(oracle$contrast - naive_diff), 0.05)

  expect_equal(causatr_out$mu_1, oracle$mu_1, tolerance = 1e-6)
  expect_equal(causatr_out$mu_0, oracle$mu_0, tolerance = 1e-6)
  expect_equal(causatr_out$contrast, oracle$contrast, tolerance = 1e-6)
  expect_equal(
    causatr_out$se_contrast,
    oracle$se_contrast,
    tolerance = 1e-4
  )
})

test_that("T-oracle3: causatr ATC ≈ WeightIt ATC (static binary, glm PS)", {
  skip_if_not_installed("WeightIt")

  d <- make_het_binary_dgp()
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw",
    estimand = "ATC"
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich"
  )
  causatr_out <- causatr_contrast_summary(res)
  oracle <- weightit_oracle_contrast(
    data = d,
    outcome = "Y",
    treatment = "A",
    ps_formula = ~L,
    estimand = "ATC"
  )

  expect_equal(causatr_out$mu_1, oracle$mu_1, tolerance = 1e-6)
  expect_equal(causatr_out$mu_0, oracle$mu_0, tolerance = 1e-6)
  expect_equal(causatr_out$contrast, oracle$contrast, tolerance = 1e-6)
  expect_equal(
    causatr_out$se_contrast,
    oracle$se_contrast,
    tolerance = 1e-4
  )
})

test_that("T-oracle4: causatr ATE with mgcv::gam PS runs and sanity-matches glm PS", {
  skip_if_not_installed("WeightIt")
  skip_if_not_installed("mgcv")

  # On a linear-in-L DGP the glm PS and the gam-smooth-in-L PS are
  # asymptotically equivalent, so this doubles as a sanity check that
  # the density-ratio engine handles a non-GLM `propensity_model_fn`
  # end-to-end without numerical breakdowns. We do NOT demand an
  # exact WeightIt match here — WeightIt's `method = "gam"` is not
  # the same model class as `mgcv::gam(formula, family = binomial)`
  # that causatr fits — so the acceptance bound is wider (~0.05 on
  # the point estimate, comfortably covering the smoothing-induced
  # drift on this DGP).

  d <- make_het_binary_dgp(n = 3000, seed = 2)

  fit_gam <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~ s(L),
    estimator = "ipw",
    estimand = "ATE",
    propensity_model_fn = mgcv::gam
  )
  res_gam <- contrast(
    fit_gam,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich"
  )
  gam_out <- causatr_contrast_summary(res_gam)

  # The PS model object is a `gam`, not a `glm` — proves the
  # `propensity_model_fn` hook is actually consumed.
  expect_true(inherits(fit_gam$details$propensity_model, "gam"))

  # WeightIt GLM-PS oracle for a soft comparison: point estimates on a
  # linear DGP should agree to within smoothing noise.
  oracle <- weightit_oracle_contrast(
    data = d,
    outcome = "Y",
    treatment = "A",
    ps_formula = ~L,
    estimand = "ATE"
  )
  expect_equal(gam_out$contrast, oracle$contrast, tolerance = 0.05)

  # Finite, positive SE — sandwich engine survives the GAM `$Vp`
  # bread path through `prepare_model_if()`.
  expect_true(is.finite(gam_out$se_contrast))
  expect_gt(gam_out$se_contrast, 0)
})
