# Monte-Carlo coverage study for causat_mice() pooling (Tier 2, skip_on_cran).
#
# Theory (Meng 1994; von Hippel 2020; Bartlett & Hughes 2020): Rubin's rules
# and Boot MI both deliver ~nominal-coverage intervals *provided the point
# estimator is consistent*. Under uncongeniality Rubin's variance "can be
# biased upwards or downwards depending on the specific situation" (Bartlett &
# Hughes 2020) -- it is conservative only for certain kinds of uncongeniality
# (Meng 1994), not in general. Crucially, a bootstrap "addresses variance
# estimation and coverage, not bias in the point estimate itself": if the
# imputation model is misspecified enough to make the causal estimator
# inconsistent, no pooling rule restores coverage.
#
# These DGPs (helper-dgp.R) make the regime explicit:
#  * simulate_mi_covariate()    -- DGP-MI1: linear outcome, Y available to the
#    imputer, so the gcomp estimator stays consistent. Both methods nominal.
#  * simulate_mi_uncongenial()  -- DGP-MI5: an A*L interaction whose imputation
#    requires P(L | A, Y); excluding Y from the predictor matrix misspecifies
#    the imputation and *biases* the causal estimate (~3.5 vs truth 3), so
#    coverage collapses for both methods. This is the cautionary case behind
#    the "include Y as a predictor" guidance, not a variance-conservatism demo.
# True marginal ATE = 3 throughout.

# Whether a pooled difference interval covers the (signed) truth. contrast()
# may report a1 - a0 or a0 - a1 depending on which intervention it picks as the
# reference, so the sign of the point estimate fixes which signed truth (+3 or
# -3) the interval should bracket.
mi_covers <- function(res, truth = 3) {
  est <- res$contrasts$estimate[1]
  signed <- sign(est) * truth
  res$contrasts$ci_lower[1] <= signed && signed <= res$contrasts$ci_upper[1]
}

# Muffle only the *expected* conditions of the misspecified-imputation arm
# (excluding Y always trips the missing-predictor warning; a modest bootstrap
# can floor a Boot MI variance component) so the loop stays quiet without
# blanket suppression. Neither condition is an assertion of this study.
mi_quiet <- function(expr) {
  withCallingHandlers(
    expr,
    causatr_mi_missing_predictors = function(w) invokeRestart("muffleWarning"),
    causatr_mi_boot_floor = function(w) invokeRestart("muffleWarning")
  )
}

test_that("Rubin pooling attains ~nominal coverage when the estimator is consistent", {
  skip_on_cran()
  skip_if_not_installed("mice")
  ints <- list(a1 = static(1), a0 = static(0))

  # DGP-MI1 has no A*L interaction and mice's default predictor matrix carries
  # the outcome, so PMM imputation of L from (A, Y) keeps the gcomp ATE
  # consistent. Rubin's rules then cover at the nominal level.
  R <- 80L
  cover <- logical(R)
  est <- numeric(R)
  for (r in seq_len(R)) {
    dat <- simulate_mi_covariate(n = 800, seed = 1000L + r)
    imp <- mice::mice(dat[, c("Y", "A", "L")], m = 10, printFlag = FALSE)
    res <- causat_mice(imp, "Y", "A", ~L, ints, estimator = "gcomp")
    est[r] <- abs(res$contrasts$estimate[1])
    cover[r] <- mi_covers(res)
  }

  # The strong, two-sided check: the pooled estimator is consistent, so the
  # mean estimate over reps pins the truth tightly (bias's MC SE ~ 0.012, so
  # |mean - 3| should be well under 0.045). This is what actually validates the
  # method; the coverage line below is a one-sided calibration floor.
  expect_equal(mean(est), 3, tolerance = 0.015)
  # Coverage floor: a Monte-Carlo coverage can only cheaply guard against gross
  # *under*-coverage (MC SE ~ sqrt(.95*.05/80) ~ 0.024 at R = 80; ruling out
  # mild *over*-coverage would need a far larger R). >= 0.90 is ~2 MC SE below
  # the 0.95 target; the observed coverage sits at ~0.96.
  expect_gte(mean(cover), 0.90)
})

test_that("Boot MI covers nominally and is no wider than Rubin in the consistent regime", {
  skip_on_cran()
  skip_if_not_installed("mice")
  ints <- list(a1 = static(1), a0 = static(0))

  # Same consistent DGP-MI1. Boot MI's resampling variance should also reach
  # nominal coverage, with intervals on par with (here marginally tighter than)
  # Rubin's -- not the dramatic over/under split the omit-Y misspecification
  # would (wrongly) suggest. B is kept modest for runtime; von Hippel's nominal
  # guarantee is asymptotic in B (he uses B ~ 1000), so this asserts the
  # *direction* of the SE comparison, not a tight coverage equality.
  R <- 25L
  rubin_cov <- logical(R)
  boot_cov <- logical(R)
  rubin_est <- numeric(R)
  boot_est <- numeric(R)
  rubin_se <- numeric(R)
  boot_se <- numeric(R)
  for (r in seq_len(R)) {
    dat <- simulate_mi_covariate(n = 700, seed = 5000L + r)
    imp <- mice::mice(dat[, c("Y", "A", "L")], m = 10, printFlag = FALSE)
    rr <- causat_mice(imp, "Y", "A", ~L, ints, estimator = "gcomp")
    bb <- mi_quiet(causat_mice(
      imp,
      "Y",
      "A",
      ~L,
      ints,
      estimator = "gcomp",
      pool_method = "boot_mi",
      B = 40,
      M = 2,
      seed = r
    ))
    rubin_cov[r] <- mi_covers(rr)
    boot_cov[r] <- mi_covers(bb)
    rubin_est[r] <- abs(rr$contrasts$estimate[1])
    boot_est[r] <- abs(bb$contrasts$estimate[1])
    rubin_se[r] <- rr$contrasts$se[1]
    boot_se[r] <- bb$contrasts$se[1]
  }

  # Strong two-sided checks: both pooled estimators are consistent here, so
  # their mean estimates pin the truth. (Boot MI's mean over reps is the
  # average of B*M fits, hence also ~unbiased.)
  expect_equal(mean(rubin_est), 3, tolerance = 0.02)
  expect_equal(mean(boot_est), 3, tolerance = 0.02)
  # Calibration floors (one-sided; see the note in the first test on why an MC
  # coverage study can only cheaply rule out under-coverage at this R).
  expect_gte(mean(rubin_cov), 0.90)
  expect_gte(mean(boot_cov), 0.85)
  # The two variances are on par in the consistent regime (here Rubin is
  # marginally larger, ~5%). Approximate parity, not a strict ordering: the gap
  # is small and a mice/B change could flip its sign, but the methods should
  # never disagree by more than a few percent.
  expect_gte(mean(rubin_se), 0.95 * mean(boot_se))
})

test_that("Omitting Y biases the causal estimate, collapsing coverage for both methods", {
  skip_on_cran()
  skip_if_not_installed("mice")
  ints <- list(a1 = static(1), a0 = static(0))

  # DGP-MI5 imputed two ways on the *same* dataset each replicate: Y included
  # (correctly specified imputation, consistent estimator) vs Y excluded
  # (misspecified, inconsistent estimator). Bartlett & Hughes' consistency
  # prerequisite predicts that excluding Y biases the point estimate, so its
  # coverage falls well below the Y-included arm -- a bias problem no variance
  # rule can repair. This is the test behind the package's "include Y as a
  # predictor" warning.
  R <- 30L
  est_incl <- numeric(R)
  est_excl <- numeric(R)
  cover_incl <- logical(R)
  cover_excl <- logical(R)
  for (r in seq_len(R)) {
    dat <- simulate_mi_uncongenial(n = 700, seed = 3000L + r)
    pred <- mice::make.predictorMatrix(dat[, c("Y", "A", "L")])
    pred[, "Y"] <- 0 # exclude Y -> misspecified imputation of L
    imp_excl <- mice::mice(
      dat[, c("Y", "A", "L")],
      m = 10,
      predictorMatrix = pred,
      printFlag = FALSE
    )
    imp_incl <- mice::mice(dat[, c("Y", "A", "L")], m = 10, printFlag = FALSE)
    rr_excl <- mi_quiet(causat_mice(
      imp_excl,
      "Y",
      "A",
      ~L,
      ints,
      estimator = "gcomp"
    ))
    rr_incl <- causat_mice(imp_incl, "Y", "A", ~L, ints, estimator = "gcomp")
    est_excl[r] <- abs(rr_excl$contrasts$estimate[1])
    est_incl[r] <- abs(rr_incl$contrasts$estimate[1])
    cover_excl[r] <- mi_covers(rr_excl)
    cover_incl[r] <- mi_covers(rr_incl)
  }

  # Strong two-sided pins on the estimates: with Y the estimator is consistent
  # (mean ~ 3); without Y it is biased clearly above the truth (~3.5). These
  # tie down the actual quantities, not just an ordering.
  expect_equal(mean(est_incl), 3, tolerance = 0.02)
  expect_gt(mean(est_excl), 3.2)
  # Consequence of the bias: the misspecified arm's coverage falls below the
  # consistent arm's. (Supporting directional checks, anchored by the pins.)
  expect_gt(abs(mean(est_excl) - 3), abs(mean(est_incl) - 3))
  expect_gt(mean(cover_incl), mean(cover_excl))
})
