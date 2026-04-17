# Tests for the effect-modification infrastructure (Phase 6, chunk 6a).
#
# This file tests:
#   - parse_effect_mod() — formula parsing and classification
#   - build_ipw_msm_formula() — IPW MSM construction
#   - build_matching_msm_formula() — matching MSM construction
#   - check_em_compat() — bare-treatment rejection
#   - Integration with fit_ipw() and fit_matching() rejection paths

# parse_effect_mod() -------------------------------------------------------

test_that("parse_effect_mod detects A:sex as a single EM term", {
  em <- parse_effect_mod(~ L + sex + A:sex, "A")
  expect_true(em$has_em)
  expect_length(em$em_terms, 1L)
  expect_equal(em$em_terms[[1]]$treatment_var, "A")
  expect_equal(em$em_terms[[1]]$modifier_vars, "sex")
  expect_equal(em$confounder_terms, c("L", "sex"))
  expect_equal(em$modifier_vars, "sex")
})

test_that("parse_effect_mod handles reversed order sex:A", {
  em <- parse_effect_mod(~ L + sex + sex:A, "A")
  expect_true(em$has_em)
  expect_length(em$em_terms, 1L)
  expect_equal(em$em_terms[[1]]$treatment_var, "A")
  expect_equal(em$em_terms[[1]]$modifier_vars, "sex")
})

test_that("parse_effect_mod handles multiple EM terms", {
  em <- parse_effect_mod(~ L + sex + race + A:sex + A:race, "A")
  expect_true(em$has_em)
  expect_length(em$em_terms, 2L)
  expect_equal(em$modifier_vars, c("sex", "race"))
  expect_equal(em$confounder_terms, c("L", "sex", "race"))
})

test_that("parse_effect_mod returns has_em = FALSE when no EM terms", {
  em <- parse_effect_mod(~ L + sex + L:sex, "A")
  expect_false(em$has_em)
  expect_length(em$em_terms, 0L)
  expect_equal(em$confounder_terms, c("L", "sex", "L:sex"))
  expect_equal(em$modifier_vars, character(0L))
})

test_that("parse_effect_mod handles I() wrapped modifier", {
  em <- parse_effect_mod(~ L + A:I(age > 65), "A")
  expect_true(em$has_em)
  expect_length(em$em_terms, 1L)
  expect_equal(em$em_terms[[1]]$modifier_vars, "age")
})

test_that("parse_effect_mod detects bare treatment as EM term with no modifier", {
  # `~ L + A` puts A in confounders — detected as an EM term with

  # empty modifier_vars (downstream check_em_compat rejects it).
  em <- parse_effect_mod(~ L + A, "A")
  expect_true(em$has_em)
  expect_length(em$em_terms, 1L)
  expect_equal(em$em_terms[[1]]$modifier_vars, character(0L))
})

test_that("parse_effect_mod handles star expansion (A * sex)", {
  # `terms(~ A * sex)` expands to main effects + interaction:
  # term.labels = c("A", "sex", "A:sex")
  em <- parse_effect_mod(~ L + A * sex, "A")
  expect_true(em$has_em)
  # "A" (bare) and "A:sex" (interaction) both involve treatment
  bare <- Filter(function(x) length(x$modifier_vars) == 0L, em$em_terms)
  interactions <- Filter(function(x) length(x$modifier_vars) > 0L, em$em_terms)
  expect_length(bare, 1L)
  expect_length(interactions, 1L)
})

test_that("parse_effect_mod handles multivariate treatment", {
  em <- parse_effect_mod(~ L + sex + A1:sex + A2:sex, c("A1", "A2"))
  expect_true(em$has_em)
  expect_length(em$em_terms, 2L)
  expect_equal(em$em_terms[[1]]$treatment_var, "A1")
  expect_equal(em$em_terms[[2]]$treatment_var, "A2")
  expect_equal(em$confounder_terms, c("L", "sex"))
})

test_that("parse_effect_mod handles empty formula", {
  em <- parse_effect_mod(~1, "A")
  expect_false(em$has_em)
  expect_length(em$em_terms, 0L)
  expect_equal(em$confounder_terms, character(0L))
})

test_that("parse_effect_mod returns causatr_em_info class", {
  em <- parse_effect_mod(~ L + A:sex, "A")
  expect_s3_class(em, "causatr_em_info")
})

test_that("parse_effect_mod handles three-way interaction", {
  em <- parse_effect_mod(~ L + A:sex:race, "A")
  expect_true(em$has_em)
  expect_length(em$em_terms, 1L)
  expect_setequal(em$em_terms[[1]]$modifier_vars, c("sex", "race"))
})

# build_ipw_msm_formula() --------------------------------------------------

test_that("build_ipw_msm_formula returns Y ~ 1 without EM", {
  em <- parse_effect_mod(~ L + sex, "A")
  f <- build_ipw_msm_formula("Y", em)
  expect_equal(deparse(f), "Y ~ 1")
})

test_that("build_ipw_msm_formula includes modifier main effects with EM", {
  em <- parse_effect_mod(~ L + sex + A:sex, "A")
  f <- build_ipw_msm_formula("Y", em)
  # Should produce Y ~ 1 + sex (modifier main effect, no treatment)
  f_terms <- attr(terms(f), "term.labels")
  expect_true("sex" %in% f_terms)
  expect_false("A" %in% f_terms)
  expect_false("A:sex" %in% f_terms)
})

test_that("build_ipw_msm_formula handles multiple modifiers", {
  em <- parse_effect_mod(~ L + sex + race + A:sex + A:race, "A")
  f <- build_ipw_msm_formula("Y", em)
  f_terms <- attr(terms(f), "term.labels")
  expect_true("sex" %in% f_terms)
  expect_true("race" %in% f_terms)
  expect_false("A" %in% f_terms)
})

# build_matching_msm_formula() ----------------------------------------------

test_that("build_matching_msm_formula returns Y ~ A without EM", {
  em <- parse_effect_mod(~ L + sex, "A")
  f <- build_matching_msm_formula("Y", "A", em)
  f_terms <- attr(terms(f), "term.labels")
  expect_equal(f_terms, "A")
})

test_that("build_matching_msm_formula includes A + modifier + A:modifier with EM", {
  em <- parse_effect_mod(~ L + sex + A:sex, "A")
  f <- build_matching_msm_formula("Y", "A", em)
  f_terms <- attr(terms(f), "term.labels")
  expect_true("A" %in% f_terms)
  expect_true("sex" %in% f_terms)
  # R may canonicalize the interaction order
  expect_true(any(grepl("A.*sex|sex.*A", f_terms)))
})

# check_em_compat() ---------------------------------------------------------

test_that("check_em_compat rejects bare treatment in confounders", {
  em <- parse_effect_mod(~ L + A, "A")
  expect_error(
    check_em_compat(em, "A", "ipw"),
    class = "causatr_bare_treatment_in_confounders"
  )
})

test_that("check_em_compat allows true EM terms (A:sex)", {
  em <- parse_effect_mod(~ L + sex + A:sex, "A")
  expect_no_error(check_em_compat(em, "A", "ipw"))
})

test_that("check_em_compat allows formulas with no EM terms", {
  em <- parse_effect_mod(~ L + sex, "A")
  expect_no_error(check_em_compat(em, "A", "ipw"))
})

# build_ps_formula integration with EM terms --------------------------------

test_that("build_ps_formula strips EM terms from propensity RHS", {
  f <- build_ps_formula(~ L + sex + A:sex, "A")
  f_terms <- attr(terms(f), "term.labels")
  expect_true("L" %in% f_terms)
  expect_true("sex" %in% f_terms)
  expect_false(any(grepl("A", f_terms)))
  # Response should be the treatment
  expect_equal(as.character(f[[2]]), "A")
})

test_that("build_ps_formula works without EM terms (unchanged behavior)", {
  f <- build_ps_formula(~ L + sex + I(age^2), "A")
  f_terms <- attr(terms(f), "term.labels")
  expect_true("L" %in% f_terms)
  expect_true("sex" %in% f_terms)
  expect_true("I(age^2)" %in% f_terms)
})

test_that("build_ps_formula handles EM-only confounders", {
  # Edge case: confounders = ~ A:sex (only EM terms, no pure confounders).
  # After stripping EM terms, propensity formula should be intercept-only.
  f <- build_ps_formula(~ A:sex, "A")
  # Should get A ~ 1 (intercept-only)
  f_terms <- attr(terms(f), "term.labels")
  expect_length(f_terms, 0L)
})

# IPW effect modification (chunk 6b) ----------------------------------------

test_that("IPW rejects bare treatment in confounders", {
  d <- simulate_binary_continuous(n = 200, seed = 2)
  expect_error(
    causat(
      d,
      outcome = "Y",
      treatment = "A",
      confounders = ~ L + A,
      estimator = "ipw"
    ),
    class = "causatr_bare_treatment_in_confounders"
  )
})

test_that("IPW accepts A:modifier and stores em_info", {
  d <- simulate_effect_mod(n = 500, seed = 1)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L + sex + A:sex,
    estimator = "ipw"
  )
  expect_s3_class(fit, "causatr_fit")
  expect_true(fit$details$em_info$has_em)
  expect_equal(fit$details$em_info$modifier_vars, "sex")
})

# Truth-based: IPW sandwich with binary treatment × binary modifier.
#
# DGP 4: Y = 2 + 3*A + 1.5*L + 1.2*sex*A + N(0,1)
#   True ATE|sex=0 = 3
#   True ATE|sex=1 = 4.2
#
# The propensity model is correctly specified (A ~ L), and the MSM
# expands to `Y ~ 1 + sex` per intervention. Under `by = "sex"`,
# `contrast()` averages the modifier-aware predictions per stratum.
test_that("IPW EM sandwich recovers stratum-specific ATEs (DGP 4)", {
  d <- simulate_effect_mod(n = 5000, seed = 42)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L + sex + A:sex,
    estimator = "ipw"
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    reference = "a0",
    ci_method = "sandwich",
    by = "sex"
  )

  ate_sex0 <- res$contrasts$estimate[res$contrasts$by == 0]
  ate_sex1 <- res$contrasts$estimate[res$contrasts$by == 1]

  # Point estimates within 0.3 of truth (the DGP is linear with
  # correctly specified propensity, so the Hajek estimator is
  # consistent).
  expect_lt(abs(ate_sex0 - 3.0), 0.3)
  expect_lt(abs(ate_sex1 - 4.2), 0.3)

  # CIs cover the truth.
  ci_sex0 <- res$contrasts[res$contrasts$by == 0, ]
  ci_sex1 <- res$contrasts[res$contrasts$by == 1, ]
  expect_lte(ci_sex0$ci_lower, 3.0)
  expect_gte(ci_sex0$ci_upper, 3.0)
  expect_lte(ci_sex1$ci_lower, 4.2)
  expect_gte(ci_sex1$ci_upper, 4.2)
})

# Cross-check: IPW EM point estimates agree with gcomp EM.
test_that("IPW EM agrees with gcomp EM on same DGP (DGP 4)", {
  d <- simulate_effect_mod(n = 5000, seed = 42)

  fit_gc <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L + sex + A:sex
  )
  res_gc <- contrast(
    fit_gc,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    reference = "a0",
    ci_method = "sandwich",
    by = "sex"
  )

  fit_ipw <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L + sex + A:sex,
    estimator = "ipw"
  )
  res_ipw <- contrast(
    fit_ipw,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    reference = "a0",
    ci_method = "sandwich",
    by = "sex"
  )

  # Stratum-specific ATEs should agree within cross-method tolerance.
  # Both are consistent under the linear DGP; differences come from
  # finite-sample variance only.
  gc_sex0 <- res_gc$contrasts$estimate[res_gc$contrasts$by == 0]
  gc_sex1 <- res_gc$contrasts$estimate[res_gc$contrasts$by == 1]
  ipw_sex0 <- res_ipw$contrasts$estimate[res_ipw$contrasts$by == 0]
  ipw_sex1 <- res_ipw$contrasts$estimate[res_ipw$contrasts$by == 1]
  expect_lt(abs(gc_sex0 - ipw_sex0), 0.5)
  expect_lt(abs(gc_sex1 - ipw_sex1), 0.5)
})

# Bootstrap: IPW EM bootstrap CIs cover the truth.
test_that("IPW EM bootstrap covers stratum-specific ATEs (DGP 4)", {
  d <- simulate_effect_mod(n = 3000, seed = 42)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L + sex + A:sex,
    estimator = "ipw"
  )
  res <- suppressWarnings(contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    reference = "a0",
    ci_method = "bootstrap",
    n_boot = 100,
    by = "sex"
  ))

  ci_sex0 <- res$contrasts[res$contrasts$by == 0, ]
  ci_sex1 <- res$contrasts[res$contrasts$by == 1, ]
  expect_lte(ci_sex0$ci_lower, 3.0)
  expect_gte(ci_sex0$ci_upper, 3.0)
  expect_lte(ci_sex1$ci_lower, 4.2)
  expect_gte(ci_sex1$ci_upper, 4.2)
})

# Regression guard: IPW without EM terms produces identical results.
test_that("IPW without EM terms is unaffected by EM infrastructure", {
  d <- simulate_binary_continuous(n = 2000, seed = 42)

  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw"
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    reference = "a0",
    ci_method = "sandwich"
  )

  # The standard DGP 1 has ATE = 3.
  expect_lt(abs(res$contrasts$estimate[1] - 3.0), 0.3)
  expect_lte(res$contrasts$ci_lower[1], 3.0)
  expect_gte(res$contrasts$ci_upper[1], 3.0)

  # em_info should report no EM terms.
  expect_false(fit$details$em_info$has_em)
})

# Matching effect modification (chunk 6c) ------------------------------------

test_that("matching rejects bare treatment in confounders", {
  skip_if_not_installed("MatchIt")
  d <- simulate_binary_continuous(n = 200, seed = 2)
  expect_error(
    causat(
      d,
      outcome = "Y",
      treatment = "A",
      confounders = ~ L + A,
      estimator = "matching"
    ),
    class = "causatr_bare_treatment_in_confounders"
  )
})

test_that("matching accepts A:modifier and stores em_info", {
  skip_if_not_installed("MatchIt")
  d <- simulate_effect_mod(n = 500, seed = 1)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L + sex + A:sex,
    estimator = "matching"
  )
  expect_s3_class(fit, "causatr_fit")
  expect_true(fit$details$em_info$has_em)
  expect_equal(fit$details$em_info$modifier_vars, "sex")
})

# Truth-based: Matching sandwich with binary treatment x binary modifier.
#
# DGP 4: Y = 2 + 3*A + 1.5*L + 1.2*sex*A + N(0,1)
#   True ATE|sex=0 = 3
#   True ATE|sex=1 = 4.2
#
# The matching MSM expands to `Y ~ A + sex + A:sex`. Under `by = "sex"`,
# `contrast()` averages the modifier-aware predictions per stratum.
# MatchIt::matchit uses full matching for ATE.
test_that("matching EM sandwich recovers stratum-specific ATEs (DGP 4)", {
  skip_if_not_installed("MatchIt")
  skip_if_not_installed("optmatch")
  d <- simulate_effect_mod(n = 5000, seed = 42)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L + sex + A:sex,
    estimator = "matching"
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    reference = "a0",
    ci_method = "sandwich",
    by = "sex"
  )

  ate_sex0 <- res$contrasts$estimate[res$contrasts$by == 0]
  ate_sex1 <- res$contrasts$estimate[res$contrasts$by == 1]

  # Point estimates within 0.5 of truth. Full matching has finite-sample
  # bias on the order of O(1/n) (Abadie & Imbens 2011), so the CI may
  # not always cover the exact truth — assert the point estimate is in
  # the right ballpark.
  expect_lt(abs(ate_sex0 - 3.0), 0.5)
  expect_lt(abs(ate_sex1 - 4.2), 0.5)

  # SEs are finite and positive.
  expect_true(all(res$contrasts$se > 0))
  expect_true(all(is.finite(res$contrasts$se)))
})

# Cross-check: matching EM agrees with gcomp EM on same DGP.
test_that("matching EM agrees with gcomp EM on same DGP (DGP 4)", {
  skip_if_not_installed("MatchIt")
  skip_if_not_installed("optmatch")
  d <- simulate_effect_mod(n = 5000, seed = 42)

  fit_gc <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L + sex + A:sex
  )
  res_gc <- contrast(
    fit_gc,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    reference = "a0",
    ci_method = "sandwich",
    by = "sex"
  )

  fit_m <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L + sex + A:sex,
    estimator = "matching"
  )
  res_m <- contrast(
    fit_m,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    reference = "a0",
    ci_method = "sandwich",
    by = "sex"
  )

  gc_sex0 <- res_gc$contrasts$estimate[res_gc$contrasts$by == 0]
  gc_sex1 <- res_gc$contrasts$estimate[res_gc$contrasts$by == 1]
  m_sex0 <- res_m$contrasts$estimate[res_m$contrasts$by == 0]
  m_sex1 <- res_m$contrasts$estimate[res_m$contrasts$by == 1]
  # Cross-method tolerance: matching is less efficient, but both are
  # consistent, so estimates should agree within ~0.5 on n=5000.
  expect_lt(abs(gc_sex0 - m_sex0), 0.5)
  expect_lt(abs(gc_sex1 - m_sex1), 0.5)
})

# Bootstrap: matching EM bootstrap CIs cover the truth.
test_that("matching EM bootstrap covers stratum-specific ATEs (DGP 4)", {
  skip_if_not_installed("MatchIt")
  skip_if_not_installed("optmatch")
  d <- simulate_effect_mod(n = 3000, seed = 42)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L + sex + A:sex,
    estimator = "matching"
  )
  res <- suppressWarnings(contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    reference = "a0",
    ci_method = "bootstrap",
    n_boot = 100,
    by = "sex"
  ))

  ci_sex0 <- res$contrasts[res$contrasts$by == 0, ]
  ci_sex1 <- res$contrasts[res$contrasts$by == 1, ]
  expect_lte(ci_sex0$ci_lower, 3.0)
  expect_gte(ci_sex0$ci_upper, 3.0)
  expect_lte(ci_sex1$ci_lower, 4.2)
  expect_gte(ci_sex1$ci_upper, 4.2)
})

# Regression guard: matching without EM terms is unaffected.
test_that("matching without EM terms is unaffected by EM infrastructure", {
  skip_if_not_installed("MatchIt")
  d <- simulate_binary_continuous(n = 2000, seed = 42)

  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "matching"
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    reference = "a0",
    ci_method = "sandwich"
  )

  # DGP 1 has ATE = 3.
  expect_lt(abs(res$contrasts$estimate[1] - 3.0), 0.5)
  expect_lte(res$contrasts$ci_lower[1], 3.0)
  expect_gte(res$contrasts$ci_upper[1], 3.0)

  # em_info should report no EM terms.
  expect_false(fit$details$em_info$has_em)
})
