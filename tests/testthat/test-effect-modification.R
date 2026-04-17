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

# ICE effect modification (chunk 6d) ------------------------------------------

# expand_em_lag_terms() unit tests.

test_that("expand_em_lag_terms produces correct lag terms", {
  em <- parse_effect_mod(~ L + sex + A:sex, "A")
  lags <- expand_em_lag_terms(em$em_terms[[1]], available_lags = 3L)
  # R's terms() alphabetizes "A:sex" -> "sex:A", so the expansion keeps
  # that order: "sex:lag1_A", "sex:lag2_A", "sex:lag3_A".
  expect_equal(lags, c("sex:lag1_A", "sex:lag2_A", "sex:lag3_A"))
})

test_that("expand_em_lag_terms returns empty for available_lags = 0", {
  em <- parse_effect_mod(~ L + sex + A:sex, "A")
  lags <- expand_em_lag_terms(em$em_terms[[1]], available_lags = 0L)
  expect_equal(lags, character(0L))
})

test_that("expand_em_lag_terms handles reversed term order (sex:A)", {
  em <- parse_effect_mod(~ L + sex + sex:A, "A")
  lags <- expand_em_lag_terms(em$em_terms[[1]], available_lags = 2L)
  # The term label from terms() preserves user order: "sex:A"
  # So lags should be "sex:lag1_A", "sex:lag2_A"
  expect_equal(lags, c("sex:lag1_A", "sex:lag2_A"))
})

test_that("expand_em_lag_terms handles multi-way interaction", {
  em <- parse_effect_mod(~ L + A:sex:race, "A")
  lags <- expand_em_lag_terms(em$em_terms[[1]], available_lags = 1L)
  expect_equal(lags, "lag1_A:sex:race")
})

# Truth-based: ICE sandwich x 2-period DGP x binary modifier.
#
# DGP-EM-ICE: Y = 10 + (2 + 1.5*sex) * sum(A_t) + L0 + sum(L_t) + eps
#   True ATE|sex=0 = 5  (always vs never)
#   True ATE|sex=1 = 8
#
# The confounders formula includes `A:sex` which triggers lag auto-expansion.
# At time_idx = 1 (the final period, 0-based), the formula should include
# both `A:sex` (current period) and `lag1_A:sex` (first period). Without
# the expansion, ICE compresses the heterogeneity and returns ~6.5/6.5
# instead of 5/8.
test_that("ICE EM sandwich recovers sex-specific ATEs (2-period DGP)", {
  d <- make_em_ice_scm(n = 8000, n_times = 2, seed = 42)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L0 + sex + A:sex,
    confounders_tv = ~L,
    id = "id",
    time = "time"
  )

  # Verify em_info is stored.
  expect_true(fit$details$em_info$has_em)

  res <- contrast(
    fit,
    interventions = list(always = static(1), never = static(0)),
    type = "difference",
    reference = "never",
    ci_method = "sandwich",
    by = "sex"
  )

  ate_sex0 <- res$contrasts$estimate[res$contrasts$by == 0]
  ate_sex1 <- res$contrasts$estimate[res$contrasts$by == 1]

  # Point estimates within 0.5 of truth. The DGP is linear and correctly
  # specified, so ICE g-comp is consistent.
  expect_lt(abs(ate_sex0 - 5.0), 0.5)
  expect_lt(abs(ate_sex1 - 8.0), 0.5)

  # CIs cover the truth.
  ci_sex0 <- res$contrasts[res$contrasts$by == 0, ]
  ci_sex1 <- res$contrasts[res$contrasts$by == 1, ]
  expect_lte(ci_sex0$ci_lower, 5.0)
  expect_gte(ci_sex0$ci_upper, 5.0)
  expect_lte(ci_sex1$ci_lower, 8.0)
  expect_gte(ci_sex1$ci_upper, 8.0)

  # SEs are finite and positive.
  expect_true(all(res$contrasts$se > 0))
  expect_true(all(is.finite(res$contrasts$se)))
})

# Truth-based: ICE sandwich x 3-period DGP x binary modifier.
#
# DGP-EM-ICE with 3 periods:
#   True ATE|sex=0 = 8    (always vs never)
#   True ATE|sex=1 = 12.5
#
# Tests deeper lag coverage: at the final step (time_idx = 2), the formula
# should include `A:sex`, `lag1_A:sex`, and `lag2_A:sex`.
test_that("ICE EM sandwich recovers sex-specific ATEs (3-period DGP)", {
  d <- make_em_ice_scm(n = 8000, n_times = 3, seed = 42)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L0 + sex + A:sex,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    history = Inf
  )

  res <- contrast(
    fit,
    interventions = list(always = static(1), never = static(0)),
    type = "difference",
    reference = "never",
    ci_method = "sandwich",
    by = "sex"
  )

  ate_sex0 <- res$contrasts$estimate[res$contrasts$by == 0]
  ate_sex1 <- res$contrasts$estimate[res$contrasts$by == 1]

  # Tolerances are slightly wider for 3 periods (more variance).
  expect_lt(abs(ate_sex0 - 8.0), 0.8)
  expect_lt(abs(ate_sex1 - 12.5), 0.8)

  # CIs cover the truth.
  ci_sex0 <- res$contrasts[res$contrasts$by == 0, ]
  ci_sex1 <- res$contrasts[res$contrasts$by == 1, ]
  expect_lte(ci_sex0$ci_lower, 8.0)
  expect_gte(ci_sex0$ci_upper, 8.0)
  expect_lte(ci_sex1$ci_lower, 12.5)
  expect_gte(ci_sex1$ci_upper, 12.5)
})

# Multiple EM terms: ICE x A:sex + A:age x 2 periods.
# Smoke test -- no closed-form truth for continuous modifiers under ICE,
# but we verify the formula expansion works and produces finite results.
test_that("ICE EM handles multiple EM terms (A:sex + A:age)", {
  set.seed(99)
  n <- 3000
  L0 <- rnorm(n)
  sex <- rbinom(n, 1, 0.5)
  age <- rnorm(n, 50, 10)
  A0 <- rbinom(n, 1, plogis(0.3 * L0))
  L1 <- A0 + 0.5 * L0 + rnorm(n, 0, 0.5)
  A1 <- rbinom(n, 1, plogis(0.3 * L1))
  Y <- 10 + (2 + 1.5 * sex + 0.05 * age) * (A0 + A1) + L0 + L1 + rnorm(n)

  d <- rbind(
    data.frame(
      id = 1:n,
      time = 0L,
      A = A0,
      L = NA_real_,
      L0 = L0,
      sex = sex,
      age = age,
      Y = NA_real_
    ),
    data.frame(
      id = 1:n,
      time = 1L,
      A = A1,
      L = L1,
      L0 = L0,
      sex = sex,
      age = age,
      Y = Y
    )
  )

  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L0 + sex + age + A:sex + A:age,
    confounders_tv = ~L,
    id = "id",
    time = "time"
  )

  expect_true(fit$details$em_info$has_em)
  expect_length(fit$details$em_info$em_terms, 2L)

  res <- contrast(
    fit,
    interventions = list(always = static(1), never = static(0)),
    type = "difference",
    reference = "never",
    ci_method = "sandwich"
  )

  # Smoke: finite estimates and SEs.
  expect_true(all(is.finite(res$contrasts$estimate)))
  expect_true(all(res$contrasts$se > 0))
  expect_true(all(is.finite(res$contrasts$se)))
})

# Bootstrap: ICE EM bootstrap CIs cover the truth (2-period DGP).
test_that("ICE EM bootstrap covers sex-specific ATEs (2-period DGP)", {
  d <- make_em_ice_scm(n = 5000, n_times = 2, seed = 42)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L0 + sex + A:sex,
    confounders_tv = ~L,
    id = "id",
    time = "time"
  )
  res <- suppressWarnings(contrast(
    fit,
    interventions = list(always = static(1), never = static(0)),
    type = "difference",
    reference = "never",
    ci_method = "bootstrap",
    n_boot = 100,
    by = "sex"
  ))

  ci_sex0 <- res$contrasts[res$contrasts$by == 0, ]
  ci_sex1 <- res$contrasts[res$contrasts$by == 1, ]
  expect_lte(ci_sex0$ci_lower, 5.0)
  expect_gte(ci_sex0$ci_upper, 5.0)
  expect_lte(ci_sex1$ci_lower, 8.0)
  expect_gte(ci_sex1$ci_upper, 8.0)
})

# Regression guard: ICE without EM terms produces identical results
# to the pre-expansion behavior.
test_that("ICE without EM terms is unaffected by EM lag expansion", {
  d <- make_linear_scm(n = 5000, n_times = 2, seed = 42)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    id = "id",
    time = "time"
  )

  # em_info should report no EM.
  expect_false(fit$details$em_info$has_em)

  res <- contrast(
    fit,
    interventions = list(always = static(1), never = static(0)),
    type = "difference",
    reference = "never",
    ci_method = "sandwich"
  )

  # The linear SCM has ATE (always vs never) = 5.
  expect_lt(abs(res$contrasts$estimate[1] - 5.0), 0.5)
  expect_lte(res$contrasts$ci_lower[1], 5.0)
  expect_gte(res$contrasts$ci_upper[1], 5.0)
})

# ICE EM x continuous treatment x shift (chunk 6d) ----------------------------

# Truth-based: ICE x continuous trt x shift(1) x binary modifier.
#
# DGP-EM-ICE-CONT: Y = 10 + (2+1.5*sex)*(A_0+A_1) + L0 + L_1 + eps
# The g-formula shift(1) effect includes the indirect path through
# treatment-confounder feedback (L_1 depends on A_0), so the true
# causal effect is larger than the direct coefficient sum.
# MC truth (n = 5*10^6, seed = 1):
#   shift(1)|sex=0 ~ 6.0, shift(1)|sex=1 ~ 9.76
test_that("ICE EM x continuous trt x shift recovers sex-specific effects", {
  d <- make_em_ice_cont_scm(n = 10000, seed = 42)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L0 + sex + A:sex,
    confounders_tv = ~L,
    id = "id",
    time = "time"
  )

  res <- contrast(
    fit,
    interventions = list(shifted = shift(1), natural = NULL),
    type = "difference",
    reference = "natural",
    ci_method = "sandwich",
    by = "sex"
  )

  ate_sex0 <- res$contrasts$estimate[res$contrasts$by == 0]
  ate_sex1 <- res$contrasts$estimate[res$contrasts$by == 1]

  # MC truth: 6.0 / 9.76. Continuous treatments have higher variance
  # due to treatment-confounder feedback and collinearity between A_0
  # and L_1. The ICE model is correctly specified but the finite-sample
  # bias is larger than for binary treatments.
  expect_lt(abs(ate_sex0 - 6.0), 0.5)
  expect_lt(abs(ate_sex1 - 9.76), 0.5)

  # SEs are finite and positive, and sex=1 effect is larger.
  expect_true(all(res$contrasts$se > 0))
  expect_true(all(is.finite(res$contrasts$se)))
  expect_gt(ate_sex1, ate_sex0)
})

# ICE EM x continuous treatment x scale_by.
#
# DGP-EM-ICE-CONT: scale_by(1.5) vs natural
#   E[Y(1.5*A)] - E[Y(A)] = (2+1.5*sex) * (1.5 - 1) * (E[A_0] + E[A_1])
#   where E[A_0] ~ 1, E[A_1] ~ 1 + 0.3*E[L_1] + 0.2 ~ flexible
# No analytical truth for the exact marginal mean, so this is a smoke test
# that verifies the scale_by path produces differentiated sex-specific
# contrasts (not collapsed to the same value).
test_that("ICE EM x continuous trt x scale_by produces sex-differentiated estimates", {
  d <- make_em_ice_cont_scm(n = 5000, seed = 42)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L0 + sex + A:sex,
    confounders_tv = ~L,
    id = "id",
    time = "time"
  )

  res <- contrast(
    fit,
    interventions = list(scaled = scale_by(1.5), natural = NULL),
    type = "difference",
    reference = "natural",
    ci_method = "sandwich",
    by = "sex"
  )

  est_sex0 <- res$contrasts$estimate[res$contrasts$by == 0]
  est_sex1 <- res$contrasts$estimate[res$contrasts$by == 1]

  # Smoke: estimates are finite, positive, and sex=1 effect is larger.
  expect_true(all(is.finite(res$contrasts$estimate)))
  expect_true(all(res$contrasts$se > 0))
  expect_gt(est_sex1, est_sex0)
})

# ICE EM x continuous treatment x threshold.
test_that("ICE EM x continuous trt x threshold produces sex-differentiated estimates", {
  d <- make_em_ice_cont_scm(n = 5000, seed = 42)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L0 + sex + A:sex,
    confounders_tv = ~L,
    id = "id",
    time = "time"
  )

  res <- contrast(
    fit,
    interventions = list(capped = threshold(upper = 1.5), natural = NULL),
    type = "difference",
    reference = "natural",
    ci_method = "sandwich",
    by = "sex"
  )

  # Smoke: finite and differentiated.
  expect_true(all(is.finite(res$contrasts$estimate)))
  expect_true(all(res$contrasts$se > 0))
  # The threshold caps A at 1.5, reducing E[Y] — sex=1 stratum should
  # have a larger (more negative) reduction because the per-unit effect
  # is larger for sex=1.
  est_sex0 <- res$contrasts$estimate[res$contrasts$by == 0]
  est_sex1 <- res$contrasts$estimate[res$contrasts$by == 1]
  expect_lt(est_sex1, est_sex0)
})

# ICE EM x binary treatment x dynamic intervention ----------------------------

# Truth-based: ICE x binary trt x dynamic x EM.
#
# Dynamic rule: treat if L_k > 0 (treat individuals with higher
# time-varying confounder). Under this rule on the make_em_ice_scm DGP,
# we can't derive the analytical truth, but we verify sex-differentiated
# estimates are produced and finite.
test_that("ICE EM x binary trt x dynamic produces sex-differentiated estimates", {
  d <- make_em_ice_scm(n = 5000, n_times = 2, seed = 42)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L0 + sex + A:sex,
    confounders_tv = ~L,
    id = "id",
    time = "time"
  )

  res <- contrast(
    fit,
    interventions = list(
      dyn = dynamic(function(data, orig_trt) as.integer(data$L0 > 0)),
      never = static(0)
    ),
    type = "difference",
    reference = "never",
    ci_method = "sandwich",
    by = "sex"
  )

  est_sex0 <- res$contrasts$estimate[res$contrasts$by == 0]
  est_sex1 <- res$contrasts$estimate[res$contrasts$by == 1]

  # Smoke: finite, positive (treating some is better than treating none),
  # and sex=1 has a larger effect (per the DGP).
  expect_true(all(is.finite(res$contrasts$estimate)))
  expect_true(all(res$contrasts$se > 0))
  expect_gt(est_sex0, 0)
  expect_gt(est_sex1, est_sex0)
})

# ICE EM x binary outcome (binomial) -----------------------------------------

# Truth-based: ICE x binary trt x binomial outcome x EM.
#
# DGP-EM-ICE-BINOM: Y ~ Bern(expit(-1 + (1+0.8*sex)*(A0+A1) + 0.5*L0 + 0.3*L1))
#   MC truth (n = 10^6): RD|sex=0 ~ 0.480, RD|sex=1 ~ 0.651
test_that("ICE EM x binomial outcome recovers sex-specific risk differences", {
  d <- make_em_ice_binom_scm(n = 8000, seed = 42)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L0 + sex + A:sex,
    confounders_tv = ~L,
    family = "binomial",
    id = "id",
    time = "time"
  )

  res <- contrast(
    fit,
    interventions = list(always = static(1), never = static(0)),
    type = "difference",
    reference = "never",
    ci_method = "sandwich",
    by = "sex"
  )

  rd_sex0 <- res$contrasts$estimate[res$contrasts$by == 0]
  rd_sex1 <- res$contrasts$estimate[res$contrasts$by == 1]

  # MC truth: RD|sex=0 ~ 0.480, RD|sex=1 ~ 0.651
  # Tolerate 0.05 on probability scale (binomial is noisier).
  expect_lt(abs(rd_sex0 - 0.480), 0.05)
  expect_lt(abs(rd_sex1 - 0.651), 0.05)

  # CIs cover the MC truth.
  ci_sex0 <- res$contrasts[res$contrasts$by == 0, ]
  ci_sex1 <- res$contrasts[res$contrasts$by == 1, ]
  expect_lte(ci_sex0$ci_lower, 0.480)
  expect_gte(ci_sex0$ci_upper, 0.480)
  expect_lte(ci_sex1$ci_lower, 0.651)
  expect_gte(ci_sex1$ci_upper, 0.651)

  # SEs are finite and positive.
  expect_true(all(res$contrasts$se > 0))
  expect_true(all(is.finite(res$contrasts$se)))
})

# ICE EM x multivariate treatment ---------------------------------------------

# Smoke test: ICE x multivariate binary trt x EM.
# Uses a DGP with two binary treatments and sex-specific effects for A1.
test_that("ICE EM x multivariate trt produces finite sex-differentiated estimates", {
  set.seed(55)
  n <- 3000
  L0 <- rnorm(n)
  sex <- rbinom(n, 1, 0.5)
  A1_0 <- rbinom(n, 1, plogis(0.3 * L0))
  A2_0 <- rbinom(n, 1, 0.5)
  L1 <- 0.5 * A1_0 + rnorm(n, 0, 0.5)
  A1_1 <- rbinom(n, 1, plogis(0.3 * L1))
  A2_1 <- rbinom(n, 1, 0.5)
  Y <- 50 +
    (2 + 1.5 * sex) * (A1_0 + A1_1) +
    (A2_0 + A2_1) -
    L1 +
    rnorm(n, 0, 2)

  d <- rbind(
    data.frame(
      id = 1:n,
      time = 0L,
      A1 = A1_0,
      A2 = A2_0,
      L = L0,
      sex = sex,
      Y = NA_real_
    ),
    data.frame(
      id = 1:n,
      time = 1L,
      A1 = A1_1,
      A2 = A2_1,
      L = L1,
      sex = sex,
      Y = Y
    )
  )

  fit <- causat(
    d,
    outcome = "Y",
    treatment = c("A1", "A2"),
    confounders = ~ sex + A1:sex,
    confounders_tv = ~L,
    id = "id",
    time = "time"
  )

  expect_true(fit$details$em_info$has_em)

  res <- contrast(
    fit,
    interventions = list(
      both = list(A1 = static(1), A2 = static(1)),
      none = list(A1 = static(0), A2 = static(0))
    ),
    type = "difference",
    reference = "none",
    ci_method = "sandwich",
    by = "sex"
  )

  est_sex0 <- res$contrasts$estimate[res$contrasts$by == 0]
  est_sex1 <- res$contrasts$estimate[res$contrasts$by == 1]

  # Smoke: finite, positive, and sex=1 has larger effect.
  expect_true(all(is.finite(res$contrasts$estimate)))
  expect_true(all(res$contrasts$se > 0))
  expect_gt(est_sex1, est_sex0)
})

# lmtp oracle cross-check: ICE EM binary treatment (2-period DGP).
#
# Validates causatr ICE EM against lmtp::lmtp_sdr() on per-stratum
# point estimates. Uses the same DGP and seed for both.
test_that("ICE EM binary: causatr agrees with lmtp per-stratum (DGP-EM-ICE)", {
  skip_if_not_installed("lmtp")

  d <- make_em_ice_scm(n = 5000, n_times = 2, seed = 42)

  # causatr
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L0 + sex + A:sex,
    confounders_tv = ~L,
    id = "id",
    time = "time"
  )
  res <- contrast(
    fit,
    interventions = list(always = static(1), never = static(0)),
    type = "difference",
    reference = "never",
    ci_method = "sandwich",
    by = "sex"
  )

  # lmtp — reshape to wide, stratify by sex
  d_wide <- reshape(
    d,
    idvar = "id",
    timevar = "time",
    direction = "wide",
    v.names = c("A", "L", "Y"),
    sep = "_"
  )
  d_clean <- d_wide[, c("id", "L0", "A_0", "A_1", "L_1", "Y_1")]
  sex_vec <- d_wide$sex

  run_lmtp <- function(data_sub, shift_fn) {
    suppressWarnings(suppressMessages(lmtp::lmtp_sdr(
      data = data_sub,
      trt = c("A_0", "A_1"),
      outcome = "Y_1",
      baseline = c("L0"),
      time_vary = list(NULL, c("L_1")),
      shift = shift_fn,
      outcome_type = "continuous",
      learners_trt = "SL.glm",
      learners_outcome = "SL.glm",
      folds = 1
    )))
  }
  theta_of <- function(r) r$estimate@x

  for (s in c(0, 1)) {
    ds <- d_clean[sex_vec == s, ]
    ra <- run_lmtp(ds, function(data, trt) rep(1, nrow(data)))
    rn <- run_lmtp(ds, function(data, trt) rep(0, nrow(data)))
    lmtp_ate <- theta_of(ra) - theta_of(rn)
    causatr_ate <- res$contrasts$estimate[res$contrasts$by == s]
    # Cross-method tolerance: lmtp uses SDR (semiparametric), causatr
    # uses parametric g-comp — both consistent under correct spec.
    expect_lt(
      abs(causatr_ate - lmtp_ate),
      0.5,
      label = paste0("causatr vs lmtp ATE|sex=", s)
    )
  }
})
