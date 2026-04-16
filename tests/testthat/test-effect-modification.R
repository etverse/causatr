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

# IPW rejection of EM terms (temporary, until chunk 6b) ---------------------

test_that("IPW rejects A:modifier interaction terms in confounders", {
  d <- simulate_binary_continuous(n = 200, seed = 1)
  d$sex <- rbinom(nrow(d), 1, 0.5)
  expect_error(
    causat(
      d,
      outcome = "Y",
      treatment = "A",
      confounders = ~ L + sex + A:sex,
      estimator = "ipw"
    ),
    class = "causatr_em_unsupported"
  )
})

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

# Matching rejection of EM terms (temporary, until chunk 6c) ----------------

test_that("matching rejects A:modifier interaction terms in confounders", {
  skip_if_not_installed("MatchIt")
  d <- simulate_binary_continuous(n = 200, seed = 1)
  d$sex <- rbinom(nrow(d), 1, 0.5)
  expect_error(
    causat(
      d,
      outcome = "Y",
      treatment = "A",
      confounders = ~ L + sex + A:sex,
      estimator = "matching"
    ),
    class = "causatr_em_unsupported"
  )
})

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
