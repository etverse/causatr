# Replication of Hernán & Robins "Causal Inference: What If" (2025)
# Chapters 12–13 using the bundled NHEFS dataset.
#
# Published estimates for qsmk (quit smoking) on wt82_71 (weight change):
#   IPW (Hájek, stabilized):  ATE ≈ 3.44 kg
#   G-computation:            ATE ≈ 3.5 kg
#   AIPW (delicatessen ref):  ATE ≈ 3.46 kg
#
# Complete-case analysis (censored == 0, no education NAs).

# Complete-case subsample: drop censored rows and education NAs.
make_nhefs_cc <- function() {
  data("nhefs", package = "causatr", envir = environment())
  nhefs[nhefs$censored == 0 & !is.na(nhefs$education), ]
}

# Shared confounder formula from Hernán & Robins, Chapters 12–13.
nhefs_confounders <- ~ sex +
  age +
  I(age^2) +
  race +
  factor(education) +
  smokeintensity +
  I(smokeintensity^2) +
  smokeyrs +
  I(smokeyrs^2) +
  factor(exercise) +
  factor(active) +
  wt71 +
  I(wt71^2) +
  qsmk:smokeintensity


# -- G-computation (Chapter 13) -------------------------------------------

test_that("NHEFS g-comp ATE ≈ 3.5 kg (Hernán & Robins Ch 13)", {
  nhefs_cc <- make_nhefs_cc()

  fit <- causat(
    nhefs_cc,
    outcome = "wt82_71",
    treatment = "qsmk",
    confounders = nhefs_confounders
  )

  res <- contrast(
    fit,
    interventions = list(quit = static(1), cont = static(0)),
    reference = "cont",
    ci_method = "sandwich"
  )

  ate <- res$contrasts$estimate[1]
  se <- res$contrasts$se[1]

  expect_equal(ate, 3.5, tolerance = 0.1)
  expect_true(is.finite(se) && se > 0)
  expect_lt(res$contrasts$ci_lower[1], 3.5)
  expect_gt(res$contrasts$ci_upper[1], 3.5)
})


# -- IPW (Chapter 12) -----------------------------------------------------

test_that("NHEFS IPW ATE ≈ 3.44 kg (Hernán & Robins Ch 12)", {
  nhefs_cc <- make_nhefs_cc()

  fit <- causat(
    nhefs_cc,
    outcome = "wt82_71",
    treatment = "qsmk",
    confounders = nhefs_confounders,
    estimator = "ipw",
    propensity_model_fn = stats::glm
  )

  res <- contrast(
    fit,
    interventions = list(quit = static(1), cont = static(0)),
    reference = "cont",
    ci_method = "sandwich"
  )

  ate <- res$contrasts$estimate[1]
  se <- res$contrasts$se[1]

  expect_equal(ate, 3.44, tolerance = 0.15)
  expect_true(is.finite(se) && se > 0)
  expect_lt(res$contrasts$ci_lower[1], 3.5)
  expect_gt(res$contrasts$ci_upper[1], 3.5)
})


# -- AIPW -----------------------------------------------------------------

test_that("NHEFS AIPW ATE ≈ 3.46 kg", {
  nhefs_cc <- make_nhefs_cc()

  fit <- causat(
    nhefs_cc,
    outcome = "wt82_71",
    treatment = "qsmk",
    confounders = nhefs_confounders,
    estimator = "aipw",
    propensity_model_fn = stats::glm
  )

  res <- contrast(
    fit,
    interventions = list(quit = static(1), cont = static(0)),
    reference = "cont",
    ci_method = "sandwich"
  )

  ate <- res$contrasts$estimate[1]
  se <- res$contrasts$se[1]

  expect_equal(ate, 3.46, tolerance = 0.1)
  expect_true(is.finite(se) && se > 0)
  expect_lt(res$contrasts$ci_lower[1], 3.5)
  expect_gt(res$contrasts$ci_upper[1], 3.5)
})


# -- Cross-estimator agreement on NHEFS -----------------------------------

test_that("NHEFS: g-comp, IPW, AIPW agree within 0.1 kg", {
  nhefs_cc <- make_nhefs_cc()

  ivs <- list(quit = static(1), cont = static(0))

  fit_gc <- causat(
    nhefs_cc,
    outcome = "wt82_71",
    treatment = "qsmk",
    confounders = nhefs_confounders
  )
  fit_ipw <- causat(
    nhefs_cc,
    outcome = "wt82_71",
    treatment = "qsmk",
    confounders = nhefs_confounders,
    estimator = "ipw",
    propensity_model_fn = stats::glm
  )
  fit_aipw <- causat(
    nhefs_cc,
    outcome = "wt82_71",
    treatment = "qsmk",
    confounders = nhefs_confounders,
    estimator = "aipw",
    propensity_model_fn = stats::glm
  )

  res_gc <- contrast(fit_gc, interventions = ivs, reference = "cont")
  res_ipw <- contrast(fit_ipw, interventions = ivs, reference = "cont")
  res_aipw <- contrast(fit_aipw, interventions = ivs, reference = "cont")

  ate_gc <- res_gc$contrasts$estimate[1]
  ate_ipw <- res_ipw$contrasts$estimate[1]
  ate_aipw <- res_aipw$contrasts$estimate[1]

  expect_lt(abs(ate_gc - ate_ipw), 0.1)
  expect_lt(abs(ate_gc - ate_aipw), 0.1)
  expect_lt(abs(ate_ipw - ate_aipw), 0.1)
})


# -- Matching on NHEFS ----------------------------------------------------

test_that("NHEFS matching ATT ≈ g-comp ATT (within 0.5 kg)", {
  nhefs_cc <- make_nhefs_cc()

  ivs <- list(quit = static(1), cont = static(0))

  fit_match <- causat(
    nhefs_cc,
    outcome = "wt82_71",
    treatment = "qsmk",
    confounders = nhefs_confounders,
    estimator = "matching",
    estimand = "ATT"
  )
  fit_gc <- causat(
    nhefs_cc,
    outcome = "wt82_71",
    treatment = "qsmk",
    confounders = nhefs_confounders
  )

  res_match <- contrast(
    fit_match,
    interventions = ivs,
    reference = "cont"
  )
  res_gc <- contrast(
    fit_gc,
    interventions = ivs,
    reference = "cont",
    estimand = "ATT",
    ci_method = "sandwich"
  )

  att_match <- res_match$contrasts$estimate[1]
  att_gc <- res_gc$contrasts$estimate[1]

  expect_gt(att_match, 2)
  expect_lt(att_match, 5)
  expect_lt(abs(att_match - att_gc), 0.5)
})


# -- IPW + IPCW on full NHEFS (including censored rows) -------------------

test_that("NHEFS IPW + IPCW ≈ complete-case IPW (within 0.5 kg)", {
  data("nhefs", package = "causatr", envir = environment())
  nhefs_full <- nhefs[!is.na(nhefs$education), ]

  ivs <- list(quit = static(1), cont = static(0))

  fit_ipcw <- causat(
    nhefs_full,
    outcome = "wt82_71",
    treatment = "qsmk",
    confounders = nhefs_confounders,
    estimator = "ipw",
    propensity_model_fn = stats::glm,
    censoring = "censored",
    ipcw = TRUE
  )
  res_ipcw <- contrast(
    fit_ipcw,
    interventions = ivs,
    reference = "cont",
    ci_method = "sandwich"
  )

  nhefs_cc <- make_nhefs_cc()
  fit_cc <- causat(
    nhefs_cc,
    outcome = "wt82_71",
    treatment = "qsmk",
    confounders = nhefs_confounders,
    estimator = "ipw",
    propensity_model_fn = stats::glm
  )
  res_cc <- contrast(
    fit_cc,
    interventions = ivs,
    reference = "cont",
    ci_method = "sandwich"
  )

  ate_ipcw <- res_ipcw$contrasts$estimate[1]
  ate_cc <- res_cc$contrasts$estimate[1]

  expect_true(is.finite(ate_ipcw))
  expect_true(is.finite(res_ipcw$contrasts$se[1]) && res_ipcw$contrasts$se[1] > 0)
  expect_lt(abs(ate_ipcw - ate_cc), 0.5)
})
