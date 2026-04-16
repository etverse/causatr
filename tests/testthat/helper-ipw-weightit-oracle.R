# WeightIt contrast-level oracle helpers for the static-binary IPW path.
#
# These helpers power the T-oracle1..4 tests in test-ipw-weightit-oracle.R
# and are the sole reason `WeightIt` remains as a test-only dependency
# after chunk 3d removes it from the runtime path. Everything below is
# called behind `testthat::skip_if_not_installed("WeightIt")` so the
# package can install without `WeightIt` available.
#
# Design: contrast-level, not IF-level. PHASE_4_INTERVENTIONS_SELF_IPW.md
# §9 spells out why — causatr and WeightIt use different internal
# decompositions of the same M-estimation system (causatr: HT indicator
# weights with `Y ~ 1` Hájek per arm; WeightIt: full IPW weights with
# saturated `Y ~ A` in a single `glm_weightit`). The per-observation
# influence functions are not elementwise comparable, but the final
# contrast point estimate and standard error agree asymptotically and
# numerically to ~1e-6 on the same data. Comparing at the contrast level
# lets the test be agnostic to the internal decomposition and robust to
# future refactoring of the IF engine.

#' Fit the WeightIt contrast oracle for a static binary IPW comparison
#'
#' @description
#' Runs the canonical `WeightIt::weightit() -> WeightIt::glm_weightit()`
#' pipeline on `data` and returns the marginal means + contrast + SE
#' in a shape that mirrors `contrast()`'s output, so the oracle test
#' can compare `causatr$estimates$estimate` / `causatr$contrasts$se`
#' directly against `oracle$mu_hat` / `oracle$se_contrast`.
#'
#' @param data A data frame with the outcome, treatment, and
#'   confounder columns.
#' @param outcome Character scalar, the outcome column name.
#' @param treatment Character scalar, the treatment column name. Must
#'   be a 0/1 binary integer / numeric.
#' @param ps_formula A one-sided propensity-score formula
#'   (e.g. `~ L + sex`).
#' @param estimand Character scalar in `c("ATE", "ATT", "ATC")`.
#' @param y_family A `stats::family` object for the outcome model
#'   (defaults to `gaussian()`).
#'
#' @return A list with elements:
#'   - `mu_1`, `mu_0` — counterfactual marginal means.
#'   - `contrast` — `mu_1 - mu_0`.
#'   - `se_contrast` — sandwich SE from `glm_weightit`'s M-estimation
#'     variance, accounting for propensity-score uncertainty.
#'   - `weightit_obj`, `glm_obj` — raw fit objects kept for
#'     introspection when a test fails.
weightit_oracle_contrast <- function(
  data,
  outcome,
  treatment,
  ps_formula,
  estimand,
  y_family = stats::gaussian()
) {
  # `weightit` takes a two-sided formula with treatment on the LHS.
  # The RHS is identical to the RHS of `ps_formula` — we rebuild it
  # to keep the helper robust to whatever form the caller passed
  # (`~ L + sex`, `~ s(L)`, etc.; we strip the intercept/env).
  ps_rhs <- paste(deparse(ps_formula[[length(ps_formula)]]), collapse = " ")
  ps_full <- stats::as.formula(paste(treatment, "~", ps_rhs))

  wobj <- WeightIt::weightit(
    ps_full,
    data = data,
    method = "glm",
    estimand = estimand
  )

  # `glm_weightit` fits the saturated MSM `Y ~ A` under the WeightIt
  # weights and returns an M-estimation sandwich vcov that propagates
  # propensity-score uncertainty through the stacked (alpha, beta)
  # system. This is the variance object we want to compare against —
  # it uses exactly the same M-estimation argument that causatr's
  # `compute_ipw_if_self_contained_one()` implements, just through a
  # different weight decomposition.
  y_formula <- stats::as.formula(paste(outcome, "~", treatment))
  gfit <- WeightIt::glm_weightit(
    y_formula,
    data = data,
    weightit = wobj,
    family = y_family
  )

  co <- stats::coef(gfit)
  V <- stats::vcov(gfit)

  # For a 0/1 treatment the saturated `Y ~ A` coefficients decode
  # the marginal means as:
  #   mu_0 = (Intercept)
  #   mu_1 = (Intercept) + coef[treatment]
  # and the contrast mu_1 - mu_0 is simply coef[treatment]. This only
  # holds for a Gaussian / identity-link outcome — for a logit link,
  # `coef[A]` is a log-OR, not a risk difference, and the contrast
  # SE on the response scale is a delta-method transform of the log
  # coefficient. T-oracle1..4 stick to continuous outcomes so this
  # linear decoding is exact.
  mu_0 <- unname(co[["(Intercept)"]])
  mu_1 <- mu_0 + unname(co[[treatment]])
  contrast_val <- unname(co[[treatment]])
  se_contrast <- sqrt(V[treatment, treatment])

  list(
    mu_0 = mu_0,
    mu_1 = mu_1,
    contrast = contrast_val,
    se_contrast = unname(se_contrast),
    weightit_obj = wobj,
    glm_obj = gfit
  )
}


#' Extract the scalar contrast + SE from a causatr_result
#'
#' @description
#' `contrast()` returns a `causatr_result` with `$estimates` (per
#' intervention) and `$contrasts` (the contrast rows). The oracle
#' tests compare a single `static(1) vs static(0)` contrast, so this
#' helper flattens the object into the same `(mu_0, mu_1, contrast,
#' se_contrast)` shape as `weightit_oracle_contrast()`.
#'
#' @param res A `causatr_result` from
#'   `contrast(fit, list(a1 = static(1), a0 = static(0)), reference = "a0")`.
#'
#' @return A list with elements `mu_0`, `mu_1`, `contrast`,
#'   `se_contrast` — structurally parallel to the WeightIt oracle
#'   output so the test bodies can call `expect_equal()` one level
#'   at a time.
causatr_contrast_summary <- function(res) {
  est <- res$estimates
  ct <- res$contrasts
  # `estimates` carries per-intervention marginal means; `contrasts`
  # carries the pairwise rows. Pull the scalar entries out in a
  # consistent order so the caller can use `$mu_1` etc.
  mu_1 <- est$estimate[est$intervention == "a1"]
  mu_0 <- est$estimate[est$intervention == "a0"]
  list(
    mu_0 = mu_0,
    mu_1 = mu_1,
    contrast = ct$estimate[1],
    se_contrast = ct$se[1]
  )
}
