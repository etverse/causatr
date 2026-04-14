#' Fit an IPW model for causal estimation (point treatment)
#'
#' @description
#' Implements inverse probability weighting (Hernán & Robins Ch. 12) by
#' delegating weight estimation to `WeightIt::weightit()` and fitting a
#' weighted outcome model (marginal structural model) via
#' `WeightIt::glm_weightit()`.
#'
#' ## Algorithm
#'
#' 1. Estimate propensity-score weights via `WeightIt::weightit()` using the
#'    treatment model `A ~ confounders`.
#' 2. Fit the weighted outcome model (saturated MSM for binary treatment:
#'    `Y ~ A`) via `WeightIt::glm_weightit()`, which provides proper
#'    M-estimation sandwich SEs that account for weight-estimation uncertainty.
#' 3. Store the `weightit` and `glm_weightit` objects in the returned
#'    `causatr_fit` for use by `contrast()`.
#'
#' @param data data.table from `prepare_data()`.
#' @param outcome Character. Outcome column name.
#' @param treatment Character scalar or vector. Treatment column name(s).
#' @param confounders One-sided formula of baseline confounders.
#' @param confounders_tv One-sided formula of time-varying confounders or `NULL`.
#' @param estimand Character. `"ATE"`, `"ATT"`, or `"ATC"`.
#' @param type Character. `"point"` or `"longitudinal"`.
#' @param history Positive integer or `Inf`. Markov order for longitudinal.
#' @param numerator One-sided formula or `NULL`. Numerator formula for
#'   stabilized weights in longitudinal models.
#' @param weights Numeric vector or `NULL`. External observation weights
#'   (e.g. survey weights), multiplied with the estimated IPW weights.
#' @param call The original `causat()` call.
#' @param ... Passed to `WeightIt::weightit()` (e.g. `method = "glm"`).
#'
#' @return A `causatr_fit` object with `weights_obj` (the `weightit` object)
#'   and `model` (the `glm_weightit` fit).
#'
#' @noRd
fit_ipw <- function(
  data,
  outcome,
  treatment,
  confounders,
  confounders_tv,
  family,
  estimand,
  type,
  history,
  numerator,
  weights,
  call,
  ...
) {
  if (type == "longitudinal") {
    rlang::abort(
      "Longitudinal IPW (via WeightIt::weightitMSM) is not yet implemented.",
      .call = FALSE
    )
  }

  # Reject A-touching terms in `confounders` before building the PS
  # formula. IPW's saturated MSM `Y ~ A` has nowhere to put an
  # `A:modifier` interaction, so any such term would be silently
  # dropped from the effect estimate. Abort early with a Phase-8
  # pointer rather than returning a wrong answer. See
  # `check_confounders_no_treatment()` in `R/utils.R`.
  check_confounders_no_treatment(confounders, treatment, estimator = "ipw")

  # Build the treatment model formula: A ~ confounders.
  # This is used by WeightIt to estimate propensity scores.
  ps_formula <- build_ps_formula(confounders, treatment)

  # Fit rows: exclude missing outcomes for the MSM.
  fit_rows <- get_fit_rows(data, outcome)
  fit_data <- data[fit_rows]

  # Step 1: Estimate propensity-score weights via WeightIt.
  # WeightIt computes stabilized weights by default.
  # Additional arguments (e.g. method = "glm", "gbm", "bart") go via `...`.
  w <- WeightIt::weightit(
    ps_formula,
    data = fit_data,
    estimand = estimand,
    ...
  )

  # `get_fit_rows()` only excludes rows with missing outcome, but
  # WeightIt drops rows with missing PS-formula covariates as well.
  # If those row sets differ, `w$weights` is shorter than `fit_rows`
  # and the external-weight multiply on the next line would either
  # recycle silently or misalign weights to the wrong individuals.
  # Require the row counts to match and ask the user to clean or
  # impute the data first.
  if (length(w$weights) != sum(fit_rows)) {
    rlang::abort(
      paste0(
        "WeightIt returned ",
        length(w$weights),
        " weights but causatr selected ",
        sum(fit_rows),
        " fitting rows (rows with non-missing outcome). This usually ",
        "means a confounder column has missing values that WeightIt ",
        "dropped. Drop those rows manually before calling `causat()` ",
        "or impute them so the row counts agree."
      ),
      .call = FALSE
    )
  }

  # If external weights are provided (e.g. survey weights), multiply them
  # with the estimated IPW weights.
  if (!is.null(weights)) {
    w$weights <- w$weights * weights[fit_rows]
  }

  # Mparts guardrail: WeightIt methods that do not support the
  # M-estimation correction (gbm, super, bart, optweight, energy, npcbps,
  # cfd) silently treat weights as fixed inside glm_weightit(). Sandwich SEs
  # then ignore propensity-estimation uncertainty. Warn at fit time so
  # users learn about the limitation before they call contrast() or
  # diagnose().
  if (is.null(attr(w, "Mparts"))) {
    rlang::warn(
      paste0(
        "WeightIt method '",
        if (is.null(w$method)) "?" else as.character(w$method),
        "' does not implement the M-estimation correction (no `Mparts` ",
        "attribute). Sandwich SEs from `ci_method = 'sandwich'` will treat ",
        "the weights as fixed and may underestimate the variance. Use ",
        "`ci_method = 'bootstrap'` for valid inference, or switch to a ",
        "method that supports Mparts (`glm`, `cbps`, `ipt`, `ebal`)."
      )
    )
  }

  # Step 2: Fit the weighted outcome model (marginal structural model).
  # For binary treatment the MSM is saturated: Y ~ A.
  # glm_weightit() provides M-estimation sandwich SEs that correctly
  # account for weight-estimation uncertainty.
  msm_formula <- stats::reformulate(treatment, response = outcome)

  msm_fit <- WeightIt::glm_weightit(
    msm_formula,
    data = fit_data,
    weightit = w,
    family = resolve_family(family)
  )

  new_causatr_fit(
    model = msm_fit,
    data = data,
    treatment = treatment,
    outcome = outcome,
    confounders = confounders,
    confounders_tv = confounders_tv,
    family = family,
    estimator = "ipw",
    type = "point",
    estimand = estimand,
    id = NULL,
    time = NULL,
    censoring = NULL,
    history = history,
    numerator = numerator,
    weights_obj = w,
    match_obj = NULL,
    call = call,
    details = list(
      fit_rows = fit_rows,
      n_fit = sum(fit_rows),
      n_total = nrow(data),
      weights = weights
    )
  )
}
