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

  # Build the treatment model formula: A ~ confounders.
  # This is used by WeightIt to estimate propensity scores.
  confounder_terms <- attr(stats::terms(confounders), "term.labels")
  ps_formula <- stats::reformulate(confounder_terms, response = treatment)

  # Fit rows: exclude missing outcomes for the MSM.
  fit_rows <- !is.na(data[[outcome]])
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

  # If external weights are provided (e.g. survey weights), multiply them
  # with the estimated IPW weights.
  if (!is.null(weights)) {
    w$weights <- w$weights * weights[fit_rows]
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
    method = "ipw",
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
      n_total = nrow(data)
    )
  )
}
