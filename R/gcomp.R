#' Dispatch g-computation to the point-treatment or longitudinal (ICE) engine
#'
#' @param data data.table from `prepare_data()`.
#' @param outcome Character. Outcome column name.
#' @param treatment Character scalar or vector. Treatment column name(s).
#' @param confounders One-sided formula of baseline confounders.
#' @param confounders_tv One-sided formula of time-varying confounders or `NULL`.
#' @param family Family object or character string for the outcome model.
#' @param estimand Character. `"ATE"`, `"ATT"`, or `"ATC"`.
#' @param type Character. `"point"` or `"longitudinal"`.
#' @param history Positive integer or `Inf`. Markov order for ICE.
#' @param censoring Character or `NULL`. Name of the censoring indicator column.
#' @param weights Numeric vector or `NULL`. Pre-computed observation weights.
#' @param model_fn Function. The fitting function to use, e.g. `stats::glm`,
#'   `mgcv::gam`, `MASS::glm.nb`. Must accept
#'   `(formula, data, family, weights, ...)`.
#' @param call The original `causat()` call (for error messages).
#' @param ... Passed to `model_fn`.
#'
#' @return A `causatr_fit` object.
#'
#' @noRd
fit_gcomp <- function(
  data,
  outcome,
  treatment,
  confounders,
  confounders_tv,
  family,
  estimand,
  type,
  history,
  censoring,
  weights,
  model_fn,
  id,
  time,
  call,
  ...
) {
  if (type == "longitudinal") {
    fit_ice(
      data,
      outcome,
      treatment,
      confounders,
      confounders_tv,
      family,
      estimand,
      history,
      censoring,
      weights,
      model_fn,
      id,
      time,
      call,
      ...
    )
  } else {
    fit_gcomp_point(
      data,
      outcome,
      treatment,
      confounders,
      family,
      estimand,
      censoring,
      weights,
      model_fn,
      call,
      ...
    )
  }
}

#' Fit the outcome model for g-computation with a point treatment
#'
#' @description
#' Implements Step 1 of the parametric g-formula (Hernán & Robins Ch. 13):
#' fit \eqn{E[Y | A, L]} on the uncensored, outcome-observed rows using the
#' user-supplied `model_fn` (default `stats::glm`). The fitted model is
#' stored in the returned `causatr_fit` and used by `contrast()` to predict
#' outcomes under counterfactual interventions.
#'
#' @param data data.table (all rows, including censored / missing outcome).
#' @param outcome Character. Outcome column name.
#' @param treatment Character scalar or vector. Treatment column name(s).
#' @param confounders One-sided formula of baseline confounders.
#' @param family Family object or character string for the outcome model.
#' @param estimand Character. `"ATE"`, `"ATT"`, or `"ATC"`.
#' @param censoring Character or `NULL`. Name of the censoring indicator column
#'   (1 = censored, 0 = uncensored). Model is fit only where `censoring == 0`.
#' @param weights Numeric vector or `NULL`. Observation weights for the full
#'   dataset (subsetted internally to match fitting rows).
#' @param model_fn Function with signature `function(formula, data, family,
#'   weights, ...)`. Defaults to `stats::glm`; pass `mgcv::gam` for GAMs,
#'   `MASS::glm.nb` for negative-binomial, etc.
#' @param call The original `causat()` call.
#' @param ... Passed to `model_fn`.
#'
#' @return A `causatr_fit` object with:
#'   - `model`: the fitted model object.
#'   - `data`: the full (unfitted) dataset for use in `contrast()`.
#'   - `details$fit_rows`: logical vector (length `nrow(data)`) flagging which
#'     rows were used to fit the model.  Needed by `variance_sandwich()` to
#'     map the `n_fit`-length score matrix back to the `n`-length dataset.
#'
#' @noRd
fit_gcomp_point <- function(
  data,
  outcome,
  treatment,
  confounders,
  family,
  estimand,
  censoring,
  weights,
  model_fn,
  call,
  ...
) {
  # Build outcome model formula: Y ~ A [+ A2 ...] + confounder_terms.
  # `term.labels` extracts the RHS terms, preserving interactions and
  # transformations the user wrote in `confounders` (e.g. ns(age, 4), L1*L2).
  confounder_terms <- attr(stats::terms(confounders), "term.labels")
  rhs <- c(treatment, confounder_terms)
  model_formula <- stats::reformulate(rhs, response = outcome)

  # Identify fitting rows: uncensored (C == 0) AND non-missing outcome.
  # contrast() will predict for ALL rows — the g-formula standardises over
  # the full target population regardless of censoring status.
  fit_rows <- is_uncensored(data, censoring) & !is.na(data[[outcome]])
  fit_data <- data[fit_rows]

  # Subset the user-supplied observation weights to match the fitting rows.
  model_weights <- if (!is.null(weights)) weights[fit_rows] else NULL

  # Fit E[Y | A, L] using the caller-supplied fitting function.
  model <- model_fn(
    model_formula,
    data = fit_data,
    family = family,
    weights = model_weights,
    ...
  )

  new_causatr_fit(
    model = model,
    data = data,
    treatment = treatment,
    outcome = outcome,
    confounders = confounders,
    confounders_tv = NULL,
    family = family,
    method = "gcomp",
    type = "point",
    estimand = estimand,
    id = NULL,
    time = NULL,
    censoring = censoring,
    history = 1L,
    numerator = NULL,
    weights_obj = NULL,
    match_obj = NULL,
    call = call,
    details = list(
      fit_rows = fit_rows,
      n_fit = sum(fit_rows),
      n_total = nrow(data),
      model_fn = model_fn
    )
  )
}
