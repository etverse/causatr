#' Construct a `causatr_fit` object
#'
#' @param model Fitted model object (glm, gam, glm_weightit, etc.) or `NULL`
#'   (ICE defers fitting to `contrast()`).
#' @param data data.table of the full dataset.
#' @param treatment Character treatment column name(s).
#' @param outcome Character outcome column name.
#' @param confounders One-sided formula of baseline confounders (soft-deprecated
#'   fallback; prefer the per-component formulas).
#' @param confounders_tv One-sided formula of time-varying confounders or `NULL`
#'   (soft-deprecated fallback; prefer the per-component formulas).
#' @param confounders_outcome One-sided formula of outcome-model confounders, or
#'   `NULL` to fall back to \code{confounders}.
#' @param confounders_treatment One-sided formula of treatment-model confounders,
#'   or `NULL` to fall back to \code{confounders}.
#' @param confounders_censoring One-sided formula of censoring-model confounders,
#'   or `NULL` to fall back to \code{confounders}.
#' @param confounders_sampling One-sided formula of sampling-model confounders,
#'   or `NULL` to fall back to \code{confounders}.
#' @param confounders_tv_outcome One-sided formula of time-varying outcome-model
#'   confounders, or `NULL` to fall back to \code{confounders_tv}.
#' @param confounders_tv_treatment One-sided formula of time-varying
#'   treatment-model confounders, or `NULL` to fall back to
#'   \code{confounders_tv}.
#' @param family Character or family object describing the outcome distribution.
#' @param estimator Character causal estimator (`"gcomp"`, `"ipw"`, `"matching"`).
#' @param type `"point"` or `"longitudinal"`.
#' @param estimand Character estimand (`"ATE"`, `"ATT"`, `"ATC"`).
#' @param id Character ID column name or `NULL`.
#' @param time Character time column name or `NULL`.
#' @param censoring Character censoring column name or `NULL`.
#' @param history Integer Markov order for longitudinal.
#' @param numerator One-sided formula for stabilised weights or `NULL`.
#' @param weights_obj Legacy slot retained for `causatr_fit` backward
#'   compatibility; always `NULL` now that IPW is self-contained.
#' @param match_obj A `matchit` object (matching) or `NULL`.
#' @param call The original `causat()` call environment.
#' @param details Named list of estimator-specific metadata.
#' @param target Character column name of the sampling indicator (S = 1 study,
#'   S = 0 target) or `NULL` when transportability is inactive.
#' @param target_subset `"target"` or `"all"` — which rows form the target
#'   population. `NULL` when `target` is `NULL`.
#' @return A list with class `"causatr_fit"`.
#' @noRd
new_causatr_fit <- function(
  model,
  data,
  treatment,
  outcome,
  confounders,
  confounders_tv,
  confounders_outcome = NULL,
  confounders_treatment = NULL,
  confounders_censoring = NULL,
  confounders_sampling = NULL,
  confounders_tv_outcome = NULL,
  confounders_tv_treatment = NULL,
  family,
  estimator,
  type,
  estimand,
  id,
  time,
  censoring,
  history,
  numerator,
  weights_obj,
  match_obj,
  call,
  details,
  target = NULL,
  target_subset = NULL
) {
  structure(
    list(
      model = model,
      data = data,
      treatment = treatment,
      outcome = outcome,
      confounders = confounders,
      confounders_tv = confounders_tv,
      confounders_outcome = confounders_outcome,
      confounders_treatment = confounders_treatment,
      confounders_censoring = confounders_censoring,
      confounders_sampling = confounders_sampling,
      confounders_tv_outcome = confounders_tv_outcome,
      confounders_tv_treatment = confounders_tv_treatment,
      family = family,
      estimator = estimator,
      type = type,
      estimand = estimand,
      id = id,
      time = time,
      censoring = censoring,
      history = history,
      numerator = numerator,
      weights_obj = weights_obj,
      match_obj = match_obj,
      call = call,
      details = details,
      target = target,
      target_subset = target_subset
    ),
    class = "causatr_fit"
  )
}

#' Construct a `causatr_result` object
#'
#' @param estimates data.table of intervention-specific marginal means.
#' @param contrasts data.table of pairwise contrasts with SEs and CIs.
#' @param type Character contrast type (`"difference"`, `"ratio"`, `"or"`).
#' @param estimand Character estimand used.
#' @param ci_method Character CI method (`"sandwich"` or `"bootstrap"`).
#' @param reference Character name of the reference intervention or `NULL`.
#' @param interventions Named list of `causatr_intervention` objects.
#' @param n Integer sample size used for estimation.
#' @param estimator Character causal estimator.
#' @param vcov Variance-covariance matrix of marginal means.
#' @param call The original `contrast()` call environment.
#' @return A list with class `"causatr_result"`.
#' @noRd
new_causatr_result <- function(
  estimates,
  contrasts,
  type,
  estimand,
  ci_method,
  reference,
  interventions,
  n,
  estimator,
  family,
  fit_type,
  vcov,
  boot_t = NULL,
  boot_info = NULL,
  call
) {
  structure(
    list(
      estimates = estimates,
      contrasts = contrasts,
      type = type,
      estimand = estimand,
      ci_method = ci_method,
      reference = reference,
      interventions = interventions,
      n = n,
      estimator = estimator,
      family = family,
      fit_type = fit_type,
      vcov = vcov,
      boot_t = boot_t,
      # NULL when ci_method = "sandwich"; otherwise a 3-element list of
      # `n_requested`, `n_ok`, `n_fail` carried up from
      # `process_boot_results()` so `print.causatr_result()` and
      # downstream consumers can surface bootstrap failure rates without
      # re-deriving them from `boot_t`.
      boot_info = boot_info,
      call = call
    ),
    class = "causatr_result"
  )
}

#' Construct a `causatr_diag` object
#'
#' @description
#' Builds the nested per-intervention diagnostic container. Each entry in
#' `per_intervention` is a named list with `positivity`, `balance`, and
#' `weights` slots (any of which may be `NULL` when the underlying estimator
#' or treatment family does not support that diagnostic). The top-level
#' `positivity` / `balance` / `weights` slots are populated from the first
#' panel for backward compatibility with the flat shape that downstream
#' callers (`print.causatr_diag()`, `plot.causatr_diag()`, existing tests)
#' originally consumed.
#'
#' @param per_intervention Named list, one entry per intervention key, each
#'   itself a list with slots `positivity`, `balance`, `weights`.
#' @param match_quality data.table of match quality metrics or `NULL`. Lives
#'   at the top level rather than per-intervention because matching is done
#'   once at fit time and is intervention-agnostic.
#' @param estimator Character causal estimator.
#' @param fit_info Named list with summary metadata about the fit
#'   (`treatment_type`, `estimand`, `type`, `has_em`). Used by the print
#'   method to render the panel header.
#' @param fit The original `causatr_fit` (stored for `plot()` method).
#' @return A list with class `"causatr_diag"`.
#' @noRd
new_causatr_diag <- function(
  per_intervention,
  match_quality = NULL,
  estimator,
  fit_info = list(),
  fit = NULL
) {
  if (!is.list(per_intervention) || length(per_intervention) == 0L) {
    rlang::abort(
      "Internal error: `per_intervention` must be a non-empty named list."
    )
  }
  if (
    is.null(names(per_intervention)) || any(!nzchar(names(per_intervention)))
  ) {
    rlang::abort(
      "Internal error: every panel in `per_intervention` must be named."
    )
  }
  # First panel feeds the backward-compat top-level slots. The flat
  # access pattern (`diag$positivity`, `diag$balance`, `diag$weights`)
  # was the public API before the per-intervention redesign; preserving it lets every
  # existing test, print-tests, and downstream user expression keep
  # working unchanged when no `interventions =` argument was passed.
  first <- per_intervention[[1L]]
  structure(
    list(
      per_intervention = per_intervention,
      interventions = names(per_intervention),
      positivity = first$positivity,
      balance = first$balance,
      weights = first$weights,
      censoring = first$censoring,
      sampling = first$sampling,
      match_quality = match_quality,
      estimator = estimator,
      fit_info = fit_info,
      fit = fit
    ),
    class = "causatr_diag"
  )
}

#' Construct a `causatr_intervention` object
#'
#' @param type Character intervention type (e.g. `"static"`, `"shift"`).
#' @param params Named list of intervention parameters.
#' @return A list with class `"causatr_intervention"`.
#' @noRd
new_causatr_intervention <- function(type, params) {
  structure(
    c(list(type = type), params),
    class = "causatr_intervention"
  )
}
