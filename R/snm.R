#' Fit a Structural Nested Mean Model (SNMM)
#'
#' Internal workhorse called by [causat()] when `estimator = "snm"`.
#' Fits the treatment model \eqn{\hat{E}[A \mid L]} via
#' [fit_treatment_model()] and stores the blip parameterisation derived
#' from effect-modification terms in the confounders formula. No point
#' estimation happens here — the g-estimating equation is solved by
#' [contrast()] (analogous to how ICE defers model fitting).
#'
#' @param data data.table of the (prepared) dataset.
#' @param outcome Character outcome column name.
#' @param treatment Character treatment column name.
#' @param confounders One-sided formula of confounders (resolved:
#'   `confounders_treatment` or the unified `confounders`).
#' @param confounders_tv One-sided formula of time-varying confounders or
#'   `NULL`.
#' @param family Character or family object for the outcome distribution.
#' @param estimand Character estimand (`"ATE"` only for now; ATT/ATC
#'   gated downstream).
#' @param type `"point"` or `"longitudinal"`.
#' @param history Integer Markov order for longitudinal lag expansion.
#' @param weights Numeric vector of external weights or `NULL`.
#' @param propensity_model_fn Model fitting function for the treatment
#'   model (default `stats::glm`).
#' @param propensity_family Family override for the treatment model or
#'   `NULL`.
#' @param id Character ID column name or `NULL`.
#' @param time Character time column name or `NULL`.
#' @param call The original `causat()` call.
#' @param target Character column name of the sampling indicator or `NULL`.
#' @param confounders_outcome One-sided formula for outcome confounders or
#'   `NULL`.
#' @param confounders_outcome_raw Raw user-supplied `confounders_outcome`.
#' @param confounders_treatment_raw Raw user-supplied `confounders_treatment`.
#' @param confounders_censoring_raw Raw user-supplied `confounders_censoring`.
#' @param confounders_sampling_raw Raw user-supplied `confounders_sampling`.
#' @param confounders_tv_outcome_raw Raw user-supplied `confounders_tv_outcome`.
#' @param confounders_tv_treatment_raw Raw user-supplied
#'   `confounders_tv_treatment`.
#' @param ... Extra arguments forwarded to the treatment model fitter.
#' @return A `causatr_fit` with `estimator = "snm"`.
#' @noRd
fit_snm <- function(
  data,
  outcome,
  treatment,
  confounders,
  confounders_tv,
  family,
  estimand,
  type,
  history,
  weights,
  propensity_model_fn,
  propensity_family = NULL,
  id = NULL,
  time = NULL,
  call,
  target = NULL,
  confounders_outcome = NULL,
  confounders_outcome_raw = NULL,
  confounders_treatment_raw = NULL,
  confounders_censoring_raw = NULL,
  confounders_sampling_raw = NULL,
  confounders_tv_outcome_raw = NULL,
  confounders_tv_treatment_raw = NULL,
  ...
) {
  # --- Rejection: multivariate treatment ---
  if (length(treatment) > 1L) {
    rlang::abort(
      c(
        "Multivariate treatment is not supported for `estimator = \"snm\"`.",
        i = paste0(
          "Multivariate SNMMs are out of scope. Use `estimator = \"ipw\"` ",
          "or `\"gcomp\"` for joint interventions on multiple treatments."
        )
      ),
      class = "causatr_snm_multivariate"
    )
  }

  # Parse effect-modification terms from confounders. The EM terms
  # (`A:modifier`) define the blip parameterisation: each interaction
  # contributes a modifier column whose coefficient in the blip
  # function is estimated by the g-estimating equation. Pure confounder
  # terms enter the treatment model only.
  em_conf <- confounders_outcome %||% confounders
  em_info <- parse_effect_mod(em_conf, treatment)

  # Inform when no EM terms are present — the SNMM reduces to a single
  # ATE parameter, which is fine but worth flagging because the main
  # motivation for SNMs is effect-modification identification.
  if (!em_info$has_em) {
    rlang::inform(
      c(
        "No effect-modification terms (e.g. `A:modifier`) found in confounders.",
        i = paste0(
          "The SNMM blip reduces to a single constant-effect parameter ",
          "(equivalent to an ATE). Add `A:modifier` terms to the confounders ",
          "formula to estimate effect modification."
        )
      ),
      class = "causatr_snm_no_em"
    )
  }

  # Build blip specification from parsed EM terms. Each EM term
  # contributes modifier columns; the intercept (`A * psi_0`) is
  # always included.
  blip_spec <- build_blip_spec(em_info, treatment)

  # --- Fit rows: exclude missing outcomes ---
  fit_rows <- get_fit_rows(data, outcome, target = target)
  fit_data <- data[fit_rows]

  # Resolve propensity fitter
  prop_fn <- propensity_model_fn %||% stats::glm

  if (type == "point") {
    # Fit a single treatment model E[A | L]. `fit_treatment_model()`
    # calls `build_ps_formula()` internally, which strips any EM terms
    # from `confounders` so the propensity model is A ~ L (no A:modifier).
    trt_model <- fit_treatment_model(
      data = fit_data,
      treatment = treatment,
      confounders = confounders,
      model_fn = prop_fn,
      propensity_family = propensity_family,
      weights = if (is.null(weights)) NULL else weights[fit_rows],
      ...
    )
  } else {
    # Longitudinal: fit per-period treatment models (deferred to 18d)
    rlang::abort(
      c(
        "Longitudinal SNMMs are not yet supported.",
        i = "Longitudinal g-estimation will be added in a future update."
      ),
      class = "causatr_snm_longitudinal_pending"
    )
  }

  # Capture user dots for bootstrap replay
  dots <- list(...)

  new_causatr_fit(
    model = NULL,
    data = data,
    treatment = treatment,
    outcome = outcome,
    confounders = confounders,
    confounders_tv = confounders_tv,
    confounders_outcome = confounders_outcome_raw,
    confounders_treatment = confounders_treatment_raw,
    confounders_censoring = confounders_censoring_raw,
    confounders_sampling = confounders_sampling_raw,
    confounders_tv_outcome = confounders_tv_outcome_raw,
    confounders_tv_treatment = confounders_tv_treatment_raw,
    family = family,
    estimator = "snm",
    type = type,
    estimand = estimand,
    id = id,
    time = time,
    censoring = NULL,
    history = history,
    numerator = NULL,
    weights_obj = NULL,
    match_obj = NULL,
    call = call,
    details = list(
      treatment_model = trt_model,
      blip_spec = blip_spec,
      blip_type = "linear",
      em_info = em_info,
      fit_rows = fit_rows,
      weights = weights,
      dots = dots,
      propensity_model_fn = prop_fn,
      propensity_family = propensity_family
    ),
    target = target
  )
}


#' Build blip specification from parsed effect-modification terms
#'
#' Translates the `causatr_em_info` object into a blip specification
#' that describes the linear blip function
#' \eqn{\gamma(a, l; \psi) = a \cdot (\psi_0 + \sum_j \psi_j \cdot m_j)}.
#' The returned list enumerates modifier columns and names the blip
#' parameters.
#'
#' @param em_info A `causatr_em_info` from [parse_effect_mod()].
#' @param treatment Character scalar treatment column name.
#' @return A list with:
#'   \describe{
#'     \item{modifier_cols}{Character vector of modifier column names
#'       (empty if no EM terms).}
#'     \item{param_names}{Character vector of blip parameter names,
#'       always starting with `"psi_intercept"`.}
#'     \item{n_params}{Integer number of blip parameters.}
#'   }
#' @noRd
build_blip_spec <- function(em_info, treatment) {
  modifier_cols <- em_info$modifier_vars
  param_names <- "psi_intercept"

  if (length(modifier_cols) > 0L) {
    param_names <- c(param_names, paste0("psi_", modifier_cols))
  }

  list(
    modifier_cols = modifier_cols,
    param_names = param_names,
    n_params = length(param_names)
  )
}
