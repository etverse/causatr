#' Fit a Structural Nested Mean Model (SNMM)
#'
#' Internal workhorse called by [causat()] when `estimator = "snm"`.
#' Fits the treatment model \eqn{\hat{E}[A \mid L]} via
#' `fit_treatment_model()` and stores the blip parameterisation derived
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
#' @param treatment_free One-sided formula or `NULL`. Specifies the
#'   treatment-free outcome model \eqn{E[Y \mid L]} for efficiency
#'   augmentation. When non-`NULL`, the model's predictions are
#'   subtracted from Y before g-estimation, projecting out the
#'   \eqn{L \to Y} association and reducing variance.
#' @param model_fn Model fitting function for the treatment-free model
#'   (default `stats::glm`).
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
  treatment_free = NULL,
  model_fn = stats::glm,
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

  # Default to GLM for the propensity model if no fitter is supplied
  prop_fn <- propensity_model_fn %||% stats::glm

  # Capture user dots early — needed by both point and longitudinal
  # branches for treatment model fitting and bootstrap replay.
  dots <- list(...)

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
    # --- Longitudinal: per-period treatment models ---
    time_points <- sort(unique(data[[time]]))
    n_times <- length(time_points)

    if (n_times < 2L) {
      rlang::abort(
        c(
          paste0(
            "`type = 'longitudinal'` requires at least 2 distinct time ",
            "points (got ",
            n_times,
            ")."
          ),
          i = "Use `type = 'point'` for single-period data."
        ),
        class = "causatr_longitudinal_too_few_times"
      )
    }

    max_lag <- if (is.infinite(history)) {
      n_times - 1L
    } else {
      as.integer(history)
    }

    # Parse baseline confounder terms (strip A:modifier interactions so
    # the propensity model is A ~ L, matching the point-treatment path).
    # When confounders is NULL (user supplied only confounders_tv), there
    # are no baseline terms — the PS formula is built from tv_vars only.
    if (!is.null(confounders)) {
      em_info_trt <- parse_effect_mod(confounders, treatment)
      baseline_terms <- em_info_trt$confounder_terms
      if (length(baseline_terms) == 0L) {
        baseline_terms <- character(0L)
      }
    } else {
      baseline_terms <- character(0L)
    }

    tv_vars <- if (!is.null(confounders_tv)) {
      all.vars(confounders_tv)
    } else {
      character(0)
    }

    treatment_models_by_time <- vector("list", n_times)
    fit_data_by_time <- vector("list", n_times)
    per_period_formula <- vector("list", n_times)

    for (k in seq_len(n_times)) {
      rows_k <- data[[time]] == time_points[k]
      data_k <- data[rows_k]
      available_lags <- min(k - 1L, max_lag)

      ps_formula <- build_longitudinal_ps_formula(
        treatment = treatment,
        baseline_terms = baseline_terms,
        tv_vars = tv_vars,
        available_lags = available_lags,
        data_at_time = data_k
      )

      tm_args <- list(
        data = data_k,
        treatment = treatment,
        confounders = remove_response(ps_formula),
        model_fn = prop_fn,
        propensity_family = propensity_family
      )
      if (!is.null(weights)) {
        tm_args$weights <- weights[rows_k]
      }
      tm_k <- do.call(fit_treatment_model, c(tm_args, dots))

      treatment_models_by_time[[k]] <- tm_k
      fit_data_by_time[[k]] <- data_k
      per_period_formula[[k]] <- ps_formula
    }
    names(treatment_models_by_time) <- as.character(time_points)
    names(fit_data_by_time) <- as.character(time_points)
    names(per_period_formula) <- as.character(time_points)

    trt_model <- list(
      models_by_time = treatment_models_by_time,
      time_points = time_points,
      n_times = n_times
    )
  }

  # --- Treatment-free model specification ---
  # When treatment_free is specified, store the formula. The actual model
  # is fit jointly with the blip parameters in compute_snm_blip_point(),
  # following the joint EE approach of Vansteelandt & Joffe (2014) and
  # DTRreg: (beta, psi) are solved simultaneously from
  #   E[Z_i (Y_i - Z_i beta - gamma(A_i, L_i; psi))] = 0
  #   E[R_i m_i (Y_i - Z_i beta - gamma(A_i, L_i; psi))] = 0
  if (!is.null(treatment_free)) {
    if (!inherits(treatment_free, "formula")) {
      rlang::abort(
        "`treatment_free` must be a one-sided formula (e.g. `~ L`).",
        class = "causatr_snm_bad_treatment_free"
      )
    }
  }

  details <- list(
    treatment_model = trt_model,
    blip_spec = blip_spec,
    blip_type = "linear",
    em_info = em_info,
    fit_rows = fit_rows,
    weights = weights,
    dots = dots,
    propensity_model_fn = prop_fn,
    propensity_family = propensity_family,
    treatment_free = treatment_free,
    model_fn = model_fn
  )

  # Longitudinal-specific details for the variance engine and contrast
  if (type == "longitudinal") {
    details$treatment_models_by_time <- trt_model$models_by_time
    details$fit_data_by_time <- fit_data_by_time
    details$per_period_formula <- per_period_formula
    details$time_points <- trt_model$time_points
    details$n_times <- trt_model$n_times
    details$max_lag <- max_lag
    details$baseline_terms <- baseline_terms
    details$tv_vars <- tv_vars
  }

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
    details = details,
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
#' @param em_info A `causatr_em_info` from `parse_effect_mod()`.
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
  # Extract modifier variables from the EM info (parsed from A:M terms
  # in the confounders formula by parse_effect_mod())
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
