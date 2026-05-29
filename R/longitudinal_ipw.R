#' Fit a longitudinal IPW model (per-period treatment density chain)
#'
#' @description
#' Longitudinal companion to `fit_ipw()` for time-varying treatments.
#' Fits the sequential factorisation
#' \deqn{f(A_0, A_1, \ldots, A_K \mid L) =
#'       \prod_{k=0}^{K} f(A_k \mid \bar A_{k-1}, \bar L_k)}
#' as `K + 1` independent per-period density models, one per
#' `time_points[k]`. Each per-period model is fit on the subset of
#' rows at `time == time_points[k]` with a formula that conditions on
#' baseline confounders, current-period time-varying confounders, and
#' lagged treatment / lagged time-varying confounders up to `history`.
#' The lag columns are already materialised by `prepare_data()` /
#' `create_lag_vars()`, so the per-period formula just references them
#' by name.
#'
#' Per-period density models are fit independently via
#' `fit_treatment_model()`, which makes the bread of the stacked
#' propensity system block-diagonal -- the variance engine sums `K`
#' single-model `apply_model_correction()` calls (same shape as the
#' multivariate IPW propensity correction).
#'
#' ## What is rejected at this fit boundary
#'
#' Longitudinal IPW supports:
#'
#' - univariate and multivariate treatment (binary / continuous / count
#'   via `propensity_family`); multivariate factorises the joint density
#'   over both the time and component axes.
#' - `static()` / `shift()` / `scale_by()` / `dynamic()` interventions
#' - intercept-only Hajek MSM, plus `Y ~ 1 + modifier` under
#'   baseline effect modification (univariate only)
#' - unstabilized and `stabilize = "marginal"` weights, including the
#'   multivariate per-period per-component numerator chain
#'
#' Effect modification combined with multivariate treatment and
#' `ipsi()` are rejected upfront with classed errors. `threshold()`
#' continues to be rejected per period by
#' `check_intervention_family_compat()`.
#'
#' @param data data.table from `prepare_data()` with lag columns and
#'   sorted by `(id, time)`.
#' @param outcome Character. Outcome column name (observed at the
#'   final time point only; earlier periods may carry NA).
#' @param treatment Character scalar. Treatment column name.
#' @param confounders One-sided formula of baseline (time-invariant)
#'   confounders.
#' @param confounders_tv One-sided formula of time-varying confounders
#'   or `NULL`.
#' @param family Character or family object for the final-period MSM.
#' @param estimand Character. Must be `"ATE"` for longitudinal data
#'   (`check_estimand_trt_compat()` already enforces this).
#' @param history Non-negative integer or `Inf`. Markov lag order for
#'   the per-period propensity formula (`0` = no lags).
#' @param numerator One-sided formula or `NULL`. Reserved for chunk
#'   10b (stabilized weights); rejected here when non-`NULL`.
#' @param weights Numeric vector or `NULL`. External observation
#'   weights (e.g. survey weights, IPCW), aligned to `data` rows.
#' @param model_fn Function. Outcome MSM fitter; default `stats::glm`.
#' @param propensity_model_fn Function or `NULL`. Per-period
#'   propensity fitter. `NULL` reuses `model_fn`.
#' @param propensity_family Character or `NULL`. Explicit propensity
#'   family override (`"poisson"` / `"negbin"`); broadcast across all
#'   periods.
#' @param id Character. Individual ID column.
#' @param time Character. Time column.
#' @param call The original `causat()` call.
#' @param ... Forwarded to the per-period propensity fitter.
#'
#' @return A `causatr_fit` of `estimator = "ipw"`, `type =
#'   "longitudinal"` carrying:
#'   \describe{
#'     \item{`details$treatment_models_by_time`}{Length-`K+1` named
#'       list of `causatr_treatment_model` objects, named by string
#'       representation of `time_points[k]`.}
#'     \item{`details$time_points`}{Sorted unique time values.}
#'     \item{`details$n_times`}{Number of distinct time points.}
#'     \item{`details$max_lag`}{Resolved Markov order.}
#'     \item{`details$baseline_terms` / `tv_vars`}{Formula bookkeeping
#'       used to rebuild per-period propensity formulas at variance /
#'       bootstrap time.}
#'   }
#'
#' @noRd
fit_longitudinal_ipw <- function(
  data,
  outcome,
  treatment,
  confounders,
  confounders_tv,
  family,
  estimand,
  history,
  numerator,
  weights,
  model_fn,
  propensity_model_fn,
  propensity_family = NULL,
  stabilize = "none",
  id,
  time,
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
  is_multivariate <- length(treatment) > 1L

  # `numerator =` is the user-facing knob for a custom numerator
  # formula. Under `stabilize = "none"` it has no meaning -- the
  # numerator is the same density as the denominator -- so reject the
  # combination upfront.
  if (!is.null(numerator)) {
    if (stabilize == "none") {
      rlang::abort(
        c(
          "`numerator =` requires `stabilize = 'marginal'`.",
          i = "Drop `numerator` or set `stabilize = 'marginal'`."
        ),
        class = "causatr_longitudinal_numerator_without_stabilize"
      )
    }
    check_formula(numerator)
  }

  # Effect modification under longitudinal IPW reuses the same parser
  # the point IPW path uses (`parse_effect_mod()`). EM terms appear in
  # `confounders` as `A:modifier`; the per-period propensity formulas
  # strip every treatment-touching term via `em_info$confounder_terms`
  # (so the propensity is `A ~ baseline + tv + lags` with no `A:`
  # interactions), and the final-period MSM expands from `Y ~ 1` to
  # `Y ~ 1 + modifier` via `build_ipw_msm_formula()`.
  #
  # Known limitation (Robins 2000): the modifier MUST be a
  # **baseline** covariate. A time-varying modifier in an MSM
  # conditions on a post-treatment variable, biasing the estimand
  # via mediator + collider paths. We don't enforce baseline-ness at
  # runtime because time-varying status isn't inferable from the data;
  # the constraint is doc-level.
  # When per-component confounders are supplied, EM terms live in the
  # outcome confounders; validate treatment confounders separately.
  em_info <- check_confounders_treatment(
    confounders_outcome %||% confounders,
    treatment,
    estimator = "ipw"
  )
  if (!is.null(confounders_outcome)) {
    check_confounders_treatment(confounders, treatment, estimator = "ipw")
  }

  # EM + multivariate longitudinal not yet implemented.
  if (is_multivariate && em_info$has_em) {
    rlang::abort(
      c(
        "Effect modification is not yet supported for multivariate longitudinal IPW.",
        i = "Remove `A:modifier` terms from confounders, or switch to a single treatment."
      ),
      class = "causatr_longitudinal_mv_em_pending"
    )
  }

  # Sorted unique time points. Per-period propensity models index into
  # `time_points[k]`; lag columns (`lag1_A`, ...) created by
  # `prepare_data()` are aligned to this ordering because
  # `create_lag_vars()` keys on `(id, time)` before shifting.
  # Sorting here guarantees `treatment_models_by_time` is keyed in
  # chronological order even when the data arrives unsorted.
  time_points <- sort(unique(data[[time]]))
  n_times <- length(time_points)

  if (n_times < 2L) {
    rlang::abort(
      c(
        paste0(
          "`type = 'longitudinal'` requires at least 2 distinct time points (got ",
          n_times,
          ")."
        ),
        i = "Use `type = 'point'` if the data really has only one period per individual."
      ),
      class = "causatr_longitudinal_too_few_times"
    )
  }

  # Resolve the Markov order. Mirrors `fit_ice()`'s convention.
  max_lag <- if (is.infinite(history)) {
    n_times - 1L
  } else {
    as.integer(history)
  }

  # Baseline confounders enter every per-period propensity. Strip the
  # treatment-touching (EM) terms for parity with point IPW
  # `build_ps_formula()` -- here EM is rejected upfront, but using
  # `em_info$confounder_terms` keeps the call self-consistent if the
  # rejection is later relaxed.
  baseline_terms <- em_info$confounder_terms
  if (length(baseline_terms) == 0L) {
    baseline_terms <- character(0L)
  }

  # Time-varying confounder names (plain identifiers). The per-period
  # formula references the current-time `L` and the lag columns
  # `lag1_L`, `lag2_L`, ... by name.
  tv_vars <- if (!is.null(confounders_tv)) {
    all.vars(confounders_tv)
  } else {
    character(0)
  }

  # Capture user `...` so the bootstrap can replay the propensity
  # fitter configuration on every replicate via `refit_ipw()`.
  dots <- list(...)

  # Resolve the propensity fitter. For univariate longitudinal IPW the
  # default fitter is `model_fn` (typically `stats::glm`); a user can
  # override with `propensity_model_fn` (e.g. `mgcv::gam`). NB / count
  # families opt in via `propensity_family`. For multivariate,
  # `fit_treatment_models()` handles per-component dispatch internally.
  if (!is_multivariate) {
    trt_family_first <- detect_treatment_family(data[[treatment]])
    prop_model_fn <- if (!is.null(propensity_model_fn)) {
      propensity_model_fn
    } else if (trt_family_first == "categorical") {
      nnet::multinom
    } else if (identical(propensity_family, "negbin")) {
      check_pkg("MASS")
      MASS::glm.nb
    } else {
      model_fn
    }
  } else {
    prop_model_fn <- propensity_model_fn %||% model_fn
  }

  # Per-period propensity fits.
  #
  # Univariate: `fit_treatment_model()` per period.
  # Multivariate: `fit_treatment_models()` per period — returns a
  #   `causatr_treatment_models` list of K component models with
  #   sequential conditioning (A_k ~ A_{1..k-1} + confounders).
  #
  # `treatment_models_by_time` is named by the string representation
  # of each time point (e.g. "0", "1", "2") so callers can retrieve
  # the k-th period model by name without positional indexing.
  treatment_models_by_time <- vector("list", n_times)
  fit_data_by_time <- vector("list", n_times)
  per_period_formula <- vector("list", n_times)

  for (k in seq_len(n_times)) {
    rows_k <- data[[time]] == time_points[k]
    data_k <- data[rows_k]

    available_lags <- min(k - 1L, max_lag)

    # Build the confounders formula with baseline + TV + lags.
    # For multivariate, `build_longitudinal_ps_formula()` creates a
    # formula with all treatment lags in the RHS; the response is
    # stripped by `remove_response()` and `fit_treatment_models()`
    # builds per-component formulas with sequential conditioning.
    ps_formula <- build_longitudinal_ps_formula(
      treatment = treatment,
      baseline_terms = baseline_terms,
      tv_vars = tv_vars,
      available_lags = available_lags,
      data_at_time = data_k
    )

    if (is_multivariate) {
      # `fit_treatment_models()` handles sequential conditioning:
      # A_k ~ A_{1..k-1} + confounders, stripping EM terms.
      # The confounders formula includes baseline + TV + lag terms.
      tm_args <- list(
        data = data_k,
        treatment = treatment,
        confounders = remove_response(ps_formula),
        model_fn = prop_model_fn,
        propensity_family = propensity_family
      )
      if (!is.null(weights)) {
        tm_args$weights <- weights[rows_k]
      }
      tm_k <- do.call(fit_treatment_models, c(tm_args, dots))
    } else {
      tm_args <- list(
        data = data_k,
        treatment = treatment,
        confounders = remove_response(ps_formula),
        model_fn = prop_model_fn,
        propensity_family = propensity_family
      )
      if (!is.null(weights)) {
        tm_args$weights <- weights[rows_k]
      }
      tm_k <- do.call(fit_treatment_model, c(tm_args, dots))
    }

    treatment_models_by_time[[k]] <- tm_k
    fit_data_by_time[[k]] <- data_k
    per_period_formula[[k]] <- ps_formula
  }
  names(treatment_models_by_time) <- as.character(time_points)
  names(fit_data_by_time) <- as.character(time_points)
  names(per_period_formula) <- as.character(time_points)

  # Stabilized weights: per-period numerator models g_k(A_k | A_{k-1},
  # ..., A_0) with L stripped from the conditioning. The full
  # density-ratio weight per period becomes
  #   w_k = g_k(d(A_k, L_k) | A_{1..k-1}^obs) * |Jac d^{-1}|
  #         / f_k(A_k | A_{1..k-1}^obs, L_{1..k}^obs).
  # Numerator parameters gamma are held fixed under the variance
  # engine's `numDeriv` perturbation of the denominator alpha (same
  # nuisance-fixed convention as multivariate IPW and
  # sigma / theta for Gaussian / negbin densities). Bootstrap refits
  # both gamma and alpha, capturing full variance.
  #
  # Custom `numerator` formula overrides the default "drop L" rule:
  # `numerator = ~ baseline` keeps only baseline confounders;
  # `numerator = ~ 1` reduces to a marginal-by-period density.
  numerator_models_by_time <- NULL
  if (stabilize == "marginal") {
    numerator_models_by_time <- vector("list", n_times)
    for (k in seq_len(n_times)) {
      data_k <- fit_data_by_time[[k]]
      available_lags <- min(k - 1L, max_lag)

      # Per-period numerator formula. Default behaviour (no `numerator`
      # argument): drop L (baseline + tv + lag_L); keep treatment lags
      # only. With a user-supplied `numerator` formula we honour it
      # verbatim as the conditioning RHS, plus the available treatment
      # lags so the per-period factorisation g_k(A_k | A_{1..k-1}, V)
      # still respects the chain rule. `V` is the user-supplied set.
      num_ps_formula <- build_longitudinal_numerator_ps_formula(
        treatment = treatment,
        numerator = numerator,
        available_lags = available_lags,
        data_at_time = data_k
      )

      num_args <- list(
        data = data_k,
        treatment = treatment,
        confounders = remove_response(num_ps_formula),
        model_fn = prop_model_fn,
        propensity_family = propensity_family
      )
      if (!is.null(weights)) {
        num_args$weights <- weights[data[[time]] == time_points[k]]
      }
      if (is_multivariate) {
        # Per-component numerator chain. `build_longitudinal_numerator_ps_formula()`
        # vectorises over `treatment`, so its RHS already carries the
        # prior-period treatment lags for every component (Ā_{t-1}) plus
        # any user `numerator = ~ V`. Passing it as `confounders` to
        # `fit_treatment_models()` adds the within-period prior components
        # A_{t,1..k-1} via the chain rule, giving the numerator
        # g_{t,k}(A_{t,k} | A_{t,1..k-1}, Ā_{t-1}, V) with L dropped.
        # We fit these as a plain (`stabilize = "none"`) component list:
        # the numerator factorisation needs only the denominator-style
        # per-component models, not a second nested numerator.
        numerator_models_by_time[[k]] <- do.call(
          fit_treatment_models,
          c(num_args, dots)
        )
      } else {
        numerator_models_by_time[[k]] <- do.call(
          fit_treatment_model,
          c(num_args, dots)
        )
      }
    }
    names(numerator_models_by_time) <- as.character(time_points)
  }

  # Final-period rows define the MSM scope and the marginal-mean
  # average. The Hajek intercept of `Y ~ 1` weighted by the cumulative
  # product weight reads `E[Y^d]` directly off the (sole) coefficient.
  # Only the final period carries a non-missing outcome in a typical
  # longitudinal study; earlier periods have `NA` in the outcome column.
  last_time <- time_points[n_times]
  rows_final <- data[[time]] == last_time
  fit_rows_final <- rows_final & !is.na(data[[outcome]])
  # Transport: restrict the final-period MSM to study rows (S=1). Target
  # rows have no observed outcome; the sampling odds weights reweight the
  # study population to the target population within the weighted MSM.
  if (!is.null(target)) {
    fit_rows_final <- fit_rows_final & data[[target]] == 1L
  }

  if (sum(fit_rows_final) == 0L) {
    rlang::abort(
      c(
        paste0(
          "No final-period (`time == ",
          last_time,
          "`) rows have non-missing outcome -- the longitudinal IPW MSM has nothing to fit."
        ),
        i = "Check that the outcome column is populated at the last time point."
      ),
      class = "causatr_longitudinal_no_final_outcome"
    )
  }

  # Placeholder MSM model. Same role as `fit_ipw()`'s placeholder
  # `Y ~ A` -- gives `print()` / `summary()` something to display
  # before any contrast is requested. The real per-intervention MSM
  # is refit inside `compute_ipw_contrast_longitudinal()`.
  # `Y ~ 1` (unweighted) rather than `Y ~ A` because in longitudinal
  # data the treatment varies by period, so `Y ~ A` would be ambiguous.
  fam_obj <- resolve_family(family)
  placeholder_args <- list(
    formula = stats::reformulate("1", response = outcome),
    data = data[fit_rows_final],
    family = fam_obj
  )
  if (!is.null(weights)) {
    placeholder_args$weights <- weights[fit_rows_final]
  }
  placeholder_model <- do.call(model_fn, placeholder_args)

  new_causatr_fit(
    model = placeholder_model,
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
    estimator = "ipw",
    type = "longitudinal",
    estimand = estimand,
    id = id,
    time = time,
    censoring = NULL,
    history = history,
    numerator = numerator,
    weights_obj = NULL,
    match_obj = NULL,
    call = call,
    details = list(
      time_points = time_points,
      n_times = n_times,
      max_lag = max_lag,
      baseline_terms = baseline_terms,
      tv_vars = tv_vars,
      em_info = em_info,
      treatment_models_by_time = treatment_models_by_time,
      numerator_models_by_time = numerator_models_by_time,
      fit_data_by_time = fit_data_by_time,
      per_period_formula = per_period_formula,
      fit_rows_final = fit_rows_final,
      n_fit_final = sum(fit_rows_final),
      n_total = nrow(data),
      weights = weights,
      dots = dots,
      model_fn = model_fn,
      propensity_model_fn = prop_model_fn,
      propensity_family = propensity_family,
      stabilize = stabilize,
      is_multivariate = is_multivariate,
      confounders_outcome = confounders_outcome
    )
  )
}
