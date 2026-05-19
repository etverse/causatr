#' Fit both nuisance models for longitudinal AIPW (ICE-AIPW)
#'
#' @description
#' Composes the ICE outcome-model metadata with the
#' longitudinal IPW per-period propensity models into a
#' single `causatr_fit` that carries both for downstream use in
#' `ice_aipw_iterate()`.
#'
#' Propensity models are fit at this stage (they are intervention-
#' independent). Outcome models are deferred to contrast time (they
#' are intervention-dependent through the ICE backward iteration),
#' following the same convention as `fit_ice()`.
#'
#' @inheritParams fit_aipw
#' @noRd
fit_aipw_longitudinal <- function(
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
  propensity_model_fn,
  propensity_family,
  id,
  time,
  call,
  confounders_treatment = NULL,
  confounders_tv_treatment = NULL,
  confounders_outcome_raw = NULL,
  confounders_treatment_raw = NULL,
  confounders_censoring_raw = NULL,
  confounders_sampling_raw = NULL,
  confounders_tv_outcome_raw = NULL,
  confounders_tv_treatment_raw = NULL,
  ...
) {
  if (length(treatment) > 1L) {
    rlang::abort(
      c(
        "Multivariate longitudinal AIPW is not supported.",
        i = "Use `estimator = 'gcomp'` (ICE) for joint interventions."
      ),
      class = "causatr_longitudinal_multivariate_pending"
    )
  }

  check_reserved_cols(data)

  # -- Outcome-side metadata (same resolution as fit_ice) ---------
  time_points <- sort(unique(data[[time]]))
  n_times <- length(time_points)

  baseline_terms <- attr(
    stats::terms(confounders),
    "term.labels"
  )
  em_info <- parse_effect_mod(confounders, treatment)

  tv_vars <- if (!is.null(confounders_tv)) {
    all.vars(confounders_tv)
  } else {
    character(0)
  }
  tv_terms <- if (!is.null(confounders_tv)) {
    attr(stats::terms(confounders_tv), "term.labels")
  } else {
    character(0)
  }

  max_lag <- if (is.infinite(history)) {
    n_times - 1L
  } else {
    as.integer(history)
  }

  family_obj <- if (is.character(family)) {
    get(family, mode = "function")()
  } else if (is.function(family)) {
    family()
  } else {
    family
  }
  # AIPW pseudo-outcomes can fall outside the original family's support
  # (e.g. negative for Poisson, outside [0,1] for binomial) because
  # the augmented value phi = m(d,L) + W*(Y - m(A,L)) is not a raw
  # outcome but an EIF contribution. Fitting the backward pseudo-
  # regression with the original family (e.g. binomial) on such
  # fractional / out-of-range values is incoherent; gaussian identity
  # regression is consistent regardless of the outcome family. Same
  # approach as lmtp::lmtp_sdr.
  family_pseudo <- if (family_obj$family == "gaussian") {
    family_obj
  } else {
    stats::gaussian()
  }

  dots <- list(...)

  # -- Propensity-side models (same loop as fit_longitudinal_ipw) --
  # EM terms strip from propensity formulas; the confounder_terms
  # slot gives propensity-safe baseline terms (no A:modifier).
  # When per-component confounders are supplied, use treatment
  # confounders for propensity formulas.
  ps_confounders <- confounders_treatment %||% confounders
  em_info_ipw <- check_confounders_treatment(
    ps_confounders,
    treatment,
    estimator = "ipw"
  )
  ps_baseline_terms <- em_info_ipw$confounder_terms
  if (length(ps_baseline_terms) == 0L) {
    ps_baseline_terms <- character(0L)
  }

  # Time-varying confounders for propensity formulas. When
  # per-component TV confounders are supplied, the propensity side
  # uses `confounders_tv_treatment`; otherwise falls back to `confounders_tv`.
  ps_tv <- confounders_tv_treatment %||% confounders_tv
  ps_tv_vars <- if (!is.null(ps_tv)) {
    all.vars(ps_tv)
  } else {
    character(0)
  }

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

  treatment_models_by_time <- vector("list", n_times)
  fit_data_by_time <- vector("list", n_times)
  per_period_formula <- vector("list", n_times)

  for (k in seq_len(n_times)) {
    rows_k <- data[[time]] == time_points[k]
    data_k <- data[rows_k]
    available_lags <- min(k - 1L, max_lag)

    ps_formula <- build_longitudinal_ps_formula(
      treatment = treatment,
      baseline_terms = ps_baseline_terms,
      tv_vars = ps_tv_vars,
      available_lags = available_lags,
      data_at_time = data_k
    )

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

    treatment_models_by_time[[k]] <- tm_k
    fit_data_by_time[[k]] <- data_k
    per_period_formula[[k]] <- ps_formula
  }
  names(treatment_models_by_time) <- as.character(time_points)
  names(fit_data_by_time) <- as.character(time_points)
  names(per_period_formula) <- as.character(time_points)

  # Final-period rows (for bootstrap and downstream alignment)
  last_time <- time_points[n_times]
  rows_final <- data[[time]] == last_time
  fit_rows_final <- rows_final & !is.na(data[[outcome]])

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
    estimator = "aipw",
    type = "longitudinal",
    estimand = estimand,
    id = id,
    time = time,
    censoring = censoring,
    history = history,
    numerator = NULL,
    weights_obj = NULL,
    match_obj = NULL,
    call = call,
    details = list(
      # ICE-side metadata
      time_points = time_points,
      n_times = n_times,
      baseline_terms = baseline_terms,
      tv_vars = tv_vars,
      tv_terms = tv_terms,
      max_lag = max_lag,
      em_info = em_info,
      model_fn = model_fn,
      family_outcome = family_obj,
      family_pseudo = family_pseudo,
      weights = weights,
      dots = dots,
      # IPW-side metadata
      treatment_models_by_time = treatment_models_by_time,
      fit_data_by_time = fit_data_by_time,
      per_period_formula = per_period_formula,
      fit_rows_final = fit_rows_final,
      propensity_model_fn = prop_model_fn,
      propensity_family = propensity_family,
      n_total = nrow(data)
    )
  )
}


#' Augmented ICE backward iteration (ICE-AIPW, Bang & Robins 2005)
#'
#' @description
#' Mirrors `ice_iterate()` but augments the pseudo-outcome at each
#' backward step with the density-ratio-weighted residual:
#' \deqn{\tilde{Y}_k(i) = \hat{m}_k(d_k, \bar{L}_k) +
#'   w_k(i) \cdot (\tilde{Y}_{k+1}(i) - \hat{m}_k(A_k, \bar{L}_k))}
#' where \eqn{w_k(i)} is the **single-period** density-ratio weight
#' at time k and \eqn{\tilde{Y}_{K+1} = Y}.
#'
#' @param fit A `causatr_fit` with `estimator = "aipw"`,
#'   `type = "longitudinal"`.
#' @param intervention A `causatr_intervention` or `NULL`.
#'
#' @return A list with `pseudo_final`, `models`, `data_iv`, `fit_ids`,
#'   `intervention`, and `cumulative_weights`.
#'
#' @noRd
ice_aipw_iterate <- function(fit, intervention) {
  data <- fit$data
  details <- fit$details
  outcome <- fit$outcome
  treatment <- fit$treatment
  id_col <- fit$id
  time_col <- fit$time
  censoring <- fit$censoring

  time_points <- details$time_points
  n_times <- details$n_times
  baseline_terms <- details$baseline_terms
  tv_vars <- details$tv_vars
  tv_terms <- details$tv_terms
  max_lag <- details$max_lag
  model_fn <- details$model_fn
  family_outcome <- details$family_outcome
  family_pseudo <- details$family_pseudo
  external_weights <- details$weights
  model_fn_dots <- details$dots
  em_info <- details$em_info
  treatment_models_by_time <- details$treatment_models_by_time
  fit_data_by_time <- details$fit_data_by_time

  binary_outcome <- is_binary_family(family_outcome)

  # IPSI shifts the propensity score, not the treatment value — longitudinal
  # AIPW needs a counterfactual treatment stream to iterate the sequential
  # outcome regressions.  Point AIPW handles IPSI as a special case, but the
  # longitudinal path does not support it.
  is_ipsi <- inherits(intervention, "causatr_intervention") &&
    intervention$type == "ipsi"
  if (is_ipsi) {
    rlang::abort(
      paste0(
        "`ipsi()` interventions are not supported under longitudinal AIPW. ",
        "Use `estimator = 'ipw'` for longitudinal IPSI, or rewrite the ",
        "intervention as `shift()` / `scale_by()` for AIPW."
      ),
      class = "causatr_longitudinal_ipsi_pending"
    )
  }

  # Build intervention-modified person-period frame (same as ICE).
  data_iv <- ice_apply_intervention_long(
    data,
    treatment,
    intervention,
    id_col,
    time_col
  )

  all_ids <- unique(data[[id_col]])
  id_chr <- as.character(all_ids)
  n_id <- length(all_ids)
  pseudo <- stats::setNames(rep(NA_real_, n_id), id_chr)
  # For non-gaussian outcomes: separate clipped vector used only as
  # the response in the gaussian backward pseudo-regression.  The
  # unclipped `pseudo` is used for residuals, the final mean, and
  # the sandwich.  Clipping EIF pseudos destroys variance.
  pseudo_reg <- pseudo

  models <- vector("list", n_times)
  names(models) <- as.character(time_points)
  fit_ids <- vector("list", n_times)
  names(fit_ids) <- as.character(time_points)

  # -- Precompute per-period density-ratio weights (forward) ------
  # W_period[i, k] = g_k(d_k | H_k) / f_k(A_k | H_k), the single-
  # period density ratio at time k. The Bang & Robins (2005) sequential
  # DR recursion augments with single-period weights at each backward
  # step:
  #   pseudo_k(i) = m_k(d_k, H_k) + w_k(i) * (pseudo_{k+1}(i) - m_k(A_k, H_k))
  # Using the CUMULATIVE product prod_{j<=k} w_j inside the recursion
  # would over-weight; the product emerges only when the full recursion
  # is expanded algebraically (Robins et al. 2004).
  #
  # Cumulative weights W_cumul are stored for the sandwich variance
  # engine, which needs them for the propensity score Jacobian.
  id_to_idx <- stats::setNames(seq_len(n_id), id_chr)
  W_period <- matrix(1, nrow = n_id, ncol = n_times)
  W_cumul <- matrix(1, nrow = n_id, ncol = n_times)

  for (k in seq_len(n_times)) {
    tm_k <- treatment_models_by_time[[k]]
    data_k <- fit_data_by_time[[k]]
    ids_k <- as.character(data_k[[id_col]])

    w_k <- compute_density_ratio_weights(
      tm_k,
      data_k,
      intervention,
      estimand = "ATE"
    )

    idx_k <- id_to_idx[ids_k]
    W_period[idx_k, k] <- w_k

    if (k == 1L) {
      W_cumul[idx_k, k] <- w_k
    } else {
      W_cumul[idx_k, k] <- W_cumul[idx_k, k - 1L] * w_k
    }
    if (k > 1L) {
      missing_k <- is.na(match(id_chr, ids_k))
      W_cumul[missing_k, k] <- W_cumul[missing_k, k - 1L]
    }
  }

  uncens <- is_uncensored(data, censoring)

  final_time <- time_points[n_times]
  final_idx <- n_times - 1L

  mask_final <- data[[time_col]] == final_time
  fit_mask <- mask_final & uncens & !is.na(data[[outcome]])

  fit_data_k <- data[fit_mask]
  fit_ids[[n_times]] <- as.character(fit_data_k[[id_col]])

  formula_k <- ice_build_formula(
    outcome,
    treatment,
    baseline_terms,
    tv_vars,
    tv_terms,
    final_idx,
    max_lag,
    fit_data_k,
    em_info
  )

  model_args <- list(formula = formula_k, data = fit_data_k)
  if (fn_accepts_family(model_fn)) {
    model_args$family <- family_outcome
  }
  if (!is.null(external_weights)) {
    model_args$weights <- external_weights[fit_mask]
  }
  models[[n_times]] <- replay_fit(
    model_fn,
    model_args,
    model_fn_dots
  )

  # Predict for all uncensored at final time
  pred_mask <- mask_final & uncens
  pred_ids <- as.character(data[pred_mask][[id_col]])
  pred_idx <- id_to_idx[pred_ids]

  preds_iv <- stats::predict(
    models[[n_times]],
    newdata = data_iv[pred_mask],
    type = "response"
  )
  # Observed-treatment predictions (for the residual)
  preds_obs <- stats::predict(
    models[[n_times]],
    newdata = data[pred_mask],
    type = "response"
  )
  resid_k <- data[pred_mask][[outcome]] - preds_obs

  # Initialize pseudo at the final time: phi_K(i) = m_K(d_K, H_K) + w_K * (Y - m_K(A_K, H_K))
  # This is the first step of the ICE-AIPW backward recursion (k = K).
  pseudo[pred_ids] <- preds_iv + W_period[pred_idx, n_times] * resid_k
  if (binary_outcome) {
    # `pseudo_reg` is the clipped version fed as the response to the
    # next backward pseudo-regression. The unclipped `pseudo` is used
    # for residuals and the final mean -- clipping it would bias the EIF.
    pseudo_reg[pred_ids] <- pmax(
      pmin(pseudo[pred_ids], 1 - 1e-5),
      1e-5
    )
  } else {
    pseudo_reg[pred_ids] <- pseudo[pred_ids]
  }

  if (n_times > 1L && binary_outcome) {
    warn_ice_boundary_saturation(preds_iv, final_time)
  }

  # -- Steps 2+: backward iteration (time K-1 down to 0) ---------
  # At each step k, augment the pseudo from the previous (later) step:
  #   pseudo_k(i) = m_k(d_k, H_k) + w_k(i) * (pseudo_{k+1}(i) - m_k(A_k, H_k))
  # This is the sequential DR condition: consistent if EITHER the k-th
  # outcome regression m_k OR the k-th propensity f_k is correct.
  for (step_i in seq(n_times - 1L, 1L, by = -1L)) {
    current_time <- time_points[step_i]
    time_idx <- step_i - 1L

    mask_current <- data[[time_col]] == current_time
    mask_uncens <- mask_current & uncens

    current_ids <- as.character(data[mask_uncens][[id_col]])
    pseudo_y <- pseudo_reg[current_ids]

    has_pseudo <- !is.na(pseudo_y)
    if (sum(has_pseudo) == 0L) {
      rlang::abort(
        paste0(
          "No valid pseudo-outcomes at time ",
          current_time,
          " for ICE-AIPW backward iteration."
        ),
        .call = FALSE
      )
    }

    fit_data_k <- data.table::copy(data[mask_uncens][has_pseudo])
    fit_data_k[, .pseudo_y := pseudo_y[has_pseudo]]
    fit_ids[[step_i]] <- as.character(fit_data_k[[id_col]])

    formula_k <- ice_build_formula(
      ".pseudo_y",
      treatment,
      baseline_terms,
      tv_vars,
      tv_terms,
      time_idx,
      max_lag,
      fit_data_k,
      em_info
    )

    model_args_k <- list(formula = formula_k, data = fit_data_k)
    if (fn_accepts_family(model_fn)) {
      model_args_k$family <- family_pseudo
    }
    if (!is.null(external_weights)) {
      model_args_k$weights <- external_weights[mask_uncens][has_pseudo]
    }
    models[[step_i]] <- replay_fit(
      model_fn,
      model_args_k,
      model_fn_dots
    )

    # Predict for ALL individuals at the current time
    pred_ids_all <- as.character(data[mask_current][[id_col]])
    pred_idx_all <- id_to_idx[pred_ids_all]

    # Intervention predictions
    pred_all_iv <- data_iv[mask_current]
    preds_iv <- stats::predict(
      models[[step_i]],
      newdata = pred_all_iv,
      type = "response"
    )

    # Observed-treatment predictions (for the residual)
    preds_obs <- stats::predict(
      models[[step_i]],
      newdata = data[mask_current],
      type = "response"
    )

    # Residual = pseudo_{k+1}(i) - m_k(A_k, H_k): uses the UNCLIPPED
    # pseudo from the previous (later) step as the "outcome" and
    # subtracts the observed-treatment prediction. Using `pseudo_reg`
    # (clipped) here would bias the residual and break double robustness.
    resid_k <- pseudo[pred_ids_all] - preds_obs

    # AIPW augmentation at step k:
    #   pseudo_k(i) = m_k(d_k, H_k) + w_k(i) * resid_k(i)
    # Where pseudo_{k+1} is NA (individual was censored at a later time
    # and never received a forward prediction), fall back to the vanilla
    # ICE prediction m_k(d_k, H_k) -- resid would be NA and propagate.
    has_prev_pseudo <- !is.na(pseudo[pred_ids_all])
    aipw_pseudo <- preds_iv
    aipw_pseudo[has_prev_pseudo] <- preds_iv[has_prev_pseudo] +
      W_period[pred_idx_all[has_prev_pseudo], step_i] *
        resid_k[has_prev_pseudo]
    pseudo[pred_ids_all] <- aipw_pseudo
    if (binary_outcome) {
      # Clip only the version used as a regression response (pseudo_reg);
      # keep `pseudo` unclipped for use in the next step's residual.
      pseudo_reg[pred_ids_all] <- pmax(
        pmin(aipw_pseudo, 1 - 1e-5),
        1e-5
      )
    } else {
      pseudo_reg[pred_ids_all] <- aipw_pseudo
    }

    if (step_i > 1L && binary_outcome) {
      warn_ice_boundary_saturation(preds_iv, current_time)
    }
  }

  # After the full backward pass, pseudo[i] at the first time point is
  # the individual-level ICE-AIPW EIF value phi_1(H_{i,1}). Its mean
  # over the target population is the doubly-robust estimate mu_hat(d).
  first_time <- time_points[1]
  rows_first <- data[[time_col]] == first_time
  first_ids <- as.character(data[rows_first][[id_col]])
  pseudo_final <- unname(pseudo[first_ids])

  list(
    pseudo_final = pseudo_final,
    models = models,
    data_iv = data_iv,
    fit_ids = fit_ids,
    intervention = intervention,
    period_weights = W_period,
    cumulative_weights = W_cumul
  )
}


#' Compute longitudinal AIPW per-intervention contrasts
#'
#' @description
#' For each intervention, runs the ICE-AIPW backward iteration via
#' `ice_aipw_iterate()` and returns the marginal mean of the baseline
#' pseudo-outcomes over the target population.
#'
#' @param fit A `causatr_fit` with `estimator = "aipw"`,
#'   `type = "longitudinal"`.
#' @param interventions Named list of interventions.
#' @param target_within_first Logical vector over first-time-point
#'   individuals flagging the target population.
#'
#' @return A list with `ice_aipw_results` (named list of
#'   `ice_aipw_iterate()` outputs) and `mu_hat` (named numeric vector).
#'
#' @noRd
compute_aipw_contrast_longitudinal <- function(
  fit,
  interventions,
  target_within_first
) {
  int_names <- names(interventions)
  ext_w <- fit$details$weights
  time_col <- fit$time
  first_time <- fit$details$time_points[1]
  rows_first <- fit$data[[time_col]] == first_time

  ice_aipw_results <- stats::setNames(
    lapply(interventions, function(iv) ice_aipw_iterate(fit, iv)),
    int_names
  )

  if (!is.null(ext_w)) {
    w_first <- ext_w[rows_first]
    w_target <- w_first[target_within_first]
  } else {
    w_target <- NULL
  }
  mu_hat <- vapply(
    ice_aipw_results,
    function(res) {
      maybe_weighted_mean(
        res$pseudo_final[target_within_first],
        w_target
      )
    },
    numeric(1)
  )
  names(mu_hat) <- int_names

  list(
    ice_aipw_results = ice_aipw_results,
    mu_hat = mu_hat
  )
}
