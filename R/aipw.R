#' Fit an AIPW model for causal estimation (point treatment)
#'
#' @description
#' Implements the augmented inverse probability weighting estimator
#' (Robins, Rotnitzky & Zhao 1994). AIPW combines the outcome model
#' from g-computation with the density-ratio weights from IPW, giving
#' an estimator that is consistent if **either** nuisance model is
#' correctly specified ("doubly robust") and semiparametrically
#' efficient when both are correct.
#'
#' The AIPW functional for intervention \eqn{g} is:
#' \deqn{\hat\psi(g) = \frac{1}{n_t}\sum_{i \in \mathrm{target}}
#'   \bigl[\hat{m}(g(L_i, A_i), L_i) +
#'   W_i(g) \cdot (Y_i - \hat{m}(A_i, L_i))\bigr]}
#'
#' where \eqn{\hat{m}} is the outcome model and \eqn{W_i(g)} is the
#' density-ratio weight produced by the self-contained IPW engine.
#'
#' @param data data.table from `prepare_data()`.
#' @param outcome Character. Outcome column name.
#' @param treatment Character scalar or vector. Treatment column name(s).
#' @param confounders One-sided formula of baseline confounders.
#' @param confounders_tv One-sided formula of time-varying confounders or `NULL`.
#' @param family Character or family object for the outcome model.
#' @param estimand Character. `"ATE"`, `"ATT"`, or `"ATC"`.
#' @param type Character. `"point"` or `"longitudinal"`.
#' @param history Positive integer or `Inf`. Markov order for longitudinal.
#' @param censoring Character column name for censoring indicator, or `NULL`.
#' @param weights Numeric vector or `NULL`. External observation weights.
#' @param model_fn Function. Outcome model fitting function.
#' @param propensity_model_fn Function or `NULL`. Treatment density fitter.
#' @param propensity_family Character or `NULL`. Treatment family override.
#' @param id Character or `NULL`. ID column for longitudinal data.
#' @param time Character or `NULL`. Time column for longitudinal data.
#' @param call The original `causat()` call.
#' @param ... Passed to `propensity_model_fn` via `fit_treatment_model()`.
#'
#' @return A `causatr_fit` object with `estimator = "aipw"`, the outcome
#'   model in `$model`, and both treatment and outcome model details in
#'   `$details`.
#'
#' @noRd
fit_aipw <- function(
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
  propensity_model_fn,
  propensity_family,
  id = NULL,
  time = NULL,
  call,
  ...
) {
  if (type == "longitudinal") {
    return(fit_aipw_longitudinal(
      data = data,
      outcome = outcome,
      treatment = treatment,
      confounders = confounders,
      confounders_tv = confounders_tv,
      family = family,
      estimand = estimand,
      history = history,
      censoring = censoring,
      weights = weights,
      model_fn = model_fn,
      propensity_model_fn = propensity_model_fn,
      propensity_family = propensity_family,
      id = id,
      time = time,
      call = call,
      ...
    ))
  }

  if (length(treatment) > 1L) {
    rlang::abort(
      c(
        "Multivariate treatment is not yet supported under AIPW.",
        i = "Use `estimator = 'gcomp'` or `estimator = 'ipw'` for multivariate treatments."
      ),
      class = "causatr_aipw_multivariate_pending"
    )
  }

  fit_aipw_point(
    data = data,
    outcome = outcome,
    treatment = treatment,
    confounders = confounders,
    family = family,
    estimand = estimand,
    censoring = censoring,
    weights = weights,
    model_fn = model_fn,
    propensity_model_fn = propensity_model_fn,
    propensity_family = propensity_family,
    call = call,
    ...
  )
}


#' Fit both nuisance models for point-treatment AIPW
#'
#' @description
#' Fits the outcome model \eqn{E[Y \mid A, L]} (same as g-comp) and
#' the treatment density model \eqn{f(A \mid L)} (same as IPW) on the
#' same row set, returning a `causatr_fit` that carries both for
#' downstream use in `compute_aipw_contrast_point()`.
#'
#' @noRd
fit_aipw_point <- function(
  data,
  outcome,
  treatment,
  confounders,
  family,
  estimand,
  censoring,
  weights,
  model_fn,
  propensity_model_fn,
  propensity_family,
  call,
  ...
) {
  # --- Effect modification parsing -----------------------------------------
  em_info <- check_confounders_treatment(
    confounders,
    treatment,
    estimator = "ipw"
  )

  # --- Shared fit rows (outcome-clean + uncensored) ------------------------
  fit_rows <- get_fit_rows(data, outcome, censoring)
  fit_data <- data[fit_rows]

  # --- Outcome model E[Y | A, L] ------------------------------------------
  confounder_terms <- attr(stats::terms(confounders), "term.labels")
  rhs <- c(treatment, confounder_terms)
  model_formula <- stats::reformulate(rhs, response = outcome)

  model_weights <- if (!is.null(weights)) weights[fit_rows] else NULL
  resolved_family <- resolve_family(family)

  dots <- list(...)

  model_args <- list(model_formula, data = fit_data)
  if (fn_accepts_family(model_fn)) {
    model_args$family <- resolved_family
  }
  if (!is.null(model_weights)) {
    model_args$weights <- model_weights
  }
  outcome_model <- replay_fit(model_fn, model_args, dots)

  # --- Treatment density model f(A | L) -----------------------------------
  trt_family <- detect_treatment_family(data[[treatment]])
  prop_model_fn <- if (!is.null(propensity_model_fn)) {
    propensity_model_fn
  } else if (trt_family == "categorical") {
    nnet::multinom
  } else if (identical(propensity_family, "negbin")) {
    check_pkg("MASS")
    MASS::glm.nb
  } else {
    model_fn
  }

  tm_args <- list(
    data = fit_data,
    treatment = treatment,
    confounders = confounders,
    model_fn = prop_model_fn,
    propensity_family = propensity_family
  )
  if (!is.null(model_weights)) {
    tm_args$weights <- model_weights
  }
  tm <- do.call(fit_treatment_model, c(tm_args, dots))

  # Row-alignment invariant: both models must operate on the same rows.
  n_fit_outcome <- sum(fit_rows)
  n_fit_ps <- sum(tm$fit_rows)
  if (n_fit_ps != n_fit_outcome) {
    rlang::abort(
      paste0(
        "Treatment density model kept ",
        n_fit_ps,
        " rows but the outcome-non-missing mask has ",
        n_fit_outcome,
        " rows. A confounder column has missing values the outcome ",
        "mask does not exclude. Drop or impute those rows before ",
        "calling `causat()` so the outcome and propensity fits agree."
      )
    )
  }

  new_causatr_fit(
    model = outcome_model,
    data = data,
    treatment = treatment,
    outcome = outcome,
    confounders = confounders,
    confounders_tv = NULL,
    family = family,
    estimator = "aipw",
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
      n_fit = n_fit_outcome,
      n_total = nrow(data),
      model_fn = model_fn,
      weights = weights,
      dots = dots,
      treatment_model = tm,
      propensity_model = tm$model,
      propensity_model_fn = prop_model_fn,
      propensity_family = propensity_family,
      em_info = em_info
    )
  )
}


#' Compute AIPW per-intervention contrasts for point treatment
#'
#' @description
#' For each intervention, computes the AIPW functional:
#' \deqn{\hat\mu(g) = \frac{1}{n_t}\sum_{i \in \mathrm{target}}
#'   \bigl[\hat{m}(g(L_i, A_i), L_i) +
#'   W_i(g) \cdot (Y_i - \hat{m}(A_i, L_i))\bigr]}
#'
#' @param fit A `causatr_fit` with `estimator = "aipw"`.
#' @param interventions Named list of interventions.
#' @param target_idx Logical vector (length n) flagging the target population.
#'
#' @return A list with components:
#'   - `mu_hat`: named numeric vector of marginal means.
#'   - `bundles`: per-intervention list carrying preds_g, preds_obs, resid, w_iv.
#'   - `fit_idx`: integer vector of fit rows.
#'   - `n_total`: total row count.
#'
#' @noRd
compute_aipw_contrast_point <- function(fit, interventions, target_idx) {
  data <- fit$data
  treatment <- fit$treatment
  outcome <- fit$outcome
  outcome_model <- fit$model
  tm <- fit$details$treatment_model
  fit_rows <- fit$details$fit_rows
  fit_data <- data[fit_rows]
  ext_w <- fit$details$weights
  estimand <- fit$estimand
  n_total <- nrow(data)
  fit_idx <- which(fit_rows)
  int_names <- names(interventions)

  # Observed-treatment predictions and residuals (shared across interventions)
  preds_obs <- stats::predict(
    outcome_model,
    newdata = fit_data,
    type = "response"
  )
  y_obs <- fit_data[[outcome]]
  resid_obs <- y_obs - preds_obs

  # Target population within fit rows
  target_fit <- target_idx[fit_rows]
  ext_w_fit <- if (!is.null(ext_w)) ext_w[fit_rows] else NULL

  bundles <- lapply(int_names, function(nm) {
    iv <- interventions[[nm]]

    # Density-ratio weights W_i(g) — computed first so that
    # check_intervention_family_compat() fires rejection errors
    # (stochastic, threshold on continuous, etc.) before we attempt
    # apply_intervention() which would fail differently.
    w_iv <- compute_density_ratio_weights(tm, fit_data, iv, estimand = estimand)

    # IPSI shifts the propensity, not the treatment value — there is no
    # counterfactual A to predict at. The outcome-model augmentation
    # evaluates at observed A, so preds_g = preds_obs and the AIPW
    # functional collapses to the IPW Hajek mean with an outcome-model
    # augmentation that adds (W-1)*resid (since preds_g + W*resid - preds_g = W*resid
    # when preds_g = preds_obs; but we keep the general formula).
    is_ipsi <- inherits(iv, "causatr_intervention") && iv$type == "ipsi"
    if (is_ipsi) {
      data_a <- fit_data
      preds_g <- preds_obs
    } else {
      data_a <- apply_intervention(fit_data, treatment, iv)
      preds_g <- stats::predict(
        outcome_model,
        newdata = data_a,
        type = "response"
      )
    }

    # AIPW individual-level contributions on fit rows
    aipw_contrib <- preds_g + w_iv * resid_obs

    # Marginal mean over target population (with optional external weights)
    if (!is.null(ext_w_fit)) {
      w_target <- ext_w_fit * target_fit
      mu_j <- sum(w_target * aipw_contrib) / sum(w_target)
    } else {
      mu_j <- mean(aipw_contrib[target_fit])
    }

    list(
      intervention = iv,
      preds_g = preds_g,
      w_iv = w_iv,
      mu_hat = mu_j,
      data_a = data_a
    )
  })
  names(bundles) <- int_names

  mu_hat <- vapply(bundles, `[[`, numeric(1), "mu_hat")
  names(mu_hat) <- int_names

  list(
    mu_hat = mu_hat,
    bundles = bundles,
    fit_idx = fit_idx,
    n_total = n_total,
    preds_obs = preds_obs,
    resid_obs = resid_obs
  )
}


#' Fit both nuisance models for longitudinal AIPW (ICE-AIPW)
#'
#' @description
#' Composes the ICE outcome-model metadata (Phase 5) with the
#' longitudinal IPW per-period propensity models (Phase 10) into a
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
  ...
) {
  if (length(treatment) > 1L) {
    rlang::abort(
      c(
        "Multivariate longitudinal AIPW is not yet supported.",
        i = "Use `estimator = 'gcomp'` (ICE) for joint interventions.",
        i = "Multivariate longitudinal support ships in Phase 19."
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
  # For AIPW, pseudo-outcomes at non-terminal steps can exceed [0,1]
  # due to augmentation. Use gaussian for the backward pseudo-regression
  # (same approach as lmtp::lmtp_sdr). ICE uses quasibinomial because
  # its pseudo-outcomes are bounded predictions; AIPW's are not.
  family_pseudo <- if (family_obj$family == "binomial") {
    stats::gaussian()
  } else {
    family_obj
  }

  dots <- list(...)

  # -- Propensity-side models (same loop as fit_longitudinal_ipw) --
  # EM terms strip from propensity formulas; the confounder_terms
  # slot gives propensity-safe baseline terms (no A:modifier).
  em_info_ipw <- check_confounders_treatment(
    confounders,
    treatment,
    estimator = "ipw"
  )
  ps_baseline_terms <- em_info_ipw$confounder_terms
  if (length(ps_baseline_terms) == 0L) {
    ps_baseline_terms <- character(0L)
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
      tv_vars = tv_vars,
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

  models <- vector("list", n_times)
  names(models) <- as.character(time_points)
  fit_ids <- vector("list", n_times)
  names(fit_ids) <- as.character(time_points)

  # -- Precompute per-period density-ratio weights (forward) ------
  # W_period[i, k] = w_k(i), the single-period density ratio at time
  # k. The Bang & Robins (2005) sequential DR recursion uses
  # single-period weights at each backward step:
  #   pseudo_k = m_k(d_k, L) + w_k * (pseudo_{k+1} - m_k(A_obs, L))
  # NOT cumulative products. The cumulative product arises from
  # expanding the recursion but must not be used inside it.
  #
  # We also store cumulative weights (for the sandwich variance
  # engine, which needs them for the propensity Jacobian).
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

  # -- Step 1: final time point (real outcome Y) ------------------
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

  # AIPW pseudo = intervention prediction + per-period weight * residual
  pseudo[pred_ids] <- preds_iv + W_period[pred_idx, n_times] * resid_k
  if (binary_outcome) {
    pseudo[pred_ids] <- pmax(pmin(pseudo[pred_ids], 1 - 1e-5), 1e-5)
  }

  if (n_times > 1L && binary_outcome) {
    warn_ice_boundary_saturation(preds_iv, final_time)
  }

  # -- Steps 2+: backward iteration (time K-1 down to 0) ---------
  for (step_i in seq(n_times - 1L, 1L, by = -1L)) {
    current_time <- time_points[step_i]
    time_idx <- step_i - 1L

    mask_current <- data[[time_col]] == current_time
    mask_uncens <- mask_current & uncens

    current_ids <- as.character(data[mask_uncens][[id_col]])
    pseudo_y <- pseudo[current_ids]

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

    # Residual: pseudo from previous backward step minus observed pred
    resid_k <- pseudo[pred_ids_all] - preds_obs

    # AIPW pseudo: augmented with cumulative weight * residual.
    # Where pseudo from the previous step is NA (censored at a later
    # time), fall back to the vanilla ICE prediction (resid = NA
    # would propagate). This mirrors ICE's has_pseudo filtering.
    has_prev_pseudo <- !is.na(pseudo[pred_ids_all])
    aipw_pseudo <- preds_iv
    aipw_pseudo[has_prev_pseudo] <- preds_iv[has_prev_pseudo] +
      W_period[pred_idx_all[has_prev_pseudo], step_i] *
        resid_k[has_prev_pseudo]
    if (binary_outcome) {
      aipw_pseudo <- pmax(pmin(aipw_pseudo, 1 - 1e-5), 1e-5)
    }
    pseudo[pred_ids_all] <- aipw_pseudo

    if (step_i > 1L && binary_outcome) {
      warn_ice_boundary_saturation(preds_iv, current_time)
    }
  }

  # Return baseline pseudo-outcomes
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
