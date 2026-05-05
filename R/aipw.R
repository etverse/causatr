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
    rlang::abort(
      c(
        "Longitudinal AIPW (ICE-AIPW) is not yet implemented.",
        i = "Use `estimator = 'gcomp'` (ICE) or `estimator = 'ipw'` (longitudinal IPW) for now."
      ),
      class = "causatr_aipw_longitudinal_pending"
    )
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
