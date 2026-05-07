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
  stabilize = "none",
  id = NULL,
  time = NULL,
  call,
  target = NULL,
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
    stabilize = stabilize,
    call = call,
    target = target,
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
  stabilize = "none",
  call,
  target = NULL,
  ...
) {
  # --- Effect modification parsing -----------------------------------------
  em_info <- check_confounders_treatment(
    confounders,
    treatment,
    estimator = "ipw"
  )

  is_multivariate <- length(treatment) > 1L

  if (!is_multivariate && stabilize != "none") {
    rlang::abort(
      c(
        paste0(
          "`stabilize = '",
          stabilize,
          "'` is only supported for multivariate AIPW."
        ),
        i = paste0(
          "Univariate AIPW already uses Hajek normalization; ",
          "a separate numerator model is not implemented ",
          "for single-treatment weights."
        )
      ),
      class = "causatr_stabilize_univariate"
    )
  }

  # --- Shared fit rows (outcome-clean + uncensored + S=1 if transport) -----
  fit_rows <- get_fit_rows(data, outcome, censoring, target = target)
  fit_data <- data[fit_rows]

  # --- Outcome model E[Y | A, L] ------------------------------------------
  model_formula <- build_outcome_formula(outcome, treatment, confounders)

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

  # --- Treatment density model(s) f(A | L) --------------------------------
  n_fit_outcome <- sum(fit_rows)

  # Resolve propensity fitter. In AIPW, model_fn is the outcome fitter;
  # do not silently reuse it for propensity. When the user omits
  # propensity_model_fn, default to stats::glm with a warning.
  if (is.null(propensity_model_fn)) {
    rlang::warn(
      c(
        paste0(
          "No `propensity_model_fn` specified; defaulting to ",
          "`stats::glm` for the treatment density model(s)."
        ),
        i = paste0(
          "Set `propensity_model_fn` explicitly if you need a ",
          "different fitter (e.g. `mgcv::gam`)."
        )
      ),
      class = "causatr_propensity_fn_default"
    )
    propensity_model_fn <- stats::glm
  }

  if (is_multivariate) {
    tm_args <- list(
      data = fit_data,
      treatment = treatment,
      confounders = confounders,
      model_fn = stats::glm,
      propensity_model_fn = propensity_model_fn,
      propensity_family = propensity_family,
      stabilize = stabilize
    )
    if (!is.null(model_weights)) {
      tm_args$weights <- model_weights
    }
    treatment_models <- do.call(
      fit_treatment_models,
      c(tm_args, dots)
    )

    for (k in seq_along(treatment_models)) {
      n_fit_k <- sum(treatment_models[[k]]$fit_rows)
      if (n_fit_k != n_fit_outcome) {
        rlang::abort(
          paste0(
            "Treatment density model for component '",
            treatment[k],
            "' kept ",
            n_fit_k,
            " rows but the outcome-non-missing mask has ",
            n_fit_outcome,
            " rows. Drop or impute the offending rows ",
            "before calling `causat()`."
          )
        )
      }
    }

    tm <- NULL
    propensity_model <- NULL
    prop_model_fn <- propensity_model_fn
  } else {
    trt_family <- detect_treatment_family(data[[treatment]])
    # Categorical and negbin need specialized fitters regardless
    # of what the user passed.
    prop_model_fn <- if (trt_family == "categorical") {
      nnet::multinom
    } else if (identical(propensity_family, "negbin")) {
      check_pkg("MASS")
      MASS::glm.nb
    } else {
      propensity_model_fn
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

    treatment_models <- structure(
      list(tm),
      class = c("causatr_treatment_models", "list"),
      names = treatment
    )
    propensity_model <- tm$model
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
      treatment_models = treatment_models,
      propensity_model = propensity_model,
      propensity_model_fn = prop_model_fn,
      propensity_family = propensity_family,
      em_info = em_info,
      is_multivariate = is_multivariate,
      stabilize = stabilize
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
  is_mv <- isTRUE(fit$details$is_multivariate)
  tms <- fit$details$treatment_models
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

  # Transport: compute sampling weights for the augmentation term.
  # Under transport the AIPW functional splits into two sums over
  # different populations (Dahabreh et al. 2020, Section 4.2):
  #   Term 1: (1/n_target) sum_{target} m(d(A,L), L)
  #   Term 2: (1/n_target) sum_{S=1} w_S * W_A * (Y - m(A,L))
  # The sampling weights w_S reweight study rows to the target.
  is_transport <- isTRUE(fit$details$transport)
  w_S_fit <- NULL
  if (is_transport) {
    w_S <- compute_sampling_weights(
      fit$details$sampling_model,
      data,
      fit$target,
      fit$target_subset
    )
    w_S_fit <- w_S[fit_rows]
  }

  # Target population within fit rows (for non-transport) or across
  # all rows (for transport Term 1, which sums over target rows that
  # may be outside fit_rows since target rows have no Y).
  target_fit <- target_idx[fit_rows]
  ext_w_fit <- if (!is.null(ext_w)) ext_w[fit_rows] else NULL

  # Transport Term 1: predictions on target rows (may include S=0
  # rows not in fit_data). Predict on all target rows using the
  # outcome model trained on S=1.
  if (is_transport) {
    target_data <- data[target_idx]
    n_target <- if (!is.null(ext_w)) {
      sum(ext_w[target_idx])
    } else {
      sum(target_idx)
    }
  }

  bundles <- lapply(int_names, function(nm) {
    iv <- interventions[[nm]]

    # Density-ratio weights W_i(g): univariate uses the single propensity
    # model; multivariate builds the joint weight as a product of K
    # per-component density ratios under the chain-rule factorisation.
    w_iv <- if (is_mv) {
      tms_local <- tms
      for (k in seq_along(tms_local)) {
        tms_local[[k]]$fit_rows <- rep(TRUE, nrow(fit_data))
      }
      class(tms_local) <- c("causatr_treatment_models", "list")
      compute_density_ratio_weights_mv(
        tms_local,
        fit_data,
        iv,
        estimand = estimand
      )
    } else {
      compute_density_ratio_weights(tm, fit_data, iv, estimand = estimand)
    }

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

    if (is_transport) {
      # AIPW transport (Dahabreh et al. 2020, Section 4.2):
      #   mu(d) = (1/n_T) sum_{target} m(d,L) + (1/n_T) sum_{S=1} w_S * W_A * resid
      # Term 1 sums over target rows; Term 2 sums over study (fit) rows.

      # Term 1: outcome-model predictions on target rows
      if (is_ipsi) {
        preds_target <- stats::predict(
          outcome_model,
          newdata = target_data,
          type = "response"
        )
      } else {
        target_data_a <- apply_intervention(target_data, treatment, iv)
        preds_target <- stats::predict(
          outcome_model,
          newdata = target_data_a,
          type = "response"
        )
      }
      if (!is.null(ext_w)) {
        term1 <- sum(ext_w[target_idx] * preds_target) / n_target
      } else {
        term1 <- mean(preds_target)
      }

      # Term 2: sampling-weighted augmentation on study rows
      aug_contrib <- w_S_fit * w_iv * resid_obs
      if (!is.null(ext_w_fit)) {
        term2 <- sum(ext_w_fit * aug_contrib) / n_target
      } else {
        term2 <- sum(aug_contrib) / n_target
      }

      mu_j <- term1 + term2

      # Per-row AIPW contributions on fit_rows (for variance engine).
      # The sandwich needs individual phi_i values on the study rows.
      aipw_contrib <- preds_g + w_S_fit * w_iv * resid_obs
    } else {
      # Standard AIPW (no transport):
      #   phi_i = m(d(A_i,L_i), L_i) + W_i(d) * (Y_i - m(A_i, L_i))
      aipw_contrib <- preds_g + w_iv * resid_obs

      if (!is.null(ext_w_fit)) {
        w_target <- ext_w_fit * target_fit
        mu_j <- sum(w_target * aipw_contrib) / sum(w_target)
      } else {
        mu_j <- mean(aipw_contrib[target_fit])
      }
    }

    list(
      intervention = iv,
      preds_g = preds_g,
      w_iv = w_iv,
      mu_hat = mu_j,
      data_a = data_a,
      w_S_fit = w_S_fit,
      preds_target = if (is_transport) preds_target else NULL
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
    resid_obs = resid_obs,
    is_transport = is_transport,
    target_idx = target_idx
  )
}
