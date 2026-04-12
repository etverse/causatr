#' Sandwich (robust) variance–covariance matrix for marginal means
#'
#' @description
#' Computes Var(μ̂) by propagating the robust parameter vcov to the marginal
#' means via the delta method (J V_β Jᵀ).  The parameter vcov V_β depends
#' on the estimation method:
#'
#' - **g-computation**: `sandwich::sandwich(model)` — standard Huber–White.
#' - **IPW**: `vcov(glm_weightit_model)` — M-estimation sandwich from
#'   `WeightIt` that accounts for weight-estimation uncertainty.
#' - **Matching**: `sandwich::vcovCL(model, cluster = ~subclass)` —
#'   cluster-robust SE on the matched-pair subclass.
#'
#' @param fit A `causatr_fit` object.
#' @param data_a_list Named list of counterfactual data.tables (one per
#'   intervention).
#' @param preds_list Named list of length-n prediction vectors (one per
#'   intervention).  Kept for interface consistency.
#' @param target_idx Logical vector (length n) flagging target-population rows.
#'
#' @return A named k × k variance–covariance matrix.
#'
#' @references
#' Arel-Bundock V (2024). marginaleffects: Predictions, Comparisons,
#' Slopes, Marginal Means, and Hypothesis Tests.
#' \url{https://marginaleffects.com/}
#'
#' @noRd
variance_sandwich <- function(fit, data_a_list, preds_list, target_idx) {
  model <- fit$model

  has_family <- !is.null(model$family) &&
    !is.null(model$family$mu.eta) &&
    !is.null(model$family$variance)

  if (fit$method == "gcomp" && has_family) {
    return(gcomp_variance_sandwich(fit, data_a_list, preds_list, target_idx))
  }

  # Choose the appropriate robust parameter vcov based on the estimation method.
  if (fit$method == "ipw") {
    # glm_weightit objects have a vcov() method that provides M-estimation
    # sandwich SEs accounting for weight-estimation uncertainty.
    V_beta <- stats::vcov(model)
  } else if (fit$method == "matching") {
    # Cluster-robust SE on the matched-pair subclass.
    matched_data <- fit$details$matched_data
    V_beta <- sandwich::vcovCL(
      model,
      cluster = matched_data[["subclass"]]
    )
  } else {
    # g-computation with model lacking family (rare): Huber–White + J V_β Jᵀ.
    V_beta <- sandwich::sandwich(model)
  }

  # Propagate parameter uncertainty: J V_β Jᵀ.
  ext_w <- fit$details$weights
  w_target <- if (!is.null(ext_w)) ext_w[target_idx] else NULL
  compute_vcov_marginal(
    model,
    data_a_list,
    target_idx,
    V_beta,
    weights = w_target
  )
}


#' Full IF sandwich for point-treatment g-computation
#'
#' @description
#' Uses the same per-individual influence function approach as
#' `ice_compute_influence()` but for a single outcome model (K = 1).
#' Accounts for both parameter uncertainty and covariate sampling
#' variability. Works with GLMs and GAMs.
#'
#' @noRd
gcomp_variance_sandwich <- function(fit, data_a_list, preds_list, target_idx) {
  model <- fit$model
  int_names <- names(data_a_list)
  k <- length(int_names)
  n <- nrow(fit$data)

  ext_w <- fit$details$weights
  has_weights <- !is.null(ext_w)

  if (has_weights) {
    w_target <- ifelse(target_idx, ext_w, 0)
    sum_w_target <- sum(w_target[target_idx])
  } else {
    n_target <- sum(target_idx)
  }

  fit_idx <- which(fit$details$fit_rows)
  na_action <- model$na.action
  if (!is.null(na_action)) {
    fit_idx <- fit_idx[-na_action]
  }
  beta_hat <- stats::coef(model)

  data_a_frames <- lapply(data_a_list, function(da) {
    as.data.frame(da)[target_idx, , drop = FALSE]
  })

  IF_list <- lapply(seq_len(k), function(j) {
    p <- preds_list[[j]]
    if (has_weights) {
      mu_hat <- sum(w_target * p) / sum_w_target
      IF_vec <- n * (w_target / sum_w_target) * (p - mu_hat)
    } else {
      mu_hat <- mean(p[target_idx])
      IF_vec <- (n / n_target) * ifelse(target_idx, p - mu_hat, 0)
    }

    X_star <- intervention_design_matrix(model, data_a_frames[[j]])
    eta_star <- as.numeric(X_star %*% beta_hat)
    mu_eta_star <- model$family$mu.eta(eta_star)

    if (has_weights) {
      w_t <- ext_w[target_idx]
      g <- as.numeric(crossprod(X_star, w_t * mu_eta_star)) / sum_w_target
    } else {
      g <- as.numeric(crossprod(X_star, mu_eta_star)) / n_target
    }

    IF_vec <- IF_vec + model_if_correction(model, g, fit_idx, n)

    IF_vec
  })

  vcov_from_if(IF_list, n, int_names)
}


#' Compute per-individual model correction for the influence function
#'
#' @description
#' Shared workhorse used by both `gcomp_variance_sandwich()` (point treatment,
#' K = 1) and `ice_compute_influence()` (longitudinal, K models in backward
#' chain). Given a sensitivity gradient `g` (how the marginal mean depends on
#' this model's parameters), computes the n × d × r correction term.
#'
#' @param model A fitted model with a `family` object (GLM or GAM).
#' @param g Sensitivity gradient (p-vector): ∂μ̂/∂β for this model.
#' @param fit_idx Integer vector of row indices (in 1..n) used for fitting.
#' @param n Total number of individuals.
#'
#' @return A length-n numeric vector of IF corrections (zero for non-fit rows).
#'
#' @noRd
model_if_correction <- function(model, g, fit_idx, n) {
  X_fit <- stats::model.matrix(model)
  r_fit <- stats::residuals(model, type = "response")
  bread_inv <- model_bread_inv(model, X_fit)

  h <- as.numeric(bread_inv %*% g)
  d_fit <- as.numeric(X_fit %*% h)

  correction <- rep(0, n)
  correction[fit_idx] <- n * d_fit * r_fit
  correction
}


#' Inverse of the model bread matrix
#'
#' @description
#' For standard GLMs, computes (X'WX)⁻¹ from the design matrix and IWLS
#' working weights. For GAMs, returns the penalized Bayesian covariance
#' `model$Vp` = (X'WX + λS)⁻¹.
#'
#' @param model A fitted model with a `family` object.
#' @param X_fit Design matrix from `model.matrix(model)`.
#'
#' @return A p × p matrix.
#'
#' @noRd
model_bread_inv <- function(model, X_fit) {
  if (inherits(model, "gam") && !is.null(model$Vp)) {
    return(model$Vp)
  }

  eta <- model$linear.predictors
  mu_eta <- model$family$mu.eta(eta)
  var_mu <- model$family$variance(stats::fitted(model))
  w <- mu_eta^2 / var_mu

  XtWX <- crossprod(X_fit, X_fit * w)
  tryCatch(
    solve(XtWX),
    error = function(e) MASS::ginv(XtWX)
  )
}


#' Design matrix for counterfactual (intervention) data
#'
#' @description
#' Constructs the model matrix for intervention data, handling differences
#' between model classes:
#' - **GAMs**: uses `predict(model, type = "lpmatrix")` for the smooth
#'   basis matrix.
#' - **GLMs** and other models: uses `model.matrix(terms, data, xlev)`.
#'
#' @param model A fitted model object.
#' @param newdata Data frame of counterfactual observations.
#'
#' @return A design matrix (n_new × p).
#'
#' @noRd
intervention_design_matrix <- function(model, newdata) {
  if (inherits(model, "gam")) {
    return(stats::predict(model, newdata = newdata, type = "lpmatrix"))
  }
  pred_terms <- stats::delete.response(stats::terms(model))
  xlev <- model$xlevels
  stats::model.matrix(pred_terms, data = newdata, xlev = xlev)
}


#' Variance–covariance matrix from per-individual influence functions
#'
#' @description
#' Computes Cov(μ̂_a, μ̂_b) = (1/n²) Σ IF_{a,i} × IF_{b,i} for all pairs
#' of interventions. Used by both `gcomp_variance_sandwich()` and
#' `ice_variance_sandwich()`.
#'
#' @param IF_list Named list of length-n IF vectors (one per intervention).
#' @param n Number of individuals.
#' @param int_names Character vector of intervention names.
#'
#' @return A named k × k variance–covariance matrix.
#'
#' @noRd
vcov_from_if <- function(IF_list, n, int_names) {
  k <- length(IF_list)
  vcov_mat <- matrix(0, k, k)
  for (j in seq_len(k)) {
    for (l in j:k) {
      vcov_mat[j, l] <- sum(IF_list[[j]] * IF_list[[l]]) / n^2
      if (j != l) vcov_mat[l, j] <- vcov_mat[j, l]
    }
  }
  rownames(vcov_mat) <- int_names
  colnames(vcov_mat) <- int_names
  vcov_mat
}


#' Sandwich variance for longitudinal ICE g-computation
#'
#' @description
#' Computes Var(μ̂) via the empirical sandwich estimator for stacked
#' estimating equations (Zivich et al. 2024). The ICE procedure — K
#' sequential outcome models plus the final mean — is expressed as a
#' single M-estimation problem. The per-individual efficient influence
#' function (EIF) propagates uncertainty from all K models to the marginal
#' mean via backward recursion through the model chain.
#'
#' For multiple interventions, the covariance between μ̂_a and μ̂_b is
#' computed from the cross-products of their influence functions (the same
#' individuals contribute to both ICE procedures).
#'
#' @param fit A `causatr_fit` of type `"longitudinal"`.
#' @param ice_results Named list of ICE iteration results (one per
#'   intervention, from `ice_iterate()`).
#' @param target_within_first Logical vector (length = number of
#'   individuals at first time) flagging the target population.
#'
#' @return A k × k variance–covariance matrix (k = number of
#'   interventions).
#'
#' @references
#' Zivich PN, Ross RK, Shook-Sa BE, Cole SR, Edwards JK (2024). Empirical
#' sandwich variance estimator for iterated conditional expectation
#' g-computation. *Statistics in Medicine* 43:5562–5572.
#'
#' @noRd
ice_variance_sandwich <- function(fit, ice_results, target_within_first) {
  int_names <- names(ice_results)

  data <- fit$data
  first_time <- fit$details$time_points[1]
  rows_first <- data[[fit$time]] == first_time
  n <- sum(rows_first)

  IF_list <- lapply(ice_results, function(res) {
    ice_compute_influence(fit, res, target_within_first)
  })

  vcov_from_if(IF_list, n, int_names)
}


#' Compute per-individual influence functions for ICE
#'
#' @description
#' The efficient influence function for the ICE g-comp estimator
#' decomposes into a **direct** term and **model correction** terms:
#'
#' \deqn{IF_i = \frac{target_i}{n_t}(\hat{Y}^*_{0,i} - \hat\mu)
#'   + \sum_k d_{k,i} \cdot r_{k,i}}
#'
#' where:
#' - \eqn{d_{k,i}} is the **sensitivity** of μ̂ to individual i's
#'   pseudo-outcome at model k (computed via backward recursion through
#'   the model chain).
#' - \eqn{r_{k,i}} is i's residual at model k (observed or pseudo
#'   response minus fitted value).
#'
#' The backward recursion for sensitivities:
#' 1. At model 0 (first time): g_0 = (1/n_t) Σ_{target} dμ/dη × X^*_{0,j}
#' 2. h_0 = (X_0^T W_0 X_0)^{-1} g_0
#' 3. d_{0,i} = h_0^T X_{0,i} (for individuals in model 0's fitting set)
#' 4. For k = 1, ..., K: g_k = Σ_{j in model_{k-1}} d_{k-1,j} × dμ/dη
#'    × X^*_{k,j}, then h_k and d_{k,i} as above.
#'
#' This recursion is the closed-form solution of the block-lower-triangular
#' bread matrix inversion in the stacked estimating equations.
#'
#' @param fit A `causatr_fit` of type `"longitudinal"`.
#' @param ice_result A single ICE iteration result (from `ice_iterate()`).
#' @param target Logical vector (length = individuals at first time)
#'   flagging the target population.
#'
#' @return A numeric vector of length n (one IF per individual).
#'
#' @noRd
ice_compute_influence <- function(fit, ice_result, target) {
  data <- fit$data
  details <- fit$details
  models <- ice_result$models
  data_iv <- ice_result$data_iv
  fit_ids_list <- ice_result$fit_ids
  pseudo_final <- ice_result$pseudo_final

  time_points <- details$time_points
  n_times <- details$n_times
  id_col <- fit$id
  time_col <- fit$time

  # Map individual IDs to integer indices (1..n).
  first_time <- time_points[1]
  rows_first <- data[[time_col]] == first_time
  all_ids <- as.character(data[rows_first][[id_col]])
  n <- length(all_ids)
  id_to_idx <- stats::setNames(seq_len(n), all_ids)

  # Direct IF term: IF_i = -A_mu^{-1} psi_mu,i.
  # The n/n_target (or n/sum_w) factor arises from A_mu^{-1} = -n/n_t.
  # This ensures the direct term has the same O(1) scaling as the model
  # correction terms (which carry the factor via g_k and the n * d * r
  # multiplication below). Without it, variance is underestimated by
  # (n_target/n)^2 when target is a proper subset.
  ext_w <- fit$details$weights
  has_weights <- !is.null(ext_w)
  if (has_weights) {
    w_first <- ext_w[rows_first]
    w_target <- ifelse(target, w_first, 0)
    sum_w_target <- sum(w_target)
    mu_hat <- sum(w_target * pseudo_final) / sum_w_target
    IF_vec <- n * (w_target / sum_w_target) * (pseudo_final - mu_hat)
  } else {
    n_target <- sum(target)
    mu_hat <- mean(pseudo_final[target], na.rm = TRUE)
    IF_vec <- (n / n_target) * ifelse(target, pseudo_final - mu_hat, 0)
  }

  # Sensitivity vector d (length n), reused across steps.
  d_vec <- rep(0, n)

  # Backward recursion for model correction terms.
  # Process models from step 1 (first time) to step n_times (final time).
  # At each step, compute the sensitivity d_{k,i} of μ̂ to individual i's

  # pseudo-outcome at model k, then add d_{k,i} × r_{k,i} to the IF.

  for (step_i in seq_len(n_times)) {
    model_k <- models[[step_i]]
    if (is.null(model_k)) {
      next
    }

    current_time <- time_points[step_i]
    fit_ids_k <- fit_ids_list[[step_i]]
    n_fit <- length(fit_ids_k)
    if (n_fit == 0L) {
      next
    }

    X_k <- stats::model.matrix(model_k)
    p_k <- ncol(X_k)
    residuals_k <- stats::residuals(model_k, type = "response")
    XtWX_inv <- model_bread_inv(model_k, X_k)

    # Compute g_k: the sensitivity gradient.
    #
    # g_k summarises how μ̂ depends on model k's parameters β_k.  The
    # derivation differs for model 0 (where μ̂ depends on β_0 through
    # the final predictions) vs. later models (where β_k affects
    # pseudo-outcomes for the model at k-1).

    # Get the intervention model matrix X^*_k (design matrix rows from
    # intervention data at the current time step).
    rows_iv_current <- data_iv[[time_col]] == current_time
    iv_data_current <- data_iv[rows_iv_current]
    iv_ids_current <- as.character(iv_data_current[[id_col]])

    # Build the intervention design matrix using model k's terms.
    pred_terms <- stats::delete.response(stats::terms(model_k))
    X_star_k <- tryCatch(
      stats::model.matrix(pred_terms, data = iv_data_current),
      error = function(e) {
        stats::model.matrix(pred_terms, data = as.data.frame(iv_data_current))
      }
    )

    # dμ/dη evaluated at the intervention predictions.
    eta_star <- as.numeric(X_star_k %*% stats::coef(model_k))
    mu_eta_star <- model_k$family$mu.eta(eta_star)

    if (step_i == 1L) {
      # Model 0 (first time): g_0 = (1/Σw_t) Σ_{target} w_j dμ/dη × X^*_{0,j}.
      target_in_iv <- match(all_ids[target], iv_ids_current)
      valid_target <- !is.na(target_in_iv)
      target_in_iv <- target_in_iv[valid_target]
      if (has_weights) {
        target_w <- w_target[target][valid_target]
        g_k <- as.numeric(
          crossprod(
            X_star_k[target_in_iv, , drop = FALSE],
            target_w * mu_eta_star[target_in_iv]
          )
        ) /
          sum_w_target
      } else {
        g_k <- as.numeric(
          crossprod(
            X_star_k[target_in_iv, , drop = FALSE],
            mu_eta_star[target_in_iv]
          )
        ) /
          n_target
      }
    } else {
      # Model k > 0: β_k affects μ̂ through pseudo-outcomes of model k-1.
      # g_k = Σ_{j ∈ model_{k-1}} d_{k-1,j} × dμ/dη × X^*_{k,j}.
      prev_fit_ids <- fit_ids_list[[step_i - 1L]]

      idx_in_all <- id_to_idx[prev_fit_ids]
      rows_in_iv <- match(prev_fit_ids, iv_ids_current)
      d_prev <- d_vec[idx_in_all]

      keep <- !is.na(idx_in_all) & !is.na(rows_in_iv) & d_prev != 0
      if (any(keep)) {
        weights_g <- d_prev[keep] * mu_eta_star[rows_in_iv[keep]]
        g_k <- as.numeric(
          crossprod(
            X_star_k[rows_in_iv[keep], , drop = FALSE],
            weights_g
          )
        )
      } else {
        g_k <- rep(0, p_k)
      }
    }

    # Solve for h_k and compute d_{k,i}.
    h_k <- as.numeric(XtWX_inv %*% g_k)

    # d_{k,i} = h_k^T X_{k,i} for each individual in model k's fitting set.
    d_vec <- rep(0, n)
    fit_id_idx <- id_to_idx[fit_ids_k]
    d_vec[fit_id_idx] <- as.numeric(X_k %*% h_k)

    # Add model correction: n × d_{k,i} × r_{k,i}.
    #
    # The factor of n arises because d_{k,i} = h_k^T X_{k,i} where
    # h_k = (X^T W X)^{-1} g_k uses the raw Hessian X^T W X = O(n),
    # not the averaged A_{kk} = (1/n) X^T W X. In the M-estimation
    # bread matrix, A_{kk}^{-1} = -n (X^T W X)^{-1}, so the correct
    # sensitivity is n × d_{k,i}. Without this factor, model correction
    # terms are O(1/n) per individual instead of O(1), which makes the
    # contrast variance O(1/n^3) instead of the correct O(1/n).
    r_vec <- rep(0, n)
    r_vec[fit_id_idx] <- residuals_k

    IF_vec <- IF_vec + n * d_vec * r_vec
  }

  IF_vec
}

#' Jacobian-based propagation of parameter uncertainty to marginal means
#'
#' @description
#' Shared workhorse for both `variance_sandwich()` and `variance_delta()`.
#' Given a parameter variance–covariance matrix `V_beta`, computes
#' `J V_beta Jᵀ` where J is the Jacobian of the marginal means w.r.t. β.
#'
#' For GLMs with canonical links, the Jacobian is computed analytically:
#' `J[j,] = (1/Σw) Σ w_i dμ/dη(η_i) × X_{j,i}`. This is exact and avoids
#' the O(p) model refits that `numDeriv` requires. For non-GLMs (GAMs, custom
#' models), falls back to `numDeriv::jacobian()` (Richardson extrapolation).
#'
#' @param model A fitted model object (`glm`, `gam`, `glm_weightit`, etc.).
#' @param data_a_list Named list of counterfactual data.tables.
#' @param target_idx Logical vector (length n) flagging target rows.
#' @param V_beta p × p variance–covariance matrix of the model coefficients.
#' @param weights Numeric vector of target-population weights (already
#'   subsetted to target rows) or `NULL` for unweighted.
#'
#' @return A named k × k variance–covariance matrix of marginal means.
#'
#' @noRd
compute_vcov_marginal <- function(
  model,
  data_a_list,
  target_idx,
  V_beta,
  weights = NULL
) {
  int_names <- names(data_a_list)
  beta_hat <- stats::coef(model)

  # Convert intervention datasets to data.frames once (predict.glm requires
  # data.frame, not data.table, for reliable model.matrix dispatch).
  data_a_frames <- lapply(data_a_list, function(da) {
    as.data.frame(da)[target_idx, , drop = FALSE]
  })

  # J is a k × p Jacobian: J[j, ] = ∂μ̂_j / ∂β.
  # For GLMs, compute analytically; for non-GLMs, fall back to numDeriv.
  is_glm <- !is.null(model$family) &&
    !is.null(model$family$mu.eta) &&
    inherits(model, "glm") &&
    !inherits(model, "gam")
  if (is_glm) {
    pred_terms <- stats::delete.response(stats::terms(model))
    # xlev preserves original factor levels so model.matrix() works even when
    # counterfactual data has a single treatment level (e.g. static("low")).
    xlev <- model$xlevels
    J <- t(vapply(
      data_a_frames,
      function(df) {
        X_j <- stats::model.matrix(pred_terms, data = df, xlev = xlev)
        eta_j <- as.numeric(X_j %*% beta_hat)
        mu_eta_j <- model$family$mu.eta(eta_j)
        if (!is.null(weights)) {
          as.numeric(crossprod(X_j, weights * mu_eta_j)) / sum(weights)
        } else {
          as.numeric(crossprod(X_j, mu_eta_j)) / nrow(X_j)
        }
      },
      numeric(length(beta_hat))
    ))
  } else {
    J <- numDeriv::jacobian(
      func = function(beta) {
        model_tmp <- model
        model_tmp$coefficients <- beta
        vapply(
          data_a_frames,
          function(df) {
            preds <- predict(model_tmp, newdata = df, type = "response")
            maybe_weighted_mean(preds, weights)
          },
          numeric(1)
        )
      },
      x = beta_hat
    )
  }

  # Propagate: Var(μ̂) = J V_β Jᵀ  (multivariate delta method).
  vcov_mat <- J %*% V_beta %*% t(J)
  rownames(vcov_mat) <- int_names
  colnames(vcov_mat) <- int_names
  vcov_mat
}
