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
    # g-computation: standard Huber–White sandwich.
    V_beta <- sandwich::sandwich(model)
  }

  # Propagate parameter uncertainty: J V_β Jᵀ.
  compute_vcov_marginal(model, data_a_list, target_idx, V_beta)
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
  n_int <- length(int_names)

  # Number of individuals (= rows at first time).
  data <- fit$data
  first_time <- fit$details$time_points[1]
  rows_first <- data[[fit$time]] == first_time
  n <- sum(rows_first)

  # Compute per-individual influence functions for each intervention.
  IF_list <- lapply(ice_results, function(res) {
    ice_compute_influence(fit, res, target_within_first)
  })

  # Vcov from cross-products: Cov(μ̂_a, μ̂_b) = (1/n²) Σ IF_a,i × IF_b,i.
  vcov_mat <- matrix(0, n_int, n_int)
  for (j in seq_len(n_int)) {
    for (k in j:n_int) {
      vcov_mat[j, k] <- sum(IF_list[[j]] * IF_list[[k]]) / n^2
      if (j != k) vcov_mat[k, j] <- vcov_mat[j, k]
    }
  }

  vcov_mat
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
  n_target <- sum(target)
  mu_hat <- mean(pseudo_final[target], na.rm = TRUE)

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

  # Direct contribution: target_i × (Ŷ*_0,i - μ̂).
  # Non-target individuals contribute 0 (they are not in the average).
  IF_vec <- ifelse(target, pseudo_final - mu_hat, 0)

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

    # Extract model components.
    # X_k: design matrix for fitting observations (n_fit × p_k).
    # w_k: IWLS working weights = (dμ/dη)² / V(μ).
    # r_k: response residuals (Y_k - μ̂_k) for score contributions.
    X_k <- stats::model.matrix(model_k)
    p_k <- ncol(X_k)
    fitted_k <- stats::fitted(model_k)
    residuals_k <- stats::residuals(model_k, type = "response")

    eta_k <- model_k$linear.predictors
    mu_eta_k <- model_k$family$mu.eta(eta_k)
    var_k <- model_k$family$variance(fitted_k)
    w_k <- mu_eta_k^2 / var_k

    # Weighted Hessian (X^T W X) and its inverse.
    XtWX <- crossprod(X_k, X_k * w_k)
    XtWX_inv <- tryCatch(
      solve(XtWX),
      error = function(e) MASS::ginv(XtWX)
    )

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
      # Model 0 (first time): μ̂ = (1/n_t) Σ_{target} predict(model_0, x^*).
      # g_0 = (1/n_t) Σ_{j ∈ target} dμ/dη × X^*_{0,j}.
      target_in_iv <- match(all_ids[target], iv_ids_current)
      g_k <- as.numeric(
        crossprod(
          X_star_k[target_in_iv, , drop = FALSE],
          mu_eta_star[target_in_iv]
        )
      ) /
        n_target
    } else {
      # Model k > 0: β_k affects μ̂ through pseudo-outcomes of model k-1.
      # g_k = Σ_{j ∈ model_{k-1}} d_{k-1,j} × dμ/dη × X^*_{k,j}.
      prev_fit_ids <- fit_ids_list[[step_i - 1L]]

      g_k <- rep(0, p_k)
      for (pid in prev_fit_ids) {
        idx_in_all <- id_to_idx[pid]
        if (is.na(idx_in_all) || d_vec[idx_in_all] == 0) {
          next
        }
        row_in_iv <- match(pid, iv_ids_current)
        if (is.na(row_in_iv)) {
          next
        }
        g_k <- g_k +
          d_vec[idx_in_all] *
            mu_eta_star[row_in_iv] *
            X_star_k[row_in_iv, ]
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
#' The Jacobian is computed numerically via `numDeriv::jacobian()`, which
#' works for any model type (GLM, GAM, glm_weightit, etc.) without needing
#' analytical gradient formulae.
#'
#' @param model A fitted model object (`glm`, `gam`, `glm_weightit`, etc.).
#' @param data_a_list Named list of counterfactual data.tables.
#' @param target_idx Logical vector (length n) flagging target rows.
#' @param V_beta p × p variance–covariance matrix of the model coefficients.
#'
#' @return A named k × k variance–covariance matrix of marginal means.
#'
#' @noRd
compute_vcov_marginal <- function(model, data_a_list, target_idx, V_beta) {
  int_names <- names(data_a_list)
  beta_hat <- stats::coef(model)

  # Convert intervention datasets to data.frames once (predict.glm requires
  # data.frame, not data.table, for reliable model.matrix dispatch).
  data_a_frames <- lapply(data_a_list, function(da) {
    as.data.frame(da)[target_idx, , drop = FALSE]
  })

  # J is a k × p Jacobian: J[j, ] = ∂μ̂_j / ∂β.
  # Each marginal mean μ̂_j = mean(predict(model, newdata_j)) is a smooth
  # function of β, so numDeriv differentiates it via Richardson extrapolation.
  J <- numDeriv::jacobian(
    func = function(beta) {
      model_tmp <- model
      model_tmp$coefficients <- beta
      vapply(
        data_a_frames,
        function(df) {
          mean(predict(model_tmp, newdata = df, type = "response"))
        },
        numeric(1)
      )
    },
    x = beta_hat
  )

  # Propagate: Var(μ̂) = J V_β Jᵀ  (multivariate delta method).
  vcov_mat <- J %*% V_beta %*% t(J)
  rownames(vcov_mat) <- int_names
  colnames(vcov_mat) <- int_names
  vcov_mat
}
