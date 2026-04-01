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
