#' Sampling-model correction for IPW transport sandwich
#'
#' @description
#' Adds the sampling-model block to the stacked sandwich for
#' IPW transport. The sampling model \eqn{P(S = 1 \mid L; \gamma)}
#' produces weights \eqn{w_S(\gamma)} that enter the MSM score. The
#' cross-derivative \eqn{A_{\beta\gamma}} captures how the MSM score
#' depends on \eqn{\gamma}; the sampling-model's own score
#' \eqn{\psi_{\gamma,i}} propagates uncertainty in \eqn{\hat\gamma}
#' into the sandwich.
#'
#' Same pattern as `compute_ipw_ipcw_correction()`:
#' 1. Decompose MSM prior weights: `other_w = pw / w_S_hat`
#' 2. Build `phi_bar_samp(gamma)` closure
#' 3. Cross-derivative via `numDeriv::jacobian`
#' 4. Sampling model correction via `apply_model_correction`
#' 5. Subset to `fit_rows`
#'
#' @param fit A `causatr_fit` with `details$transport == TRUE`.
#' @param msm_model Fitted weighted MSM for one intervention.
#' @param J Marginal-mean Jacobian.
#' @param fit_rows Logical vector (length `nrow(fit$data)`).
#' @param n_sub Integer. Number of MSM fit rows.
#'
#' @return Numeric vector of length `n_sub`.
#' @noRd
compute_ipw_sampling_correction <- function(
  fit,
  msm_model,
  J,
  fit_rows,
  n_sub
) {
  samp_model <- fit$details$sampling_model
  data <- fit$data
  n_total <- nrow(data)

  msm_prep <- prepare_model_if(msm_model, seq_len(n_sub), n_sub)
  msm_res <- apply_model_correction(msm_prep, J)
  h_msm <- n_sub * msm_res$h

  beta_hat <- coef_clean(msm_model)
  X_msm <- msm_prep$X_fit
  fam <- msm_model$family
  eta_msm <- as.numeric(X_msm %*% beta_hat)
  mu_msm <- fam$linkinv(eta_msm)
  mu_eta_msm <- fam$mu.eta(eta_msm)
  var_mu_msm <- fam$variance(mu_msm)
  y_msm <- stats::model.response(stats::model.frame(msm_model))
  r_msm <- y_msm - mu_msm

  # Decompose MSM prior weights: total_w = ext_w * DR_w * w_S.
  # Hold everything except w_S fixed for the gamma-varying closure.
  pw <- msm_model$prior.weights
  if (is.null(pw)) {
    pw <- rep(1, n_sub)
  }

  w_S_fit <- compute_sampling_weights(
    samp_model,
    data,
    fit$target,
    fit$target_subset
  )[fit_rows]
  other_w <- ifelse(w_S_fit > 0, pw / w_S_fit, 0)

  samp_wfn <- make_sampling_weight_fn(
    samp_model,
    data,
    fit$target,
    fit$target_subset
  )

  phi_bar_samp <- function(gamma) {
    w_S_full <- samp_wfn(gamma)
    w_S_sub <- w_S_full[fit_rows]
    s_per_i <- other_w * w_S_sub * mu_eta_msm * r_msm / var_mu_msm
    as.numeric(crossprod(X_msm, s_per_i)) / n_sub
  }

  gamma_hat <- samp_model$gamma_hat
  A_beta_gamma <- -numDeriv::jacobian(phi_bar_samp, x = gamma_hat)
  g_samp <- as.numeric(crossprod(A_beta_gamma, h_msm))

  samp_prep <- prepare_model_if(
    samp_model$model,
    which(samp_model$fit_rows),
    n_total
  )
  samp_res <- apply_model_correction(samp_prep, g_samp)
  samp_res$correction[fit_rows]
}


#' Sampling-model correction for longitudinal IPW transport sandwich
#'
#' @description
#' Longitudinal companion to `compute_ipw_sampling_correction()`. Adds the
#' sampling-model block to the id-level stacked sandwich for longitudinal
#' IPW transport. The sampling model \eqn{P(S = 1 \mid L; \gamma)} was fit
#' on cross-sectional (first-period) data, so the correction lives at id
#' scale (length `n_id`), consistent with the rest of the longitudinal IF.
#'
#' Steps mirror the point-treatment correction:
#' 1. Decompose MSM prior weights: `other_w = total_w / w_S_hat`
#' 2. Build `phi_bar_samp(gamma)` on final-period rows
#' 3. Cross-derivative via `numDeriv::jacobian`
#' 4. Sampling-model correction at id scale via `apply_model_correction`
#'
#' @param fit A `causatr_fit` with `details$transport == TRUE`.
#' @param msm_model Fitted weighted final-period MSM for one intervention.
#' @param J Marginal-mean Jacobian (scalar for `Y ~ 1`).
#' @param data_final data.table of final-period rows (all study rows).
#' @param n_final Integer. Number of final-period rows.
#' @param n_id Integer. Number of unique individuals.
#'
#' @return Numeric vector of length `n_id`.
#'
#' @noRd
compute_ipw_sampling_correction_longitudinal <- function(
  fit,
  msm_model,
  J,
  data_final,
  n_final,
  n_id
) {
  samp_model <- fit$details$sampling_model

  msm_prep <- prepare_model_if(msm_model, seq_len(n_final), n_final)
  msm_res <- apply_model_correction(msm_prep, J)
  h_msm <- n_final * msm_res$h

  beta_hat <- coef_clean(msm_model)
  X_msm <- msm_prep$X_fit
  fam <- msm_model$family
  eta_msm <- as.numeric(X_msm %*% beta_hat)
  mu_msm <- fam$linkinv(eta_msm)
  mu_eta_msm <- fam$mu.eta(eta_msm)
  var_mu_msm <- fam$variance(mu_msm)
  y_msm <- stats::model.response(stats::model.frame(msm_model))
  r_msm <- y_msm - mu_msm

  # Decompose MSM prior weights: total_w = w_A * w_S.
  # Hold w_A fixed and vary w_S with gamma.
  pw <- msm_model$prior.weights
  if (is.null(pw)) {
    pw <- rep(1, n_final)
  }

  w_S_final <- compute_sampling_weights(
    samp_model,
    data_final,
    fit$target,
    fit$target_subset
  )
  other_w <- ifelse(w_S_final > 0, pw / w_S_final, 0)

  samp_wfn_final <- make_sampling_weight_fn(
    samp_model,
    data_final,
    fit$target,
    fit$target_subset
  )

  phi_bar_samp <- function(gamma) {
    w_S <- samp_wfn_final(gamma)
    s_per_i <- other_w * w_S * mu_eta_msm * r_msm / var_mu_msm
    as.numeric(crossprod(X_msm, s_per_i)) / n_final
  }

  gamma_hat <- samp_model$gamma_hat
  A_beta_gamma <- -numDeriv::jacobian(phi_bar_samp, x = gamma_hat)
  g_samp <- as.numeric(crossprod(A_beta_gamma, h_msm))

  # The sampling model was fit on first-period (cross-sectional) data,
  # so fit_rows is relative to that data, which has length n_id.
  # n_total = n_id ensures the correction is on the id scale.
  samp_prep <- prepare_model_if(
    samp_model$model,
    which(samp_model$fit_rows),
    n_id
  )
  samp_res <- apply_model_correction(samp_prep, g_samp)
  samp_res$correction
}
