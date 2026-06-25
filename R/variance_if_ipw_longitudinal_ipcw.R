#' Per-period censoring (IPCW) correction for the longitudinal IPW sandwich
#'
#' @description
#' Longitudinal companion to `compute_ipw_ipcw_correction()` (point) and the
#' structural twin of `compute_ipw_sampling_correction_longitudinal()`
#' (transport). Adds the censoring-model block to the id-level stacked sandwich
#' when `details$ipcw == TRUE`.
#'
#' Under built-in IPCW the final-period MSM is fit weighted by
#' \eqn{w_i = \mathrm{DR}_i \cdot w^{\mathrm{ext}}_i \cdot w^{C}_i(\gamma)},
#' where the censoring weight \eqn{w^{C}_i(\gamma)} is the stabilized IPCW
#' weight produced by the per-period censoring models. Treating
#' \eqn{w^{C}(\hat\gamma)} as a fixed external weight (the previous behaviour)
#' yields the conservative known-weights variance. Propagating the censoring
#' model's estimation uncertainty subtracts the projection of the MSM score
#' onto the censoring score, recovering the Robins-Rotnitzky-Zhao (1994)
#' efficiency gain. For longitudinal IPW this cross-term is **not** negligible
#' (unlike the doubly-robust AIPW case, where double-robust orthogonality makes
#' it near-zero): when censoring depends on treatment, the treated arm's
#' density-ratio reweighting concentrates on heavily IPCW-upweighted rows, so
#' the cross-term materially narrows that arm's SE.
#'
#' @details
#' The correction reuses the same M-estimation primitives as the propensity
#' correction in `compute_ipw_if_self_contained_long_one()`:
#'
#' \enumerate{
#'   \item Build the MSM bread \eqn{h = A_{\beta\beta}^{-1} J} and the MSM score
#'     ingredients (residual, link derivative, variance) at \eqn{\hat\beta} on
#'     the final-period rows.
#'   \item Decompose the MSM prior weight into its \eqn{\gamma}-free part
#'     \eqn{\mathrm{other}_i = w_i / w^{C}_i(\hat\gamma)}, so the closure
#'     \eqn{\bar\Phi_\beta(\gamma)} varies only the censoring weight as
#'     \eqn{\gamma} moves. `make_ipcw_weight_fn_longitudinal()`'s full-length
#'     weight closure supplies \eqn{w^{C}(\gamma)} and reproduces
#'     `details$ipcw_weights` at \eqn{\hat\gamma}.
#'   \item Form the cross-derivative
#'     \eqn{A_{\beta\gamma} = -\partial\bar\Phi_\beta/\partial\gamma} over the
#'     full stacked \eqn{\gamma} via `numDeriv::jacobian()`, then project onto
#'     the censoring-parameter space: \eqn{g = A_{\beta\gamma}^{\mathsf T} h}.
#'   \item Apply each per-period censoring model's correction
#'     (`apply_model_correction()`) on its period-row subset and project to the
#'     id ordering, exactly as the per-period propensity blocks do. Periods
#'     whose \eqn{\gamma} does not enter the final-period weight contribute a
#'     numerically-zero cross-derivative and hence a zero correction.
#' }
#'
#' The returned id-level vector is **subtracted** from the per-intervention IF
#' (same sign convention as the propensity, transport, and point-IPCW
#' corrections).
#'
#' @param fit A `causatr_fit` of estimator `"ipw"`, type `"longitudinal"`, with
#'   `details$ipcw == TRUE` and per-period `details$censoring_models`.
#' @param msm_model The fitted weighted final-period MSM for one intervention.
#' @param J Numeric `p`-vector. The marginal-mean Jacobian
#'   \eqn{\partial\hat\mu/\partial\beta} (scalar for `Y ~ 1`; longer with
#'   effect modification).
#' @param ids_first Character vector of the canonical first-time-point id
#'   ordering (length `n_id`).
#' @param n_final Integer. Number of final-period MSM rows.
#' @param n_id Integer. Number of unique individuals.
#'
#' @returns Numeric vector of length `n_id`, the id-level censoring-model
#'   correction to subtract from the intervention's influence function. All
#'   zeros when no period carries a censoring model.
#'
#' @seealso `compute_ipw_ipcw_correction()` (point analog),
#'   `compute_ipw_sampling_correction_longitudinal()` (transport twin),
#'   `make_ipcw_weight_fn_longitudinal()` (the shared per-period γ block).
#' @family variance
#' @noRd
compute_ipw_ipcw_correction_longitudinal <- function(
  fit,
  msm_model,
  J,
  ids_first,
  n_final,
  n_id
) {
  data <- fit$data
  id_col <- fit$id

  # Per-period censoring EE block (shared with the AIPW sandwich). `blocks`
  # carries each period's model, global rows, score rows, and stacked-gamma
  # columns; `weight_full_fn(gamma)` reproduces `details$ipcw_weights` at
  # gamma_hat and is the only gamma-dependent piece of the MSM weight.
  cb <- make_ipcw_weight_fn_longitudinal(fit)
  if (cb$n_gamma == 0L) {
    # No censoring variation in any period: nothing to propagate.
    return(numeric(n_id))
  }

  # ---- MSM bread + score ingredients on final-period rows --------
  # Same construction as compute_ipw_if_self_contained_long_one(): the MSM
  # lives on the n_final final-period rows. h_msm = A_bb^{-1} J in the
  # M-estimation scaling (n_final * the (X'WX)^{-1} J of apply_model_correction).
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

  # Final-period global row indices, aligned to the MSM (data_final) row order.
  # `fit_rows_final` selects the same rows the MSM was fit on (final period,
  # non-missing outcome, target subset), so `which()` is the data_final order.
  final_global_rows <- which(fit$details$fit_rows_final)

  # Decompose the MSM prior weight: total_w = DR_w * ext_w * ipcw_w(gamma).
  # other_w holds the density-ratio and external pieces fixed while gamma
  # varies, so the perturbed score only moves through the censoring weight.
  pw <- msm_model$prior.weights
  if (is.null(pw)) {
    pw <- rep(1, n_final)
  }
  ipcw_at_final_hat <- cb$weight_full_fn(cb$gamma_hat)[final_global_rows]
  # Censored rows carry ipcw_w = 0; they are excluded from the MSM fit rows so
  # this guard only protects against a defensive 0/0, contributing 0 either way.
  other_w <- ifelse(ipcw_at_final_hat > 0, pw / ipcw_at_final_hat, 0)

  # phi_bar(gamma) = (1/n_final) sum_i X_i * other_w_i * ipcw_w_i(gamma) *
  #                  mu'(eta_i) * (Y_i - mu_i) / V(mu_i)
  # At gamma_hat this reproduces the fitted (vanishing) MSM score; the
  # numerical jacobian gives A_{beta,gamma}.
  phi_bar_cens <- function(gamma) {
    w_ipcw_final <- cb$weight_full_fn(gamma)[final_global_rows]
    s_per_i <- other_w * w_ipcw_final * mu_eta_msm * r_msm / var_mu_msm
    as.numeric(crossprod(X_msm, s_per_i)) / n_final
  }

  # A_{beta,gamma} = -d phi_bar / d gamma (negative-Hessian convention,
  # Stefanski & Boos 2002 eq. 8), then g = A_{beta,gamma}^T h_msm projects the
  # MSM sensitivity onto the stacked censoring-parameter space.
  A_beta_gamma <- -numDeriv::jacobian(phi_bar_cens, x = cb$gamma_hat)
  g_cens <- as.numeric(crossprod(A_beta_gamma, h_msm))

  # ---- Per-period censoring corrections, projected to id-level ----
  # Each period's censoring model was fit on its own disjoint row subset, so the
  # stacked censoring bread is block-diagonal across periods (same structure as
  # the per-period propensity blocks). Slice g into each period's columns,
  # apply the model correction on the period subset, and scatter to id order.
  total_cens_correction_id <- numeric(n_id)

  for (b in cb$blocks) {
    g_cens_k <- g_cens[b$gamma_cols]

    n_period_k <- length(b$rows_global)
    # Local fit-row indices within the period subset (1..n_period_k): the rows
    # the censoring GLM actually fit on, in data_k order. `rows_global` is the
    # ascending period-k row order, so match() recovers `which(cm_k$fit_rows)`.
    local_fit_idx <- match(b$score_rows_global, b$rows_global)

    cens_prep_k <- prepare_model_if(b$model, local_fit_idx, n_period_k)
    cens_res_k <- apply_model_correction(cens_prep_k, g_cens_k)

    # Project the period-row correction to the first-time-point id ordering.
    ids_k <- as.character(data[[id_col]][b$rows_global])
    pos_k <- match(ids_k, ids_first)
    correction_id <- numeric(n_id)
    correction_id[pos_k] <- cens_res_k$correction
    # Rescale from the period's own n_total to the id scale (= 1 when every id
    # has a period-k row, e.g. final-period censoring; < 1 under prior dropout).
    if (n_period_k != n_id) {
      correction_id <- correction_id * (n_id / n_period_k)
    }
    total_cens_correction_id <- total_cens_correction_id + correction_id
  }

  total_cens_correction_id
}
