#' Longitudinal IPW branch of variance_if()
#'
#' @description
#' Per-intervention loop over the longitudinal stacked sandwich. The
#' IF lives at the **id** level (one IF per individual), built on the
#' final-period MSM rows: each id's contribution combines Channel 1
#' (final-period prediction minus the marginal mean), the MSM
#' correction (final-period weighted intercept-only Hajek), and the
#' per-period propensity correction (sum over `K` per-period blocks
#' with `apply_model_correction()` -- block-diagonal bread).
#'
#' Same scale convention as `variance_if_ipw()`:
#'   - all channels are `n_id`-scaled before aggregation;
#'   - `vcov_from_if()` divides by `n_id^2` to get the Hajek-mean
#'     variance.
#'
#' Cluster aggregation is **id-level by construction**. Each
#' individual contributes one IF row, so there's no cross-period
#' double-counting; an explicit `cluster_vec` is collapsed to the
#' first-time-point rows (mirroring `variance_if_ice()`'s convention)
#' so cluster-robust aggregation works the same way as ICE.
#'
#' @param fit A `causatr_fit` of estimator `"ipw"`, type `"longitudinal"`.
#' @param bundles Named list of per-intervention bundles built by
#'   `compute_ipw_contrast_longitudinal()`.
#' @param target_within_first Logical vector flagging the target
#'   population at the first-time-point id ordering.
#' @param cluster_vec Optional cluster vector, length `nrow(fit$data)`
#'   (full person-period). Subset to first-time-point rows for id-level
#'   aggregation.
#'
#' @return A `k x k` variance-covariance matrix.
#'
#' @noRd
variance_if_ipw_longitudinal <- function(
  fit,
  bundles,
  target_within_first,
  cluster_vec = NULL,
  trim = 1
) {
  int_names <- names(bundles)
  details <- fit$details
  data <- fit$data
  treatment_models_by_time <- details$treatment_models_by_time
  numerator_models_by_time <- details$numerator_models_by_time
  fit_data_by_time <- details$fit_data_by_time
  time_points <- details$time_points
  n_times <- details$n_times
  id_col <- fit$id
  time_col <- fit$time
  ext_w <- details$weights

  first_time <- time_points[1]
  rows_first <- data[[time_col]] == first_time
  ids_first <- as.character(data[rows_first][[id_col]])
  n_id <- length(ids_first)
  id_to_first_idx <- stats::setNames(seq_len(n_id), ids_first)

  # Subset cluster vector to first-time-point rows for id-level
  # aggregation. Mirrors `variance_if_ice()`'s convention.
  cluster_first <- if (is.null(cluster_vec)) NULL else cluster_vec[rows_first]

  # Final-period rows + ids: where the MSM is fit and where the IF
  # lives. We assume one row per id at the final period (no MAR dropout);
  # under that, n_final == n_id and the id_to_first mapping is invertible.
  fit_rows_final <- details$fit_rows_final
  data_final <- data[fit_rows_final]
  ids_final <- as.character(data_final[[id_col]])
  final_to_first <- id_to_first_idx[ids_final]
  n_final <- length(ids_final)
  ext_w_final <- if (is.null(ext_w)) NULL else ext_w[fit_rows_final]

  IF_list <- lapply(int_names, function(nm) {
    b <- bundles[[nm]]
    msm_model <- b$msm_model
    mu_hat_j <- b$mu_hat
    intervention <- b$intervention

    # ---- Channel 1 (id-level on final-period rows) -----------------
    # Build the id-indexed Channel-1 vector. For each id present at
    # the final period, the contribution is
    #   (n_id / sum_w_target) * w_i * (pred_i - mu_hat_j)   if target
    #   0                                                    otherwise
    # The Hajek-style scaling matches `variance_if_ipw()` and ensures
    # `vcov_from_if(n = n_id)` reproduces the classical sandwich.
    target_final <- b$target_ids_final
    valid_final <- b$valid_final
    preds_final <- b$preds_final

    # Same `w_target_vec` convention as the point-IPW branch: zero off
    # target, identity 1 on target (unweighted) or external weight on
    # target (weighted). Used to compute both Channel 1 and the
    # marginal-mean Jacobian J.
    if (is.null(ext_w_final)) {
      w_target_vec <- rep(1, n_final)
      w_target_vec[!target_final] <- 0
      sum_w_target <- sum(target_final)
    } else {
      w_target_vec <- ext_w_final
      w_target_vec[!target_final] <- 0
      sum_w_target <- sum(ext_w_final[target_final])
    }
    if (sum_w_target <= 0) {
      rlang::abort(
        "variance_if_ipw_longitudinal(): target-population weights sum to 0.",
        class = "causatr_empty_target"
      )
    }

    # Ch1_i = n_id * (w_i / sum_w) * (\hat{mu}^g_i - \hat{mu}^g),
    # same Hajek scaling as the point-IPW branch but denominator is n_id.
    Ch1_final <- n_id * (w_target_vec / sum_w_target) * (preds_final - mu_hat_j)
    Ch1_final[!valid_final] <- 0

    # Map final-period contributions to id-level (length n_id).
    # Under the complete-case assumption `final_to_first` is a bijection,
    # so this is a permutation; the explicit projection makes
    # loss-of-final-period-rows behavior visible if the assumption is
    # ever relaxed.
    Ch1_i <- numeric(n_id)
    Ch1_i[final_to_first] <- Ch1_final

    # ---- Marginal-mean Jacobian J ----------------------------------
    # MSM is `Y ~ 1` -- the design matrix is a single column of 1s on
    # final-period rows. So `iv_design_matrix(msm_model, data_final)`
    # returns an n_final x 1 matrix and the Jacobian is the
    # final-period-target-weighted average of mu_eta:
    #   J = (1 / sum_w_target) * sum_target (X_star_i * mu_eta_i * w_i)
    # Same shape as `compute_ipw_if_self_contained_one()`.
    X_star <- iv_design_matrix(msm_model, data_final)
    beta_hat <- coef_clean(msm_model)
    eta_star <- as.numeric(X_star %*% beta_hat)
    mu_eta_star <- msm_model$family$mu.eta(eta_star)
    J <- as.numeric(crossprod(X_star, w_target_vec * mu_eta_star)) /
      sum_w_target

    # ---- Stacked-alpha closure for cross-derivative ---------------
    mv_closure <- make_weight_fn_longitudinal(
      treatment_models_by_time = treatment_models_by_time,
      fit_data_by_time = fit_data_by_time,
      ids_first = ids_first,
      id_col = id_col,
      intervention = intervention,
      estimand = "ATE",
      numerator_models_by_time = numerator_models_by_time,
      trim = trim
    )

    # The longitudinal weight closure operates at id-level (length n_id),
    # but phi_bar() inside compute_ipw_if_self_contained_long_one() needs
    # weights aligned to final-period MSM rows (length n_final).
    # `[final_to_first]` re-indexes from id-space to final-period-row order
    # without re-sorting: row j of the MSM corresponds to id final_to_first[j].
    base_wfn <- mv_closure$weight_fn
    if (is.null(ext_w_final)) {
      wfn_final <- function(alpha) base_wfn(alpha)[final_to_first]
    } else {
      ext_w_closure <- ext_w_final
      wfn_final <- function(alpha) {
        base_wfn(alpha)[final_to_first] * ext_w_closure
      }
    }

    if_vec <- compute_ipw_if_self_contained_long_one(
      msm_model = msm_model,
      treatment_models_by_time = treatment_models_by_time,
      fit_data_by_time = fit_data_by_time,
      ids_first = ids_first,
      id_col = id_col,
      weight_fn = wfn_final,
      alpha_offsets = mv_closure$offsets,
      alpha_hat_stacked = mv_closure$alpha_hat,
      J = J,
      Ch1_i = Ch1_i,
      n_id = n_id,
      n_final = n_final,
      final_to_first = final_to_first
    )

    if (isTRUE(fit$details$transport)) {
      if_vec <- if_vec -
        compute_ipw_sampling_correction_longitudinal(
          fit = fit,
          msm_model = msm_model,
          J = J,
          data_final = data_final,
          n_final = n_final,
          n_id = n_id
        )
    }

    if_vec
  })
  names(IF_list) <- int_names

  vcov_from_if(IF_list, n_id, int_names, cluster = cluster_first)
}


#' Per-individual IF for one longitudinal IPW intervention
#'
#' @description
#' Longitudinal companion to `compute_ipw_if_self_contained_one()` and
#' `compute_ipw_if_self_contained_mv_one()`. Same shape: builds the
#' MSM correction once, then sums `K` block-diagonal propensity
#' corrections (one per period) plus the Channel-1 vector.
#'
#' The K propensity models are fit on **disjoint** row subsets (one
#' per period), so the bread of the stacked propensity system is
#' block-diagonal -- the same structural property the multivariate-IPW
#' primitive exploits. The cross-derivative
#' \eqn{[A_{\beta\alpha_1}, \ldots, A_{\beta\alpha_K}]} is computed
#' once via `numDeriv::jacobian(phi_bar, alpha_hat_stacked)` on the
#' full stacked closure; we then slice into K column blocks and call
#' `apply_model_correction()` once per period.
#'
#' Each per-period correction lives in the period-k row space
#' (length n_period_k). To aggregate at the id level, we project each
#' period's correction back to the first-time-point id ordering via
#' the `align_idx` map already used by
#' `compute_longitudinal_weights()`.
#'
#' @param msm_model Fitted weighted final-period MSM (`Y ~ 1`).
#' @param treatment_models_by_time Per-period density models.
#' @param fit_data_by_time Per-period data subsets.
#' @param ids_first Character vector of canonical id ordering.
#' @param id_col Character. Id column name.
#' @param weight_fn Closure `alpha -> w(alpha)` that returns a length
#'   `n_final` weight vector aligned to the MSM's row order. Built
#'   from the longitudinal stacked closure with the final-row
#'   projection already baked in.
#' @param alpha_offsets Integer (K+1)-vector of stacked-alpha block
#'   boundaries.
#' @param alpha_hat_stacked Numeric stacked alpha vector.
#' @param J Numeric scalar/vector. Marginal-mean Jacobian.
#' @param Ch1_i Numeric vector of length `n_id`. Channel 1, id-indexed.
#' @param n_id Integer. Number of unique individuals.
#' @param n_final Integer. Number of final-period rows.
#' @param final_to_first Integer mapping final-period rows to first-
#'   time-point ids.
#'
#' @return Numeric vector of length `n_id`.
#'
#' @noRd
compute_ipw_if_self_contained_long_one <- function(
  msm_model,
  treatment_models_by_time,
  fit_data_by_time,
  ids_first,
  id_col,
  weight_fn,
  alpha_offsets,
  alpha_hat_stacked,
  J,
  Ch1_i,
  n_id,
  n_final,
  final_to_first
) {
  # ---- MSM correction (final-period row space) -------------------
  # MSM lives on final-period rows (length n_final). `prepare_model_if`
  # builds the bread on the MSM's `model.matrix(model)` rows; the
  # correction is then `n_final * d * r_score` per final-period row.
  # We project to id space at the end.
  msm_prep <- prepare_model_if(msm_model, seq_len(n_final), n_final)
  msm_res <- apply_model_correction(msm_prep, J)
  n_fit <- nrow(msm_prep$X_fit)

  # Natural-course short-circuit: W_i = \prod_k g_k/f_k = 1 identically,
  # so d(W_i)/d(alpha_k) = 0 and the propensity correction vanishes.
  # `alpha_hat_stacked` is length 0; numDeriv errors on empty alpha,
  # so skip to Channel 1 + MSM correction.
  if (length(alpha_hat_stacked) == 0L) {
    msm_correction_id <- numeric(n_id)
    msm_correction_id[final_to_first] <- msm_res$correction
    if (n_final != n_id) {
      msm_correction_id <- msm_correction_id * (n_id / n_final)
    }
    return(Ch1_i + msm_correction_id)
  }

  # ---- Stacked cross-derivative via numDeriv ---------------------
  # Same `phi_bar(alpha)` shape as the univariate / multivariate
  # IPW primitives: psi_beta_i = X * w(alpha) * (Y - mu) * mu_eta /
  # var_mu, averaged over the n_fit MSM rows.
  beta_hat <- coef_clean(msm_model)
  X_msm <- msm_prep$X_fit
  y_fit <- stats::model.response(stats::model.frame(msm_model))
  fam <- msm_model$family
  eta <- as.numeric(X_msm %*% beta_hat)
  mu <- fam$linkinv(eta)
  mu_eta <- fam$mu.eta(eta)
  var_mu <- fam$variance(mu)
  r_fit <- y_fit - mu

  phi_bar <- function(alpha) {
    w_alpha <- weight_fn(alpha)
    s_per_i <- w_alpha * mu_eta * r_fit / var_mu
    as.numeric(crossprod(X_msm, s_per_i)) / n_fit
  }
  # Negative-Hessian convention (matches univariate / multivariate).
  A_beta_alpha <- -numDeriv::jacobian(phi_bar, x = alpha_hat_stacked)

  # `h_msm_true = A_bb^{-1} J = n_fit * msm_res$h` per the
  # `apply_model_correction()` scaling convention.
  h_msm_true <- n_fit * msm_res$h

  # MSM correction is on final-period row scale (length n_final);
  # project to id-level via the `final_to_first` index map.
  msm_correction_id <- numeric(n_id)
  msm_correction_id[final_to_first] <- msm_res$correction

  # Rescale id-level Channel 1 + MSM correction so the stacked
  # sandwich uses the n_id scale. `apply_model_correction()` returns
  # `n_final * d_fit * r_score` for the MSM piece, but we want
  # `n_id * d_id * r_score`. Multiplying by `n_id / n_final` (= 1
  # under the bijection assumption) keeps the convention exact even
  # if the bijection ever breaks.
  if (n_final != n_id) {
    msm_correction_id <- msm_correction_id * (n_id / n_final)
  }

  # ---- Per-period propensity corrections -------------------------
  # The propensity models are fit on disjoint row subsets (one per
  # period), so the stacked propensity bread is block-diagonal.
  # For univariate: each period contributes one propensity block.
  # For multivariate: each period contributes K_comp sub-blocks (one
  # per treatment component), all fit on the same period-k rows. We
  # iterate over periods, then over components within each period.
  n_periods <- length(treatment_models_by_time)
  is_mv <- inherits(
    treatment_models_by_time[[1L]],
    "causatr_treatment_models"
  )
  total_prop_correction_id <- rep(0, n_id)

  for (k in seq_len(n_periods)) {
    period_idx <- alpha_offsets[k]:(alpha_offsets[k + 1L] - 1L)
    A_block_period <- A_beta_alpha[, period_idx, drop = FALSE]

    tm_k <- treatment_models_by_time[[k]]
    data_k <- fit_data_by_time[[k]]
    ids_k <- as.character(data_k[[id_col]])

    if (is_mv) {
      # Multivariate: tm_k is a causatr_treatment_models list of
      # K_comp per-component models. Sub-slice the period block into
      # per-component blocks and apply corrections independently.
      K_comp <- length(tm_k)
      period_ids <- ids_k[tm_k[[1L]]$fit_rows]
      pos_k <- match(period_ids, ids_first)
      n_period_k <- length(period_ids)

      comp_lens <- vapply(
        tm_k,
        function(m) length(m$alpha_hat),
        integer(1)
      )
      comp_offsets <- c(1L, cumsum(comp_lens) + 1L)

      for (j in seq_len(K_comp)) {
        comp_idx <- comp_offsets[j]:(comp_offsets[j + 1L] - 1L)
        A_block_j <- A_block_period[, comp_idx, drop = FALSE]
        g_prop_j <- as.numeric(crossprod(A_block_j, h_msm_true))

        prop_model_j <- tm_k[[j]]$model
        prop_prep_j <- if (inherits(prop_model_j, "multinom")) {
          prepare_model_if_multinom(
            prop_model_j,
            seq_len(n_period_k),
            n_period_k
          )
        } else {
          prepare_model_if(
            prop_model_j,
            seq_len(n_period_k),
            n_period_k
          )
        }
        prop_res_j <- apply_model_correction(prop_prep_j, g_prop_j)

        correction_id <- numeric(n_id)
        correction_id[pos_k] <- prop_res_j$correction
        if (n_period_k != n_id) {
          correction_id <- correction_id * (n_id / n_period_k)
        }
        total_prop_correction_id <-
          total_prop_correction_id + correction_id
      }
    } else {
      # Univariate: one propensity model per period.
      g_prop_k <- as.numeric(
        crossprod(A_block_period, h_msm_true)
      )

      period_ids <- ids_k[tm_k$fit_rows]
      pos_k <- match(period_ids, ids_first)
      n_period_k <- length(period_ids)

      prop_model_k <- tm_k$model
      prop_prep_k <- if (inherits(prop_model_k, "multinom")) {
        prepare_model_if_multinom(
          prop_model_k,
          seq_len(n_period_k),
          n_period_k
        )
      } else {
        prepare_model_if(
          prop_model_k,
          seq_len(n_period_k),
          n_period_k
        )
      }
      prop_res_k <- apply_model_correction(prop_prep_k, g_prop_k)

      correction_id <- numeric(n_id)
      correction_id[pos_k] <- prop_res_k$correction
      if (n_period_k != n_id) {
        correction_id <- correction_id * (n_id / n_period_k)
      }
      total_prop_correction_id <-
        total_prop_correction_id + correction_id
    }
  }

  # Block-lower-triangular M-estimation:
  #   IF_beta_i = A_bb^{-1}(psi_beta_i - sum_k A_{beta, alpha_k} A_{aa,k}^{-1} psi_alpha_{k,i})
  # i.e. the propensity correction is SUBTRACTED. Same sign convention
  # as the univariate / multivariate IPW primitives.
  Ch1_i + msm_correction_id - total_prop_correction_id
}
