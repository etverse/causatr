#' ICE branch of variance_if()
#'
#' @description
#' Computes the per-individual IF for chained ICE g-computation by
#' iterating `correct_model()` once per outcome model in the chain. The
#' sensitivity vector `d` is propagated forward (model 0 -> model K),
#' even though the models themselves are fit backward in time
#' (model K -> model 0) -- both directions correspond to the back-
#' substitution of the block-triangular bread for the stacked
#' \eqn{K{+}1}-model M-estimation system (vignette Section 5.4; D2).
#'
#' @param fit A `causatr_fit` of type `"longitudinal"`.
#' @param ice_results Named list of `ice_iterate()` results, one per
#'   intervention.
#' @param target_within_first Logical vector over first-time individuals
#'   flagging the target population.
#'
#' @noRd
variance_if_ice <- function(
  fit,
  ice_results,
  target_within_first,
  cluster_vec = NULL
) {
  int_names <- names(ice_results)
  data <- fit$data
  first_time <- fit$details$time_points[1]
  rows_first <- data[[fit$time]] == first_time
  n <- sum(rows_first)

  # Per-individual IFs in the ICE branch are indexed by unique id at
  # the first time point (one IF per person). `cluster_vec` is supplied
  # at full person-period length `nrow(fit$data)`; subset to first-time
  # rows and align by id order so it matches the IF vector ordering
  # produced by `variance_if_ice_one()`. If multiple person-period rows
  # for the same id disagree on cluster, the first-time row wins — the
  # cluster is expected to be an id-level attribute anyway.
  cluster_first <- if (is.null(cluster_vec)) {
    NULL
  } else {
    cluster_vec[rows_first]
  }

  IF_list <- lapply(ice_results, function(res) {
    variance_if_ice_one(fit, res, target_within_first)
  })

  vcov_from_if(IF_list, n, int_names, cluster = cluster_first)
}


#' Per-individual IF for one ICE intervention
#'
#' @description
#' Workhorse called once per intervention by `variance_if_ice()`. Assembles
#' the IF as
#' \deqn{\mathrm{IF}_i = \frac{n}{n_t}\,t_i\,(\tilde Y_{0,i} - \hat\mu)
#'       + \sum_{k=0}^{K} n_k\,d_{k,i}\,r^{\mathrm{score}}_{k,i}}
#' by looping over the \eqn{K{+}1} outcome models. At each step, the
#' previous step's per-individual sensitivity `d_{k-1,j}` weights the
#' construction of `g_k`, then `correct_model(model_k, g_k, ...)` returns
#' both the model-k correction and the new `d_k` for the next iteration.
#'
#' Mirrors the closed-form back-substitution of the stacked-EE bread
#' from vignette Section 5.4-5.6.
#'
#' @noRd
variance_if_ice_one <- function(fit, ice_result, target) {
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

  first_time <- time_points[1]
  rows_first <- data[[time_col]] == first_time
  all_ids <- as.character(data[rows_first][[id_col]])
  n <- length(all_ids)
  id_to_idx <- stats::setNames(seq_len(n), all_ids)

  # `pseudo_final` may carry NA on first-time rows that were dropped
  # during the backward ICE iteration (e.g. intermediate covariates
  # missing at a later time point). Indexing by `target` rather than
  # multiplying the full-length vector by a zero mask avoids `0 * NA =
  # NA` propagating into the IF.
  #
  # B7 (2026-04-15 review): unify the weighted and unweighted IF on a
  # single formula. Previously the weighted branch used
  #   IF_i = n * (w_i / sum_w_target) * (Y_i - mu_hat)
  # and the unweighted branch used
  #   IF_i = (n / n_target) * (Y_i - mu_hat)
  # which agree only when sum(w) == n_target. For arbitrary external
  # weights they drift, and the Channel-2 cross-term (which uses the
  # n-scaled gradient d_i) is mis-scaled relative to Channel 1.
  # Setting w_i = 1 in the weighted form recovers the unweighted case
  # exactly (sum_w_target = n_target, so n/n_target drops out), so
  # there is one formula to reason about and one to truth-test.
  ext_w <- details$weights
  has_weights <- !is.null(ext_w)
  if (has_weights) {
    w_first <- ext_w[rows_first]
    w_t <- w_first[target]
  } else {
    w_t <- rep(1, sum(target))
  }
  sum_w_target <- sum(w_t)
  # \hat{mu} = (1/sum_w) sum_{target} w_i \tilde{Y}_{0,i}: Hajek-weighted
  # average of the backward-iterated pseudo-outcome over the target population.
  mu_hat <- sum(w_t * pseudo_final[target]) / sum_w_target
  # Channel 1: IF_i = n * (w_i / sum_w) * (\tilde{Y}_{0,i} - \hat{mu}).
  # Scaled by n so vcov_from_if(n) gives the sandwich (1/n^2) sum IF_i^2.
  IF_vec <- numeric(n)
  IF_vec[target] <- n *
    (w_t / sum_w_target) *
    (pseudo_final[target] - mu_hat)

  # Per-time-step id-to-weight lookup for the cascade gradient.
  # At step k > 1, g_k = sum_{j: prev fit} w_{k-1,j} d_{k-1,j} X*_{k,j} \mu'_{k,j}.
  # The w_{k-1} factor arises from A_{k-1,k} = (1/n) sum \partial s_{k-1}/\partial \beta_k
  # in the block-triangular bread (see variance-theory vignette Section 5.4).
  # Omitting it silently under-scales the off-diagonal correction by ext_w,
  # causing ~2x SE underestimation under heterogeneous external weights.
  w_at_step <- if (has_weights) {
    lapply(seq_len(n_times), function(k) {
      rows_k <- data[[time_col]] == time_points[k]
      ids_k <- as.character(data[rows_k][[id_col]])
      stats::setNames(ext_w[rows_k], ids_k)
    })
  } else {
    NULL
  }

  d_vec <- rep(0, n)

  for (step_i in seq_len(n_times)) {
    model_k <- models[[step_i]]
    if (is.null(model_k)) {
      next
    }

    current_time <- time_points[step_i]
    fit_ids_k <- fit_ids_list[[step_i]]
    if (length(fit_ids_k) == 0L) {
      next
    }

    # Subset the (already intervention-modified) longitudinal data to
    # the rows at the current time step, and record which individuals
    # are present (some may have been censored out before this step).
    intervention <- ice_result$intervention
    is_stoch <- has_stochastic_component(intervention)
    rows_iv_current <- data_iv[[time_col]] == current_time
    iv_data_current <- data_iv[rows_iv_current]
    iv_ids_current <- as.character(iv_data_current[[id_col]])

    if (step_i == 1L) {
      # Step 1 (earliest time): g_1 = (1/sum_w) sum_{target} w_i X*_{1,i} \mu'(\eta*_{1,i}).
      # This is \partial \hat{mu} / \partial \beta_1, the gradient of the ICE
      # estimand w.r.t. the first outcome model's parameters.
      target_in_iv <- match(all_ids[target], iv_ids_current)
      valid_target <- !is.na(target_in_iv)
      target_in_iv <- target_in_iv[valid_target]
      target_w <- w_t[valid_target]

      if (is_stoch) {
        n_mc_ice <- get_stochastic_n_mc(intervention)
        p_coef <- length(stats::coef(model_k))
        g_k_sum <- numeric(p_coef)
        for (m in seq_len(n_mc_ice)) {
          dv_m <- ice_apply_intervention_long(
            data,
            fit$treatment,
            intervention,
            id_col,
            time_col
          )
          iv_cur_m <- dv_m[dv_m[[time_col]] == current_time]
          X_m <- iv_design_matrix(model_k, iv_cur_m)
          eta_m <- as.numeric(X_m %*% stats::coef(model_k))
          mu_eta_m <- model_k$family$mu.eta(eta_m)
          ids_m <- as.character(iv_cur_m[[id_col]])
          tgt_in_m <- match(all_ids[target], ids_m)
          vt_m <- !is.na(tgt_in_m)
          tgt_in_m <- tgt_in_m[vt_m]
          # target_w is length sum(valid_target); vt_m is length
          # sum(target). Index into the already-subsetted target_w
          # via vt_m[valid_target] so both operands align.
          vt_in_valid <- vt_m[valid_target]
          g_k_sum <- g_k_sum +
            as.numeric(
              crossprod(
                X_m[tgt_in_m, , drop = FALSE],
                target_w[vt_in_valid] * mu_eta_m[tgt_in_m]
              )
            )
        }
        g_k <- g_k_sum / (n_mc_ice * sum_w_target)
      } else {
        X_star_k <- iv_design_matrix(model_k, iv_data_current)
        eta_star <- as.numeric(
          X_star_k %*% stats::coef(model_k)
        )
        mu_eta_star <- model_k$family$mu.eta(eta_star)
        # g_k = (1/sum_w) X*^T diag(w_target) \mu'(\eta*): the raw gradient
        # before apply_model_correction scales it through the model bread.
        g_k <- as.numeric(
          crossprod(
            X_star_k[target_in_iv, , drop = FALSE],
            target_w * mu_eta_star[target_in_iv]
          )
        ) /
          sum_w_target
      }
    } else {
      # Step k > 1: chain rule through the ICE recursion.
      # g_k = sum_{j: prev fit} w_{k-1,j} d_{k-1,j} X*_{k,j} \mu'_{k,j}.
      # d_{k-1,j} is the sensitivity of \hat{mu} to Q_{k-1,j}; it comes
      # from the previous correct_model() call. w_{k-1,j} is the external
      # weight at the previous time step (see w_at_step lookup above).
      prev_fit_ids <- fit_ids_list[[step_i - 1L]]
      idx_in_all <- id_to_idx[prev_fit_ids]
      rows_in_iv <- match(prev_fit_ids, iv_ids_current)
      keep <- !is.na(idx_in_all) & !is.na(rows_in_iv)

      if (any(keep)) {
        d_prev <- d_vec[idx_in_all[keep]]
        w_prev <- if (has_weights) {
          unname(w_at_step[[step_i - 1L]][prev_fit_ids[keep]])
        } else {
          rep(1, sum(keep))
        }

        if (is_stoch) {
          n_mc_ice <- get_stochastic_n_mc(intervention)
          p_coef <- length(stats::coef(model_k))
          g_k_sum <- numeric(p_coef)
          for (m in seq_len(n_mc_ice)) {
            dv_m <- ice_apply_intervention_long(
              data,
              fit$treatment,
              intervention,
              id_col,
              time_col
            )
            iv_cur_m <- dv_m[dv_m[[time_col]] == current_time]
            X_m <- iv_design_matrix(model_k, iv_cur_m)
            eta_m <- as.numeric(X_m %*% stats::coef(model_k))
            mu_eta_m <- model_k$family$mu.eta(eta_m)
            ids_m <- as.character(iv_cur_m[[id_col]])
            riv_m <- match(prev_fit_ids, ids_m)
            keep_m <- !is.na(idx_in_all) & !is.na(riv_m)
            # keep_m_in_keep: which of the already-kept indices also
            # survive in keep_m (both are length P, but w_prev/d_prev
            # are length sum(keep), so index into them via keep[keep_m]).
            keep_m_in_keep <- keep_m[keep]
            if (any(keep_m)) {
              wg_m <- w_prev[keep_m_in_keep] *
                d_prev[keep_m_in_keep] *
                mu_eta_m[riv_m[keep_m]]
              g_k_sum <- g_k_sum +
                as.numeric(
                  crossprod(
                    X_m[riv_m[keep_m], , drop = FALSE],
                    wg_m
                  )
                )
            }
          }
          g_k <- g_k_sum / n_mc_ice
        } else {
          X_star_k <- iv_design_matrix(model_k, iv_data_current)
          eta_star <- as.numeric(
            X_star_k %*% stats::coef(model_k)
          )
          mu_eta_star <- model_k$family$mu.eta(eta_star)
          weights_g <- w_prev * d_prev * mu_eta_star[rows_in_iv[keep]]
          g_k <- as.numeric(
            crossprod(
              X_star_k[rows_in_iv[keep], , drop = FALSE],
              weights_g
            )
          )
        }
      } else {
        p_coef <- length(stats::coef(model_k))
        g_k <- rep(0, p_coef)
      }
    }

    fit_id_idx <- id_to_idx[fit_ids_k]
    na_act_k <- model_k$na.action
    if (!is.null(na_act_k)) fit_id_idx <- fit_id_idx[-na_act_k]
    # correct_model() computes two outputs:
    #   $correction: per-individual model-k IF contribution (n-scaled)
    #   $d:          updated sensitivity vector d_k = d_{k-1} * (\partial Q_{k-1}/\partial Q_k)
    # The correction is added to IF_vec; d_vec is passed to the next step.
    res <- correct_model(model_k, g_k, fit_id_idx, n)
    IF_vec <- IF_vec + res$correction
    d_vec <- res$d
  }

  IF_vec
}
