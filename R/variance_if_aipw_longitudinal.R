#' Longitudinal AIPW sandwich variance (ICE-AIPW)
#'
#' @description
#' Analytic sandwich for the ICE-AIPW estimator (Bang & Robins 2005).
#' Derives from composition of Zivich et al. (2024, *Stat. Med.* 43:5562-5572,
#' arXiv:2306.10976) for the ICE outcome model chain and Shook-Sa et al.
#' (2025, *Biometrics* 81(2):ujaf054, arXiv:2404.16166) for the point-treatment
#' AIPW propensity correction.
#'
#' **DR caveat:** This sandwich SE is consistent only when **both**
#' nuisance models (outcome + propensity) are correctly specified.
#' Under one-model misspecification the AIPW point estimate remains
#' consistent (DR property), but the sandwich SE is not. For DR-valid
#' variance under misspecification, use `ci_method = "bootstrap"`
#' (Rotnitzky et al. 2012, *Biometrika*).
#'
#' @param fit `causatr_fit` with `estimator = "aipw"`,
#'   `type = "longitudinal"`.
#' @param ice_aipw_results Named list of `ice_aipw_iterate()` outputs.
#' @param target_within_first Logical vector over first-time-point
#'   individuals flagging the target population.
#' @param cluster_vec Optional cluster vector for clustered SE.
#'
#' @return A named `k x k` variance-covariance matrix.
#'
#' @noRd
variance_if_aipw_longitudinal <- function(
  fit,
  ice_aipw_results,
  target_within_first,
  cluster_vec = NULL,
  trim = 1
) {
  data <- fit$data
  details <- fit$details
  int_names <- names(ice_aipw_results)
  time_points <- details$time_points
  n_times <- details$n_times
  id_col <- fit$id
  time_col <- fit$time
  ext_w <- details$weights
  treatment_models_by_time <- details$treatment_models_by_time
  fit_data_by_time <- details$fit_data_by_time

  first_time <- time_points[1]
  rows_first <- data[[time_col]] == first_time
  all_ids <- as.character(data[rows_first][[id_col]])
  n <- length(all_ids)
  id_to_idx <- stats::setNames(seq_len(n), all_ids)

  # Unbalanced-panel sandwich caveat. When some individuals are not
  # observed at every period (monotone dropout, censoring row-filter), a
  # period's propensity/outcome models are fit on fewer than `n` rows and
  # the per-period IF contributions are embedded into id-space with the
  # `n / n_period_k` rescale below. A Monte-Carlo study shows this
  # rescaled sandwich underestimates the true SE by ~15% under informative
  # dropout, because dropped ids contribute exactly zero to later-period
  # channels rather than their (unobserved) counterfactual contribution —
  # a limitation of the row-filtering IF that a constant rescale cannot
  # repair. The bootstrap reproduces the truth. Warn and steer users to
  # `ci_method = "bootstrap"` rather than return a silently-low SE.
  n_by_period <- vapply(fit_data_by_time, nrow, integer(1))
  if (any(n_by_period != n)) {
    rlang::warn(
      c(
        "Longitudinal AIPW sandwich SE is unreliable on an unbalanced panel.",
        i = paste0(
          "Some individuals are not observed at every period, so the ",
          "sandwich SE can underestimate the truth by ~15%."
        ),
        i = "Use `ci_method = 'bootstrap'` for valid inference here."
      ),
      class = "causatr_longitudinal_aipw_unbalanced_sandwich",
      .frequency = "once",
      .frequency_id = "causatr_longitudinal_aipw_unbalanced_sandwich"
    )
  }

  cluster_first <- if (is.null(cluster_vec)) {
    NULL
  } else {
    cluster_vec[rows_first]
  }

  IF_list <- lapply(int_names, function(nm) {
    res <- ice_aipw_results[[nm]]
    variance_if_aipw_long_one(
      fit,
      res,
      target_within_first,
      all_ids,
      n,
      id_to_idx,
      rows_first,
      trim = trim
    )
  })
  names(IF_list) <- int_names

  vcov_from_if(IF_list, n, int_names, cluster = cluster_first)
}


#' Per-intervention IF for longitudinal AIPW
#'
#' @description
#' Assembles the influence function for one intervention under the
#' ICE-AIPW estimator. Three channels:
#' \enumerate{
#'   \item Direct: Hajek-scaled deviation of individual pseudo-outcomes.
#'   \item Outcome model corrections: forward sensitivity cascade
#'     (same structure as ICE) with augmented gradient accounting for
#'     the \eqn{-W_k X_{obs} \mu'_{obs}} term from the residual.
#'   \item Propensity model corrections: per-period numerical Jacobian
#'     of the augmented mean w.r.t. stacked propensity parameters.
#' }
#'
#' @param fit `causatr_fit`.
#' @param aipw_result Output from `ice_aipw_iterate()`.
#' @param target Logical vector for target population.
#' @param all_ids Character vector of all individual IDs.
#' @param n Integer. Number of individuals.
#' @param id_to_idx Named integer mapping IDs to indices.
#' @param rows_first Logical vector for first-time-point rows.
#'
#' @return Numeric vector of length `n`.
#'
#' @noRd
variance_if_aipw_long_one <- function(
  fit,
  aipw_result,
  target,
  all_ids,
  n,
  id_to_idx,
  rows_first,
  trim = 1
) {
  data <- fit$data
  details <- fit$details
  models <- aipw_result$models
  data_iv <- aipw_result$data_iv
  fit_ids_list <- aipw_result$fit_ids
  pseudo_final <- aipw_result$pseudo_final
  W_period <- aipw_result$period_weights
  W_cumul <- aipw_result$cumulative_weights
  intervention <- aipw_result$intervention

  time_points <- details$time_points
  n_times <- details$n_times
  id_col <- fit$id
  time_col <- fit$time
  model_fn <- details$model_fn
  treatment_models_by_time <- details$treatment_models_by_time
  fit_data_by_time <- details$fit_data_by_time
  ext_w <- details$weights
  binary_outcome <- is_binary_family(details$family_outcome)

  # ---- Weights and mu_hat -----------------------------------------
  has_weights <- !is.null(ext_w)
  if (has_weights) {
    w_first <- ext_w[rows_first]
    w_t <- w_first[target]
  } else {
    w_t <- rep(1, sum(target))
  }
  sum_w_target <- sum(w_t)
  mu_hat <- sum(w_t * pseudo_final[target]) / sum_w_target

  # ---- Channel 1: direct contribution ----------------------------
  # Ch1_i = n * (w_i / sum_w) * (\tilde{Y}_{0,i} - \hat{mu}) for target ids,
  # where \tilde{Y}_{0,i} is the ICE-AIPW pseudo-outcome propagated backward
  # to time 0 (pseudo_final). Zero for non-target ids.
  IF_vec <- numeric(n)
  IF_vec[target] <- n *
    (w_t / sum_w_target) *
    (pseudo_final[target] - mu_hat)

  # Per-time-step external weight lookup (same as ICE sandwich)
  w_at_step <- if (has_weights) {
    lapply(seq_len(n_times), function(k) {
      rows_k <- data[[time_col]] == time_points[k]
      ids_k <- as.character(data[rows_k][[id_col]])
      stats::setNames(ext_w[rows_k], ids_k)
    })
  } else {
    NULL
  }

  # ---- Channel 2a: outcome model corrections ---------------------
  # Forward sensitivity cascade (ICE-style) with augmented gradient.
  # At each step k the gradient gains an additional term from the
  # AIPW residual: -W_cumul[i,k] * X_obs_i * mu_eta_obs_i.
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

    rows_iv_current <- data_iv[[time_col]] == current_time
    iv_data_current <- data_iv[rows_iv_current]
    iv_ids_current <- as.character(iv_data_current[[id_col]])

    # Observed data at this time step
    rows_obs_current <- data[[time_col]] == current_time
    obs_data_current <- data[rows_obs_current]
    obs_ids_current <- as.character(obs_data_current[[id_col]])

    fam_k <- model_k$family
    beta_k <- coef_clean(model_k)

    if (step_i == 1L) {
      # First time step: gradient from target population
      target_in_iv <- match(all_ids[target], iv_ids_current)
      valid_target <- !is.na(target_in_iv)
      target_in_iv <- target_in_iv[valid_target]
      target_w <- w_t[valid_target]

      # Term (a): intervention predictions d/dbeta m_k(d_k, L)
      X_star_k <- iv_design_matrix(model_k, iv_data_current)
      eta_star <- as.numeric(X_star_k %*% beta_k)
      mu_eta_star <- fam_k$mu.eta(eta_star)

      grad_a <- as.numeric(
        crossprod(
          X_star_k[target_in_iv, , drop = FALSE],
          target_w * mu_eta_star[target_in_iv]
        )
      ) /
        sum_w_target

      # Term (b): gradient from the augmentation residual at step 1.
      # -W_{1,i} * X_obs_i * \mu'(\eta_obs_i) averaged over target.
      # W_period[i, 1] is the per-period (not cumulative) weight,
      # matching the backward recursion's residual structure.
      target_in_obs <- match(
        all_ids[target],
        obs_ids_current
      )
      valid_obs <- !is.na(target_in_obs)
      target_in_obs <- target_in_obs[valid_obs]
      target_w_obs <- w_t[valid_obs]

      X_obs_k <- iv_design_matrix(model_k, obs_data_current)
      eta_obs <- as.numeric(X_obs_k %*% beta_k)
      mu_eta_obs <- fam_k$mu.eta(eta_obs)

      target_ids_obs <- all_ids[target][valid_obs]
      w_period_target <- W_period[
        id_to_idx[target_ids_obs],
        step_i
      ]

      grad_b <- -as.numeric(
        crossprod(
          X_obs_k[target_in_obs, , drop = FALSE],
          target_w_obs * w_period_target * mu_eta_obs[target_in_obs]
        )
      ) /
        sum_w_target

      g_k <- grad_a + grad_b
    } else {
      # Later steps (k > 1): sensitivity d_vec cascades forward from the
      # previous model -- same chain-rule mechanism as variance_if_ice_one().
      # grad_a uses intervention predictions; grad_b uses observed predictions
      # weighted by W_period[i, k] (per-period, not cumulative density ratio).
      prev_fit_ids <- fit_ids_list[[step_i - 1L]]
      idx_in_all <- id_to_idx[prev_fit_ids]
      rows_in_iv <- match(prev_fit_ids, iv_ids_current)
      rows_in_obs <- match(prev_fit_ids, obs_ids_current)
      keep <- !is.na(idx_in_all) &
        !is.na(rows_in_iv) &
        !is.na(rows_in_obs)

      if (any(keep)) {
        d_prev <- d_vec[idx_in_all[keep]]
        w_prev <- if (has_weights) {
          unname(
            w_at_step[[step_i - 1L]][prev_fit_ids[keep]]
          )
        } else {
          rep(1, sum(keep))
        }

        # Intervention gradient (same as ICE cascade)
        X_star_k <- iv_design_matrix(
          model_k,
          iv_data_current
        )
        eta_star <- as.numeric(X_star_k %*% beta_k)
        mu_eta_star <- fam_k$mu.eta(eta_star)

        weights_g_iv <- w_prev * d_prev * mu_eta_star[rows_in_iv[keep]]
        grad_a <- as.numeric(
          crossprod(
            X_star_k[rows_in_iv[keep], , drop = FALSE],
            weights_g_iv
          )
        )

        # Observed gradient with cumulative weights
        X_obs_k <- iv_design_matrix(
          model_k,
          obs_data_current
        )
        eta_obs <- as.numeric(X_obs_k %*% beta_k)
        mu_eta_obs <- fam_k$mu.eta(eta_obs)

        prev_ids_keep <- prev_fit_ids[keep]
        w_period_prev <- W_period[
          id_to_idx[prev_ids_keep],
          step_i
        ]

        weights_g_obs <- w_prev *
          d_prev *
          w_period_prev *
          mu_eta_obs[rows_in_obs[keep]]
        grad_b <- -as.numeric(
          crossprod(
            X_obs_k[rows_in_obs[keep], , drop = FALSE],
            weights_g_obs
          )
        )

        g_k <- grad_a + grad_b
      } else {
        p_coef <- length(beta_k)
        g_k <- rep(0, p_coef)
      }
    }

    fit_id_idx <- id_to_idx[fit_ids_k]
    na_act_k <- model_k$na.action
    if (!is.null(na_act_k)) {
      fit_id_idx <- fit_id_idx[-na_act_k]
    }
    res <- correct_model(model_k, g_k, fit_id_idx, n)
    IF_vec <- IF_vec + res$correction
    d_vec <- res$d
  }

  # ---- Channel 2b: propensity model corrections ------------------
  # Per-period density-ratio corrections via numerical Jacobian on
  # the AIPW augmented mean as a function of the stacked propensity
  # parameters alpha. Same block-diagonal structure as longitudinal
  # IPW: K disjoint per-period propensity models.
  #
  # We build per-period weight sub-closures (same way
  # make_weight_fn_longitudinal does internally) so we can
  # reconstruct per-period cumulative products under perturbed alpha.
  ids_first <- all_ids

  # Multivariate longitudinal AIPW: each period's propensity is a
  # K-component `causatr_treatment_models` chain rather than a single
  # model. The per-period sub-closure is then `make_weight_fn_mv()`
  # (which stacks the K component alphas internally) and the per-period
  # alpha block carries all K component parameters concatenated. Net:
  # the stacked alpha is T×K blocks in period-major order, mirroring
  # `make_weight_fn_longitudinal()`'s MV delegation.
  is_mv <- inherits(
    treatment_models_by_time[[1L]],
    "causatr_treatment_models"
  )

  K <- length(treatment_models_by_time)
  sub_fns <- vector("list", K)
  block_lens <- integer(K)
  alpha_blocks <- vector("list", K)
  align_idx_list <- vector("list", K)

  for (k in seq_len(K)) {
    tm_k <- treatment_models_by_time[[k]]
    data_k <- fit_data_by_time[[k]]
    ids_k <- as.character(data_k[[id_col]])

    if (is_mv) {
      # Reset fit_rows to all-TRUE so the per-period MV closure computes a
      # weight for every row of the (already-filtered) period-k data. We do
      # NOT re-attach a stabilized numerator chain: `fit_aipw()` rejects
      # `stabilize` for every longitudinal AIPW fit upstream
      # (`causatr_stabilize_longitudinal_aipw`), so no period model ever
      # carries a marginal numerator here.
      tms_local <- tm_k
      for (j in seq_along(tms_local)) {
        tms_local[[j]]$fit_rows <- rep(TRUE, nrow(data_k))
      }
      class(tms_local) <- c("causatr_treatment_models", "list")
      mv_closure_k <- make_weight_fn_mv(
        tms_local,
        data_k,
        intervention,
        estimand = "ATE",
        trim = trim
      )
      sub_fns[[k]] <- mv_closure_k$weight_fn
      # The closure spans ALL period-k rows (fit_rows reset above), so the
      # broadcast target is every period-k id, not a per-component fit-row
      # subset. Using a single component's fit_rows here would mis-length the
      # broadcast whenever components drop different rows.
      period_ids <- ids_k
      align_idx_list[[k]] <- match(period_ids, ids_first)
      alpha_blocks[[k]] <- mv_closure_k$alpha_hat
      block_lens[k] <- length(mv_closure_k$alpha_hat)
    } else {
      sub_fns[[k]] <- make_weight_fn(
        treatment_model = tm_k,
        data = data_k,
        intervention = intervention,
        estimand = "ATE",
        trim = trim
      )
      period_ids <- ids_k[tm_k$fit_rows]
      align_idx_list[[k]] <- match(period_ids, ids_first)
      alpha_k <- tm_k$alpha_hat
      alpha_blocks[[k]] <- alpha_k
      block_lens[k] <- length(alpha_k)
    }
  }

  alpha_offsets <- c(1L, cumsum(block_lens) + 1L)
  alpha_hat_stacked <- unlist(alpha_blocks, use.names = FALSE)

  # Skip propensity correction for natural course
  if (length(alpha_hat_stacked) > 0L) {
    # Build the augmented-mean closure for numDeriv. Perturbing alpha changes
    # only the period weights W_new; outcome models are held at their fitted
    # values. The backward recursion re-runs under the perturbed weights so
    # numDeriv captures d(aug_mean)/d(alpha) through the full weight product.
    models_fixed <- models
    data_iv_fixed <- data_iv

    aug_mean <- function(alpha) {
      # Reconstruct per-period weights from perturbed alpha
      W_new <- matrix(1, nrow = n, ncol = K)
      for (kk in seq_len(K)) {
        idx <- alpha_offsets[kk]:(alpha_offsets[kk + 1L] - 1L)
        w_k_raw <- sub_fns[[kk]](alpha[idx])

        # Broadcast period-k weights from the subset of fitted ids to all n ids.
        # Unobserved ids (not in align_idx_list[[kk]]) retain w_k = 1 (no
        # augmentation contribution for missing periods).
        w_k <- rep(1, n)
        w_k[align_idx_list[[kk]]] <- w_k_raw
        W_new[, kk] <- w_k
      }

      # Backward AIPW pseudo-outcome recursion with fixed models
      pseudo_a <- rep(NA_real_, n)

      for (step_i in rev(seq_len(n_times))) {
        model_k <- models_fixed[[step_i]]
        if (is.null(model_k)) {
          next
        }

        current_time <- time_points[step_i]
        rows_iv <- data_iv_fixed[[time_col]] == current_time
        rows_obs <- data[[time_col]] == current_time
        obs_ids <- as.character(data[rows_obs][[id_col]])
        obs_idx <- id_to_idx[obs_ids]

        preds_iv_k <- stats::predict(
          model_k,
          newdata = data_iv_fixed[rows_iv],
          type = "response"
        )
        preds_obs_k <- stats::predict(
          model_k,
          newdata = data[rows_obs],
          type = "response"
        )

        if (step_i == n_times) {
          # At the final period the residual is Y_i - Q_K(A_obs, H_K).
          y_k <- data[rows_obs][[fit$outcome]]
          resid_k <- y_k - preds_obs_k
        } else {
          # At earlier periods the residual is \tilde{Y}_{k+1,i} - Q_k(A_obs, H_k),
          # i.e. the backward pseudo-outcome minus the observed-treatment prediction.
          resid_k <- pseudo_a[obs_idx] - preds_obs_k
        }

        has_prev <- !is.na(resid_k)
        pseudo_k <- preds_iv_k
        pseudo_k[has_prev] <- preds_iv_k[has_prev] +
          W_new[obs_idx[has_prev], step_i] *
            resid_k[has_prev]
        if (binary_outcome) {
          pseudo_k <- pmax(pmin(pseudo_k, 1 - 1e-5), 1e-5)
        }
        pseudo_a[obs_idx] <- pseudo_k
      }

      sum(w_t * pseudo_a[which(target)]) / sum_w_target
    }

    # J_alpha = \partial aug_mean / \partial \alpha (stacked over K periods).
    # Block-diagonal bread means each period's slice J_alpha_k is independent.
    J_alpha <- as.numeric(
      numDeriv::jacobian(aug_mean, x = alpha_hat_stacked)
    )

    # Per-period propensity corrections (block-diagonal bread).
    # Slice the k-th alpha block and pass it as the sensitivity gradient
    # to apply_model_correction, which computes sum_i J_alpha_k^T A_{kk}^{-1} psi_{k,i}.
    for (k in seq_len(K)) {
      idx <- alpha_offsets[k]:(alpha_offsets[k + 1L] - 1L)
      J_alpha_k <- J_alpha[idx]

      tm_k <- treatment_models_by_time[[k]]
      data_k <- fit_data_by_time[[k]]
      ids_k <- as.character(data_k[[id_col]])

      if (is_mv) {
        # Multivariate: the period block is the concatenation of K
        # per-component alpha sub-blocks (one Bernoulli/Gaussian/etc.
        # density per treatment component). The within-period chain rule
        # makes the components fit independently, so the period bread is
        # block-diagonal across components — slice J_alpha_k into K
        # sub-blocks and apply a correction per component model.
        comp_lens <- vapply(
          tm_k,
          function(m) length(m$alpha_hat),
          integer(1)
        )
        comp_offsets <- c(1L, cumsum(comp_lens) + 1L)

        for (j in seq_along(tm_k)) {
          comp_idx <- comp_offsets[j]:(comp_offsets[j + 1L] - 1L)
          J_alpha_j <- J_alpha_k[comp_idx]

          # Each component model carries its own fit_rows (NA handling can
          # differ across components), so align the component score to the
          # rows that component was actually fit on rather than assuming a
          # shared per-period fit-row set.
          comp_ids_j <- ids_k[tm_k[[j]]$fit_rows]
          pos_j <- match(comp_ids_j, ids_first)
          n_period_j <- length(pos_j)

          prop_model_j <- tm_k[[j]]$model
          prop_prep_j <- if (inherits(prop_model_j, "multinom")) {
            prepare_model_if_multinom(
              prop_model_j,
              seq_len(n_period_j),
              n_period_j
            )
          } else {
            prepare_model_if(
              prop_model_j,
              seq_len(n_period_j),
              n_period_j
            )
          }
          prop_res_j <- apply_model_correction(prop_prep_j, J_alpha_j)

          correction_id <- numeric(n)
          correction_id[pos_j] <- prop_res_j$correction
          if (n_period_j != n) {
            correction_id <- correction_id * (n / n_period_j)
          }
          IF_vec <- IF_vec + correction_id
        }
        next
      }

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
      prop_res_k <- apply_model_correction(
        prop_prep_k,
        J_alpha_k
      )

      # Project period-k correction (n_period_k-scaled) to id-space (n-scaled).
      # n/n_period_k = 1 under the complete-case bijection; kept for robustness.
      correction_id <- numeric(n)
      correction_id[pos_k] <- prop_res_k$correction
      if (n_period_k != n) {
        correction_id <- correction_id * (n / n_period_k)
      }
      # Propensity correction added: J_alpha already carries the correct
      # sign from d(mu)/d(alpha), matching the standard M-estimation
      # block-inverse identity (Stefanski & Boos 2002).
      IF_vec <- IF_vec + correction_id
    }
  }

  IF_vec
}
