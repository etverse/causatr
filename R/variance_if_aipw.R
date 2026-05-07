#' AIPW branch of variance_if()
#'
#' @description
#' Computes the sandwich variance for the AIPW estimator via a
#' three-block stacked estimating-equation system:
#' 1. Outcome model score equations (defines \eqn{\beta})
#' 2. Propensity model score equations (defines \eqn{\alpha})
#' 3. AIPW plug-in equation (defines \eqn{\mu})
#'
#' The per-individual IF for \eqn{\hat\mu(g)} is:
#' \deqn{\mathrm{IF}_i = \mathrm{Ch1}_i +
#'   \mathrm{outcome\_correction}_i -
#'   \mathrm{propensity\_correction}_i}
#'
#' where Ch1 is the direct AIPW functional evaluated at individual i
#' minus the estimate, and the corrections propagate nuisance-model
#' parameter uncertainty through the plug-in functional.
#'
#' @noRd
variance_if_aipw <- function(
  fit,
  aipw_bundles,
  target_idx,
  mu_hat,
  aipw_fit_idx,
  aipw_n_total,
  cluster_vec = NULL
) {
  data <- fit$data
  outcome_model <- fit$model
  tm <- fit$details$treatment_model
  is_mv <- isTRUE(fit$details$is_multivariate)
  tms <- fit$details$treatment_models
  propensity_model <- if (is_mv) NULL else tm$model
  fit_rows <- fit$details$fit_rows
  fit_data <- data[fit_rows]
  ext_w <- fit$details$weights
  estimand <- fit$estimand
  int_names <- names(aipw_bundles)
  n_sub <- sum(fit_rows)
  n_total <- aipw_n_total
  is_transport <- isTRUE(fit$details$transport)

  # For transport the IF spans all n_total rows (target + study contribute
  # to different terms). For non-transport all work is on fit_rows = n_sub.
  n_if <- if (is_transport) n_total else n_sub

  # Target population indexing depends on transport mode:
  # - Non-transport: target_fit within fit_rows (ATE: all TRUE; ATT/ATC: subset)
  # - Transport: target_idx spans full data; fit_rows are the study (S=1) rows
  if (is_transport) {
    n_target_raw <- sum(target_idx)
    if (!is.null(ext_w)) {
      sum_w_target <- sum(ext_w[target_idx])
    } else {
      sum_w_target <- n_target_raw
    }
  } else {
    target_fit <- target_idx[fit_rows]
    ext_w_fit <- if (!is.null(ext_w)) ext_w[fit_rows] else NULL
    if (!is.null(ext_w_fit)) {
      w_target_vec <- ext_w_fit * as.numeric(target_fit)
      sum_w_target <- sum(w_target_vec[target_fit])
    } else {
      w_target_vec <- as.numeric(target_fit)
      sum_w_target <- sum(target_fit)
    }
  }

  # Observed-treatment predictions and residuals on study (fit) rows
  preds_obs <- stats::predict(
    outcome_model,
    newdata = fit_data,
    type = "response"
  )
  y_obs <- fit_data[[fit$outcome]]
  resid_obs <- y_obs - preds_obs

  # Outcome model prep (shared across interventions)
  outcome_prep <- prepare_model_if(outcome_model, aipw_fit_idx, n_total)

  # Univariate: single propensity model prep shared across interventions.
  # Multivariate: per-component prep deferred to inside the loop.
  prop_prep <- NULL
  if (!is_mv) {
    fit_idx_ps <- which(tm$fit_rows)
    if (inherits(propensity_model, "multinom")) {
      prop_prep <- prepare_model_if_multinom(
        propensity_model,
        fit_idx_ps,
        n_total
      )
    } else {
      prop_prep <- prepare_model_if(propensity_model, fit_idx_ps, n_total)
    }
  }

  # Outcome model family objects for Jacobian computation
  fam <- outcome_model$family
  beta_hat <- stats::coef(outcome_model)
  X_fit <- outcome_prep$X_fit

  # Transport: sampling weights and model prep
  w_S_fit <- NULL
  samp_prep <- NULL
  samp_wfn <- NULL
  gamma_hat <- NULL
  if (is_transport) {
    w_S <- compute_sampling_weights(
      fit$details$sampling_model,
      data,
      fit$target,
      fit$target_subset
    )
    w_S_fit <- w_S[fit_rows]

    samp_model <- fit$details$sampling_model
    samp_prep <- prepare_model_if(
      samp_model$model,
      which(samp_model$fit_rows),
      n_total
    )
    samp_wfn <- make_sampling_weight_fn(
      samp_model,
      data,
      fit$target,
      fit$target_subset
    )
    gamma_hat <- samp_model$gamma_hat
  }

  IF_list <- lapply(int_names, function(nm) {
    b <- aipw_bundles[[nm]]
    preds_g <- b$preds_g
    w_iv <- b$w_iv
    mu_hat_j <- mu_hat[[nm]]

    # Effective augmentation weight on study rows
    w_aug <- if (is_transport) w_S_fit * w_iv else w_iv

    # --- Channel 1: direct contribution ------------------------------------
    if (is_transport) {
      # Transport AIPW: IF contributions from two disjoint populations.
      # mu = (1/n_T) sum_{target} m(d,L) + (1/n_T) sum_{study} w_S * W_A * resid
      # Ch1_i = (n/n_T) * xi_i - mu_hat, where:
      #   xi_i = m(d,L_i) for target rows
      #   xi_i = w_S * W_A * resid for study rows
      #   xi_i = 0 otherwise
      Ch1_i <- rep(-mu_hat_j, n_total)

      # Target rows: prediction contribution
      target_which <- which(target_idx)
      preds_target <- b$preds_target
      if (!is.null(ext_w)) {
        Ch1_i[target_which] <- Ch1_i[target_which] +
          (n_total / sum_w_target) * ext_w[target_idx] * preds_target
      } else {
        Ch1_i[target_which] <- Ch1_i[target_which] +
          (n_total / sum_w_target) * preds_target
      }

      # Study rows: augmentation contribution
      study_which <- which(fit_rows)
      aug_study <- w_S_fit * w_iv * resid_obs
      if (!is.null(ext_w)) {
        Ch1_i[study_which] <- Ch1_i[study_which] +
          (n_total / sum_w_target) * ext_w[fit_rows] * aug_study
      } else {
        Ch1_i[study_which] <- Ch1_i[study_which] +
          (n_total / sum_w_target) * aug_study
      }

      # Rows in neither target nor study contribute only -mu_hat
      # (they are not in any estimating equation, so their net
      # contribution is zero after the corrections cancel).
      neither <- !target_idx & !fit_rows
      Ch1_i[neither] <- 0
    } else {
      # Non-transport: standard Hajek-scaled AIPW.
      # phi_i = Q_i(g) + W_i * (Y_i - Q_i(obs))
      # Ch1_i = n * (w_i / sum_w) * (phi_i - mu_hat)
      aipw_contrib <- preds_g + w_aug * resid_obs
      Ch1_i <- n_sub * (w_target_vec / sum_w_target) * (aipw_contrib - mu_hat_j)
      Ch1_i[!target_fit] <- 0
    }

    # --- Channel 2a: outcome model correction ------------------------------
    # J_beta = d mu_hat / d beta.
    data_a <- b$data_a
    X_star <- iv_design_matrix(outcome_model, data_a)
    eta_star <- as.numeric(X_star %*% beta_hat)
    mu_eta_star <- fam$mu.eta(eta_star)

    eta_obs <- as.numeric(X_fit %*% beta_hat)
    mu_eta_obs <- fam$mu.eta(eta_obs)

    if (is_transport) {
      # Term (a): d/dbeta of (1/n_T) sum_{target} m(d,L). The outcome
      # model was fit on study rows but we predict on target rows.
      # Use target predictions' design matrix.
      target_data <- data[target_idx]
      if (
        inherits(b$intervention, "causatr_intervention") &&
          b$intervention$type == "ipsi"
      ) {
        data_a_target <- target_data
      } else {
        data_a_target <- apply_intervention(
          target_data,
          fit$treatment,
          b$intervention
        )
      }
      X_star_target <- iv_design_matrix(outcome_model, data_a_target)
      eta_star_target <- as.numeric(X_star_target %*% beta_hat)
      mu_eta_star_target <- fam$mu.eta(eta_star_target)

      if (!is.null(ext_w)) {
        grad_a <- as.numeric(
          crossprod(X_star_target, ext_w[target_idx] * mu_eta_star_target)
        ) /
          sum_w_target
      } else {
        grad_a <- as.numeric(
          crossprod(X_star_target, mu_eta_star_target)
        ) /
          sum_w_target
      }

      # Term (b): d/dbeta of -(1/n_T) sum_{study} w_S * W_A * m(A,L)
      if (!is.null(ext_w)) {
        grad_b <- -as.numeric(
          crossprod(X_fit, ext_w[fit_rows] * w_aug * mu_eta_obs)
        ) /
          sum_w_target
      } else {
        grad_b <- -as.numeric(
          crossprod(X_fit, w_aug * mu_eta_obs)
        ) /
          sum_w_target
      }
    } else {
      # Non-transport: grad over target within fit_rows
      grad_a <- as.numeric(
        crossprod(X_star, w_target_vec * mu_eta_star)
      ) /
        sum_w_target

      grad_b <- -as.numeric(
        crossprod(X_fit, w_target_vec * w_aug * mu_eta_obs)
      ) /
        sum_w_target
    }

    J_beta <- grad_a + grad_b
    outcome_res <- apply_model_correction(outcome_prep, J_beta)

    # --- Channel 2b: propensity model correction ---------------------------
    # Natural-course intervention (NULL) has constant w_iv = 1, so
    # d(aug_mean)/d(alpha) = 0 and propensity correction vanishes.
    is_natural_course <- is.null(b$intervention)

    if (is_mv && !is_natural_course) {
      # Stacked-alpha closure across K independent propensity models.
      # The block-diagonal bread means each component's IF correction
      # depends only on its own alpha block.
      tms_local <- tms
      for (kk in seq_along(tms_local)) {
        tms_local[[kk]]$fit_rows <- rep(TRUE, nrow(fit_data))
      }
      class(tms_local) <- c("causatr_treatment_models", "list")

      mv_closure <- make_weight_fn_mv(
        treatment_models = tms_local,
        data = fit_data,
        interventions = b$intervention,
        estimand = "ATE"
      )
      mv_weight_fn <- mv_closure$weight_fn

      if (is_transport) {
        aug_mean_alpha <- function(alpha) {
          w_alpha <- mv_weight_fn(alpha)
          w_aug_alpha <- w_S_fit * w_alpha
          if (!is.null(ext_w)) {
            sum(ext_w[fit_rows] * w_aug_alpha * resid_obs) / sum_w_target
          } else {
            sum(w_aug_alpha * resid_obs) / sum_w_target
          }
        }
      } else {
        aug_mean_alpha <- function(alpha) {
          w_alpha <- mv_weight_fn(alpha)
          sum(w_target_vec * w_alpha * resid_obs) / sum_w_target
        }
      }

      J_alpha <- as.numeric(
        numDeriv::jacobian(aug_mean_alpha, x = mv_closure$alpha_hat)
      )

      # Per-component propensity corrections summed over K blocks.
      K <- length(tms_local)
      prop_correction <- rep(0, n_total)
      for (kk in seq_len(K)) {
        idx <- mv_closure$offsets[kk]:(mv_closure$offsets[kk + 1L] - 1L)
        J_alpha_k <- J_alpha[idx]
        prop_model_k <- tms_local[[kk]]$model
        prop_prep_k <- if (inherits(prop_model_k, "multinom")) {
          prepare_model_if_multinom(prop_model_k, aipw_fit_idx, n_total)
        } else {
          prepare_model_if(prop_model_k, aipw_fit_idx, n_total)
        }
        prop_res_k <- apply_model_correction(prop_prep_k, J_alpha_k)
        prop_correction <- prop_correction + prop_res_k$correction
      }
    } else if (is_natural_course) {
      prop_correction <- rep(0, n_total)
    } else {
      weight_fn <- make_weight_fn(
        tm,
        fit_data,
        b$intervention,
        estimand = estimand
      )

      if (is_transport) {
        aug_mean_alpha <- function(alpha) {
          w_alpha <- weight_fn(alpha)
          w_aug_alpha <- w_S_fit * w_alpha
          if (!is.null(ext_w)) {
            sum(ext_w[fit_rows] * w_aug_alpha * resid_obs) / sum_w_target
          } else {
            sum(w_aug_alpha * resid_obs) / sum_w_target
          }
        }
      } else {
        aug_mean_alpha <- function(alpha) {
          w_alpha <- weight_fn(alpha)
          sum(w_target_vec * w_alpha * resid_obs) / sum_w_target
        }
      }

      alpha_hat_raw <- stats::coef(propensity_model)
      if (!is.null(dim(alpha_hat_raw))) {
        alpha_hat <- as.vector(t(alpha_hat_raw))
      } else {
        alpha_hat <- alpha_hat_raw
      }

      J_alpha <- as.numeric(
        numDeriv::jacobian(aug_mean_alpha, x = alpha_hat)
      )
      prop_res <- apply_model_correction(prop_prep, J_alpha)
      prop_correction <- prop_res$correction
    }

    # --- Channel 2c: sampling model correction (transport only) ------------
    samp_correction <- rep(0, n_total)
    if (is_transport) {
      aug_mean_gamma <- function(gamma) {
        w_S_gamma <- samp_wfn(gamma)[fit_rows]
        w_aug_gamma <- w_S_gamma * w_iv
        if (!is.null(ext_w)) {
          sum(ext_w[fit_rows] * w_aug_gamma * resid_obs) / sum_w_target
        } else {
          sum(w_aug_gamma * resid_obs) / sum_w_target
        }
      }

      J_gamma <- as.numeric(numDeriv::jacobian(aug_mean_gamma, x = gamma_hat))
      samp_res <- apply_model_correction(samp_prep, J_gamma)
      samp_correction <- samp_res$correction
    }

    # --- Assembly ----------------------------------------------------------
    Ch1_i + outcome_res$correction - prop_correction - samp_correction
  })
  names(IF_list) <- int_names

  vcov_from_if(IF_list, n_if, int_names, cluster = cluster_vec)
}


#' Longitudinal AIPW sandwich variance (ICE-AIPW)
#'
#' @description
#' Analytic sandwich for the ICE-AIPW estimator (Bang & Robins 2005).
#' Derives from composition of Zivich et al. (2024, *Stat in Med*,
#' arXiv:2306.10976) for the ICE outcome model chain and Zivich et al.
#' (2024, *Biometrics* 81(2), arXiv:2404.16166) for the point-treatment
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
  cluster_vec = NULL
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
      rows_first
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
  rows_first
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
    beta_k <- stats::coef(model_k)

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

  K <- length(treatment_models_by_time)
  sub_fns <- vector("list", K)
  block_lens <- integer(K)
  alpha_blocks <- vector("list", K)
  align_idx_list <- vector("list", K)

  for (k in seq_len(K)) {
    tm_k <- treatment_models_by_time[[k]]
    data_k <- fit_data_by_time[[k]]
    ids_k <- as.character(data_k[[id_col]])

    sub_fns[[k]] <- make_weight_fn(
      treatment_model = tm_k,
      data = data_k,
      intervention = intervention,
      estimand = "ATE"
    )

    period_ids <- ids_k[tm_k$fit_rows]
    align_idx_list[[k]] <- match(period_ids, ids_first)
    alpha_k <- tm_k$alpha_hat
    alpha_blocks[[k]] <- alpha_k
    block_lens[k] <- length(alpha_k)
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
      # Propensity correction subtracted: block-lower-triangular M-estimation sign.
      IF_vec <- IF_vec - correction_id
    }
  }

  IF_vec
}
