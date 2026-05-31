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
#'   \mathrm{propensity\_correction}_i -
#'   \mathrm{sampling\_correction}_i -
#'   \mathrm{censoring\_correction}_i}
#'
#' where Ch1 is the direct AIPW functional evaluated at individual i
#' minus the estimate, and the corrections propagate nuisance-model
#' parameter uncertainty through the plug-in functional. The censoring
#' correction (Channel 2d) accounts for IPCW weight uncertainty via
#' three cross-derivative sub-components: direct Hajek effect,
#' outcome-model cross-block, and propensity-model cross-block.
#'
#' @param fit A `causatr_fit` object (AIPW estimator).
#' @param aipw_bundles Named list of per-intervention AIPW bundles from
#'   `compute_aipw_contrast_point()`.
#' @param target_idx Logical vector identifying target-population rows.
#' @param mu_hat Named numeric vector of point estimates per intervention.
#' @param aipw_fit_idx Logical vector of rows used for model fitting.
#' @param aipw_n_total Integer. Total number of target observations.
#' @param cluster_vec Integer or factor vector for cluster-robust SE, or `NULL`.
#'
#' @return A variance-covariance matrix (interventions x interventions).
#'
#' @noRd
variance_if_aipw <- function(
  fit,
  aipw_bundles,
  target_idx,
  mu_hat,
  aipw_fit_idx,
  aipw_n_total,
  cluster_vec = NULL,
  trim = 1
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
    if (inherits(propensity_model, "multinom")) {
      prop_prep <- prepare_model_if_multinom(
        propensity_model,
        aipw_fit_idx,
        n_total
      )
    } else {
      prop_prep <- prepare_model_if(propensity_model, aipw_fit_idx, n_total)
    }
  }

  # Outcome model family objects for Jacobian computation
  fam <- outcome_model$family
  beta_hat <- coef_clean(outcome_model)
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

  # IPCW: censoring model prep (shared across interventions)
  is_ipcw <- isTRUE(fit$details$ipcw)
  cens_prep <- NULL
  cens_wfn <- NULL
  cens_gamma_hat <- NULL
  w_pre_ipcw_fit <- NULL
  ipcw_w_fit <- NULL
  if (is_ipcw) {
    cens_model <- fit$details$censoring_model
    cens_prep <- prepare_model_if(
      cens_model$model,
      which(cens_model$fit_rows),
      n_total
    )
    cens_wfn <- make_ipcw_weight_fn(
      cens_model,
      n_total = n_total,
      censoring_col = as.integer(data[[fit$censoring]]),
      stabilize = TRUE
    )
    cens_gamma_hat <- cens_model$alpha_hat
    ipcw_w_fit <- fit$details$ipcw_weights[fit_rows]

    w_pre <- fit$details$weights_pre_ipcw
    w_pre_ipcw_fit <- if (!is.null(w_pre)) w_pre[fit_rows] else rep(1, n_sub)
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
      Ch1_fit <- n_sub *
        (w_target_vec / sum_w_target) *
        (aipw_contrib - mu_hat_j)
      Ch1_fit[!target_fit] <- 0
      # Pad to n_total so Ch1_i aligns with correction vectors from
      # apply_model_correction(), which are always length n_total.
      Ch1_i <- rep(0, n_total)
      Ch1_i[fit_rows] <- Ch1_fit
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
        estimand = "ATE",
        trim = trim
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
        estimand = estimand,
        trim = trim
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
        alpha_hat <- alpha_hat_raw[!is.na(alpha_hat_raw)]
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

    # --- Channel 2d: censoring model correction (IPCW only) ---------------
    # Three sub-components from the A_{*,gamma} cross-blocks:
    # (i)   Direct: d(AIPW_mean)/d(gamma) — Hajek weights depend on IPCW.
    # (ii)  Outcome: A_{beta,gamma}^T h_outcome — IPCW fitting weights
    #       enter the outcome model score.
    # (iii) Propensity: A_{alpha,gamma}^T h_prop — IPCW fitting weights
    #       enter the propensity model score.
    cens_correction <- rep(0, n_total)
    if (is_ipcw) {
      # (i) Direct effect: vary gamma in the AIPW augmented mean
      # (ext_w = w_pre_ipcw * w_ipcw, so varying gamma changes ext_w).
      # Transport and non-transport functionals differ: transport sums
      # over study rows with fixed denominator n_T; non-transport uses
      # the Hajek ratio where both numerator and denominator vary.
      if (is_transport) {
        aug_mean_cens <- function(gamma_c) {
          w_ipcw_c <- cens_wfn(gamma_c)[fit_rows]
          w_c <- w_pre_ipcw_fit * w_ipcw_c
          sum(w_c * w_S_fit * w_iv * resid_obs) / sum_w_target
        }
      } else {
        aipw_phi <- preds_g + w_iv * resid_obs
        aug_mean_cens <- function(gamma_c) {
          w_ipcw_c <- cens_wfn(gamma_c)[fit_rows]
          w_total_c <- w_pre_ipcw_fit * w_ipcw_c
          w_tgt_c <- w_total_c * as.numeric(target_fit)
          sw_c <- sum(w_tgt_c[target_fit])
          if (sw_c <= 0) {
            return(0)
          }
          sum(w_tgt_c * aipw_phi) / sw_c
        }
      }
      J_gamma_direct <- as.numeric(
        numDeriv::jacobian(aug_mean_cens, x = cens_gamma_hat)
      )

      # (ii) Outcome cross-block: the outcome model was fit with IPCW
      # weights, so its score psi_beta depends on gamma through the
      # fitting weights.
      # h = A_{bb}^{-1} J in M-estimation scaling: A_{bb} = (1/n)X'WX,
      # so A_{bb}^{-1} = n(X'WX)^{-1} = n * res$h.
      h_outcome <- n_sub * outcome_res$h
      out_score_factor <- mu_eta_obs * resid_obs / fam$variance(preds_obs)
      phi_bar_out_cens <- function(gamma_c) {
        w_ipcw_c <- cens_wfn(gamma_c)[fit_rows]
        w_c <- w_pre_ipcw_fit * w_ipcw_c
        as.numeric(crossprod(X_fit, w_c * out_score_factor)) / n_sub
      }
      A_beta_gamma <- -numDeriv::jacobian(phi_bar_out_cens, x = cens_gamma_hat)
      g_out_cens <- as.numeric(crossprod(A_beta_gamma, h_outcome))

      # (iii) Propensity cross-block: the propensity model was also fit
      # with IPCW weights, so its score depends on gamma too.
      g_prop_cens <- rep(0, length(cens_gamma_hat))
      if (!is_natural_course) {
        if (!is_mv) {
          X_ps <- stats::model.matrix(propensity_model)
          fam_ps <- propensity_model$family
          eta_ps <- propensity_model$linear.predictors
          mu_ps <- fam_ps$linkinv(eta_ps)
          ps_score_factor <- fam_ps$mu.eta(eta_ps) *
            (stats::model.response(stats::model.frame(propensity_model)) -
              mu_ps) /
            fam_ps$variance(mu_ps)
          h_prop <- n_sub * prop_res$h

          phi_bar_ps_cens <- function(gamma_c) {
            w_ipcw_c <- cens_wfn(gamma_c)[fit_rows]
            w_c <- w_pre_ipcw_fit * w_ipcw_c
            as.numeric(crossprod(X_ps, w_c * ps_score_factor)) / n_sub
          }
          A_alpha_gamma <- -numDeriv::jacobian(
            phi_bar_ps_cens,
            x = cens_gamma_hat
          )
          g_prop_cens <- as.numeric(crossprod(A_alpha_gamma, h_prop))
        }
        if (is_mv) {
          for (kk in seq_along(tms)) {
            prop_model_kk <- tms[[kk]]$model
            if (inherits(prop_model_kk, "multinom")) {
              next
            }
            X_ps_k <- stats::model.matrix(prop_model_kk)
            fam_ps_k <- prop_model_kk$family
            eta_ps_k <- prop_model_kk$linear.predictors
            mu_ps_k <- fam_ps_k$linkinv(eta_ps_k)
            ps_sf_k <- fam_ps_k$mu.eta(eta_ps_k) *
              (stats::model.response(stats::model.frame(prop_model_kk)) -
                mu_ps_k) /
              fam_ps_k$variance(mu_ps_k)

            idx_k <- mv_closure$offsets[kk]:(mv_closure$offsets[kk + 1L] - 1L)
            prop_prep_kk <- prepare_model_if(
              prop_model_kk,
              aipw_fit_idx,
              n_total
            )
            prop_res_kk <- apply_model_correction(prop_prep_kk, J_alpha[idx_k])
            h_prop_k <- n_sub * prop_res_kk$h

            phi_bar_ps_k <- function(gamma_c) {
              w_ipcw_c <- cens_wfn(gamma_c)[fit_rows]
              w_c <- w_pre_ipcw_fit * w_ipcw_c
              as.numeric(crossprod(X_ps_k, w_c * ps_sf_k)) / n_sub
            }
            A_ak_gamma <- -numDeriv::jacobian(
              phi_bar_ps_k,
              x = cens_gamma_hat
            )
            g_prop_cens <- g_prop_cens +
              as.numeric(crossprod(A_ak_gamma, h_prop_k))
          }
        }
      }

      g_cens_total <- J_gamma_direct + g_out_cens + g_prop_cens
      cens_res <- apply_model_correction(cens_prep, g_cens_total)
      cens_correction <- cens_res$correction
    }

    # --- Assembly ----------------------------------------------------------
    # All corrections are additive: each propagates nuisance-model
    # parameter uncertainty through J_k A_k^{-1} psi_{k,i}.  The
    # Jacobian J_k = d mu / d theta_k already carries the correct sign,
    # so the correction is added, not subtracted. This follows the
    # standard M-estimation block-inverse identity (Stefanski & Boos
    # 2002, eq. 9): B^{-1}[target, nuisance] = -A_{22}^{-1} A_{21}
    # A_{11}^{-1}, where A_{21} = -d(psi_target)/d(theta_nuisance)
    # incorporates a sign flip from the bread definition.
    Ch1_i +
      outcome_res$correction +
      prop_correction +
      samp_correction +
      cens_correction
  })
  names(IF_list) <- int_names

  vcov_from_if(IF_list, n_if, int_names, cluster = cluster_vec)
}
