#' Numerical fallback for the IF variance engine
#'
#' @description
#' Two-tier fallback used by `variance_if()` when the model exposes neither
#' `family$mu.eta` / `family$variance` (no analytic GLM bread) nor `$Vp`
#' (no GAM bread). Both tiers compute Channel 1 directly from the
#' predictions and use `numDeriv::jacobian()` on `predict()` to get the
#' marginal-mean Jacobian \eqn{J} that would otherwise come from
#' `iv_design_matrix()`.
#'
#' ## Tier 1 \(\(`sandwich::estfun` available\)\)
#'
#' If `sandwich::estfun(fit$model)` returns the per-observation score
#' matrix \eqn{\psi_i}, combine it with `sandwich::bread(fit$model)/nobs`
#' to get \eqn{\mathrm{IF}_{\beta,i} = A^{-1}\psi_i} per observation, then
#' add the Channel 1 vector and aggregate via `vcov_from_if()`. **This
#' recovers the full IF (Channel 1 + Channel 2 + cross-term) exactly**,
#' matching the analytic main path whenever `sandwich::estfun` has a
#' method for the model class.
#'
#' ## Tier 2 \(plain delta shortcut, with a warning\)
#'
#' If no `estfun` method exists, fall back to
#' \eqn{V_1 + V_2 = (1/n^2)\mathrm{Ch.1}^T\mathrm{Ch.1} + J V_\beta J^T}
#' and return their sum, where \eqn{\mathrm{Ch.1}} is the \eqn{n \times k}
#' matrix stacking Channel-1 vectors column-wise so that \eqn{V_1} carries
#' full off-diagonal covariance, not just diagonal variances. This
#' **drops the cross-term**
#' \eqn{(1/n^2)(\mathrm{Ch.1}^T\mathrm{Ch.2} + \mathrm{Ch.2}^T\mathrm{Ch.1})}
#' from the decomposition in vignette Section 3.3, so the variance is
#' slightly miscalibrated. A `rlang::warn()` is emitted so the user knows.
#'
#' @param fit A `causatr_fit` object.
#' @param data_a_list Named list of counterfactual data.tables.
#' @param preds_list Named list of length-`n_total` prediction vectors.
#' @param mu_hat Named numeric vector of marginal-mean point estimates.
#' @param target_idx Logical vector (length `n_total`) flagging target rows.
#' @param weights Numeric vector of target-row weights, or `NULL`.
#'
#' @return A named `k x k` variance-covariance matrix.
#'
#' @noRd
variance_if_numeric <- function(
  fit,
  data_a_list,
  preds_list,
  mu_hat,
  target_idx,
  weights = NULL,
  cluster_vec = NULL
) {
  model <- fit$model
  int_names <- names(data_a_list)
  k <- length(int_names)
  n_total <- length(target_idx)
  beta_hat <- stats::coef(model)

  # Row-subset each counterfactual frame to the target population.
  # `da[target_idx, , drop = FALSE]` is polymorphic: it does row-selection
  # on both data.frame and data.table (the trailing empty `j` forces
  # data.frame's `[` into row-index mode instead of column-select).
  # Avoids the previous `as.data.frame(da)[...]` roundtrip.
  data_a_frames <- lapply(
    data_a_list,
    function(da) da[target_idx, , drop = FALSE]
  )

  # Normalise weights to sum to 1 on target rows; non-target rows stay 0.
  target_w <- numeric(n_total)
  if (!is.null(weights)) {
    sum_w <- sum(weights)
    target_w[target_idx] <- weights / sum_w
  } else {
    n_t <- sum(target_idx)
    target_w[target_idx] <- 1 / n_t
  }

  # Ch1_i = n * w_i * (pred_i - mu_hat). The n_total prefactor converts
  # the normalised weight into the IF scale: Var = (1/n^2) sum Ch1_i^2.
  Ch1_list <- lapply(seq_len(k), function(j) {
    p <- preds_list[[j]]
    n_total * target_w * (p - mu_hat[j])
  })

  # Marginal-mean Jacobian J (k x p), via numDeriv on predict().
  is_betareg <- inherits(model, "betareg")
  pred_fun <- function(beta) {
    model_tmp <- model
    if (is_betareg) {
      # betareg stores coefficients as list($mean, $precision).
      n_mean <- length(model$coefficients$mean)
      model_tmp$coefficients$mean <- beta[seq_len(n_mean)]
      model_tmp$coefficients$precision <- beta[(n_mean + 1L):length(beta)]
    } else {
      model_tmp$coefficients <- beta
    }
    vapply(
      data_a_frames,
      function(df) {
        preds <- stats::predict(model_tmp, newdata = df, type = "response")
        if (!is.null(weights)) {
          stats::weighted.mean(preds, weights, na.rm = TRUE)
        } else {
          mean(preds, na.rm = TRUE)
        }
      },
      numeric(1)
    )
  }
  # J[j, l] = d mu_hat_j / d beta_l, the (k x p) Jacobian of the k
  # counterfactual means w.r.t. the outcome model parameters.
  # numDeriv::jacobian() returns a (k x p) matrix when pred_fun maps
  # R^p -> R^k; coerce in case k=1 collapses it to a vector.
  J <- numDeriv::jacobian(pred_fun, x = beta_hat)
  if (!is.matrix(J)) {
    J <- matrix(J, nrow = k)
  }

  # Tier 1: IF_beta = psi_i A^{-1} (row i of score times bread inverse).
  # `sandwich::bread(model)` returns n * (X'WX)^{-1}, so dividing by nobs
  # gives (X'WX)^{-1} = A^{-1} on the per-observation scale. The product
  # psi %*% A_inv then yields a (n_fit x k) matrix of per-obs IF
  # contributions to the k counterfactual means, recovering Channel 2
  # exactly for any model class that has an estfun() method.
  estfun_ok <- tryCatch(
    {
      ef <- sandwich::estfun(model)
      is.matrix(ef) && nrow(ef) > 0L
    },
    error = function(e) FALSE
  )

  if (estfun_ok) {
    psi <- sandwich::estfun(model)
    # sandwich::bread() = n * (X'WX)^{-1}; dividing by nobs normalises to A^{-1}.
    A_inv <- sandwich::bread(model) / stats::nobs(model)
    # IF_beta[i, ] = psi_i A^{-1}: the per-obs nuisance IF, an (n_fit x p) matrix.
    IF_beta <- psi %*% A_inv
    fit_idx <- resolve_fit_idx(fit, model)
    if (length(fit_idx) != nrow(IF_beta)) {
      rlang::warn(
        paste0(
          "variance_if_numeric() Tier 1: estfun() row count (",
          nrow(IF_beta),
          ") does not match fit_idx (",
          length(fit_idx),
          "). Falling back to Tier 2."
        )
      )
    } else {
      # Ch2_fit[i, j] = n * IF_beta[i, ] %*% J[j, ]: the Channel-2 contribution
      # of obs i to intervention j. The n_total factor converts the score-scale
      # IF_beta (which uses the normalised A^{-1}) back to the IF scale where
      # Var = (1/n^2) sum IF_i^2. Transposing J gives a (n_fit x k) result
      # in one BLAS call rather than k separate matrix-vector products.
      Ch2_fit <- (IF_beta %*% t(J)) * n_total
      IF_list <- lapply(seq_len(k), function(j) {
        IF <- Ch1_list[[j]]
        IF[fit_idx] <- IF[fit_idx] + Ch2_fit[, j]
        IF
      })
      names(IF_list) <- int_names
      return(vcov_from_if(
        IF_list,
        n_total,
        int_names,
        cluster = cluster_vec
      ))
    }
  }

  # Tier 2: drop the cross-term, warn. We attach a classed warning
  # (`causatr_tier2_fallback`) so batch / CI pipelines can grep for
  # this exact condition via `withCallingHandlers(..., causatr_tier2_fallback = ...)`.
  # Plain `rlang::warn()` with a string is not discoverable once it
  # lands in a knit log, and this fallback is a real calibration risk
  # the user should be able to react to programmatically.
  rlang::warn(
    paste0(
      "Model class '",
      class(model)[1],
      "' exposes neither analytic bread (family$mu.eta / $Vp) nor a ",
      "sandwich::estfun() method. Falling back to V1 + J V_beta J^T, ",
      "which drops the IF cross-term (vignette Section 3.3) and may ",
      "slightly miscalibrate the variance. Use ci_method = 'bootstrap' ",
      "for an exact answer."
    ),
    class = "causatr_tier2_fallback"
  )

  Ch1_mat <- do.call(cbind, Ch1_list)
  # Under Tier 2 we only own the V1 piece fully; V2 = J V_beta J^T
  # ships whatever vcov the model exposes and we cannot cluster-adjust
  # it without plumbing through the model's score. V1 is clusterable
  # directly -- sum Ch1 within cluster before squaring -- so we apply
  # the adjustment there and warn that V2 remains model-level.
  if (!is.null(cluster_vec)) {
    cluster_f <- factor(cluster_vec, levels = unique(cluster_vec))
    Ch1_mat <- rowsum(Ch1_mat, cluster_f, reorder = FALSE)
    rlang::warn(
      paste0(
        "Tier 2 numerical fallback: cluster-robust aggregation applied ",
        "only to the Channel-1 variance block. The J V_beta J^T piece ",
        "uses the model's standard vcov. Use `ci_method = \"bootstrap\"` ",
        "for a fully cluster-robust variance."
      ),
      class = "causatr_tier2_cluster_partial"
    )
  }
  V1 <- crossprod(Ch1_mat) / n_total^2
  V_beta <- tryCatch(stats::vcov(model), error = function(e) NULL)
  if (is.null(V_beta)) {
    V_beta <- diag(0, length(beta_hat))
  }
  V2 <- J %*% V_beta %*% t(J)

  vcov_mat <- V2 + V1
  rownames(vcov_mat) <- int_names
  colnames(vcov_mat) <- int_names
  # Tag the returned vcov so downstream code can detect Tier-2 fallback
  # post-hoc without parsing warning output. The classed
  # `causatr_tier2_fallback` warning is only visible during the variance
  # calculation; the attribute survives as long as the vcov is stored,
  # which lets users test `attr(result$vcov, "tier2_approximate")`.
  attr(vcov_mat, "tier2_approximate") <- TRUE
  vcov_mat
}


#' Influence-function variance dispatcher
#'
#' @description
#' Single entry point for sandwich/IF variance estimation across all four
#' methods. Dispatches to a method-specific branch and returns the
#' \eqn{k \times k} marginal-mean variance-covariance matrix.
#'
#' Branches (vignette `variance-theory.qmd`, Sections 4 and 5):
#'
#' - **g-comp** (`variance_if_gcomp`): one model correction via
#'   `correct_model()`; standard (independent) IF aggregation.
#' - **matching** (`variance_if_matching`): one model correction; the
#'   matched sample defines the IF scope (`n = nrow(matched_data)`,
#'   `fit_idx = 1..n_match`); cluster-robust aggregation on
#'   `matched$subclass` per D1.
#' - **ICE** (`variance_if_ice`): \eqn{K{+}1} model corrections in a
#'   forward sensitivity recursion. Models are fit backward in time
#'   (\eqn{K \to 0}) inside `ice_iterate()`, but sensitivities propagate
#'   **forward** (\eqn{0 \to K}) through this loop because the
#'   block-triangular bread is inverted by back-substitution from the
#'   mean equation downward (vignette Section 5.4 / D2).
#' - **IPW** (`variance_if_ipw`): same Channel 1 shape as g-comp;
#'   Channel 2 is the density-ratio stacked sandwich -- the MSM
#'   correction plus the propensity correction tied together by a
#'   numerical cross-derivative \eqn{A_{\beta\alpha}} from
#'   `numDeriv::jacobian()`. Per-intervention assembly lives in
#'   `compute_ipw_if_self_contained_one()`.
#'
#' Non-GLM models without `family$mu.eta` and without `$Vp` are routed to
#' `variance_if_numeric()`.
#'
#' @param fit A `causatr_fit` object.
#' @param interventions Named list of `causatr_intervention` objects (or
#'   `NULL` natural-course entries). Required for the matching branch
#'   because matching rebuilds counterfactual data on the matched sample;
#'   may be `NULL` for g-comp / ICE branches that already have
#'   `data_a_list` / `ice_results`.
#' @param data_a_list Named list of counterfactual data.tables (g-comp).
#' @param preds_list Named list of length-`n_total` prediction vectors
#'   (g-comp).
#' @param mu_hat Named numeric vector of marginal-mean point estimates.
#' @param target_idx Logical vector (length `n_total`) flagging target rows
#'   (g-comp).
#' @param ice_results Named list of `ice_iterate()` results (ICE only).
#' @param ice_aipw_results Named list of `ice_aipw_iterate()` results
#'   (longitudinal AIPW only).
#' @param target_within_first Logical vector over first-time individuals
#'   (ICE only).
#'
#' @return A named `k x k` variance-covariance matrix.
#'
#' @noRd
variance_if <- function(
  fit,
  interventions = NULL,
  data_a_list = NULL,
  preds_list = NULL,
  mu_hat = NULL,
  target_idx = NULL,
  ice_results = NULL,
  ice_aipw_results = NULL,
  target_within_first = NULL,
  ipw_bundles = NULL,
  ipw_fit_idx = NULL,
  ipw_n_total = NULL,
  aipw_bundles = NULL,
  aipw_fit_idx = NULL,
  aipw_n_total = NULL,
  cluster_vec = NULL
) {
  if (fit$type == "longitudinal") {
    # Longitudinal IPW dispatch. The IPW path needs the
    # per-intervention bundles + the `ids_first` mapping (carried in
    # `ipw_bundles` / `ipw_long`) that `compute_ipw_contrast_longitudinal()`
    # built; ICE keeps its own `ice_results` route. Both go through
    # `vcov_from_if()` after the per-intervention IF assembly.
    if (fit$estimator == "ipw") {
      return(variance_if_ipw_longitudinal(
        fit,
        ipw_bundles,
        target_within_first,
        cluster_vec = cluster_vec
      ))
    }
    if (fit$estimator == "aipw") {
      return(variance_if_aipw_longitudinal(
        fit,
        ice_aipw_results,
        target_within_first,
        cluster_vec = cluster_vec
      ))
    }
    return(variance_if_ice(
      fit,
      ice_results,
      target_within_first,
      cluster_vec = cluster_vec
    ))
  }

  estimator <- fit$estimator

  if (estimator == "gcomp") {
    return(variance_if_gcomp(
      fit,
      data_a_list,
      preds_list,
      mu_hat,
      target_idx,
      interventions = interventions,
      cluster_vec = cluster_vec
    ))
  }

  if (estimator == "matching") {
    # Matching's cluster-robust aggregation on `subclass` is structural
    # and mutually exclusive with a user-supplied `cluster`. The
    # `resolve_cluster()` guard upstream already hard-aborts in that
    # case, so by the time we land here `cluster_vec` must be NULL.
    return(variance_if_matching(fit, interventions))
  }

  if (estimator == "ipw") {
    return(variance_if_ipw(
      fit,
      ipw_bundles,
      target_idx,
      mu_hat,
      ipw_fit_idx,
      ipw_n_total,
      cluster_vec = cluster_vec
    ))
  }

  if (estimator == "aipw") {
    return(variance_if_aipw(
      fit,
      aipw_bundles,
      target_idx,
      mu_hat,
      aipw_fit_idx,
      aipw_n_total,
      cluster_vec = cluster_vec
    ))
  }

  rlang::abort(paste0("Unknown estimator '", estimator, "' in variance_if()."))
}


#' Shared Channel-1 + Jacobian builder for point-treatment branches
#'
#' @description
#' Collects the per-intervention pieces that g-comp and IPW both need
#' before their Channel-2 dispatch:
#'
#' - `Ch1_list`: length-`n` Channel-1 IF vectors (one per intervention),
#'   zero off target rows; weighted or unweighted as determined by
#'   `fit$details$weights`.
#' - `grad_list`: length-`p` marginal-mean Jacobians
#'   \eqn{g_j = \partial\hat\mu_j/\partial\beta} (one per intervention),
#'   via \eqn{X^{*T}(d\mu/d\eta)} averaged over the target population.
#' - `fit_idx`, `n`: the row-alignment scope the Channel-2 call will use.
#'
#' Used by `variance_if_gcomp()` to build Channel 1 and the
#' per-intervention marginal-mean Jacobian off a single fitted
#' outcome model, so the Channel-2 step is just a loop over
#' `apply_model_correction(prep, grad_list[[j]])`. `variance_if_ipw()`
#' does not share this plumbing because its MSM is refit per
#' intervention on a density-ratio-weighted dataset -- it builds its
#' own Channel 1 / Jacobian pieces row-locally on the MSM fit-row
#' subset.
#'
#' @param fit A `causatr_fit` object.
#' @param data_a_list Named list of counterfactual data.tables.
#' @param preds_list Named list of length-`n` prediction vectors.
#' @param mu_hat Named numeric vector of marginal-mean point estimates.
#' @param target_idx Logical vector (length `n`) flagging target rows.
#'
#' @return A list with components `Ch1_list`, `grad_list`, `fit_idx`, `n`.
#'
#' @noRd
build_point_channel_pieces <- function(
  fit,
  data_a_list,
  preds_list,
  mu_hat,
  target_idx,
  interventions = NULL
) {
  model <- fit$model
  int_names <- names(data_a_list)
  k <- length(int_names)
  n <- nrow(fit$data)

  ext_w <- fit$details$weights
  has_weights <- !is.null(ext_w)

  # Guard against an empty target population. `compute_contrast()`
  # normally filters this out upstream, but a defensive check here
  # converts a cryptic NaN-propagated vcov into a clear abort at the
  # variance-engine boundary.
  if (has_weights) {
    sum_w_target <- sum(ext_w[target_idx])
    if (sum_w_target <= 0) {
      # Classed abort so `compute_contrast()`'s `by`-skip path can match
      # on class, not on the English message text. See C3 in the
      # 2026-04-15 second-round critical review.
      rlang::abort(
        "build_point_channel_pieces(): target-population weights sum to 0.",
        class = "causatr_empty_target"
      )
    }
  } else {
    n_target <- sum(target_idx)
    if (n_target == 0L) {
      rlang::abort(
        "build_point_channel_pieces(): target population is empty.",
        class = "causatr_empty_target"
      )
    }
  }

  fit_idx <- resolve_fit_idx(fit, model)
  beta_hat <- stats::coef(model)

  # Polymorphic row-subset: works on both data.frame and data.table.
  data_a_frames <- lapply(
    data_a_list,
    function(da) da[target_idx, , drop = FALSE]
  )

  Ch1_list <- vector("list", k)
  grad_list <- vector("list", k)

  # `preds_list[[j]]` is length-`n` (full data). `compute_contrast()`
  # upstream has already masked `target_idx` so that all TRUE rows have
  # non-NA predictions, but non-target rows may still carry NA. We must
  # therefore index `p` by `target_idx` rather than multiplying the
  # full-length `p` by a zero vector -- `0 * NA = NA` in R, which would
  # silently corrupt the Channel-1 contribution on non-target rows.
  for (j in seq_len(k)) {
    p <- preds_list[[j]]
    ch1 <- numeric(n)
    if (has_weights) {
      # Weighted Ch1_i = n * (w_i / sum_w) * (pred_i - mu_hat_j).
      # The (w_i / sum_w) ratio normalises so the weighted average of
      # weights equals 1/n_target; multiplying by n puts the result on
      # the IF scale where Var = (1/n^2) sum Ch1_i^2.
      ch1[target_idx] <- n *
        (ext_w[target_idx] / sum_w_target) *
        (p[target_idx] - mu_hat[j])
    } else {
      # Unweighted: n / n_target accounts for the fact that the target
      # subpopulation (e.g. ATT treated rows) may be smaller than n.
      ch1[target_idx] <- (n / n_target) * (p[target_idx] - mu_hat[j])
    }
    Ch1_list[[j]] <- ch1

    iv_j <- if (!is.null(interventions)) interventions[[j]] else NULL
    if (has_stochastic_component(iv_j)) {
      # For stochastic interventions, E[mu'(eta^*)] must be estimated by
      # MC since the counterfactual distribution has no closed form. Each
      # draw samples a new A* ~ g(A | L), producing a different X^* and
      # hence a different mu'(eta^*). Averaging over n_mc draws approximates
      # the analytic gradient g_j = E_{g}[X^{*T} mu'(eta^*)] / n_target.
      n_mc <- get_stochastic_n_mc(iv_j)
      p_coef <- length(beta_hat)
      grad_sum <- numeric(p_coef)
      data_full <- fit$data
      treatment <- fit$treatment
      for (m in seq_len(n_mc)) {
        data_a_m <- apply_intervention(data_full, treatment, iv_j)
        da_m_target <- data_a_m[target_idx, , drop = FALSE]
        X_star_m <- iv_design_matrix(model, da_m_target)
        eta_star_m <- as.numeric(X_star_m %*% beta_hat)
        mu_eta_star_m <- model$family$mu.eta(eta_star_m)
        if (has_weights) {
          w_t <- ext_w[target_idx]
          grad_sum <- grad_sum +
            as.numeric(crossprod(X_star_m, w_t * mu_eta_star_m))
        } else {
          grad_sum <- grad_sum +
            as.numeric(crossprod(X_star_m, mu_eta_star_m))
        }
      }
      denom <- if (has_weights) sum_w_target * n_mc else n_target * n_mc
      grad_list[[j]] <- grad_sum / denom
    } else {
      # Analytic gradient: g_j = (X^{*T} diag(mu'(eta^*)) 1) / n_target.
      # This is d mu_hat_j / d beta via the chain rule:
      #   d/d beta E[mu(X^* beta)] = E[mu'(X^* beta) X^{*T}]
      # evaluated at beta_hat over the target population.
      X_star <- iv_design_matrix(model, data_a_frames[[j]])
      eta_star <- as.numeric(X_star %*% beta_hat)
      mu_eta_star <- model$family$mu.eta(eta_star)

      if (has_weights) {
        w_t <- ext_w[target_idx]
        grad_list[[j]] <- as.numeric(
          crossprod(X_star, w_t * mu_eta_star)
        ) /
          sum_w_target
      } else {
        grad_list[[j]] <- as.numeric(
          crossprod(X_star, mu_eta_star)
        ) /
          n_target
      }
    }
  }

  names(Ch1_list) <- int_names
  names(grad_list) <- int_names

  list(
    Ch1_list = Ch1_list,
    grad_list = grad_list,
    fit_idx = fit_idx,
    n = n
  )
}


#' Bundle `build_point_channel_pieces()` + `prepare_model_if()` together
#'
#' @description
#' G-comp's branch follows a fixed prepare sequence: build the Channel-1
#' / Jacobian pieces, then prepare the outcome model for repeated
#' `apply_model_correction()` calls. Passing `pieces$fit_idx` /
#' `pieces$n` to `prepare_model_if()` by hand leaves two ways to get out
#' of sync -- a caller could pass `resolve_fit_idx(fit, model)` and
#' `nrow(fit$data)` independently, which would silently drift if one of
#' them were computed differently. This helper bundles the two into one
#' call so the alignment is structural, not by convention.
#'
#' The IPW branch does not use this bundler because it routes its
#' Channel 1 / Jacobian construction through
#' `variance_if_ipw()` directly on the MSM fit-row subset -- the MSM is
#' refit per intervention on a density-ratio-weighted dataset, so the
#' shared point-treatment plumbing does not apply.
#'
#' @inheritParams build_point_channel_pieces
#' @param model The fitted outcome (or MSM) model to prepare.
#'
#' @return A list with components `pieces` (output of
#'   `build_point_channel_pieces()`) and `prep` (output of
#'   `prepare_model_if()`) sharing the same `fit_idx` / `n`.
#'
#' @noRd
prepare_point_variance <- function(
  fit,
  model,
  data_a_list,
  preds_list,
  mu_hat,
  target_idx,
  interventions = NULL
) {
  pieces <- build_point_channel_pieces(
    fit,
    data_a_list,
    preds_list,
    mu_hat,
    target_idx,
    interventions = interventions
  )
  prep <- prepare_model_if(model, pieces$fit_idx, pieces$n)
  list(pieces = pieces, prep = prep)
}


#' IPCW Channel 2 correction for the outcome/MSM model
#'
#' @description
#' Computes the censoring model's contribution to the stacked
#' M-estimation IF via the cross-derivative
#' \eqn{A_{\beta,\gamma} = -(1/n) \sum_i \partial\psi_\beta/\partial\gamma}.
#'
#' The closure \eqn{\bar\Phi(\gamma)} evaluates the outcome/MSM score
#' with IPCW weights recomputed at candidate \eqn{\gamma}:
#' \deqn{\bar\Phi_j(\gamma) = \frac{1}{n} \sum_i X_i \, w_{\mathrm{ipcw}}(\gamma)_i
#'       \, w_{\mathrm{ext},i} \, r_i \, \mu'(\eta_i) / V(\mu_i)}
#' where \eqn{r_i}, \eqn{\mu'}, \eqn{V} are held fixed at
#' \eqn{\hat\beta}. Only the IPCW weight factor varies with \eqn{\gamma}.
#'
#' The per-individual IF correction is then
#' \deqn{\mathrm{cens\_correction}_i = g_\gamma^T A_{\gamma\gamma}^{-1}
#'       \psi_{\gamma,i}}
#' where \eqn{g_\gamma = A_{\beta,\gamma}^T h} and
#' \eqn{h = A_{\beta\beta}^{-1} J} from the outcome model.
#'
#' @param fit A `causatr_fit` with `details$ipcw == TRUE`.
#' @param outcome_model The fitted outcome (or MSM) GLM.
#' @param outcome_prep Output of `prepare_model_if()` for the outcome model.
#' @param h_outcome Numeric p-vector. The bread-projected gradient
#'   \eqn{A_{\beta\beta}^{-1} J}. For gcomp, this is
#'   `apply_model_correction(prep, grad)$h`; for IPW, `n_fit * msm_res$h`.
#' @param n_fit Integer. Number of rows in the outcome model fit.
#'
#' @return A named list:
#'   \describe{
#'     \item{`correction`}{Numeric vector of length `n_total`.}
#'   }
#'
#' @noRd
compute_ipcw_if_correction <- function(
  fit,
  outcome_model,
  outcome_prep,
  h_outcome,
  n_fit
) {
  cens_model <- fit$details$censoring_model
  n_total <- outcome_prep$n_total

  # Build the phi_bar_cens(gamma) closure. The outcome model's GLM
  # score for observation i is
  #   psi_{beta,i} = X_i * w_total_i * (y_i - mu_i) * mu'(eta_i) / V(mu_i)
  # where w_total_i = w_ext_i * w_ipcw_i. Varying gamma only changes
  # w_ipcw_i; everything else is fixed at beta_hat.
  beta_hat <- stats::coef(outcome_model)
  X_fit <- outcome_prep$X_fit
  fam <- outcome_model$family
  eta_fit <- as.numeric(X_fit %*% beta_hat)
  mu_fit <- fam$linkinv(eta_fit)
  mu_eta_fit <- fam$mu.eta(eta_fit)
  var_mu_fit <- fam$variance(mu_fit)
  y_fit <- stats::model.response(stats::model.frame(outcome_model))
  r_fit <- y_fit - mu_fit

  # Pre-IPCW weights: the weight factor that doesn't depend on gamma.
  # This is the external weights (survey, etc.) that were stashed
  # before IPCW composition.
  w_ext <- fit$details$weights_pre_ipcw
  if (is.null(w_ext)) {
    w_ext_fit <- rep(1, n_fit)
  } else {
    fit_idx <- outcome_prep$fit_idx
    w_ext_fit <- w_ext[fit_idx]
  }

  # IPCW weight closure: maps gamma -> IPCW weight vector (full length)
  ipcw_wfn <- make_ipcw_weight_fn(
    cens_model,
    n_total = n_total,
    censoring_col = as.integer(fit$data[[fit$censoring]]),
    stabilize = TRUE
  )

  fit_idx <- outcome_prep$fit_idx

  # phi_bar(gamma) = (1/n) sum_i X_i * w_ext_i * w_ipcw_i(gamma) * r_i * mu'(eta_i) / V(mu_i)
  # This is the outcome/MSM score as a function of the censoring parameter gamma,
  # holding beta fixed at beta_hat. Its negative Jacobian is the cross-block
  # A_{beta,gamma} of the stacked M-estimation bread.
  phi_bar_cens <- function(gamma) {
    w_ipcw <- ipcw_wfn(gamma)
    w_ipcw_fit <- w_ipcw[fit_idx]
    s_per_i <- w_ext_fit * w_ipcw_fit * mu_eta_fit * r_fit / var_mu_fit
    as.numeric(crossprod(X_fit, s_per_i)) / n_fit
  }

  gamma_hat <- cens_model$alpha_hat
  # A_{beta,gamma} = -d phi_bar / d gamma evaluated at (beta_hat, gamma_hat).
  # The minus sign comes from the convention A = -E[d psi / d theta] for the
  # M-estimation bread (Stefanski & Boos 2002, eq. 8).
  A_beta_gamma <- -numDeriv::jacobian(phi_bar_cens, x = gamma_hat)

  # g_cens = A_{beta,gamma}^T h: projects the outcome-model sensitivity h
  # onto the censoring-parameter space, giving the gradient of mu_hat
  # w.r.t. gamma via the chain rule through the stacked system.
  g_cens <- as.numeric(crossprod(A_beta_gamma, h_outcome))

  cens_prep <- prepare_model_if(
    cens_model$model,
    which(cens_model$fit_rows),
    n_total
  )
  cens_res <- apply_model_correction(cens_prep, g_cens)

  list(correction = cens_res$correction)
}


#' IPCW Channel 2 correction for the IPW MSM
#'
#' @description
#' IPW-specific IPCW correction. Unlike `compute_ipcw_if_correction()`
#' (which works on the full data for gcomp), the IPW path operates on a
#' subset `data[fit_rows]` of length `n_sub`. The censoring model was
#' fit on the full data, so we need to:
#' 1. Build the phi_bar closure on the n_sub-sized MSM
#' 2. Compute the cross-derivative A_{beta,gamma}
#' 3. Apply the censoring model correction on n_sub rows (subsetting
#'    the censoring model's score to the fit_rows subset)
#'
#' @param fit A `causatr_fit` with `details$ipcw == TRUE`.
#' @param msm_model The per-intervention MSM model.
#' @param J Numeric p-vector. Marginal-mean Jacobian.
#' @param fit_rows Logical vector (length n_total). MSM fit rows.
#' @param n_sub Integer. Number of MSM fit rows.
#'
#' @return Numeric vector of length `n_sub`.
#'
#' @noRd
compute_ipw_ipcw_correction <- function(
  fit,
  msm_model,
  J,
  fit_rows,
  n_sub
) {
  cens_model <- fit$details$censoring_model
  n_total <- nrow(fit$data)

  # fit_idx for the MSM is 1..n_sub because the MSM was fit on the
  # already-subsetted n_sub rows, not on the full n_total dataset.
  msm_prep <- prepare_model_if(msm_model, seq_len(n_sub), n_sub)
  msm_res <- apply_model_correction(msm_prep, J)
  # h_msm = A_{bb}^{-1} J in M-estimation scaling (same n-rescaling as
  # the gcomp IPCW path in compute_ipcw_if_correction).
  h_msm <- n_sub * msm_res$h

  # MSM score ingredients (fixed at beta_hat)
  beta_hat <- stats::coef(msm_model)
  X_msm <- msm_prep$X_fit
  fam <- msm_model$family
  eta_msm <- as.numeric(X_msm %*% beta_hat)
  mu_msm <- fam$linkinv(eta_msm)
  mu_eta_msm <- fam$mu.eta(eta_msm)
  var_mu_msm <- fam$variance(mu_msm)
  y_msm <- stats::model.response(stats::model.frame(msm_model))
  r_msm <- y_msm - mu_msm

  # Decompose the MSM prior weights: total_w = ext_w * DR_w * ipcw_w.
  # The phi_bar_cens closure must vary only ipcw_w as gamma changes while
  # holding DR_w and ext_w fixed. So other_w = total_w / ipcw_w_hat.
  pw <- msm_model$prior.weights
  if (is.null(pw)) {
    pw <- rep(1, n_sub)
  }

  ipcw_w_fit <- fit$details$ipcw_weights[fit_rows]
  # Censored rows have ipcw_w = 0 and pw = 0; the ratio is 0/0, so
  # use ifelse to avoid NaN propagation (the contribution is 0 either way).
  other_w <- ifelse(ipcw_w_fit > 0, pw / ipcw_w_fit, 0)

  # IPCW weight closure on full data, then subset to fit_rows
  ipcw_wfn <- make_ipcw_weight_fn(
    cens_model,
    n_total = n_total,
    censoring_col = as.integer(fit$data[[fit$censoring]]),
    stabilize = TRUE
  )

  phi_bar_cens <- function(gamma) {
    w_ipcw_full <- ipcw_wfn(gamma)
    w_ipcw_sub <- w_ipcw_full[fit_rows]
    s_per_i <- other_w * w_ipcw_sub * mu_eta_msm * r_msm / var_mu_msm
    as.numeric(crossprod(X_msm, s_per_i)) / n_sub
  }

  gamma_hat <- cens_model$alpha_hat
  # A_{beta,gamma} = -d phi_bar / d gamma, the cross-block of the stacked bread.
  A_beta_gamma <- -numDeriv::jacobian(phi_bar_cens, x = gamma_hat)
  g_cens <- as.numeric(crossprod(A_beta_gamma, h_msm))

  # The censoring model was fit on the full data (n_total), so its prep
  # uses n_total as the denominator; we then extract only the fit_rows
  # slice to stay aligned with the IPW MSM's n_sub-length IF vectors.
  cens_prep <- prepare_model_if(
    cens_model$model,
    which(cens_model$fit_rows),
    n_total
  )
  cens_res <- apply_model_correction(cens_prep, g_cens)
  cens_res$correction[fit_rows]
}


#' G-computation branch of variance_if()
#'
#' @noRd
variance_if_gcomp <- function(
  fit,
  data_a_list,
  preds_list,
  mu_hat,
  target_idx,
  interventions = NULL,
  cluster_vec = NULL
) {
  model <- fit$model
  int_names <- names(data_a_list)

  has_family <- !is.null(model$family) &&
    !is.null(model$family$mu.eta) &&
    !is.null(model$family$variance)
  is_gam <- inherits(model, "gam") && !is.null(model$Vp)

  if (!has_family && !is_gam) {
    return(variance_if_numeric(
      fit,
      data_a_list,
      preds_list,
      mu_hat,
      target_idx,
      weights = fit$details$weights,
      cluster_vec = cluster_vec
    ))
  }

  bundle <- prepare_point_variance(
    fit,
    model,
    data_a_list,
    preds_list,
    mu_hat,
    target_idx,
    interventions = interventions
  )
  pieces <- bundle$pieces
  prep <- bundle$prep

  has_ipcw <- isTRUE(fit$details$ipcw)
  n_fit <- nrow(prep$X_fit)

  IF_list <- lapply(seq_along(pieces$grad_list), function(j) {
    res <- apply_model_correction(prep, pieces$grad_list[[j]])
    if_j <- pieces$Ch1_list[[j]] + res$correction

    if (has_ipcw) {
      # The IPCW cross-term needs h = A_{bb}^{-1} J in M-estimation
      # scaling. `res$h` = (X'WX)^{-1} J is the bread-projected gradient
      # in GLM scaling; A_{bb} = (1/n) X'WX so A_{bb}^{-1} = n (X'WX)^{-1},
      # hence h_outcome = n_fit * res$h.
      ipcw_corr <- compute_ipcw_if_correction(
        fit,
        model,
        prep,
        h_outcome = n_fit * res$h,
        n_fit = n_fit
      )
      if_j <- if_j - ipcw_corr$correction
    }

    if_j
  })
  names(IF_list) <- int_names

  # `cluster_vec` has length `n = nrow(fit$data)` (validated upstream by
  # `resolve_cluster()`), matching the length of each IF vector; passing
  # it to `vcov_from_if()` switches aggregation from independent-sum to
  # sum-within-cluster-then-square (Liang & Zeger 1986).
  vcov_from_if(IF_list, pieces$n, int_names, cluster = cluster_vec)
}


#' Matching branch of variance_if()
#'
#' @description
#' D1: the IF for matching is computed **on the matched sample** (not on
#' the full data). `n` is the matched sample size, `fit_idx = 1..n_match`,
#' and IFs are aggregated cluster-robustly with `cluster = subclass`. Match
#' weights enter twice: as `prior_weights` of the outcome glm (so they
#' propagate through `correct_model()` via the IWLS working weights) and
#' as the target-population weights for Channel 1.
#'
#' For a saturated MSM (`Y ~ A`) with static interventions, predictions
#' are constant within each treatment level so Channel 1 is identically
#' zero and the result reduces to the Channel-2-only formula
#' `J vcovCL(model) J^T`, matching the legacy implementation.
#'
#' **Why this branch doesn't use `build_point_channel_pieces()`.** G-comp
#' and IPW share that helper because they both operate on the full
#' `fit$data` with optional external weights. Matching operates on
#' `fit$details$matched_data` (a strict subset with its own row count
#' `n_m`), interventions are re-applied on `matched_dt` rather than
#' `data`, and match weights play a dual role -- they enter both as the
#' outcome-model `prior.weights` (propagated through `prepare_model_if()`
#' via IWLS working weights) and as the target-population weights for
#' Channel 1. Forcing this into `build_point_channel_pieces()` would
#' require a `scope = c("full", "matched")` switch whose branches would
#' be larger than the inline construction below. The duplication here
#' is intentional -- any future Channel-1 fix must be mirrored into both
#' places.
#'
#' @noRd
variance_if_matching <- function(fit, interventions) {
  model <- fit$model
  int_names <- names(interventions)
  k <- length(int_names)

  # `fit$details$matched_data` is stashed as a data.table by
  # `fit_matching()` (R/matching.R:168). Use it directly -- no round-trip
  # through data.frame. `$subclass` / `$weights` column access and
  # `nrow()` all work natively on data.table, and downstream
  # `stats::predict()` + `iv_design_matrix()` accept data.table inputs.
  matched_dt <- fit$details$matched_data
  n_m <- nrow(matched_dt)
  cluster <- matched_dt$subclass
  match_w <- matched_dt$weights
  if (is.null(match_w)) {
    match_w <- rep(1, n_m)
  }
  sum_w <- sum(match_w)

  fit_idx <- seq_len(n_m)
  beta_hat <- stats::coef(model)

  treatment <- fit$treatment

  prep <- prepare_model_if(model, fit_idx, n_m)

  IF_list <- lapply(seq_len(k), function(j) {
    iv <- interventions[[j]]
    da <- apply_intervention(matched_dt, treatment, iv)
    p <- as.numeric(stats::predict(model, newdata = da, type = "response"))
    mu_hat_j <- sum(match_w * p) / sum_w

    # Channel 1 on the matched sample: Ch1_i = n_m * (w_i / sum_w) * (pred_i - mu_hat_j).
    IF_vec <- n_m * (match_w / sum_w) * (p - mu_hat_j)

    X_star <- iv_design_matrix(model, da)
    eta_star <- as.numeric(X_star %*% beta_hat)
    mu_eta_star <- model$family$mu.eta(eta_star)
    # Gradient g_j = (X^{*T} diag(w_i * mu'(eta^*_i)) 1) / sum_w, the
    # match-weight-adjusted marginal-mean Jacobian for intervention j.
    g <- as.numeric(crossprod(X_star, match_w * mu_eta_star)) / sum_w

    # Total IF_i = Ch1_i + Ch2_i; cluster-robust aggregation done by
    # vcov_from_if() on subclass, not here.
    IF_vec + apply_model_correction(prep, g)$correction
  })
  names(IF_list) <- int_names

  vcov_from_if(IF_list, n_m, int_names, cluster = cluster)
}
