#' Stacked estimating-equation sandwich for longitudinal AIPW
#'
#' @description
#' Machinery for the full M-estimation sandwich variance of the
#' ICE-AIPW estimator (Bang & Robins 2005), valid on **both balanced
#' and unbalanced** panels. The estimator solves a stacked system of
#' estimating equations \eqn{\sum_i \psi(O_i; \theta) = 0} where
#' \eqn{\theta = (\alpha, \beta, \gamma, \mu)} concatenates (the
#' \eqn{\gamma} block is present only under IPCW):
#'
#' \itemize{
#'   \item \eqn{\alpha}: the per-period propensity coefficients (one
#'     block per time point; for multivariate treatment, one sub-block
#'     per within-period component). Estimated by the model score gated
#'     by each period's fit mask -- the GLM score \eqn{X_i (A_i -
#'     \pi_i(\alpha))} for a GLM treatment model, the softmax residual
#'     score for a multinomial (categorical) treatment model.
#'   \item \eqn{\beta}: the per-step ICE outcome / pseudo-outcome
#'     coefficients (one block per backward step). Estimated by the GLM
#'     score \eqn{X_i (\mathrm{resp}_{k,i} - m_k(\beta_k))} gated by
#'     step \eqn{k}'s observed/uncensored fit mask, with
#'     \eqn{\mathrm{resp}_k = Y} at the final step and the
#'     \emph{previous step's pseudo-outcome} otherwise.
#'   \item \eqn{\gamma} (IPCW only): the per-period censoring-model
#'     coefficients (one block per period with censoring variation).
#'     Estimated by the logistic score \eqn{X_{\mathrm{cens},i}\,
#'     (C^{0}_{i,k} - \pi_{C,k}(\gamma_k))}. The per-period stabilized IPCW
#'     weight \eqn{w_{C,k}(\gamma_k)} is threaded into each outcome score's
#'     prior weight, so the numerical bread captures the censoring
#'     cross-terms (\eqn{\gamma \to \beta \to \mu}) mechanically -- no
#'     hand-derived cross-derivative. By the double-robust orthogonality of
#'     the augmented estimator this cross-term is near-zero, but it is
#'     carried for an exactly correct sandwich.
#'   \item \eqn{\mu}: the marginal mean, estimated by
#'     \eqn{w_i t_i (\tilde Y_{0,i} - \mu)} over the target population.
#' }
#'
#' The augmented backward recursion is, for \eqn{k} from the last step
#' to the first,
#' \deqn{\tilde Y_k = m_k(d_k; \beta_k) + W_k(\alpha)\,
#'   (\tilde Y_{k+1} - m_k(A_k; \beta_k)),}
#' with \eqn{\tilde Y_{K+1} = Y}, NA propagated to prediction-only for
#' ids absent at later steps, and \eqn{\hat\mu} the (weighted) mean of
#' \eqn{\tilde Y_0} over the target.
#'
#' The bread \eqn{B = -\frac1n \partial(\sum_i\psi)/\partial\theta} is a
#' numerical Jacobian (`numDeriv::jacobian`) of the summed stacked score;
#' it captures **every** block-triangular cross-term -- including the
#' baseline-pseudo-regression block \eqn{B^{-1}[\mu,\beta_1]} that the
#' previous forward-cascade assembly dropped under selection. The meat is
#' \eqn{M = \frac1n\sum_i\psi_i\psi_i^\top}, and the per-id influence
#' function for the marginal mean is the \eqn{\mu}-row of \eqn{B^{-1}}
#' applied to each \eqn{\psi_i}:
#' \deqn{\mathrm{IF}_i = e_\mu^\top B^{-1}\psi_i.}
#'
#' Designs, weights and fit masks are extracted from causatr's own
#' fitted models (the `ice_aipw_iterate()` result and
#' `fit$details$treatment_models_by_time`) so the estimating equation is
#' faithful by construction: \eqn{\sum_i\psi(\hat\theta)} vanishes (to GLM
#' convergence) and the recursion reproduces the stored `mu_hat` exactly.
#'
#' **Supported nuisance models.** The reconstructed score is the
#' (maximum-likelihood) GLM / multinomial score, whose Jacobian is the
#' model information -- so the analytic sandwich covers GLM-family outcome
#' and propensity models (gaussian, binomial, poisson, Gamma, quasi-*,
#' inverse-gaussian, `MASS::glm.nb`) and multinomial (categorical)
#' propensities. Penalised / non-likelihood fitters (`mgcv::gam`,
#' `betareg::betareg`) do not have a vanishing GLM score whose Jacobian is
#' the sandwich bread, so they are routed to the bootstrap (consistent
#' with the longitudinal ICE betareg path).
#'
#' Composition follows Zivich et al. (2024, *Stat. Med.* 43:5562-5572)
#' for the ICE outcome-model chain and Shook-Sa et al. (2025,
#' *Biometrics* 81(2):ujaf054) for the propensity correction.
#'
#' @name variance_if_aipw_longitudinal_ee
#' @keywords internal
NULL


#' Is a nuisance model supported by the longitudinal-AIPW analytic sandwich?
#'
#' @description
#' The stacked-EE sandwich reconstructs each nuisance model's score from
#' the maximum-likelihood GLM / multinomial form. GLMs (including
#' `MASS::glm.nb`, whose family carries the estimated dispersion) and
#' `nnet::multinom` qualify. `mgcv::gam` is excluded even though it
#' inherits `"glm"`: its penalised fit's bread is `Vp`, not the score
#' Jacobian. `betareg` and other non-GLM fitters are excluded too.
#'
#' @param model A fitted nuisance model.
#'
#' @returns `TRUE` if the analytic sandwich can score `model`.
#'
#' @examples
#' m <- stats::glm(c(0, 1, 0, 1) ~ c(1, 2, 3, 4), family = stats::binomial())
#' causatr:::aipw_long_model_supported(m)
#'
#' @noRd
aipw_long_model_supported <- function(model) {
  if (inherits(model, "gam")) {
    return(FALSE)
  }
  inherits(model, "glm") || inherits(model, "multinom")
}


#' Per-observation GLM score matrix as a function of coefficients
#'
#' @description
#' Returns the \eqn{n_{\mathrm{fit}} \times p} matrix whose row \eqn{i}
#' is the GLM log-likelihood gradient
#' \eqn{X_i\, w_i\, \mu'(\eta_i)\, (y_i - \mu_i) / V(\mu_i)} evaluated at
#' the supplied coefficient vector. Column-sums vanish at the fitted MLE
#' (the GLM normal equations), which is exactly the faithfulness the
#' stacked sandwich relies on for both the propensity scores (treatment
#' response, fixed) and -- after substituting the recursion-dependent
#' response -- the outcome scores.
#'
#' For continuous-treatment propensities the gaussian GLM score
#' \eqn{X(A - X\alpha)} differs from the density score
#' \eqn{X(A - X\alpha)/\sigma^2} by the constant \eqn{\sigma^2}; that
#' factor cancels in \eqn{A_{\alpha\alpha}^{-1}\psi_\alpha}, so the
#' influence function is invariant to it.
#'
#' @param model A fitted GLM (`stats::glm` or compatible) with a
#'   `family` object.
#' @param alpha Numeric coefficient vector aligned to the model's
#'   non-aliased columns.
#'
#' @returns A numeric matrix (`n_fit x p`), one score row per fit
#'   observation.
#'
#' @examples
#' m <- stats::glm(c(0, 1, 1, 0) ~ c(-1, 0, 1, 2), family = stats::binomial())
#' s <- causatr:::glm_score_matrix(m, stats::coef(m))
#' # Column sums vanish at the MLE.
#' max(abs(colSums(s))) < 1e-6
#'
#' @noRd
glm_score_matrix <- function(model, alpha) {
  X <- stats::model.matrix(model)
  aliased <- is.na(stats::coef(model))
  if (any(aliased)) {
    X <- X[, !aliased, drop = FALSE]
  }
  y <- stats::model.response(stats::model.frame(model))
  fam <- model$family
  eta <- as.numeric(X %*% alpha)
  mu <- fam$linkinv(eta)
  mu_eta <- fam$mu.eta(eta)
  v <- fam$variance(mu)
  pw <- model$prior.weights
  if (is.null(pw)) {
    pw <- rep(1, length(eta))
  }
  fac <- pw * mu_eta * (y - mu) / v
  X * fac
}


#' Per-observation multinomial-logit score matrix as a function of coefficients
#'
#' @description
#' Softmax residual score for a `nnet::multinom` treatment model, used for
#' a categorical (k > 2) propensity in the stacked sandwich. For each
#' non-reference class \eqn{k} the score block is \eqn{X_i\,(I(A_i = k) -
#' p_{ik}(\alpha))}, stacked class-major to match `prepare_model_if_multinom()`
#' and `nnet::multinom`'s coefficient layout (so the column-sums vanish at
#' the fitted coefficients).
#'
#' @param model A fitted `nnet::multinom` model.
#' @param alpha Numeric coefficient vector, length `(K-1) * p`, laid out
#'   class-major (`c(t(coef(model)))`).
#'
#' @returns A numeric matrix (`n_fit x (K-1)*p`), one stacked score row
#'   per fit observation.
#'
#' @examples
#' if (requireNamespace("nnet", quietly = TRUE)) {
#'   d <- data.frame(y = factor(rep(c("a", "b", "c"), 20)), x = rnorm(60))
#'   m <- nnet::multinom(y ~ x, data = d, trace = FALSE)
#'   s <- causatr:::multinom_score_matrix(m, as.numeric(t(stats::coef(m))))
#'   max(abs(colSums(s))) < 1e-4
#' }
#'
#' @noRd
multinom_score_matrix <- function(model, alpha) {
  X <- stats::model.matrix(model)
  n <- nrow(X)
  p <- ncol(X)
  lev <- model$lev
  km1 <- length(lev) - 1L

  # alpha is class-major: [class2 coefs (p), class3 coefs (p), ...].
  # matrix(alpha, p, km1) fills column k with class (k+1)'s coefficients.
  Amat <- matrix(alpha, nrow = p, ncol = km1)
  eta <- X %*% Amat
  # Reference class has linear predictor 0; softmax over all K classes.
  expe <- cbind(1, exp(eta))
  probs <- expe / rowSums(expe)
  p_non_ref <- probs[, -1, drop = FALSE]

  resp <- as.character(stats::model.response(stats::model.frame(model)))
  non_ref <- lev[-1]

  s <- matrix(0, nrow = n, ncol = km1 * p)
  for (k in seq_len(km1)) {
    r_k <- as.numeric(resp == non_ref[k]) - p_non_ref[, k]
    cols <- ((k - 1L) * p + 1L):(k * p)
    s[, cols] <- X * r_k
  }
  s
}


#' Per-observation score matrix for a nuisance model (GLM or multinomial)
#'
#' @description
#' Dispatches to `glm_score_matrix()` or `multinom_score_matrix()` so the
#' stacked sandwich treats GLM and categorical propensities uniformly.
#'
#' @param model A fitted GLM or `nnet::multinom` model.
#' @param alpha Numeric coefficient vector for `model`.
#'
#' @returns A numeric score matrix (`n_fit x p`).
#'
#' @noRd
nuisance_score_matrix <- function(model, alpha) {
  if (inherits(model, "multinom")) {
    multinom_score_matrix(model, alpha)
  } else {
    glm_score_matrix(model, alpha)
  }
}


#' Assemble the stacked-EE pieces for one longitudinal-AIPW intervention
#'
#' @description
#' Extracts fixed designs / masks / responses from causatr's fitted
#' propensity and outcome models and returns a closure `psi_fn(theta)`
#' that recomputes the per-id stacked score matrix as a function of
#' \eqn{\theta = (\alpha, \beta, \mu)}, together with the fitted
#' \eqn{\hat\theta} and the index of the \eqn{\mu} component. Designs are
#' coefficient-independent and so are precomputed once; `psi_fn` is then
#' pure linear algebra plus the per-period weight closures, which keeps
#' the numerical Jacobian fast.
#'
#' For the natural-course arm (`intervention = NULL`) the per-period
#' weights are identically 1 and carry no propensity uncertainty, so the
#' propensity block is omitted entirely.
#'
#' @param fit `causatr_fit` (`estimator = "aipw"`, longitudinal).
#' @param aipw_result One element of `ice_aipw_iterate()` outputs.
#' @param all_ids Character vector of first-time-point ids (the id order).
#' @param n Integer. Number of individuals.
#' @param id_to_idx Named integer map from id to `1..n`.
#' @param target Logical vector (length `n`) flagging the target.
#' @param w_t Numeric weights aligned to `which(target)` (external
#'   weights on the target, or all ones when unweighted).
#' @param sum_w_target Numeric. `sum(w_t)`.
#' @param trim Numeric upper-quantile weight truncation (1 = none).
#'
#' @returns A list with `psi_fn`, `theta_hat`, `mu_index`,
#'   `total_alpha`, `total_beta`, `total_gamma` (the censoring block size, 0
#'   without IPCW), and `unsupported` (NULL, or a string naming an unsupported
#'   nuisance model class that should route to the bootstrap).
#'
#' @noRd
build_aipw_long_psi <- function(
  fit,
  aipw_result,
  all_ids,
  n,
  id_to_idx,
  target,
  w_t,
  sum_w_target,
  trim = 1
) {
  data <- fit$data
  details <- fit$details
  outcome <- fit$outcome
  id_col <- fit$id
  time_col <- fit$time
  censoring <- fit$censoring
  time_points <- details$time_points
  n_times <- details$n_times
  binary_outcome <- is_binary_family(details$family_outcome)
  treatment_models_by_time <- details$treatment_models_by_time
  fit_data_by_time <- details$fit_data_by_time

  models <- aipw_result$models
  data_iv <- aipw_result$data_iv
  fit_ids_list <- aipw_result$fit_ids
  intervention <- aipw_result$intervention

  uncens <- is_uncensored(data, censoring)
  target_idx <- which(target)

  unsupported <- NULL

  # ---- Propensity side: per-period weight closures + score terms ----
  # Reuse the same per-period closure construction as the longitudinal
  # IPW sandwich so the perturbed weights W(alpha) reproduce causatr's
  # density ratios exactly. `prop_terms` carries one (model, design,
  # id-positions, alpha column range) record per period (univariate) or
  # per within-period component (multivariate); its scores supply the
  # alpha block of psi. The natural-course arm has W == 1 identically, so
  # it carries no propensity parameters.
  K <- length(treatment_models_by_time)
  is_mv <- inherits(treatment_models_by_time[[1L]], "causatr_treatment_models")

  sub_fns <- vector("list", K)
  align_idx_list <- vector("list", K)
  alpha_blocks <- vector("list", K)
  block_lens <- integer(K)
  prop_terms <- list()
  has_prop <- !is.null(intervention)

  if (has_prop) {
    for (k in seq_len(K)) {
      tm_k <- treatment_models_by_time[[k]]
      data_k <- fit_data_by_time[[k]]
      ids_k <- as.character(data_k[[id_col]])

      if (is_mv) {
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
        align_idx_list[[k]] <- match(ids_k, all_ids)
        alpha_blocks[[k]] <- mv_closure_k$alpha_hat
        block_lens[k] <- length(mv_closure_k$alpha_hat)

        # Per-component score terms from the original (un-reset) component
        # models, whose fit rows carry the true NA handling per component.
        # Only build them when this period actually contributes alpha
        # (a fully natural-course component leaves block_lens == 0).
        if (block_lens[k] > 0L) {
          comp_lens <- vapply(tm_k, function(m) length(m$alpha_hat), integer(1))
          comp_offsets <- c(0L, cumsum(comp_lens))
          for (j in seq_along(tm_k)) {
            if (comp_lens[j] == 0L) {
              next
            }
            model_j <- tm_k[[j]]$model
            if (!aipw_long_model_supported(model_j)) {
              unsupported <- class(model_j)[1L]
            }
            comp_ids_j <- ids_k[tm_k[[j]]$fit_rows]
            prop_terms[[length(prop_terms) + 1L]] <- list(
              model = model_j,
              pos = match(comp_ids_j, all_ids),
              local_cols = (comp_offsets[j] + 1L):comp_offsets[j + 1L],
              period = k
            )
          }
        }
      } else {
        sub_fns[[k]] <- make_weight_fn(
          treatment_model = tm_k,
          data = data_k,
          intervention = intervention,
          estimand = "ATE",
          trim = trim
        )
        period_ids <- ids_k[tm_k$fit_rows]
        align_idx_list[[k]] <- match(period_ids, all_ids)
        alpha_blocks[[k]] <- tm_k$alpha_hat
        block_lens[k] <- length(tm_k$alpha_hat)
        if (block_lens[k] > 0L) {
          if (!aipw_long_model_supported(tm_k$model)) {
            unsupported <- class(tm_k$model)[1L]
          }
          prop_terms[[length(prop_terms) + 1L]] <- list(
            model = tm_k$model,
            pos = match(period_ids, all_ids),
            local_cols = seq_len(block_lens[k]),
            period = k
          )
        }
      }
    }
  }

  alpha_offsets <- c(1L, cumsum(block_lens) + 1L)
  alpha_hat_stacked <- unlist(alpha_blocks, use.names = FALSE)
  if (is.null(alpha_hat_stacked)) {
    alpha_hat_stacked <- numeric(0)
  }
  total_alpha <- length(alpha_hat_stacked)

  # Resolve each propensity term's global alpha column indices from its
  # period offset and within-period local columns.
  for (t_i in seq_along(prop_terms)) {
    pk <- prop_terms[[t_i]]$period
    prop_terms[[t_i]]$alpha_cols <-
      alpha_offsets[pk] - 1L + prop_terms[[t_i]]$local_cols
  }

  # ---- Censoring side (IPCW only): per-period gamma block ------------
  # Under IPCW the per-step outcome models are fit weighted by the per-period
  # stabilized IPCW weight (carried in `details$weights`). Including the
  # censoring coefficients gamma in theta -- with the per-period logistic score
  # pinning gamma and the IPCW weight threaded into each outcome score as a
  # function of gamma -- lets the numerical bread capture the censoring
  # cross-terms (gamma -> outcome beta -> mu) mechanically. `cens_block`'s
  # gamma columns are appended after beta; `cens_score_terms` map each period's
  # logistic score to first-period id positions, mirroring `prop_terms`.
  has_ipcw <- isTRUE(details$ipcw)
  cens_block <- if (has_ipcw) make_ipcw_weight_fn_longitudinal(fit) else NULL
  total_gamma <- if (has_ipcw) cens_block$n_gamma else 0L
  weights_pre_ipcw <- details$weights_pre_ipcw
  cens_score_terms <- list()
  if (has_ipcw) {
    for (b in cens_block$blocks) {
      ids_b <- as.character(data[[id_col]][b$score_rows_global])
      cens_score_terms[[length(cens_score_terms) + 1L]] <- list(
        model = b$model,
        pos = match(ids_b, all_ids),
        # gamma columns are offset past the alpha and beta blocks below once
        # total_beta is known; store the within-gamma local columns for now.
        local_cols = b$gamma_cols
      )
    }
  }

  # ---- Outcome side: per-step fixed designs and fit masks -----------
  # Step ordering for the theta vector is ascending (step 1 .. n_times);
  # the recursion below processes descending. NULL models contribute no
  # block.
  beta_blocks <- vector("list", n_times)
  steps <- vector("list", n_times)
  beta_offset <- 0L

  for (step_i in seq_len(n_times)) {
    model_k <- models[[step_i]]
    if (is.null(model_k)) {
      next
    }
    if (!aipw_long_model_supported(model_k)) {
      unsupported <- class(model_k)[1L]
    }
    beta_k <- coef_clean(model_k)
    p_k <- length(beta_k)
    beta_blocks[[step_i]] <- beta_k

    current_time <- time_points[step_i]
    is_final <- step_i == n_times

    # Observed-treatment fit design: model.matrix gives the exact rows
    # (in fit order) the step was estimated on; aligns with prior.weights
    # and the recursion-dependent response.
    X_obs_fit <- stats::model.matrix(model_k)
    aliased <- is.na(stats::coef(model_k))
    if (any(aliased)) {
      X_obs_fit <- X_obs_fit[, !aliased, drop = FALSE]
    }
    prior_w_fit <- model_k$prior.weights
    if (is.null(prior_w_fit)) {
      prior_w_fit <- rep(1, nrow(X_obs_fit))
    }
    fit_ids_k <- fit_ids_list[[step_i]]
    fit_idx <- unname(id_to_idx[fit_ids_k])

    # IPCW weight threading: the step's fit-row global person-period indices
    # (in model.matrix / fit_ids_k order) and the survey-weight part of
    # `prior_w_fit` that does NOT vary with gamma. `prior_w_fit = external x
    # ipcw(gamma)`, so the perturbed score uses `external x ipcw_full(gamma)`
    # at these rows. With no survey weights `external` is all ones.
    fit_global_rows <- NULL
    external_part <- NULL
    if (has_ipcw) {
      period_rows <- which(data[[time_col]] == current_time)
      period_ids <- as.character(data[[id_col]][period_rows])
      fit_global_rows <- period_rows[match(fit_ids_k, period_ids)]
      external_part <- if (is.null(weights_pre_ipcw)) {
        rep(1, length(fit_global_rows))
      } else {
        weights_pre_ipcw[fit_global_rows]
      }
    }

    if (is_final) {
      pred_mask <- (data[[time_col]] == current_time) & uncens
      y_fit <- stats::model.response(stats::model.frame(model_k))
    } else {
      pred_mask <- data[[time_col]] == current_time
      y_fit <- NULL
    }
    pred_data <- data[pred_mask]
    pred_data_iv <- data_iv[pred_mask]
    pred_idx <- unname(id_to_idx[as.character(pred_data[[id_col]])])

    X_obs_pred <- iv_design_matrix(model_k, pred_data)
    X_iv_pred <- iv_design_matrix(model_k, pred_data_iv)
    y_pred <- if (is_final) pred_data[[outcome]] else NULL

    steps[[step_i]] <- list(
      beta_cols = (beta_offset + 1L):(beta_offset + p_k),
      family = model_k$family,
      is_final = is_final,
      X_obs_fit = X_obs_fit,
      prior_w_fit = prior_w_fit,
      fit_idx = fit_idx,
      fit_global_rows = fit_global_rows,
      external_part = external_part,
      y_fit = y_fit,
      X_obs_pred = X_obs_pred,
      X_iv_pred = X_iv_pred,
      pred_idx = pred_idx,
      y_pred = y_pred,
      period = step_i
    )
    beta_offset <- beta_offset + p_k
  }
  total_beta <- beta_offset
  beta_hat_stacked <- unlist(beta_blocks, use.names = FALSE)
  if (is.null(beta_hat_stacked)) {
    beta_hat_stacked <- numeric(0)
  }

  gamma_hat <- if (has_ipcw) cens_block$gamma_hat else numeric(0)

  # theta = (alpha, beta, gamma, mu); the gamma block sits between the outcome
  # coefficients and the marginal mean so the bread's block-triangular cross
  # terms (gamma -> beta -> mu) are captured by the numerical Jacobian.
  mu_index <- total_alpha + total_beta + total_gamma + 1L
  mu_hat <- sum(w_t * aipw_result$pseudo_final[target_idx]) / sum_w_target
  theta_hat <- c(alpha_hat_stacked, beta_hat_stacked, gamma_hat, mu_hat)

  clip_lo <- 1e-5
  clip_hi <- 1 - 1e-5

  step_order <- rev(which(!vapply(steps, is.null, logical(1))))

  psi_fn <- function(theta) {
    alpha <- theta[seq_len(total_alpha)]
    beta_all <- theta[total_alpha + seq_len(total_beta)]
    gamma <- theta[total_alpha + total_beta + seq_len(total_gamma)]
    mu <- theta[mu_index]

    # Full-length per-row IPCW weight under the perturbed gamma. Each outcome
    # step's prior weight is then `external x ipcw(gamma)` at its fit rows, so
    # the bread differentiates the outcome score through the censoring model.
    w_ipcw_full <- if (has_ipcw) cens_block$weight_full_fn(gamma) else NULL

    # Per-period density-ratio weights under the perturbed alpha.
    W_new <- matrix(1, nrow = n, ncol = K)
    for (kk in seq_len(K)) {
      if (block_lens[kk] == 0L) {
        next
      }
      idx <- alpha_offsets[kk]:(alpha_offsets[kk + 1L] - 1L)
      w_raw <- sub_fns[[kk]](alpha[idx])
      w_full <- rep(1, n)
      w_full[align_idx_list[[kk]]] <- w_raw
      W_new[, kk] <- w_full
    }

    # Backward augmented recursion; accumulate per-step outcome scores.
    pseudo <- rep(NA_real_, n)
    pseudo_reg <- rep(NA_real_, n)
    psi_beta <- matrix(0, nrow = n, ncol = total_beta)

    for (step_i in step_order) {
      st <- steps[[step_i]]
      fam <- st$family
      beta_k <- beta_all[st$beta_cols]

      # Outcome score on the step's fit rows FIRST, while `pseudo_reg`
      # still holds the later step's pseudo-outcome -- that (clipped)
      # value is exactly the response the step was fit on. The response
      # is Y at the final step. Computing the score before the pseudo
      # update is essential: the fit rows are a subset of this step's
      # prediction rows, so updating `pseudo_reg` first would overwrite
      # the regression response. This block-triangular coupling is what
      # the bread's numerical Jacobian differentiates through.
      eta_fit <- as.numeric(st$X_obs_fit %*% beta_k)
      mu_fit <- fam$linkinv(eta_fit)
      mu_eta_fit <- fam$mu.eta(eta_fit)
      v_fit <- fam$variance(mu_fit)
      resp_fit <- if (st$is_final) st$y_fit else pseudo_reg[st$fit_idx]
      # Recompute the outcome model's prior weight at the perturbed gamma:
      # `external x ipcw(gamma)` on the step's fit rows. At gamma_hat this is
      # exactly the fitted `prior_w_fit`, so the score still vanishes (the
      # faithfulness gate holds). With no IPCW the stored weight is used.
      pw_fit <- if (has_ipcw) {
        st$external_part * w_ipcw_full[st$fit_global_rows]
      } else {
        st$prior_w_fit
      }
      fac <- pw_fit * mu_eta_fit * (resp_fit - mu_fit) / v_fit
      psi_beta[st$fit_idx, st$beta_cols] <- st$X_obs_fit * fac

      # Augmented pseudo-outcome update on the prediction rows, producing
      # the response for the next (earlier) step.
      eta_iv <- as.numeric(st$X_iv_pred %*% beta_k)
      mu_iv <- fam$linkinv(eta_iv)
      eta_obs <- as.numeric(st$X_obs_pred %*% beta_k)
      mu_obs <- fam$linkinv(eta_obs)
      w_step <- W_new[st$pred_idx, st$period]

      if (st$is_final) {
        resid <- st$y_pred - mu_obs
        new_pseudo <- mu_iv + w_step * resid
      } else {
        prev_pseudo <- pseudo[st$pred_idx]
        has_prev <- !is.na(prev_pseudo)
        new_pseudo <- mu_iv
        new_pseudo[has_prev] <- mu_iv[has_prev] +
          w_step[has_prev] * (prev_pseudo[has_prev] - mu_obs[has_prev])
      }
      pseudo[st$pred_idx] <- new_pseudo
      pseudo_reg[st$pred_idx] <- if (binary_outcome) {
        pmax(pmin(new_pseudo, clip_hi), clip_lo)
      } else {
        new_pseudo
      }
    }

    # Propensity scores (treatment response fixed; depend on alpha only).
    psi_alpha <- matrix(0, nrow = n, ncol = total_alpha)
    for (term in prop_terms) {
      s <- nuisance_score_matrix(term$model, alpha[term$alpha_cols])
      psi_alpha[term$pos, term$alpha_cols] <- s
    }

    # Per-period censoring scores (uncensoring response fixed; depend on gamma
    # only). The logistic score pins gamma in the stacked system so the bread
    # picks up the gamma cross-terms.
    psi_gamma <- matrix(0, nrow = n, ncol = total_gamma)
    for (term in cens_score_terms) {
      s <- glm_score_matrix(term$model, gamma[term$local_cols])
      psi_gamma[term$pos, term$local_cols] <- s
    }

    # Marginal-mean estimating equation over the target population.
    psi_mu <- numeric(n)
    psi_mu[target_idx] <- w_t * (pseudo[target_idx] - mu)

    cbind(psi_alpha, psi_beta, psi_gamma, psi_mu)
  }

  list(
    psi_fn = psi_fn,
    theta_hat = theta_hat,
    mu_index = mu_index,
    total_alpha = total_alpha,
    total_beta = total_beta,
    total_gamma = total_gamma,
    unsupported = unsupported
  )
}


#' Per-id influence function for one longitudinal-AIPW intervention
#'
#' @description
#' Builds the stacked estimating equation via `build_aipw_long_psi()`,
#' checks faithfulness (the summed score must vanish at \eqn{\hat\theta}),
#' forms the numerical bread, and returns the per-id influence function
#' \eqn{\mathrm{IF}_i = e_\mu^\top B^{-1}\psi_i} for the marginal mean.
#'
#' @inheritParams build_aipw_long_psi
#'
#' @returns A numeric vector of length `n` (the per-id IF), or signals a
#'   classed error if a nuisance model is unsupported or faithfulness
#'   fails (both steer to the bootstrap).
#'
#' @noRd
aipw_long_stacked_if <- function(
  fit,
  aipw_result,
  all_ids,
  n,
  id_to_idx,
  target,
  w_t,
  sum_w_target,
  trim = 1
) {
  built <- build_aipw_long_psi(
    fit,
    aipw_result,
    all_ids,
    n,
    id_to_idx,
    target,
    w_t,
    sum_w_target,
    trim = trim
  )

  if (!is.null(built$unsupported)) {
    rlang::abort(
      c(
        paste0(
          "Longitudinal AIPW sandwich variance does not support a `",
          built$unsupported,
          "` nuisance model (penalised / non-likelihood fit)."
        ),
        i = "Use `ci_method = \"bootstrap\"` for valid inference here."
      ),
      class = "causatr_longitudinal_aipw_sandwich_model"
    )
  }

  psi_fn <- built$psi_fn
  theta_hat <- built$theta_hat
  mu_index <- built$mu_index

  psi_hat <- psi_fn(theta_hat)

  # Faithfulness gate: the summed estimating equation must vanish at the
  # fitted theta (to model convergence). A gross violation means the
  # reconstructed EE does not match the fitted nuisances (an unsupported
  # configuration), where returning a sandwich SE would be silently
  # wrong -- abort and steer to the bootstrap instead. The 1e-3 tolerance is
  # set by the loosest-converging supported block: GLM/glm.nb scores vanish to
  # ~1e-10 (IRLS), but a `nnet::multinom` propensity only reaches ~1e-4
  # (its BFGS `reltol`), so the threshold leaves that block ~7x of headroom
  # while still catching an order-of-magnitude reconstruction mismatch.
  col_sums <- colSums(psi_hat)
  col_scale <- sqrt(colSums(psi_hat^2))
  ratio <- abs(col_sums) / pmax(col_scale, 1e-8)
  if (max(ratio) > 1e-3) {
    rlang::abort(
      c(
        paste0(
          "Longitudinal AIPW stacked estimating equation did not vanish ",
          "at the fitted estimate (max relative score ",
          format(max(ratio), digits = 3),
          ")."
        ),
        i = "Use `ci_method = \"bootstrap\"` for valid inference here."
      ),
      class = "causatr_longitudinal_aipw_sandwich_unfaithful"
    )
  }

  # Bread B = -(1/n) d(sum psi)/d(theta) via numerical Jacobian. The
  # block-triangular cross-terms (alpha -> beta -> mu) are captured
  # automatically, so the assembly is consistent on balanced and
  # unbalanced panels alike.
  psi_sum <- function(theta) colSums(psi_fn(theta))
  B <- -numDeriv::jacobian(psi_sum, theta_hat) / n

  # mu-row of B^{-1}: r = B^{-T} e_mu, so IF_i = r^T psi_i.
  e_mu <- numeric(length(theta_hat))
  e_mu[mu_index] <- 1
  r <- tryCatch(
    solve(t(B), e_mu),
    error = function(e) {
      rlang::abort(
        c(
          "Longitudinal AIPW sandwich bread matrix is singular.",
          i = "Use `ci_method = \"bootstrap\"` for valid inference here."
        ),
        class = "causatr_longitudinal_aipw_sandwich_singular"
      )
    }
  )

  as.numeric(psi_hat %*% r)
}
