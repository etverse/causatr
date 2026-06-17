#' Analytic IF sandwich for a multinomial-outcome point g-computation
#'
#' @description
#' Per-class influence-function variance for \eqn{P(Y = k \mid do(A = a))} from a
#' fitted `nnet::multinom` outcome model. Returns the same per-class structure as
#' `variance_bootstrap_multinom()` -- a named list (one entry per outcome class)
#' of `k_int x k_int` matrices over the interventions -- so the contrast / S3
#' layer is reused unchanged; only the variance source differs.
#'
#' The estimand for class `k` under intervention `a` is the marginal class
#' probability
#' \deqn{\mu_{k,a} = \frac{\sum_{i \in \mathrm{target}} w_i\, p_{k,i}(a)}
#'                        {\sum_{i \in \mathrm{target}} w_i},}
#' where \eqn{p_{k,i}(a) = \mathrm{softmax}_k(X_i^* \beta)} is the predicted
#' probability at the counterfactual design \eqn{X_i^*} and \eqn{w_i} is the
#' survey / external weight (all 1 for the complete-case fit, which recovers the
#' \eqn{1/n_t} average). Its influence function decomposes into the usual two
#' channels:
#' \itemize{
#'   \item \strong{Channel 1} (sampling):
#'     \eqn{(n/\sum w)\, w_i\,(p_{k,i}(a) - \mu_{k,a})} on target rows, zero
#'     elsewhere -- the per-class analogue of the Channel-1 term of
#'     `build_point_channel_pieces()`.
#'   \item \strong{Channel 2} (nuisance): the weighted multinomial-logit score
#'     correction \eqn{g_{k,a}^\top H^{-1} s_i}, where \eqn{s_i} is the stacked
#'     per-observation weighted score and \eqn{H} the weighted stacked
#'     information, both supplied by `prepare_model_if_multinom()`, and
#'     \eqn{g_{k,a} = \partial\mu_{k,a}/\partial\beta} is the softmax
#'     marginal-mean gradient (the one new derivation; see Details).
#' }
#'
#' @details
#' Stack the K-1 non-reference linear predictors as \eqn{\beta = (\beta_2,
#' \dots, \beta_K)} (the reference class `l = 1` has \eqn{\eta_1 = 0}). The
#' softmax marginal-mean gradient stacks, over the non-reference classes
#' \eqn{l = 1,\dots,K-1} (level `class_labels[l+1]`),
#' \deqn{g_{k,a}^{(l)} = \frac{1}{\sum w}\sum_{i \in \mathrm{target}}
#'        w_i\, p_{k,i}(a)\,(\delta_{k,\,l+1} - p_{l+1,i}(a))\,X_i^*,}
#' the softmax Jacobian \eqn{\partial p_k/\partial\eta_l =
#' p_k(\delta_{kl}-p_l)} chained through \eqn{\partial\eta_l/\partial\beta_l =
#' X_i^*}, weighted by the survey weight \eqn{w_i} (the denominator
#' \eqn{\sum w} is the constant total weight, so it does not contribute a
#' gradient term). The reference class (\eqn{k = 1}) is the \eqn{\delta = 0}
#' special case and falls out of the same expression, so no class is handled
#' separately. This mirrors the scalar gradient \eqn{g_j = E[\mu'(\eta) X]} of
#' `build_point_channel_pieces()`, with the \eqn{1/\sum w} factor carried in `g`
#' (not in `apply_model_correction()`, which applies the M-estimation
#' \eqn{\times n} rescale; Stefanski & Boos 2002). With `weights = NULL` every
#' \eqn{w_i = 1}, \eqn{\sum w = n_t}, and the whole derivation collapses to the
#' complete-case sandwich byte-for-byte.
#'
#' IPCW multinomial sandwiches are gated upstream in
#' `compute_contrast_multinom()` (raising
#' `causatr_categorical_outcome_sandwich`) until that chunk lands.
#'
#' @param fit A `causatr_fit` with a `nnet::multinom` outcome model.
#' @param data_a_list Named list of counterfactual data.tables (one per
#'   intervention; `apply_intervention()` already applied).
#' @param preds_list Named list of `n x K` predicted class-probability matrices,
#'   columns in `class_labels` order.
#' @param mu_mat `k_int x K` matrix of marginal class probabilities, rows named
#'   by intervention, columns by class. Already weighted by `weights` when
#'   present (the caller standardises with `maybe_weighted_mean()`).
#' @param target_idx Logical vector (length `n`) flagging target-population rows
#'   (already NA-pruned by the caller).
#' @param class_labels Character vector of the K outcome class labels (the column
#'   order of `preds_list`).
#' @param weights Numeric vector (length `n`) of survey / external weights, or
#'   `NULL` for the unweighted complete-case sandwich.
#'
#' @return A named list (one entry per class label) of `k_int x k_int`
#'   variance-covariance matrices over the interventions.
#'
#' @seealso `variance_if_gcomp()` for the scalar analogue,
#'   `prepare_model_if_multinom()` for the reused (now weighted) score/bread
#'   template.
#' @family variance
#' @noRd
variance_if_gcomp_multinom <- function(
  fit,
  data_a_list,
  preds_list,
  mu_mat,
  target_idx,
  class_labels,
  weights = NULL
) {
  model <- fit$model
  int_names <- names(data_a_list)
  k_int <- length(int_names)
  K <- length(class_labels)
  Km1 <- K - 1L
  n <- nrow(fit$data)
  n_target <- sum(target_idx)

  # Defensive: an empty target population would propagate NaNs into the vcov.
  # `compute_contrast_multinom()` prunes NA-prediction rows upstream, but a
  # classed abort here keeps the failure at the variance-engine boundary.
  if (n_target == 0L) {
    rlang::abort(
      "variance_if_gcomp_multinom(): target population is empty.",
      class = "causatr_empty_target"
    )
  }

  # Survey / external weights. The complete-case path passes `weights = NULL`,
  # which collapses every weighted expression below to the unweighted form:
  # `w_target` is all-ones and `sum_w` is `n_target`, so Channel 1's `(n/sum_w)`
  # becomes `n/n_target` and the gradient denominator becomes `n_target`.
  has_weights <- !is.null(weights)
  w_target <- if (has_weights) weights[target_idx] else rep(1, n_target)
  sum_w <- sum(w_target)
  if (sum_w <= 0) {
    rlang::abort(
      "variance_if_gcomp_multinom(): target-population weights sum to 0.",
      class = "causatr_empty_target"
    )
  }

  # Stacked weighted multinomial-logit score + inverse information. The model
  # class is identical to the multinomial *propensity* model, so the propensity
  # IF template is reused verbatim; the survey weights enter via the fit-row
  # prior weight vector (`prepare_model_if_multinom()` weights both the bread
  # and the score residual).
  fit_idx <- resolve_fit_idx(fit, model)
  w_fit <- if (has_weights) weights[fit_idx] else NULL
  prep <- prepare_model_if_multinom(model, fit_idx, n, weights = w_fit)
  # Design columns per non-reference class block. `prep$X_fit` is the
  # n_fit x ((K-1) * p) stacked design, so dividing by K-1 recovers p.
  p <- ncol(prep$X_fit) %/% Km1

  # Counterfactual design matrices are built directly from the model terms
  # (not via `iv_design_matrix()`): a multinom `coef()` is a matrix, so that
  # helper's vector-coef aliased-column drop does not apply here.
  pred_terms <- stats::delete.response(stats::terms(model))
  xlev <- model$xlevels

  # Per intervention, cache the target-row counterfactual design X^* and the
  # target-row probability matrix P (n_t x K); both are reused across all K
  # class gradients, so building them once avoids K-fold recomputation.
  Xstar_list <- vector("list", k_int)
  P_list <- vector("list", k_int)
  for (a in seq_len(k_int)) {
    da_t <- data_a_list[[a]][target_idx, , drop = FALSE]
    X_star <- stats::model.matrix(pred_terms, data = da_t, xlev = xlev)
    if (ncol(X_star) != p) {
      # A counterfactual design with a different column count than the fitted
      # model breaks the stacked-gradient alignment; surface it as a classed
      # variance error rather than a cryptic dimension mismatch downstream.
      rlang::abort(
        paste0(
          "variance_if_gcomp_multinom(): counterfactual design has ",
          ncol(X_star),
          " columns but the fitted model has ",
          p,
          "."
        ),
        class = "causatr_variance_row_mismatch"
      )
    }
    Xstar_list[[a]] <- X_star
    P_list[[a]] <- preds_list[[a]][target_idx, , drop = FALSE]
  }

  vcov_by_class <- vector("list", K)
  names(vcov_by_class) <- class_labels

  for (ci in seq_len(K)) {
    IF_list <- vector("list", k_int)
    for (a in seq_len(k_int)) {
      X_star <- Xstar_list[[a]]
      P <- P_list[[a]]
      p_k <- P[, ci] # n_t vector: predicted prob of class ci
      P_nonref <- P[, -1L, drop = FALSE] # n_t x (K-1), non-reference classes

      # Softmax marginal-mean gradient for class ci under intervention a:
      # block l (non-reference class l+1) is
      #   (1/sum_w) sum_i w_i p_{k,i} (delta_{k,l+1} - p_{l+1,i}) X_i^*.
      # With unit weights sum_w = n_t, recovering the complete-case gradient.
      g <- numeric(Km1 * p)
      for (l in seq_len(Km1)) {
        delta <- as.numeric((l + 1L) == ci)
        # n_t softmax-Jacobian weights, scaled by the per-row survey weight.
        wfac <- w_target * p_k * (delta - P_nonref[, l])
        block <- as.numeric(crossprod(X_star, wfac)) / sum_w
        g[((l - 1L) * p + 1L):(l * p)] <- block
      }

      # Channel 2: nuisance correction. `apply_model_correction()` applies the
      # M-estimation x n rescale; prep$r_score carries the prior weight w_i
      # (1 in the complete-case fit), so the stacked score is weighted once.
      corr <- apply_model_correction(prep, g)$correction

      # Channel 1: empirical-distribution sampling term on target rows.
      # Weighted form (n / sum_w) w_i (p_{k,i} - mu); with unit weights this is
      # (n / n_t)(p_{k,i} - mu), the complete-case Channel 1.
      ch1 <- numeric(n)
      ch1[target_idx] <- (n / sum_w) * w_target * (p_k - mu_mat[a, ci])

      IF_list[[a]] <- ch1 + corr
    }
    names(IF_list) <- int_names
    vcov_by_class[[ci]] <- vcov_from_if(IF_list, n, int_names)
  }

  vcov_by_class
}
