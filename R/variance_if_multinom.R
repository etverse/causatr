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
#' Under IPCW (missing Y), a third channel is added: the censoring model's
#' estimation uncertainty propagates into \eqn{\hat\mu} through the IPCW
#' weights via two paths, both carried through the censoring model's own
#' influence function. (1) \emph{Indirect} -- the stacked multinomial score
#' \eqn{s_i} depends on \eqn{\gamma} via the weight \eqn{w_i(\gamma)}, so the
#' bread cross-block \eqn{A_{\beta,\gamma} = -\partial\bar\Phi/\partial\gamma}
#' (with \eqn{\bar\Phi(\gamma) = n^{-1}\sum_i w_{\mathrm{ext},i}\,w_i(\gamma)\,
#' s_i^{(0)}} the weighted score at \eqn{\hat\beta}) projects each
#' class/intervention gradient onto the censoring-parameter space. (2)
#' \emph{Direct} -- the IPCW-weighted marginal average itself carries
#' \eqn{\gamma} in its weights, contributing
#' \eqn{\partial\mu_{k,a}/\partial\gamma = (\sum w)^{-1}\sum_i
#' (dw_i/d\gamma)(p_{k,i}(a) - \mu_{k,a})}. The total sensitivity
#' \eqn{d\mu/d\gamma} stacks both before the projection; the direct path
#' recovers the efficiency gain from estimating the censoring model (Robins,
#' Rotnitzky & Zhao 1994), cancels in contrasts, and matters for the per-class
#' means. With `ipcw = FALSE` this channel is skipped and the sandwich is
#' byte-identical to the complete-case / weighted form.
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
  n_fit <- nrow(prep$X_fit)

  # IPCW (missing-Y) censoring cross-term. The outcome multinom is fit on the
  # uncensored rows with the total weight (survey x IPCW), so Channel 1 / 2
  # above are already correct for the *fitted* nuisance; what is missing is the
  # third term -- the censoring model's own estimation uncertainty, which feeds
  # mu_hat through the IPCW weights. The bread cross-block A_{beta,gamma} and the
  # censoring-model prep do not depend on the outcome class or the intervention,
  # so build them once and project per (class, intervention) inside the loop.
  has_ipcw <- isTRUE(fit$details$ipcw)
  ipcw_cross <- if (has_ipcw) {
    multinom_ipcw_cross_setup(fit, prep, target_idx)
  } else {
    NULL
  }

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
      res <- apply_model_correction(prep, g)
      corr <- res$correction

      # Channel 3 (IPCW only): subtract the censoring-model estimation
      # cross-term. The IPCW-weighted marginal mean mu_{k,a}(gamma, beta)
      # depends on the censoring parameter gamma through two paths, so the
      # total sensitivity d mu / d gamma stacks both before projecting through
      # the censoring-model IF:
      #   * Indirect (gamma -> outcome beta -> mu): g_cens = A_{beta,gamma}^T h,
      #     with h = A_{bb}^{-1} g = n_fit * res$h in M-estimation scaling
      #     (A_{bb}^{-1} = n_fit * B_inv). This mirrors the scalar gcomp IPCW
      #     path (`compute_ipcw_if_correction()`), generalised to the stacked
      #     multinomial score.
      #   * Direct (gamma -> IPCW weights -> covariate average): the weighted
      #     mean carries gamma in its weights, so
      #       d mu_{k,a} / d gamma = (1/sum_w) sum_i (dw_i/dgamma)(p_{k,i}(a) - mu),
      #     with dw_i/dgamma = -w_i (1 - p_unc_i) X_cens_i. Omitting this term
      #     leaves the SE at its known-weights (conservative) value; including
      #     it recovers the efficiency gain from estimating the censoring model
      #     (Robins, Rotnitzky & Zhao 1994). It cancels in contrasts but
      #     matters for the per-class marginal means.
      if (has_ipcw) {
        g_cens <- as.numeric(crossprod(
          ipcw_cross$A_beta_gamma,
          n_fit * res$h
        ))
        # Direct path: the IPCW-weighted average's own gamma-dependence, summed
        # over the target rows that carry a positive weight (`ipcw_direct_grad()`
        # restricts to the target so ATT / by / subset stay aligned). The
        # per-class column `preds_list[[a]][, ci]` plays the role of p_{ij}.
        grad_mu_gamma <- ipcw_direct_grad(
          ipcw_cross$direct,
          preds_list[[a]][, ci],
          mu_mat[a, ci]
        )
        corr <- corr -
          apply_model_correction(
            ipcw_cross$cens_prep,
            g_cens - grad_mu_gamma
          )$correction
      }

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


#' Precompute the IPCW censoring cross-block for the multinomial sandwich
#'
#' @description
#' Builds the intervention- and class-independent ingredients of the
#' multinomial IPCW Channel-3 correction: the bread cross-block
#' \eqn{A_{\beta,\gamma}} (indirect path), the censoring-model
#' `prepare_model_if()` output, and the direct-path pieces (the censoring
#' design \eqn{X_{\mathrm{cens}}} and the weight derivative
#' \eqn{dw_i/d\gamma} on the outcome fit rows). All are reused across every
#' outcome class and every intervention, so they are computed once and the
#' per-class loop only assembles the total sensitivity
#' \eqn{d\mu/d\gamma = A_{\beta,\gamma}^\top h - \partial\mu/\partial\gamma}
#' and projects it through `apply_model_correction()`.
#'
#' @details
#' The fitted multinomial-logit score for observation `i` is
#' \eqn{s_i = w_{\mathrm{ext},i}\, w_i(\gamma)\, s_i^{(0)}}, where
#' \eqn{s_i^{(0)} = (I(Y_i = l) - p_{l,i})_l \otimes X_i} is the unweighted
#' stacked residual score (the rows of `prep$X_fit`) held fixed at
#' \eqn{\hat\beta}, \eqn{w_{\mathrm{ext},i}} is the survey weight that does not
#' depend on \eqn{\gamma} (`weights_pre_ipcw`, 1 when absent), and
#' \eqn{w_i(\gamma)} is the stabilized IPCW weight from `make_ipcw_weight_fn()`.
#' Only \eqn{w_i(\gamma)} varies with the censoring parameter, so
#' \deqn{\bar\Phi(\gamma) = \frac{1}{n_{\mathrm{fit}}}\sum_i
#'       w_{\mathrm{ext},i}\, w_i(\gamma)\, s_i^{(0)}, \qquad
#'       A_{\beta,\gamma} = -\,\partial\bar\Phi/\partial\gamma}
#' is a \eqn{((K-1)p) \times q} cross-block (numerically differentiated with
#' `numDeriv::jacobian()`). At \eqn{\hat\gamma} the per-row scale
#' \eqn{w_{\mathrm{ext},i} w_i(\hat\gamma)} equals the total fit weight
#' `prep$r_score`, so \eqn{\bar\Phi(\hat\gamma)} is the vanishing MLE score --
#' the closure is faithful to the fitted model by construction.
#'
#' @param fit A `causatr_fit` with `details$ipcw == TRUE` and a point censoring
#'   model in `details$censoring_model`.
#' @param prep Output of `prepare_model_if_multinom()` for the outcome model
#'   (carries the stacked score `X_fit`, the prior weight `r_score`, the fit
#'   indices `fit_idx`, and the total denominator `n_total`).
#' @param target_idx Logical vector (length `n`) flagging the target-population
#'   rows the marginal mean is averaged over (already NA-pruned by the caller).
#'   The direct path is summed over the target rows with a positive IPCW weight.
#'
#' @return A named list with `A_beta_gamma` (the \eqn{((K-1)p) \times q}
#'   indirect-path cross-block), `cens_prep` (the censoring-model
#'   `prepare_model_if()` output, reused for every projection), and `direct`
#'   (the shared `ipcw_direct_grad_setup()` output feeding `ipcw_direct_grad()`
#'   for the per-class direct gradient).
#'
#' @seealso `compute_ipcw_if_correction()` for the scalar-outcome analogue,
#'   `make_ipcw_weight_fn()` for the IPCW weight closure.
#' @family variance
#' @noRd
multinom_ipcw_cross_setup <- function(fit, prep, target_idx) {
  cens_model <- fit$details$censoring_model
  n_total <- prep$n_total
  fit_idx <- prep$fit_idx
  n_fit <- nrow(prep$X_fit)

  # Survey weight that does NOT vary with gamma. `weights_pre_ipcw` is the
  # external weight stashed before IPCW composition (NULL when there are no
  # survey weights), so the gamma-varying factor is exactly the IPCW weight.
  w_ext <- fit$details$weights_pre_ipcw
  w_ext_fit <- if (is.null(w_ext)) rep(1, n_fit) else w_ext[fit_idx]

  # Rows of the stacked design ARE the unweighted per-obs multinomial score
  # (residual tensor X_i), held fixed at beta_hat; the weight enters only as
  # the scalar `w_ext_i * w_ipcw_i(gamma)` per row.
  X_stacked <- prep$X_fit

  # IPCW weight closure: gamma -> full-length stabilized weight vector.
  ipcw_wfn <- make_ipcw_weight_fn(
    cens_model,
    n_total = n_total,
    censoring_col = as.integer(fit$data[[fit$censoring]]),
    stabilize = TRUE
  )

  # phi_bar(gamma): the weighted stacked score as a function of the censoring
  # parameter, holding beta fixed. crossprod(X_stacked, scale) sums
  # scale_i * X_stacked[i, ] over the fit rows, giving the ((K-1)*p) score.
  phi_bar_cens <- function(gamma) {
    w_ipcw_fit <- ipcw_wfn(gamma)[fit_idx]
    as.numeric(crossprod(X_stacked, w_ext_fit * w_ipcw_fit)) / n_fit
  }

  gamma_hat <- cens_model$alpha_hat
  # A_{beta,gamma} = -d phi_bar / d gamma (the M-estimation bread convention
  # A = -E[d psi / d theta]; Stefanski & Boos 2002).
  A_beta_gamma <- -numDeriv::jacobian(phi_bar_cens, x = gamma_hat)

  cens_prep <- prepare_model_if(
    cens_model$model,
    which(cens_model$fit_rows),
    n_total
  )

  # Direct-path ingredients (censoring design + dw_i/dgamma on the target rows
  # that carry a positive weight, plus the weight total) are shared with the
  # scalar gcomp IPCW path, so they come from the common helper. The per-class
  # gradient is then `ipcw_direct_grad()` in the loop.
  direct <- ipcw_direct_grad_setup(fit, fit_idx, target_idx)

  list(
    A_beta_gamma = A_beta_gamma,
    cens_prep = cens_prep,
    direct = direct
  )
}
