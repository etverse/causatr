#' IPW branch of variance_if()
#'
#' @description
#' Straight per-intervention loop over the density-ratio stacked
#' sandwich. For each intervention, builds the three pieces
#' `compute_ipw_if_self_contained_one()` needs -- Channel 1, the
#' marginal-mean Jacobian `J`, and the `weight_fn` closure -- from the
#' per-intervention bundle, then aggregates via `vcov_from_if()`.
#'
#' ## Scale convention
#'
#' Variance is computed **on the MSM fit-row subset** (`n = n_fit`).
#' `compute_ipw_if_self_contained_one()` asserts `n_fit == n_ps ==
#' n_total` internally, so we feed it the same `n_fit` for all three
#' and align Channel 1 / the target mask to the fit-row subset via
#' `target_idx[fit_rows]`. This matches the hand-derivation in
#' `test-ipw-branch-b.R`, where the analytic stacked sandwich is also
#' computed on the full `n = nrow(d)` (all rows had non-missing
#' outcomes).
#'
#' @param fit A `causatr_fit` of estimator `"ipw"`.
#' @param bundles Named list of per-intervention bundles built by
#'   `compute_ipw_contrast_point()`. Each carries `intervention`,
#'   `msm_model`, `weights_final`, `mu_hat`.
#' @param target_idx Logical vector (length `nrow(fit$data)`) flagging
#'   target-population rows.
#' @param mu_hat Named numeric vector of marginal-mean point estimates.
#' @param fit_idx_full Integer vector of MSM fit rows in `1..nrow(fit$data)`.
#' @param n_total Integer. `nrow(fit$data)`; used only for the
#'   final vcov scaling step (`vcov_from_if` divides by `n_sub^2`).
#'
#' @return A `k x k` variance-covariance matrix.
#'
#' @noRd
variance_if_ipw <- function(
  fit,
  bundles,
  target_idx,
  mu_hat,
  fit_idx_full,
  n_total,
  cluster_vec = NULL,
  trim = 1
) {
  int_names <- names(bundles)
  data <- fit$data
  is_mv <- isTRUE(fit$details$is_multivariate)
  tm <- fit$details$treatment_model
  tms <- fit$details$treatment_models
  propensity_model <- if (is_mv) NULL else tm$model
  ext_w <- fit$details$weights
  # Same estimand threading as `compute_ipw_contrast_point()` -- the
  # weight closure must match the weights the MSM was fit with, or the
  # cross-derivative `A_{beta, alpha}` is computed for the wrong weight
  # formula and the propensity correction drifts. For ATT / ATC the
  # closure's f*(p) term picks up the propensity-score dependence that
  # makes the sandwich SE agree with `WeightIt::glm_weightit`.
  estimand <- fit$estimand

  # Subset to the MSM fit rows. Everything below operates in the
  # length-`n_sub` space.
  fit_rows <- fit$details$fit_rows
  n_sub <- length(fit_idx_full)
  target_sub <- target_idx[fit_rows]
  # For transport, the MSM with sampling weights already targets the
  # target population. All study rows contribute to Ch1 / Jacobian;
  # target_idx selects S=0 rows which aren't in fit_rows.
  if (isTRUE(fit$details$transport)) {
    target_sub <- rep(TRUE, n_sub)
  }
  ext_w_sub <- if (is.null(ext_w)) NULL else ext_w[fit_rows]
  # Align the cluster vector to the MSM fit-row subset. `cluster_vec`
  # comes from `resolve_cluster()` at length `n_total = nrow(fit$data)`,
  # and IFs in this branch live on `fit_rows` (outcome-filtered subset),
  # so we subset in the same order the IFs are assembled.
  cluster_sub <- if (is.null(cluster_vec)) NULL else cluster_vec[fit_rows]

  # Target-population weights for the Channel-1 / Jacobian averaging.
  # `sum_w_target` plays the role of the denominator in a Hajek mean
  # over the target population; unweighted case uses a uniform weight
  # so the formula degenerates to `/n_target` -- matching the
  # `build_point_channel_pieces()` recipe.
  if (is.null(ext_w_sub)) {
    w_target_vec <- rep(1, n_sub)
    w_target_vec[!target_sub] <- 0
    sum_w_target <- sum(target_sub)
  } else {
    w_target_vec <- ext_w_sub
    w_target_vec[!target_sub] <- 0
    sum_w_target <- sum(ext_w_sub[target_sub])
  }
  if (sum_w_target <= 0) {
    rlang::abort(
      "variance_if_ipw(): target-population weights sum to 0.",
      class = "causatr_empty_target"
    )
  }

  IF_list <- lapply(int_names, function(nm) {
    b <- bundles[[nm]]
    msm_model <- b$msm_model
    mu_hat_j <- mu_hat[[nm]]

    # Counterfactual design matrix on the MSM fit rows under the
    # per-intervention rule. Without EM, the MSM is `Y ~ 1` (single
    # column of ones). With EM, the MSM is `Y ~ 1 + modifier` (extra
    # columns for modifier main effects). `iv_design_matrix()` handles
    # both. `apply_intervention()` runs on the full data; we subset to
    # the fit rows to keep lengths aligned with the MSM row count.
    #
    # IPSI does not materialize a counterfactual treatment value -- the
    # intervention acts on the propensity, not on A itself. The MSM
    # design matrix depends only on modifier columns (unchanged by
    # IPSI), so we skip `apply_intervention()` and use the original
    # data.
    iv_type_v <- if (inherits(b$intervention, "causatr_intervention")) {
      b$intervention$type
    } else {
      NULL
    }
    data_a_full <- if (identical(iv_type_v, "ipsi")) {
      data
    } else {
      apply_intervention(data, fit$treatment, b$intervention)
    }
    data_a_sub <- data_a_full[fit_rows]
    X_star <- iv_design_matrix(msm_model, data_a_sub)

    beta_hat <- coef_clean(msm_model)
    eta_star <- as.numeric(X_star %*% beta_hat)
    mu_eta_star <- msm_model$family$mu.eta(eta_star)
    preds_sub <- msm_model$family$linkinv(eta_star)

    # Channel 1: Ch1_i = n * (w_i / sum_w) * (\hat{mu}^g_i - \hat{mu}^g),
    # i.e. the Hajek residual scaled so vcov_from_if(n_sub) gives 1/n^2 * sum IF^2.
    # Off-target rows contribute zero; preds_sub is \hat{mu}^g_i from the MSM.
    Ch1_i <- n_sub * (w_target_vec / sum_w_target) * (preds_sub - mu_hat_j)
    Ch1_i[!target_sub] <- 0

    # Marginal-mean Jacobian J = \partial \hat{mu} / \partial \beta.
    # Hajek mean over target: J = (1/sum_w) sum_{i: target} w_i X*_i mu'(\eta*_i).
    # Chain rule on g(\hat\beta) = (1/sum_w) sum w_i g(X*_i \hat\beta).
    w_vec_j <- w_target_vec # zero off target already
    J <- as.numeric(crossprod(X_star, w_vec_j * mu_eta_star)) /
      sum_w_target

    # Per-intervention weight closure. `make_weight_fn()` captures
    # the treatment model and the intervention by value, then
    # returns a `function(alpha)` that `compute_ipw_if_self_contained_one()`
    # feeds into `numDeriv::jacobian()` for the cross-derivative
    # `A_{beta, alpha}`.
    # `tm$fit_rows` is relative to the outcome-filtered subset, so
    # the weight closure must receive `data[fit_rows]` (not the full
    # `data`). Without this, outcome NAs create a length mismatch.
    fit_data_local <- data[fit_rows]

    # Multivariate dispatch: build a stacked-alpha closure across the
    # K propensity models, then route to the multi-model variance
    # primitive that sums per-component propensity corrections.
    if (is_mv) {
      tms_local <- tms
      for (k in seq_along(tms_local)) {
        tms_local[[k]]$fit_rows <- rep(TRUE, nrow(fit_data_local))
      }
      class(tms_local) <- c("causatr_treatment_models", "list")

      mv_closure <- make_weight_fn_mv(
        treatment_models = tms_local,
        data = fit_data_local,
        interventions = b$intervention,
        estimand = "ATE",
        trim = trim
      )
      wfn <- mv_closure$weight_fn

      if (!is.null(ext_w_sub)) {
        ext_w_closure <- ext_w_sub
        base_wfn <- wfn
        wfn <- function(alpha) base_wfn(alpha) * ext_w_closure
      }

      if_i_mv <- compute_ipw_if_self_contained_mv_one(
        msm_model = msm_model,
        treatment_models = tms_local,
        weight_fn = wfn,
        alpha_offsets = mv_closure$offsets,
        alpha_hat_stacked = mv_closure$alpha_hat,
        J = J,
        Ch1_i = Ch1_i,
        fit_idx = seq_len(n_sub),
        n_total = n_sub
      )
      # Add the sampling-model correction block when transport is active.
      # The multivariate propensity correction is already handled inside
      # compute_ipw_if_self_contained_mv_one() (K per-component blocks).
      # The gamma block for the sampling model is independent of the
      # number of propensity components: same formula as the univariate
      # transport path, only the MSM structure changes.
      if (isTRUE(fit$details$transport)) {
        if_i_mv <- if_i_mv -
          compute_ipw_sampling_correction(
            fit,
            msm_model,
            J,
            fit_rows = fit_rows,
            n_sub = n_sub
          )
      }
      return(if_i_mv)
    }

    wfn <- make_weight_fn(
      tm,
      fit_data_local,
      b$intervention,
      estimand = estimand,
      trim = trim
    )

    # The base weight closure covers only the density-ratio piece.
    # The MSM score is psi_beta_i = (ext_w * DR_w) * X * r * mu_eta / var_mu,
    # so the cross-derivative A_{beta,alpha} = -(1/n) sum d psi_beta / d alpha
    # needs ext_w as a constant multiplier of the alpha-varying DR_w term.
    # Without ext_w in the closure, phi_bar() inside
    # compute_ipw_if_self_contained_one() would compute A_{beta,alpha}
    # on the wrong (unweighted) score, and the propensity correction
    # would be under-scaled by ext_w. The MSM bread A_{beta,beta}
    # already has ext_w baked in (via the MSM model's IWLS weights),
    # and A_{beta,alpha} must be computed on the same scale.
    if (!is.null(ext_w_sub)) {
      ext_w_closure <- ext_w_sub
      base_wfn <- wfn
      wfn <- function(alpha) base_wfn(alpha) * ext_w_closure
    }

    if_i <- compute_ipw_if_self_contained_one(
      msm_model = msm_model,
      propensity_model = propensity_model,
      weight_fn = wfn,
      J = J,
      Ch1_i = Ch1_i,
      fit_idx = seq_len(n_sub),
      fit_idx_ps = seq_len(n_sub),
      n_total = n_sub
    )

    if (isTRUE(fit$details$ipcw)) {
      if_i <- if_i -
        compute_ipw_ipcw_correction(
          fit,
          msm_model,
          J,
          fit_rows = fit_rows,
          n_sub = n_sub
        )
    }

    if (isTRUE(fit$details$transport)) {
      if_i <- if_i -
        compute_ipw_sampling_correction(
          fit,
          msm_model,
          J,
          fit_rows = fit_rows,
          n_sub = n_sub
        )
    }

    if_i
  })
  names(IF_list) <- int_names

  vcov_from_if(IF_list, n_sub, int_names, cluster = cluster_sub)
}


#' Per-individual IF for one IPW intervention
#'
#' @description
#' Returns the length-`n_total` per-individual influence function for
#' ONE intervention, assembled from three channels following
#' the IF decomposition in the variance theory vignette (Section 4):
#'
#' \deqn{\mathrm{IF}_i = \underbrace{\mathrm{Ch.1}_i}_{\text{direct}}
#'       + \underbrace{J^{\mathsf T}\,A_{\beta\beta}^{-1}\,\psi_{\beta,i}}_{\text{MSM correction}}
#'       + \underbrace{J^{\mathsf T}\,A_{\beta\beta}^{-1}\,A_{\beta\alpha}\,A_{\alpha\alpha}^{-1}\,\psi_{\alpha,i}}_{\text{propensity correction}}.}
#'
#' The interesting piece is the cross-derivative
#' \eqn{A_{\beta\alpha} = -\partial\bar\psi_\beta/\partial\alpha}.
#' We compute it numerically via `numDeriv::jacobian()` on a closure
#' that recomputes the average weighted MSM score as a function of
#' candidate propensity parameters \eqn{\alpha}. That closure is the
#' `weight_fn` built by `make_weight_fn()` for this specific
#' intervention, wrapped with the fixed-at-`beta_hat` residual /
#' linkinv machinery.
#'
#' ## Why this function uses `apply_model_correction()` twice
#'
#' The user's treatment density model is a plain `glm`, so using
#' `sandwich::bread.glm` here would land in the zero-weight-dispersion
#' sharp edge (`summary.glm$dispersion` computed on the effective-
#' weight subset, vs `sandwich::bread.glm`'s use of nominal degrees of
#' freedom). The two disagree by a ratio of effective-to-nominal n,
#' which is non-trivial on HT-weighted fits where half the rows have
#' `w = 0`. Rather than patch around that, we use causatr's own
#' `apply_model_correction()`, which computes the same quantity from
#' the GLM's working residuals + working weights directly -- no
#' `summary.glm` intermediary -- and handles zero-weight rows by
#' simply letting `r_score_i = (Y_i - mu_i) * 0 = 0` for those rows.
#' That matches the gcomp / ICE / matching convention exactly.
#'
#' Both channels go through the same primitive:
#'
#' - **MSM correction.** `apply_model_correction(msm_prep, J)`. The
#'   gradient `J` is the per-intervention Jacobian
#'   \eqn{\partial\hat\mu/\partial\beta}. Returns `$correction`
#'   (per-individual MSM-side IF contribution) and `$h = B_{\text{inv}} J
#'   = (X^{\mathsf T}WX)^{-1} J`.
#'
#' - **Propensity correction.** `apply_model_correction(prop_prep, g_prop)`.
#'   The gradient `g_prop` is the cross-derivative-projected MSM
#'   Jacobian, \eqn{A_{\beta\alpha}^{\mathsf T}\,h_{\text{msm}}}. Under
#'   causatr's convention `msm_res$h = A_{\beta\beta}^{-1} J / n_{\text{fit}}`
#'   (because `bread_inv` returns the raw `(X^{\mathsf T}WX)^{-1}`, not
#'   `n \cdot (X^{\mathsf T}WX)^{-1}`), so we multiply by `n_fit` once
#'   to recover the "true" \eqn{A_{\beta\beta}^{-1} J}.
#'
#' ## Sign of the propensity correction
#'
#' The full IF for the beta block of a stacked M-estimator with
#' block-lower-triangular bread
#' \eqn{A = \begin{pmatrix} A_{\alpha\alpha} & 0 \\ A_{\beta\alpha} & A_{\beta\beta} \end{pmatrix}}
#' is
#' \deqn{\mathrm{IF}_i(\beta) = A_{\beta\beta}^{-1}\,\bigl(\psi_{\beta,i}
#'       - A_{\beta\alpha}\,A_{\alpha\alpha}^{-1}\,\psi_{\alpha,i}\bigr),}
#' i.e. the propensity correction is **subtracted** from the MSM
#' term. This comes out of the block inversion
#' \eqn{A^{-1} = \begin{pmatrix} A_{\alpha\alpha}^{-1} & 0 \\ -A_{\beta\beta}^{-1} A_{\beta\alpha} A_{\alpha\alpha}^{-1} & A_{\beta\beta}^{-1} \end{pmatrix}}
#' applied to \eqn{\psi_i = (\psi_{\alpha,i}, \psi_{\beta,i})}. We
#' return `msm_res$correction - prop_res$correction` at the end.
#'
#' (The existing `variance-theory.qmd Sec.4.2` writes the IF with a
#' `+` sign because it uses Wooldridge's convention
#' \eqn{A_{\beta\alpha}^{\text{W}} = +(1/n)\sum \partial\psi_\beta/\partial\alpha},
#' whereas we use the negative-Hessian convention
#' \eqn{A_{\beta\alpha} = -(1/n)\sum \partial\psi_\beta/\partial\alpha}
#' that `numDeriv::jacobian(phi_bar)` naturally gives. The two
#' conventions differ by a sign on `A_{beta,alpha}`, which flips the
#' composition sign. The numerical result is the same; the code
#' must just be self-consistent.)
#'
#' ## Channel scaling
#'
#' - **Channel 1** comes from `build_point_channel_pieces()` and is
#'   \eqn{n_{\text{total}}}-scaled
#'   (`Ch1_list[[j]] = n_total * target_w * (p - mu_hat[j])`). For a
#'   saturated MSM with static intervention Ch1 is identically zero.
#' - **MSM correction** from `apply_model_correction()` is also
#'   \eqn{n_{\text{total}}}-scaled (`correction[fit_idx] <- n_total *
#'   d_fit * r_score`).
#' - **Propensity correction** same story.
#'
#' All three channels are thus in the same scaling convention. Under
#' the invariant `n_fit == n_ps == n_total` (enforced upstream by
#' `fit_ipw()`), `vcov_from_if(IF_list, n_total)` reproduces the
#' classical stacked sandwich variance.
#'
#' ## Scope (what this function does NOT do)
#'
#' - It does not fit any model. The MSM and the propensity model are
#'   passed in fully fitted.
#' - It does not build `Ch1_i`. That comes from
#'   `build_point_channel_pieces()` in the usual `variance_if_*` flow.
#' - It does not aggregate IFs into a vcov. The caller hands the
#'   returned vector to `vcov_from_if()` alongside every other
#'   intervention's IF.
#' - It does not iterate over interventions. One call per intervention.
#'
#' The caller owns the per-intervention loop.
#'
#' @param msm_model Fitted weighted GLM for this intervention. Must
#'   support `stats::coef()`, `stats::model.matrix()`,
#'   `stats::model.response(stats::model.frame(.))`, and the GLM
#'   working-residual / working-weight accessors used by
#'   `prepare_model_if()`.
#' @param propensity_model Fitted treatment density model (typically
#'   `tm$model` from `fit_treatment_model()`). Same GLM accessor
#'   contract as `msm_model`.
#' @param weight_fn Closure `alpha -> w_i(alpha)` built by
#'   `make_weight_fn()` for this intervention. The closure must
#'   return a length-`n_fit` vector for any candidate `alpha` of the
#'   same length as `stats::coef(propensity_model)`.
#' @param J Numeric `p_beta`-vector. The per-intervention Jacobian
#'   \eqn{J = \partial\hat\mu_a/\partial\beta}. For a saturated MSM
#'   with static binary treatment, `J` is the counterfactual design
#'   row averaged over the target population. For an intercept-only
#'   MSM under a non-saturated intervention (shift, MTP, IPSI), `J =
#'   dmu/dbeta_0` which for an identity-link Gaussian is `1` and for
#'   a logit-link binomial is `mu_hat * (1 - mu_hat)`.
#' @param Ch1_i Numeric vector of length `n_total`. The Channel 1
#'   contribution from `build_point_channel_pieces()`. Zero on rows
#'   outside the target population and (for saturated static MSMs
#'   with predictions constant per treatment level) zero everywhere.
#' @param fit_idx Integer vector. Indices of MSM fit rows in
#'   `1..n_total`. Typically `seq_len(n_total)` under the
#'   "same row set for both models" invariant.
#' @param fit_idx_ps Integer vector. Indices of propensity fit rows
#'   in `1..n_total`. Same invariant as `fit_idx`.
#' @param n_total Integer. Length of the returned IF vector.
#'
#' @return Numeric vector of length `n_total`, the per-individual
#'   influence function for one intervention under the IPW engine.
#'
#' @noRd
compute_ipw_if_self_contained_one <- function(
  msm_model,
  propensity_model,
  weight_fn,
  J,
  Ch1_i,
  fit_idx,
  fit_idx_ps,
  n_total
) {
  # ---- MSM correction via causatr's primitive --------------------
  # `apply_model_correction()` returns:
  #   $correction: per-individual MSM-side contribution
  #                (n_total-scaled, zero off `fit_idx`)
  #   $h:          bread_inv %*% J = (X'WX)^{-1} J
  # The "true" A_bb^{-1} J is n_fit * $h because `bread_inv` returns
  # the raw (X'WX)^{-1} without the factor of n that the M-estimation
  # definition of A_bb carries.
  msm_prep <- prepare_model_if(msm_model, fit_idx, n_total)
  msm_res <- apply_model_correction(msm_prep, J)
  n_fit <- nrow(msm_prep$X_fit)

  # ---- Cross-derivative A_{beta, alpha} via numDeriv -------------
  # A_{\beta\alpha} = -(1/n) sum_i \partial\psi_\beta_i / \partial\alpha.
  # psi_\beta_i = X_i w_i(\alpha) (Y_i - \mu_i) mu'_i / V(\mu_i);
  # \mu_i, mu'_i, V(\mu_i), Y_i - \mu_i all fixed at \hat\beta.
  # Only w_i(\alpha) varies, so numDeriv re-evaluates only the weight formula.
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
    # \bar\psi_\beta(\alpha) = (1/n) X' diag(w(alpha) * mu'(eta) * r / V(\mu)) 1
    # Only w(alpha) changes with alpha; beta-dependent terms are pre-computed.
    w_alpha <- weight_fn(alpha)
    s_per_i <- w_alpha * mu_eta * r_fit / var_mu
    as.numeric(crossprod(X_msm, s_per_i)) / n_fit
  }
  # For multinomial models `coef()` returns a matrix; flatten to match
  # `make_weight_fn()`'s convention (row-major: `as.vector(t(coef_mat))`).
  # Drop aliased (NA) coefficients — they correspond to collinear columns
  # already removed from `X_prop` at fit time (see `fit_treatment_model()`).
  alpha_hat_raw <- stats::coef(propensity_model)
  if (!is.null(dim(alpha_hat_raw))) {
    alpha_hat <- as.vector(t(alpha_hat_raw))
  } else {
    alpha_hat <- alpha_hat_raw[!is.na(alpha_hat_raw)]
  }
  # Negative-Hessian convention: A_{beta, alpha} = -(1/n) sum d psi/d alpha.
  # numDeriv::jacobian(phi_bar, alpha) = d phi_bar/d alpha = +(1/n) sum d psi/d alpha.
  # Flipping the sign gives A_{beta, alpha}.
  A_beta_alpha <- -numDeriv::jacobian(phi_bar, x = alpha_hat)

  # ---- Propensity correction via causatr's primitive --------------
  # h_msm_true = A_{\beta\beta}^{-1} J. msm_res$h = (X'WX)^{-1} J,
  # which equals A_{\beta\beta}^{-1} J / n_fit under the M-estimation
  # definition A_{\beta\beta} = (1/n) X'WX. Multiply by n_fit to recover.
  h_msm_true <- n_fit * msm_res$h
  # g_prop = A_{\beta\alpha}^T h_msm_true (p_alpha-vector). Passed as
  # the sensitivity gradient to apply_model_correction on the propensity
  # model, which returns sum_i g_prop^T A_{\alpha\alpha}^{-1} \psi_{\alpha,i}
  # -- the per-individual propensity correction.
  g_prop <- as.numeric(crossprod(A_beta_alpha, h_msm_true))

  # Route to the multinomial-specific prep when the propensity model
  # is from `nnet::multinom` (which is not a GLM and lacks `$family`,
  # working weights, and the other internals that `prepare_model_if()`
  # relies on). The multinomial prep computes bread and estfun from the
  # multinomial log-likelihood score equations directly.
  if (inherits(propensity_model, "multinom")) {
    prop_prep <- prepare_model_if_multinom(
      propensity_model,
      fit_idx_ps,
      n_total
    )
  } else {
    prop_prep <- prepare_model_if(propensity_model, fit_idx_ps, n_total)
  }
  prop_res <- apply_model_correction(prop_prep, g_prop)

  # ---- Assemble ---------------------------------------------------
  # Invariant check: every channel must be in n_total scaling. That
  # holds as long as `n_fit == n_ps == n_total`, which is the
  # "same row set for both models" invariant. Violating it would
  # silently mis-scale the cross-model composition.
  if (nrow(X_msm) != n_total) {
    rlang::abort(
      paste0(
        "compute_ipw_if_self_contained_one(): n_fit (",
        nrow(X_msm),
        ") != n_total (",
        n_total,
        "). The IPW sandwich engine assumes the MSM fits on the same ",
        "row set as the full data. Drop NA rows in `causat()` before ",
        "the IPW path builds the MSM."
      ),
      class = "causatr_variance_row_mismatch"
    )
  }
  if (nrow(prop_prep$X_fit) != n_total) {
    rlang::abort(
      paste0(
        "compute_ipw_if_self_contained_one(): n_ps (",
        nrow(prop_prep$X_fit),
        ") != n_total (",
        n_total,
        "). Same row-alignment invariant as the MSM above."
      ),
      class = "causatr_variance_row_mismatch"
    )
  }

  # Block-lower-triangular M-estimation result:
  #   IF_beta_i = A_bb^{-1}(psi_beta_i - A_{beta, alpha} A_aa^{-1} psi_alpha_i)
  # i.e. the propensity correction is SUBTRACTED. See the docstring
  # note on the sign convention.
  Ch1_i + msm_res$correction - prop_res$correction
}


#' Per-individual IF for one multivariate IPW intervention
#'
#' @description
#' Multivariate companion to `compute_ipw_if_self_contained_one()`.
#' For a joint treatment `A = (A_1, ..., A_K)` factorised as
#' \eqn{f(A_1, \ldots, A_K \mid L) = \prod_k f(A_k \mid A_{1..k-1}, L)},
#' the K propensity models are fit independently (sequentially) on the
#' same row set, so the bread of the stacked propensity system is
#' **block-diagonal** in the per-component alpha blocks. The IF for
#' the beta block of the full M-estimation system therefore decomposes
#' as
#' \deqn{\mathrm{IF}_i(\beta) = A_{\beta\beta}^{-1}\bigl(\psi_{\beta,i}
#'       - \sum_{k=1}^K A_{\beta\alpha_k}\,A_{\alpha_k\alpha_k}^{-1}\,\psi_{\alpha_k,i}\bigr).}
#'
#' Implementation:
#'
#' 1. Compute the MSM correction once (`apply_model_correction(msm_prep, J)`).
#' 2. Compute the full cross-derivative
#'    \eqn{[A_{\beta\alpha_1}, \ldots, A_{\beta\alpha_K}]} via
#'    `numDeriv::jacobian(phi_bar, alpha_hat_stacked)` on the stacked
#'    weight closure. The closure splits `alpha` into K blocks and
#'    multiplies the K per-component sub-closures (built by
#'    `make_weight_fn_mv()`).
#' 3. Slice the cross-derivative into K column-blocks and call
#'    `apply_model_correction(prop_prep_k, g_prop_k)` for each
#'    component, where `g_prop_k = A_{\beta\alpha_k}^T h_msm_true`.
#' 4. Return `Ch1_i + msm_res$correction - sum_k prop_res_k$correction`.
#'
#' Block-diagonal bread relies on the K propensity models being fit
#' on disjoint parameter blocks. `fit_treatment_models()` fits each
#' component's univariate density independently via
#' `fit_treatment_model()`, so this invariant holds by construction.
#' If a future enhancement introduces shared structure across
#' components (e.g. joint MLE on a multivariate normal), this
#' primitive would need to be replaced by a true stacked-bread
#' computation including the off-diagonal `A_{\alpha_j\alpha_k}` blocks.
#'
#' @param msm_model Fitted weighted MSM for this intervention.
#' @param treatment_models A `causatr_treatment_models` list, K models.
#' @param weight_fn Stacked-alpha closure from `make_weight_fn_mv()`.
#' @param alpha_offsets Integer (K+1)-vector of block boundaries.
#' @param alpha_hat_stacked Numeric vector of stacked propensity
#'   coefficients (concatenation of per-component `alpha_hat`).
#' @param J Numeric `p_beta`-vector. Marginal-mean Jacobian as in the
#'   univariate primitive.
#' @param Ch1_i Numeric vector of length `n_total`. Channel 1
#'   contribution.
#' @param fit_idx Integer vector. Indices of MSM fit rows.
#' @param n_total Integer. Length of returned IF vector.
#'
#' @return Numeric vector of length `n_total`.
#'
#' @noRd
compute_ipw_if_self_contained_mv_one <- function(
  msm_model,
  treatment_models,
  weight_fn,
  alpha_offsets,
  alpha_hat_stacked,
  J,
  Ch1_i,
  fit_idx,
  n_total
) {
  # ---- MSM correction --------------------------------------------
  msm_prep <- prepare_model_if(msm_model, fit_idx, n_total)
  msm_res <- apply_model_correction(msm_prep, J)
  n_fit <- nrow(msm_prep$X_fit)

  # Natural course (NULL) for all K components: weight function is

  # constant (all ones), so no propensity parameters exist and the
  # cross-derivative A_{beta,alpha} is vacuous. Skip the jacobian
  # and return only the MSM correction + Ch1 (sampling channel).
  if (length(alpha_hat_stacked) == 0L) {
    return(Ch1_i + msm_res$correction)
  }

  # ---- Stacked cross-derivative via numDeriv ---------------------
  # phi_bar(alpha) = (1/n_fit) sum_i psi_beta_i(alpha, beta_hat).
  # Same shape as the univariate case; only `weight_fn` is the
  # stacked product closure here.
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
  # Negative-Hessian convention, same as the univariate primitive.
  A_beta_alpha <- -numDeriv::jacobian(phi_bar, x = alpha_hat_stacked)

  # `h_msm_true = A_bb^{-1} J = n_fit * msm_res$h` per the
  # `apply_model_correction()` scaling convention.
  h_msm_true <- n_fit * msm_res$h

  # ---- Per-component propensity correction -----------------------
  # Block-diagonal bread: each component's IF correction depends only
  # on its own alpha block, so the propensity correction is the SUM
  # over components of `apply_model_correction(prop_prep_k, g_prop_k)`.
  K <- length(treatment_models)
  total_prop_correction <- rep(0, n_total)

  for (k in seq_len(K)) {
    idx <- alpha_offsets[k]:(alpha_offsets[k + 1L] - 1L)
    # Slice the k-th column-block of A_{beta, alpha}; project onto
    # the MSM bread to get the gradient for the k-th propensity bread.
    A_block_k <- A_beta_alpha[, idx, drop = FALSE]
    g_prop_k <- as.numeric(crossprod(A_block_k, h_msm_true))

    # Dispatch on multinomial vs GLM (same convention as the
    # univariate `compute_ipw_if_self_contained_one()`). Multinomial
    # propensity models are not GLMs and need their own
    # bread/score primitive.
    prop_model_k <- treatment_models[[k]]$model
    prop_prep_k <- if (inherits(prop_model_k, "multinom")) {
      prepare_model_if_multinom(prop_model_k, fit_idx, n_total)
    } else {
      prepare_model_if(prop_model_k, fit_idx, n_total)
    }
    prop_res_k <- apply_model_correction(prop_prep_k, g_prop_k)
    total_prop_correction <- total_prop_correction + prop_res_k$correction
  }

  if (nrow(X_msm) != n_total) {
    rlang::abort(
      paste0(
        "compute_ipw_if_self_contained_mv_one(): n_fit (",
        nrow(X_msm),
        ") != n_total (",
        n_total,
        "). Multivariate IPW assumes the MSM and every propensity model fit on the same rows."
      ),
      class = "causatr_variance_row_mismatch"
    )
  }

  Ch1_i + msm_res$correction - total_prop_correction
}
