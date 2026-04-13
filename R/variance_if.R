#' Influence-function variance engine
#'
#' @description
#' This file contains the unified variance engine for causatr. Every causal
#' estimand \eqn{\hat\mu} is defined by a system of estimating equations
#' (model scores + mean equation), and its variance is obtained by computing
#' a per-individual influence function (IF) and summing the squares:
#'
#' \deqn{\widehat{\mathrm{Var}}(\hat\mu) = \frac{1}{n^2}\sum_i \mathrm{IF}_i^2}
#'
#' The IF decomposes into two channels (vignette `variance-theory.qmd`,
#' Sections 3 and 4):
#'
#' - **Channel 1** (universal): \eqn{(n/n_t) t_i (\mathrm{pred}_i - \hat\mu)},
#'   the direct contribution of individual \eqn{i}'s prediction.
#' - **Channel 2** (one term per nuisance model \eqn{k}):
#'   \eqn{J_k A_k^{-1}\psi_{k,i}}, the propagation of model parameter
#'   uncertainty into \eqn{\hat\mu}.
#'
#' This file holds:
#'
#' - The Channel-2 primitives — `bread_inv()`, `iv_design_matrix()`,
#'   `prepare_model_if()` / `apply_model_correction()` (the prep/apply
#'   split that lets g-comp, matching, and IPW pay the `p x p` bread solve
#'   once per model and reuse it across interventions), plus the thin
#'   single-gradient wrapper `correct_model()` that ICE uses.
#' - The `resolve_fit_idx()` helper (maps `fit$details$fit_rows` to
#'   model-local row indices, asserting the `na.action` invariant).
#' - The `vcov_from_if()` aggregator and the `variance_if_numeric()`
#'   two-tier numerical fallback.
#' - The dispatcher `variance_if()` with branches for g-computation,
#'   matching (cluster-robust), ICE (chained forward sensitivity
#'   recursion), and IPW.
#' - The IPW-specific `prepare_propensity_if()` /
#'   `apply_propensity_correction()` pair (Branch A: WeightIt shortcut;
#'   Branch B: self-contained, deferred to Phase 4).
#'
#' @name variance_if
#' @keywords internal
NULL


#' Inverse of a model's bread matrix
#'
#' @description
#' For standard GLMs returns \eqn{(X^T W X)^{-1}} where \eqn{W} are the
#' IWLS working weights \eqn{(d\mu/d\eta)^2 / V(\mu)}. For GAMs returns the
#' penalized Bayesian covariance `model$Vp = (X'WX + \lambda S)^{-1}`.
#'
#' This is the building block for `correct_model()` and the `IF_\beta`
#' construction inside `prepare_propensity_if()` Branch B.
#'
#' @param model A fitted model with a `family` object (GLM or GAM).
#' @param X_fit Design matrix from `model.matrix(model)`.
#'
#' @return A `p x p` matrix.
#'
#' @noRd
bread_inv <- function(model, X_fit) {
  if (inherits(model, "gam") && !is.null(model$Vp)) {
    return(model$Vp)
  }

  # Use model$weights when available — these are the IWLS working weights
  # from the converged fit, *including* prior weights, and match exactly
  # what sandwich::bread() and summary(glm)$cov.unscaled use. Recomputing
  # them from family$mu.eta and family$variance is consistent in theory
  # but loses 5–6 digits of precision in practice (different evaluation
  # paths through the family functions vs. the IWLS internals).
  w_iwls <- model$weights
  if (is.null(w_iwls)) {
    eta <- model$linear.predictors
    mu_eta <- model$family$mu.eta(eta)
    var_mu <- model$family$variance(stats::fitted(model))
    w_iwls <- mu_eta^2 / var_mu
    prior_w <- model$prior.weights
    if (!is.null(prior_w)) {
      w_iwls <- w_iwls * prior_w
    }
  }

  XtWX <- crossprod(X_fit, X_fit * w_iwls)
  tryCatch(
    solve(XtWX),
    error = function(e) {
      # Rate-limited: bootstrap loops can hit this error once per
      # replicate — without throttling the console fills up. The
      # underlying rank deficiency is the same each time.
      rlang::warn(
        c(
          "`X'WX` is singular; using `MASS::ginv()` as a fallback.",
          i = paste0(
            "This usually means the outcome model is rank-deficient ",
            "(collinear covariates, aliased factor levels, or a saturated ",
            "design with too few observations per cell). The pseudo-inverse ",
            "gives a minimum-norm solution but the resulting sandwich ",
            "variance may be miscalibrated. Inspect `summary(fit$model)` ",
            "for NA coefficients."
          )
        ),
        .frequency = "once",
        .frequency_id = "causatr_bread_inv_singular"
      )
      MASS::ginv(XtWX)
    }
  )
}


#' Design matrix for counterfactual (intervention) data
#'
#' @description
#' Builds the design matrix \eqn{X^*} for counterfactual data, handling
#' the GLM / GAM split:
#'
#' - **GAMs**: `predict(model, newdata, type = "lpmatrix")` returns the
#'   smooth basis matrix evaluated at the counterfactual covariates.
#' - **GLMs and other models**: `model.matrix()` on the model's terms,
#'   reusing `model$xlevels` so factor levels survive interventions
#'   that produce single-level treatment columns (e.g. `static("low")`).
#'
#' `newdata` is coerced to `data.frame` up front to avoid
#' `model.matrix()` failing on `data.table` inputs — `model.matrix()`
#' sometimes chokes on `data.table` because of how it dispatches the
#' `[` subset used internally when evaluating terms. The coercion is
#' cheap (shared column pointers) and side-steps a class of
#' environment-related edge cases that previously required callers to
#' wrap the call in a `tryCatch`.
#'
#' @param model A fitted model object.
#' @param newdata Data frame or data.table of counterfactual observations.
#'
#' @return A design matrix (`n_new x p`).
#'
#' @noRd
iv_design_matrix <- function(model, newdata) {
  # `predict(gam, ..., type = "lpmatrix")` accepts data.table directly.
  if (inherits(model, "gam")) {
    return(stats::predict(model, newdata = newdata, type = "lpmatrix"))
  }
  # GLM path: force data.frame so `model.matrix()` doesn't trip on
  # data.table's `[` dispatch. Using `delete.response()` on the model's
  # terms drops the LHS (Y) so we don't need the response column in
  # `newdata` — important because the counterfactual frame is built
  # from covariates alone.
  if (data.table::is.data.table(newdata)) {
    newdata <- as.data.frame(newdata)
  }
  pred_terms <- stats::delete.response(stats::terms(model))
  xlev <- model$xlevels
  stats::model.matrix(pred_terms, data = newdata, xlev = xlev)
}


#' Resolve a model's fit rows relative to the full data
#'
#' @description
#' Returns the integer row indices (in `1..n_total`) of the rows `model`
#' was actually fit on, starting from `fit$details$fit_rows` and removing
#' `model$na.action` rows dropped during fitting.
#'
#' `model$na.action` carries indices that are local to the subset passed
#' to `glm()` (not the full `fit$data`), so `fit_idx[-na_action]` is the
#' correct removal as long as the pipeline upstream passes pre-subsetted
#' data. This helper wraps the calculation and asserts the invariant so a
#' future regression in the upstream `fit_*` code paths surfaces loudly
#' rather than silently corrupting the IF alignment.
#'
#' @param fit A `causatr_fit` with a valid `$details$fit_rows` logical.
#' @param model The fitted outcome model.
#'
#' @return Integer vector of row indices in `1..nrow(fit$data)`.
#'
#' @noRd
resolve_fit_idx <- function(fit, model) {
  fit_idx <- which(fit$details$fit_rows)
  na_action <- model$na.action
  if (is.null(na_action)) {
    return(fit_idx)
  }
  if (max(na_action, 0L) > length(fit_idx)) {
    rlang::abort(
      paste0(
        "resolve_fit_idx(): `model$na.action` max index (",
        max(na_action),
        ") exceeds `sum(fit$details$fit_rows)` (",
        length(fit_idx),
        "); upstream fit pipeline is not pre-subsetting data."
      )
    )
  }
  fit_idx[-na_action]
}


#' Precompute the gradient-independent ingredients of `correct_model()`
#'
#' @description
#' `correct_model()` needs three things that depend only on the model and
#' the fit data, not on the intervention-specific gradient: the design
#' matrix `X_fit`, the inverse bread `B_inv`, and the GLM score residual
#' `r_score`. Callers that apply the correction for many interventions
#' (g-comp, IPW, matching) should call `prepare_model_if()` **once** per
#' model and then `apply_model_correction(prep, gradient)` per intervention,
#' avoiding `O(k)` recomputation of a `p \times p` `solve()`.
#'
#' The score residual
#' \eqn{r^{\mathrm{score}}_i = (Y_i - \hat\mu_i)(d\mu/d\eta)/V(\mu)} is
#' obtained from the GLM internals as `residuals(m, "working") *
#' weights(m, "working")` so it matches `sandwich::estfun()` row-for-row.
#' For canonical links it collapses to the response residual; for
#' non-canonical links (probit, cloglog, Gamma-log, ...) the link-scale
#' factor is essential or the IF is miscalibrated.
#'
#' @param model A fitted model with a `family` object (GLM or GAM).
#' @param fit_idx Integer vector. Row indices in `1..n_total` corresponding
#'   to the rows the model was fit on.
#' @param n_total Integer. The total denominator used to scale the IF.
#'
#' @return A list with components `model`, `X_fit`, `B_inv`, `r_score`,
#'   `fit_idx`, `n_total`, suitable for passing to
#'   `apply_model_correction()`.
#'
#' @noRd
prepare_model_if <- function(model, fit_idx, n_total) {
  X_fit <- stats::model.matrix(model)
  B_inv <- bread_inv(model, X_fit)

  r_score <- tryCatch(
    stats::residuals(model, type = "working") *
      stats::weights(model, type = "working"),
    error = function(e) {
      eta <- model$linear.predictors
      mu_eta <- model$family$mu.eta(eta)
      var_mu <- model$family$variance(stats::fitted(model))
      r <- stats::residuals(model, type = "response") * mu_eta / var_mu
      pw <- model$prior.weights
      if (!is.null(pw)) {
        r <- r * pw
      }
      r
    }
  )

  list(
    model = model,
    X_fit = X_fit,
    B_inv = B_inv,
    r_score = r_score,
    fit_idx = fit_idx,
    n_total = n_total
  )
}


#' Apply a prepared model correction to a single gradient
#'
#' @description
#' Gradient-specific half of `correct_model()`. Consumes a `prep` list
#' from `prepare_model_if()` and a sensitivity gradient
#' \eqn{g = \partial\hat\mu/\partial\beta}, and returns the per-individual
#' correction vector along with `d` and `h` (needed by ICE and by
#' `prepare_propensity_if()` Branch B).
#'
#' Following the three-line template of vignette Section 4:
#'
#' \deqn{h = A^{-1} g, \qquad d_i = X_i^T h, \qquad
#'       \mathrm{correction}_i = n \cdot d_i \cdot r^{\mathrm{score}}_i}
#'
#' The factor of \eqn{n} arises because \eqn{A^{-1} = n(X^TWX)^{-1}}
#' (see vignette Section 5.5).
#'
#' @param prep Output of `prepare_model_if()`.
#' @param gradient Numeric `p`-vector. The sensitivity gradient.
#'
#' @return A list with `correction`, `d`, `h` as in `correct_model()`.
#'
#' @noRd
apply_model_correction <- function(prep, gradient) {
  h <- as.numeric(prep$B_inv %*% gradient)
  d_fit <- as.numeric(prep$X_fit %*% h)

  n_total <- prep$n_total
  fit_idx <- prep$fit_idx

  d_full <- rep(0, n_total)
  d_full[fit_idx] <- d_fit

  correction <- rep(0, n_total)
  correction[fit_idx] <- n_total * d_fit * prep$r_score

  list(correction = correction, d = d_full, h = h)
}


#' Channel 2 correction for one model
#'
#' @description
#' Convenience wrapper around `prepare_model_if()` + `apply_model_correction()`
#' for callers that only need the correction for a single gradient. ICE uses
#' this entry point because its outer loop cycles over distinct models, one
#' per time step; g-comp / IPW / matching cycle over interventions at a
#' fixed model and should call the two halves directly.
#'
#' @inheritParams prepare_model_if
#' @param gradient Numeric `p`-vector. The sensitivity gradient
#'   \eqn{g = \partial\hat\mu/\partial\beta}.
#'
#' @return A list with components:
#'   \describe{
#'     \item{`correction`}{Numeric vector of length `n_total`. Zero for
#'       rows outside `fit_idx`. `n_total * d_i * r_score_i` on fit rows.}
#'     \item{`d`}{Numeric vector of length `n_total`. The per-individual
#'       sensitivity \eqn{d_i = X_i^T h}; zero off `fit_idx`.}
#'     \item{`h`}{Numeric `p`-vector. The bread-projected gradient
#'       \eqn{h = A^{-1} g}.}
#'   }
#'
#' @noRd
correct_model <- function(model, gradient, fit_idx, n_total) {
  prep <- prepare_model_if(model, fit_idx, n_total)
  apply_model_correction(prep, gradient)
}


#' Aggregate per-individual influence functions into a vcov matrix
#'
#' @description
#' Computes the \eqn{k \times k} variance-covariance matrix
#' \deqn{\widehat{\mathrm{Cov}}(\hat\mu_a, \hat\mu_b)
#'       = \frac{1}{n^2}\sum_i \mathrm{IF}_{a,i}\,\mathrm{IF}_{b,i}}
#' from a list of length-`n` IF vectors. With `cluster = NULL` (the
#' standard case for g-comp, IPW, and ICE) IFs are squared individually.
#' With a non-`NULL` cluster vector (matching), IFs are first **summed
#' within each cluster** and then squared, implementing the cluster-robust
#' aggregation from vignette Section 4.3:
#'
#' \deqn{\widehat{\mathrm{Var}}(\hat\mu_a)
#'       = \frac{1}{n^2}\sum_{c=1}^{C}\Bigl(\sum_{i \in c}\mathrm{IF}_{a,i}\Bigr)^2}
#'
#' Singletons (clusters of size 1) reduce trivially to the standard
#' formula, so unmatched rows mixed in with matched ones are handled
#' correctly without special-casing.
#'
#' @param IF_list Named list of length-`n` numeric vectors, one per
#'   intervention.
#' @param n Integer. The denominator for the \eqn{1/n^2} scaling. For
#'   matching this is the matched-sample size, not the original `nrow(data)`.
#' @param int_names Character vector of intervention names (length `k`).
#' @param cluster Optional vector of length `n` identifying cluster
#'   membership (e.g. matched subclass ids). `NULL` for independent
#'   aggregation.
#'
#' @return A named `k x k` variance-covariance matrix.
#'
#' @noRd
vcov_from_if <- function(IF_list, n, int_names, cluster = NULL) {
  IF_mat <- do.call(cbind, IF_list)

  if (!is.null(cluster)) {
    if (length(cluster) != nrow(IF_mat)) {
      rlang::abort(
        "`cluster` length must match the IF vector length in `vcov_from_if()`."
      )
    }
    IF_mat <- rowsum(IF_mat, as.factor(cluster), reorder = FALSE)
  }

  vcov_mat <- crossprod(IF_mat) / n^2
  rownames(vcov_mat) <- int_names
  colnames(vcov_mat) <- int_names
  vcov_mat
}


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
  weights = NULL
) {
  model <- fit$model
  int_names <- names(data_a_list)
  k <- length(int_names)
  n_total <- length(target_idx)
  beta_hat <- stats::coef(model)

  data_a_frames <- lapply(data_a_list, function(da) {
    as.data.frame(da)[target_idx, , drop = FALSE]
  })

  target_w <- numeric(n_total)
  if (!is.null(weights)) {
    sum_w <- sum(weights)
    target_w[target_idx] <- weights / sum_w
  } else {
    n_t <- sum(target_idx)
    target_w[target_idx] <- 1 / n_t
  }

  # Channel 1 vectors per intervention (length n_total). The
  # (n_total / n_t) scaling is folded into target_w * n_total below.
  Ch1_list <- lapply(seq_len(k), function(j) {
    p <- preds_list[[j]]
    n_total * target_w * (p - mu_hat[j])
  })

  # Marginal-mean Jacobian J (k x p), via numDeriv on predict().
  pred_fun <- function(beta) {
    model_tmp <- model
    model_tmp$coefficients <- beta
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
  J <- numDeriv::jacobian(pred_fun, x = beta_hat)
  if (!is.matrix(J)) {
    J <- matrix(J, nrow = k)
  }

  # Tier 1: sandwich::estfun() available -> recover full IF.
  estfun_ok <- tryCatch(
    {
      ef <- sandwich::estfun(model)
      is.matrix(ef) && nrow(ef) > 0L
    },
    error = function(e) FALSE
  )

  if (estfun_ok) {
    psi <- sandwich::estfun(model)
    A_inv <- sandwich::bread(model) / stats::nobs(model)
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
      # Batched Channel-2: one matrix multiply gives the (n_fit x k)
      # per-observation Ch2 contributions across all interventions.
      Ch2_fit <- (IF_beta %*% t(J)) * n_total
      IF_list <- lapply(seq_len(k), function(j) {
        IF <- Ch1_list[[j]]
        IF[fit_idx] <- IF[fit_idx] + Ch2_fit[, j]
        IF
      })
      names(IF_list) <- int_names
      return(vcov_from_if(IF_list, n_total, int_names))
    }
  }

  # Tier 2: drop the cross-term, warn.
  rlang::warn(
    paste0(
      "Model class '",
      class(model)[1],
      "' exposes neither analytic bread (family$mu.eta / $Vp) nor a ",
      "sandwich::estfun() method. Falling back to V1 + J V_beta J^T, ",
      "which drops the IF cross-term (vignette Section 3.3) and may ",
      "slightly miscalibrate the variance. Use ci_method = 'bootstrap' ",
      "for an exact answer."
    )
  )

  Ch1_mat <- do.call(cbind, Ch1_list)
  V1 <- crossprod(Ch1_mat) / n_total^2
  V_beta <- tryCatch(stats::vcov(model), error = function(e) NULL)
  if (is.null(V_beta)) {
    V_beta <- diag(0, length(beta_hat))
  }
  V2 <- J %*% V_beta %*% t(J)

  vcov_mat <- V2 + V1
  rownames(vcov_mat) <- int_names
  colnames(vcov_mat) <- int_names
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
#'   Channel 2 (the combined MSM correction + propensity correction)
#'   is delegated to `prepare_propensity_if()` / `apply_propensity_correction()`,
#'   which dispatch between Branch A (WeightIt shortcut via
#'   `sandwich::estfun(asympt = TRUE)`) and Branch B (self-contained,
#'   deferred to Phase 4).
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
  target_within_first = NULL
) {
  if (fit$type == "longitudinal") {
    return(variance_if_ice(fit, ice_results, target_within_first))
  }

  method <- fit$method

  if (method == "gcomp") {
    return(variance_if_gcomp(
      fit,
      data_a_list,
      preds_list,
      mu_hat,
      target_idx
    ))
  }

  if (method == "matching") {
    return(variance_if_matching(fit, interventions))
  }

  if (method == "ipw") {
    return(variance_if_ipw(
      fit,
      data_a_list,
      preds_list,
      mu_hat,
      target_idx
    ))
  }

  rlang::abort(paste0("Unknown method '", method, "' in variance_if()."))
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
#' Factoring this out keeps `variance_if_gcomp()` and `variance_if_ipw()`
#' identical up to the final Channel-2 call — they both call
#' `build_point_channel_pieces()` and then loop over
#' `apply_model_correction(prep, grad_list[[j]])` (g-comp) or
#' `apply_propensity_correction(prop_prep, grad_list[[j]])` (IPW).
#' Any future Channel-1 fix lands in one place instead of drifting
#' between two near-duplicate loops.
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
  target_idx
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
      rlang::abort(
        "build_point_channel_pieces(): target-population weights sum to 0."
      )
    }
  } else {
    n_target <- sum(target_idx)
    if (n_target == 0L) {
      rlang::abort(
        "build_point_channel_pieces(): target population is empty."
      )
    }
  }

  fit_idx <- resolve_fit_idx(fit, model)
  beta_hat <- stats::coef(model)

  data_a_frames <- lapply(data_a_list, function(da) {
    as.data.frame(da)[target_idx, , drop = FALSE]
  })

  Ch1_list <- vector("list", k)
  grad_list <- vector("list", k)

  # `preds_list[[j]]` is length-`n` (full data). `compute_contrast()`
  # upstream has already masked `target_idx` so that all TRUE rows have
  # non-NA predictions, but non-target rows may still carry NA. We must
  # therefore index `p` by `target_idx` rather than multiplying the
  # full-length `p` by a zero vector — `0 * NA = NA` in R, which would
  # silently corrupt the Channel-1 contribution on non-target rows.
  for (j in seq_len(k)) {
    p <- preds_list[[j]]
    ch1 <- numeric(n)
    if (has_weights) {
      ch1[target_idx] <- n *
        (ext_w[target_idx] / sum_w_target) *
        (p[target_idx] - mu_hat[j])
    } else {
      ch1[target_idx] <- (n / n_target) * (p[target_idx] - mu_hat[j])
    }
    Ch1_list[[j]] <- ch1

    X_star <- iv_design_matrix(model, data_a_frames[[j]])
    eta_star <- as.numeric(X_star %*% beta_hat)
    mu_eta_star <- model$family$mu.eta(eta_star)

    if (has_weights) {
      w_t <- ext_w[target_idx]
      grad_list[[j]] <- as.numeric(crossprod(X_star, w_t * mu_eta_star)) /
        sum_w_target
    } else {
      grad_list[[j]] <- as.numeric(crossprod(X_star, mu_eta_star)) / n_target
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
#' of sync — a caller could pass `resolve_fit_idx(fit, model)` and
#' `nrow(fit$data)` independently, which would silently drift if one of
#' them were computed differently. This helper bundles the two into one
#' call so the alignment is structural, not by convention.
#'
#' **Phase 3 IPW** does *not* use this bundler because it routes its
#' Channel 2 through `prepare_propensity_if()` (WeightIt shortcut), not
#' through `prepare_model_if()`. **Phase 4's self-contained IPW** — which
#' fits a plain weighted GLM as the MSM — will reuse this bundler for
#' the MSM term alongside a parallel `prepare_propensity_if()` call for
#' the propensity term; see `PHASE_4_INTERVENTIONS_SELF_IPW.md`.
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
  target_idx
) {
  pieces <- build_point_channel_pieces(
    fit,
    data_a_list,
    preds_list,
    mu_hat,
    target_idx
  )
  prep <- prepare_model_if(model, pieces$fit_idx, pieces$n)
  list(pieces = pieces, prep = prep)
}


#' G-computation branch of variance_if()
#'
#' @noRd
variance_if_gcomp <- function(
  fit,
  data_a_list,
  preds_list,
  mu_hat,
  target_idx
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
      weights = fit$details$weights
    ))
  }

  bundle <- prepare_point_variance(
    fit,
    model,
    data_a_list,
    preds_list,
    mu_hat,
    target_idx
  )
  pieces <- bundle$pieces
  prep <- bundle$prep

  IF_list <- lapply(seq_along(pieces$grad_list), function(j) {
    pieces$Ch1_list[[j]] +
      apply_model_correction(prep, pieces$grad_list[[j]])$correction
  })
  names(IF_list) <- int_names

  vcov_from_if(IF_list, pieces$n, int_names)
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
#' `data`, and match weights play a dual role — they enter both as the
#' outcome-model `prior.weights` (propagated through `prepare_model_if()`
#' via IWLS working weights) and as the target-population weights for
#' Channel 1. Forcing this into `build_point_channel_pieces()` would
#' require a `scope = c("full", "matched")` switch whose branches would
#' be larger than the inline construction below. The duplication here
#' is intentional — any future Channel-1 fix must be mirrored into both
#' places.
#'
#' @noRd
variance_if_matching <- function(fit, interventions) {
  model <- fit$model
  int_names <- names(interventions)
  k <- length(int_names)

  matched <- as.data.frame(fit$details$matched_data)
  n_m <- nrow(matched)
  cluster <- matched$subclass
  match_w <- matched$weights
  if (is.null(match_w)) {
    match_w <- rep(1, n_m)
  }
  sum_w <- sum(match_w)

  fit_idx <- seq_len(n_m)
  beta_hat <- stats::coef(model)

  treatment <- fit$treatment
  matched_dt <- data.table::as.data.table(matched)

  prep <- prepare_model_if(model, fit_idx, n_m)

  IF_list <- lapply(seq_len(k), function(j) {
    iv <- interventions[[j]]
    da <- apply_intervention(matched_dt, treatment, iv)
    da_df <- as.data.frame(da)
    p <- as.numeric(stats::predict(model, newdata = da_df, type = "response"))
    mu_hat_j <- sum(match_w * p) / sum_w

    IF_vec <- n_m * (match_w / sum_w) * (p - mu_hat_j)

    X_star <- iv_design_matrix(model, da_df)
    eta_star <- as.numeric(X_star %*% beta_hat)
    mu_eta_star <- model$family$mu.eta(eta_star)
    g <- as.numeric(crossprod(X_star, match_w * mu_eta_star)) / sum_w

    IF_vec + apply_model_correction(prep, g)$correction
  })
  names(IF_list) <- int_names

  vcov_from_if(IF_list, n_m, int_names, cluster = cluster)
}


#' ICE branch of variance_if()
#'
#' @description
#' Computes the per-individual IF for chained ICE g-computation by
#' iterating `correct_model()` once per outcome model in the chain. The
#' sensitivity vector `d` is propagated forward (model 0 -> model K),
#' even though the models themselves are fit backward in time
#' (model K -> model 0) — both directions correspond to the back-
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
variance_if_ice <- function(fit, ice_results, target_within_first) {
  int_names <- names(ice_results)
  data <- fit$data
  first_time <- fit$details$time_points[1]
  rows_first <- data[[fit$time]] == first_time
  n <- sum(rows_first)

  IF_list <- lapply(ice_results, function(res) {
    variance_if_ice_one(fit, res, target_within_first)
  })

  vcov_from_if(IF_list, n, int_names)
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
#' from vignette Section 5.4–5.6.
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
  ext_w <- details$weights
  has_weights <- !is.null(ext_w)
  if (has_weights) {
    w_first <- ext_w[rows_first]
    w_t <- w_first[target]
    sum_w_target <- sum(w_t)
    mu_hat <- sum(w_t * pseudo_final[target]) / sum_w_target
    IF_vec <- numeric(n)
    IF_vec[target] <- n *
      (w_t / sum_w_target) *
      (pseudo_final[target] - mu_hat)
  } else {
    n_target <- sum(target)
    mu_hat <- mean(pseudo_final[target], na.rm = TRUE)
    IF_vec <- numeric(n)
    IF_vec[target] <- (n / n_target) * (pseudo_final[target] - mu_hat)
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
    rows_iv_current <- data_iv[[time_col]] == current_time
    iv_data_current <- data_iv[rows_iv_current]
    iv_ids_current <- as.character(iv_data_current[[id_col]])

    # Build the counterfactual design matrix for this time step's
    # model. `iv_design_matrix()` handles both GLMs (via `model.matrix`
    # with stored xlevels) and GAMs (via `predict(..., type =
    # "lpmatrix")`), and accepts data.table input — so ICE can run the
    # same sandwich-variance pipeline on GAM outcome models without any
    # special casing here.
    X_star_k <- iv_design_matrix(model_k, iv_data_current)

    # Counterfactual linear predictor and its derivative w.r.t. eta.
    # `mu_eta_star = dmu/deta` is needed to convert the gradient of the
    # marginal mean w.r.t. beta into the gradient of the linear
    # predictor, which is what the Channel-2 correction operates on.
    eta_star <- as.numeric(X_star_k %*% stats::coef(model_k))
    mu_eta_star <- model_k$family$mu.eta(eta_star)

    if (step_i == 1L) {
      target_in_iv <- match(all_ids[target], iv_ids_current)
      valid_target <- !is.na(target_in_iv)
      target_in_iv <- target_in_iv[valid_target]
      if (has_weights) {
        # `w_t` is already the length-`sum(target)` vector of
        # target-population weights from the IF setup above; align it
        # with `valid_target` to cover cases where some target individuals
        # are missing from the current iv frame.
        target_w <- w_t[valid_target]
        g_k <- as.numeric(
          crossprod(
            X_star_k[target_in_iv, , drop = FALSE],
            target_w * mu_eta_star[target_in_iv]
          )
        ) /
          sum_w_target
      } else {
        g_k <- as.numeric(
          crossprod(
            X_star_k[target_in_iv, , drop = FALSE],
            mu_eta_star[target_in_iv]
          )
        ) /
          n_target
      }
    } else {
      prev_fit_ids <- fit_ids_list[[step_i - 1L]]
      idx_in_all <- id_to_idx[prev_fit_ids]
      rows_in_iv <- match(prev_fit_ids, iv_ids_current)

      keep <- !is.na(idx_in_all) & !is.na(rows_in_iv)
      if (any(keep)) {
        d_prev <- d_vec[idx_in_all[keep]]
        weights_g <- d_prev * mu_eta_star[rows_in_iv[keep]]
        g_k <- as.numeric(
          crossprod(
            X_star_k[rows_in_iv[keep], , drop = FALSE],
            weights_g
          )
        )
      } else {
        g_k <- rep(0, ncol(X_star_k))
      }
    }

    fit_id_idx <- id_to_idx[fit_ids_k]
    res <- correct_model(model_k, g_k, fit_id_idx, n)
    IF_vec <- IF_vec + res$correction
    d_vec <- res$d
  }

  IF_vec
}


#' IPW branch of variance_if()
#'
#' @description
#' Channel 1 + Channel 2 IF for IPW. Shares the Channel-1 / Jacobian
#' construction with the g-comp branch via `build_point_channel_pieces()`
#' — the two branches differ only in whether Channel 2 routes through
#' `apply_model_correction()` (g-comp) or `apply_propensity_correction()`
#' (IPW, combined MSM + propensity correction). The prep step is hoisted
#' out of the intervention loop so `sandwich::estfun()` / `sandwich::bread()`
#' run once regardless of `k`.
#'
#' For a saturated MSM (`Y ~ A`) with a static intervention — the only
#' shape supported in Phase 3 — predictions are constant within each
#' treatment level so Channel 1 is identically zero. The
#' `build_point_channel_pieces()` call still computes `X_star` and `J_j`
#' via `iv_design_matrix()` and `crossprod()`, which is technically
#' redundant in the saturated case but architecturally consistent with
#' Phase 4's non-static interventions (shift, MTP, IPSI, dynamic), where
#' Channel 1 is nonzero and `J_j` carries real information. Keeping one
#' code path handles both shapes without special-casing.
#'
#' @noRd
variance_if_ipw <- function(
  fit,
  data_a_list,
  preds_list,
  mu_hat,
  target_idx
) {
  int_names <- names(data_a_list)

  pieces <- build_point_channel_pieces(
    fit,
    data_a_list,
    preds_list,
    mu_hat,
    target_idx
  )
  prop_prep <- prepare_propensity_if(fit, pieces$fit_idx, pieces$n)

  IF_list <- lapply(seq_along(pieces$grad_list), function(j) {
    pieces$Ch1_list[[j]] +
      apply_propensity_correction(prop_prep, pieces$grad_list[[j]])
  })
  names(IF_list) <- int_names

  vcov_from_if(IF_list, pieces$n, int_names)
}


#' Precompute the gradient-independent IPW Channel 2 ingredients
#'
#' @description
#' Gradient-independent half of the IPW Channel 2 computation (MSM
#' correction + propensity correction). Builds the per-observation
#' adjusted-score IF matrix \eqn{\mathrm{IF}_{\beta,i} = A_{\beta\beta}^{-1}
#' (\psi_{\beta,i} + A_{\beta\alpha}A_{\alpha\alpha}^{-1}\psi_{\alpha,i})},
#' pays the `sandwich::estfun()` / `sandwich::bread()` cost once, and
#' returns an opaque `prop_prep` list that `apply_propensity_correction()`
#' consumes per intervention.
#'
#' Two branches:
#'
#' - **Branch A — WeightIt shortcut.** When `fit$model` is a
#'   `glm_weightit`, `WeightIt::glm_weightit()` already implements
#'   Wooldridge's (2010, eq. 12.41) stacked M-estimation correction
#'   internally. `sandwich::estfun(model, asympt = TRUE)` returns the
#'   **adjusted** per-observation score matrix in a single call, and
#'   `sandwich::bread(model)` returns \eqn{A_{\beta\beta}^{-1}}. Their
#'   product is the full per-observation IF matrix. Covers binary /
#'   multinomial / continuous static IPW with any Mparts-supporting
#'   WeightIt method (`glm`, `cbps`, `ipt`, `ebal`).
#'
#' - **Branch B — self-contained, deferred to Phase 4.** When
#'   `fit$model` is a plain weighted GLM (no `glm_weightit` class),
#'   handles the non-static interventions (shift, MTP, IPSI, dynamic
#'   rules) that WeightIt does not support. See
#'   `PHASE_4_INTERVENTIONS_SELF_IPW.md` for the implementation plan.
#'   Currently aborts.
#'
#' @param fit A `causatr_fit` of method `"ipw"`.
#' @param fit_idx Integer vector of row indices used by the MSM (in
#'   `1..n_total`).
#' @param n_total Integer. Total length of the IF vector.
#'
#' @return A list with `IF_beta` (n_fit x p matrix), `fit_idx`, `n_total`.
#'
#' @noRd
prepare_propensity_if <- function(fit, fit_idx, n_total) {
  if (inherits(fit$model, "glm_weightit")) {
    return(prepare_propensity_if_weightit(fit, fit_idx, n_total))
  }
  prepare_propensity_if_self_contained(fit, fit_idx, n_total)
}


#' Apply a prepared IPW Channel 2 correction to a single gradient
#'
#' @description
#' Returns the per-individual IPW Channel 2 contribution for a single
#' marginal-mean Jacobian `J`:
#'
#' \deqn{\mathrm{Ch.2}_i =
#'   J\,A_{\beta\beta}^{-1}\,\psi_{\beta,i}
#'   + J\,A_{\beta\beta}^{-1}\,A_{\beta\alpha}\,A_{\alpha\alpha}^{-1}\,
#'     \psi_{\alpha,i}}
#'
#' Both terms are carried inside `prep$IF_beta` via the adjusted-score
#' identity (see `prepare_propensity_if()`).
#'
#' @param prep Output of `prepare_propensity_if()`.
#' @param J Numeric `p`-vector. The marginal-mean Jacobian
#'   \eqn{J = \partial\hat\mu/\partial\beta} for one intervention.
#'
#' @return A numeric vector of length `prep$n_total`. Zero outside
#'   `prep$fit_idx`.
#'
#' @noRd
apply_propensity_correction <- function(prep, J) {
  Ch2_fit <- as.numeric(prep$IF_beta %*% J)
  Ch2 <- rep(0, prep$n_total)
  Ch2[prep$fit_idx] <- Ch2_fit
  Ch2
}


#' Branch A: WeightIt shortcut for `prepare_propensity_if()`
#'
#' @noRd
prepare_propensity_if_weightit <- function(fit, fit_idx, n_total) {
  model <- fit$model

  # `sandwich::bread()` follows the sandwich-package convention and
  # returns A^{-1} (NOT n*A^{-1}), so multiplying estfun by bread directly
  # yields per-observation IF_beta,i = A^{-1} psi_i. Dividing by nobs(m)
  # would shrink the IF by a factor of n and the variance by n^2.
  E <- sandwich::estfun(model, asympt = TRUE)
  IF_beta <- E %*% sandwich::bread(model)

  if (nrow(IF_beta) != length(fit_idx)) {
    rlang::abort(
      paste0(
        "prepare_propensity_if() Branch A: estfun() returned ",
        nrow(IF_beta),
        " rows but fit_idx has ",
        length(fit_idx),
        ". This indicates an alignment bug between the MSM fit data ",
        "and fit$details$fit_rows."
      )
    )
  }

  list(IF_beta = IF_beta, fit_idx = fit_idx, n_total = n_total)
}


#' Branch B scaffold: self-contained IPW IF for Phase 4
#'
#' @description
#' Phase 4 will fill in this body. The dispatcher in
#' `prepare_propensity_if()` already routes any non-`glm_weightit` MSM
#' here, so once the body is implemented Phase 4's self-contained IPW
#' (for shift/MTP/IPSI/dynamic interventions) plugs in without touching
#' `variance_if()` itself.
#'
#' Implementation steps (summarised from `PHASE_4_INTERVENTIONS_SELF_IPW.md`):
#'
#' 1. Pull the MSM, the propensity model, and the closure
#'    \eqn{\alpha \mapsto w_i(\alpha)} from `fit$details`.
#' 2. **MSM term.** Reuse `prepare_model_if(msm)` + a per-intervention
#'    `apply_model_correction()` for the MSM correction, keeping
#'    `h_msm = A_{ββ}^{-1} J`.
#' 3. **Cross-derivative.** Compute `A_beta_alpha` numerically via
#'    `numDeriv::jacobian()` on the averaged weighted score, generalising
#'    to any intervention shape.
#' 4. **Propensity term.** Form `g_prop = t(A_beta_alpha) %*% h_msm` and
#'    call `apply_model_correction(prep_prop, g_prop)`.
#' 5. Assemble a `prop_prep` list whose `apply_propensity_correction()`
#'    returns the combined MSM + propensity correction per individual.
#'
#' @noRd
prepare_propensity_if_self_contained <- function(fit, fit_idx, n_total) {
  rlang::abort(
    paste0(
      "Self-contained IPW IF (for non-static interventions: shift, MTP, ",
      "IPSI, dynamic) is planned for Phase 4. Use `ci_method = 'bootstrap'` ",
      "in the meantime, or fit with a WeightIt method that supports the ",
      "Mparts correction (`glm`, `cbps`, `ipt`, `ebal`)."
    )
  )
}
