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
#' - The Channel-2 primitives -- `bread_inv()`, `iv_design_matrix()`,
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
#'   recursion), and IPW (density-ratio stacked sandwich).
#'
#' @section Layout:
#' Functions appear in this file in the order below. Grep for the name
#' to jump; line numbers are deliberately omitted so this index does not
#' rot as the file grows.
#'
#' **Primitives (shared across methods).**
#' - `bread_inv()` -- `(X'WX)^{-1}` for GLMs, `Vp` for GAMs, with a warned
#'   `MASS::ginv()` fallback on singular bread.
#' - `iv_design_matrix()` -- counterfactual model matrix under an
#'   intervention (GLM/GAM split).
#' - `resolve_fit_idx()` -- maps `fit$details$fit_rows` to model-local row
#'   indices and asserts the `na.action` invariant.
#'
#' **Channel-2 per-model correction (prep/apply split).**
#' - `prepare_model_if()` -- one bread solve per model, reused across
#'   interventions (O(1) instead of O(k)).
#' - `apply_model_correction()` -- per-intervention application step.
#' - `correct_model()` -- thin single-gradient wrapper used by ICE.
#'
#' **Aggregation and fallbacks.**
#' - `vcov_from_if()` -- `crossprod`-based aggregator with optional
#'   `cluster =` sum-then-square for matching.
#' - `variance_if_numeric()` -- two-tier numerical fallback (Tier 1
#'   recovers the full IF via `sandwich::estfun + bread`; Tier 2 falls
#'   back to `V1 + J V_\beta J^T`).
#'
#' **Dispatcher and per-method branches.**
#' - `variance_if()` -- single entry point, routes to one of the four.
#' - `build_point_channel_pieces()`, `prepare_point_variance()` --
#'   shared point-treatment plumbing for gcomp / matching.
#' - `variance_if_gcomp()` -- g-computation branch.
#' - `variance_if_matching()` -- matching branch; cluster-robust on
#'   `subclass`.
#' - `variance_if_ice()` / `variance_if_ice_one()` -- ICE forward
#'   sensitivity recursion; cycles models rather than interventions so
#'   prep-hoisting does not apply.
#' - `variance_if_ipw()` / `compute_ipw_if_self_contained_one()` -- IPW
#'   density-ratio stacked sandwich; per-intervention workhorse assembles
#'   Channel 1 + MSM correction + propensity correction via
#'   `apply_model_correction()` and a numerical cross-derivative
#'   \eqn{A_{\beta\alpha}} from `numDeriv::jacobian()`.
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
#' This is the building block for `correct_model()` and for the MSM
#' plus propensity Channel-2 pair inside
#' `compute_ipw_if_self_contained_one()`.
#'
#' @param model A fitted model with a `family` object (GLM or GAM).
#' @param X_fit Design matrix from `model.matrix(model)`.
#'
#' @return A `p x p` matrix.
#'
#' @noRd
bread_inv <- function(model, X_fit) {
  if (inherits(model, "gam")) {
    # Properly fitted `mgcv::gam` objects always carry `$Vp` (the
    # Bayesian posterior covariance of the smooth coefficients,
    # which plays the role of the inverse bread for the penalised
    # IWLS fit). If `$Vp` is absent we cannot fall back to a GLM-style
    # bread on `model.matrix(model)`: for a GAM that matrix is the
    # linear-predictor design (including basis-expanded smooth
    # columns) but the penalty has warped the IWLS weights in a way
    # the naive `X'WX` solve cannot recover. Silently falling through
    # would silently miscompute the sandwich variance. Abort loudly.
    if (is.null(model$Vp)) {
      rlang::abort(
        paste0(
          "`bread_inv()`: GAM fit object is missing `$Vp`; cannot ",
          "compute sandwich bread. Rebuild the fit with `mgcv::gam()` ",
          "(the default path produces `$Vp`), or switch to ",
          "`ci_method = 'bootstrap'`."
        ),
        class = "causatr_gam_missing_vp"
      )
    }
    return(model$Vp)
  }

  # Prefer `stats::weights(model, type = "working")` when it returns a
  # non-empty vector -- this is the public accessor that sandwich::bread()
  # and summary(glm)$cov.unscaled both use, so routing through it keeps
  # us aligned with sandwich-ecosystem conventions even for GLM subclasses
  # (e.g. glm_weightit, glmnet's glm wrapper) that override $weights to
  # carry something other than the IWLS working weights. See R7 in the
  # 2026-04-15 critical review. Fall back to $weights and finally to a
  # family-based recomputation.
  w_iwls <- tryCatch(
    stats::weights(model, type = "working"),
    error = function(e) NULL
  )
  if (is.null(w_iwls) || length(w_iwls) == 0L) {
    w_iwls <- model$weights
  }
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
      # replicate -- without throttling the console fills up. The
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
        class = "causatr_singular_bread",
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
#' Accepts `data.table` directly on both branches -- no coercion to
#' `data.frame`. Verified against GLM, GLM-with-interaction, and
#' `na.action`-triggering inputs; `stats::model.matrix()` dispatches
#' through `model.frame()` which handles `data.table` without going
#' through the `[.data.table` subset that can evaluate bare symbols in
#' the frame's environment. `delete.response()` strips the LHS from
#' `terms(model)` so the counterfactual frame does not need the
#' response column.
#'
#' @param model A fitted model object.
#' @param newdata Data frame or data.table of counterfactual observations.
#'
#' @return A design matrix (`n_new x p`).
#'
#' @noRd
iv_design_matrix <- function(model, newdata) {
  if (inherits(model, "gam")) {
    return(stats::predict(model, newdata = newdata, type = "lpmatrix"))
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


#' Precompute the IF ingredients for a multinomial propensity model
#'
#' @description
#' Parallel to `prepare_model_if()` but for `nnet::multinom` models,
#' which are not GLMs and lack `$family`, `$linear.predictors`, and
#' IWLS working weights. Computes bread and score from the multinomial
#' log-likelihood directly.
#'
#' The multinomial logit score for observation i and non-reference
#' class k is
#' \deqn{s_{ik} = (I(A_i = k) - p_{ik}) X_i,}
#' stacked into a (K-1)*p vector. The expected information (bread) is
#' \deqn{H = \sum_i \mathrm{diag}(p_i) - p_i p_i^T) \otimes X_i X_i^T,}
#' where the Kronecker product is over the K-1 non-reference classes.
#'
#' The return value has the same shape as `prepare_model_if()` so
#' `apply_model_correction()` can consume it transparently.
#'
#' @param model A fitted `nnet::multinom` model.
#' @param fit_idx Integer vector. Row indices in `1..n_total`.
#' @param n_total Integer. Total row count for scaling.
#'
#' @return A list with `model`, `X_fit`, `B_inv`, `r_score`, `fit_idx`,
#'   `n_total`. `X_fit` is the n x ((K-1)*p) stacked design matrix so
#'   the standard `apply_model_correction()` algebra works.
#'
#' @noRd
prepare_model_if_multinom <- function(model, fit_idx, n_total) {
  X_base <- stats::model.matrix(model)
  n <- nrow(X_base)
  p <- ncol(X_base)

  # Predicted probabilities: n x K matrix.
  prob_raw <- stats::predict(model, type = "probs")
  cc <- stats::coef(model)
  if (is.null(dim(cc))) {
    Km1 <- 1L
    # 2-level: prob_raw is P(second level), need full matrix
    trt_levels <- model$lev
    prob_mat <- cbind(1 - prob_raw, prob_raw)
  } else {
    Km1 <- nrow(cc)
    prob_mat <- prob_raw
    trt_levels <- model$lev
  }

  # Response indicators: n x (K-1) matrix, I(A_i = level_k) for each
  # non-reference class. Reference level is `model$lev[1]`.
  response <- stats::model.response(stats::model.frame(model))
  response_char <- as.character(response)
  # Non-reference levels are trt_levels[2:K]
  non_ref <- trt_levels[-1]
  Y_mat <- matrix(0, nrow = n, ncol = Km1)
  for (k in seq_len(Km1)) {
    Y_mat[, k] <- as.numeric(response_char == non_ref[k])
  }

  # Non-reference probabilities: n x (K-1)
  P_non_ref <- prob_mat[, -1, drop = FALSE]

  # Score residual matrix: n x (K-1), each column is (I(A=k) - p_k)
  R_mat <- Y_mat - P_non_ref

  # Stacked score: n x ((K-1)*p). Row i of R_score is the Kronecker
  # product of the K-1 residuals with X_i. We stack column-major
  # within each row to match the row-major alpha flattening:
  # (r_{i,1}*X_i, r_{i,2}*X_i, ...).
  X_stacked <- matrix(0, nrow = n, ncol = Km1 * p)
  for (k in seq_len(Km1)) {
    cols <- ((k - 1L) * p + 1L):(k * p)
    X_stacked[, cols] <- X_base * R_mat[, k]
  }

  # Bread: information matrix of the multinomial logit. For the (j,k)
  # block (p x p each):
  #   H_{jk} = -sum_i (delta_{jk} * p_{ij} - p_{ij} * p_{ik}) * X_i X_i'
  # where j, k are 1-indexed non-reference classes.
  # We build H as a (Km1*p) x (Km1*p) matrix.
  H <- matrix(0, nrow = Km1 * p, ncol = Km1 * p)
  for (j in seq_len(Km1)) {
    for (k in seq_len(Km1)) {
      j_cols <- ((j - 1L) * p + 1L):(j * p)
      k_cols <- ((k - 1L) * p + 1L):(k * p)
      if (j == k) {
        w_jk <- P_non_ref[, j] * (1 - P_non_ref[, j])
      } else {
        w_jk <- -P_non_ref[, j] * P_non_ref[, k]
      }
      H[j_cols, k_cols] <- crossprod(X_base, X_base * w_jk)
    }
  }

  B_inv <- tryCatch(
    solve(H),
    error = function(e) {
      rlang::warn(
        c(
          "Multinomial information matrix is singular; using `MASS::ginv()` fallback.",
          i = "Inspect the propensity model for collinear confounders."
        ),
        class = "causatr_singular_bread",
        .frequency = "once",
        .frequency_id = "causatr_multinom_bread_singular"
      )
      MASS::ginv(H)
    }
  )

  # `r_score` in the standard prep is a length-n vector of per-obs
  # scores. For multinomial, the score is (K-1)*p-dimensional per obs.
  # `apply_model_correction()` computes `d_fit = X_fit %*% h` then
  # correction = (d_fit * r_score) summed over fit_idx. For the
  # stacked system, `X_fit` = X_stacked (n x (Km1*p)), `r_score` = 1
  # (a scalar) because the score is already embedded in X_stacked.
  # This is equivalent to saying: each row of X_stacked IS the per-obs
  # score vector (the estimating equation). The prep/apply split needs
  # `r_score * X_fit` = the n x (Km1*p) score matrix. We achieve this
  # by setting r_score = rep(1, n).
  r_score <- rep(1, n)

  list(
    model = model,
    X_fit = X_stacked,
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
#' `compute_ipw_if_self_contained_one()` for the MSM-to-propensity
#' cross-term).
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
    # R8 (2026-04-15 review): `as.factor(cluster)` would reorder the
    # levels via `sort(unique(x))`, which for an integer cluster vector
    # uses numeric sort and for a character cluster vector uses
    # lexicographic sort. Either way it stops tracking the row order of
    # `IF_mat`, and `reorder = FALSE` cannot save us once the factor
    # levels themselves have been permuted. Use a first-seen factor so
    # `rowsum(..., reorder = FALSE)` groups consistently with IF_mat's
    # row ordering.
    cluster_f <- factor(cluster, levels = unique(cluster))
    IF_mat <- rowsum(IF_mat, cluster_f, reorder = FALSE)
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

  # Row-subset each counterfactual frame to the target population.
  # `da[target_idx, , drop = FALSE]` is polymorphic: it does row-selection
  # on both data.frame and data.table (the trailing empty `j` forces
  # data.frame's `[` into row-index mode instead of column-select).
  # Avoids the previous `as.data.frame(da)[...]` roundtrip.
  data_a_frames <- lapply(
    data_a_list,
    function(da) da[target_idx, , drop = FALSE]
  )

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
  target_within_first = NULL,
  ipw_bundles = NULL,
  ipw_fit_idx = NULL,
  ipw_n_total = NULL
) {
  if (fit$type == "longitudinal") {
    return(variance_if_ice(fit, ice_results, target_within_first))
  }

  estimator <- fit$estimator

  if (estimator == "gcomp") {
    return(variance_if_gcomp(
      fit,
      data_a_list,
      preds_list,
      mu_hat,
      target_idx
    ))
  }

  if (estimator == "matching") {
    return(variance_if_matching(fit, interventions))
  }

  if (estimator == "ipw") {
    return(variance_if_ipw(
      fit,
      ipw_bundles,
      target_idx,
      mu_hat,
      ipw_fit_idx,
      ipw_n_total
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

    IF_vec <- n_m * (match_w / sum_w) * (p - mu_hat_j)

    X_star <- iv_design_matrix(model, da)
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
  mu_hat <- sum(w_t * pseudo_final[target]) / sum_w_target
  IF_vec <- numeric(n)
  IF_vec[target] <- n *
    (w_t / sum_w_target) *
    (pseudo_final[target] - mu_hat)

  # Per-time-step id -> external-weight lookup. Needed for the
  # step > 1 cascade gradient, which must multiply by the weight of
  # the PREVIOUS model's fit row (w_{k-1} enters through A_{k-1,k} =
  # E[\partial s_{k-1}/\partial \beta_k]; see variance-theory vignette
  # Section 5.4). Without this factor the sandwich silently drops a
  # term and underestimates the ICE SE under non-uniform weights by
  # ~2x on heterogeneous designs. Only built when external weights
  # are present -- unweighted ICE is unchanged.
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
    rows_iv_current <- data_iv[[time_col]] == current_time
    iv_data_current <- data_iv[rows_iv_current]
    iv_ids_current <- as.character(iv_data_current[[id_col]])

    # Build the counterfactual design matrix for this time step's
    # model. `iv_design_matrix()` handles both GLMs (via `model.matrix`
    # with stored xlevels) and GAMs (via `predict(..., type =
    # "lpmatrix")`), and accepts data.table input -- so ICE can run the
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
      # Unified gradient using the same weighted form as the IF above:
      # `w_t` is rep(1, n_target) in the unweighted branch (constructed
      # at the top of this function), so dividing by `sum_w_target`
      # recovers the unweighted `/ n_target` scale exactly. See B7 in
      # the 2026-04-15 review.
      target_w <- w_t[valid_target]
      g_k <- as.numeric(
        crossprod(
          X_star_k[target_in_iv, , drop = FALSE],
          target_w * mu_eta_star[target_in_iv]
        )
      ) /
        sum_w_target
    } else {
      prev_fit_ids <- fit_ids_list[[step_i - 1L]]
      idx_in_all <- id_to_idx[prev_fit_ids]
      rows_in_iv <- match(prev_fit_ids, iv_ids_current)

      keep <- !is.na(idx_in_all) & !is.na(rows_in_iv)
      if (any(keep)) {
        d_prev <- d_vec[idx_in_all[keep]]
        # Cascade gradient g_k^eff = A_{k-1,k}^T h_{k-1}
        #   = sum_j w_{k-1,j} * d_{k-1,j} * X^*_{k,j} * mu_eta_{k,j}
        # The w_{k-1,j} factor comes from \partial s_{k-1,j}/\partial
        # \beta_k, which carries the prior weights from the (k-1)-th
        # fit because the pseudo-outcome model at step k-1 is weighted.
        # Unweighted ICE has w == 1 so this collapses to the previous
        # formula; see variance-theory vignette Section 5.4.
        w_prev <- if (has_weights) {
          unname(w_at_step[[step_i - 1L]][prev_fit_ids[keep]])
        } else {
          rep(1, sum(keep))
        }
        weights_g <- w_prev * d_prev * mu_eta_star[rows_in_iv[keep]]
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
  n_total
) {
  int_names <- names(bundles)
  data <- fit$data
  tm <- fit$details$treatment_model
  propensity_model <- tm$model
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
  ext_w_sub <- if (is.null(ext_w)) NULL else ext_w[fit_rows]

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

    beta_hat <- stats::coef(msm_model)
    eta_star <- as.numeric(X_star %*% beta_hat)
    mu_eta_star <- msm_model$family$mu.eta(eta_star)
    preds_sub <- msm_model$family$linkinv(eta_star)

    # Channel 1: n_sub * (w_i / sum_w_target) * (pred_i - mu_hat_j),
    # zero off target. `target_sub` masks the contribution; the
    # unweighted branch uses `w_target_vec = 1{target}`, and the
    # weighted branch uses `ext_w_i * 1{target}`. Scaled by
    # `n_sub` (the `n` passed to `vcov_from_if`).
    Ch1_i <- n_sub * (w_target_vec / sum_w_target) * (preds_sub - mu_hat_j)
    Ch1_i[!target_sub] <- 0

    # Marginal-mean Jacobian J = d mu_hat_j / d beta.
    # For a weighted Hajek marginal mean over the target population,
    # J = (1/sum_w_target) * sum_target (X_star_i * mu_eta_i * w_i).
    # This matches `build_point_channel_pieces()`'s gradient.
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
    wfn <- make_weight_fn(
      tm,
      fit_data_local,
      b$intervention,
      estimand = estimand
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

    compute_ipw_if_self_contained_one(
      msm_model = msm_model,
      propensity_model = propensity_model,
      weight_fn = wfn,
      J = J,
      Ch1_i = Ch1_i,
      fit_idx = seq_len(n_sub),
      fit_idx_ps = seq_len(n_sub),
      n_total = n_sub
    )
  })
  names(IF_list) <- int_names

  vcov_from_if(IF_list, n_sub, int_names)
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
  # Phi_bar(alpha) = (1/n_fit) sum_i psi_beta_i(alpha, beta_hat).
  # For a GLM with canonical / non-canonical link,
  #   psi_beta_i = X_i * w_i(alpha) * (Y_i - mu_i) * mu_eta_i / var_mu_i
  # where mu, mu_eta, var_mu, Y - mu are all functions of beta_hat
  # (fixed inside the closure). Only `w_i(alpha)` varies with alpha,
  # so numDeriv only has to re-run the weight formula per perturbation.
  beta_hat <- stats::coef(msm_model)
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
  # For multinomial models `coef()` returns a matrix; flatten to match
  # `make_weight_fn()`'s convention (row-major: `as.vector(t(coef_mat))`).
  alpha_hat_raw <- stats::coef(propensity_model)
  if (!is.null(dim(alpha_hat_raw))) {
    alpha_hat <- as.vector(t(alpha_hat_raw))
  } else {
    alpha_hat <- alpha_hat_raw
  }
  # Negative-Hessian convention: A_{beta, alpha} = -(1/n) sum d psi/d alpha.
  # numDeriv::jacobian(phi_bar, alpha) = d phi_bar/d alpha = +(1/n) sum d psi/d alpha.
  # Flipping the sign gives A_{beta, alpha}.
  A_beta_alpha <- -numDeriv::jacobian(phi_bar, x = alpha_hat)

  # ---- Propensity correction via causatr's primitive --------------
  # h_msm_true = A_{beta, beta}^{-1} J. msm_res$h holds
  # (X'WX)^{-1} J = A_bb^{-1} J / n_fit, so multiply by n_fit to
  # recover the "true" h.
  h_msm_true <- n_fit * msm_res$h
  # g_prop = A_{beta, alpha}^T h_msm_true is a p_alpha-vector. Feeding
  # it as the "sensitivity gradient" to `apply_model_correction()` on
  # the propensity model returns
  #   prop_res$correction_i = g_prop^T A_{alpha, alpha}^{-1} psi_{alpha, i}
  # which is exactly the quantity we need (up to sign; see below).
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
      )
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
      )
    )
  }

  # Block-lower-triangular M-estimation result:
  #   IF_beta_i = A_bb^{-1}(psi_beta_i - A_{beta, alpha} A_aa^{-1} psi_alpha_i)
  # i.e. the propensity correction is SUBTRACTED. See the docstring
  # note on the sign convention.
  Ch1_i + msm_res$correction - prop_res$correction
}
