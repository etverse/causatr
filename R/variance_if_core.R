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
