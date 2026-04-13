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
#' - The five primitives — `bread_inv()`, `iv_design_matrix()`,
#'   `correct_model()`, `vcov_from_if()`, `variance_if_numeric()`.
#' - The dispatcher `variance_if()` with branches for g-computation,
#'   matching (cluster-robust), and ICE (chained backward recursion).
#' - The IPW-specific `correct_propensity()` (added in Phase C).
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
#' construction inside `correct_propensity()` Branch B.
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
    error = function(e) MASS::ginv(XtWX)
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
#' @param model A fitted model object.
#' @param newdata Data frame of counterfactual observations.
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


#' Channel 2 correction for one model
#'
#' @description
#' Computes the per-individual contribution of a single nuisance model
#' \eqn{k} to the influence function of \eqn{\hat\mu}, given the gradient
#' \eqn{g_k = \partial\hat\mu/\partial\beta_k}. Following the three-line
#' template of vignette Section 4:
#'
#' \deqn{h = A_k^{-1} g_k, \qquad d_i = X_i^T h, \qquad
#'       \mathrm{correction}_i = n \cdot d_i \cdot r^{\mathrm{score}}_i}
#'
#' The factor of \eqn{n} arises because \eqn{A_k^{-1} = n(X^TWX)^{-1}}
#' (see vignette Section 5.5). The score residual
#' \eqn{r^{\mathrm{score}}_i = (Y_i - \hat\mu_i) \cdot (d\mu/d\eta)/V(\mu)}
#' is the one that appears in the GLM score \eqn{\psi_i = X_i\,r^{\mathrm{score}}_i}.
#' For canonical links \eqn{(d\mu/d\eta)/V(\mu) = 1} and the score
#' residual collapses to the response residual; for non-canonical links
#' (probit, cloglog, Gamma-log, Poisson-sqrt, ...) the link-scale factor
#' must be applied or the IF is miscalibrated.
#'
#' Returns the correction vector along with `d` and `h` because:
#'
#' - **ICE** consumes `d` to build the next model's gradient via
#'   \eqn{g_{k+1} = \sum_j d_{k,j}\,(d\mu/d\eta)\,X^*_{k+1,j}}.
#' - **`correct_propensity()` Branch B** consumes `h` to construct the
#'   derived gradient \eqn{g^{\mathrm{prop}} = A_{\beta\alpha}^T h_{\mathrm{msm}}}
#'   for the propensity correction term.
#'
#' G-comp and matching only read `$correction` and ignore the other slots.
#'
#' @param model A fitted model with a `family` object (GLM or GAM).
#' @param gradient Numeric `p`-vector. The sensitivity gradient
#'   \eqn{g = \partial\hat\mu/\partial\beta}.
#' @param fit_idx Integer vector of length \eqn{n_{\mathrm{fit}}}. Row indices
#'   in `1..n_total` corresponding to the rows the model was fit on.
#' @param n_total Integer. The total denominator used to scale the IF
#'   (matches the length of the IF vector, not necessarily `nobs(model)`).
#'
#' @return A list with components:
#'   \describe{
#'     \item{`correction`}{Numeric vector of length `n_total`. Zero for
#'       rows outside `fit_idx`. `n_total * d_i * r_score_i` on fit rows.}
#'     \item{`d`}{Numeric vector of length `n_total`. The per-individual
#'       sensitivity \eqn{d_i = X_i^T h}; zero off `fit_idx`.}
#'     \item{`h`}{Numeric `p`-vector. The bread-projected gradient
#'       \eqn{h = A_k^{-1} g_k}.}
#'   }
#'
#' @noRd
correct_model <- function(model, gradient, fit_idx, n_total) {
  X_fit <- stats::model.matrix(model)
  B_inv <- bread_inv(model, X_fit)

  h <- as.numeric(B_inv %*% gradient)
  d_fit <- as.numeric(X_fit %*% h)

  # Score residual = response residual * (dmu/deta) / V(mu) * prior_w.
  # For canonical links the factor (dmu/deta)/V is 1; for non-canonical
  # (probit, cloglog, Gamma-log, Poisson-sqrt, ...) it is not, and using
  # the bare response residual gives a miscalibrated IF.
  #
  # We obtain it directly from the GLM internals via
  #   working_residual * working_weight = ((y-mu)/mu_eta) * (mu_eta^2/V * pw)
  #                                     = (y-mu) * mu_eta/V * pw
  # which matches sandwich::estfun(glm) row-for-row. Recomputing through
  # family$mu.eta / family$variance is correct in theory but introduces
  # ~1e-6 numerical drift relative to the IWLS internals.
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

  d_full <- rep(0, n_total)
  d_full[fit_idx] <- d_fit

  correction <- rep(0, n_total)
  correction[fit_idx] <- n_total * d_fit * r_score

  list(correction = correction, d = d_full, h = h)
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
  k <- length(IF_list)

  if (!is.null(cluster)) {
    if (length(cluster) != length(IF_list[[1]])) {
      rlang::abort(
        "`cluster` length must match the IF vector length in `vcov_from_if()`."
      )
    }
    cl <- as.factor(cluster)
    IF_list <- lapply(IF_list, function(IF) {
      as.numeric(rowsum(IF, cl, reorder = FALSE))
    })
  }

  vcov_mat <- matrix(0, k, k)
  for (j in seq_len(k)) {
    for (l in j:k) {
      vcov_mat[j, l] <- sum(IF_list[[j]] * IF_list[[l]]) / n^2
      if (j != l) {
        vcov_mat[l, j] <- vcov_mat[j, l]
      }
    }
  }
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
#' \eqn{V_1 + V_2 = (1/n^2)\sum_i \mathrm{Ch.1}_i^2 + J V_\beta J^T} and
#' return their sum. This **drops the cross-term**
#' \eqn{(2/n^2)\sum_i \mathrm{Ch.1}_i\,\mathrm{Ch.2}_i} from the
#' decomposition in vignette Section 3.3, so the variance is slightly
#' miscalibrated. A `rlang::warn()` is emitted so the user knows.
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

  if (!is.null(weights)) {
    sum_w <- sum(weights)
    target_w <- ifelse(target_idx, 0, 0)
    target_w[target_idx] <- weights / sum_w
  } else {
    n_t <- sum(target_idx)
    target_w <- ifelse(target_idx, 1 / n_t, 0)
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
    fit_idx <- which(fit$details$fit_rows)
    na_action <- model$na.action
    if (!is.null(na_action)) {
      fit_idx <- fit_idx[-na_action]
    }
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
      IF_list <- lapply(seq_len(k), function(j) {
        Ch2_fit <- as.numeric(IF_beta %*% J[j, ]) * n_total
        IF <- Ch1_list[[j]]
        IF[fit_idx] <- IF[fit_idx] + Ch2_fit
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

  V1 <- vapply(
    Ch1_list,
    function(IF) sum(IF^2) / n_total^2,
    numeric(1)
  )
  V_beta <- tryCatch(stats::vcov(model), error = function(e) NULL)
  if (is.null(V_beta)) {
    V_beta <- diag(0, length(beta_hat))
  }
  V2 <- J %*% V_beta %*% t(J)

  vcov_mat <- V2 + diag(V1, nrow = k)
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
#'   is delegated to `correct_propensity()`, which dispatches between
#'   Branch A (WeightIt shortcut via `sandwich::estfun(asympt = TRUE)`)
#'   and Branch B (self-contained, deferred to Phase 4).
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
  k <- length(int_names)
  n <- nrow(fit$data)

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
      target_idx
    ))
  }

  ext_w <- fit$details$weights
  has_weights <- !is.null(ext_w)

  if (has_weights) {
    w_target <- ifelse(target_idx, ext_w, 0)
    sum_w_target <- sum(w_target[target_idx])
  } else {
    n_target <- sum(target_idx)
  }

  fit_idx <- which(fit$details$fit_rows)
  na_action <- model$na.action
  if (!is.null(na_action)) {
    fit_idx <- fit_idx[-na_action]
  }
  beta_hat <- stats::coef(model)

  data_a_frames <- lapply(data_a_list, function(da) {
    as.data.frame(da)[target_idx, , drop = FALSE]
  })

  IF_list <- lapply(seq_len(k), function(j) {
    p <- preds_list[[j]]
    if (has_weights) {
      IF_vec <- n * (w_target / sum_w_target) * (p - mu_hat[j])
    } else {
      IF_vec <- (n / n_target) * ifelse(target_idx, p - mu_hat[j], 0)
    }

    X_star <- iv_design_matrix(model, data_a_frames[[j]])
    eta_star <- as.numeric(X_star %*% beta_hat)
    mu_eta_star <- model$family$mu.eta(eta_star)

    if (has_weights) {
      w_t <- ext_w[target_idx]
      g <- as.numeric(crossprod(X_star, w_t * mu_eta_star)) / sum_w_target
    } else {
      g <- as.numeric(crossprod(X_star, mu_eta_star)) / n_target
    }

    IF_vec + correct_model(model, g, fit_idx, n)$correction
  })
  names(IF_list) <- int_names

  vcov_from_if(IF_list, n, int_names)
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

    IF_vec + correct_model(model, g, fit_idx, n_m)$correction
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

  ext_w <- details$weights
  has_weights <- !is.null(ext_w)
  if (has_weights) {
    w_first <- ext_w[rows_first]
    w_target <- ifelse(target, w_first, 0)
    sum_w_target <- sum(w_target)
    mu_hat <- sum(w_target * pseudo_final) / sum_w_target
    IF_vec <- n * (w_target / sum_w_target) * (pseudo_final - mu_hat)
  } else {
    n_target <- sum(target)
    mu_hat <- mean(pseudo_final[target], na.rm = TRUE)
    IF_vec <- (n / n_target) * ifelse(target, pseudo_final - mu_hat, 0)
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

    rows_iv_current <- data_iv[[time_col]] == current_time
    iv_data_current <- data_iv[rows_iv_current]
    iv_ids_current <- as.character(iv_data_current[[id_col]])

    pred_terms <- stats::delete.response(stats::terms(model_k))
    X_star_k <- tryCatch(
      stats::model.matrix(pred_terms, data = iv_data_current),
      error = function(e) {
        stats::model.matrix(pred_terms, data = as.data.frame(iv_data_current))
      }
    )

    eta_star <- as.numeric(X_star_k %*% stats::coef(model_k))
    mu_eta_star <- model_k$family$mu.eta(eta_star)

    if (step_i == 1L) {
      target_in_iv <- match(all_ids[target], iv_ids_current)
      valid_target <- !is.na(target_in_iv)
      target_in_iv <- target_in_iv[valid_target]
      if (has_weights) {
        target_w <- w_target[target][valid_target]
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
      d_prev <- d_vec[idx_in_all]

      keep <- !is.na(idx_in_all) & !is.na(rows_in_iv) & d_prev != 0
      if (any(keep)) {
        weights_g <- d_prev[keep] * mu_eta_star[rows_in_iv[keep]]
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
#' Channel 1 + Channel 2 IF for IPW. Channel 1 is identical in form to the
#' g-comp branch (a marginal-mean residual). Channel 2 is the **combined**
#' MSM correction + propensity correction, delegated entirely to
#' `correct_propensity()` — see its docstring for the two branches and
#' the WeightIt vs self-contained dispatch.
#'
#' For a saturated MSM (`Y ~ A`) with a static intervention — the only
#' shape supported in Phase 3 — predictions are constant within each
#' treatment level so Channel 1 is identically zero and the entire
#' variance comes from `correct_propensity()`. For Phase 4's non-static
#' interventions (shift, MTP, IPSI, dynamic) Channel 1 will be nonzero
#' and the same branch handles both.
#'
#' @noRd
variance_if_ipw <- function(
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

  if (has_weights) {
    w_target <- ifelse(target_idx, ext_w, 0)
    sum_w_target <- sum(w_target[target_idx])
  } else {
    n_target <- sum(target_idx)
  }

  fit_idx <- which(fit$details$fit_rows)
  na_action <- model$na.action
  if (!is.null(na_action)) {
    fit_idx <- fit_idx[-na_action]
  }
  beta_hat <- stats::coef(model)

  data_a_frames <- lapply(data_a_list, function(da) {
    as.data.frame(da)[target_idx, , drop = FALSE]
  })

  IF_list <- lapply(seq_len(k), function(j) {
    p <- preds_list[[j]]
    if (has_weights) {
      IF_vec <- n * (w_target / sum_w_target) * (p - mu_hat[j])
    } else {
      IF_vec <- (n / n_target) * ifelse(target_idx, p - mu_hat[j], 0)
    }

    X_star <- iv_design_matrix(model, data_a_frames[[j]])
    eta_star <- as.numeric(X_star %*% beta_hat)
    mu_eta_star <- model$family$mu.eta(eta_star)

    if (has_weights) {
      w_t <- ext_w[target_idx]
      J_j <- as.numeric(crossprod(X_star, w_t * mu_eta_star)) / sum_w_target
    } else {
      J_j <- as.numeric(crossprod(X_star, mu_eta_star)) / n_target
    }

    IF_vec + correct_propensity(fit, J_j, fit_idx, n)
  })
  names(IF_list) <- int_names

  vcov_from_if(IF_list, n, int_names)
}


#' IPW Channel 2 (MSM + propensity correction)
#'
#' @description
#' Returns the **combined** Channel 2 contribution for IPW per individual:
#' the MSM correction term plus the propensity correction term from
#' vignette Section 4.2:
#'
#' \deqn{\mathrm{Ch.2}_i =
#'   J\,A_{\beta\beta}^{-1}\,\psi_{\beta,i}
#'   + J\,A_{\beta\beta}^{-1}\,A_{\beta\alpha}\,A_{\alpha\alpha}^{-1}\,
#'     \psi_{\alpha,i}}
#'
#' Two branches:
#'
#' - **Branch A — WeightIt shortcut.** When `fit$model` is a
#'   `glm_weightit`, `WeightIt::glm_weightit()` already implements
#'   Wooldridge's (2010, eq. 12.41) stacked M-estimation correction
#'   internally. `sandwich::estfun(model, asympt = TRUE)` returns the
#'   **adjusted** score matrix
#'   \eqn{\psi_{\beta,i} + A_{\beta\alpha}A_{\alpha\alpha}^{-1}\psi_{\alpha,i}}
#'   in a single matrix; `sandwich::bread(model)/nobs(model)` returns
#'   \eqn{A_{\beta\beta}^{-1}}. Their product times \eqn{J} yields both
#'   correction terms in five lines, with no need to compute
#'   \eqn{A_{\beta\alpha}} ourselves. Covers binary / multinomial /
#'   continuous static IPW with any Mparts-supporting WeightIt method
#'   (`glm`, `cbps`, `ipt`, `ebal`).
#'
#' - **Branch B — self-contained, deferred to Phase 4.** When
#'   `fit$model` is a plain weighted GLM (no `glm_weightit` class),
#'   handles the non-static interventions (shift, MTP, IPSI, dynamic
#'   rules) that WeightIt does not support. Plan: (1) call
#'   `correct_model(msm_model, J, ...)` for the MSM term; (2) compute
#'   \eqn{A_{\beta\alpha}} numerically via `numDeriv::jacobian()` on the
#'   averaged weighted score (so the same code generalises across all
#'   intervention types); (3) build the derived gradient
#'   \eqn{g^{\mathrm{prop}} = A_{\beta\alpha}^T h_{\mathrm{msm}}} and
#'   call `correct_model(propensity_model, g_prop, ...)` for the
#'   propensity correction. See `PHASE_4_INTERVENTIONS_SELF_IPW.md`
#'   (Section "Self-contained IPW influence function") for the
#'   step-by-step plan and `VARIANCE_REFACTOR.qmd` for the derivation.
#'   The dispatcher is in place; the body currently aborts.
#'
#' @param fit A `causatr_fit` of method `"ipw"`.
#' @param J Numeric `p`-vector. The marginal-mean Jacobian
#'   \eqn{J = \partial\hat\mu/\partial\beta} for one intervention.
#' @param fit_idx Integer vector of row indices used by the MSM (in
#'   `1..n_total`).
#' @param n_total Integer. Total length of the IF vector.
#'
#' @return A numeric vector of length `n_total`. Zero outside `fit_idx`.
#'
#' @noRd
correct_propensity <- function(fit, J, fit_idx, n_total) {
  if (inherits(fit$model, "glm_weightit")) {
    return(correct_propensity_weightit(fit, J, fit_idx, n_total))
  }
  correct_propensity_self_contained(fit, J, fit_idx, n_total)
}


#' Branch A: WeightIt shortcut for `correct_propensity()`
#'
#' @noRd
correct_propensity_weightit <- function(fit, J, fit_idx, n_total) {
  model <- fit$model

  # `sandwich::bread()` follows the sandwich-package convention and
  # returns A^{-1} (NOT n*A^{-1}), so multiplying estfun by bread directly
  # yields per-observation IF_beta,i = A^{-1} psi_i. Dividing by nobs(m)
  # would shrink the IF by a factor of n and the variance by n^2.
  E <- sandwich::estfun(model, asympt = TRUE)
  IF_beta <- E %*% sandwich::bread(model)
  Ch2_fit <- as.numeric(IF_beta %*% J)

  if (length(Ch2_fit) != length(fit_idx)) {
    rlang::abort(
      paste0(
        "correct_propensity() Branch A: estfun() returned ",
        length(Ch2_fit),
        " rows but fit_idx has ",
        length(fit_idx),
        ". This indicates an alignment bug between the MSM fit data ",
        "and fit$details$fit_rows."
      )
    )
  }

  Ch2 <- rep(0, n_total)
  Ch2[fit_idx] <- Ch2_fit
  Ch2
}


#' Branch B scaffold: self-contained IPW IF for Phase 4
#'
#' @description
#' Phase 4 will fill in this body. The dispatcher in `correct_propensity()`
#' already routes any non-`glm_weightit` MSM here, so once the body is
#' implemented Phase 4's self-contained IPW (for shift/MTP/IPSI/dynamic
#' interventions) plugs in without touching `variance_if()` itself.
#'
#' Implementation steps (summarised from `VARIANCE_REFACTOR.qmd` and
#' `PHASE_4_INTERVENTIONS_SELF_IPW.md`):
#'
#' 1. Pull the MSM, the propensity model, and the closure
#'    \eqn{\alpha \mapsto w_i(\alpha)} from `fit$details`. (Phase 4 must
#'    add `propensity_model` and `weight_fn` slots in `fit_ipw()`.)
#' 2. **MSM term.** Call `correct_model(msm, J, fit_idx, n_total)` to get
#'    the MSM correction *and* `h_msm = A_{ββ}^{-1} J`.
#' 3. **Cross-derivative.** Define
#'    `psi_beta_bar(alpha)` = `(1/n) sum_i psi_β,i(alpha, beta_hat)`, the
#'    averaged weighted MSM score as a function of the propensity
#'    parameters with `beta` fixed. Compute
#'    `A_beta_alpha = -numDeriv::jacobian(psi_beta_bar, x = alpha_hat)`.
#'    The numerical jacobian generalises to any intervention shape
#'    (analytic derivatives differ for static vs shift vs IPSI vs dynamic).
#' 4. **Propensity term.** Form the derived gradient
#'    `g_prop = t(A_beta_alpha) %*% h_msm` and call
#'    `correct_model(propensity_model, g_prop, fit_idx_prop, n_total)`.
#'    Add its `$correction` to the MSM correction; return the sum.
#'
#' @noRd
correct_propensity_self_contained <- function(fit, J, fit_idx, n_total) {
  rlang::abort(
    paste0(
      "correct_propensity() Branch B (self-contained IPW for non-static ",
      "interventions: shift, MTP, IPSI, dynamic) is planned for Phase 4. ",
      "See PHASE_4_INTERVENTIONS_SELF_IPW.md, section 'Self-contained IPW ",
      "influence function', for the implementation steps. The dispatcher ",
      "in correct_propensity() already routes here when fit$model is not ",
      "a glm_weightit object."
    )
  )
}
