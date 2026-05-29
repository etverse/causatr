#' Per-component MV closure for HT branch (static / dynamic on discrete)
#'
#' @description
#' Builds the alpha-closure for a per-component multivariate weight
#' factor under static or dynamic intervention on a discrete
#' (Bernoulli) treatment:
#' \deqn{w_k(\alpha_k) = \frac{\mathbb 1\{A_k = \mathrm{target}_k\}}{f_k(A_k \mid \mathrm{cond\_obs}, L; \alpha_k)}.}
#' Multivariate is ATE-only, so the Bayes numerator is 1. The
#' indicator is invariant under alpha perturbation and is captured at
#' closure-creation time.
#'
#' @param fam_tag Character family.
#' @param X_prop_obs Denominator-side propensity design matrix.
#' @param a_obs_k Observed treatment vector.
#' @param ind Length-n 0/1 indicator vector `I(a_obs_k == target_k)`.
#' @return `function(alpha)` returning the per-component weight.
#' @noRd
mv_ht_closure <- function(
  fam_tag,
  X_prop_obs,
  a_obs_k,
  ind,
  trt_levels = NULL
) {
  # `force()` ensures each captured variable is evaluated at closure-creation
  # time, preventing lazy-evaluation bugs when called inside a loop over k.
  force(X_prop_obs)
  force(a_obs_k)
  force(ind)
  force(trt_levels)
  if (fam_tag == "bernoulli") {
    return(function(alpha) {
      p_obs <- stats::plogis(as.numeric(X_prop_obs %*% alpha))
      # f_obs = P(A_k = a_obs_k | cond_obs, L; alpha) for a Bernoulli(p_obs).
      f_obs <- ifelse(a_obs_k == 1, p_obs, 1 - p_obs)
      ind / f_obs
    })
  }
  if (fam_tag == "categorical") {
    # Multinomial HT closure. The flattened `alpha` is row-major
    # `as.vector(t(coef_mat))` per `fit_treatment_model()`'s convention,
    # so reshape to (K-1) x p to recover per-non-reference log-odds.
    # Probabilities for non-reference levels follow the softmax;
    # reference level is `1 / (1 + sum(exp(eta)))`. We index per-row
    # by the observed treatment value to recover f_obs.
    if (is.null(trt_levels)) {
      rlang::abort(
        "mv_ht_closure: categorical branch requires `trt_levels`."
      )
    }
    K_lev <- length(trt_levels)
    Km1 <- K_lev - 1L
    p_cols <- ncol(X_prop_obs)
    n_obs <- length(a_obs_k)
    a_obs_char <- as.character(a_obs_k)
    # col_idx[i] is the column of prob_mat that holds P(A = a_obs[i] | L).
    # Precomputed at closure-creation time; invariant under alpha perturbation.
    col_idx <- match(a_obs_char, trt_levels)
    return(function(alpha) {
      alpha_mat <- matrix(alpha, nrow = Km1, ncol = p_cols, byrow = TRUE)
      # eta: n x (K-1) log-odds vs the reference level.
      eta <- X_prop_obs %*% t(alpha_mat)
      exp_eta <- exp(eta)
      denom <- 1 + rowSums(exp_eta)
      # prob_mat: n x K; column 1 = reference level P = 1/denom.
      prob_mat <- cbind(1 / denom, exp_eta / denom)
      # Two-column matrix index selects prob_mat[i, col_idx[i]] per row.
      f_obs <- prob_mat[cbind(seq_len(n_obs), col_idx)]
      ind / f_obs
    })
  }
  rlang::abort(
    paste0("mv_ht_closure: unsupported family '", fam_tag, "' for HT branch.")
  )
}


#' Per-component MV closure for shift / scale (smooth pushforward)
#'
#' @description
#' Builds the alpha-closure for a per-component multivariate weight
#' factor under shift or scale intervention on a continuous (Gaussian
#' or count) treatment under sequential MTP semantics:
#' \deqn{w_k(\alpha_k) = \frac{f_k(d_k^{-1}(A_k) \mid A_{1..k-1}^{\mathrm{obs}}, L; \alpha_k) \cdot |\mathrm{Jac}|}{f_k(A_k \mid A_{1..k-1}^{\mathrm{obs}}, L; \alpha_k)}.}
#' Both numerator and denominator condition on observed upstream
#' treatments (Diaz et al. 2023 Sec 2). No intervened-conditioning
#' substitution is needed; only the k-th argument (A_k vs
#' d_k^{-1}(A_k)) differs between numerator and denominator.
#'
#' @param fam_tag Character family.
#' @param X_prop Propensity design matrix at OBSERVED upstream
#'   conditioning.
#' @param a_obs_k Observed treatment vector.
#' @param a_eval Inverse-map value of A_k (`A_k - delta` for shift,
#'   `A_k / factor` for scale).
#' @param jac_abs Absolute Jacobian of the inverse map.
#' @param sigma Residual SD (Gaussian only).
#' @param theta NB dispersion (negbin only).
#' @return `function(alpha)` returning the per-component weight.
#' @noRd
mv_pushforward_closure <- function(
  fam_tag,
  X_prop,
  a_obs_k,
  a_eval,
  jac_abs,
  sigma = NULL,
  theta = NULL
) {
  # `force()` prevents late-binding bugs inside the k-loop.
  force(X_prop)
  force(a_obs_k)
  force(a_eval)
  force(jac_abs)
  force(sigma)
  force(theta)
  # Sequential MTP semantics: both numerator f_k(d^{-1}(A_k) | obs_hist, L)
  # and denominator f_k(A_k | obs_hist, L) evaluate at the SAME
  # conditioning linear predictor (X_prop %*% alpha). `a_eval` = d^{-1}(A_k)
  # is fixed at closure-creation time; only the mean parameter mu varies.
  if (fam_tag == "gaussian") {
    return(function(alpha) {
      mu <- as.numeric(X_prop %*% alpha)
      f_num <- stats::dnorm(a_eval, mean = mu, sd = sigma)
      f_den <- stats::dnorm(a_obs_k, mean = mu, sd = sigma)
      (f_num / f_den) * jac_abs
    })
  }
  if (fam_tag == "poisson") {
    return(function(alpha) {
      lam <- as.numeric(exp(X_prop %*% alpha))
      f_num <- stats::dpois(a_eval, lam)
      f_den <- stats::dpois(a_obs_k, lam)
      (f_num / f_den) * jac_abs
    })
  }
  if (fam_tag == "negbin") {
    return(function(alpha) {
      lam <- as.numeric(exp(X_prop %*% alpha))
      f_num <- stats::dnbinom(a_eval, mu = lam, size = theta)
      f_den <- stats::dnbinom(a_obs_k, mu = lam, size = theta)
      (f_num / f_den) * jac_abs
    })
  }
  rlang::abort(
    paste0(
      "mv_pushforward_closure: unsupported family '",
      fam_tag,
      "' for pushforward branch."
    )
  )
}


#' Per-component MV closure under stabilize = "marginal"
#'
#' @description
#' Stabilized per-component multivariate weight factor, with a FIXED
#' numerator density vector. The numerator `f_num_fixed` was
#' precomputed from the separate numerator model
#' `g_k(A_k \mid A_{1..k-1}; \hat\gamma)` at closure-creation time;
#' its parameters are held fixed under the variance engine's numDeriv
#' perturbation of alpha (same "nuisance-fixed" convention as sigma /
#' theta). Only the denominator density depends on alpha:
#' \deqn{w_k(\alpha_k) = \mathrm{ind\_or\_jac}_i \cdot \frac{f^{\mathrm{num}}_{\mathrm{fixed}, i}}{f_k(A_{k,i} \mid A_{1..k-1}^{\mathrm{obs}}, L_i; \alpha_k)}.}
#'
#' `ind_or_jac` is a length-`n` vector that carries either:
#' - the HT indicator `I(A_k = target)` for static / dynamic
#'   interventions (zero outside the target set);
#' - a constant `|Jac|` per row for pushforward (shift / scale)
#'   interventions;
#' - all-`1`s for natural course.
#'
#' This one helper covers all three branches because the
#' alpha-dependence always lives in the denominator only under the
#' "nuisance-fixed" convention.
#'
#' Bootstrap variance correctly captures the full uncertainty
#' (including gamma) because `refit_ipw()` re-fits both the
#' denominator and numerator models on each replicate.
#'
#' @param fam_tag Character. Treatment density family tag.
#' @param X_prop Design matrix of the full-L propensity model at
#'   observed upstream conditioning.
#' @param a_obs_k Observed treatment vector for component k.
#' @param f_num_fixed Precomputed numerator density vector of length
#'   `n` (evaluated on the numerator model at the appropriate A_k
#'   point for the intervention).
#' @param ind_or_jac Length-`n` vector carrying the intervention-
#'   specific indicator / Jacobian (see description).
#' @param sigma Residual SD (Gaussian only).
#' @param theta NB dispersion (negbin only).
#' @param trt_levels Character vector of factor levels (categorical
#'   only).
#' @return `function(alpha)` returning a length-`n` per-component
#'   weight vector.
#' @noRd
mv_stabilized_closure <- function(
  fam_tag,
  X_prop,
  a_obs_k,
  f_num_fixed,
  ind_or_jac,
  sigma = NULL,
  theta = NULL,
  trt_levels = NULL
) {
  # `force()` prevents late-binding bugs inside the k-loop.
  force(X_prop)
  force(a_obs_k)
  force(f_num_fixed)
  force(ind_or_jac)
  force(sigma)
  force(theta)
  force(trt_levels)
  # `f_num_fixed` is precomputed from the numerator model at closure-creation
  # time and held constant under numDeriv perturbation of alpha. Only the
  # denominator f_obs(alpha) varies; the formula is:
  #   w_k = ind_or_jac * f_num_fixed / f_k(A_k | ..., L; alpha)
  if (fam_tag == "bernoulli") {
    return(function(alpha) {
      p <- stats::plogis(as.numeric(X_prop %*% alpha))
      f_obs <- ifelse(a_obs_k == 1, p, 1 - p)
      ind_or_jac * f_num_fixed / f_obs
    })
  }
  if (fam_tag == "categorical") {
    if (is.null(trt_levels)) {
      rlang::abort(
        "mv_stabilized_closure: categorical branch requires `trt_levels`."
      )
    }
    K_lev <- length(trt_levels)
    Km1 <- K_lev - 1L
    p_cols <- ncol(X_prop)
    n_obs <- length(a_obs_k)
    a_obs_char <- as.character(a_obs_k)
    col_idx <- match(a_obs_char, trt_levels)
    return(function(alpha) {
      alpha_mat <- matrix(alpha, nrow = Km1, ncol = p_cols, byrow = TRUE)
      eta <- X_prop %*% t(alpha_mat)
      exp_eta <- exp(eta)
      denom <- 1 + rowSums(exp_eta)
      prob_mat <- cbind(1 / denom, exp_eta / denom)
      f_obs <- prob_mat[cbind(seq_len(n_obs), col_idx)]
      ind_or_jac * f_num_fixed / f_obs
    })
  }
  if (fam_tag == "gaussian") {
    return(function(alpha) {
      mu <- as.numeric(X_prop %*% alpha)
      f_obs <- stats::dnorm(a_obs_k, mean = mu, sd = sigma)
      ind_or_jac * f_num_fixed / f_obs
    })
  }
  if (fam_tag == "poisson") {
    return(function(alpha) {
      lam <- as.numeric(exp(X_prop %*% alpha))
      f_obs <- stats::dpois(a_obs_k, lam)
      ind_or_jac * f_num_fixed / f_obs
    })
  }
  if (fam_tag == "negbin") {
    return(function(alpha) {
      lam <- as.numeric(exp(X_prop %*% alpha))
      f_obs <- stats::dnbinom(a_obs_k, mu = lam, size = theta)
      ind_or_jac * f_num_fixed / f_obs
    })
  }
  rlang::abort(
    paste0(
      "mv_stabilized_closure: unsupported family '",
      fam_tag,
      "'."
    )
  )
}
