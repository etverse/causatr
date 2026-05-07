#' Compute joint density-ratio weights for one multivariate intervention
#'
#' @description
#' Multivariate companion to `compute_density_ratio_weights()`. Builds
#' the joint density-ratio weight under sequential factorisation
#' \eqn{f(A_1, \ldots, A_K \mid L) = \prod_k f(A_k \mid A_{1..k-1}, L)}.
#'
#' For each component `k = 1..K`:
#'
#' 1. Build a per-component `newdata` where the upstream treatment
#'    columns `A_1, ..., A_{k-1}` are set to their **intervened**
#'    evaluation values `a_{j,eval}` (for shift / scale this is
#'    `d_j^{-1}(A_j_obs)`; for static it is the static value). This is
#'    the conditioning the k-th factor of the chain-rule numerator
#'    requires:
#'    \deqn{f^*_k(A_k \mid a_{1,eval}, \ldots, a_{k-1,eval}, L).}
#'    For all-static interventions the I-indicator collapses surviving
#'    rows so observed = intervened on those columns and the extra
#'    machinery is a no-op; the routine does it uniformly so the same
#'    code path covers shift-then-anything correctly.
#' 2. Compute the per-component weight by reusing the univariate
#'    `compute_density_ratio_weights()` against a temporarily
#'    rebound `treatment_model` whose `fit_rows` is full (we have
#'    already filtered to the analysis subset).
#'
#' The joint weight is the product across `k`. IPSI components are
#' rejected upstream by `check_intervention_family_compat_mv()` because
#' Kennedy's closed form assumes a single binary treatment with no
#' upstream conditioning.
#'
#' @param treatment_models A `causatr_treatment_models` list from
#'   `fit_treatment_models()`.
#' @param data data.table on the IPW analysis fit-rows (same row set
#'   the per-component models were fit on).
#' @param interventions Named list with one `causatr_intervention` per
#'   treatment component, ordered to match
#'   `names(treatment_models)`. `NULL` entries (per-component natural
#'   course) are allowed.
#' @param estimand Character. ATE only for multivariate; ATT/ATC are
#'   rejected upstream by `check_estimand_trt_compat()`.
#'
#' @return Numeric weight vector of length `nrow(data)`.
#'
#' @noRd
compute_density_ratio_weights_mv <- function(
  treatment_models,
  data,
  interventions,
  estimand = "ATE"
) {
  if (!inherits(treatment_models, "causatr_treatment_models")) {
    rlang::abort(
      "`treatment_models` must be a `causatr_treatment_models` list."
    )
  }
  if (estimand != "ATE") {
    rlang::abort(
      "Multivariate IPW only supports estimand = 'ATE'.",
      class = "causatr_multivariate_estimand"
    )
  }
  K <- length(treatment_models)

  # Natural course (no intervention list): unit weights.
  if (is.null(interventions)) {
    return(rep(1, nrow(data)))
  }

  if (!is.list(interventions) || length(interventions) != K) {
    rlang::abort(
      paste0(
        "Multivariate intervention must be a list with one entry per ",
        "treatment component (expected ",
        K,
        ", got ",
        length(interventions),
        ")."
      )
    )
  }

  # Reorder per-component interventions to match the treatment-model
  # order. `apply_intervention()` already validated the names earlier.
  iv_names <- names(treatment_models)
  if (!setequal(names(interventions), iv_names)) {
    rlang::abort(
      paste0(
        "Multivariate intervention names must match the `treatment` ",
        "argument exactly. Expected: ",
        paste(shQuote(iv_names), collapse = ", "),
        ". Got: ",
        paste(shQuote(names(interventions)), collapse = ", "),
        "."
      )
    )
  }
  interventions <- interventions[iv_names]

  # Reject IPSI components in any slot. IPSI's closed-form weight
  # (Kennedy 2019) presumes the single-treatment Bernoulli case with
  # no upstream conditioning, which the chain-rule factorisation does
  # not preserve.
  for (k in seq_len(K)) {
    iv_k <- interventions[[k]]
    if (!is.null(iv_k) && iv_k$type == "ipsi") {
      rlang::abort(
        c(
          "`ipsi()` is not supported for multivariate IPW.",
          x = paste0(
            "Intervention on component '",
            iv_names[k],
            "' is `ipsi()`."
          ),
          i = "Use `static()` / `shift()` / `scale_by()` / `dynamic()` per component, or fit a single-treatment IPW model."
        ),
        class = "causatr_multivariate_ipsi"
      )
    }
  }

  # Per-component density ratio under the sequential MTP semantics
  # of Diaz, Williams, Hoffman, Schenck (2023). The k-th factor is
  #   w_k = f^num_k(d_k^{-1}(A_k_obs) | A_{1..k-1}_obs, L?) * |Jac d_k^{-1}|
  #         / f_k(A_k_obs | A_{1..k-1}_obs, L)
  # The DENOMINATOR always uses the full-L density model `f_k` fit on
  # baseline confounders + prior treatments. Under `stabilize = "none"`
  # the NUMERATOR evaluates on the same full-L model; under
  # `stabilize = "marginal"` the numerator evaluates on a separately
  # fit `g_k(A_k | A_{1..k-1})` that drops L -- `numerator_models[[k]]`
  # from `fit_treatment_models(stabilize = "marginal")`. Both halves
  # condition on OBSERVED upstream treatments.
  #
  # For natural course on component k under stabilize = "none", the
  # ratio simplifies to 1. Under stabilize = "marginal" the ratio is
  # g_k(A_k | ..._obs) / f_k(A_k | ..._obs, L) -- not 1, because the
  # two models have different conditioning sets.
  num_models <- attr(treatment_models, "numerator_models")
  stabilize <- attr(treatment_models, "stabilize") %||% "none"

  # Running product across k. Each factor is the per-component ratio
  # w_k = f^num_k(A_k | ...) / f_k(A_k | A_{1..k-1}_obs, L).
  # Multiplying sequentially is equivalent to the joint density ratio
  # prod_k w_k = g(d(A) | L) / f(A | L) by the chain-rule factorisation.
  joint_w <- rep(1, nrow(data))
  for (k in seq_len(K)) {
    iv_k <- interventions[[k]]
    tm_k <- treatment_models[[k]]
    a_obs_k <- data[[tm_k$treatment]]
    # Pick numerator model: marginal g_k if stabilize; otherwise same
    # full-L f_k as denominator.
    tm_num_k <- if (stabilize == "marginal") num_models[[k]] else tm_k

    # Denominator: observed density of A_k at observed conditioning
    # under the full-L model f_k.
    f_obs <- evaluate_density(tm_k, a_obs_k, data)
    check_density_positivity(
      f_obs,
      paste0("multivariate density (component ", k, ")")
    )

    if (is.null(iv_k)) {
      # Natural course. Under stabilize = "none" the ratio is 1
      # (identical models); under stabilize = "marginal" it's the
      # density-ratio between the marginal and full-L models.
      if (stabilize == "none") {
        w_k <- rep(1, nrow(data))
      } else {
        f_num <- evaluate_density(tm_num_k, a_obs_k, data)
        w_k <- f_num / f_obs
      }
    } else {
      # Reuse the univariate compatibility check (rejects threshold,
      # bad family combos, etc.); IPSI was already rejected upstream
      # in this function.
      check_intervention_family_compat(iv_k, tm_k)
      iv_type <- iv_k$type

      if (iv_type %in% c("static", "dynamic")) {
        # HT branch on discrete treatment. Numerator density at the
        # target value (evaluated under the numerator model when
        # stabilized).
        # w_k = I(A_k = target) * f^num_k(A_k | ...) / f_k(A_k | ..., L)
        target <- apply_intervention_to_values(iv_k, data, a_obs_k)
        ind <- as.numeric(a_obs_k == target)
        if (stabilize == "none") {
          w_k <- ind / f_obs
        } else {
          # Stabilized HT: indicator * g_k(target | ...) / f_k(A_obs | ..., L)
          # The g evaluated at `target` (the intervention point)
          # equals g evaluated at `a_obs_k` on surviving rows (where
          # the indicator is 1), so we can reuse f_num at a_obs_k.
          f_num <- evaluate_density(tm_num_k, a_obs_k, data)
          w_k <- ind * f_num / f_obs
        }
      } else if (iv_type == "shift") {
        delta <- iv_k$delta
        # a_eval = d^{-1}(A_k) = A_k - delta; no Jacobian for shift.
        a_eval <- a_obs_k - delta
        f_num <- evaluate_density(tm_num_k, a_eval, data)
        warn_intervened_density_near_zero(
          f_num,
          paste0("multivariate shift (component ", k, ")")
        )
        w_k <- f_num / f_obs
      } else if (iv_type == "scale") {
        fct <- iv_k$factor
        if (fct == 0) {
          rlang::abort(
            "`scale_by(0)` collapses the treatment support; not a valid MTP."
          )
        }
        # a_eval = A_k / c; |Jac d^{-1}| = 1/|c|.
        # w_k = f^num_k(A_k / c | ...) / (|c| * f_k(A_k | ..., L))
        a_eval <- a_obs_k / fct
        f_num <- evaluate_density(tm_num_k, a_eval, data)
        warn_intervened_density_near_zero(
          f_num,
          paste0("multivariate scale (component ", k, ")")
        )
        w_k <- (f_num / f_obs) / abs(fct)
      } else {
        rlang::abort(
          paste0(
            "Internal error: multivariate density-ratio has no branch for iv_type='",
            iv_type,
            "'."
          )
        )
      }
    }

    joint_w <- joint_w * w_k
  }

  joint_w
}


#' Construct the stacked alpha closure for multivariate IPW variance
#'
#' @description
#' Multivariate analogue of `make_weight_fn()`. Returns a closure
#' \eqn{w(\alpha) = \prod_k w_k(\alpha_k)} where `alpha` is the
#' concatenation of per-component propensity coefficient vectors. The
#' variance engine feeds this closure to `numDeriv::jacobian()` to
#' compute the cross-derivative \eqn{A_{\beta\alpha}} that the
#' propensity-uncertainty correction needs.
#'
#' Each per-component sub-closure is built by `make_weight_fn()` so it
#' inherits the same per-family branching (Bernoulli HT, Gaussian
#' pushforward, count pushforward, IPSI -- though IPSI is rejected
#' upstream for the multivariate path). Upstream conditioning is
#' baked in by passing a per-component newdata where the prior
#' treatment columns hold their intervened evaluation values; the
#' univariate closure's `X_prop` and `a_obs` are computed from that
#' newdata at closure-creation time.
#'
#' Because each component's score equations depend only on its own
#' alpha block (the propensity models are fit independently), the
#' bread of the stacked propensity system is block-diagonal -- the
#' variance engine handles the propensity correction as a sum of K
#' single-model corrections. This closure is only used for the
#' cross-derivative; the per-block bread is computed by
#' `apply_model_correction()` per component.
#'
#' @inheritParams compute_density_ratio_weights_mv
#'
#' @return A list with components:
#'   \describe{
#'     \item{`weight_fn`}{`function(alpha)` returning a length-n joint
#'       weight vector.}
#'     \item{`offsets`}{Integer (K+1)-vector. `offsets[k]:offsets[k+1]-1`
#'       gives the alpha-block indices for component `k`. Used by the
#'       variance engine to project the cross-derivative onto each
#'       component's bread.}
#'     \item{`alpha_hat`}{Stacked initial alpha (concatenated
#'       per-component `alpha_hat`).}
#'   }
#'
#' @noRd
make_weight_fn_mv <- function(
  treatment_models,
  data,
  interventions,
  estimand = "ATE"
) {
  if (!inherits(treatment_models, "causatr_treatment_models")) {
    rlang::abort(
      "`treatment_models` must be a `causatr_treatment_models` list."
    )
  }
  if (estimand != "ATE") {
    rlang::abort(
      "Multivariate IPW only supports estimand = 'ATE'.",
      class = "causatr_multivariate_estimand"
    )
  }

  K <- length(treatment_models)
  iv_names <- names(treatment_models)

  if (is.null(interventions)) {
    n <- nrow(data)
    return(list(
      weight_fn = function(alpha) rep(1, n),
      offsets = c(1L, 1L),
      alpha_hat = numeric(0)
    ))
  }

  interventions <- interventions[iv_names]

  # Per-component sub-closures + alpha block lengths. Each sub-closure
  # produces the per-component MV weight
  #   w_k(alpha_k) = f_num / f_k(A_k | A_{1..k-1}_obs, L; alpha_k)
  # under the sequential MTP semantics of Diaz et al. 2023. The
  # DENOMINATOR is always the full-L density model `f_k` perturbed
  # through alpha_k. The NUMERATOR depends on `stabilize`:
  #   - "none": f_num is evaluated on the same `f_k` at d_k^{-1}(A_k)
  #     (pushforward) or via the indicator (HT) or equals f_obs
  #     (natural course, ratio = 1).
  #   - "marginal": f_num is a FIXED vector evaluated once from the
  #     numerator model `g_k(A_k | A_{1..k-1}; gamma_hat)` at
  #     closure-creation time. The numerator model's parameters are
  #     held fixed under the variance engine's numDeriv perturbation
  #     of alpha -- same "nuisance-fixed" convention as sigma
  #     (Gaussian) and theta (negbin). Bootstrap refits both gamma
  #     and alpha so the full-variance path is unaffected.
  sub_fns <- vector("list", K)
  block_lens <- integer(K)
  alpha_blocks <- vector("list", K)

  stabilize <- attr(treatment_models, "stabilize") %||% "none"
  num_models <- attr(treatment_models, "numerator_models")
  n_rows <- nrow(data)

  for (k in seq_len(K)) {
    iv_k <- interventions[[k]]
    tm_k <- treatment_models[[k]]
    a_obs_k <- data[[tm_k$treatment]]

    # Propensity design matrix at OBSERVED upstream conditioning. Use
    # the model's actual formula (not ps_formula) so that custom
    # propensity_model_fn that override the formula stay conformable
    # with alpha_hat.
    X_prop <- stats::model.matrix(
      stats::terms(tm_k$model), data = data
    )

    fam_tag <- tm_k$family
    sigma <- tm_k$sigma

    # Precompute the fixed-nuisance numerator vector under stabilization.
    # NULL when not stabilized (closure uses alpha-dependent numerator).
    f_num_fixed <- NULL
    if (stabilize == "marginal") {
      tm_num_k <- num_models[[k]]
      # Pick the evaluation point of A_k in the numerator density:
      # - natural course, HT (surviving rows): a_obs_k
      # - shift: a_obs_k - delta
      # - scale: a_obs_k / factor
      a_eval_num <- if (is.null(iv_k)) {
        a_obs_k
      } else if (iv_k$type %in% c("static", "dynamic")) {
        a_obs_k
      } else if (iv_k$type == "shift") {
        a_obs_k - iv_k$delta
      } else if (iv_k$type == "scale") {
        fct_tmp <- iv_k$factor
        if (fct_tmp == 0) {
          rlang::abort(
            "`scale_by(0)` collapses the treatment support; not a valid MTP."
          )
        }
        a_obs_k / fct_tmp
      } else {
        rlang::abort(
          paste0(
            "Internal error: make_weight_fn_mv stabilize branch has no numerator for iv_type='",
            iv_k$type,
            "'."
          )
        )
      }
      f_num_fixed <- evaluate_density(tm_num_k, a_eval_num, data)
    }

    # Build the per-component closure based on (intervention, family).
    # Natural course -> constant 1 (unstabilized) or f_num_fixed / f_obs(α)
    # (stabilized); static/dynamic -> HT indicator * (f_num_fixed? /) f_obs;
    # shift/scale -> pushforward ratio with |Jac|, using f_num_fixed if
    # stabilized, else the alpha-dependent f_num.
    if (is.null(iv_k)) {
      if (stabilize == "none") {
        # Natural course: w_k = f_k(A_k | obs_hist, L) / f_k(A_k | obs_hist, L) = 1.
        n <- n_rows
        sub_fns[[k]] <- function(alpha) rep(1, n)
      } else {
        # Stabilized natural course: w_k = g_k(A_k | ...)/f_k(A_k | ..., L).
        # Numerator is fixed; denominator via the alpha-dependent
        # full-L density closure (reuse mv_natural_denominator_closure).
        sub_fns[[k]] <- mv_stabilized_closure(
          fam_tag = fam_tag,
          X_prop = X_prop,
          a_obs_k = a_obs_k,
          f_num_fixed = f_num_fixed,
          ind_or_jac = rep(1, n_rows),
          sigma = sigma,
          theta = tm_k$theta,
          trt_levels = tm_k$levels
        )
      }
    } else {
      iv_type <- iv_k$type
      if (iv_type %in% c("static", "dynamic")) {
        target <- apply_intervention_to_values(iv_k, data, a_obs_k)
        ind <- as.numeric(a_obs_k == target)
        if (stabilize == "none") {
          sub_fns[[k]] <- mv_ht_closure(
            fam_tag,
            X_prop,
            a_obs_k,
            ind,
            trt_levels = tm_k$levels
          )
        } else {
          # Stabilized HT: indicator * f_num_fixed / f_obs(alpha).
          sub_fns[[k]] <- mv_stabilized_closure(
            fam_tag = fam_tag,
            X_prop = X_prop,
            a_obs_k = a_obs_k,
            f_num_fixed = f_num_fixed,
            ind_or_jac = ind,
            sigma = sigma,
            theta = tm_k$theta,
            trt_levels = tm_k$levels
          )
        }
      } else if (iv_type == "shift") {
        delta <- iv_k$delta
        a_eval <- a_obs_k - delta
        if (stabilize == "none") {
          sub_fns[[k]] <- mv_pushforward_closure(
            fam_tag,
            X_prop,
            a_obs_k,
            a_eval,
            jac_abs = 1,
            sigma = sigma,
            theta = tm_k$theta
          )
        } else {
          # Stabilized pushforward: f_num_fixed * |Jac| / f_obs(alpha).
          sub_fns[[k]] <- mv_stabilized_closure(
            fam_tag = fam_tag,
            X_prop = X_prop,
            a_obs_k = a_obs_k,
            f_num_fixed = f_num_fixed,
            ind_or_jac = rep(1, n_rows),
            sigma = sigma,
            theta = tm_k$theta
          )
        }
      } else if (iv_type == "scale") {
        fct <- iv_k$factor
        if (fct == 0) {
          rlang::abort(
            "`scale_by(0)` collapses the treatment support; not a valid MTP."
          )
        }
        a_eval <- a_obs_k / fct
        if (stabilize == "none") {
          sub_fns[[k]] <- mv_pushforward_closure(
            fam_tag,
            X_prop,
            a_obs_k,
            a_eval,
            jac_abs = abs(1 / fct),
            sigma = sigma,
            theta = tm_k$theta
          )
        } else {
          sub_fns[[k]] <- mv_stabilized_closure(
            fam_tag = fam_tag,
            X_prop = X_prop,
            a_obs_k = a_obs_k,
            f_num_fixed = f_num_fixed,
            ind_or_jac = rep(abs(1 / fct), n_rows),
            sigma = sigma,
            theta = tm_k$theta
          )
        }
      } else {
        rlang::abort(
          paste0(
            "Internal error: make_weight_fn_mv has no branch for iv_type='",
            iv_type,
            "'."
          )
        )
      }
    }

    alpha_k <- tm_k$alpha_hat
    alpha_blocks[[k]] <- alpha_k
    block_lens[k] <- length(alpha_k)
  }

  # offsets[k]:offsets[k+1]-1 slices the stacked alpha for component k.
  # e.g. three components with p=3,2,4 coefficients -> offsets = c(1,4,6,10).
  offsets <- c(1L, cumsum(block_lens) + 1L)
  alpha_hat <- unlist(alpha_blocks, use.names = FALSE)

  weight_fn <- function(alpha) {
    # Joint weight = prod_k w_k(alpha_k). Each sub-closure captures its own
    # X_prop and data, so they evaluate independently and the joint weight
    # is correct even though alpha is perturbed as a whole by numDeriv.
    w <- 1
    for (k in seq_len(K)) {
      idx <- offsets[k]:(offsets[k + 1L] - 1L)
      alpha_k <- alpha[idx]
      w <- w * sub_fns[[k]](alpha_k)
    }
    w
  }

  list(
    weight_fn = weight_fn,
    offsets = offsets,
    alpha_hat = alpha_hat
  )
}


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


#' Apply an intervention to a treatment vector (no data-frame mutation)
#'
#' @description
#' Lightweight companion to `apply_single_intervention()` in
#' `R/interventions.R`. The difference: this function takes a
#' treatment *vector* and returns a new vector, without copying or
#' mutating any `data.table`. It's what `compute_density_ratio_weights()`
#' and `make_weight_fn()` call to get `a_int` for density evaluation;
#' they cannot use `apply_single_intervention()` directly because that
#' function mutates a `data.table` column in place, which is wasteful
#' inside a `numDeriv::jacobian` loop (we rebuild the closure's `a_int`
#' once at creation time and then the closure itself never touches the
#' data.table).
#'
#' @param intervention A `causatr_intervention` object.
#' @param data `data.table` passed to `dynamic()` rules (only used by
#'   the dynamic branch).
#' @param a_obs Numeric / character / factor vector of observed
#'   treatment values.
#'
#' @return A vector of the same length and type (with the caveats
#'   documented below for dynamic / character interventions).
#'
#' @noRd
apply_intervention_to_values <- function(intervention, data, a_obs) {
  iv <- intervention
  switch(
    iv$type,
    static = {
      # `rep(value, length)` preserves scalar type -- `static("chemo")`
      # returns a character vector, `static(1)` returns a double.
      rep(iv$value, length(a_obs))
    },
    shift = {
      # Only meaningful for numeric treatments. The upstream
      # `check_intervention_family_compat()` has already rejected
      # non-numeric cases (shift on a factor treatment would
      # produce NaN silently).
      a_obs + iv$delta
    },
    scale = a_obs * iv$factor,
    threshold = pmax(pmin(a_obs, iv$upper), iv$lower),
    dynamic = {
      # User-supplied rule. We mirror the type-preservation checks
      # from `apply_single_intervention()` so the density evaluator
      # downstream doesn't get a wrong-type vector. Replicating the
      # checks (rather than calling `apply_single_intervention()`)
      # avoids the unnecessary data.table copy. The error messages
      # point at the same sharp edges -- factor coercion, unknown
      # levels, wrong length.
      new_trt <- iv$rule(data, a_obs)
      if (length(new_trt) != length(a_obs)) {
        rlang::abort(
          paste0(
            "`dynamic()` rule must return a vector of length ",
            length(a_obs),
            " (got ",
            length(new_trt),
            ")."
          )
        )
      }
      if (is.numeric(a_obs) && !is.numeric(new_trt)) {
        rlang::abort(
          "`dynamic()` rule returned non-numeric for a numeric treatment column."
        )
      }
      if (is.factor(a_obs)) {
        if (is.character(new_trt)) {
          unknown <- setdiff(unique(new_trt), levels(a_obs))
          if (length(unknown) > 0L) {
            rlang::abort(
              paste0(
                "`dynamic()` rule returned level(s) not in the factor ",
                "treatment: ",
                paste(shQuote(unknown), collapse = ", ")
              )
            )
          }
          new_trt <- factor(new_trt, levels = levels(a_obs))
        } else if (is.factor(new_trt)) {
          if (!identical(levels(new_trt), levels(a_obs))) {
            rlang::abort(
              "`dynamic()` rule returned a factor with mismatched levels."
            )
          }
        } else {
          rlang::abort(
            "`dynamic()` rule returned a non-factor, non-character value for a factor treatment."
          )
        }
      }
      new_trt
    },
    ipsi = {
      # IPSI callers should go through the closed-form branch in
      # `compute_density_ratio_weights()` / `make_weight_fn()`
      # directly -- reaching here is a bug in the caller.
      rlang::abort(
        "Internal error: `apply_intervention_to_values()` should not be called with an IPSI intervention."
      )
    },
    rlang::abort(paste0("Unknown intervention type: '", iv$type, "'."))
  )
}


#' IPSI closed-form density-ratio weight
#'
#' @description
#' Kennedy (2019) incremental propensity score intervention weight:
#' \deqn{w_i = \frac{\delta A_i + (1 - A_i)}{\delta p_i + (1 - p_i)}.}
#'
#' Derivation. Under IPSI the intervened density is
#' \eqn{f_\delta(a \mid L) = (\delta p)^a (1 - p)^{1-a} / (1 - p + \delta p)}.
#' The density ratio simplifies to
#' \eqn{f_\delta(A_i \mid L_i) / f(A_i \mid L_i) = (\delta A_i + (1 - A_i))
#'  / (\delta p_i + (1 - p_i))}. No density evaluation at a
#' counterfactual treatment value is needed -- the intervention acts
#' directly on the propensity.
#'
#' Positivity: IPSI always preserves positivity (finite `delta`
#' cannot push `p` to 0 or 1), so the denominator is strictly between
#' `min(1, delta)` and `max(1, delta)`. No guard needed.
#'
#' @param a_obs Numeric binary 0/1 treatment vector.
#' @param p Numeric vector of predicted treatment probabilities on the
#'   same rows.
#' @param delta Positive scalar odds multiplier.
#'
#' @return Numeric weight vector, same length as `a_obs`.
#'
#' @noRd
ipsi_weight_formula <- function(a_obs, p, delta) {
  (delta * a_obs + (1 - a_obs)) / (delta * p + (1 - p))
}


#' Bayes-rule numerator f*(L) for the HT estimand weight
#'
#' @description
#' Returns the per-individual multiplier `f_star_i = f(A* | L_i)` that
#' converts the ATE density-ratio weight into an ATT or ATC weight:
#'
#' \deqn{w_i = \mathbb 1\{A_i = a\} \cdot f^\*_i / f(a \mid L_i).}
#'
#' The derivation is the standard Bayes-rule rewrite of
#' \eqn{E[Y^a \mid A = A^\*]} (Imbens 2004; Hernan & Robins Ch. 12).
#' For ATE the target is the whole population and \eqn{f^\* \equiv 1};
#' for ATT the target is the treated and \eqn{f^\*_i = p(L_i)}; for
#' ATC the target is the controls and \eqn{f^\*_i = 1 - p(L_i)}.
#'
#' Only the Bernoulli treatment family is supported because
#' `check_estimand_intervention_compat()` has already rejected ATT /
#' ATC for non-binary static interventions -- hitting the fallback
#' `rlang::abort` here would indicate a missed upstream guard.
#'
#' @param estimand Character scalar in `c("ATE", "ATT", "ATC")`.
#' @param treatment_model A `causatr_treatment_model`.
#' @param fit_data The `fit_rows`-subset `data.table` the caller is
#'   building weights for.
#' @param family_tag Character. The treatment family tag from
#'   `treatment_model$family`.
#'
#' @return Numeric vector of length `nrow(fit_data)` (or a length-1
#'   vector of `1` for ATE, which `R`'s recycling promotes to the
#'   correct per-row constant without allocating).
#'
#' @noRd
ht_bayes_numerator <- function(
  estimand,
  treatment_model,
  fit_data,
  family_tag
) {
  if (estimand == "ATE") {
    return(1)
  }
  if (family_tag != "bernoulli") {
    rlang::abort(
      paste0(
        "Internal error: ATT / ATC Bayes numerator requested for a non-",
        "Bernoulli treatment family ('",
        family_tag,
        "'). `check_estimand_intervention_compat()` should have ",
        "rejected this upstream."
      )
    )
  }
  p <- unname(stats::predict(
    treatment_model$model,
    newdata = fit_data,
    type = "response"
  ))
  if (estimand == "ATT") {
    return(p)
  }
  if (estimand == "ATC") {
    return(1 - p)
  }
  rlang::abort(
    paste0("Internal error: unknown estimand '", estimand, "'.")
  )
}


#' Check that an intervention is compatible with the treatment family
#'
#' @description
#' Rejects nonsensical (intervention, treatment family) combinations
#' up front, with actionable error messages pointing at the right
#' alternative. The compatibility matrix the self-contained IPW
#' engine supports is:
#'
#' | intervention    | bernoulli | gaussian         | categorical | poisson / negbin    |
#' |-----------------|-----------|------------------|-------------|---------------------|
#' | `static()`      | yes HT      | no                | yes HT        | no                   |
#' | `shift()`       | no         | yes smooth ratio   | no           | yes integer only      |
#' | `scale_by()`    | no         | yes smooth ratio   | no           | yes integer-preserving|
#' | `threshold()`   | no         | no (use gcomp)    | no           | no (use gcomp)       |
#' | `dynamic()`     | yes HT      | no (use gcomp)    | yes HT        | no (use gcomp)       |
#' | `ipsi()`        | yes Kennedy | no                | no           | no                   |
#'
#' @param intervention A `causatr_intervention`.
#' @param treatment_model A `causatr_treatment_model`.
#'
#' @return `NULL` invisibly; aborts on invalid combinations.
#'
#' @noRd
check_intervention_family_compat <- function(
  intervention,
  treatment_model
) {
  iv_type <- intervention$type
  fam <- treatment_model$family
  is_count <- fam %in% c("poisson", "negbin")

  if (iv_type == "stochastic") {
    rlang::abort(
      c(
        paste0(
          "`stochastic()` interventions are only supported ",
          "under `estimator = 'gcomp'`."
        ),
        i = paste0(
          "Stochastic interventions require Monte Carlo ",
          "integration over the counterfactual treatment ",
          "distribution, which is only available via the ",
          "g-formula (outcome-model prediction)."
        ),
        i = paste0(
          "Use `causat(..., estimator = 'gcomp')`, or ",
          "rewrite the intervention as a deterministic ",
          "`dynamic()` rule if the policy is actually ",
          "deterministic."
        )
      )
    )
  }

  if (iv_type == "ipsi" && fam != "bernoulli") {
    rlang::abort(
      c(
        "`ipsi()` interventions are only defined for binary (0/1) treatments.",
        i = paste0(
          "The treatment column is classified as '",
          fam,
          "'. Use `shift()` or `scale_by()` for a continuous MTP instead."
        )
      )
    )
  }

  # shift / scale_by require a numeric treatment with a well-defined
  # density: gaussian or count (poisson / negbin).
  if (
    iv_type %in%
      c("shift", "scale") &&
      !fam %in% c("gaussian", "poisson", "negbin")
  ) {
    rlang::abort(
      c(
        paste0(
          "`",
          iv_type,
          "()` interventions require a numeric continuous treatment."
        ),
        i = paste0(
          "The treatment column is classified as '",
          fam,
          "'. Use `static()` for binary / categorical treatments, ",
          "or pass `propensity_family = 'poisson'` / `'negbin'` for count treatments."
        )
      )
    )
  }

  # threshold on any non-gcomp-compatible family
  if (iv_type == "threshold" && fam != "gaussian" && !is_count) {
    rlang::abort(
      c(
        paste0(
          "`threshold()` interventions require a numeric continuous treatment."
        ),
        i = paste0(
          "The treatment column is classified as '",
          fam,
          "'. Use `static()` for binary / categorical treatments."
        )
      )
    )
  }

  # Count-specific guards: shift must be integer, scale must preserve
  # integer support, threshold and dynamic are rejected.
  if (is_count) {
    if (iv_type == "static") {
      rlang::abort(
        c(
          paste0(
            "`static(v)` on a count treatment is degenerate for IPW."
          ),
          i = "The Horvitz-Thompson indicator weight is zero for almost all observations.",
          i = "Use `shift()` for integer shifts, or switch to `estimator = 'gcomp'`."
        )
      )
    }
    if (iv_type == "threshold") {
      rlang::abort(
        c(
          "`threshold()` on a count treatment is not supported by the IPW engine.",
          i = "The pushforward of a count density under a boundary clamp is a mixed measure.",
          i = "Use `estimator = 'gcomp'`."
        )
      )
    }
    if (iv_type == "dynamic") {
      rlang::abort(
        c(
          "`dynamic()` rules on count treatments are not supported by the IPW engine.",
          i = "Use `shift()` for integer shifts, or switch to `estimator = 'gcomp'` for deterministic rules."
        )
      )
    }
    if (iv_type == "shift") {
      delta <- intervention$delta
      if (delta != round(delta)) {
        rlang::abort(
          c(
            paste0(
              "`shift(",
              delta,
              ")` is not integer-valued."
            ),
            i = paste0(
              "Count treatments (",
              fam,
              ") require integer shift deltas because `dpois()` / `dnbinom()` return 0 at non-integer arguments."
            ),
            i = "Use `estimator = 'gcomp'` for fractional shifts on count data."
          )
        )
      }
    }
    if (iv_type == "scale") {
      fct <- intervention$factor
      a_obs <- treatment_model$model$data[[treatment_model$treatment]]
      if (is.null(a_obs)) {
        a_obs <- stats::model.response(stats::model.frame(
          treatment_model$model
        ))
      }
      # The density ratio evaluates dpois(a_obs / factor, lambda).
      # For dpois to return a non-zero value, a_obs / factor must be
      # a non-negative integer for every observed treatment value.
      inv_scaled <- a_obs / fct
      if (!all(inv_scaled == round(inv_scaled) & inv_scaled >= 0)) {
        rlang::abort(
          c(
            paste0(
              "`scale_by(",
              fct,
              ")` does not produce integer inverse values (A / ",
              fct,
              ") for all observed treatment values."
            ),
            i = paste0(
              "Count treatments (",
              fam,
              ") require that A / factor is a non-negative integer for every observation, because the density ratio evaluates the pmf at A / factor."
            ),
            i = "Use `estimator = 'gcomp'` for non-integer-preserving scales on count data."
          )
        )
      }
    }
  }

  if (iv_type == "static" && fam == "gaussian") {
    rlang::abort(
      c(
        "`static(v)` on a continuous treatment is degenerate for IPW.",
        i = "No observations lie exactly at `v`, so the Horvitz-Thompson weight is zero almost surely.",
        i = "Use `shift()` or `scale_by()` to move the whole treatment distribution, or switch to `estimator = 'gcomp'`."
      )
    )
  }

  if (iv_type == "threshold" && fam == "gaussian") {
    rlang::abort(
      c(
        "`threshold(lo, hi)` on a continuous treatment is not supported by the IPW engine.",
        i = "The pushforward of a continuous density under a boundary clamp has point masses at `lo` and `hi`, so the density ratio w.r.t. the fitted `f(a|l)` is not well-defined.",
        i = "Use `estimator = 'gcomp'` -- the clamped-treatment counterfactual is well-defined there via predict-then-average on the outcome model."
      )
    )
  }

  if (iv_type == "dynamic" && fam == "gaussian") {
    rlang::abort(
      c(
        "`dynamic()` rules on continuous treatments are not supported by the IPW engine.",
        i = "A deterministic per-individual target is a Dirac per individual; the density ratio is degenerate.",
        i = "Use `shift()` / `scale_by()` for smooth MTPs, or switch to `estimator = 'gcomp'` for deterministic rules."
      )
    )
  }

  invisible(NULL)
}
