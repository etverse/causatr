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
  estimand = "ATE",
  trim = 1
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

    # Per-component truncation before the cross-component product.
    # Truncating individual density ratios at the source prevents a
    # single extreme component from dominating the joint weight
    # (Cole & Hernán 2008; lmtp uses per-component 0.999 by default).
    w_k <- truncate_weights(w_k, trim)
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
  estimand = "ATE",
  trim = 1
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

  # Precompute per-component truncation thresholds at alpha_hat.
  # Per-component truncation before the cross-component product matches
  # the semantics of compute_density_ratio_weights_mv(), which
  # truncates each component's density ratio via truncate_weights()
  # before multiplying across k. Fixing each threshold under numDeriv
  # perturbation ensures the sandwich SE reflects the truncated weight
  # surface, not a shifting truncation boundary.
  comp_thresholds <- vector("list", K)
  if (trim < 1) {
    for (kk in seq_len(K)) {
      iv_kk <- interventions[[kk]]
      tm_kk <- treatment_models[[kk]]
      if (is.null(iv_kk)) {
        next
      }
      w_kk <- compute_density_ratio_weights(
        treatment_model = tm_kk,
        data = data,
        intervention = iv_kk,
        estimand = estimand,
        trim = 1
      )
      comp_thresholds[[kk]] <- stats::quantile(w_kk, trim, names = FALSE)
    }
  }

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
      stats::terms(tm_k$model),
      data = data
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

    # Wrap the sub-closure with per-component truncation at the
    # precomputed threshold. Matches compute_density_ratio_weights_mv()
    # semantics: each component's density ratio is truncated before
    # the cross-component product.
    if (trim < 1 && !is.null(comp_thresholds[[k]])) {
      sub_fns[[k]] <- wrap_closure_with_trim(sub_fns[[k]], comp_thresholds[[k]])
    }
  }

  # offsets[k]:offsets[k+1]-1 slices the stacked alpha for component k.
  # e.g. three components with p=3,2,4 coefficients -> offsets = c(1,4,6,10).
  offsets <- c(1L, cumsum(block_lens) + 1L)
  alpha_hat <- unlist(alpha_blocks, use.names = FALSE)

  weight_fn <- function(alpha) {
    # Joint weight = prod_k w_k(alpha_k). Each sub-closure captures its own
    # X_prop and data, so they evaluate independently and the joint weight
    # is correct even though alpha is perturbed as a whole by numDeriv.
    # Per-component truncation is baked into each sub-closure (via
    # wrap_closure_with_trim), so no post-product truncation is needed.
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
