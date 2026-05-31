#' Construct the stacked-alpha closure for longitudinal IPW variance
#'
#' @description
#' Longitudinal companion to `make_weight_fn()` /
#' `make_weight_fn_mv()`. For a fixed intervention applied at every
#' period, builds a closure
#' \deqn{W_i(\alpha) = \prod_{t=1}^{T} w_t(\alpha_t)}
#' where `\alpha = (\alpha_1, \ldots, \alpha_T)` is the concatenation
#' of per-period propensity coefficients. Each per-period sub-closure
#' is built by `make_weight_fn()` (univariate) or
#' `make_weight_fn_mv()` (multivariate) applied to the period-t
#' subset and density model(s).
#'
#' For multivariate treatment, each period's alpha block contains
#' K component sub-blocks (from `make_weight_fn_mv()`), so the
#' total stacked alpha has T x K blocks in period-major order.
#'
#' The T propensity model groups are fit on **disjoint** row subsets
#' (one subset per period), so the bread of the stacked propensity
#' system is block-diagonal -- the variance engine handles the
#' propensity correction as a sum of T (x K) single-model
#' corrections. The closure itself is what `numDeriv::jacobian()`
#' consumes to compute the cross-derivative `A_{beta, alpha}`.
#'
#' Each sub-closure returns a length-`n_id` vector ordered to match
#' the canonical `ids_first` ordering -- per-period ids are mapped
#' back to that index space before multiplication, mirroring the
#' alignment in `compute_longitudinal_weights()`.
#'
#' @param treatment_models_by_time Named list of T per-period models.
#'   Each element is either a `causatr_treatment_model` (univariate)
#'   or a `causatr_treatment_models` list of K components (MV).
#' @param fit_data_by_time Named list of K per-period data subsets.
#' @param ids_first Character vector of individual ids in the
#'   canonical first-period order.
#' @param id_col Character. ID column name.
#' @param intervention A `causatr_intervention` object or `NULL`.
#' @param estimand Character. ATE only for longitudinal IPW.
#'
#' @return A list with components:
#'   \describe{
#'     \item{`weight_fn`}{`function(alpha)` returning a length-`n_id`
#'       cumulative product weight vector.}
#'     \item{`offsets`}{Integer (K+1)-vector. `offsets[k]:offsets[k+1]
#'       - 1` gives the alpha-block indices for period k.}
#'     \item{`alpha_hat`}{Stacked initial alpha (concatenation of
#'       per-period `alpha_hat`).}
#'   }
#'
#' @noRd
make_weight_fn_longitudinal <- function(
  treatment_models_by_time,
  fit_data_by_time,
  ids_first,
  id_col,
  intervention,
  estimand = "ATE",
  numerator_models_by_time = NULL,
  trim = 1
) {
  if (estimand != "ATE") {
    rlang::abort(
      "Longitudinal IPW only supports estimand = 'ATE'.",
      class = "causatr_longitudinal_estimand"
    )
  }

  n_periods <- length(treatment_models_by_time)
  n_id <- length(ids_first)
  stabilize <- !is.null(numerator_models_by_time)

  # Detect multivariate: per-period models are `causatr_treatment_models`
  # (plural) rather than scalar `causatr_treatment_model` objects.
  is_mv <- inherits(
    treatment_models_by_time[[1L]],
    "causatr_treatment_models"
  )

  # Natural course unstabilized: weight is identically 1 regardless of
  # alpha. Still return a closure so the variance engine's loop has a
  # uniform contract. `numDeriv::jacobian` on a constant function
  # returns zero, which is correct -- natural course carries no
  # propensity-uncertainty.
  #
  # Natural course **stabilized** is non-trivial: the per-period
  # weight is `g_k(A_obs | ...) / f_k(A_obs | ..., L; alpha_k)`, which
  # depends on alpha through the denominator. We fall through to the
  # general per-period closure builder below.
  if (is.null(intervention) && !stabilize) {
    return(list(
      weight_fn = function(alpha) rep(1, n_id),
      offsets = c(1L, 1L),
      alpha_hat = numeric(0)
    ))
  }

  # Precompute per-period truncation thresholds at alpha_hat.
  # Per-period truncation before the cumulative product matches the
  # semantics of compute_longitudinal_weights(), which truncates each
  # period's density ratio via compute_density_ratio_weights(trim=).
  # Fixing the threshold under numDeriv perturbation ensures the
  # sandwich SE reflects the truncated weight surface, not a shifting
  # truncation boundary.
  period_thresholds <- vector("list", n_periods)
  if (trim < 1) {
    for (kk in seq_len(n_periods)) {
      tm_kk <- treatment_models_by_time[[kk]]
      data_kk <- fit_data_by_time[[kk]]
      if (is_mv) {
        tms_local <- tm_kk
        for (j in seq_along(tms_local)) {
          tms_local[[j]]$fit_rows <- rep(TRUE, nrow(data_kk))
        }
        class(tms_local) <- c("causatr_treatment_models", "list")
        if (stabilize) {
          attr(tms_local, "numerator_models") <-
            numerator_models_by_time[[kk]]
          attr(tms_local, "stabilize") <- "marginal"
        }
        w_kk <- compute_density_ratio_weights_mv(
          tms_local,
          data_kk,
          intervention,
          estimand = estimand,
          trim = 1
        )
      } else if (stabilize) {
        w_kk <- compute_stabilized_period_weight(
          tm_denom = tm_kk,
          tm_num = numerator_models_by_time[[kk]],
          data = data_kk,
          intervention = intervention,
          trim = 1
        )
      } else {
        w_kk <- compute_density_ratio_weights(
          treatment_model = tm_kk,
          data = data_kk,
          intervention = intervention,
          estimand = estimand,
          trim = 1
        )
      }
      period_thresholds[[kk]] <- stats::quantile(
        w_kk,
        trim,
        names = FALSE
      )
    }
  }

  # Per-period sub-closures + alpha block lengths + alignment maps
  # back to the canonical id ordering.
  # For univariate: each period contributes one propensity block.
  # For multivariate: each period contributes K component blocks
  # (via make_weight_fn_mv()), so the total stacked alpha has
  # T×K blocks in period-major order.
  sub_fns <- vector("list", n_periods)
  block_lens <- integer(n_periods)
  alpha_blocks <- vector("list", n_periods)
  align_idx <- vector("list", n_periods)

  for (k in seq_len(n_periods)) {
    tm_k <- treatment_models_by_time[[k]]
    data_k <- fit_data_by_time[[k]]
    ids_k <- as.character(data_k[[id_col]])

    if (is_mv) {
      # Multivariate: delegate to make_weight_fn_mv() which handles the
      # K-component stacking internally. The returned closure takes a
      # stacked alpha of length sum_k(p_k) and returns a per-row weight.
      # Reset fit_rows to all-TRUE (period-k data is already filtered).
      tms_local <- tm_k
      for (j in seq_along(tms_local)) {
        tms_local[[j]]$fit_rows <- rep(TRUE, nrow(data_k))
      }
      class(tms_local) <- c("causatr_treatment_models", "list")
      # Stabilized multivariate: re-attach the per-period numerator
      # chain + flag so `make_weight_fn_mv()` builds the fixed-numerator
      # closures (numerator gamma held fixed under numDeriv perturbation,
      # only the denominator alpha varies). Stripped by the `class<-` /
      # `[[<-` edits above, so re-attach here.
      if (stabilize) {
        attr(tms_local, "numerator_models") <- numerator_models_by_time[[k]]
        attr(tms_local, "stabilize") <- "marginal"
      }
      mv_closure_k <- make_weight_fn_mv(
        tms_local,
        data_k,
        intervention,
        estimand = estimand,
        trim = trim
      )
      sub_fn_k <- mv_closure_k$weight_fn
    } else if (stabilize) {
      # Stabilized closure: numerator density g_k is precomputed once
      # at closure-creation time (gamma fixed under numDeriv
      # perturbation; same nuisance-fixed convention as multivariate IPW).
      # Only the denominator f_k(A | ..., L; alpha) varies with alpha.
      tm_num_k <- numerator_models_by_time[[k]]
      raw_fn_k <- make_long_stabilized_period_closure(
        tm_denom = tm_k,
        tm_num = tm_num_k,
        data = data_k,
        intervention = intervention
      )
      # Wrap with per-period truncation at the precomputed threshold.
      if (trim < 1) {
        sub_fn_k <- wrap_closure_with_trim(
          raw_fn_k,
          period_thresholds[[k]]
        )
      } else {
        sub_fn_k <- raw_fn_k
      }
    } else {
      # Unstabilized univariate: reuse the existing per-period builder.
      # Per-period truncation applied here via trim + precomputed
      # threshold — matches compute_longitudinal_weights() semantics.
      sub_fn_k <- make_weight_fn(
        treatment_model = tm_k,
        data = data_k,
        intervention = intervention,
        estimand = estimand,
        trim = trim,
        trim_threshold = period_thresholds[[k]]
      )
    }

    if (is_mv) {
      period_ids <- ids_k[tm_k[[1L]]$fit_rows]
    } else {
      period_ids <- ids_k[tm_k$fit_rows]
    }
    pos <- match(period_ids, ids_first)
    align_idx[[k]] <- pos

    sub_fns[[k]] <- sub_fn_k
    if (is_mv) {
      alpha_blocks[[k]] <- mv_closure_k$alpha_hat
      block_lens[k] <- length(mv_closure_k$alpha_hat)
    } else {
      alpha_blocks[[k]] <- tm_k$alpha_hat
      block_lens[k] <- length(tm_k$alpha_hat)
    }
  }

  # offsets[k]:offsets[k+1]-1 slices the stacked alpha for period k.
  offsets <- c(1L, cumsum(block_lens) + 1L)
  alpha_hat <- unlist(alpha_blocks, use.names = FALSE)

  weight_fn <- function(alpha) {
    # Cumulative product W_i = prod_k w_{ik}(alpha_k). Each sub-closure
    # evaluates independently (period-k data is disjoint from period-k' data),
    # so numDeriv can perturb the full stacked alpha safely. Per-period
    # truncation is baked into each sub-closure (via make_weight_fn's
    # maybe_trim wrapper), so no post-product truncation is needed.
    W <- rep(1, n_id)
    for (k in seq_len(n_periods)) {
      idx <- offsets[k]:(offsets[k + 1L] - 1L)
      alpha_k <- alpha[idx]
      w_period <- sub_fns[[k]](alpha_k)

      # Project per-period weight (length n_period_k) onto the
      # canonical id ordering. Periods missing an id contribute
      # weight 0 (drops the id from the Hajek mean). Under the
      # complete-case assumption, alignment is the identity.
      # `align_idx[[k]]` was precomputed at closure-creation time to
      # avoid re-calling `match()` on every numDeriv perturbation.
      w_aligned <- rep(0, n_id)
      w_aligned[align_idx[[k]]] <- w_period
      W <- W * w_aligned
    }
    W
  }

  list(
    weight_fn = weight_fn,
    offsets = offsets,
    alpha_hat = alpha_hat
  )
}


#' Per-period stabilized closure for longitudinal IPW variance
#'
#' @description
#' Wrapper around `mv_stabilized_closure()` for one period of
#' longitudinal IPW. Builds the per-period inputs (precomputed
#' numerator density, indicator / Jacobian) from `tm_num`,
#' `tm_denom`, and `intervention`, then delegates to
#' `mv_stabilized_closure()` -- both helpers compute the same
#' "fixed numerator / alpha-dependent denominator" weight closure.
#'
#' Reusing `mv_stabilized_closure()` keeps a single
#' implementation of the Bernoulli / Gaussian / categorical / count
#' branches; the only difference between the multivariate and
#' longitudinal stabilized paths is what the per-component / per-period
#' upstream conditioning is (component-prior vs time-prior), and that's
#' baked into `tm_denom$X_prop` / `tm_num` at fit time.
#'
#' @param tm_denom Per-period denominator density model (full-L `f_k`).
#' @param tm_num Per-period numerator density model (`g_k`,
#'   no-L conditioning).
#' @param data data.table for the period-k subset.
#' @param intervention `causatr_intervention` object or `NULL`.
#' @return `function(alpha)` returning a length `nrow(data)` weight
#'   vector, with `alpha` indexing the denominator's flattened
#'   coefficient vector.
#' @noRd
make_long_stabilized_period_closure <- function(
  tm_denom,
  tm_num,
  data,
  intervention
) {
  fam_tag <- tm_denom$family
  X_prop <- tm_denom$X_prop
  trt_col <- tm_denom$treatment
  fit_rows <- tm_denom$fit_rows
  fit_data <- data[fit_rows]
  a_obs <- fit_data[[trt_col]]
  n_rows <- length(a_obs)

  # Per-intervention build-up of `f_num_fixed` (the precomputed
  # numerator density vector) and `ind_or_jac` (the
  # indicator-or-Jacobian multiplier). Mirrors the construction in
  # `make_weight_fn_mv()`'s stabilized branch.
  # `f_num_fixed` is fixed at closure-creation time; it never varies
  # with alpha during the numDeriv perturbation loop.
  if (is.null(intervention)) {
    # Natural course stabilized: numerator is g_k(A_obs | ..., V), no Jac.
    f_num_fixed <- evaluate_density(tm_num, a_obs, fit_data)
    ind_or_jac <- rep(1, n_rows)
  } else {
    check_intervention_family_compat(intervention, tm_denom)
    iv_type <- intervention$type
    if (iv_type == "ipsi") {
      # IPSI + `stabilize = "marginal"` is rejected upstream
      # (`causatr_longitudinal_ipsi_stabilize`). Unstabilized IPSI never
      # reaches this stabilized closure — its variance closure is built
      # by `make_weight_fn()`'s IPSI branch.
      rlang::abort(
        "Internal error: stabilized longitudinal closure reached IPSI branch."
      )
    }
    if (iv_type %in% c("static", "dynamic")) {
      target <- apply_intervention_to_values(intervention, fit_data, a_obs)
      ind_or_jac <- as.numeric(a_obs == target)
      # Same surviving-rows shortcut as `mv_stabilized_closure()`'s HT
      # branch: at rows where ind = 1, target = a_obs, so
      # g_k(target | ...) = g_k(a_obs | ...).
      f_num_fixed <- evaluate_density(tm_num, a_obs, fit_data)
    } else if (iv_type == "shift") {
      delta <- intervention$delta
      # a_eval = d^{-1}(A_obs) = A_obs - delta; |Jac| = 1.
      a_eval <- a_obs - delta
      f_num_fixed <- evaluate_density(tm_num, a_eval, fit_data)
      ind_or_jac <- rep(1, n_rows)
    } else if (iv_type == "scale") {
      fct <- intervention$factor
      if (fct == 0) {
        rlang::abort(
          "`scale_by(0)` collapses the treatment support; not a valid MTP."
        )
      }
      # a_eval = A_obs / c; |Jac d^{-1}| = 1/|c| baked into ind_or_jac.
      a_eval <- a_obs / fct
      f_num_fixed <- evaluate_density(tm_num, a_eval, fit_data)
      ind_or_jac <- rep(abs(1 / fct), n_rows)
    } else {
      rlang::abort(
        paste0(
          "Internal error: stabilized longitudinal closure has no branch for iv_type='",
          iv_type,
          "'."
        )
      )
    }
  }

  mv_stabilized_closure(
    fam_tag = fam_tag,
    X_prop = X_prop,
    a_obs_k = a_obs,
    f_num_fixed = f_num_fixed,
    ind_or_jac = ind_or_jac,
    sigma = tm_denom$sigma,
    theta = tm_denom$theta,
    trt_levels = tm_denom$levels
  )
}
