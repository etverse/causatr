#' Compute cumulative density-ratio weights for one longitudinal intervention
#'
#' @description
#' Longitudinal companion to `compute_density_ratio_weights()` and
#' `compute_density_ratio_weights_mv()`. For a fixed intervention `d`
#' applied at every period `k = 1..K`, builds the per-individual joint
#' weight
#' \deqn{W_i = \prod_{k=1}^{K} w_k(A_{k,i}, L_{k,i})
#'      = \prod_{k=1}^{K} \frac{f_k(d(A_{k,i}, L_{k,i}) \mid \bar A_{k-1,i}^{\mathrm{obs}}, \bar L_{k,i})}
#'                              {f_k(A_{k,i} \mid \bar A_{k-1,i}^{\mathrm{obs}}, \bar L_{k,i})}}
#' under the sequential-MTP semantics of Robins, Hernán, Brumback (2000) and
#' Díaz et al. (2023). The k-th per-period factor is computed by
#' reusing `compute_density_ratio_weights()` against the period-k
#' subset `data_at_time = fit_data_by_time[[k]]` and the period-k
#' density model `tm_k`. The product is taken **across time** within
#' each individual.
#'
#' Issues a sequential-positivity warning when any per-period weight
#' has a heavy right tail (`max > positivity_threshold`); the per-period
#' breakdown helps the user identify which period is unstable, which
#' the cumulative product alone would obscure.
#'
#' Natural course (NULL intervention) collapses to a constant 1 vector
#' -- every per-period factor equals 1 because numerator and
#' denominator densities coincide.
#'
#' Under `stabilize = "marginal"` the per-period numerator density is
#' swapped for `g_k(d(A_k) | A_{1..k-1}, V)` from a separately fit
#' numerator model that drops time-varying confounders `L_k` (with
#' optional baseline conditioning `V` if `numerator =` was supplied).
#' The denominator stays at the full-L density `f_k`. Numerator
#' parameters are held fixed in the variance engine; bootstrap refits
#' both numerator and denominator and captures the full uncertainty.
#'
#' @param treatment_models_by_time Named list of T per-period models
#'   from `fit_longitudinal_ipw()`. Each element is either a scalar
#'   `causatr_treatment_model` (univariate) or a
#'   `causatr_treatment_models` list of K components (multivariate).
#' @param fit_data_by_time Named list of T per-period data subsets
#'   (one row per individual at each period).
#' @param ids_first Character vector of individual ids in the order
#'   the returned weight vector is indexed by.
#' @param id_col Character. Name of the id column.
#' @param intervention A `causatr_intervention` object, a named list
#'   of per-component interventions (multivariate), or `NULL`.
#' @param estimand Character. ATE only for longitudinal IPW.
#' @param positivity_threshold Positive scalar. Per-period max weight
#'   above this triggers the sequential-positivity warning.
#' @param numerator_models_by_time Optional named list of K per-period
#'   numerator density models (from
#'   `fit_longitudinal_ipw(stabilize = "marginal")`). When `NULL`
#'   (the default) the routine returns unstabilized weights identical
#'   to `stabilize = "none"`.
#'
#' @return Numeric vector of length `length(ids_first)` -- the
#'   cumulative product weight per individual, ordered by `ids_first`.
#'
#' @noRd
compute_longitudinal_weights <- function(
  treatment_models_by_time,
  fit_data_by_time,
  ids_first,
  id_col,
  intervention,
  estimand = "ATE",
  positivity_threshold = 100,
  numerator_models_by_time = NULL,
  trim = 1
) {
  if (estimand != "ATE") {
    rlang::abort(
      "Longitudinal IPW only supports estimand = 'ATE'.",
      class = "causatr_longitudinal_estimand"
    )
  }

  n_times <- length(treatment_models_by_time)
  n_id <- length(ids_first)

  # Detect multivariate: per-period models are `causatr_treatment_models`
  # (plural, from `fit_treatment_models()`) rather than scalar
  # `causatr_treatment_model` objects.
  is_mv <- inherits(
    treatment_models_by_time[[1L]],
    "causatr_treatment_models"
  )

  # Natural course: w = 1 for every individual at every period; the
  # product is identically 1.
  if (is.null(intervention)) {
    return(rep(1, n_id))
  }

  # Per-period weights stacked into an `n_id × T` matrix so we can
  # diagnose positivity per period and then take row-products to get
  # the per-id cumulative weight. Each column is built by reusing the
  # univariate / multivariate density-ratio engine against the
  # period-k model and subset.
  # W[i, k] = w_{ik}; cumulative weight W_i = prod_k W[i, k] (row-product).
  W_per_period <- matrix(NA_real_, nrow = n_id, ncol = n_times)

  stabilize <- !is.null(numerator_models_by_time)

  for (k in seq_len(n_times)) {
    tm_k <- treatment_models_by_time[[k]]
    data_k <- fit_data_by_time[[k]]
    ids_k <- as.character(data_k[[id_col]])

    if (stabilize) {
      tm_num_k <- numerator_models_by_time[[k]]
      w_k <- compute_stabilized_period_weight(
        tm_denom = tm_k,
        tm_num = tm_num_k,
        data = data_k,
        intervention = intervention,
        trim = trim
      )
    } else if (is_mv) {
      # Multivariate: per-period model is a list of K component models.
      # Reset fit_rows to all-TRUE (period-k subset is already filtered),
      # matching the convention in compute_ipw_contrast_point().
      tms_local <- tm_k
      for (j in seq_along(tms_local)) {
        tms_local[[j]]$fit_rows <- rep(TRUE, nrow(data_k))
      }
      class(tms_local) <- c("causatr_treatment_models", "list")
      w_k <- compute_density_ratio_weights_mv(
        tms_local,
        data_k,
        intervention,
        estimand = estimand,
        trim = trim
      )
    } else {
      w_k <- compute_density_ratio_weights(
        treatment_model = tm_k,
        data = data_k,
        intervention = intervention,
        estimand = estimand,
        trim = trim
      )
    }

    # Align per-period weight to the canonical id ordering. Some ids
    # may have been dropped at period k by `fit_treatment_model()`'s
    # NA mask; we expect a complete-case design so this
    # is identity. Defensive: rows that weren't fit at period k get
    # weight 0 (absorbed into the cumulative product, drops them
    # from the Hajek mean).
    # `match()` returns NA for ids absent from ids_first; the
    # assignment `w_aligned[NA] <- ...` is silently ignored, so
    # unmatched ids stay 0 -- correct.
    #
    # For MV, `fit_rows` alignment uses the first component's
    # `fit_rows` (all components share the same row mask after the
    # reset above).
    if (is_mv) {
      period_ids <- ids_k[tm_k[[1L]]$fit_rows]
    } else {
      period_ids <- ids_k[tm_k$fit_rows]
    }
    w_aligned <- rep(0, n_id)
    pos <- match(period_ids, ids_first)
    w_aligned[pos] <- w_k
    W_per_period[, k] <- w_aligned
  }

  # Sequential positivity warning. We surface the per-period max weight
  # (and the time index) so users can spot which period is the
  # bottleneck. The default threshold of 100 mirrors a common rule of
  # thumb (`cobalt::col_w_dist`) for IPW weights; a single per-period
  # excursion above this is enough to inflate the cumulative product
  # noticeably.
  warn_seq_positivity(
    W_per_period,
    time_points = names(treatment_models_by_time),
    threshold = positivity_threshold
  )

  # Row-product across time. `apply(W, 1, prod)` is fine here because
  # K is typically small (< 10) and the matrix is dense.
  apply(W_per_period, 1L, prod)
}


#' Compute one period's stabilized density-ratio weight
#'
#' @description
#' Helper for `compute_longitudinal_weights()` under
#' `stabilize = "marginal"`. The per-period weight is
#' \deqn{w_k = \frac{g_k(d(A_{k}, L_k) \mid A_{1..k-1}^{\mathrm{obs}}, V) \cdot |\mathrm{Jac}\,d^{-1}|}{f_k(A_k \mid A_{1..k-1}^{\mathrm{obs}}, \bar L_k)},}
#' where the **numerator** density `g_k` (no L conditioning) replaces
#' the denominator `f_k` (full-L conditioning) used in the
#' unstabilized form. This dampens the multiplicative L-dependence
#' across K factors and typically reduces the cumulative product's
#' variance.
#'
#' Implementation reuses the existing per-family branches of
#' `compute_density_ratio_weights()` for the **denominator**
#' (which still depends on L via the full-L model `f_k`) and for the
#' **shape of the numerator evaluation**: HT indicator on discrete,
#' pushforward on continuous / count. Only the density-evaluation
#' model is swapped from `tm_denom` to `tm_num`.
#'
#' @param tm_denom Per-period denominator density model (full-L `f_k`).
#' @param tm_num Per-period numerator density model (`g_k`,
#'   no-L conditioning under default; user-customised via
#'   `numerator =`).
#' @param data data.table for the period-k subset.
#' @param intervention `causatr_intervention` object or `NULL`
#'   (natural course; unstabilized branch -> ratio is 1, stabilized
#'   -> ratio is `g_k(A_k | ...) / f_k(A_k | ..., L)`).
#' @return Numeric weight vector, length `sum(tm_denom$fit_rows)`.
#' @noRd
compute_stabilized_period_weight <- function(
  tm_denom,
  tm_num,
  data,
  intervention,
  trim = 1
) {
  fit_rows <- tm_denom$fit_rows
  fit_data <- data[fit_rows]
  trt_col <- tm_denom$treatment
  a_obs <- fit_data[[trt_col]]
  family_tag <- tm_denom$family

  # Denominator: f_k(A_obs | A_{1..k-1}_obs, L_obs).
  f_obs <- evaluate_density(tm_denom, a_obs, fit_data)
  check_density_positivity(
    f_obs,
    "stabilized longitudinal weight (denominator)"
  )

  if (is.null(intervention)) {
    # Stabilized natural course: g_k(A_obs | ...) / f_k(A_obs | ..., L).
    f_num <- evaluate_density(tm_num, a_obs, fit_data)
    return(truncate_weights(f_num / f_obs, trim))
  }

  check_intervention_family_compat(intervention, tm_denom)
  iv_type <- intervention$type

  if (iv_type == "ipsi") {
    # IPSI under longitudinal IPW is rejected upstream
    # (`causatr_longitudinal_ipsi_pending`); reaching this branch
    # means the rejection was bypassed.
    rlang::abort(
      "Internal error: stabilized longitudinal weight reached IPSI branch."
    )
  }

  # Discrete HT branch (static / dynamic on Bernoulli / categorical).
  # Numerator evaluates `g_k` at the target value; surviving rows
  # (where the indicator is 1) have observed A == target, so the
  # numerator density at those rows equals `g_k(A_obs | ...)` -- no
  # need to substitute the target into newdata. Multivariate IPW's
  # `mv_stabilized_closure()` HT branch uses the same shortcut.
  # w_k = I(A_k = target) * g_k(A_k | ...) / f_k(A_k | ..., L)
  is_ht <- iv_type %in%
    c("static", "dynamic") &&
    family_tag %in% c("bernoulli", "categorical")
  if (is_ht) {
    target <- apply_intervention_to_values(intervention, fit_data, a_obs)
    ind <- as.numeric(a_obs == target)
    f_num <- evaluate_density(tm_num, a_obs, fit_data)
    return(truncate_weights(ind * f_num / f_obs, trim))
  }

  # Smooth pushforward branch (shift / scale_by on Gaussian / count).
  if (iv_type == "shift") {
    delta <- intervention$delta
    # a_eval = d^{-1}(A_obs) = A_obs - delta; numerator evaluated there.
    a_eval <- a_obs - delta
    f_num <- evaluate_density(tm_num, a_eval, fit_data)
    warn_intervened_density_near_zero(
      f_num,
      "stabilized longitudinal shift"
    )
    return(truncate_weights(f_num / f_obs, trim))
  }
  if (iv_type == "scale") {
    fct <- intervention$factor
    if (fct == 0) {
      rlang::abort(
        "`scale_by(0)` collapses the treatment support; not a valid MTP."
      )
    }
    # a_eval = A_obs / c; divide by |c| for the inverse-map Jacobian.
    a_eval <- a_obs / fct
    f_num <- evaluate_density(tm_num, a_eval, fit_data)
    warn_intervened_density_near_zero(
      f_num,
      "stabilized longitudinal scale"
    )
    return(truncate_weights((f_num / f_obs) / abs(fct), trim))
  }

  rlang::abort(
    paste0(
      "Internal error: stabilized longitudinal weight has no branch for iv_type='",
      iv_type,
      "'."
    )
  )
}


#' Warn when any per-period weight has a heavy right tail
#'
#' @description
#' Sequential positivity check for longitudinal IPW. Issues a single
#' rate-limited `rlang::warn()` listing every period whose maximum
#' weight exceeds `threshold`, along with the per-period max weight
#' (sortable by user). No-op when every per-period weight is below
#' threshold.
#'
#' Detection lives here rather than in `diagnose()` because the
#' per-period decomposition is not preserved past
#' `compute_longitudinal_weights()` -- once the product is taken,
#' a single bad period is indistinguishable from a uniform mild
#' inflation across periods. The full diagnostic integration is
#' planned for the `diagnose()` rewrite.
#'
#' @param W_per_period Numeric matrix `n_id × K` of per-period
#'   weights.
#' @param time_points Character vector (length K) labelling the time
#'   points -- typically `as.character(time_points)` from the fit.
#' @param threshold Positive scalar. Per-period max above which a
#'   warning fires.
#'
#' @return `NULL` invisibly.
#'
#' @noRd
warn_seq_positivity <- function(W_per_period, time_points, threshold = 100) {
  per_period_max <- apply(W_per_period, 2L, max, na.rm = TRUE)
  bad_idx <- which(per_period_max > threshold)
  if (length(bad_idx) == 0L) {
    return(invisible(NULL))
  }

  bad_summary <- vapply(
    bad_idx,
    function(k) {
      paste0(
        "time = ",
        time_points[k],
        " (max weight ",
        format(per_period_max[k], digits = 3),
        ")"
      )
    },
    character(1)
  )

  rlang::warn(
    c(
      paste0(
        "Sequential positivity warning: ",
        length(bad_idx),
        " period(s) have a per-period density-ratio weight > ",
        threshold,
        "."
      ),
      i = paste0(
        "Affected periods: ",
        paste(bad_summary, collapse = "; "),
        "."
      ),
      i = paste0(
        "The cumulative product weight may be unstable. Consider a ",
        "richer propensity model (`propensity_model_fn = mgcv::gam`), ",
        "tighter intervention (smaller shift), or an MTP that stays ",
        "closer to the natural treatment distribution."
      )
    ),
    class = "causatr_longitudinal_seq_positivity",
    .frequency = "regularly",
    .frequency_id = "causatr_longitudinal_seq_positivity"
  )
  invisible(NULL)
}


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
      if (stabilize) {
        w_kk <- compute_stabilized_period_weight(
          tm_denom = tm_kk,
          tm_num = numerator_models_by_time[[kk]],
          data = data_kk,
          intervention = intervention,
          trim = 1
        )
      } else if (is_mv) {
        tms_local <- tm_kk
        for (j in seq_along(tms_local)) {
          tms_local[[j]]$fit_rows <- rep(TRUE, nrow(data_kk))
        }
        class(tms_local) <- c("causatr_treatment_models", "list")
        w_kk <- compute_density_ratio_weights_mv(
          tms_local,
          data_kk,
          intervention,
          estimand = estimand,
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

    if (stabilize) {
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
    } else if (is_mv) {
      # Multivariate: delegate to make_weight_fn_mv() which handles the
      # K-component stacking internally. The returned closure takes a
      # stacked alpha of length sum_k(p_k) and returns a per-row weight.
      # Reset fit_rows to all-TRUE (period-k data is already filtered).
      tms_local <- tm_k
      for (j in seq_along(tms_local)) {
        tms_local[[j]]$fit_rows <- rep(TRUE, nrow(data_k))
      }
      class(tms_local) <- c("causatr_treatment_models", "list")
      mv_closure_k <- make_weight_fn_mv(
        tms_local,
        data_k,
        intervention,
        estimand = estimand,
        trim = trim
      )
      sub_fn_k <- mv_closure_k$weight_fn
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
