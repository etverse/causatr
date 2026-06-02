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
#' When `fit$details$stratified` is set, the per-step model is a list of
#' per-stratum models and the assembly is delegated to
#' `variance_if_ice_one_stratified()`, which runs the same Channel-2 chain
#' once per stratum (block-diagonal across strata; see
#' `R/variance_if_ice_stratified.R`).
#'
#' @param fit A `causatr_fit` of type `"longitudinal"`.
#' @param ice_results Named list of `ice_iterate()` results, one per
#'   intervention.
#' @param target_within_first Logical vector over first-time individuals
#'   flagging the target population.
#'
#' @noRd
variance_if_ice <- function(
  fit,
  ice_results,
  target_within_first,
  cluster_vec = NULL
) {
  int_names <- names(ice_results)
  data <- fit$data
  first_time <- fit$details$time_points[1]
  rows_first <- data[[fit$time]] == first_time
  n <- sum(rows_first)

  # Per-individual IFs in the ICE branch are indexed by unique id at
  # the first time point (one IF per person). `cluster_vec` is supplied
  # at full person-period length `nrow(fit$data)`; subset to first-time
  # rows and align by id order so it matches the IF vector ordering
  # produced by `variance_if_ice_one()`. If multiple person-period rows
  # for the same id disagree on cluster, the first-time row wins — the
  # cluster is expected to be an id-level attribute anyway.
  cluster_first <- if (is.null(cluster_vec)) {
    NULL
  } else {
    cluster_vec[rows_first]
  }

  # Pooled vs stratified ICE select different per-intervention IF
  # assemblers. Both return a length-n per-individual IF; the downstream
  # `vcov_from_if()` aggregation is identical.
  is_stratified <- !is.null(fit$details$stratified)
  IF_list <- lapply(ice_results, function(res) {
    if (is_stratified) {
      variance_if_ice_one_stratified(fit, res, target_within_first)
    } else {
      variance_if_ice_one(fit, res, target_within_first)
    }
  })

  vcov_from_if(IF_list, n, int_names, cluster = cluster_first)
}


#' Per-individual IF for one (pooled) ICE intervention
#'
#' @description
#' Workhorse called once per intervention by `variance_if_ice()` for the
#' pooled (non-stratified) case. Assembles the IF as
#' \deqn{\mathrm{IF}_i = \frac{n}{n_t}\,t_i\,(\tilde Y_{0,i} - \hat\mu)
#'       + \sum_{k=0}^{K} n_k\,d_{k,i}\,r^{\mathrm{score}}_{k,i}}
#' from a Channel-1 sampling term (`ice_if_setup()`) and a Channel-2
#' nuisance-correction chain (`variance_if_ice_chain()`). The chain loops
#' over the \eqn{K{+}1} outcome models; at each step the previous step's
#' per-individual sensitivity `d_{k-1,j}` weights the construction of
#' `g_k`, then `correct_model(model_k, g_k, ...)` returns both the model-k
#' correction and the new `d_k` for the next iteration.
#'
#' Mirrors the closed-form back-substitution of the stacked-EE bread
#' from vignette Section 5.4-5.6.
#'
#' @param fit A `causatr_fit` object (ICE estimator).
#' @param ice_result Per-intervention ICE result from `ice_iterate()`.
#' @param target Logical vector identifying target-population rows.
#'
#' @return Numeric vector of length `n` (individuals), the per-individual IF.
#'
#' @noRd
variance_if_ice_one <- function(fit, ice_result, target) {
  ctx <- ice_if_setup(fit, ice_result, target)
  # Pooled path: one Channel-2 chain over all rows / all models, with no
  # stratum restriction. Channel 1 is already in `ctx$IF_vec`.
  ctx$IF_vec +
    variance_if_ice_chain(ctx, ice_result$models, ice_result$fit_ids)
}


#' Assemble the Channel-1 IF and shared state for an ICE intervention
#'
#' @description
#' Splits the intervention-independent bookkeeping out of
#' `variance_if_ice_one()` so the pooled assembler and the stratified
#' assembler (`variance_if_ice_one_stratified()`) share identical Channel-1
#' construction, mean estimation, and weight lookups. The returned context
#' carries everything `variance_if_ice_chain()` needs.
#'
#' Channel 1 is the sampling term
#' \eqn{n\,(w_i / \sum w)\,(\tilde Y_{0,i} - \hat\mu)} over target rows. It
#' is identical for pooled and stratified fits: the estimand
#' \eqn{\hat\mu = E[Y^{\bar d}]} is marginal over all strata, so the
#' centring mean and the target weighting never partition by stratum.
#'
#' @param fit A `causatr_fit` object (ICE estimator).
#' @param ice_result Per-intervention ICE result from `ice_iterate()`.
#' @param target Logical vector identifying target-population rows.
#'
#' @return A list (context) with the Channel-1 `IF_vec` plus all shared
#'   state consumed by `variance_if_ice_chain()`: `fit`, `data`,
#'   `data_iv`, `intervention`, `time_points`, `n_times`, `id_col`,
#'   `time_col`, `all_ids`, `n`, `id_to_idx`, `target`, `w_t`,
#'   `sum_w_target`, `has_weights`, `w_at_step`.
#'
#' @noRd
ice_if_setup <- function(fit, ice_result, target) {
  data <- fit$data
  details <- fit$details
  data_iv <- ice_result$data_iv
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
  # \hat{mu} = (1/sum_w) sum_{target} w_i \tilde{Y}_{0,i}: Hajek-weighted
  # average of the backward-iterated pseudo-outcome over the target population.
  mu_hat <- sum(w_t * pseudo_final[target]) / sum_w_target
  # Channel 1: IF_i = n * (w_i / sum_w) * (\tilde{Y}_{0,i} - \hat{mu}).
  # Scaled by n so vcov_from_if(n) gives the sandwich (1/n^2) sum IF_i^2.
  IF_vec <- numeric(n)
  IF_vec[target] <- n *
    (w_t / sum_w_target) *
    (pseudo_final[target] - mu_hat)

  # Per-time-step id-to-weight lookup for the cascade gradient.
  # At step k > 1, g_k = sum_{j: prev fit} w_{k-1,j} d_{k-1,j} X*_{k,j} \mu'_{k,j}.
  # The w_{k-1} factor arises from A_{k-1,k} = (1/n) sum \partial s_{k-1}/\partial \beta_k
  # in the block-triangular bread (see variance-theory vignette Section 5.4).
  # Omitting it silently under-scales the off-diagonal correction by ext_w,
  # causing ~2x SE underestimation under heterogeneous external weights.
  w_at_step <- if (has_weights) {
    lapply(seq_len(n_times), function(k) {
      rows_k <- data[[time_col]] == time_points[k]
      ids_k <- as.character(data[rows_k][[id_col]])
      stats::setNames(ext_w[rows_k], ids_k)
    })
  } else {
    NULL
  }

  list(
    fit = fit,
    data = data,
    data_iv = data_iv,
    intervention = ice_result$intervention,
    time_points = time_points,
    n_times = n_times,
    id_col = id_col,
    time_col = time_col,
    all_ids = all_ids,
    n = n,
    id_to_idx = id_to_idx,
    target = target,
    w_t = w_t,
    sum_w_target = sum_w_target,
    has_weights = has_weights,
    w_at_step = w_at_step,
    IF_vec = IF_vec
  )
}


#' Channel-2 nuisance-correction chain for an ICE model sequence
#'
#' @description
#' Runs the backward per-step model-correction cascade for one ordered
#' sequence of outcome models, returning the accumulated length-`n`
#' Channel-2 IF (Channel 1 is added by the caller). At step `k > 1` the
#' previous step's per-individual sensitivity `d_{k-1,j}` weights the
#' construction of the gradient `g_k`; `correct_model()` then returns both
#' the model-k correction (added to the IF) and the updated `d_k` (fed to
#' step `k-1`). The sensitivity vector `d_vec` is local to the call, so a
#' caller running several disjoint chains (one per stratum) keeps their
#' cascades independent.
#'
#' `restrict` selects a row subset for stratified ICE. When non-`NULL` the
#' current-time prediction frame is filtered to rows whose `restrict$col`
#' equals `restrict$val`; the supplied `fit_ids_by_step` must already be
#' the stratum's fit ids. With `restrict = NULL` the chain runs over all
#' rows (the pooled path), behaviourally identical to the historical
#' single-function implementation.
#'
#' @param ctx Context list from `ice_if_setup()`.
#' @param models_by_step List of fitted models indexed by time step (one
#'   pooled model per step, or a single stratum's per-step models).
#' @param fit_ids_by_step List of character id vectors indexed by time
#'   step -- the individuals in each step's fitting set (already restricted
#'   to the stratum for the stratified path).
#' @param restrict `NULL` (pooled) or `list(col = <stratum column>,
#'   val = <stratum value as character>)` (stratified).
#'
#' @return Numeric vector of length `ctx$n`, the accumulated Channel-2 IF.
#'
#' @noRd
variance_if_ice_chain <- function(
  ctx,
  models_by_step,
  fit_ids_by_step,
  restrict = NULL
) {
  fit <- ctx$fit
  data <- ctx$data
  data_iv <- ctx$data_iv
  intervention <- ctx$intervention
  time_points <- ctx$time_points
  n_times <- ctx$n_times
  id_col <- ctx$id_col
  time_col <- ctx$time_col
  all_ids <- ctx$all_ids
  n <- ctx$n
  id_to_idx <- ctx$id_to_idx
  target <- ctx$target
  w_t <- ctx$w_t
  sum_w_target <- ctx$sum_w_target
  has_weights <- ctx$has_weights
  w_at_step <- ctx$w_at_step

  # Stratum membership predicate over an arbitrary person-period frame.
  # For the pooled path every row passes; for the stratified path only
  # rows in the active stratum pass (the baseline stratum column is
  # carried verbatim into `data_iv` and the re-drawn stochastic frames).
  strat_ok <- function(frame) {
    if (is.null(restrict)) {
      rep(TRUE, nrow(frame))
    } else {
      as.character(frame[[restrict$col]]) == restrict$val
    }
  }

  is_stoch <- has_stochastic_component(intervention)
  IF_acc <- numeric(n)
  d_vec <- rep(0, n)

  for (step_i in seq_len(n_times)) {
    model_k <- models_by_step[[step_i]]
    if (is.null(model_k)) {
      next
    }

    current_time <- time_points[step_i]
    fit_ids_k <- fit_ids_by_step[[step_i]]
    if (length(fit_ids_k) == 0L) {
      next
    }

    # Subset the (already intervention-modified) longitudinal data to
    # the rows at the current time step (and, when stratified, the active
    # stratum), and record which individuals are present (some may have
    # been censored out before this step).
    rows_iv_current <- data_iv[[time_col]] == current_time
    iv_data_current <- data_iv[rows_iv_current]
    iv_data_current <- iv_data_current[strat_ok(iv_data_current)]
    iv_ids_current <- as.character(iv_data_current[[id_col]])

    if (step_i == 1L) {
      # Step 1 (earliest time): g_1 = (1/sum_w) sum_{target} w_i X*_{1,i} \mu'(\eta*_{1,i}).
      # This is \partial \hat{mu} / \partial \beta_1, the gradient of the ICE
      # estimand w.r.t. the first outcome model's parameters. For stratified
      # ICE the iv frame is already restricted to the stratum, so the sum
      # naturally covers only target rows in this stratum -- exactly
      # \partial \hat{mu} / \partial \beta_{g,1} -- while the normaliser
      # sum_w_target stays global (the 1/n_t factor of the marginal mean).
      target_in_iv <- match(all_ids[target], iv_ids_current)
      valid_target <- !is.na(target_in_iv)
      target_in_iv <- target_in_iv[valid_target]
      target_w <- w_t[valid_target]

      if (is_stoch) {
        n_mc_ice <- get_stochastic_n_mc(intervention)
        p_coef <- length(coef_clean(model_k))
        g_k_sum <- numeric(p_coef)
        for (m in seq_len(n_mc_ice)) {
          dv_m <- ice_apply_intervention_long(
            data,
            fit$treatment,
            intervention,
            id_col,
            time_col
          )
          iv_cur_m <- dv_m[dv_m[[time_col]] == current_time]
          iv_cur_m <- iv_cur_m[strat_ok(iv_cur_m)]
          X_m <- iv_design_matrix(model_k, iv_cur_m)
          eta_m <- as.numeric(X_m %*% coef_clean(model_k))
          mu_eta_m <- model_k$family$mu.eta(eta_m)
          ids_m <- as.character(iv_cur_m[[id_col]])
          tgt_in_m <- match(all_ids[target], ids_m)
          vt_m <- !is.na(tgt_in_m)
          tgt_in_m <- tgt_in_m[vt_m]
          # target_w is length sum(valid_target); vt_m is length
          # sum(target). Index into the already-subsetted target_w
          # via vt_m[valid_target] so both operands align.
          vt_in_valid <- vt_m[valid_target]
          g_k_sum <- g_k_sum +
            as.numeric(
              crossprod(
                X_m[tgt_in_m, , drop = FALSE],
                target_w[vt_in_valid] * mu_eta_m[tgt_in_m]
              )
            )
        }
        g_k <- g_k_sum / (n_mc_ice * sum_w_target)
      } else {
        X_star_k <- iv_design_matrix(model_k, iv_data_current)
        eta_star <- as.numeric(
          X_star_k %*% coef_clean(model_k)
        )
        mu_eta_star <- model_k$family$mu.eta(eta_star)
        # g_k = (1/sum_w) X*^T diag(w_target) \mu'(\eta*): the raw gradient
        # before apply_model_correction scales it through the model bread.
        g_k <- as.numeric(
          crossprod(
            X_star_k[target_in_iv, , drop = FALSE],
            target_w * mu_eta_star[target_in_iv]
          )
        ) /
          sum_w_target
      }
    } else {
      # Step k > 1: chain rule through the ICE recursion.
      # g_k = sum_{j: prev fit} w_{k-1,j} d_{k-1,j} X*_{k,j} \mu'_{k,j}.
      # d_{k-1,j} is the sensitivity of \hat{mu} to Q_{k-1,j}; it comes
      # from the previous correct_model() call. w_{k-1,j} is the external
      # weight at the previous time step (see w_at_step lookup above).
      prev_fit_ids <- fit_ids_by_step[[step_i - 1L]]
      idx_in_all <- id_to_idx[prev_fit_ids]
      rows_in_iv <- match(prev_fit_ids, iv_ids_current)
      keep <- !is.na(idx_in_all) & !is.na(rows_in_iv)

      if (any(keep)) {
        d_prev <- d_vec[idx_in_all[keep]]
        w_prev <- if (has_weights) {
          unname(w_at_step[[step_i - 1L]][prev_fit_ids[keep]])
        } else {
          rep(1, sum(keep))
        }

        if (is_stoch) {
          n_mc_ice <- get_stochastic_n_mc(intervention)
          p_coef <- length(coef_clean(model_k))
          g_k_sum <- numeric(p_coef)
          for (m in seq_len(n_mc_ice)) {
            dv_m <- ice_apply_intervention_long(
              data,
              fit$treatment,
              intervention,
              id_col,
              time_col
            )
            iv_cur_m <- dv_m[dv_m[[time_col]] == current_time]
            iv_cur_m <- iv_cur_m[strat_ok(iv_cur_m)]
            X_m <- iv_design_matrix(model_k, iv_cur_m)
            eta_m <- as.numeric(X_m %*% coef_clean(model_k))
            mu_eta_m <- model_k$family$mu.eta(eta_m)
            ids_m <- as.character(iv_cur_m[[id_col]])
            riv_m <- match(prev_fit_ids, ids_m)
            keep_m <- !is.na(idx_in_all) & !is.na(riv_m)
            # keep_m_in_keep: which of the already-kept indices also
            # survive in keep_m (both are length P, but w_prev/d_prev
            # are length sum(keep), so index into them via keep[keep_m]).
            keep_m_in_keep <- keep_m[keep]
            if (any(keep_m)) {
              wg_m <- w_prev[keep_m_in_keep] *
                d_prev[keep_m_in_keep] *
                mu_eta_m[riv_m[keep_m]]
              g_k_sum <- g_k_sum +
                as.numeric(
                  crossprod(
                    X_m[riv_m[keep_m], , drop = FALSE],
                    wg_m
                  )
                )
            }
          }
          g_k <- g_k_sum / n_mc_ice
        } else {
          X_star_k <- iv_design_matrix(model_k, iv_data_current)
          eta_star <- as.numeric(
            X_star_k %*% coef_clean(model_k)
          )
          mu_eta_star <- model_k$family$mu.eta(eta_star)
          weights_g <- w_prev * d_prev * mu_eta_star[rows_in_iv[keep]]
          g_k <- as.numeric(
            crossprod(
              X_star_k[rows_in_iv[keep], , drop = FALSE],
              weights_g
            )
          )
        }
      } else {
        p_coef <- length(coef_clean(model_k))
        g_k <- rep(0, p_coef)
      }
    }

    fit_id_idx <- id_to_idx[fit_ids_k]
    na_act_k <- model_k$na.action
    if (!is.null(na_act_k)) {
      fit_id_idx <- fit_id_idx[-na_act_k]
    }
    # correct_model() computes two outputs:
    #   $correction: per-individual model-k IF contribution (n-scaled)
    #   $d:          updated sensitivity vector d_k = d_{k-1} * (\partial Q_{k-1}/\partial Q_k)
    # The correction is added to the accumulator; d_vec is passed to the
    # next step.
    res <- correct_model(model_k, g_k, fit_id_idx, n)
    IF_acc <- IF_acc + res$correction
    d_vec <- res$d
  }

  IF_acc
}
