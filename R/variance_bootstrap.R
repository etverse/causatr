#' Process bootstrap results into vcov matrix and boot_t
#'
#' Shared post-processing for `variance_bootstrap()` and
#' `ice_variance_bootstrap()`: extracts complete replicates, warns on
#' failures, and computes the variance-covariance matrix.
#'
#' @param boot_res A `boot::boot` result object.
#' @param int_names Character vector of intervention names.
#' @param n_boot Integer. Total number of requested replicates.
#' @return A list with `vcov` (k x k matrix), `boot_t` (matrix of
#'   successful replicates), and `boot_info` (a 3-element list of
#'   `n_requested`, `n_ok`, `n_fail`) so downstream consumers
#'   (`contrast()`, `print.causatr_result()`) can surface the bootstrap
#'   failure rate without re-deriving it from `boot_t`.
#' @noRd
process_boot_results <- function(boot_res, int_names, n_boot) {
  # `boot_res$t` is an (R x k) matrix of statistics. Failed replicates
  # are flagged by NA rows (each per-rep function wraps its body in
  # tryCatch and returns `rep(NA_real_, k)` on error), so
  # `complete.cases()` identifies the usable ones.
  #
  # We also track per-intervention failure counts in the returned
  # `boot_info$n_fail_by_int`. This is the only way to see whether
  # failures cluster on one intervention (e.g. a rare static value
  # that triggers separation in every resample) versus are spread
  # across the whole replicate (e.g. a bad factor-level sample): the
  # whole-row drop biases the vcov more in the former case. Downstream
  # `print.causatr_result()` surfaces the vector so users spot the
  # pattern.
  t_mat <- boot_res$t
  colnames(t_mat) <- int_names
  complete_rows <- stats::complete.cases(t_mat)
  n_ok <- sum(complete_rows)
  n_fail <- n_boot - n_ok
  # Per-intervention failure counts: for each column, count NA rows.
  # This is strictly `>= n_fail / k` (a failed replicate contributes
  # an NA to every column), but is a better diagnostic because a
  # heterogeneous failure pattern reveals which intervention caused
  # the failure.
  n_fail_by_int <- vapply(
    seq_along(int_names),
    function(j) sum(is.na(t_mat[, j])),
    integer(1)
  )
  names(n_fail_by_int) <- int_names

  if (n_ok < 2L) {
    # With fewer than 2 successful replicates the sample variance is
    # undefined. Return a NA vcov so downstream `sqrt(diag(.))` yields
    # NA SE rather than a misleading zero.
    rlang::warn(
      "Bootstrap produced fewer than 2 non-NA replicates; SE estimates will be NA."
    )
    vcov_mat <- matrix(
      NA_real_,
      nrow = length(int_names),
      ncol = length(int_names)
    )
  } else {
    if (n_fail > 0L) {
      pct <- round(100 * n_fail / n_boot, 1)
      # Two different warnings for the same problem — the 20%
      # threshold is an ad-hoc "this is getting bad" line. High
      # failure rates usually mean a fragile pipeline (e.g. factor
      # levels absent from some resamples), and the resulting SEs
      # are biased because the successful replicates aren't a
      # random sample of the sampling distribution.
      if (pct > 20) {
        rlang::warn(paste0(
          n_fail,
          " of ",
          n_boot,
          " bootstrap replicates (",
          pct,
          "%) failed. High failure rate may indicate model ",
          "instability; variance estimates may be unreliable."
        ))
      } else {
        rlang::warn(paste0(
          n_fail,
          " of ",
          n_boot,
          " bootstrap replicates (",
          pct,
          "%) failed and were discarded. ",
          "Variance is estimated from the ",
          n_ok,
          " successful replicates."
        ))
      }
    }
    # The bootstrap vcov is just the sample covariance of the
    # successful replicates (Davison & Hinkley 1997, §2.5.3).
    vcov_mat <- stats::var(t_mat[complete_rows, , drop = FALSE])
  }

  rownames(vcov_mat) <- int_names
  colnames(vcov_mat) <- int_names

  boot_t <- t_mat[complete_rows, , drop = FALSE]
  colnames(boot_t) <- int_names

  list(
    vcov = vcov_mat,
    boot_t = boot_t,
    boot_info = list(
      n_requested = n_boot,
      n_ok = n_ok,
      n_fail = n_fail,
      n_fail_by_int = n_fail_by_int
    )
  )
}

#' Bootstrap variance–covariance matrix for marginal means
#'
#' @description
#' Estimates Var(μ̂) by resampling the entire estimation pipeline `n_boot`
#' times (Hernán & Robins, Technical Point 13.1).  Works for all causal
#' estimators (g-comp, IPW, matching).
#'
#' ## Algorithm
#'
#' For each bootstrap replicate b = 1, …, B:
#' 1. Draw n individuals with replacement from the full dataset.
#' 2. Refit the full estimation pipeline on the bootstrap sample:
#'    - **g-comp**: refit the outcome model on uncensored rows.
#'    - **IPW**: re-estimate weights + refit the weighted MSM.
#'    - **matching**: re-match + refit on the matched sample.
#' 3. For each intervention, apply it and compute the marginal mean.
#' 4. Collect the k-vector of marginal means.
#'
#' The k × k bootstrap vcov is `var(boot_replicates)`.
#'
#' @param fit A `causatr_fit` object.
#' @param interventions Named list of `causatr_intervention` objects.
#' @param n_boot Positive integer. Number of bootstrap replicates.
#' @param target_idx Logical vector (length n) flagging target-population rows.
#' @param est Character. Estimand string (`"ATE"`, `"ATT"`, or `"ATC"`).
#' @param subset Quoted expression or `NULL`.
#'
#' @return A named k × k variance–covariance matrix.
#'
#' @noRd
variance_bootstrap <- function(
  fit,
  interventions,
  n_boot,
  target_idx,
  est,
  subset,
  parallel = "no",
  ncpus = 1L
) {
  data <- fit$data
  int_names <- names(interventions)
  estimator <- fit$estimator
  outcome <- fit$outcome
  treatment <- fit$treatment
  censoring <- fit$censoring

  # boot_fn: called by boot::boot() for each replicate.
  # The entire body is wrapped in tryCatch because bootstrap samples may
  # cause downstream failures (e.g. factor levels absent from uncensored
  # rows but present in prediction data).  Failed replicates return NA and
  # are excluded from the variance calculation — this is the standard
  # bootstrap approach (Davison & Hinkley, 1997, §2.5.3).
  boot_fn <- function(d, indices) {
    tryCatch(
      {
        d_b <- d[indices]
        # Resample weights alongside data rows.
        w_b <- if (!is.null(fit$details$weights)) {
          fit$details$weights[indices]
        } else {
          NULL
        }

        # Demote only the specific warnings we expect to emit on
        # nearly every replicate and that would otherwise flood the
        # console: GLM "fitted probabilities numerically 0 or 1",
        # `X'WX` near-singular on thin resamples, and MatchIt's
        # "Fewer control/treated units". Everything else — non-
        # convergence, family mismatches, factor-level surprises — is
        # allowed to surface so users can spot pipeline instability
        # instead of seeing a silent NA column in `boot_t`.
        model_b <- withCallingHandlers(
          refit_model(fit, d_b, weights = w_b),
          warning = function(w) {
            msg <- conditionMessage(w)
            if (
              grepl(
                "fitted probabilities numerically 0 or 1",
                msg,
                fixed = TRUE
              ) ||
                grepl("X'WX` is singular", msg, fixed = TRUE) ||
                grepl("Fewer (control|treated) units", msg)
            ) {
              invokeRestart("muffleWarning")
            }
          }
        )

        target_idx_b <- get_target_idx(d_b, treatment, est, subset)

        vapply(
          interventions,
          function(iv) {
            data_a_b <- apply_intervention(d_b, treatment, iv)
            pred_a_b <- predict(model_b, newdata = data_a_b, type = "response")
            valid <- target_idx_b & !is.na(pred_a_b)
            maybe_weighted_mean(pred_a_b[valid], if (!is.null(w_b)) w_b[valid])
          },
          numeric(1)
        )
      },
      error = function(e) rep(NA_real_, length(int_names))
    )
  }

  boot_res <- boot::boot(
    data = data,
    statistic = boot_fn,
    R = n_boot,
    parallel = parallel,
    ncpus = ncpus
  )

  process_boot_results(boot_res, int_names, n_boot)
}

#' Refit the full estimation pipeline on a bootstrap sample
#'
#' @description
#' Dispatches to the appropriate refitting logic based on the causal
#' estimator stored in `fit$estimator`.
#'
#' @param fit A `causatr_fit` object (original fit, for extracting formulas etc.).
#' @param d_b A data.table bootstrap sample.
#'
#' @return A fitted model object (glm, glm_weightit, etc.) or NULL on failure.
#'
#' @noRd
refit_model <- function(fit, d_b, weights = NULL) {
  if (fit$estimator == "gcomp") {
    refit_gcomp(fit, d_b, weights = weights)
  } else if (fit$estimator == "ipw") {
    refit_ipw(fit, d_b, weights = weights)
  } else if (fit$estimator == "matching") {
    refit_matching(fit, d_b, weights = weights)
  } else {
    rlang::abort(
      paste0(
        "Bootstrap not yet supported for estimator = '",
        fit$estimator,
        "'."
      ),
      .call = FALSE
    )
  }
}

#' Refit g-comp outcome model on a bootstrap sample
#'
#' @param fit A `causatr_fit` object.
#' @param d_b A data.table bootstrap sample.
#' @return A fitted model object.
#' @noRd
refit_gcomp <- function(fit, d_b, weights = NULL) {
  model_formula <- stats::formula(fit$model)
  family <- fit$model$family
  model_fn <- fit$details$model_fn
  censoring <- fit$censoring
  outcome <- fit$outcome

  fit_rows_b <- get_fit_rows(d_b, outcome, censoring)

  args <- list(model_formula, data = d_b[fit_rows_b], family = family)
  if (!is.null(weights)) {
    args$weights <- weights[fit_rows_b]
  }
  do.call(model_fn, args)
}

#' Refit IPW propensity weights and MSM on a bootstrap sample
#'
#' @param fit A `causatr_fit` object.
#' @param d_b A data.table bootstrap sample.
#' @return A `glm_weightit` model object.
#' @noRd
refit_ipw <- function(fit, d_b, weights = NULL) {
  ps_formula <- build_ps_formula(fit$confounders, fit$treatment)

  fit_rows_b <- get_fit_rows(d_b, fit$outcome)
  fit_data_b <- d_b[fit_rows_b]

  w_b <- WeightIt::weightit(
    ps_formula,
    data = fit_data_b,
    estimand = fit$estimand
  )

  # Mirror the alignment guard in `fit_ipw()` (R/ipw.R): WeightIt drops
  # rows with missing PS-formula covariates, so `w_b$weights` may be
  # shorter than `fit_rows_b`. The next-line external-weight multiply
  # would then recycle silently and misalign weights to the wrong
  # individuals. Abort the bootstrap replicate via an error that the
  # outer `tryCatch` in `boot_fn()` will catch and convert to an NA
  # row, so the failure is counted in `boot_info$n_fail` rather than
  # silently producing a wrong vcov.
  if (length(w_b$weights) != sum(fit_rows_b)) {
    rlang::abort(
      paste0(
        "Bootstrap replicate: WeightIt returned ",
        length(w_b$weights),
        " weights but causatr selected ",
        sum(fit_rows_b),
        " fitting rows. A confounder column has missing values that ",
        "WeightIt dropped; the multiply on the next line would recycle ",
        "silently. Clean or impute the data before calling `causat()`."
      ),
      .call = FALSE
    )
  }

  if (!is.null(weights)) {
    w_b$weights <- w_b$weights * weights[fit_rows_b]
  }

  msm_formula <- stats::reformulate(fit$treatment, response = fit$outcome)
  WeightIt::glm_weightit(
    msm_formula,
    data = fit_data_b,
    weightit = w_b,
    family = fit$model$family
  )
}

#' Re-match and refit outcome model on a bootstrap sample
#'
#' @param fit A `causatr_fit` object.
#' @param d_b A data.table bootstrap sample.
#' @return A `glm` model fit on the matched bootstrap data.
#' @noRd
refit_matching <- function(fit, d_b, weights = NULL) {
  ps_formula <- build_ps_formula(fit$confounders, fit$treatment)

  fit_rows_b <- get_fit_rows(d_b, fit$outcome)
  fit_data_b <- as.data.frame(d_b[fit_rows_b])

  match_args <- list(ps_formula, data = fit_data_b, estimand = fit$estimand)
  if (fit$estimand == "ATE") {
    check_pkg("optmatch")
    match_args$method <- "full"
  }
  m_b <- do.call(MatchIt::matchit, match_args)
  matched_b <- MatchIt::match.data(m_b)

  # Combine match weights with external weights via the shared helper
  # (defined in R/matching.R). Using the same code path as `fit_matching()`
  # means any row-name invariant violation fails loudly at the bootstrap
  # boundary rather than silently producing NA-tainted or misaligned
  # weights on a subset of replicates.
  matched_weights <- combine_match_and_external_weights(
    matched_b,
    external_weights = weights,
    fit_rows = fit_rows_b
  )

  msm_formula <- stats::reformulate(fit$treatment, response = fit$outcome)
  stats::glm(
    msm_formula,
    data = matched_b,
    weights = matched_weights,
    family = fit$model$family
  )
}


#' Bootstrap variance for longitudinal ICE g-computation
#'
#' @description
#' Resamples **individuals** (all their person-period rows together) and
#' re-runs the full ICE procedure on each bootstrap sample. This is the
#' standard nonparametric bootstrap for longitudinal data (Hernán & Robins,
#' Technical Point 13.1): each resample preserves complete individual
#' trajectories.
#'
#' Unlike point-treatment bootstrap (which resamples rows), ICE bootstrap
#' must resample by individual ID to maintain within-person correlation
#' structure and treatment-confounder feedback.
#'
#' @param fit A `causatr_fit` of type `"longitudinal"`.
#' @param interventions Named list of `causatr_intervention` objects.
#' @param n_boot Positive integer. Number of bootstrap replicates.
#' @param target_within_first Logical vector (length = individuals at
#'   first time) flagging the target population.
#' @param est Character. Estimand string (`"ATE"` for longitudinal).
#' @param subset Quoted expression or `NULL`.
#'
#' @return A k × k variance–covariance matrix (k = number of
#'   interventions).
#'
#' @noRd
ice_variance_bootstrap <- function(
  fit,
  interventions,
  n_boot,
  target_within_first,
  est,
  subset,
  parallel = "no",
  ncpus = 1L
) {
  data <- fit$data
  int_names <- names(interventions)
  id_col <- fit$id
  time_col <- fit$time
  treatment <- fit$treatment
  first_time <- fit$details$time_points[1]

  all_ids <- unique(data[[id_col]])
  n_ids <- length(all_ids)

  # boot_fn: resamples individual IDs (not person-period rows).
  # For each bootstrap sample, reconstcts the person-period data by
  # extracting all rows for the sampled IDs, re-runs fit_ice() +
  # ice_iterate(), and returns marginal means under each intervention.
  orig_weights <- fit$details$weights

  boot_fn <- function(ids, indices) {
    sampled_ids <- ids[indices]

    # Cluster bootstrap trickiness: when an individual is sampled
    # multiple times, we can't just include k copies of their rows
    # under the same `id` — downstream ICE code would treat them as
    # one individual with duplicate rows at each time point. Instead
    # we clone the rows AND assign fresh integer IDs so the ICE
    # recursion sees them as distinct people. Weights are cloned in
    # the same order so they align.
    id_counts <- table(sampled_ids)
    d_b_list <- vector("list", length(sampled_ids))
    w_b_list <- if (!is.null(orig_weights)) {
      vector("list", length(sampled_ids))
    }
    new_id <- 0L
    for (orig_id in names(id_counts)) {
      n_copies <- as.integer(id_counts[[orig_id]])
      orig_rows <- which(data[[id_col]] == orig_id)
      sub <- data[orig_rows]
      sub_w <- if (!is.null(orig_weights)) orig_weights[orig_rows]
      for (cc in seq_len(n_copies)) {
        new_id <- new_id + 1L
        # Deep copy + ID reassignment. Without `copy()` the mutation
        # would leak into the shared `sub` across iterations of the
        # inner loop and corrupt earlier entries of `d_b_list`.
        sub_copy <- data.table::copy(sub)
        sub_copy[, (id_col) := new_id]
        d_b_list[[new_id]] <- sub_copy
        if (!is.null(orig_weights)) w_b_list[[new_id]] <- sub_w
      }
    }
    d_b <- data.table::rbindlist(d_b_list)
    w_b <- if (!is.null(orig_weights)) unlist(w_b_list)

    # Refit the ICE object on the bootstrap sample.
    fit_b <- tryCatch(
      suppressWarnings(
        fit_ice(
          data = d_b,
          outcome = fit$outcome,
          treatment = treatment,
          confounders = fit$confounders,
          confounders_tv = fit$confounders_tv,
          family = fit$family,
          estimand = fit$estimand,
          history = fit$history,
          censoring = fit$censoring,
          weights = w_b,
          model_fn = fit$details$model_fn,
          id = id_col,
          time = time_col,
          call = fit$call
        )
      ),
      error = function(e) NULL
    )
    if (is.null(fit_b)) {
      return(rep(NA_real_, length(int_names)))
    }

    # Determine target population in the bootstrap sample. `d_b` is
    # a data.table (list-like), so `eval()` resolves column references
    # directly; `enclos = parent.frame()` preserves any session-scoped
    # variables referenced by the user's subset expression.
    rows_b_first <- d_b[[time_col]] == first_time
    if (!is.null(subset)) {
      target_b <- rows_b_first &
        as.logical(eval(subset, envir = d_b, enclos = parent.frame()))
    } else {
      target_b <- rows_b_first
    }
    target_b_within <- target_b[rows_b_first]

    # Run ICE for each intervention and compute marginal means.
    w_b_target <- if (!is.null(w_b)) w_b[rows_b_first][target_b_within]
    vapply(
      interventions,
      function(iv) {
        res_b <- tryCatch(
          ice_iterate(fit_b, iv),
          error = function(e) NULL
        )
        if (is.null(res_b)) {
          return(NA_real_)
        }
        maybe_weighted_mean(
          res_b$pseudo_final[target_b_within],
          w_b_target
        )
      },
      numeric(1)
    )
  }

  # Run bootstrap: pass individual IDs as data, boot::boot resamples them.
  boot_res <- boot::boot(
    data = all_ids,
    statistic = boot_fn,
    R = n_boot,
    parallel = parallel,
    ncpus = ncpus
  )

  process_boot_results(boot_res, int_names, n_boot)
}
