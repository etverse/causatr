#' Bootstrap variance for longitudinal IPW
#'
#' @description
#' Resamples **individuals** (all their person-period rows together)
#' and re-runs the full longitudinal-IPW pipeline on each replicate.
#' Mirrors `ice_variance_bootstrap()` exactly -- same id-cluster logic
#' (clone rows + reassign integer ids on duplicates) and same target /
#' subset / external-weight handling -- but the per-replicate body
#' refits via `fit_longitudinal_ipw()` and reads marginal means via
#' `compute_ipw_contrast_longitudinal()`.
#'
#' @inheritParams ice_variance_bootstrap
#'
#' @return A list with `vcov`, `boot_t`, `boot_info`.
#'
#' @noRd
ipw_longitudinal_variance_bootstrap <- function(
  fit,
  interventions,
  n_boot,
  target_within_first,
  est,
  subset,
  parallel = "no",
  ncpus = 1L,
  subset_env = parent.frame()
) {
  data <- fit$data
  int_names <- names(interventions)
  id_col <- fit$id
  time_col <- fit$time
  treatment <- fit$treatment
  first_time <- fit$details$time_points[1]

  all_ids <- unique(data[[id_col]])
  n_ids <- length(all_ids)
  # When IPCW is active, use pre-IPCW weights and let refit_ipw()
  # recompute censoring weights on each bootstrap resample.
  orig_weights <- if (isTRUE(fit$details$ipcw)) {
    fit$details$weights_pre_ipcw
  } else {
    fit$details$weights
  }

  boot_fn <- function(ids, indices) {
    sampled_ids <- ids[indices]

    # Cluster bootstrap clone-and-reassign (same logic as
    # `ice_variance_bootstrap()`). When an id is sampled multiple
    # times, replicate its rows with fresh integer ids so the
    # downstream pipeline treats the duplicates as distinct
    # individuals.
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
        sub_copy <- data.table::copy(sub)
        sub_copy[, (id_col) := new_id]
        d_b_list[[new_id]] <- sub_copy
        if (!is.null(orig_weights)) w_b_list[[new_id]] <- sub_w
      }
    }
    d_b <- data.table::rbindlist(d_b_list)
    w_b <- if (!is.null(orig_weights)) unlist(w_b_list)

    fit_b <- tryCatch(
      withCallingHandlers(
        suppressWarnings(refit_ipw(fit, d_b, weights = w_b)),
        warning = function(w) {
          if (inherits(w, "causatr_singular_bread")) {
            invokeRestart("muffleWarning")
          }
        }
      ),
      error = function(e) NULL
    )
    if (is.null(fit_b)) {
      return(rep(NA_real_, length(int_names)))
    }

    # Determine the bootstrap-replicate target population at the first
    # period. Same shape as `ice_variance_bootstrap()`: evaluate the
    # user's `subset` against `d_b` if provided, AND with the
    # first-period mask, then collapse to a length-n_ids logical.
    rows_b_first <- d_b[[time_col]] == first_time
    if (!is.null(subset)) {
      target_b <- rows_b_first &
        as.logical(eval(subset, envir = d_b, enclos = subset_env))
    } else {
      target_b <- rows_b_first
    }
    target_b_within <- target_b[rows_b_first]

    ipw_long_b <- tryCatch(
      compute_ipw_contrast_longitudinal(
        fit_b,
        interventions,
        target_b_within
      ),
      error = function(e) NULL
    )
    if (is.null(ipw_long_b)) {
      return(rep(NA_real_, length(int_names)))
    }

    ipw_long_b$mu_hat
  }

  boot_res <- dispatch_boot(
    data = all_ids,
    statistic = boot_fn,
    R = n_boot,
    parallel = parallel,
    ncpus = ncpus
  )

  process_boot_results(boot_res, int_names, n_boot)
}


#' Bootstrap variance for longitudinal ICE g-computation
#'
#' @description
#' Resamples **individuals** (all their person-period rows together) and
#' re-runs the full ICE procedure on each bootstrap sample. This is the
#' standard nonparametric bootstrap for longitudinal data (Hernan & Robins,
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
#' @return A k x k variance-covariance matrix (k = number of
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
  ncpus = 1L,
  subset_env = parent.frame()
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
  # For each bootstrap sample, reconstructs the person-period data by
  # extracting all rows for the sampled IDs, re-runs fit_ice() +
  # ice_iterate(), and returns marginal means under each intervention.
  # When IPCW is active, use pre-IPCW weights and refit censoring
  # models on each resample to capture estimation uncertainty.
  orig_weights <- if (isTRUE(fit$details$ipcw)) {
    fit$details$weights_pre_ipcw
  } else {
    fit$details$weights
  }

  boot_fn <- function(ids, indices) {
    sampled_ids <- ids[indices]

    # Cluster bootstrap trickiness: when an individual is sampled
    # multiple times, we can't just include k copies of their rows
    # under the same `id` -- downstream ICE code would treat them as
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

    # IPCW: refit per-period censoring models on the bootstrap
    # sample and compose weights before fitting ICE.
    if (isTRUE(fit$details$ipcw)) {
      ipcw_w_b <- refit_censoring_weights(fit, d_b)
      w_b <- if (is.null(w_b)) {
        ipcw_w_b
      } else {
        w_b * ipcw_w_b
      }
    }

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
    # directly; `subset_env` was captured at `contrast()`'s top level
    # and threaded down so session-scoped variables referenced by the
    # user's subset expression (e.g. `quote(age > cutoff)`) still
    # resolve from the user's frame, not the boot worker's.
    rows_b_first <- d_b[[time_col]] == first_time
    if (!is.null(subset)) {
      target_b <- rows_b_first &
        as.logical(eval(subset, envir = d_b, enclos = subset_env))
    } else {
      target_b <- rows_b_first
    }
    target_b_within <- target_b[rows_b_first]

    # Run ICE for each intervention and compute marginal means.
    # Each call to `ice_iterate()` refits all K period outcome models
    # for that intervention: Q_k(\bar{a}) depends on the specific
    # intervention value applied at period k, so the K models cannot be
    # shared across interventions or pre-computed outside the bootstrap.
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

  # Run bootstrap: pass individual IDs as data, the dispatcher
  # resamples them via `boot::boot()` or `future.apply` as configured.
  boot_res <- dispatch_boot(
    data = all_ids,
    statistic = boot_fn,
    R = n_boot,
    parallel = parallel,
    ncpus = ncpus
  )

  process_boot_results(boot_res, int_names, n_boot)
}


#' Bootstrap variance for longitudinal AIPW (ICE-AIPW)
#'
#' @description
#' ID-clustered nonparametric bootstrap for the longitudinal AIPW
#' estimator. Resamples individuals (complete trajectories), refits
#' both propensity and outcome models on each replicate, and runs
#' the augmented backward iteration via `ice_aipw_iterate()`.
#'
#' Same resampling logic as `ice_variance_bootstrap()`: clone rows
#' with fresh integer IDs so multiply-sampled individuals are treated
#' as distinct people in the ICE recursion.
#'
#' @param fit A `causatr_fit` with `estimator = "aipw"`,
#'   `type = "longitudinal"`.
#' @param interventions Named list of interventions.
#' @param n_boot Positive integer. Number of replicates.
#' @param target_within_first Logical vector over first-time-point
#'   rows flagging the target population.
#' @param est Character estimand label.
#' @param subset Quoted subset expression or `NULL`.
#' @param parallel Bootstrap parallelisation backend.
#' @param ncpus Number of cores.
#' @param subset_env Environment for evaluating `subset`.
#'
#' @return A list with `vcov`, `boot_t`, `boot_info` (same contract
#'   as other `*_variance_bootstrap()` functions).
#'
#' @noRd
aipw_longitudinal_variance_bootstrap <- function(
  fit,
  interventions,
  n_boot,
  target_within_first,
  est,
  subset,
  parallel = "no",
  ncpus = 1L,
  subset_env = parent.frame()
) {
  data <- fit$data
  int_names <- names(interventions)
  id_col <- fit$id
  time_col <- fit$time
  treatment <- fit$treatment
  first_time <- fit$details$time_points[1]

  all_ids <- unique(data[[id_col]])

  orig_weights <- if (isTRUE(fit$details$ipcw)) {
    fit$details$weights_pre_ipcw
  } else {
    fit$details$weights
  }

  boot_fn <- function(ids, indices) {
    sampled_ids <- ids[indices]

    # Clone individuals with fresh IDs (same as ICE bootstrap)
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
      sub_w <- if (!is.null(orig_weights)) {
        orig_weights[orig_rows]
      }
      for (cc in seq_len(n_copies)) {
        new_id <- new_id + 1L
        sub_copy <- data.table::copy(sub)
        sub_copy[, (id_col) := new_id]
        d_b_list[[new_id]] <- sub_copy
        if (!is.null(orig_weights)) {
          w_b_list[[new_id]] <- sub_w
        }
      }
    }
    d_b <- data.table::rbindlist(d_b_list)
    w_b <- if (!is.null(orig_weights)) unlist(w_b_list)

    # IPCW: refit censoring models on bootstrap sample
    if (isTRUE(fit$details$ipcw)) {
      ipcw_w_b <- refit_censoring_weights(fit, d_b)
      w_b <- if (is.null(w_b)) {
        ipcw_w_b
      } else {
        w_b * ipcw_w_b
      }
    }

    # Refit the full longitudinal AIPW object on the bootstrap sample.
    # Both the K propensity models and the K outcome models must be
    # refit here: the augmentation term in ICE-AIPW couples the two
    # (each Q_k update uses g-weights from the propensity models), so
    # fixing either set at the original estimates would understate
    # variance by ignoring nuisance-estimation uncertainty.
    fit_b <- tryCatch(
      suppressWarnings(
        fit_aipw_longitudinal(
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
          propensity_model_fn = fit$details$propensity_model_fn,
          propensity_family = fit$details$propensity_family,
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

    # Determine target in bootstrap sample
    rows_b_first <- d_b[[time_col]] == first_time
    if (!is.null(subset)) {
      target_b <- rows_b_first &
        as.logical(
          eval(
            subset,
            envir = d_b,
            enclos = subset_env
          )
        )
    } else {
      target_b <- rows_b_first
    }
    target_b_within <- target_b[rows_b_first]

    w_b_target <- if (!is.null(w_b)) {
      w_b[rows_b_first][target_b_within]
    }

    # Run ICE-AIPW for each intervention
    vapply(
      interventions,
      function(iv) {
        res_b <- tryCatch(
          ice_aipw_iterate(fit_b, iv),
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

  boot_res <- dispatch_boot(
    data = all_ids,
    statistic = boot_fn,
    R = n_boot,
    parallel = parallel,
    ncpus = ncpus
  )

  process_boot_results(boot_res, int_names, n_boot)
}
