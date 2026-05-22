#' Bootstrap variance for SNM blip parameters
#'
#' Resamples individuals (rows for point, complete trajectories for
#' longitudinal) and re-estimates the full g-estimation pipeline on
#' each replicate: treatment model refit + blip parameter solve. The
#' bootstrap statistic is the psi vector itself, not marginal means.
#'
#' For longitudinal fits, uses ID-clustered resampling (resample entire
#' individual trajectories with clone-and-reassign for duplicates),
#' matching the cluster bootstrap used by
#' `ipw_longitudinal_variance_bootstrap()`.
#'
#' @param fit A `causatr_fit` with `estimator = "snm"`.
#' @param treatment_values Numeric vector of length 2 or `NULL`. When
#'   non-`NULL`, the bootstrap statistic is the averaged blip effect
#'   rather than the raw psi vector.
#' @param n_boot Positive integer. Number of bootstrap replicates.
#' @param parallel Character. `"no"`, `"multicore"`, `"snow"`, or
#'   `"future"`.
#' @param ncpus Integer. Number of CPU cores for parallel backends.
#' @return A list with `vcov`, `boot_t`, `boot_info` (same contract
#'   as `process_boot_results()`).
#' @noRd
snm_variance_bootstrap <- function(
  fit,
  treatment_values = NULL,
  n_boot = 500L,
  parallel = "no",
  ncpus = 1L
) {
  if (fit$type == "longitudinal") {
    snm_longitudinal_variance_bootstrap(
      fit,
      n_boot = n_boot,
      parallel = parallel,
      ncpus = ncpus
    )
  } else {
    snm_point_variance_bootstrap(
      fit,
      treatment_values = treatment_values,
      n_boot = n_boot,
      parallel = parallel,
      ncpus = ncpus
    )
  }
}


#' Point-treatment SNM bootstrap
#'
#' Resamples rows, refits the treatment model and solves the
#' g-estimating equation on each bootstrap sample. The statistic is
#' either the raw psi vector (`treatment_values = NULL`) or the
#' scalar averaged blip effect.
#'
#' @param fit A `causatr_fit` with `estimator = "snm"`,
#'   `type = "point"`.
#' @param treatment_values Numeric vector of length 2 or `NULL`.
#' @param n_boot Positive integer.
#' @param parallel Character parallelisation backend.
#' @param ncpus Integer number of cores.
#' @return A list with `vcov`, `boot_t`, `boot_info`.
#' @noRd
snm_point_variance_bootstrap <- function(
  fit,
  treatment_values = NULL,
  n_boot = 500L,
  parallel = "no",
  ncpus = 1L
) {
  data <- fit$data
  blip_spec <- fit$details$blip_spec
  outcome <- fit$outcome
  treatment <- fit$treatment
  confounders <- resolve_confounders_treatment(fit)
  prop_fn <- fit$details$propensity_model_fn
  prop_family <- fit$details$propensity_family
  tf_formula <- fit$details$treatment_free
  dots <- fit$details$dots
  orig_weights <- fit$details$weights
  target <- fit$target

  # Determine psi dimension from a trial treatment design — for
  # categorical treatments, p_psi = (K-1) * p_mod, not blip_spec$n_params
  fit_rows_orig <- get_fit_rows(data, outcome, target = target)
  fit_data_orig <- data[fit_rows_orig]
  trt_model_orig <- fit$details$treatment_model
  td_orig <- snm_treatment_design(
    fit_data_orig,
    treatment,
    blip_spec,
    trt_model_orig
  )

  if (!is.null(treatment_values)) {
    stat_names <- "avg_blip_effect"
  } else {
    stat_names <- td_orig$param_names
  }

  boot_fn <- function(d, indices) {
    tryCatch(
      {
        d_b <- d[indices]
        w_b <- if (!is.null(orig_weights)) orig_weights[indices]

        fit_rows_b <- get_fit_rows(d_b, outcome, target = target)
        fit_data_b <- d_b[fit_rows_b]

        tm_args <- list(
          data = fit_data_b,
          treatment = treatment,
          confounders = confounders,
          model_fn = prop_fn,
          propensity_family = prop_family
        )
        if (!is.null(w_b)) {
          tm_args$weights <- w_b[fit_rows_b]
        }
        trt_model_b <- do.call(
          fit_treatment_model,
          c(tm_args, dots)
        )

        Y_b <- fit_data_b[[outcome]]
        td_b <- snm_treatment_design(
          fit_data_b,
          treatment,
          blip_spec,
          trt_model_b
        )
        AM_b <- td_b$AM
        RM_b <- td_b$RM
        p_psi_b <- td_b$p_psi

        if (!is.null(tf_formula)) {
          Z_b <- stats::model.matrix(tf_formula, data = fit_data_b)
          lhs <- rbind(
            cbind(crossprod(Z_b, Z_b), crossprod(Z_b, AM_b)),
            cbind(crossprod(RM_b, Z_b), crossprod(RM_b, AM_b))
          )
          rhs <- c(
            as.numeric(crossprod(Z_b, Y_b)),
            as.numeric(crossprod(RM_b, Y_b))
          )
          theta_b <- as.numeric(solve(lhs, rhs))
          p_beta <- ncol(Z_b)
          psi_b <- theta_b[p_beta + seq_len(p_psi_b)]
        } else {
          lhs <- crossprod(RM_b, AM_b)
          rhs <- crossprod(RM_b, Y_b)
          psi_b <- as.numeric(solve(lhs, rhs))
        }

        if (!is.null(treatment_values)) {
          delta_a <- treatment_values[2] - treatment_values[1]
          m_bar_b <- colMeans(td_b$M)
          return(delta_a * sum(psi_b * m_bar_b))
        }

        psi_b
      },
      error = function(e) rep(NA_real_, length(stat_names))
    )
  }

  boot_res <- dispatch_boot(
    data = data,
    statistic = boot_fn,
    R = n_boot,
    parallel = parallel,
    ncpus = ncpus
  )

  process_boot_results(boot_res, stat_names, n_boot)
}


#' Longitudinal SNM cluster bootstrap
#'
#' Resamples complete individual trajectories (ID-clustered bootstrap)
#' and re-runs the backward sequential g-estimation on each replicate.
#' Uses the same clone-and-reassign pattern as
#' `ipw_longitudinal_variance_bootstrap()`.
#'
#' @param fit A `causatr_fit` with `estimator = "snm"`,
#'   `type = "longitudinal"`.
#' @param n_boot Positive integer.
#' @param parallel Character parallelisation backend.
#' @param ncpus Integer number of cores.
#' @return A list with `vcov`, `boot_t`, `boot_info`.
#' @noRd
snm_longitudinal_variance_bootstrap <- function(
  fit,
  n_boot = 500L,
  parallel = "no",
  ncpus = 1L
) {
  data <- fit$data
  blip_spec <- fit$details$blip_spec
  outcome <- fit$outcome
  treatment <- fit$treatment
  confounders <- resolve_confounders_treatment(fit)
  confounders_tv <- resolve_confounders_tv_treatment(fit)
  confounders_outcome_raw <- fit$confounders_outcome
  confounders_treatment_raw <- fit$confounders_treatment
  confounders_tv_outcome_raw <- fit$confounders_tv_outcome
  confounders_tv_treatment_raw <- fit$confounders_tv_treatment
  prop_fn <- fit$details$propensity_model_fn
  prop_family <- fit$details$propensity_family
  tf_formula <- fit$details$treatment_free
  model_fn <- fit$details$model_fn
  dots <- fit$details$dots
  orig_weights <- fit$details$weights
  id_col <- fit$id
  time_col <- fit$time
  history <- fit$history

  all_ids <- unique(data[[id_col]])

  n_times <- fit$details$n_times

  # Derive stat_names from a trial longitudinal solve — for categorical
  # treatments, per-stage psi dimension is (K-1) * p_mod, and the
  # param_names include level prefixes.
  trial_result <- compute_snm_blip_longitudinal(fit)
  stat_names <- names(trial_result$psi_hat)

  boot_fn <- function(ids, indices) {
    sampled_ids <- ids[indices]

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
        suppressMessages(
          fit_snm(
            data = d_b,
            outcome = outcome,
            treatment = treatment,
            confounders = confounders,
            confounders_tv = confounders_tv,
            family = fit$family,
            estimand = fit$estimand,
            type = "longitudinal",
            history = history,
            weights = w_b,
            propensity_model_fn = prop_fn,
            propensity_family = prop_family,
            id = id_col,
            time = time_col,
            call = fit$call,
            confounders_outcome = confounders_outcome_raw,
            confounders_tv_outcome_raw = confounders_tv_outcome_raw,
            confounders_tv_treatment_raw = confounders_tv_treatment_raw,
            confounders_treatment_raw = confounders_treatment_raw,
            treatment_free = tf_formula,
            model_fn = model_fn
          )
        ),
        warning = function(w) {
          if (inherits(w, "causatr_singular_bread")) {
            invokeRestart("muffleWarning")
          }
        }
      ),
      error = function(e) NULL
    )
    if (is.null(fit_b)) {
      return(rep(NA_real_, length(stat_names)))
    }

    snm_b <- tryCatch(
      compute_snm_blip_longitudinal(fit_b),
      error = function(e) NULL
    )
    if (is.null(snm_b)) {
      return(rep(NA_real_, length(stat_names)))
    }

    as.numeric(snm_b$psi_hat)
  }

  boot_res <- dispatch_boot(
    data = all_ids,
    statistic = boot_fn,
    R = n_boot,
    parallel = parallel,
    ncpus = ncpus
  )

  process_boot_results(boot_res, stat_names, n_boot)
}
