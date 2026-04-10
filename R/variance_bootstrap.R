#' Bootstrap variance–covariance matrix for marginal means
#'
#' @description
#' Estimates Var(μ̂) by resampling the entire estimation pipeline `n_boot`
#' times (Hernán & Robins, Technical Point 13.1).  Works for all estimation
#' methods (g-comp, IPW, matching).
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
  method <- fit$method
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

        model_b <- suppressWarnings(refit_model(fit, d_b))

        target_idx_b <- get_target_idx(d_b, treatment, est, subset)

        vapply(
          interventions,
          function(iv) {
            data_a_b <- apply_intervention(d_b, treatment, iv)
            pred_a_b <- predict(model_b, newdata = data_a_b, type = "response")
            valid <- target_idx_b & !is.na(pred_a_b)
            mean(pred_a_b[valid])
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

  t_mat <- boot_res$t
  complete_rows <- stats::complete.cases(t_mat)
  n_ok <- sum(complete_rows)
  n_fail <- n_boot - n_ok

  if (n_ok < 2L) {
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
    vcov_mat <- stats::var(t_mat[complete_rows, , drop = FALSE])
  }

  rownames(vcov_mat) <- int_names
  colnames(vcov_mat) <- int_names
  vcov_mat
}

#' Refit the full estimation pipeline on a bootstrap sample
#'
#' @description
#' Dispatches to the appropriate refitting logic based on the estimation
#' method stored in `fit$method`.
#'
#' @param fit A `causatr_fit` object (original fit, for extracting formulas etc.).
#' @param d_b A data.table bootstrap sample.
#'
#' @return A fitted model object (glm, glm_weightit, etc.) or NULL on failure.
#'
#' @noRd
refit_model <- function(fit, d_b) {
  if (fit$method == "gcomp") {
    refit_gcomp(fit, d_b)
  } else if (fit$method == "ipw") {
    refit_ipw(fit, d_b)
  } else if (fit$method == "matching") {
    refit_matching(fit, d_b)
  } else {
    rlang::abort(
      paste0("Bootstrap not yet supported for method = '", fit$method, "'."),
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
refit_gcomp <- function(fit, d_b) {
  model_formula <- stats::formula(fit$model)
  family <- fit$model$family
  model_fn <- fit$details$model_fn
  censoring <- fit$censoring
  outcome <- fit$outcome

  fit_rows_b <- is_uncensored(d_b, censoring) & !is.na(d_b[[outcome]])

  model_fn(model_formula, data = d_b[fit_rows_b], family = family)
}

#' Refit IPW propensity weights and MSM on a bootstrap sample
#'
#' @param fit A `causatr_fit` object.
#' @param d_b A data.table bootstrap sample.
#' @return A `glm_weightit` model object.
#' @noRd
refit_ipw <- function(fit, d_b) {
  # Re-estimate propensity-score weights.
  confounder_terms <- attr(stats::terms(fit$confounders), "term.labels")
  ps_formula <- stats::reformulate(confounder_terms, response = fit$treatment)

  fit_rows_b <- !is.na(d_b[[fit$outcome]])
  fit_data_b <- d_b[fit_rows_b]

  w_b <- WeightIt::weightit(
    ps_formula,
    data = fit_data_b,
    estimand = fit$estimand
  )

  # Refit the weighted MSM.
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
refit_matching <- function(fit, d_b) {
  confounder_terms <- attr(stats::terms(fit$confounders), "term.labels")
  ps_formula <- stats::reformulate(confounder_terms, response = fit$treatment)

  fit_rows_b <- !is.na(d_b[[fit$outcome]])
  fit_data_b <- as.data.frame(d_b[fit_rows_b])

  match_args <- list(ps_formula, data = fit_data_b, estimand = fit$estimand)
  if (fit$estimand == "ATE") {
    check_pkg("optmatch")
    match_args$method <- "full"
  }
  m_b <- do.call(MatchIt::matchit, match_args)
  matched_b <- MatchIt::match.data(m_b)

  msm_formula <- stats::reformulate(fit$treatment, response = fit$outcome)
  stats::glm(
    msm_formula,
    data = matched_b,
    weights = matched_b$weights,
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
  boot_fn <- function(ids, indices) {
    sampled_ids <- ids[indices]

    # Handle duplicate IDs from sampling with replacement by assigning
    # new unique IDs to each copy.
    id_counts <- table(sampled_ids)
    d_b_list <- vector("list", length(sampled_ids))
    new_id <- 0L
    for (orig_id in names(id_counts)) {
      n_copies <- as.integer(id_counts[[orig_id]])
      sub <- data[data[[id_col]] == orig_id]
      for (cc in seq_len(n_copies)) {
        new_id <- new_id + 1L
        sub_copy <- data.table::copy(sub)
        sub_copy[, (id_col) := new_id]
        d_b_list[[new_id]] <- sub_copy
      }
    }
    d_b <- data.table::rbindlist(d_b_list)

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
          weights = NULL,
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

    # Determine target population in the bootstrap sample.
    rows_b_first <- d_b[[time_col]] == first_time
    if (!is.null(subset)) {
      target_b <- rows_b_first &
        as.logical(eval(subset, envir = as.list(d_b)))
    } else {
      target_b <- rows_b_first
    }
    target_b_within <- target_b[rows_b_first]

    # Run ICE for each intervention and compute marginal means.
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
        mean(res_b$pseudo_final[target_b_within], na.rm = TRUE)
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

  # Compute variance from the bootstrap replicates.
  t_mat <- boot_res$t
  complete_rows <- stats::complete.cases(t_mat)
  n_ok <- sum(complete_rows)
  n_fail <- n_boot - n_ok

  if (n_ok < 2L) {
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
    vcov_mat <- stats::var(t_mat[complete_rows, , drop = FALSE])
  }

  rownames(vcov_mat) <- int_names
  colnames(vcov_mat) <- int_names
  vcov_mat
}
