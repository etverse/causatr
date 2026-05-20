#' Compute covariate balance (via cobalt or simple SMD fallback)
#'
#' @param fit A `causatr_fit` object.
#' @param stats Character vector of balance statistics for cobalt.
#' @param thresholds Named list of thresholds for cobalt.
#' @param by Character scalar or `NULL`. When non-`NULL`, passed as
#'   `cluster =` to `cobalt::bal.tab()` for stratified balance.
#' @return A cobalt `bal.tab` object or a data.table of SMDs.
#' @noRd
compute_balance <- function(fit, stats, thresholds, by = NULL) {
  if (!rlang::is_installed("cobalt")) {
    return(compute_balance_simple(fit, by = by))
  }

  # cobalt::bal.tab() accepts `cluster =` for stratified balance.
  cluster_arg <- if (!is.null(by)) {
    fit_rows <- if (!is.null(fit$details$fit_rows)) {
      fit$details$fit_rows
    } else {
      get_fit_rows(fit$data, fit$outcome, fit$censoring)
    }
    fit$data[fit_rows][[by]]
  }

  if (
    fit$estimator == "ipw" &&
      (!is.null(fit$details$treatment_model) ||
        !is.null(fit$details$treatment_models))
  ) {
    trt <- if (length(fit$treatment) > 1L) {
      fit$treatment[1]
    } else {
      fit$treatment
    }
    ps_formula <- build_ps_formula(resolve_confounders_treatment(fit), trt)
    fit_rows <- fit$details$fit_rows
    bal_args <- list(
      x = ps_formula,
      data = as.data.frame(fit$data[fit_rows]),
      stats = stats,
      thresholds = thresholds,
      binary = "std"
    )
    if (!is.null(cluster_arg)) {
      bal_args$cluster <- cluster_arg
    }
    do.call(cobalt::bal.tab, bal_args)
  } else if (fit$estimator == "matching" && !is.null(fit$match_obj)) {
    bal_args <- list(
      x = fit$match_obj,
      un = TRUE,
      stats = stats,
      thresholds = thresholds,
      binary = "std"
    )
    if (!is.null(cluster_arg)) {
      bal_args$cluster <- cluster_arg
    }
    do.call(cobalt::bal.tab, bal_args)
  } else {
    ps_formula <- build_ps_formula(
      resolve_confounders_treatment(fit),
      fit$treatment
    )
    fit_rows <- get_fit_rows(fit$data, fit$outcome, fit$censoring)
    bal_args <- list(
      x = ps_formula,
      data = as.data.frame(fit$data[fit_rows]),
      stats = stats,
      thresholds = thresholds,
      binary = "std"
    )
    if (!is.null(cluster_arg)) {
      cluster_rows <- get_fit_rows(
        fit$data,
        fit$outcome,
        fit$censoring
      )
      bal_args$cluster <- fit$data[cluster_rows][[by]]
    }
    do.call(cobalt::bal.tab, bal_args)
  }
}

#' Compute simple balance table without cobalt
#'
#' @description
#' Minimal fallback when cobalt isn't installed. For binary
#' treatments, computes unadjusted SMDs per confounder. For
#' continuous treatments, computes Pearson correlations between
#' treatment and each confounder. For categorical / multivariate,
#' returns NULL (no simple fallback -- cobalt is needed).
#'
#' @param fit A `causatr_fit` object.
#' @return A data.table of balance metrics, or `NULL`.
#' @noRd
compute_balance_simple <- function(fit, by = NULL) {
  data <- fit$data
  treatment <- fit$treatment
  outcome <- fit$outcome

  if (length(treatment) > 1L) {
    return(NULL)
  }

  fit_rows <- get_fit_rows(data, outcome, fit$censoring)
  d <- data[fit_rows]
  confounder_vars <- all.vars(resolve_confounders_treatment(fit))
  confounder_vars <- intersect(confounder_vars, names(d))

  trt_vals <- unique(stats::na.omit(d[[treatment]]))
  is_binary <- all(trt_vals %in% c(0, 1))

  if (is_binary) {
    return(compute_balance_simple_binary(d, treatment, confounder_vars))
  }

  # Continuous / count: Pearson correlation between treatment and
  # each numeric confounder.
  if (is.numeric(d[[treatment]])) {
    return(compute_balance_simple_corr(d, treatment, confounder_vars))
  }

  # Categorical: no simple fallback.
  NULL
}

#' SMD balance table for binary treatments (cobalt fallback)
#'
#' @param d data.table of analysis-sample rows.
#' @param treatment Character treatment column name.
#' @param confounder_vars Character vector of confounder names.
#' @return A data.table with `variable`, `mean_treated`,
#'   `mean_control`, `smd` columns.
#' @noRd
compute_balance_simple_binary <- function(d, treatment, confounder_vars) {
  rows_1 <- d[[treatment]] == 1
  rows_0 <- d[[treatment]] == 0

  results <- lapply(confounder_vars, function(v) {
    x <- d[[v]]
    if (!is.numeric(x)) {
      return(NULL)
    }
    m1 <- mean(x[rows_1], na.rm = TRUE)
    m0 <- mean(x[rows_0], na.rm = TRUE)
    # Rosenbaum & Rubin (1985) pooled SD: average of treated and control
    # variances (not the overall SD), so the denominator is unaffected
    # by imbalance in group sizes. SMD > 0.1 is the Austin (2009)
    # convention for flagging meaningful imbalance.
    s_pooled <- sqrt(
      (stats::var(x[rows_1], na.rm = TRUE) +
        stats::var(x[rows_0], na.rm = TRUE)) /
        2
    )
    smd <- if (s_pooled > 0) (m1 - m0) / s_pooled else NA_real_
    data.table::data.table(
      variable = v,
      mean_treated = m1,
      mean_control = m0,
      smd = smd
    )
  })

  data.table::rbindlist(results[!vapply(results, is.null, logical(1))])
}

#' Correlation balance table for continuous treatments (cobalt fallback)
#'
#' @param d data.table of analysis-sample rows.
#' @param treatment Character treatment column name.
#' @param confounder_vars Character vector of confounder names.
#' @return A data.table with `variable` and `correlation` columns.
#' @noRd
compute_balance_simple_corr <- function(d, treatment, confounder_vars) {
  a <- d[[treatment]]
  results <- lapply(confounder_vars, function(v) {
    x <- d[[v]]
    if (!is.numeric(x)) {
      return(NULL)
    }
    r <- stats::cor(a, x, use = "pairwise.complete.obs")
    data.table::data.table(variable = v, correlation = r)
  })
  data.table::rbindlist(results[!vapply(results, is.null, logical(1))])
}
