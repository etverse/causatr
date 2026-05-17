#' Diagnose longitudinal (ICE / longitudinal-IPW) fits
#'
#' @description
#' Dispatches per-period positivity, balance, and weight diagnostics
#' for longitudinal fits. Each panel stores `positivity`, `balance`,
#' and `weights` as named lists keyed by time-point strings
#' (e.g. `"0"`, `"1"`). For longitudinal IPW the cumulative weight
#' summary is appended under the `"cumulative"` key.
#'
#' ICE (longitudinal g-comp) fits have no treatment model, so
#' positivity and weights are NULL; only per-period balance is
#' reported.
#'
#' @param fit A longitudinal `causatr_fit`.
#' @param interventions Named list of interventions (already resolved).
#' @param default_view Logical: TRUE when no explicit interventions
#'   were requested.
#' @param stats Character: balance statistics for cobalt.
#' @param thresholds Named numeric: balance thresholds.
#' @param ps_bounds Numeric(2): positivity bounds.
#' @return A `causatr_diag` object.
#' @noRd
diagnose_longitudinal <- function(
  fit,
  interventions,
  default_view,
  stats,
  thresholds,
  ps_bounds
) {
  is_ipw <- fit$estimator == "ipw"
  time_points <- fit$details$time_points
  tp_labels <- as.character(time_points)

  # Positivity: per-period treatment density diagnostics (IPW only).
  # Each period's positivity is computed independently from its own
  # treatment model so that a single problematic period (e.g. near-
  # deterministic treatment at time 2) is identified precisely rather
  # than masked by the cumulative product weight.
  positivity_shared <- if (is_ipw) {
    compute_positivity_longitudinal(fit, ps_bounds)
  }

  # Balance: per-period covariate balance.
  balance_shared <- compute_balance_longitudinal(
    fit,
    stats,
    thresholds
  )

  censoring_shared <- compute_censoring_diagnostics(fit)
  sampling_shared <- compute_sampling_diagnostics(fit)

  per_intervention <- lapply(seq_along(interventions), function(i) {
    iv <- interventions[[i]]
    is_default_obs <- default_view && i == 1L
    weights_panel <- if (is_ipw) {
      compute_weights_longitudinal(fit, iv, is_default_obs)
    }
    pct_panel <- if (is_default_obs) {
      NULL
    } else {
      compute_pct_intervened_longitudinal(fit, iv)
    }
    list(
      positivity = positivity_shared,
      balance = balance_shared,
      weights = weights_panel,
      censoring = censoring_shared,
      sampling = sampling_shared,
      pct_intervened = pct_panel
    )
  })
  names(per_intervention) <- names(interventions)

  fit_info <- list(
    treatment_type = detect_diag_treatment_type(fit),
    estimand = fit$estimand,
    type = fit$type,
    has_em = !is.null(fit$details$em_info) &&
      length(fit$details$em_info$em_terms) > 0L
  )

  new_causatr_diag(
    per_intervention = per_intervention,
    match_quality = NULL,
    estimator = fit$estimator,
    fit_info = fit_info,
    fit = fit
  )
}


#' Per-period positivity for longitudinal IPW
#'
#' @description
#' Loops over `treatment_models_by_time` and calls the appropriate
#' single-period positivity helper on each. Returns a named list
#' keyed by time-point string, each element a positivity data.table.
#'
#' @param fit A longitudinal IPW `causatr_fit`.
#' @param ps_bounds Numeric(2).
#' @return Named list of data.tables.
#' @noRd
compute_positivity_longitudinal <- function(fit, ps_bounds) {
  tms_by_time <- fit$details$treatment_models_by_time
  fit_data_by_time <- fit$details$fit_data_by_time

  result <- lapply(names(tms_by_time), function(tp) {
    tm_k <- tms_by_time[[tp]]
    data_k <- fit_data_by_time[[tp]]
    # Build a fake single-period fit to reuse the point-IPW positivity
    # helpers without duplicating the dispatch logic. The fake fit only
    # needs the slots that `compute_positivity_*` actually reads.
    fake_fit <- list(
      treatment = fit$treatment,
      estimator = "ipw",
      data = data_k,
      confounders = fit$confounders,
      outcome = fit$outcome,
      details = list(
        treatment_model = tm_k,
        propensity_model = tm_k$model,
        fit_rows = tm_k$fit_rows
      )
    )
    fam <- tm_k$family
    if (fam == "bernoulli") {
      compute_positivity_binary(fake_fit, ps_bounds)
    } else if (fam == "categorical") {
      compute_positivity_categorical(fake_fit, ps_bounds)
    } else {
      compute_positivity_density(fake_fit, ps_bounds)
    }
  })
  names(result) <- names(tms_by_time)
  result
}


#' Per-period covariate balance for longitudinal fits
#'
#' @description
#' For longitudinal IPW: computes unadjusted balance at each period
#' using the per-period data subset and the treatment formula for
#' that period. For ICE: computes balance at each period using the
#' time-varying confounders on the per-period data subset.
#'
#' @param fit A longitudinal `causatr_fit`.
#' @param stats Character: balance statistics for cobalt.
#' @param thresholds Named numeric: balance thresholds.
#' @return Named list of cobalt `bal.tab` objects or data.tables
#'   (one per time point), or NULL.
#' @noRd
compute_balance_longitudinal <- function(fit, stats, thresholds) {
  time_points <- fit$details$time_points
  time_col <- fit$time
  treatment <- fit$treatment
  data <- fit$data

  # Confounders for each period: baseline + time-varying.
  conf_baseline <- fit$confounders
  conf_tv <- fit$confounders_tv

  result <- lapply(as.character(time_points), function(tp) {
    # Subset to this period's rows.
    rows_k <- data[[time_col]] == as.numeric(tp)
    d_k <- data[rows_k]

    # Build the confounder set: baseline terms + time-varying terms
    # that exist as columns in the per-period data.
    baseline_vars <- all.vars(conf_baseline)
    tv_vars <- if (!is.null(conf_tv)) all.vars(conf_tv) else character(0)
    all_vars <- unique(c(baseline_vars, tv_vars))
    available <- intersect(all_vars, names(d_k))
    available <- available[
      vapply(available, function(v) !all(is.na(d_k[[v]])), logical(1))
    ]

    if (length(available) == 0L) {
      return(NULL)
    }

    # Build a formula for cobalt: treatment ~ available confounders.
    ps_formula <- stats::reformulate(available, response = treatment)

    if (rlang::is_installed("cobalt")) {
      tryCatch(
        cobalt::bal.tab(
          ps_formula,
          data = as.data.frame(d_k),
          stats = stats,
          thresholds = thresholds,
          binary = "std"
        ),
        error = function(e) NULL
      )
    } else {
      # Simple fallback.
      trt_vals <- unique(stats::na.omit(d_k[[treatment]]))
      is_binary <- all(trt_vals %in% c(0, 1))
      if (is_binary) {
        compute_balance_simple_binary(
          d_k,
          treatment,
          available
        )
      } else if (is.numeric(d_k[[treatment]])) {
        compute_balance_simple_corr(
          d_k,
          treatment,
          available
        )
      }
    }
  })
  names(result) <- as.character(time_points)
  # Drop NULL entries (periods with no confounders).
  result <- result[!vapply(result, is.null, logical(1))]
  if (length(result) == 0L) NULL else result
}


#' Per-period weight diagnostics for longitudinal IPW
#'
#' @description
#' Computes per-period density-ratio weight summaries and the
#' cumulative (product) weight summary. For the default `obs` panel,
#' each period's "observed" weight is the Horvitz-Thompson weight
#' for binary treatments or unit weight for non-binary. For
#' user-supplied interventions, per-period density-ratio weights are
#' computed and then multiplied into the cumulative weight.
#'
#' @param fit A longitudinal IPW `causatr_fit`.
#' @param intervention A `causatr_intervention` or NULL.
#' @param is_default_obs Logical.
#' @return Named list of data.tables (one per time point +
#'   `"cumulative"`).
#' @noRd
compute_weights_longitudinal <- function(
  fit,
  intervention,
  is_default_obs
) {
  tms_by_time <- fit$details$treatment_models_by_time
  fit_data_by_time <- fit$details$fit_data_by_time
  time_points <- fit$details$time_points
  tp_labels <- as.character(time_points)
  id_col <- fit$id
  stabilize <- !is.null(fit$details$numerator_models_by_time)

  # First-period ids define the canonical ordering. All per-period
  # weight vectors are then aligned to this ordering so the cumulative
  # product (row-wise product of W_per_period) is well-defined.
  first_tp <- tp_labels[1]
  ids_first <- as.character(
    fit_data_by_time[[first_tp]][[id_col]]
  )
  n_id <- length(ids_first)

  # Per-period weight matrix (n_id x K). Each column k holds w_{ik}
  # (the period-k density-ratio weight for individual i) aligned to
  # `ids_first` order. Cumulative weight W_i = prod_k w_{ik} is the
  # row product used in the longitudinal Hajek estimator.
  W_per_period <- matrix(1, nrow = n_id, ncol = length(tp_labels))
  period_summaries <- vector("list", length(tp_labels))
  names(period_summaries) <- tp_labels

  for (k in seq_along(tp_labels)) {
    tp <- tp_labels[k]
    tm_k <- tms_by_time[[tp]]
    data_k <- fit_data_by_time[[tp]]
    ids_k <- as.character(data_k[[id_col]])

    if (is_default_obs) {
      # Default "observed" view: HT weights for binary, unit for rest.
      if (tm_k$family == "bernoulli") {
        a_obs <- data_k[[fit$treatment[1]]]
        p <- as.numeric(stats::predict(
          tm_k$model,
          newdata = data_k,
          type = "response"
        ))
        w_k <- ifelse(a_obs == 1, 1 / p, 1 / (1 - p))
      } else {
        w_k <- rep(1, nrow(data_k))
      }
    } else if (is.null(intervention)) {
      w_k <- rep(1, nrow(data_k))
    } else {
      if (stabilize) {
        tm_num_k <- fit$details$numerator_models_by_time[[tp]]
        w_k <- compute_stabilized_period_weight(
          tm_denom = tm_k,
          tm_num = tm_num_k,
          data = data_k,
          intervention = intervention
        )
      } else {
        w_k <- compute_density_ratio_weights(
          treatment_model = tm_k,
          data = data_k,
          intervention = intervention,
          estimand = fit$estimand
        )
      }
    }

    # Align to canonical id ordering. `tm_k$fit_rows` is logical and
    # relative to `data_k`; applying it here gives the ids of individuals
    # for whom the period-k model was actually fit (non-missing treatment
    # and confounders). Individuals absent from `period_ids` retain the
    # default weight of 1 (neutral for the product).
    period_ids <- ids_k[tm_k$fit_rows]
    w_aligned <- rep(1, n_id)
    pos <- match(period_ids, ids_first)
    w_aligned[pos] <- w_k
    W_per_period[, k] <- w_aligned

    # Per-period summary.
    if (
      tm_k$family == "bernoulli" &&
        (is_default_obs || is.null(intervention))
    ) {
      a_obs <- data_k[[fit$treatment[1]]]
      period_summaries[[tp]] <- summarise_weights_by_arm(
        w_k,
        a_obs[tm_k$fit_rows]
      )
    } else {
      period_summaries[[tp]] <- summarise_weights_overall(w_k)
    }
  }

  # Cumulative (product) weight W_i = prod_{k=0}^{K-1} w_{ik}.
  # Row-wise product of the period weight matrix. Extreme cumulative
  # weights signal compounding instability across periods -- the
  # per-period summaries above help diagnose which period is the source.
  w_cumulative <- apply(W_per_period, 1, prod)
  period_summaries[["cumulative"]] <- summarise_weights_overall(
    w_cumulative
  )

  period_summaries
}
