#' Observed-treatment IPW weight distribution summary
#'
#' @description
#' For binary treatments, reconstructs the Horvitz-Thompson weights
#' `1/p` on treated rows and `1/(1-p)` on controls and summarises by
#' arm. For non-binary treatments (continuous, categorical, count),
#' the observed-treatment density-ratio weight is identically 1
#' (natural-course view), so the summary reports unit weights. For
#' multivariate treatments, computes the product of per-component
#' natural-course weights (also unit). Used for the auto-injected
#' `obs` panel under `diagnose(fit)`.
#'
#' @param fit A `causatr_fit` with `estimator = "ipw"`.
#' @return A `data.table` with columns `group`, `n`, `mean`, `sd`,
#'   `min`, `max`, `ess`.
#' @noRd
compute_weight_summary_observed <- function(fit) {
  # Multivariate: natural-course product weight = 1 for all.
  if (isTRUE(fit$details$is_multivariate)) {
    fit_rows <- fit$details$fit_rows
    n <- sum(fit_rows)
    return(summarise_weights_overall(rep(1, n)))
  }

  tm <- fit$details$treatment_model
  fit_rows <- fit$details$fit_rows
  fit_data <- fit$data[fit_rows]

  # Binary: Horvitz-Thompson observed-arm view. The weight formula
  # depends on the estimand; each form is the Bayes numerator f*(A|L)
  # divided by the propensity f(A|L):
  #   ATE: w_i = 1 / p(A_i|L_i)   (f* = 1, a constant marginal)
  #   ATT: numerator is p(A=1|L), so treated get 1 and controls get p/(1-p)
  #   ATC: numerator is p(A=0|L), so controls get 1 and treated get (1-p)/p
  if (tm$family == "bernoulli") {
    a_obs <- fit_data[[fit$treatment[1]]]
    p <- as.numeric(stats::predict(
      fit$details$propensity_model,
      newdata = fit_data,
      type = "response"
    ))
    estimand <- fit$estimand %||% "ATE"
    if (estimand == "ATT") {
      w <- ifelse(a_obs == 1, 1, p / (1 - p))
    } else if (estimand == "ATC") {
      w <- ifelse(a_obs == 1, (1 - p) / p, 1)
    } else {
      w <- ifelse(a_obs == 1, 1 / p, 1 / (1 - p))
    }
    return(summarise_weights_by_arm(w, a_obs))
  }

  # Non-binary (continuous / categorical / count): the "observed"
  # density-ratio weight f(A|L)/f(A|L) = 1. Report unit weights so
  # the panel is still meaningful (ESS = n, no instability).
  n <- sum(fit_rows)
  summarise_weights_overall(rep(1, n))
}

#' Per-intervention IPW weight distribution summary
#'
#' @description
#' Computes the density-ratio weight for a specific intervention and
#' summarises it. Binary treatments are split by arm (treated /
#' control / overall); non-binary treatments report overall only.
#' Multivariate treatments compute the combined product weight
#' across all components.
#'
#' Under static-arm interventions on a binary treatment, the
#' Horvitz-Thompson indicator zeroes the off-arm rows. The "control"
#' row for `static(1)` reports `mean = 0`, `ess = 0` -- this is the
#' correct M-estimation view.
#'
#' @param fit A `causatr_fit` with `estimator = "ipw"`.
#' @param intervention A `causatr_intervention` object (or `NULL` for
#'   natural course). For multivariate, a named list of per-component
#'   interventions.
#' @return A `data.table` with weight summary columns.
#' @noRd
compute_weight_summary_intervention <- function(fit, intervention) {
  # Multivariate: compute combined product weight via the multivariate
  # weight engine.
  if (isTRUE(fit$details$is_multivariate)) {
    tms <- fit$details$treatment_models
    fit_rows <- fit$details$fit_rows
    fit_data <- fit$data[fit_rows]
    w <- compute_density_ratio_weights_mv(
      tms,
      fit_data,
      intervention,
      fit$estimand
    )
    return(summarise_weights_overall(w))
  }

  tm <- fit$details$treatment_model
  fit_rows <- fit$details$fit_rows
  fit_data <- fit$data[fit_rows]

  w <- compute_density_ratio_weights(
    tm,
    fit_data,
    intervention,
    fit$estimand
  )

  # Binary: split by observed treatment arm.
  if (tm$family == "bernoulli") {
    a_obs <- fit_data[[fit$treatment[1]]]
    if (!isTRUE(all(tm$fit_rows))) {
      a_obs <- a_obs[tm$fit_rows]
    }
    return(summarise_weights_by_arm(w, a_obs))
  }

  # Non-binary: overall summary only (no arm split).
  summarise_weights_overall(w)
}

#' Aggregate a weight vector into the (treated / control / overall) summary
#' table used by `compute_weight_summary_*` helpers.
#'
#' @param w Numeric weight vector aligned with `a_obs`.
#' @param a_obs Observed binary treatment vector (length matches `w`).
#' @return A `data.table` with columns `group`, `n`, `mean`, `sd`, `min`,
#'   `max`, `ess`.
#' @noRd
summarise_weights_by_arm <- function(w, a_obs) {
  # Effective sample size (Kish 1965): ESS = (sum w)^2 / sum(w^2).
  # Equals n when all weights are 1, and decreases as weight variance
  # grows. Intuition: a weight w_i > 1 effectively "replicates" that
  # observation, inflating precision estimates if uncorrected. ESS
  # represents the nominal sample size of an unweighted study with
  # equivalent precision. The 0/0 guard handles the all-zero off-arm
  # case under static-arm density-ratio weights (e.g. controls under
  # `static(1)`); ESS = 0 is correct because no row contributes.
  ess <- function(wts) {
    s2 <- sum(wts^2)
    if (s2 == 0) {
      return(0)
    }
    sum(wts)^2 / s2
  }
  masks <- list(a_obs == 1, a_obs == 0, rep(TRUE, length(w)))
  labels <- c("treated", "control", "overall")
  rows <- lapply(seq_along(labels), function(i) {
    w_sub <- w[masks[[i]]]
    data.table::data.table(
      group = labels[i],
      n = length(w_sub),
      mean = if (length(w_sub) == 0L) NA_real_ else mean(w_sub),
      sd = if (length(w_sub) < 2L) NA_real_ else stats::sd(w_sub),
      min = if (length(w_sub) == 0L) NA_real_ else min(w_sub),
      max = if (length(w_sub) == 0L) NA_real_ else max(w_sub),
      ess = if (length(w_sub) == 0L) 0 else ess(w_sub)
    )
  })
  data.table::rbindlist(rows)
}

#' Aggregate a weight vector into a single overall row
#'
#' @description
#' Used by non-binary treatment types where splitting by "treated /
#' control" is not meaningful. Reports a single `overall` row with
#' the same column schema as `summarise_weights_by_arm()` so
#' downstream print / plot code can handle both shapes uniformly.
#'
#' @param w Numeric weight vector.
#' @return A `data.table` with a single `overall` row.
#' @noRd
summarise_weights_overall <- function(w) {
  ess_val <- if (sum(w^2) == 0) 0 else sum(w)^2 / sum(w^2)
  data.table::data.table(
    group = "overall",
    n = length(w),
    mean = mean(w),
    sd = if (length(w) < 2L) NA_real_ else stats::sd(w),
    min = min(w),
    max = max(w),
    ess = ess_val
  )
}

#' Compute matching quality metrics
#'
#' @param fit A `causatr_fit` object (matching estimator only).
#' @return A data.table with total, matched, discarded counts and retention
#'   percentage, or `NULL` for non-matching fits.
#' @noRd
compute_match_quality <- function(fit) {
  if (fit$estimator != "matching" || is.null(fit$match_obj)) {
    return(NULL)
  }

  m <- fit$match_obj
  n_total <- length(m$weights)
  n_matched <- sum(m$weights > 0)
  n_discarded <- n_total - n_matched

  data.table::data.table(
    statistic = c(
      "n_total",
      "n_matched",
      "n_discarded",
      "pct_retained"
    ),
    value = c(
      n_total,
      n_matched,
      n_discarded,
      round(100 * n_matched / n_total, 1)
    )
  )
}
