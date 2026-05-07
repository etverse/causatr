#' Compute positivity diagnostics
#'
#' @description
#' Dispatches to the appropriate positivity helper based on treatment
#' type. Binary treatments get the classic propensity-score quantile
#' table with violation counts (Crump et al. 2009). Continuous and
#' count treatments get a density-range summary: quantiles of
#' \eqn{f(A_i \mid L_i)} and a count of low-density tail
#' observations. Categorical treatments get per-level
#' \eqn{P(A = k \mid L)} quantile summaries. Multivariate
#' treatments get a per-component positivity list.
#'
#' @param fit A `causatr_fit` object.
#' @param ps_bounds Numeric vector of length 2 defining violation
#'   thresholds for binary treatments, or the density tail quantile
#'   for non-binary treatments. Default `c(0.025, 0.975)`.
#' @return A data.table (or list of data.tables for multivariate),
#'   or `NULL` when no positivity diagnostic applies.
#' @noRd
compute_positivity <- function(fit, ps_bounds) {
  treatment <- fit$treatment

  # Multivariate: per-component positivity.
  if (length(treatment) > 1L) {
    return(compute_positivity_multivariate(fit, ps_bounds))
  }

  # IPW with a treatment density model: dispatch by family.
  if (fit$estimator == "ipw" && !is.null(fit$details$treatment_model)) {
    tm <- fit$details$treatment_model
    fam <- tm$family
    if (fam == "bernoulli") {
      return(compute_positivity_binary(fit, ps_bounds))
    }
    if (fam == "gaussian") {
      return(compute_positivity_density(fit, ps_bounds))
    }
    if (fam == "categorical") {
      return(compute_positivity_categorical(fit, ps_bounds))
    }
    if (fam %in% c("poisson", "negbin")) {
      return(compute_positivity_density(fit, ps_bounds))
    }
    return(NULL)
  }

  # Non-IPW estimators (gcomp, matching): binary PS check only.
  compute_positivity_binary(fit, ps_bounds)
}

#' Compute binary propensity-score positivity diagnostics
#'
#' @description
#' The classic propensity-score positivity table: quantiles of
#' \eqn{P(A = 1 \mid L)} plus violation counts. Only applicable
#' to binary 0/1 treatments.
#'
#' @param fit A `causatr_fit` object.
#' @param ps_bounds Numeric vector of length 2.
#' @return A data.table or `NULL`.
#' @noRd
compute_positivity_binary <- function(fit, ps_bounds) {
  treatment <- fit$treatment

  if (length(treatment) > 1L) {
    return(NULL)
  }

  data <- fit$data
  trt_vals <- unique(stats::na.omit(data[[treatment]]))
  if (!all(trt_vals %in% c(0, 1))) {
    return(NULL)
  }

  # Source propensity scores from the fitted model:
  #   IPW:      bernoulli density model in fit$details$propensity_model
  #   Matching: MatchIt distance vector (logistic PS by default)
  #   Gcomp:    fit a quick logistic regression for diagnostics only
  if (fit$estimator == "ipw") {
    tm <- fit$details$treatment_model
    if (is.null(tm) || tm$family != "bernoulli") {
      return(NULL)
    }
    ps <- as.numeric(stats::predict(
      fit$details$propensity_model,
      newdata = fit$data[fit$details$fit_rows],
      type = "response"
    ))
  } else if (!is.null(fit$match_obj)) {
    ps <- fit$match_obj$distance
    if (is.null(ps) || length(ps) == 0L) return(NULL)
  } else {
    ps_formula <- build_ps_formula(fit$confounders, treatment)
    fit_rows <- get_fit_rows(data, fit$outcome, fit$censoring)
    ps_model <- stats::glm(
      ps_formula,
      data = data[fit_rows],
      family = stats::binomial()
    )
    ps <- stats::fitted(ps_model)
  }

  # Crump et al. (2009) tail thresholds. Values near 0 or near 1 mean
  # the treatment is nearly deterministic given covariates, so the
  # corresponding density-ratio weight 1/p or 1/(1-p) diverges.
  # `ps_bounds[1]` flags the lower tail (near-certain control);
  # `ps_bounds[2]` flags the upper tail (near-certain treatment).
  n_low <- sum(ps < ps_bounds[1], na.rm = TRUE)
  n_high <- sum(ps > ps_bounds[2], na.rm = TRUE)
  n_total <- length(ps)

  data.table::data.table(
    statistic = c(
      "min",
      "q25",
      "median",
      "q75",
      "max",
      "n_below_lower",
      "n_above_upper",
      "n_violations",
      "pct_violations"
    ),
    value = c(
      min(ps, na.rm = TRUE),
      stats::quantile(ps, 0.25, na.rm = TRUE),
      stats::quantile(ps, 0.50, na.rm = TRUE),
      stats::quantile(ps, 0.75, na.rm = TRUE),
      max(ps, na.rm = TRUE),
      n_low,
      n_high,
      n_low + n_high,
      round(100 * (n_low + n_high) / n_total, 2)
    )
  )
}

#' Compute density-range positivity diagnostics
#'
#' @description
#' For continuous and count treatments, "positivity" means the fitted
#' conditional density \eqn{f(A_i \mid L_i)} evaluated at the
#' observed treatment value is bounded away from zero. A very small
#' density means the observed treatment is unlikely given the
#' covariates, yielding extreme density-ratio weights. This helper
#' reports quantiles of the fitted density plus a count of low-density
#' observations (those below the 1st percentile of the density
#' distribution).
#'
#' @param fit A `causatr_fit` with an IPW treatment model.
#' @param ps_bounds Numeric vector of length 2 (used to define the
#'   low-density fraction; observations with density below the
#'   `ps_bounds[1]`-th quantile are flagged).
#' @return A data.table with density quantiles and low-density counts.
#' @noRd
compute_positivity_density <- function(fit, ps_bounds) {
  tm <- fit$details$treatment_model
  fit_rows <- fit$details$fit_rows
  fit_data <- fit$data[fit_rows]
  a_obs <- fit_data[[fit$treatment[1]]]

  f_obs <- evaluate_density(tm, a_obs, fit_data)

  # Low-density threshold: observations below the 1st percentile of the
  # density distribution. A low f(A_i|L_i) means the observed treatment
  # value is improbable given covariates, producing a large density-ratio
  # weight g/f. Using a data-adaptive quantile rather than a fixed cutoff
  # keeps the flag interpretable across different treatment scales.
  low_thresh <- stats::quantile(f_obs, 0.01, na.rm = TRUE)
  n_low <- sum(f_obs < low_thresh, na.rm = TRUE)
  n_total <- length(f_obs)

  data.table::data.table(
    statistic = c(
      "min",
      "q25",
      "median",
      "q75",
      "max",
      "n_low_density",
      "pct_low_density"
    ),
    value = c(
      min(f_obs, na.rm = TRUE),
      stats::quantile(f_obs, 0.25, na.rm = TRUE),
      stats::quantile(f_obs, 0.50, na.rm = TRUE),
      stats::quantile(f_obs, 0.75, na.rm = TRUE),
      max(f_obs, na.rm = TRUE),
      n_low,
      round(100 * n_low / n_total, 2)
    )
  )
}

#' Compute per-level positivity for categorical treatments
#'
#' @description
#' For categorical treatments with \eqn{k \ge 2} levels, reports
#' quantiles of \eqn{P(A = k \mid L)} for each level, plus the
#' count of low-probability cells (predicted probability below 0.01).
#' A low predicted probability for a level means some individuals
#' are very unlikely to receive that treatment level given their
#' covariates -- a positivity concern.
#'
#' @param fit A `causatr_fit` with a categorical IPW treatment model.
#' @param ps_bounds Numeric vector of length 2 (unused, kept for
#'   interface consistency).
#' @return A data.table with per-level probability summaries.
#' @noRd
compute_positivity_categorical <- function(fit, ps_bounds) {
  tm <- fit$details$treatment_model
  fit_rows <- fit$details$fit_rows
  fit_data <- fit$data[fit_rows]
  a_obs <- fit_data[[fit$treatment[1]]]
  trt_levels <- tm$levels

  # Get the full predicted probability matrix from the multinomial
  # model. `predict(type = "probs")` returns n x K for K > 2 or a
  # vector for K = 2.
  prob_raw <- stats::predict(
    tm$model,
    newdata = fit_data,
    type = "probs"
  )
  if (is.null(dim(prob_raw))) {
    prob_mat <- cbind(1 - prob_raw, prob_raw)
    colnames(prob_mat) <- trt_levels
  } else {
    prob_mat <- prob_raw
  }

  # Per-level summary: quantiles + low-probability count.
  rows <- lapply(trt_levels, function(lev) {
    p_lev <- prob_mat[, lev]
    n_low <- sum(p_lev < 0.01, na.rm = TRUE)
    data.table::data.table(
      level = lev,
      min = min(p_lev, na.rm = TRUE),
      q25 = stats::quantile(p_lev, 0.25, na.rm = TRUE),
      median = stats::quantile(p_lev, 0.50, na.rm = TRUE),
      q75 = stats::quantile(p_lev, 0.75, na.rm = TRUE),
      max = max(p_lev, na.rm = TRUE),
      n_low_prob = n_low,
      pct_low_prob = round(100 * n_low / length(p_lev), 2)
    )
  })
  data.table::rbindlist(rows)
}

#' Compute per-component positivity for multivariate treatments
#'
#' @description
#' Loops over each component's `causatr_treatment_model` and calls
#' the appropriate single-component positivity helper. Returns a
#' named list of positivity tables, one per treatment component.
#'
#' @param fit A `causatr_fit` with multivariate IPW.
#' @param ps_bounds Numeric vector of length 2.
#' @return A named list of data.tables (one per component), or NULL
#'   for non-IPW fits.
#' @noRd
compute_positivity_multivariate <- function(fit, ps_bounds) {
  if (fit$estimator != "ipw") {
    return(NULL)
  }
  tms <- fit$details$treatment_models
  if (is.null(tms)) {
    return(NULL)
  }

  # Build a fake single-component fit for each component, then
  # dispatch to the family-appropriate helper.
  result <- lapply(names(tms), function(trt_name) {
    tm_k <- tms[[trt_name]]
    fake_fit <- list(
      treatment = trt_name,
      estimator = "ipw",
      data = fit$data,
      confounders = fit$confounders,
      outcome = fit$outcome,
      details = list(
        treatment_model = tm_k,
        propensity_model = tm_k$model,
        fit_rows = fit$details$fit_rows
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
  names(result) <- names(tms)
  result
}
