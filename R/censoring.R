#' Censoring model fitting and IPCW weight computation
#'
#' @description
#' Fits a logistic censoring model \eqn{P(C = 0 \mid A, L)} and computes
#' inverse-probability-of-censoring weights (IPCW) for MAR outcome
#' censoring. The censoring model is always binary (C is 0/1), so the
#' implementation is substantially simpler than `fit_treatment_model()`.
#'
#' **Point treatments:** A single censoring model conditions on treatment
#' and confounders: \eqn{P(C = 0 \mid A, L)}.
#'
#' **Longitudinal treatments:** Per-period censoring models condition on
#' treatment and covariate history up to that period:
#' \eqn{P(C_k = 0 \mid \bar{A}_k, \bar{L}_k, C_{k-1} = 0)}.
#' Cumulative IPCW weight is the product over periods.
#'
#' @section Stabilization:
#' Stabilized IPCW weights use \eqn{P(C = 0)} in the numerator
#' (point) or the product of marginal per-period uncensoring
#' probabilities (longitudinal). Stabilization keeps the weights
#' centered near 1, reducing finite-sample variability.
#'
#' @name censoring
#' @keywords internal
NULL


# ── Point treatment censoring model ─────────────────────────────────

#' Fit a point-treatment censoring model
#'
#' @description
#' Fits \eqn{P(C = 0 \mid A, L)} via a logistic model. The censoring
#' column uses the convention 1 = censored, 0 = uncensored; the model
#' predicts the probability of being **uncensored**.
#'
#' @param data A data.table with columns for treatment, confounders,
#'   and the censoring indicator.
#' @param censoring Character. Column name of the censoring indicator
#'   (1 = censored, 0 = uncensored).
#' @param treatment Character vector. Treatment column name(s).
#' @param confounders One-sided formula for confounders.
#' @param model_fn Model fitting function (default `stats::glm`).
#' @param weights Optional numeric vector of external weights
#'   (e.g. survey weights). Length must equal `nrow(data)`.
#' @param ... Additional arguments passed to `model_fn`.
#'
#' @return A `causatr_censoring_model` S3 object with:
#'   - `model`: the fitted GLM
#'   - `cens_formula`: the formula used
#'   - `fit_rows`: logical vector of rows used for fitting
#'   - `alpha_hat`: coefficient vector
#'   - `X_fit`: design matrix at fit time
#'   - `p_uncensored`: fitted P(C=0|A,L) for fit rows
#'   - `p_marginal`: marginal P(C=0) across fit rows (for stabilization)
#'
#' @noRd
fit_censoring_model <- function(
  data,
  censoring,
  treatment,
  confounders,
  model_fn = stats::glm,
  weights = NULL,
  ...
) {
  check_string(censoring)
  check_col_exists(data, censoring)
  check_formula(confounders)

  cens_vals <- data[[censoring]]

  # Fit rows: all rows with non-missing censoring indicator.
  # Unlike the outcome model, we use ALL non-NA censoring rows
  # (both censored and uncensored) to fit the censoring model.
  fit_rows <- !is.na(cens_vals)
  for (v in treatment) {
    fit_rows <- fit_rows & !is.na(data[[v]])
  }
  for (v in all.vars(confounders)) {
    fit_rows <- fit_rows & !is.na(data[[v]])
  }
  fit_data <- data[fit_rows]

  # Response: uncensored indicator (1 = uncensored, 0 = censored).
  # The column convention is C=1 means censored, so we invert.
  uncens_response <- 1L - as.integer(fit_data[[censoring]])

  # Build formula: .uncens ~ A + confounders
  confounder_terms <- attr(stats::terms(confounders), "term.labels")
  rhs <- c(treatment, confounder_terms)
  cens_formula <- stats::reformulate(rhs, response = ".uncens")

  fit_data$.uncens <- uncens_response

  model_args <- list(
    formula = cens_formula,
    data = fit_data,
    family = stats::binomial()
  )
  weights_fit <- if (is.null(weights)) NULL else weights[fit_rows]
  if (!is.null(weights_fit)) {
    model_args$weights <- weights_fit
  }
  model <- replay_fit(model_fn, model_args, list(...))

  p_uncensored <- stats::fitted(model)
  p_marginal <- mean(uncens_response)

  alpha_hat <- stats::coef(model)
  X_fit <- stats::model.matrix(model)

  structure(
    list(
      model = model,
      cens_formula = cens_formula,
      fit_rows = fit_rows,
      alpha_hat = alpha_hat,
      X_fit = X_fit,
      p_uncensored = p_uncensored,
      p_marginal = p_marginal
    ),
    class = "causatr_censoring_model"
  )
}


#' Compute IPCW weights from a fitted censoring model
#'
#' @description
#' For uncensored rows: stabilized weight = P(C=0) / P(C=0|A,L).
#' For censored rows: weight = 0 (defensive; these rows are excluded
#' by `get_fit_rows()` downstream).
#' For rows not used in fitting (NA censoring): weight = 1.
#'
#' @param censoring_model A `causatr_censoring_model` object.
#' @param n_total Total number of rows in the original data.
#' @param censoring_col Integer vector. The censoring column values
#'   from the full dataset (length `n_total`).
#' @param stabilize Logical. If TRUE (default), use stabilized weights.
#'
#' @return Numeric vector of length `n_total`.
#'
#' @noRd
compute_ipcw_weights <- function(
  censoring_model,
  n_total,
  censoring_col,
  stabilize = TRUE
) {
  fit_rows <- censoring_model$fit_rows
  p_uncens <- censoring_model$p_uncensored
  p_marg <- censoring_model$p_marginal

  w <- rep(1, n_total)

  # For rows that were fit: assign IPCW weights
  fit_idx <- which(fit_rows)
  if (stabilize) {
    w[fit_idx] <- p_marg / p_uncens
  } else {
    w[fit_idx] <- 1 / p_uncens
  }

  # Zero out censored rows (C != 0)
  censored <- !is.na(censoring_col) & censoring_col != 0L
  w[censored] <- 0

  w
}


# ── Longitudinal censoring models ───────────────────────────────────

#' Fit per-period censoring models for longitudinal data
#'
#' @description
#' At each time step k, fits
#' \eqn{P(C_k = 0 \mid \bar{A}_k, \bar{L}_k, C_{k-1} = 0)}
#' on individuals still uncensored at k-1. Returns per-period models
#' and cumulative IPCW weights.
#'
#' The formula-building follows the same lag structure as
#' `build_longitudinal_ps_formula()` in `longitudinal_ipw.R`.
#'
#' @param data A person-period data.table (after `prepare_data()`).
#' @param censoring Character. Censoring indicator column.
#' @param treatment Character. Treatment column name.
#' @param confounders One-sided formula for baseline confounders.
#' @param confounders_tv Optional one-sided formula for time-varying
#'   confounders.
#' @param model_fn Model fitting function.
#' @param id Character. ID column.
#' @param time Character. Time column.
#' @param time_points Sorted unique time values.
#' @param history Integer. Maximum lag depth.
#' @param weights Optional external weight vector (full length).
#'
#' @return A list with:
#'   - `models`: list of `causatr_censoring_model` objects (one per
#'     period with censoring variation; NULL for periods with no censoring).
#'   - `cumulative_weights`: numeric vector of length `nrow(data)`
#'     holding per-row cumulative IPCW weights.
#'   - `per_period_weights`: list of per-period weight vectors.
#'
#' @noRd
fit_censoring_models_longitudinal <- function(
  data,
  censoring,
  treatment,
  confounders,
  confounders_tv = NULL,
  model_fn = stats::glm,
  id,
  time,
  time_points,
  history = 1L,
  weights = NULL
) {
  n_times <- length(time_points)
  max_lag <- if (is.infinite(history)) n_times - 1L else as.integer(history)

  baseline_terms <- attr(stats::terms(confounders), "term.labels")
  tv_vars <- if (!is.null(confounders_tv)) {
    attr(stats::terms(confounders_tv), "term.labels")
  } else {
    character(0)
  }

  models <- vector("list", n_times)
  per_period_weights <- vector("list", n_times)
  n_total <- nrow(data)

  # Track cumulative IPCW weights per row (start at 1)
  cum_w <- rep(1, n_total)

  for (k in seq_len(n_times)) {
    rows_k <- data[[time]] == time_points[k]
    data_k <- data[rows_k]
    cens_k <- data_k[[censoring]]

    # At the first time point, all individuals are present.
    # At later time points, restrict to those uncensored at k-1.
    # In person-period format, censored individuals at k-1 should
    # not have rows at k, but we filter explicitly for safety.
    if (k > 1L) {
      uncens_prev <- is_uncensored(data_k, censoring)
      # If all rows are uncensored or all censored, skip model fitting
      if (all(uncens_prev) || !any(uncens_prev)) {
        per_period_weights[[k]] <- rep(1, sum(rows_k))
        next
      }
    }

    # Skip if no censoring variation at this period
    cens_vals_k <- as.integer(cens_k)
    has_variation <- length(unique(stats::na.omit(cens_vals_k))) > 1L
    if (!has_variation) {
      per_period_weights[[k]] <- rep(1, sum(rows_k))
      next
    }

    # Build formula with lags (same structure as longitudinal IPW)
    available_lags <- min(k - 1L, max_lag)
    cens_formula_k <- build_longitudinal_censoring_formula(
      treatment = treatment,
      baseline_terms = baseline_terms,
      tv_vars = tv_vars,
      available_lags = available_lags,
      data_at_time = data_k
    )

    # Fit the per-period censoring model
    cm_k <- fit_censoring_model(
      data = data_k,
      censoring = censoring,
      treatment = treatment,
      confounders = remove_response(cens_formula_k),
      model_fn = model_fn,
      weights = if (!is.null(weights)) weights[rows_k] else NULL
    )

    models[[k]] <- cm_k

    # Per-period IPCW weights for this time step
    w_k <- compute_ipcw_weights(
      cm_k,
      n_total = sum(rows_k),
      censoring_col = cens_vals_k,
      stabilize = TRUE
    )
    per_period_weights[[k]] <- w_k

    # Accumulate into cumulative weights
    cum_w[rows_k] <- cum_w[rows_k] * w_k
  }

  names(models) <- as.character(time_points)
  names(per_period_weights) <- as.character(time_points)

  list(
    models = models,
    cumulative_weights = cum_w,
    per_period_weights = per_period_weights
  )
}


#' Build the per-period censoring formula for longitudinal IPCW
#'
#' @description
#' Mirrors `build_longitudinal_ps_formula()`: conditions on treatment,
#' baseline confounders, TV confounders, and available lags.
#'
#' @param treatment Character. Treatment column name.
#' @param baseline_terms Character vector. Baseline confounder terms.
#' @param tv_vars Character vector. Time-varying covariate names.
#' @param available_lags Integer. Number of available lags at this period.
#' @param data_at_time data.table. Period-specific data for column validation.
#'
#' @return A two-sided formula with `.uncens` as response (placeholder;
#'   callers use `remove_response()` to extract the RHS).
#'
#' @noRd
build_longitudinal_censoring_formula <- function(
  treatment,
  baseline_terms,
  tv_vars,
  available_lags,
  data_at_time
) {
  rhs_dynamic <- c(treatment, tv_vars)

  if (available_lags > 0L) {
    for (lag_k in seq_len(available_lags)) {
      rhs_dynamic <- c(rhs_dynamic, paste0("lag", lag_k, "_", treatment))
      for (v in tv_vars) {
        rhs_dynamic <- c(rhs_dynamic, paste0("lag", lag_k, "_", v))
      }
    }
  }

  # Drop columns that don't exist or are all-NA at this time
  valid <- vapply(
    rhs_dynamic,
    function(col) {
      col %in% names(data_at_time) && !all(is.na(data_at_time[[col]]))
    },
    logical(1)
  )
  rhs_dynamic <- rhs_dynamic[valid]

  rhs_terms <- c(rhs_dynamic, baseline_terms)
  if (length(rhs_terms) == 0L) {
    rhs_terms <- "1"
  }

  stats::reformulate(rhs_terms, response = ".uncens")
}


# ── Bootstrap replay ────────────────────────────────────────────────

#' Refit the censoring model on a bootstrap sample
#'
#' @description
#' Replays the censoring model fit on a resampled dataset, recomputes
#' IPCW weights. Used by the bootstrap variance engine.
#'
#' @param fit A `causatr_fit` object with `details$ipcw == TRUE`.
#' @param d_b data.table. Bootstrap sample.
#'
#' @return Numeric vector of IPCW weights for the bootstrap sample.
#'
#' @noRd
refit_censoring_weights <- function(fit, d_b) {
  cm <- fit$details$censoring_model
  censoring <- fit$censoring
  treatment <- fit$treatment
  model_fn <- fit$details$censoring_model_fn

  if (inherits(cm, "causatr_censoring_model")) {
    # Point treatment: refit single censoring model
    cm_b <- fit_censoring_model(
      data = d_b,
      censoring = censoring,
      treatment = treatment,
      confounders = remove_response(cm$cens_formula),
      model_fn = model_fn
    )
    compute_ipcw_weights(
      cm_b,
      n_total = nrow(d_b),
      censoring_col = as.integer(d_b[[censoring]]),
      stabilize = TRUE
    )
  } else {
    # Longitudinal: refit per-period models
    cens_result <- fit_censoring_models_longitudinal(
      data = d_b,
      censoring = censoring,
      treatment = treatment,
      confounders = fit$confounders,
      confounders_tv = fit$confounders_tv,
      model_fn = model_fn,
      id = fit$id,
      time = fit$time,
      time_points = fit$details$time_points,
      history = fit$details$max_lag
    )
    cens_result$cumulative_weights
  }
}


# ── Sandwich helper ─────────────────────────────────────────────────

#' Create a closure that maps censoring parameters to IPCW weights
#'
#' @description
#' Returns `function(gamma)` that computes the stabilized IPCW weight
#' vector at candidate parameters `gamma`. Used by
#' `numDeriv::jacobian()` in the sandwich variance engine to compute
#' the cross-derivative \eqn{A_{\beta,\gamma}}.
#'
#' @param censoring_model A `causatr_censoring_model` object.
#' @param n_total Integer. Total rows in the original dataset.
#' @param censoring_col Integer. Full-length censoring indicator.
#' @param stabilize Logical. Use stabilized weights.
#'
#' @return A function `function(gamma)` returning a numeric vector
#'   of IPCW weights.
#'
#' @noRd
make_ipcw_weight_fn <- function(
  censoring_model,
  n_total,
  censoring_col,
  stabilize = TRUE
) {
  X <- censoring_model$X_fit
  fit_rows <- censoring_model$fit_rows
  p_marg <- censoring_model$p_marginal
  fit_idx <- which(fit_rows)
  censored <- !is.na(censoring_col) & censoring_col != 0L

  function(gamma) {
    eta <- as.vector(X %*% gamma)
    p_uncens <- stats::plogis(eta)

    w <- rep(1, n_total)
    if (stabilize) {
      w[fit_idx] <- p_marg / p_uncens
    } else {
      w[fit_idx] <- 1 / p_uncens
    }
    w[censored] <- 0
    w
  }
}


# ── S3 print ─────────────────────────────────────────────��──────────

#' @export
print.causatr_censoring_model <- function(x, ...) {
  cat("<causatr_censoring_model>\n")
  cat("  Formula:", deparse(x$cens_formula), "\n")
  n_fit <- sum(x$fit_rows)
  n_uncens <- sum(x$p_uncensored > 0.5)
  cat("  Observations:", n_fit, "\n")
  cat("  Marginal P(uncensored):", round(x$p_marginal, 3), "\n")
  cat("  Coefficients:", length(x$alpha_hat), "\n")
  invisible(x)
}
