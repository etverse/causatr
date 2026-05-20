#' Plot a causatr result
#'
#' @description
#' Produces a forest plot of intervention means or pairwise contrasts from a
#' `causatr_result` object using the `forrest` package.
#'
#' The plot adapts to the context:
#' - **Contrast type**: reference line at 0 (difference) or 1 (ratio/OR);
#'   log scale for ratio and OR.
#' - **Outcome family**: axis label uses "risk" (binomial) or "mean" (other).
#' - **Effect modification (`by`)**: sections group subgroup-specific estimates.
#' - **Fit type**: title notes point vs longitudinal analysis.
#'
#' @param x A `causatr_result` object.
#' @param which Character. What to plot: `"contrasts"` (default) for pairwise
#'   comparisons, or `"means"` for intervention-specific marginal means.
#' @param ... Additional arguments passed to [forrest::forrest()] (e.g.
#'   `title`, `stripe`, `dodge`, `cols`, `widths`, `theme`).
#'
#' @return Invisibly returns `x`.
#'
#' @seealso [contrast()], [print.causatr_result()], [forrest::forrest()]
#' @export
plot.causatr_result <- function(x, which = c("contrasts", "means"), ...) {
  check_pkg("forrest")
  which <- rlang::arg_match(which)

  # Three contextual choices drive the plot:
  #   (a) `by` column present -> section the forest plot by subgroup
  #   (b) binary outcome      -> label axes "Risk" instead of "Mean"
  #   (c) contrast type       -> reference line at 0 (diff) or 1 (ratio/OR)
  #                              and log-scale axis for ratios/ORs
  # These are all derivable from the stored result slots rather than
  # requiring the caller to pass format flags.
  has_by <- "by" %in% names(x$estimates)
  is_binary_outcome <- is_binary_family(x$family)
  outcome_noun <- if (is_binary_outcome) "Risk" else "Mean"

  if (which == "contrasts") {
    # Contrast plot: one row per pairwise comparison. Reference
    # line at 0 for differences, 1 for ratios/ORs -- that's the
    # "no effect" value on each scale. Log axis for ratios keeps
    # symmetric multiplicative effects visually balanced.
    dt <- as.data.frame(x$contrasts)
    label_col <- "comparison"
    xlab <- switch(
      x$type,
      difference = paste0(outcome_noun, " difference (95% CI)"),
      ratio = paste0(outcome_noun, " ratio (95% CI)"),
      or = "Odds ratio (95% CI)"
    )
    ref_line <- if (x$type == "difference") 0 else 1
    log_scale <- x$type %in% c("ratio", "or")
    header <- "Contrast"
  } else {
    # Means plot: one row per intervention, on the absolute scale
    # (no reference line, no log axis).
    dt <- as.data.frame(x$estimates)
    label_col <- "intervention"
    xlab <- paste0(outcome_noun, " (95% CI)")
    ref_line <- NULL
    log_scale <- FALSE
    header <- "Intervention"
  }

  if (nrow(dt) == 0L) {
    rlang::warn("No data to plot.")
    return(invisible(x))
  }

  dt$ci_label <- format_ci(dt$estimate, dt$ci_lower, dt$ci_upper)

  estimator_label <- switch(
    x$estimator,
    gcomp = "G-computation",
    ipw = "IPW",
    matching = "Matching",
    x$estimator
  )
  type_label <- if (!is.null(x$fit_type) && x$fit_type == "longitudinal") {
    " (ICE)"
  } else {
    ""
  }
  default_title <- paste0(estimator_label, type_label, ", ", x$estimand)

  # Assemble forrest() arguments. We build a named list so callers
  # can override any default via `...` without us needing to
  # explicitly enumerate every forrest argument here.
  forrest_args <- list(
    data = dt,
    estimate = "estimate",
    lower = "ci_lower",
    upper = "ci_upper",
    label = label_col,
    xlab = xlab,
    log_scale = log_scale,
    header = header,
    title = default_title,
    cols = stats::setNames("ci_label", "Est (95% CI)"),
    stripe = nrow(dt) > 1L
  )

  # Optional params added only when relevant -- adding NULL would
  # trigger `match.arg`-style errors in some forrest versions.
  if (!is.null(ref_line)) {
    forrest_args$ref_line <- ref_line
  }

  if (has_by) {
    forrest_args$section <- "by"
  }

  # Late-bind user overrides. Anything passed in `...` overwrites
  # the defaults we just set -- this is the standard "plot with
  # sensible defaults but let the user customize" pattern.
  dots <- list(...)
  forrest_args[names(dots)] <- dots

  do.call(forrest::forrest, forrest_args)

  invisible(x)
}

#' Format a point estimate and confidence interval
#'
#' @param est Numeric vector of estimates.
#' @param lo Numeric vector of lower CI bounds.
#' @param hi Numeric vector of upper CI bounds.
#' @param digits Integer number of decimal places.
#' @return Character vector.
#' @noRd
format_ci <- function(est, lo, hi, digits = 3L) {
  # `formatC(..., format = "f")` produces fixed-width decimal strings.
  # This matters for the forest plot label column: all rows in the
  # `ci_label` column must have the same character width so the decimal
  # points align visually when displayed in a monospace table cell.
  paste0(
    formatC(est, format = "f", digits = digits),
    " (",
    formatC(lo, format = "f", digits = digits),
    ", ",
    formatC(hi, format = "f", digits = digits),
    ")"
  )
}

#' Plot diagnostics for a causatr fit
#'
#' @description
#' Produces diagnostic visualisations from a [causatr_diag][diagnose]
#' object. The `which` argument selects the plot type:
#'
#' - `"balance"` (default): a Love plot showing standardised mean
#'   differences via `cobalt::love.plot()`. Requires `cobalt`.
#' - `"weights"`: weight-distribution histograms (one panel per
#'   intervention). Binary treatments are faceted by arm (treated /
#'   control); non-binary show the overall distribution.
#'   Requires `tinyplot`.
#' - `"positivity"`: propensity-score histogram for binary treatments;
#'   conditional density histogram for continuous / count treatments;
#'   per-level bar chart for categorical treatments.
#'   Requires `tinyplot`.
#'
#' @param x A `causatr_diag` object returned by [diagnose()].
#' @param which Character. Type of diagnostic plot. One of `"balance"`
#'   (default), `"weights"`, or `"positivity"`.
#' @param log_scale Logical. For `which = "weights"`, apply log10 to
#'   the weight axis. Default `FALSE`.
#' @param stats Character. For `which = "balance"`, which balance
#'   statistic(s) to plot. Default `"m"` (standardised mean
#'   differences). See `cobalt::love.plot()` for options.
#' @param abs Logical. For `which = "balance"`, whether to plot
#'   absolute values. Default `TRUE`.
#' @param thresholds Named numeric vector. For `which = "balance"`,
#'   threshold lines to draw. Default `c(m = 0.1)`.
#' @param ... Additional arguments passed to the plotting backend
#'   (`cobalt::love.plot()` for balance; `tinyplot::plt()` for
#'   weights and positivity).
#'
#' @return Invisibly returns `x`.
#'
#' @seealso [diagnose()], [print.causatr_diag()]
#' @export
plot.causatr_diag <- function(
  x,
  which = c("balance", "weights", "positivity"),
  log_scale = FALSE,
  stats = "m",
  abs = TRUE,
  thresholds = c(m = 0.1),
  ...
) {
  which <- rlang::arg_match(which)
  if (which == "balance") {
    return(plot_diag_balance(x, stats, abs, thresholds, ...))
  }
  if (which == "weights") {
    return(plot_diag_weights(x, log_scale, ...))
  }
  plot_diag_positivity(x, ...)
}

#' Love plot (balance) for diagnostics
#'
#' @param x A `causatr_diag` object.
#' @param stats Character: balance statistics.
#' @param abs Logical: absolute values.
#' @param thresholds Named numeric: threshold lines.
#' @param ... Forwarded to `cobalt::love.plot()`.
#' @return Invisibly returns `x`.
#' @noRd
plot_diag_balance <- function(x, stats, abs, thresholds, ...) {
  check_pkg("cobalt")

  obj <- get_cobalt_object(x)
  if (is.null(obj)) {
    rlang::abort(
      "Love plot requires an IPW or matching fit with a stored weightit/matchit object.",
      .call = FALSE
    )
  }

  p <- if (is.list(obj) && !is.null(obj$formula) && !is.null(obj$data)) {
    cobalt::love.plot(
      obj$formula,
      data = obj$data,
      stats = stats,
      abs = abs,
      thresholds = thresholds,
      var.order = "unadjusted",
      binary = "std",
      ...
    )
  } else {
    cobalt::love.plot(
      obj,
      stats = stats,
      abs = abs,
      thresholds = thresholds,
      var.order = "unadjusted",
      binary = "std",
      ...
    )
  }

  print(p)
  invisible(x)
}

#' Weight-distribution plot for diagnostics
#'
#' @description
#' Histograms of the weight vectors stored in each intervention panel.
#' Binary treatments facet by arm (treated / control); non-binary
#' treatments show the overall distribution. Multi-intervention
#' diagnostics facet by intervention key.
#'
#' @param x A `causatr_diag` object.
#' @param log_scale Logical: log10 x-axis.
#' @param ... Forwarded to `tinyplot::plt()`.
#' @return Invisibly returns `x`.
#' @noRd
plot_diag_weights <- function(x, log_scale, ...) {
  check_pkg("tinyplot")

  panels <- x$per_intervention
  if (is.null(panels) || length(panels) == 0L) {
    rlang::abort(
      "No weight panels available to plot.",
      .call = FALSE
    )
  }

  # Collect weight vectors from each panel into a long data.frame.
  rows <- list()
  for (nm in names(panels)) {
    wts <- panels[[nm]]$weights
    if (is.null(wts)) {
      next
    }

    # Longitudinal: per-period named list. Plot cumulative only.
    if (is.list(wts) && !data.table::is.data.table(wts)) {
      wts <- wts[["cumulative"]]
      if (is.null(wts)) next
    }

    # Reconstruct the raw weight vector from the fit. The weight
    # summary table stores aggregates, not raw vectors. Recompute
    # from the fit object.
    fit <- x$fit
    if (is.null(fit) || fit$estimator != "ipw") {
      next
    }

    is_default_obs <- nm == "obs" && length(panels) == 1L
    iv <- if (is_default_obs) NULL else x$fit$details$interventions[[nm]]

    w_raw <- reconstruct_weights(fit, nm, is_default_obs)
    if (is.null(w_raw)) {
      next
    }

    # Binary: label by arm.
    trt_type <- x$fit_info$treatment_type
    if (trt_type == "binary" && !isTRUE(fit$details$is_multivariate)) {
      fit_rows <- fit$details$fit_rows
      a_obs <- fit$data[fit_rows][[fit$treatment[1]]]
      arm <- ifelse(a_obs == 1, "treated", "control")
      rows[[length(rows) + 1L]] <- data.frame(
        weight = w_raw,
        arm = arm,
        intervention = nm,
        stringsAsFactors = FALSE
      )
    } else {
      rows[[length(rows) + 1L]] <- data.frame(
        weight = w_raw,
        arm = "overall",
        intervention = nm,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(rows) == 0L) {
    rlang::warn("No weight data available to plot.")
    return(invisible(x))
  }

  plot_df <- do.call(rbind, rows)

  log_arg <- if (log_scale) "x" else ""
  multi_intervention <- length(unique(plot_df$intervention)) > 1L
  has_arms <- length(unique(plot_df$arm)) > 1L

  plt_args <- list(
    x = plot_df$weight,
    type = tinyplot::type_histogram(),
    xlab = if (log_scale) "Weight (log scale)" else "Weight",
    ylab = "Frequency",
    main = "Weight distribution",
    log = log_arg
  )

  if (has_arms && multi_intervention) {
    plt_args$by <- plot_df$arm
    plt_args$facet <- plot_df$intervention
  } else if (has_arms) {
    plt_args$by <- plot_df$arm
  } else if (multi_intervention) {
    plt_args$facet <- plot_df$intervention
  }

  dots <- list(...)
  plt_args[names(dots)] <- dots

  do.call(tinyplot::plt, plt_args)
  invisible(x)
}

#' Reconstruct raw weight vector from a fit for plotting
#'
#' @param fit A `causatr_fit`.
#' @param panel_name Character: panel key.
#' @param is_default_obs Logical.
#' @return Numeric weight vector, or NULL.
#' @noRd
reconstruct_weights <- function(fit, panel_name, is_default_obs) {
  if (fit$estimator != "ipw") {
    return(NULL)
  }

  # Multivariate: product weight
  if (isTRUE(fit$details$is_multivariate)) {
    if (is_default_obs) {
      n <- sum(fit$details$fit_rows)
      return(rep(1, n))
    }
    return(NULL)
  }

  tm <- fit$details$treatment_model
  if (is.null(tm)) {
    return(NULL)
  }
  fit_rows <- fit$details$fit_rows
  fit_data <- fit$data[fit_rows]

  if (is_default_obs) {
    if (tm$family == "bernoulli") {
      a_obs <- fit_data[[fit$treatment[1]]]
      p <- as.numeric(stats::predict(
        fit$details$propensity_model,
        newdata = fit_data,
        type = "response"
      ))
      estimand <- fit$estimand %||% "ATE"
      if (estimand == "ATT") {
        return(ifelse(a_obs == 1, 1, p / (1 - p)))
      }
      if (estimand == "ATC") {
        return(ifelse(a_obs == 1, (1 - p) / p, 1))
      }
      return(ifelse(a_obs == 1, 1 / p, 1 / (1 - p)))
    }
    return(rep(1, sum(fit_rows)))
  }

  # Per-intervention density-ratio weights cannot be reconstructed from
  # the stored fit: they require evaluating the counterfactual target
  # density p^*(a) under the specific intervention, which is only
  # available during `contrast()`. A future extension could cache
  # them in `fit$details$intervention_weights`.
  NULL
}

#' Positivity plot for diagnostics
#'
#' @description
#' For binary treatments: propensity-score histogram with optional
#' threshold lines. For continuous/count: conditional density
#' histogram. For categorical: per-level probability histogram.
#'
#' @param x A `causatr_diag` object.
#' @param ... Forwarded to `tinyplot::plt()`.
#' @return Invisibly returns `x`.
#' @noRd
plot_diag_positivity <- function(x, ...) {
  check_pkg("tinyplot")

  fit <- x$fit
  if (is.null(fit)) {
    rlang::abort(
      "Positivity plot requires the original fit object.",
      .call = FALSE
    )
  }

  trt_type <- x$fit_info$treatment_type

  if (trt_type == "binary") {
    return(plot_positivity_binary(x, ...))
  }
  if (trt_type %in% c("continuous", "count")) {
    return(plot_positivity_density(x, ...))
  }
  if (trt_type == "categorical") {
    return(plot_positivity_categorical(x, ...))
  }

  rlang::warn(
    paste0("Positivity plot not supported for treatment type '", trt_type, "'.")
  )
  invisible(x)
}

#' Binary propensity-score histogram
#'
#' @param x A `causatr_diag` object.
#' @param ... Forwarded to `tinyplot::plt()`.
#' @return Invisibly returns `x`.
#' @noRd
plot_positivity_binary <- function(x, ...) {
  fit <- x$fit
  fit_rows <- fit$details$fit_rows
  fit_data <- fit$data[fit_rows]
  a_obs <- fit_data[[fit$treatment[1]]]

  # Source propensity scores.
  if (fit$estimator == "ipw" && !is.null(fit$details$propensity_model)) {
    ps <- as.numeric(stats::predict(
      fit$details$propensity_model,
      newdata = fit_data,
      type = "response"
    ))
  } else if (!is.null(fit$match_obj) && !is.null(fit$match_obj$distance)) {
    ps <- fit$match_obj$distance
  } else {
    ps_formula <- build_ps_formula(
      resolve_confounders_treatment(fit),
      fit$treatment
    )
    ps_model <- stats::glm(
      ps_formula,
      data = fit_data,
      family = stats::binomial()
    )
    ps <- stats::fitted(ps_model)
  }

  arm <- ifelse(a_obs == 1, "treated", "control")
  plt_args <- list(
    x = ps,
    by = arm,
    type = tinyplot::type_density(),
    xlab = "Propensity score P(A=1|L)",
    ylab = "Density",
    main = "Propensity score distribution"
  )

  dots <- list(...)
  plt_args[names(dots)] <- dots
  do.call(tinyplot::plt, plt_args)
  invisible(x)
}

#' Continuous/count conditional density histogram
#'
#' @param x A `causatr_diag` object.
#' @param ... Forwarded to `tinyplot::plt()`.
#' @return Invisibly returns `x`.
#' @noRd
plot_positivity_density <- function(x, ...) {
  fit <- x$fit
  tm <- fit$details$treatment_model
  if (is.null(tm)) {
    rlang::warn("No treatment model available for density plot.")
    return(invisible(x))
  }

  fit_rows <- fit$details$fit_rows
  fit_data <- fit$data[fit_rows]
  a_obs <- fit_data[[fit$treatment[1]]]
  f_obs <- evaluate_density(tm, a_obs, fit_data)

  plt_args <- list(
    x = f_obs,
    type = tinyplot::type_histogram(),
    xlab = "Conditional density f(A|L)",
    ylab = "Frequency",
    main = "Treatment density distribution"
  )

  dots <- list(...)
  plt_args[names(dots)] <- dots
  do.call(tinyplot::plt, plt_args)
  invisible(x)
}

#' Categorical per-level probability histogram
#'
#' @param x A `causatr_diag` object.
#' @param ... Forwarded to `tinyplot::plt()`.
#' @return Invisibly returns `x`.
#' @noRd
plot_positivity_categorical <- function(x, ...) {
  fit <- x$fit
  tm <- fit$details$treatment_model
  if (is.null(tm)) {
    rlang::warn("No treatment model available for categorical plot.")
    return(invisible(x))
  }

  fit_rows <- fit$details$fit_rows
  fit_data <- fit$data[fit_rows]
  trt_levels <- tm$levels

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

  # Long format for tinyplot faceting.
  rows <- lapply(trt_levels, function(lev) {
    data.frame(
      probability = prob_mat[, lev],
      level = lev,
      stringsAsFactors = FALSE
    )
  })
  plot_df <- do.call(rbind, rows)

  plt_args <- list(
    x = plot_df$probability,
    facet = plot_df$level,
    type = tinyplot::type_histogram(),
    xlab = "P(A = k | L)",
    ylab = "Frequency",
    main = "Per-level treatment probability"
  )

  dots <- list(...)
  plt_args[names(dots)] <- dots
  do.call(tinyplot::plt, plt_args)
  invisible(x)
}

#' Extract the cobalt-compatible object from diagnostics
#'
#' @param diag A `causatr_diag` object.
#' @return A `weightit` or `matchit` object for cobalt, or `NULL`.
#' @noRd
get_cobalt_object <- function(diag) {
  fit <- diag$fit
  if (is.null(fit)) {
    return(NULL)
  }

  if (fit$estimator == "ipw" && !is.null(fit$details$treatment_model)) {
    # The self-contained density-ratio engine has no `weightit` object
    # to hand to `cobalt::love.plot()`. Feed it the propensity formula
    # directly -- this reproduces the "unadjusted" love plot (SMDs
    # before weighting), the most universal balance view the engine
    # can surface without committing to one intervention's
    # post-weighting balance.
    ps_formula <- build_ps_formula(
      resolve_confounders_treatment(fit),
      fit$treatment
    )
    fit_rows <- fit$details$fit_rows
    return(list(
      formula = ps_formula,
      data = as.data.frame(fit$data[fit_rows])
    ))
  }
  if (fit$estimator == "matching" && !is.null(fit$match_obj)) {
    return(fit$match_obj)
  }
  NULL
}
