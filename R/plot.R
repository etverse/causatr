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
    # line at 0 for differences, 1 for ratios/ORs — that's the
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

  # Optional params added only when relevant — adding NULL would
  # trigger `match.arg`-style errors in some forrest versions.
  if (!is.null(ref_line)) {
    forrest_args$ref_line <- ref_line
  }

  if (has_by) {
    forrest_args$section <- "by"
  }

  # Late-bind user overrides. Anything passed in `...` overwrites
  # the defaults we just set — this is the standard "plot with
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
#' Produces a Love plot showing covariate balance before and after adjustment.
#' Uses `cobalt::love.plot()` if the `cobalt` package is available.
#'
#' For IPW fits, the plot shows balance before and after weighting. For
#' matching fits, it shows balance before and after matching. For g-computation
#' fits, it shows unadjusted balance only.
#'
#' @param x A `causatr_diag` object returned by [diagnose()].
#' @param stats Character. Which balance statistic(s) to plot. Default `"m"`
#'   (standardised mean differences). See `cobalt::love.plot()` for options.
#' @param abs Logical. Whether to plot absolute values. Default `TRUE`.
#' @param thresholds Named numeric vector. Threshold lines to draw on the
#'   plot. Default `c(m = 0.1)`.
#' @param ... Additional arguments passed to `cobalt::love.plot()`.
#'
#' @return A `ggplot` object (invisibly).
#'
#' @seealso [diagnose()], [print.causatr_diag()]
#' @export
plot.causatr_diag <- function(
  x,
  stats = "m",
  abs = TRUE,
  thresholds = c(m = 0.1),
  ...
) {
  check_pkg("cobalt")

  obj <- get_cobalt_object(x)
  if (is.null(obj)) {
    rlang::abort(
      "Love plot requires an IPW or matching fit with a stored weightit/matchit object.",
      .call = FALSE
    )
  }

  p <- cobalt::love.plot(
    obj,
    stats = stats,
    abs = abs,
    thresholds = thresholds,
    var.order = "unadjusted",
    binary = "std",
    ...
  )

  print(p)
  invisible(p)
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

  if (fit$estimator == "ipw" && !is.null(fit$weights_obj)) {
    return(fit$weights_obj)
  }
  if (fit$estimator == "matching" && !is.null(fit$match_obj)) {
    return(fit$match_obj)
  }
  NULL
}
