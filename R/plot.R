#' Plot a causatr result
#'
#' @description
#' Produces a forest plot of intervention means and pairwise contrasts from a
#' `causatr_result` object (via `tinyplot` if available, otherwise base
#' graphics).
#'
#' @param x A `causatr_result` object.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns `x`.
#'
#' @seealso [contrast()], [print.causatr_result()]
#' @export
plot.causatr_result <- function(x, ...) {
  rlang::abort("plot.causatr_result() is not yet implemented.", .call = FALSE)
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

  if (fit$method == "ipw" && !is.null(fit$weights_obj)) {
    return(fit$weights_obj)
  }
  if (fit$method == "matching" && !is.null(fit$match_obj)) {
    return(fit$match_obj)
  }
  NULL
}
