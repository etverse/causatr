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
#' Produces a Love plot (covariate balance) for IPW and matching fits, or a
#' positivity plot for g-computation fits. Uses `cobalt::love.plot()` if
#' the `cobalt` package is available, otherwise `tinyplot`.
#'
#' @param x A `causatr_diag` object returned by [diagnose()].
#' @param ... Additional arguments passed to the underlying plot function.
#'
#' @return Invisibly returns `x`.
#'
#' @seealso [diagnose()], [print.causatr_diag()]
#' @export
plot.causatr_diag <- function(x, ...) {
  rlang::abort("plot.causatr_diag() is not yet implemented.", .call = FALSE)
}
