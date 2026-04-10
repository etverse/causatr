#' Extract point estimates from a causatr result
#'
#' @description
#' Returns the point estimates for each intervention mean (E\[Y^a\]) as a
#' named numeric vector.
#'
#' @param object A `causatr_result` object.
#' @param ... Currently unused.
#'
#' @return A named numeric vector with one element per intervention.
#'
#' @examples
#' \dontrun{
#' result <- contrast(fit, interventions = list(a1 = static(1), a0 = static(0)))
#' coef(result)
#' }
#'
#' @seealso [confint.causatr_result()], [contrast()]
#' @export
coef.causatr_result <- function(object, ...) {
  stats::setNames(object$estimates$estimate, object$estimates$intervention)
}

#' Variance-covariance matrix for a causatr result
#'
#' @description
#' Returns the variance-covariance matrix of the intervention-specific
#' marginal means (E\[Y^a\]) from a `causatr_result` object.
#'
#' When `by` was used in [contrast()], returns a named list of per-stratum
#' vcov matrices keyed by stratum level.
#'
#' @param object A `causatr_result` object.
#' @param ... Currently unused.
#'
#' @return A named k x k matrix (k = number of interventions), or a named
#'   list of such matrices when `by` was used.
#'
#' @examples
#' \dontrun{
#' result <- contrast(fit, interventions = list(a1 = static(1), a0 = static(0)))
#' vcov(result)
#' }
#'
#' @seealso [coef.causatr_result()], [confint.causatr_result()], [contrast()]
#' @export
vcov.causatr_result <- function(object, ...) {
  object$vcov
}
