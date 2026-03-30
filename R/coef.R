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
