#' Confidence intervals for a causatr result
#'
#' @description
#' Returns confidence intervals for each intervention mean (E\[Y^a\]) from a
#' `causatr_result` object.
#'
#' @param object A `causatr_result` object.
#' @param parm Ignored. Intervals are returned for all interventions.
#' @param level Numeric. Confidence level. Default `0.95`.
#' @param ... Currently unused.
#'
#' @return A matrix with columns `"lower"` and `"upper"` and one row per
#'   intervention.
#'
#' @examples
#' \dontrun{
#' result <- contrast(fit, interventions = list(a1 = static(1), a0 = static(0)))
#' confint(result)
#' }
#'
#' @seealso [coef.causatr_result()], [contrast()]
#' @export
confint.causatr_result <- function(object, parm, level = 0.95, ...) {
  ci <- cbind(
    lower = object$estimates$ci_lower,
    upper = object$estimates$ci_upper
  )
  rownames(ci) <- object$estimates$intervention
  ci
}
