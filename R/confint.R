#' Confidence intervals for a causatr result
#'
#' @description
#' Returns confidence intervals for each intervention mean (E\[Y^a\]) from a
#' `causatr_result` object. By default uses the stored CIs; if `level` differs
#' from the level used in [contrast()], recomputes from the vcov matrix.
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
#' confint(result, level = 0.99)
#' }
#'
#' @seealso [coef.causatr_result()], [contrast()]
#' @export
confint.causatr_result <- function(object, parm, level = 0.95, ...) {
  est <- object$estimates$estimate
  se <- object$estimates$se
  z <- stats::qnorm((1 + level) / 2)
  ci <- cbind(
    lower = est - z * se,
    upper = est + z * se
  )
  rownames(ci) <- object$estimates$intervention
  ci
}
