#' Confidence intervals for a causatr result
#'
#' @description
#' Returns confidence intervals for each intervention mean (E\[Y^a\]) from a
#' `causatr_result` object. When the result was computed with
#' `ci_method = "bootstrap"`, percentile-based CIs are returned from the
#' stored bootstrap replicates. Otherwise, Wald-type CIs are computed from
#' the standard errors.
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
  int_names <- object$estimates$intervention
  alpha <- (1 - level) / 2

  if (!is.null(object$boot_t) && !is.list(object$boot_t)) {
    ci <- t(apply(object$boot_t, 2, stats::quantile,
      probs = c(alpha, 1 - alpha), na.rm = TRUE
    ))
    colnames(ci) <- c("lower", "upper")
  } else {
    est <- object$estimates$estimate
    se <- object$estimates$se
    z <- stats::qnorm(1 - alpha)
    ci <- cbind(
      lower = est - z * se,
      upper = est + z * se
    )
  }
  rownames(ci) <- int_names
  ci
}
