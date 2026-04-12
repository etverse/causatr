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
  alpha <- (1 - level) / 2

  # Stratified bootstrap: boot_t is a named list keyed by by-level.
  # Iterate in estimates order (not boot_t insertion order) to guarantee
  # CI rows align with estimates rows.
  if (!is.null(object$boot_t) && is.list(object$boot_t)) {
    n_int <- length(unique(object$estimates$intervention))
    by_levels <- unique(object$estimates$by)
    ci_list <- lapply(by_levels, function(lev) {
      bt <- object$boot_t[[as.character(lev)]]
      if (is.null(bt) || nrow(bt) < 2L) {
        na_row <- matrix(NA_real_, nrow = n_int, ncol = 2)
        colnames(na_row) <- c("lower", "upper")
        return(na_row)
      }
      ci_lev <- t(apply(bt, 2, stats::quantile,
        probs = c(alpha, 1 - alpha), na.rm = TRUE
      ))
      colnames(ci_lev) <- c("lower", "upper")
      ci_lev
    })
    ci <- do.call(rbind, ci_list)
  } else if (!is.null(object$boot_t) && is.matrix(object$boot_t)) {
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
  rownames(ci) <- object$estimates$intervention
  ci
}
