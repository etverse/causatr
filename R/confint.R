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
  # Half-alpha for two-sided CI: 0.95 level -> alpha = 0.025 per tail.
  alpha <- (1 - level) / 2

  # Three cases based on the shape of `object$boot_t`:
  #   (a) list                 -> stratified (by-group) bootstrap,
  #                               one matrix per subgroup
  #   (b) matrix               -> unstratified bootstrap, rows = reps
  #   (c) NULL                 -> sandwich path, reconstruct Wald CI
  #                               from stored estimate/SE
  # The percentile-bootstrap branches recompute CIs here rather than
  # reading `object$estimates$ci_lower/ci_upper` because the stored
  # CIs were computed at the original `conf_level` passed to
  # `contrast()`, not the current `level` argument.
  if (!is.null(object$boot_t) && is.list(object$boot_t)) {
    # Stratified bootstrap: boot_t is a named list keyed by by-level.
    # We iterate in the order of the levels as they appear in
    # `object$estimates$by` to guarantee the returned CI rows line
    # up with the estimates table row-for-row.
    n_int <- length(unique(object$estimates$intervention))
    by_levels <- unique(object$estimates$by)
    ci_list <- lapply(by_levels, function(lev) {
      bt <- object$boot_t[[as.character(lev)]]
      if (is.null(bt) || nrow(bt) < 2L) {
        # Degenerate stratum — return NA so downstream code
        # (print/tidy) still has a well-formed matrix to work with.
        na_row <- matrix(NA_real_, nrow = n_int, ncol = 2)
        colnames(na_row) <- c("lower", "upper")
        return(na_row)
      }
      # Percentile CI: alpha and 1-alpha quantiles across
      # bootstrap replicates, per intervention column.
      ci_lev <- t(apply(
        bt,
        2,
        stats::quantile,
        probs = c(alpha, 1 - alpha),
        na.rm = TRUE
      ))
      colnames(ci_lev) <- c("lower", "upper")
      ci_lev
    })
    ci <- do.call(rbind, ci_list)
  } else if (!is.null(object$boot_t) && is.matrix(object$boot_t)) {
    # Unstratified bootstrap — same percentile formula applied once
    # to the full (R x k) replicate matrix.
    ci <- t(apply(
      object$boot_t,
      2,
      stats::quantile,
      probs = c(alpha, 1 - alpha),
      na.rm = TRUE
    ))
    colnames(ci) <- c("lower", "upper")
  } else {
    # Sandwich path: no replicates stored, so we Wald-reconstruct
    # from the stored point estimate and SE. This respects the
    # user's `level` argument, unlike the stored ci_lower/ci_upper
    # which are frozen at contrast()'s conf_level.
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
