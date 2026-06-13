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
#' @param boot_ci Character or `NULL`. Override the bootstrap CI flavour:
#'   `"percentile"` (empirical replicate quantiles) or `"normal"` (Wald from
#'   the bootstrap SE). `NULL` (default) uses the convention recorded on the
#'   result by [contrast()]. Ignored for sandwich results.
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
confint.causatr_result <- function(
  object,
  parm,
  level = 0.95,
  boot_ci = NULL,
  ...
) {
  # Half-alpha for two-sided CI: 0.95 level -> alpha = 0.025 per tail.
  alpha <- (1 - level) / 2

  # Resolve the bootstrap CI flavour: an explicit override, else the
  # convention recorded on the result, else percentile. Percentile reads
  # empirical replicate quantiles; normal falls through to the Wald branch
  # (estimate +/- z * bootstrap SE), so it shares the sandwich reconstruction.
  bc <- if (!is.null(boot_ci)) {
    rlang::arg_match(boot_ci, c("percentile", "normal"))
  } else {
    object$boot_ci %||% "percentile"
  }
  use_perc <- !is.null(object$boot_t) && bc == "percentile"

  # Multiple-imputation results (`ci_method = "rubin"` / `"boot_mi"`) carry
  # per-row degrees of freedom in the "mi_details" attribute. Reconstruct
  # t-based CIs from the pooled SE and df so the user's `level` is honored,
  # matching how the Rubin / Boot MI engines built the stored bounds.
  mi <- attr(object, "mi_details")
  if (!is.null(mi) && object$ci_method %in% c("rubin", "boot_mi")) {
    est <- object$estimates$estimate
    se <- object$estimates$se
    df <- mi$estimates$df
    crit <- stats::qt(1 - alpha, df = df)
    ci <- cbind(lower = est - crit * se, upper = est + crit * se)
    rownames(ci) <- object$estimates[[
      if ("parameter" %in% names(object$estimates)) {
        "parameter"
      } else {
        "intervention"
      }
    ]]
    return(ci)
  }

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
  if (use_perc && is.list(object$boot_t)) {
    # Stratified bootstrap: boot_t is a named list keyed by by-level.
    # We iterate in the order of the levels as they appear in
    # `object$estimates$by` to guarantee the returned CI rows line
    # up with the estimates table row-for-row.
    by_levels <- unique(object$estimates$by)
    # Rows the CI must emit per stratum. For a scalar outcome this is the
    # number of interventions; for a multinomial outcome each stratum's
    # `boot_t` is the flat (intervention x class) matrix, so it is
    # interventions x classes. Deriving it from the estimates table (which is
    # balanced across strata) covers both without special-casing.
    n_rows_per_stratum <- nrow(object$estimates) / length(by_levels)
    ci_list <- lapply(by_levels, function(lev) {
      bt <- object$boot_t[[as.character(lev)]]
      if (is.null(bt) || nrow(bt) < 2L) {
        # Degenerate stratum -- return NA so downstream code
        # (print/tidy) still has a well-formed matrix to work with.
        na_row <- matrix(NA_real_, nrow = n_rows_per_stratum, ncol = 2)
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
  } else if (use_perc && is.matrix(object$boot_t)) {
    # Unstratified bootstrap -- same percentile formula applied once
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
    # Sandwich path, or bootstrap with `boot_ci = "normal"`: Wald-reconstruct
    # from the stored point estimate and SE (the bootstrap SE under the normal
    # flavour). This respects the user's `level` argument, unlike the stored
    # ci_lower/ci_upper which are frozen at contrast()'s conf_level.
    est <- object$estimates$estimate
    se <- object$estimates$se
    z <- stats::qnorm(1 - alpha)
    ci <- cbind(
      lower = est - z * se,
      upper = est + z * se
    )
  }
  # Label CI rows by whichever key the estimates table carries: SNM uses a
  # blip `parameter` column, a multinomial outcome uses (intervention, class)
  # pairs (the bootstrap `boot_t` columns are in the same class-major order as
  # the estimates rows), otherwise the intervention name.
  rownames(ci) <- if ("parameter" %in% names(object$estimates)) {
    object$estimates$parameter
  } else if ("class" %in% names(object$estimates)) {
    paste(object$estimates$intervention, object$estimates$class, sep = ":")
  } else {
    object$estimates$intervention
  }
  ci
}
