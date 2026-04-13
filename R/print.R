#' Print a causatr fit
#'
#' @description
#' Displays a compact summary of a [causatr_fit][causat] object, showing the
#' estimation method, treatment type, outcome, treatment variable, and sample
#' size.
#'
#' @param x A `causatr_fit` object.
#' @param ... Currently unused.
#' @return Invisibly returns `x`.
#' @seealso [summary.causatr_fit()], [causat()]
#' @export
print.causatr_fit <- function(x, ...) {
  # Pretty-print the method. Defaults-through the raw value if the
  # user has somehow set an unrecognized method so print() never
  # errors on a malformed fit object.
  method_label <- switch(
    x$method,
    gcomp = "G-computation",
    ipw = "IPW (WeightIt)",
    matching = "Matching (MatchIt)",
    x$method
  )
  family_label <- format_family(x$family)
  # Collapse multivariate treatments into a comma-separated label.
  trt_label <- paste(x$treatment, collapse = ", ")

  cat("<causatr_fit>\n")
  cat(" Method:     ", method_label, "\n", sep = "")
  cat(" Type:       ", x$type, "\n", sep = "")
  cat(" Outcome:    ", x$outcome, " (", family_label, ")\n", sep = "")
  cat(" Treatment:  ", trt_label, "\n", sep = "")
  cat(" Estimand:   ", x$estimand, "\n", sep = "")
  cat(" Confounders:", deparse(x$confounders), "\n", sep = " ")
  # TV confounders and id/time are only shown when relevant — keeping
  # the display compact for point-treatment fits.
  if (!is.null(x$confounders_tv)) {
    cat(" TV conf.:  ", deparse(x$confounders_tv), "\n", sep = " ")
  }
  if (!is.null(x$id)) {
    cat(" ID:         ", x$id, "\n", sep = "")
    cat(" Time:       ", x$time, "\n", sep = "")
  }
  cat(" N:          ", nrow(x$data), "\n", sep = "")
  invisible(x)
}

#' Format a family object for display
#'
#' @param family A character string, family object, or function.
#' @return A character string.
#' @noRd
format_family <- function(family) {
  # Three-shape dispatch matching `resolve_family()` but producing
  # a display string rather than a family object. Used only for
  # print() output, so the fallback is "<custom>" rather than an
  # abort — printing a fit should never error out.
  if (is.character(family)) {
    # Already a string like "gaussian" — pass through.
    return(family)
  }
  if (is.list(family) && !is.null(family$family)) {
    # Already a family object (list with $family and $link slots):
    # format as "gaussian(identity)", matching base R's convention.
    return(paste0(family$family, "(", family$link, ")"))
  }
  if (is.function(family)) {
    # A family closure (e.g. `stats::binomial`) — try to evaluate it.
    # Wrapped in tryCatch in case the closure needs arguments; if it
    # fails we fall through to "<custom>".
    fam <- tryCatch(family(), error = function(e) NULL)
    if (!is.null(fam) && !is.null(fam$family)) {
      return(paste0(fam$family, "(", fam$link, ")"))
    }
  }
  "<custom>"
}

#' Print a causatr result
#'
#' @description
#' Displays the estimation method, contrast type, CI method, sample size,
#' intervention-specific marginal means, and pairwise contrasts (with SEs
#' and confidence intervals).
#'
#' @param x A `causatr_result` object.
#' @param ... Currently unused.
#' @return Invisibly returns `x`.
#' @seealso [summary.causatr_result()], [contrast()]
#' @export
print.causatr_result <- function(x, ...) {
  # Human-readable labels for the header block. Same defaults-through
  # pattern as `print.causatr_fit` in case of unrecognized slot values.
  method_label <- switch(
    x$method,
    gcomp = "G-computation",
    ipw = "IPW",
    matching = "Matching",
    x$method
  )
  type_label <- switch(
    x$type,
    difference = "Difference",
    ratio = "Ratio",
    or = "Odds ratio"
  )

  cat("<causatr_result>\n")
  cat(" Method:    ", method_label, "\n", sep = "")
  cat(" Estimand:  ", x$estimand, "\n", sep = "")
  cat(" Contrast:  ", type_label, "\n", sep = "")
  cat(" CI method: ", x$ci_method, "\n", sep = "")
  cat(" N:         ", x$n, "\n", sep = "")

  # When `by = ...` was used, the estimates / contrasts tables have
  # an extra `by` column. Section headers acknowledge this so the
  # reader knows to interpret each row as a per-stratum estimate.
  has_by <- "by" %in% names(x$estimates)

  if (has_by) {
    cat("\nIntervention means (by subgroup):\n")
  } else {
    cat("\nIntervention means:\n")
  }
  print(x$estimates, digits = 4)

  if (has_by) {
    cat("\nContrasts (by subgroup):\n")
  } else {
    cat("\nContrasts:\n")
  }
  print(x$contrasts, digits = 4)
  invisible(x)
}

#' Print causatr diagnostics
#'
#' @description
#' Displays positivity summaries, covariate balance tables, weight
#' distributions (IPW), and match quality metrics (matching) from a
#' [causatr_diag][diagnose] object.
#'
#' @param x A `causatr_diag` object.
#' @param ... Currently unused.
#' @return Invisibly returns `x`.
#' @seealso [summary.causatr_diag()], [diagnose()]
#' @export
print.causatr_diag <- function(x, ...) {
  cat("<causatr_diag>\n", " Method:", x$method, "\n\n", sep = "")
  # Each section is conditionally printed based on which diagnostics
  # the underlying method supports:
  #   - positivity: any method (when treatment is binary)
  #   - balance:    always (cobalt table or simple SMD fallback)
  #   - weights:    IPW only (NULL elsewhere)
  #   - match_quality: matching only (NULL elsewhere)
  # The NULL checks let us print one compact block per method
  # without a switch() on method type.
  if (!is.null(x$positivity)) {
    cat("Positivity (propensity score):\n")
    print(x$positivity, row.names = FALSE)
    cat("\n")
  }
  if (!is.null(x$balance)) {
    cat("Covariate balance:\n")
    print(x$balance)
    cat("\n")
  }
  if (!is.null(x$weights)) {
    cat("Weight distribution:\n")
    print(x$weights, row.names = FALSE)
    cat("\n")
  }
  if (!is.null(x$match_quality)) {
    cat("Match quality:\n")
    print(x$match_quality, row.names = FALSE)
    cat("\n")
  }
  invisible(x)
}
