#' Print a causatr fit
#'
#' @description
#' Displays a compact summary of a [causatr_fit][causat] object, showing the
#' causal estimator, treatment type, outcome, treatment variable, and sample
#' size.
#'
#' @param x A `causatr_fit` object.
#' @param ... Currently unused.
#' @return Invisibly returns `x`.
#' @seealso [summary.causatr_fit()], [causat()]
#' @export
print.causatr_fit <- function(x, ...) {
  # Pretty-print the estimator. Defaults-through the raw value if the
  # user has somehow set an unrecognized estimator so print() never
  # errors on a malformed fit object.
  estimator_label <- switch(
    x$estimator,
    gcomp = "G-computation",
    ipw = "IPW (WeightIt)",
    matching = "Matching (MatchIt)",
    x$estimator
  )
  family_label <- format_family(x$family)
  # Collapse multivariate treatments into a comma-separated label.
  trt_label <- paste(x$treatment, collapse = ", ")

  cat("<causatr_fit>\n")
  cat(" Estimator:  ", estimator_label, "\n", sep = "")
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
#' Displays the causal estimator, contrast type, CI method, sample size,
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
  estimator_label <- switch(
    x$estimator,
    gcomp = "G-computation",
    ipw = "IPW",
    matching = "Matching",
    x$estimator
  )
  type_label <- switch(
    x$type,
    difference = "Difference",
    ratio = "Ratio",
    or = "Odds ratio"
  )

  cat("<causatr_result>\n")
  cat(" Estimator: ", estimator_label, "\n", sep = "")
  cat(" Estimand:  ", x$estimand, "\n", sep = "")
  cat(" Contrast:  ", type_label, "\n", sep = "")
  cat(" CI method: ", x$ci_method, "\n", sep = "")
  cat(" N:         ", x$n, "\n", sep = "")

  # Bootstrap diagnostics: surface the success/failure split so users
  # see at print time when replicates were silently discarded. The
  # `boot_info` slot is only populated for `ci_method = "bootstrap"`,
  # and may be a single 3-element list (no `by`) or a list-of-lists
  # keyed by subgroup level when `by = ...` was used. In the by case
  # we reduce to totals across strata; per-stratum detail is available
  # in `x$boot_info` for users who want it.
  if (!is.null(x$boot_info)) {
    bi <- x$boot_info
    if (
      is.list(bi) && !all(c("n_requested", "n_ok", "n_fail") %in% names(bi))
    ) {
      # by-stratum list-of-lists: collapse to totals.
      n_req <- sum(vapply(bi, function(b) b$n_requested %||% 0L, integer(1)))
      n_fail <- sum(vapply(bi, function(b) b$n_fail %||% 0L, integer(1)))
    } else {
      n_req <- bi$n_requested
      n_fail <- bi$n_fail
    }
    if (!is.null(n_req) && n_req > 0L) {
      n_ok <- n_req - n_fail
      pct <- if (n_fail > 0L) {
        paste0(" (", round(100 * n_fail / n_req, 1), "% failed)")
      } else {
        ""
      }
      cat(
        " Bootstrap: ",
        n_ok,
        "/",
        n_req,
        " replicates",
        pct,
        "\n",
        sep = ""
      )
    }
  }

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
  cat("<causatr_diag>\n", " Estimator:", x$estimator, "\n\n", sep = "")
  # Each section is conditionally printed based on which diagnostics
  # the underlying estimator supports:
  #   - positivity: any estimator (when treatment is binary)
  #   - balance:    always (cobalt table or simple SMD fallback)
  #   - weights:    IPW only (NULL elsewhere)
  #   - match_quality: matching only (NULL elsewhere)
  # The NULL checks let us print one compact block per estimator
  # without a switch() on estimator type.
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
