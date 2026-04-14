#' Summarise a causatr fit
#'
#' @description
#' Provides a detailed summary of a [causatr_fit][causat] object, including
#' the causal estimator, outcome family, confounder formula, and the
#' underlying model (outcome model, propensity weights, or matching object).
#'
#' @param object A `causatr_fit` object.
#' @param ... Currently unused.
#' @return Invisibly returns `object`.
#' @seealso [print.causatr_fit()], [causat()]
#' @export
summary.causatr_fit <- function(object, ...) {
  cat(
    "causatr fit\n",
    "Estimator:   ",
    object$estimator,
    "\n",
    "Type:        ",
    object$type,
    "\n",
    "Outcome:     ",
    object$outcome,
    " (",
    object$family,
    ")\n",
    "Treatment:   ",
    object$treatment,
    "\n",
    "Confounders: ",
    deparse(object$confounders),
    "\n",
    "N:           ",
    nrow(object$data),
    "\n",
    sep = ""
  )
  if (!is.null(object$id)) {
    cat(" ID:         ", object$id, "\n", sep = "")
    cat(" Time:       ", object$time, "\n", sep = "")
  }
  if (!is.null(object$model)) {
    cat("\nOutcome model:\n")
    print(summary(object$model))
  }
  if (!is.null(object$weights_obj)) {
    cat("\nPropensity weights:\n")
    print(summary(object$weights_obj))
  }
  if (!is.null(object$match_obj)) {
    cat("\nMatching:\n")
    print(summary(object$match_obj))
  }
  invisible(object)
}

#' Summarise a causatr result
#'
#' @description
#' Displays intervention-specific marginal means and pairwise contrasts with
#' standard errors and confidence intervals from a [causatr_result][contrast]
#' object.
#'
#' @param object A `causatr_result` object.
#' @param ... Currently unused.
#' @return Invisibly returns `object`.
#' @seealso [print.causatr_result()], [contrast()]
#' @export
summary.causatr_result <- function(object, ...) {
  # Start with the standard print() header + estimates/contrasts tables.
  print(object)

  # Then append per-intervention detail lines. Three shapes to
  # handle, matching the constructors in R/interventions.R:
  #   - NULL                 -> "natural course (NULL)"
  #   - causatr_intervention -> "type, p1 = v1, p2 = v2, ..."
  # Multivariate treatment interventions (named list of sub-interventions)
  # are not expanded here — they'd clutter the output; users can
  # inspect `object$interventions` directly for that detail.
  cat("\nIntervention details:\n")
  for (nm in names(object$interventions)) {
    iv <- object$interventions[[nm]]
    if (is.null(iv)) {
      cat("  ", nm, ": natural course (NULL)\n", sep = "")
    } else if (inherits(iv, "causatr_intervention")) {
      cat("  ", nm, ": ", iv$type, sep = "")
      # Drop the type slot (already printed) and list remaining
      # params. Function-valued params (dynamic rule closures) get
      # a placeholder rather than a cryptic environment dump.
      params <- iv[names(iv) != "type"]
      for (p in names(params)) {
        if (!is.function(params[[p]])) {
          cat(", ", p, " = ", params[[p]], sep = "")
        } else {
          cat(", ", p, " = <function>", sep = "")
        }
      }
      cat("\n")
    }
  }
  # Finally, the full vcov matrix — useful for hand-computing custom
  # linear contrasts beyond the pairwise ones in `object$contrasts`.
  cat("\nVariance-covariance matrix of marginal means:\n")
  print(object$vcov, digits = 6)
  invisible(object)
}

#' Summarise causatr diagnostics
#'
#' @description
#' Prints positivity, balance, weight, and match quality diagnostics from a
#' [causatr_diag][diagnose] object.
#'
#' @param object A `causatr_diag` object.
#' @param ... Currently unused.
#' @return Invisibly returns `object`.
#' @seealso [print.causatr_diag()], [diagnose()]
#' @export
summary.causatr_diag <- function(object, ...) {
  print(object)
  if (!is.null(object$fit)) {
    cat("Underlying fit:\n")
    print(object$fit)
  }
  invisible(object)
}
