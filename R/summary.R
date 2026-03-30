#' @export
summary.causatr_fit <- function(object, ...) {
  cat(
    "causatr fit\n",
    "Method:      ",
    object$method,
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

#' @export
summary.causatr_result <- function(object, ...) {
  cat(
    "causatr result\n",
    "Method:    ",
    object$method,
    "\n",
    "Contrast:  ",
    object$type,
    "\n",
    "CI method: ",
    object$ci_method,
    "\n",
    "N:         ",
    object$n,
    "\n\n",
    sep = ""
  )
  cat("Intervention means:\n")
  print(object$estimates)
  cat("\nContrasts:\n")
  print(object$contrasts)
  invisible(object)
}

#' @export
summary.causatr_diag <- function(object, ...) {
  print(object)
  invisible(object)
}
