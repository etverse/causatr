#' Print a causatr fit
#'
#' @param x A `causatr_fit` object.
#' @param ... Currently unused.
#' @return Invisibly returns `x`.
#' @seealso [summary.causatr_fit()], [causat()]
#' @export
print.causatr_fit <- function(x, ...) {
  cat(
    "<causatr_fit>\n",
    " Method:   ",
    x$method,
    "\n",
    " Type:     ",
    x$type,
    "\n",
    " Outcome:  ",
    x$outcome,
    "\n",
    " Treatment:",
    x$treatment,
    "\n",
    " N:        ",
    nrow(x$data),
    "\n",
    sep = ""
  )
  invisible(x)
}

#' Print a causatr result
#'
#' @param x A `causatr_result` object.
#' @param ... Currently unused.
#' @return Invisibly returns `x`.
#' @seealso [summary.causatr_result()], [contrast()]
#' @export
print.causatr_result <- function(x, ...) {
  cat(
    "<causatr_result>\n",
    " Method:      ",
    x$method,
    "\n",
    " Contrast:    ",
    x$type,
    "\n",
    " CI method:   ",
    x$ci_method,
    "\n",
    " N:           ",
    x$n,
    "\n\n",
    sep = ""
  )
  cat("Intervention means:\n")
  print(x$estimates, digits = 4)
  cat("\nContrasts:\n")
  print(x$contrasts, digits = 4)
  invisible(x)
}

#' Print causatr diagnostics
#'
#' @param x A `causatr_diag` object.
#' @param ... Currently unused.
#' @return Invisibly returns `x`.
#' @seealso [summary.causatr_diag()], [diagnose()]
#' @export
print.causatr_diag <- function(x, ...) {
  cat("<causatr_diag>\n", " Method:", x$method, "\n\n", sep = "")
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
