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

#' @export
print.causatr_diag <- function(x, ...) {
  cat("<causatr_diag>\n Method:", x$method, "\n\n")
  if (!is.null(x$positivity)) {
    cat("Positivity:\n")
    print(x$positivity)
  }
  if (!is.null(x$balance)) {
    cat("\nBalance:\n")
    print(x$balance)
  }
  if (!is.null(x$weights)) {
    cat("\nWeight distribution:\n")
    print(x$weights)
  }
  if (!is.null(x$match_quality)) {
    cat("\nMatch quality:\n")
    print(x$match_quality)
  }
  invisible(x)
}
