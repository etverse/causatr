#' Package load hook: register knit_print S3 method
#'
#' @description
#' knitr is a Suggests (not Imports) dependency, so we cannot use the
#' `@method` roxygen tag + `NAMESPACE` export path to register
#' `knit_print.causatr_result()` — that would force an unconditional
#' dependency on knitr. Instead, we register it at load time if and only
#' if knitr is actually installed. This is the standard tidyverse pattern
#' for Suggests-only S3 methods (see, e.g., `vctrs:::s3_register()`).
#'
#' @keywords internal
#' @noRd
.onLoad <- function(libname, pkgname) {
  # Only register if knitr is available — otherwise users without knitr
  # installed would get an error on package load.
  if (requireNamespace("knitr", quietly = TRUE)) {
    registerS3method(
      "knit_print",
      "causatr_result",
      knit_print.causatr_result,
      envir = asNamespace("knitr")
    )
  }
  invisible()
}
