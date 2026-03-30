#' @noRd
fit_ipw <- function(
  data,
  outcome,
  treatment,
  confounders,
  type,
  weights,
  call,
  ...
) {
  rlang::abort(
    "IPW estimation is not yet implemented.",
    .call = FALSE
  )
}
