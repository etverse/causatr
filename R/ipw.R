#' @noRd
fit_ipw <- function(
  data,
  outcome,
  treatment,
  confounders,
  confounders_tv,
  estimand,
  type,
  history,
  numerator,
  weights,
  call,
  ...
) {
  rlang::abort(
    "IPW estimation is not yet implemented.",
    .call = FALSE
  )
}
