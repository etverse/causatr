#' @noRd
fit_ice <- function(
  data,
  outcome,
  treatment,
  confounders,
  confounders_tv,
  family,
  estimand,
  history,
  weights,
  call,
  ...
) {
  rlang::abort(
    "ICE g-computation for longitudinal treatments is not yet implemented.",
    .call = FALSE
  )
}
