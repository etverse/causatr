#' @noRd
fit_ice <- function(
  data,
  outcome,
  treatment,
  confounders,
  family,
  weights,
  call,
  ...
) {
  rlang::abort(
    "ICE g-computation for longitudinal treatments is not yet implemented.",
    .call = FALSE
  )
}
