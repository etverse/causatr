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
  check_pkg("WeightIt")

  if (type == "longitudinal") {
    rlang::abort(
      paste0(
        "Longitudinal IPW is not yet implemented. ",
        "For time-varying treatments, use method = 'gcomp' (ICE g-computation)."
      ),
      .call = FALSE
    )
  }

  rlang::abort(
    "IPW estimation is not yet implemented.",
    .call = FALSE
  )
}
