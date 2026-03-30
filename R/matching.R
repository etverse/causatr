#' @noRd
fit_matching <- function(
  data,
  outcome,
  treatment,
  confounders,
  type,
  weights,
  call,
  ...
) {
  check_pkg("MatchIt")

  if (type == "longitudinal") {
    rlang::abort(
      "Matching is only supported for point treatments.",
      .call = FALSE
    )
  }

  rlang::abort(
    "Matching estimation is not yet implemented.",
    .call = FALSE
  )
}
