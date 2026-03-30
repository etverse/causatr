#' @noRd
fit_gcomp <- function(
  data,
  outcome,
  treatment,
  confounders,
  family,
  type,
  weights,
  call,
  ...
) {
  if (type == "longitudinal") {
    fit_ice(data, outcome, treatment, confounders, family, weights, call, ...)
  } else {
    fit_gcomp_point(
      data,
      outcome,
      treatment,
      confounders,
      family,
      weights,
      call,
      ...
    )
  }
}

#' @noRd
fit_gcomp_point <- function(
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
    "G-computation for point treatments is not yet implemented.",
    .call = FALSE
  )
}
