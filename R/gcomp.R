#' @noRd
fit_gcomp <- function(
  data,
  outcome,
  treatment,
  confounders,
  confounders_tv,
  family,
  estimand,
  type,
  history,
  weights,
  call,
  ...
) {
  if (type == "longitudinal") {
    fit_ice(
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
    )
  } else {
    fit_gcomp_point(
      data,
      outcome,
      treatment,
      confounders,
      family,
      estimand,
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
  estimand,
  weights,
  call,
  ...
) {
  rlang::abort(
    "G-computation for point treatments is not yet implemented.",
    .call = FALSE
  )
}
