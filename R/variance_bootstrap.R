#' @noRd
variance_bootstrap <- function(
  fit,
  interventions,
  estimates,
  n_boot,
  conf_level,
  by
) {
  rlang::abort(
    "Bootstrap variance estimation is not yet implemented.",
    .call = FALSE
  )
}
