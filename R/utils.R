#' @noRd
new_causatr_fit <- function(
  model,
  data,
  treatment,
  outcome,
  confounders,
  family,
  method,
  type,
  id,
  time,
  censoring,
  weights_obj,
  match_obj,
  call,
  details
) {
  structure(
    list(
      model = model,
      data = data,
      treatment = treatment,
      outcome = outcome,
      confounders = confounders,
      family = family,
      method = method,
      type = type,
      id = id,
      time = time,
      censoring = censoring,
      weights_obj = weights_obj,
      match_obj = match_obj,
      call = call,
      details = details
    ),
    class = "causatr_fit"
  )
}

#' @noRd
new_causatr_result <- function(
  estimates,
  contrasts,
  type,
  ci_method,
  interventions,
  n,
  method,
  vcov,
  call
) {
  structure(
    list(
      estimates = estimates,
      contrasts = contrasts,
      type = type,
      ci_method = ci_method,
      interventions = interventions,
      n = n,
      method = method,
      vcov = vcov,
      call = call
    ),
    class = "causatr_result"
  )
}

#' @noRd
new_causatr_diag <- function(
  balance,
  positivity,
  weights,
  match_quality,
  method
) {
  structure(
    list(
      balance = balance,
      positivity = positivity,
      weights = weights,
      match_quality = match_quality,
      method = method
    ),
    class = "causatr_diag"
  )
}

#' @noRd
new_causatr_intervention <- function(type, params) {
  structure(
    c(list(type = type), params),
    class = "causatr_intervention"
  )
}

#' @noRd
check_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    rlang::abort(
      paste0(
        "Package '",
        pkg,
        "' is required but not installed. ",
        "Install it with: install.packages('",
        pkg,
        "')"
      ),
      .call = FALSE
    )
  }
}
