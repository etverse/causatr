#' @noRd
new_causatr_fit <- function(
  model,
  data,
  treatment,
  outcome,
  confounders,
  confounders_tv,
  family,
  method,
  type,
  estimand,
  id,
  time,
  censoring,
  history,
  numerator,
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
      confounders_tv = confounders_tv,
      family = family,
      method = method,
      type = type,
      estimand = estimand,
      id = id,
      time = time,
      censoring = censoring,
      history = history,
      numerator = numerator,
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
  estimand,
  ci_method,
  reference,
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
      estimand = estimand,
      ci_method = ci_method,
      reference = reference,
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
  method,
  fit = NULL
) {
  structure(
    list(
      balance = balance,
      positivity = positivity,
      weights = weights,
      match_quality = match_quality,
      method = method,
      fit = fit
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
