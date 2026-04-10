#' Construct a `causatr_fit` object
#'
#' @param model Fitted model object (glm, gam, glm_weightit, etc.) or `NULL`
#'   (ICE defers fitting to `contrast()`).
#' @param data data.table of the full dataset.
#' @param treatment Character treatment column name(s).
#' @param outcome Character outcome column name.
#' @param confounders One-sided formula of baseline confounders.
#' @param confounders_tv One-sided formula of time-varying confounders or `NULL`.
#' @param family Character or family object describing the outcome distribution.
#' @param method Character estimation method (`"gcomp"`, `"ipw"`, `"matching"`).
#' @param type `"point"` or `"longitudinal"`.
#' @param estimand Character estimand (`"ATE"`, `"ATT"`, `"ATC"`).
#' @param id Character ID column name or `NULL`.
#' @param time Character time column name or `NULL`.
#' @param censoring Character censoring column name or `NULL`.
#' @param history Integer Markov order for longitudinal.
#' @param numerator One-sided formula for stabilised weights or `NULL`.
#' @param weights_obj A `weightit` object (IPW) or `NULL`.
#' @param match_obj A `matchit` object (matching) or `NULL`.
#' @param call The original `causat()` call environment.
#' @param details Named list of method-specific metadata.
#' @return A list with class `"causatr_fit"`.
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

#' Construct a `causatr_result` object
#'
#' @param estimates data.table of intervention-specific marginal means.
#' @param contrasts data.table of pairwise contrasts with SEs and CIs.
#' @param type Character contrast type (`"difference"`, `"ratio"`, `"or"`).
#' @param estimand Character estimand used.
#' @param ci_method Character CI method (`"sandwich"` or `"bootstrap"`).
#' @param reference Character name of the reference intervention or `NULL`.
#' @param interventions Named list of `causatr_intervention` objects.
#' @param n Integer sample size used for estimation.
#' @param method Character estimation method.
#' @param vcov Variance-covariance matrix of marginal means.
#' @param call The original `contrast()` call environment.
#' @return A list with class `"causatr_result"`.
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
  family,
  fit_type,
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
      family = family,
      fit_type = fit_type,
      vcov = vcov,
      call = call
    ),
    class = "causatr_result"
  )
}

#' Construct a `causatr_diag` object
#'
#' @param balance Balance table (from cobalt or simple SMD computation).
#' @param positivity data.table of propensity score summaries.
#' @param weights data.table of weight distribution summaries (IPW) or `NULL`.
#' @param match_quality data.table of match quality metrics or `NULL`.
#' @param method Character estimation method.
#' @param fit The original `causatr_fit` (stored for `plot()` method).
#' @return A list with class `"causatr_diag"`.
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

#' Construct a `causatr_intervention` object
#'
#' @param type Character intervention type (e.g. `"static"`, `"shift"`).
#' @param params Named list of intervention parameters.
#' @return A list with class `"causatr_intervention"`.
#' @noRd
new_causatr_intervention <- function(type, params) {
  structure(
    c(list(type = type), params),
    class = "causatr_intervention"
  )
}

#' Check whether rows are uncensored
#'
#' A row is uncensored when the censoring column is `NA` or `0`.
#' Any other non-NA value (e.g. `1` for censored, `2` for a competing event)
#' is treated as censored.
#'
#' @param data A data.table.
#' @param censoring Character censoring column name, or `NULL`.
#' @return Logical vector of length `nrow(data)` (`TRUE` = uncensored).
#' @noRd
is_uncensored <- function(data, censoring) {
  if (is.null(censoring)) {
    return(rep(TRUE, nrow(data)))
  }
  cens <- data[[censoring]]
  is.na(cens) | cens == 0L
}

#' Resolve a family argument to a family object
#'
#' Accepts a character string (e.g. `"gaussian"`), a family function
#' (e.g. `stats::binomial`), or an already-evaluated family object
#' (e.g. `stats::binomial()`).
#'
#' @param family Character, function, or family object.
#' @return A family object (list with `$family`, `$link`, etc.).
#' @noRd
resolve_family <- function(family) {
  if (is.character(family)) {
    return(get(family, mode = "function", envir = asNamespace("stats"))())
  }
  if (is.function(family)) {
    return(family())
  }
  family
}

#' Check whether a family describes a binary outcome
#'
#' @param family A character string, family object, or function.
#' @return Logical scalar.
#' @noRd
is_binary_family <- function(family) {
  if (is.null(family)) {
    return(FALSE)
  }
  fam <- tryCatch(resolve_family(family), error = function(e) NULL)
  if (!is.null(fam) && is.character(fam$family)) {
    return(fam$family %in% c("binomial", "quasibinomial"))
  }
  FALSE
}

#' Check that an optional package is installed
#'
#' @param pkg Character package name.
#' @return `NULL` invisibly; aborts with installation instructions if absent.
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
