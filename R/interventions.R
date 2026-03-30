#' Set treatment to a fixed value
#'
#' @description
#' Creates a static intervention that sets the treatment to a fixed value for
#' all individuals at all time points. The most common intervention type,
#' corresponding to "always treat" (`static(1)`) or "never treat" (`static(0)`).
#'
#' @param value The fixed treatment value.
#'
#' @return A `causatr_intervention` object.
#'
#' @references
#' Hernán MA, Robins JM (2025). *Causal Inference: What If*. Chapman &
#' Hall/CRC. Chapter 1 (static interventions).
#'
#' @examples
#' \dontrun{
#' data("nhefs", package = "causatr")
#' fit <- causat(nhefs, outcome = "wt82_71", treatment = "qsmk",
#'               confounders = ~ sex + age + wt71)
#' contrast(fit, interventions = list(quit = static(1), continue = static(0)))
#' }
#'
#' @seealso [shift()], [dynamic()], [scale()], [threshold()], [ipsi()],
#'   [contrast()]
#' @export
static <- function(value) {
  new_causatr_intervention("static", list(value = value))
}

#' Shift treatment by a fixed amount
#'
#' @description
#' Creates a modified treatment policy (MTP) that adds a fixed `delta` to each
#' individual's observed treatment value. Useful for continuous treatments where
#' a population-level shift is the relevant intervention (e.g., "reduce
#' exposure by 10 units").
#'
#' @param delta Numeric. The amount to add to the observed treatment.
#'
#' @return A `causatr_intervention` object.
#'
#' @references
#' Díaz I, Williams N, Hoffman KL, Schenck EJ (2023). Non-parametric causal
#' effects based on longitudinal modified treatment policies. *Journal of the
#' American Statistical Association* 118:846–857.
#'
#' @examples
#' \dontrun{
#' fit <- causat(nhefs, outcome = "wt82_71", treatment = "smokeintensity",
#'               confounders = ~ sex + age + wt71)
#' contrast(fit, interventions = list(
#'   reduce10 = shift(-10),
#'   observed = shift(0)
#' ))
#' }
#'
#' @seealso [static()], [scale()], [threshold()], [dynamic()], [ipsi()]
#' @export
shift <- function(delta) {
  new_causatr_intervention("shift", list(delta = delta))
}

#' Multiply treatment by a fixed factor
#'
#' @description
#' Creates a modified treatment policy that multiplies each individual's
#' observed treatment by `factor`. Useful for proportional reductions or
#' increases in continuous treatments.
#'
#' @param factor Numeric. The multiplicative factor.
#'
#' @return A `causatr_intervention` object.
#'
#' @examples
#' \dontrun{
#' fit <- causat(nhefs, outcome = "wt82_71", treatment = "smokeintensity",
#'               confounders = ~ sex + age + wt71)
#' contrast(fit, interventions = list(
#'   halved = scale(0.5),
#'   observed = scale(1)
#' ))
#' }
#'
#' @seealso [static()], [shift()], [threshold()], [dynamic()], [ipsi()]
#' @export
scale <- function(factor) {
  new_causatr_intervention("scale", list(factor = factor))
}

#' Clamp treatment within bounds
#'
#' @description
#' Creates a modified treatment policy that clamps each individual's observed
#' treatment to lie within `[lower, upper]`. Values below `lower` are set to
#' `lower`; values above `upper` are set to `upper`.
#'
#' @param lower Numeric. Lower bound (use `-Inf` for no lower bound).
#' @param upper Numeric. Upper bound (use `Inf` for no upper bound).
#'
#' @return A `causatr_intervention` object.
#'
#' @examples
#' \dontrun{
#' fit <- causat(nhefs, outcome = "wt82_71", treatment = "smokeintensity",
#'               confounders = ~ sex + age + wt71)
#' contrast(fit, interventions = list(
#'   capped20 = threshold(0, 20),
#'   observed  = shift(0)
#' ))
#' }
#'
#' @seealso [static()], [shift()], [scale()], [dynamic()], [ipsi()]
#' @export
threshold <- function(lower = -Inf, upper = Inf) {
  new_causatr_intervention("threshold", list(lower = lower, upper = upper))
}

#' Dynamic treatment rule
#'
#' @description
#' Creates a dynamic intervention where the treatment at each time point is
#' determined by a user-supplied function of the covariate history. The
#' function `rule` receives the current data (subset to the current time point)
#' and the observed treatment vector, and returns the intervened treatment
#' vector.
#'
#' @param rule A function with signature `function(data, treatment)` that
#'   returns a vector of treatment values of the same length as `nrow(data)`.
#'
#' @return A `causatr_intervention` object.
#'
#' @references
#' Hernán MA, Robins JM (2025). *Causal Inference: What If*. Chapman &
#' Hall/CRC. Chapter 19 (dynamic treatment strategies).
#'
#' @examples
#' # Treat if CD4 count is below 200
#' cd4_rule <- dynamic(\(data, trt) ifelse(data$cd4 < 200, 1, 0))
#'
#' @seealso [static()], [shift()], [scale()], [threshold()], [ipsi()]
#' @export
dynamic <- function(rule) {
  if (!is.function(rule)) {
    rlang::abort(
      "`rule` must be a function with signature function(data, treatment)."
    )
  }
  new_causatr_intervention("dynamic", list(rule = rule))
}

#' Incremental propensity score intervention
#'
#' @description
#' Creates an incremental propensity score intervention (IPSI) that multiplies
#' each individual's odds of treatment by `delta`. Values of `delta > 1`
#' increase the probability of treatment; `delta < 1` decrease it. This is a
#' stochastic modified treatment policy indexed by a single scalar.
#'
#' @param delta Positive numeric. The odds multiplier.
#'
#' @return A `causatr_intervention` object.
#'
#' @references
#' Kennedy EH (2019). Nonparametric causal effects based on incremental
#' propensity score interventions. *Journal of the American Statistical
#' Association* 114:645–656.
#'
#' @examples
#' \dontrun{
#' fit <- causat(nhefs, outcome = "wt82_71", treatment = "qsmk",
#'               confounders = ~ sex + age + wt71)
#' contrast(fit, interventions = list(
#'   double_odds = ipsi(2),
#'   half_odds   = ipsi(0.5)
#' ))
#' }
#'
#' @seealso [static()], [shift()], [dynamic()], [scale()], [threshold()]
#' @export
ipsi <- function(delta) {
  if (!rlang::is_scalar_double(delta) && !rlang::is_scalar_integer(delta)) {
    rlang::abort("`delta` must be a single positive number.")
  }
  if (delta <= 0) {
    rlang::abort("`delta` must be positive.")
  }
  new_causatr_intervention("ipsi", list(delta = delta))
}

#' @export
print.causatr_intervention <- function(x, ...) {
  cat("<causatr_intervention: ", x$type, ">\n", sep = "")
  params <- x[names(x) != "type"]
  for (nm in names(params)) {
    if (is.function(params[[nm]])) {
      cat(" ", nm, ": <function>\n", sep = "")
    } else {
      cat(" ", nm, ": ", params[[nm]], "\n", sep = "")
    }
  }
  invisible(x)
}
