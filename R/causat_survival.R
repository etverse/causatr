#' Fit a causal survival model
#'
#' @description
#' Convenience wrapper for causal survival analysis using pooled logistic
#' regression as a discrete-time hazard model. Converts the data to
#' person-period format if needed, fits one logistic model per time point
#' (or a single model with a flexible time function), and returns a
#' `causatr_fit` that [contrast()] uses to compute survival curves and risk
#' differences under each intervention.
#'
#' Censoring is handled within the model (restrict to uncensored at each
#' time) or via inverse probability of censoring weighting (IPCW) supplied
#' through the `weights` argument. For competing risks, a separate
#' cause-specific hazard model is fitted per event type.
#'
#' @param data A data frame or data.table in either wide format (one row per
#'   individual, with separate time-to-event and event indicator columns) or
#'   long (person-period) format.
#' @param outcome Character. Name of the binary event indicator (1 = event,
#'   0 = censored or no event at that interval).
#' @param treatment Character. Name of the treatment variable.
#' @param confounders A one-sided formula specifying confounders.
#' @param id Character. Name of the individual ID variable.
#' @param time Character. Name of the time variable (interval index).
#' @param censoring Character or `NULL`. Name of the censoring indicator.
#' @param competing Character or `NULL`. Name of a variable indicating the
#'   type of competing event (for competing risks analysis).
#' @param time_formula A one-sided formula specifying how time enters the
#'   hazard model. Default `~ splines::ns(time, 4)`. Use `~ factor(time)`
#'   for a fully saturated (non-parametric) baseline hazard.
#' @param weights Numeric vector or `NULL`. Pre-computed IPCW or survey
#'   weights.
#' @param ... Additional arguments passed to `glm()`.
#'
#' @return A `causatr_fit` object (with `type = "survival"`) suitable for
#'   use with [contrast()]. When passed to `contrast()`, estimates are
#'   survival curves at each time point and risk/survival differences at
#'   user-specified times.
#'
#' @references
#' Hernán MA, Robins JM (2025). *Causal Inference: What If*. Chapman &
#' Hall/CRC. Chapter 17.
#'
#' @examples
#' \dontrun{
#' data("nhefs", package = "causatr")
#' fit_surv <- causat_survival(
#'   nhefs,
#'   outcome = "death",
#'   treatment = "qsmk",
#'   confounders = ~ sex + age + race + education +
#'     smokeintensity + smokeyrs + exercise + active + wt71,
#'   id = "seqn",
#'   time = "year"
#' )
#' result <- contrast(fit_surv,
#'   interventions = list(quit = static(1), continue = static(0)),
#'   type = "difference"
#' )
#' }
#'
#' @seealso [causat()], [contrast()], [to_person_period()]
#' @export
causat_survival <- function(
  data,
  outcome,
  treatment,
  confounders,
  id,
  time,
  censoring = NULL,
  competing = NULL,
  time_formula = ~ splines::ns(time, 4),
  weights = NULL,
  ...
) {
  rlang::abort("causat_survival() is not yet implemented.", .call = FALSE)
}
