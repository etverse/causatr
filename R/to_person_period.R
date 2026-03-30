#' Convert wide data to person-period (long) format
#'
#' @description
#' Reshapes individual-level data in wide format (one row per person, with
#' columns for each time-varying variable at each time point) into long
#' person-period format (one row per person per time interval), as required
#' by [causat()] for longitudinal analyses and by [causat_survival()] for
#' survival analyses.
#'
#' @param data A data frame or data.table in wide format.
#' @param id Character. Name of the individual ID variable.
#' @param time_varying Named list specifying how to reshape time-varying
#'   columns. Each element is a character vector of column names for one
#'   variable across time, in time order. For example:
#'   `list(A = c("A0", "A1", "A2"), L = c("L0", "L1", "L2"))`.
#' @param time_invariant Character vector of column names to carry forward
#'   unchanged (baseline covariates, outcome, etc.).
#' @param time_name Character. Name for the new time column. Default `"time"`.
#'
#' @return A data.table in person-period format with one row per individual
#'   per time point, sorted by `id` and `time_name`.
#'
#' @examples
#' \dontrun{
#' # Toy wide dataset
#' wide <- data.table::data.table(
#'   id  = 1:3,
#'   sex = c(0, 1, 0),
#'   A0  = c(1, 0, 1),
#'   A1  = c(1, 1, 0),
#'   L0  = c(5, 3, 7),
#'   L1  = c(4, 6, 8),
#'   Y   = c(0, 1, 0)
#' )
#' long <- to_person_period(
#'   wide,
#'   id = "id",
#'   time_varying = list(A = c("A0", "A1"), L = c("L0", "L1")),
#'   time_invariant = c("sex", "Y")
#' )
#' }
#'
#' @seealso [causat()], [causat_survival()]
#' @export
to_person_period <- function(
  data,
  id,
  time_varying,
  time_invariant = character(0),
  time_name = "time"
) {
  rlang::abort("to_person_period() is not yet implemented.", .call = FALSE)
}
