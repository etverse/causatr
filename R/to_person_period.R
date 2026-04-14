#' Convert wide data to person-period (long) format
#'
#' @description
#' Reshapes individual-level data in wide format (one row per person, with
#' columns for each time-varying variable at each time point) into long
#' person-period format (one row per person per time interval), as required
#' by [causat()] for longitudinal analyses and by [causat_survival()] for
#' discrete-time survival analyses.
#'
#' For survival data the event indicator \eqn{D_{k+1} = 1} means the individual
#' experienced the event during interval k+1 (Hernán & Robins Ch. 17).
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
  if (!data.table::is.data.table(data)) {
    data <- data.table::as.data.table(data)
  }

  all_tv_cols <- unlist(time_varying, use.names = FALSE)
  missing_cols <- setdiff(all_tv_cols, names(data))
  if (length(missing_cols) > 0L) {
    rlang::abort(
      paste0(
        "Column(s) not found in data: ",
        paste(missing_cols, collapse = ", ")
      )
    )
  }

  n_times <- length(time_varying[[1]])

  # Validate that each id appears exactly once in the wide input. The
  # whole reshape assumes one row per individual (since each (var, t)
  # cell is a distinct column), so a duplicated id in the wide table
  # silently triggers an n*n_times-row stack with mangled time
  # alignment. Catch this up front instead of letting downstream
  # ICE / survival code see malformed long data.
  id_counts <- data[, .N, by = c(id)]
  dup_ids <- id_counts[N > 1L][[id]]
  if (length(dup_ids) > 0L) {
    rlang::abort(
      paste0(
        "`to_person_period()` requires each id to appear exactly once in ",
        "the wide input. Found ",
        length(dup_ids),
        " duplicated id(s): ",
        paste(utils::head(dup_ids, 5), collapse = ", "),
        if (length(dup_ids) > 5L) ", ..." else "",
        "."
      ),
      .call = FALSE
    )
  }

  # Validate that all time-varying variables have the same number of columns.
  for (nm in names(time_varying)) {
    if (length(time_varying[[nm]]) != n_times) {
      rlang::abort(
        paste0(
          "All elements of `time_varying` must have the same length. '",
          nm,
          "' has ",
          length(time_varying[[nm]]),
          " columns but expected ",
          n_times,
          "."
        )
      )
    }
  }

  # Build the long table by stacking n_times copies of the data,
  # one per time point. This is the simple/explicit reshape idiom —
  # clearer than `data.table::melt()` when every variable has its
  # own column-name pattern (the `time_varying` list specifies
  # exactly which original column maps to which time point).
  #
  # For each time point k (0-indexed):
  #   1. Start from a slice containing id + baseline columns.
  #   2. Tag the slice with the time column.
  #   3. Pull each time-varying variable's k-th column into a
  #      canonical name (e.g. "A0" -> "A", "A1" -> "A" at t=1).
  #
  # The result is one row per (id, time) with a consistent column
  # schema — exactly the shape that prepare_data() and fit_ice()
  # expect.
  long_list <- lapply(seq_len(n_times), function(k) {
    row_dt <- data[, c(id, time_invariant), with = FALSE]
    row_dt[, (time_name) := k - 1L] # 0-indexed time

    for (nm in names(time_varying)) {
      col_name <- time_varying[[nm]][k]
      row_dt[, (nm) := data[[col_name]]]
    }

    row_dt
  })

  result <- data.table::rbindlist(long_list)
  data.table::setkeyv(result, c(id, time_name))
  result
}
