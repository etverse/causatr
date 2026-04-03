#' Validate inputs and prepare data for estimation
#'
#' Coerces to data.table, selects relevant columns, and creates lag variables
#' for longitudinal data.
#'
#' @param data A data.frame or data.table.
#' @param outcome Character outcome column name.
#' @param treatment Character treatment column name(s).
#' @param confounders One-sided formula of baseline confounders.
#' @param confounders_tv One-sided formula of time-varying confounders or `NULL`.
#' @param id Character ID column name or `NULL`.
#' @param time Character time column name or `NULL`.
#' @param censoring Character censoring column name or `NULL`.
#' @param history Positive integer Markov order or `Inf`.
#' @param call Caller environment for error messages.
#' @return A data.table with only the needed columns (and lag columns for
#'   longitudinal data).
#' @noRd
prepare_data <- function(
  data,
  outcome,
  treatment,
  confounders,
  confounders_tv = NULL,
  id = NULL,
  time = NULL,
  censoring = NULL,
  history = 1L,
  call = rlang::caller_env()
) {
  if (!data.table::is.data.table(data)) {
    data <- data.table::as.data.table(data)
  }

  tv_vars <- if (!is.null(confounders_tv)) {
    all.vars(confounders_tv)
  } else {
    character(0)
  }

  keep_cols <- unique(c(
    outcome,
    treatment,
    all.vars(confounders),
    tv_vars,
    id,
    time,
    censoring
  ))
  keep_cols <- intersect(keep_cols, names(data))

  data <- data[, keep_cols, with = FALSE]

  if (!is.null(id) && !is.null(time)) {
    data <- create_lag_vars(data, treatment, tv_vars, id, time, history)
    warn_confounder_variation(data, confounders, confounders_tv, id)
  }

  data
}

#' Create lagged columns for treatment and time-varying confounders
#'
#' @param data A data.table sorted by id and time.
#' @param treatment Character treatment column name(s).
#' @param tv_vars Character vector of time-varying confounder column names.
#' @param id Character ID column name.
#' @param time Character time column name.
#' @param history Positive integer Markov order or `Inf`.
#' @return The data.table with added lag columns (e.g. `lag1_A`, `lag2_L`).
#' @noRd
create_lag_vars <- function(data, treatment, tv_vars, id, time, history) {
  lag_vars <- c(treatment, tv_vars)
  max_lag <- if (is.infinite(history)) {
    time_vec <- data[[time]]
    id_vec <- data[[id]]
    max(tapply(time_vec, id_vec, function(x) max(x) - min(x)))
  } else {
    as.integer(history)
  }

  data <- data.table::copy(data)
  data.table::setkeyv(data, c(id, time))

  for (v in lag_vars) {
    for (k in seq_len(max_lag)) {
      lag_col <- paste0("lag", k, "_", v)
      data[,
        (lag_col) := data.table::shift(get(v), n = k, type = "lag"), # nolint: object_usage_linter
        by = c(id)
      ]
    }
  }

  data
}

#' Warn if confounders appear misclassified as baseline vs time-varying
#'
#' @param data A data.table of person-period data.
#' @param confounders One-sided formula of baseline confounders.
#' @param confounders_tv One-sided formula of time-varying confounders or `NULL`.
#' @param id Character ID column name.
#' @return `NULL` invisibly; issues warnings for suspected misclassifications.
#' @noRd
warn_confounder_variation <- function(data, confounders, confounders_tv, id) {
  baseline_vars <- all.vars(confounders)
  tv_vars <- if (!is.null(confounders_tv)) {
    all.vars(confounders_tv)
  } else {
    character(0)
  }

  for (v in tv_vars) {
    if (v %in% names(data)) {
      n_unique <- data[, data.table::uniqueN(get(v)), by = c(id)]$V1
      if (all(n_unique == 1L)) {
        rlang::warn(
          paste0(
            "Time-varying confounder '",
            v,
            "' does not vary within any individual. ",
            "Consider moving it to `confounders` (baseline)."
          )
        )
      }
    }
  }

  for (v in baseline_vars) {
    if (v %in% names(data)) {
      n_unique <- data[, data.table::uniqueN(get(v)), by = c(id)]$V1
      if (any(n_unique > 1L)) {
        rlang::warn(
          paste0(
            "Baseline confounder '",
            v,
            "' varies within some individuals. ",
            "Consider moving it to `confounders_tv` (time-varying)."
          )
        )
      }
    }
  }
}
