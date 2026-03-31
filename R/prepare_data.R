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

#' @noRd
create_lag_vars <- function(data, treatment, tv_vars, id, time, history) {
  lag_vars <- c(treatment, tv_vars)
  max_lag <- if (is.infinite(history)) {
    data[, max(get(time) - min(get(time))), by = id]$V1 |> max()
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
