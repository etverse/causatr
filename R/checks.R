#' @noRd
check_string <- function(
  x,
  arg = rlang::caller_arg(x),
  call = rlang::caller_env()
) {
  if (!rlang::is_string(x)) {
    rlang::abort(
      paste0("`", arg, "` must be a single character string."),
      call = call
    )
  }
}

#' @noRd
check_col_exists <- function(
  data,
  col,
  arg = rlang::caller_arg(col),
  call = rlang::caller_env()
) {
  if (!col %in% names(data)) {
    rlang::abort(
      paste0("Column `", col, "` (", arg, ") not found in `data`."),
      call = call
    )
  }
}

#' @noRd
check_formula <- function(
  x,
  arg = rlang::caller_arg(x),
  call = rlang::caller_env()
) {
  if (!inherits(x, "formula")) {
    rlang::abort(
      paste0("`", arg, "` must be a formula (e.g. `~ L1 + L2`)."),
      call = call
    )
  }
}

#' @noRd
check_intervention_list <- function(x, call = rlang::caller_env()) {
  if (!is.list(x) || length(x) == 0) {
    rlang::abort(
      "`interventions` must be a named list with at least one intervention.",
      call = call
    )
  }
  if (is.null(names(x)) || any(names(x) == "")) {
    rlang::abort(
      "All elements of `interventions` must be named.",
      call = call
    )
  }
  for (nm in names(x)) {
    if (!inherits(x[[nm]], "causatr_intervention")) {
      rlang::abort(
        paste0(
          "`interventions$",
          nm,
          "` must be a `causatr_intervention` object. ",
          "Use `static()`, `shift()`, `dynamic()`, etc."
        ),
        call = call
      )
    }
  }
}

#' @noRd
check_causat_inputs <- function(
  data,
  outcome,
  treatment,
  confounders,
  method,
  id,
  time,
  call = rlang::caller_env()
) {
  check_string(outcome, call = call)
  check_string(treatment, call = call)
  check_formula(confounders, call = call)
  check_col_exists(data, outcome, call = call)
  check_col_exists(data, treatment, call = call)

  confounder_vars <- all.vars(confounders)
  missing_vars <- setdiff(confounder_vars, names(data))
  if (length(missing_vars) > 0) {
    rlang::abort(
      paste0(
        "Confounder variable(s) not found in `data`: ",
        paste(missing_vars, collapse = ", ")
      ),
      call = call
    )
  }

  if (outcome == treatment) {
    rlang::abort(
      "`outcome` and `treatment` must be different columns.",
      call = call
    )
  }

  if (!is.null(id)) {
    check_string(id, arg = "id", call = call)
    check_col_exists(data, id, arg = "id", call = call)
  }
  if (!is.null(time)) {
    check_string(time, arg = "time", call = call)
    check_col_exists(data, time, arg = "time", call = call)
  }
  if (xor(is.null(id), is.null(time))) {
    rlang::abort(
      "Both `id` and `time` must be provided together for longitudinal data.",
      call = call
    )
  }
}
