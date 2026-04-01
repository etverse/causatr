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
    el <- x[[nm]]
    if (is.null(el)) {
      next
    }
    if (is.list(el) && !inherits(el, "causatr_intervention")) {
      if (is.null(names(el)) || any(names(el) == "")) {
        rlang::abort(
          paste0(
            "`interventions$",
            nm,
            "` is a list but not all elements are named. ",
            "For multivariate treatment, supply a named list with one ",
            "entry per treatment variable."
          ),
          call = call
        )
      }
      for (sub_nm in names(el)) {
        if (!inherits(el[[sub_nm]], "causatr_intervention")) {
          rlang::abort(
            paste0(
              "`interventions$",
              nm,
              "$",
              sub_nm,
              "` must be a `causatr_intervention` object. ",
              "Use `static()`, `shift()`, `dynamic()`, etc."
            ),
            call = call
          )
        }
      }
    } else if (!inherits(el, "causatr_intervention")) {
      rlang::abort(
        paste0(
          "`interventions$",
          nm,
          "` must be a `causatr_intervention` object or `NULL` (natural ",
          "course). Use `static()`, `shift()`, `dynamic()`, etc."
        ),
        call = call
      )
    }
  }
}

#' @noRd
check_estimand_compat <- function(
  estimand,
  fit_method,
  fit_estimand,
  call = rlang::caller_env()
) {
  if (is.null(estimand)) {
    return(invisible(NULL))
  }

  if (fit_method %in% c("ipw", "matching") && estimand != fit_estimand) {
    rlang::abort(
      paste0(
        "For method = '",
        fit_method,
        "', the estimand is fixed at fitting time because it determines the ",
        "weights. Refit with causat(estimand = '",
        estimand,
        "')."
      ),
      call = call
    )
  }
}

#' @noRd
check_estimand_trt_compat <- function(
  estimand,
  treatment,
  type,
  data = NULL,
  call = rlang::caller_env()
) {
  if (estimand == "ATE") {
    return(invisible(NULL))
  }

  msg <- paste0(
    "estimand = '",
    estimand,
    "' is only defined for binary point treatments. ",
    "Use estimand = 'ATE' or subset = quote(...) for subgroup effects."
  )

  if (type == "longitudinal") {
    rlang::abort(msg, call = call)
  }

  if (length(treatment) > 1L) {
    rlang::abort(msg, call = call)
  }

  if (!is.null(data)) {
    trt_vals <- unique(stats::na.omit(data[[treatment]]))
    if (!all(trt_vals %in% c(0, 1))) {
      rlang::abort(msg, call = call)
    }
  }
}

#' @noRd
check_treatment_nas <- function(
  data,
  treatment,
  censoring,
  call = rlang::caller_env()
) {
  trt_cols <- treatment
  for (col in trt_cols) {
    n_na <- sum(is.na(data[[col]]))
    if (n_na > 0 && is.null(censoring)) {
      rlang::abort(
        c(
          paste0(
            "Treatment variable '",
            col,
            "' has ",
            n_na,
            " missing value",
            if (n_na == 1) "" else "s",
            "."
          ),
          i = paste0(
            "Use `censoring = '...'` for inverse probability of ",
            "censoring weights."
          ),
          i = paste0(
            "Use `causat_mice()` with a mice `mids` object for ",
            "multiple imputation."
          ),
          i = "Or remove incomplete cases before calling `causat()`."
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
  confounders_tv,
  method,
  estimand,
  id,
  time,
  history,
  call = rlang::caller_env()
) {
  type <- if (!is.null(id) && !is.null(time)) "longitudinal" else "point"

  if (type == "longitudinal" && method %in% c("ipw", "matching")) {
    rlang::abort(
      paste0(
        "method = '",
        method,
        "' does not support longitudinal data. Use method = 'gcomp'."
      ),
      call = call
    )
  }

  check_estimand_trt_compat(estimand, treatment, type, data = data, call = call)

  check_string(outcome, call = call)

  if (!is.character(treatment) || length(treatment) == 0L) {
    rlang::abort(
      "`treatment` must be a character string or character vector.",
      call = call
    )
  }

  check_formula(confounders, call = call)
  check_col_exists(data, outcome, call = call)

  for (trt in treatment) {
    check_col_exists(data, trt, arg = "treatment", call = call)
  }

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

  if (any(outcome == treatment)) {
    rlang::abort(
      "`outcome` and `treatment` must be different columns.",
      call = call
    )
  }

  if (!is.null(confounders_tv)) {
    check_formula(confounders_tv, arg = "confounders_tv", call = call)
    tv_vars <- all.vars(confounders_tv)
    missing_tv <- setdiff(tv_vars, names(data))
    if (length(missing_tv) > 0) {
      rlang::abort(
        paste0(
          "Time-varying confounder variable(s) not found in `data`: ",
          paste(missing_tv, collapse = ", ")
        ),
        call = call
      )
    }
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

  if (!is.null(history)) {
    if (
      !rlang::is_scalar_double(history) &&
        !rlang::is_scalar_integer(history) &&
        !identical(history, Inf)
    ) {
      rlang::abort(
        "`history` must be a positive integer or `Inf`.",
        call = call
      )
    }
    if (!is.infinite(history) && (history < 1 || history != floor(history))) {
      rlang::abort(
        "`history` must be a positive integer or `Inf`.",
        call = call
      )
    }
  }
}
