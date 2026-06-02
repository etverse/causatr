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
#' @param confounders_outcome One-sided formula of outcome-model confounders, or
#'   `NULL`.
#' @param confounders_treatment One-sided formula of treatment-model confounders,
#'   or `NULL`.
#' @param confounders_censoring One-sided formula of censoring-model confounders,
#'   or `NULL`.
#' @param confounders_sampling One-sided formula of sampling-model confounders,
#'   or `NULL`.
#' @param confounders_tv_outcome One-sided formula of time-varying outcome-model
#'   confounders, or `NULL`.
#' @param confounders_tv_treatment One-sided formula of time-varying
#'   treatment-model confounders, or `NULL`.
#' @param id Character ID column name or `NULL`.
#' @param time Character time column name or `NULL`.
#' @param censoring Character censoring column name or `NULL`.
#' @param history Positive integer Markov order or `Inf`.
#' @param target Character sampling indicator column name or `NULL`.
#' @param stratified Character baseline column name to retain for
#'   stratified ICE, or `NULL`. Added to the keep set so the column
#'   survives column-stripping.
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
  confounders_outcome = NULL,
  confounders_treatment = NULL,
  confounders_censoring = NULL,
  confounders_sampling = NULL,
  confounders_tv_outcome = NULL,
  confounders_tv_treatment = NULL,
  id = NULL,
  time = NULL,
  censoring = NULL,
  history = 1L,
  cluster = NULL,
  target = NULL,
  stratified = NULL,
  call = rlang::caller_env()
) {
  # Coerce to data.table once, up front. Everything downstream
  # assumes data.table semantics: in-place `:=`, `.SD`, by-group
  # operations. Copying here avoids mutating the user's data frame.
  if (!data.table::is.data.table(data)) {
    data <- data.table::as.data.table(data)
  }

  # R12 (2026-04-15 review): reject user data carrying a column name
  # reserved by causatr internals (`.pseudo_y`). Otherwise the in-place
  # `:=` mutations below (or inside ice.R) would silently clobber the
  # user's column.
  check_reserved_cols(data)

  tv_vars <- if (!is.null(confounders_tv)) {
    all.vars(confounders_tv)
  } else {
    character(0)
  }

  # Collect variables from per-component confounder formulas. These may
  # reference columns not in the unified `confounders` / `confounders_tv`
  # formulas (the whole point of per-component formulas), so they must be
  # included in the retention set or the downstream model fit will fail.
  per_component_vars <- character(0)
  per_component_tv_vars <- character(0)
  for (fml in list(
    confounders_outcome,
    confounders_treatment,
    confounders_censoring,
    confounders_sampling
  )) {
    if (!is.null(fml)) {
      per_component_vars <- c(per_component_vars, all.vars(fml))
    }
  }
  for (fml in list(confounders_tv_outcome, confounders_tv_treatment)) {
    if (!is.null(fml)) {
      per_component_tv_vars <- c(per_component_tv_vars, all.vars(fml))
    }
  }
  # Time-varying per-component variables need lags too
  tv_vars <- unique(c(tv_vars, per_component_tv_vars))

  # Strip to only the columns we need. Keeping extra columns is
  # harmless for correctness but wastes memory across bootstrap
  # iterations where data is repeatedly resampled. `intersect`
  # handles the case where `confounders` references transformed
  # variables that aren't themselves column names -- those transforms
  # are re-evaluated at fit time.
  keep_cols <- unique(c(
    outcome,
    treatment,
    if (!is.null(confounders)) all.vars(confounders) else character(0),
    tv_vars,
    per_component_vars,
    id,
    time,
    censoring,
    cluster,
    target,
    stratified
  ))
  keep_cols <- intersect(keep_cols, names(data))

  data <- data[, keep_cols, with = FALSE]

  # Longitudinal-only: build `lag1_A`, `lag2_A`, ..., `lag1_L`, ...
  # columns up to `history`. Materializing lags here rather than
  # computing them on-the-fly at model-fit time means every fitting
  # function (ICE, longitudinal IPW, longitudinal AIPW) can express
  # its formula as a plain character string like `Y ~ A + lag1_A + L`
  # that works with any model_fn (glm, gam, ...) without requiring
  # model-fn-specific lag hooks. It also ensures that bootstrap
  # resamples see the same lag structure without re-running the shift.
  if (!is.null(id) && !is.null(time)) {
    data <- create_lag_vars(data, treatment, tv_vars, id, time, history)
    # Sanity check: warn about confounders that are obviously
    # misclassified (e.g. baseline that varies within person, or
    # time-varying that doesn't). These are often typos in the
    # user's formula rather than intentional. Pass `treatment` so
    # interaction terms that include the treatment (e.g.
    # `~ sex + A:sex`) don't trigger a spurious warning about `A`
    # varying within individuals -- `A` is the treatment, not a
    # baseline confounder.
    warn_confounder_variation(
      data,
      confounders,
      confounders_tv,
      id,
      treatment
    )
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
  # Lag every treatment column and every time-varying confounder.
  # Baseline confounders are time-invariant and don't need lagging.
  lag_vars <- c(treatment, tv_vars)

  # `history = Inf` means "full history": cap at the maximum number
  # of time points per individual minus one (there's no lag at t=0).
  # Using `tapply(id, id, length)` to count rows per id is faster
  # than `uniqueN` here because we only need max, not the full set.
  max_lag <- if (is.infinite(history)) {
    id_vec <- data[[id]]
    max(tapply(id_vec, id_vec, length)) - 1L
  } else {
    as.integer(history)
  }

  # Deep copy + key. `setkeyv` is load-bearing: `data.table::shift`
  # with `by = id` processes rows in their current order within each
  # group, so if an id's rows are not sorted by time the lag at
  # position k will correspond to the wrong time point. `setkeyv`
  # sorts by (id, time) in place, guaranteeing correct lag alignment.
  data <- data.table::copy(data)
  data.table::setkeyv(data, c(id, time))

  # Nested loop: (variable x lag order). Each iteration adds one
  # new column in place via `:=`. `by = c(id)` keeps the shift
  # per-individual so lags at t=0 within an id are NA, not bleeding
  # from the previous id's last row.
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
warn_confounder_variation <- function(
  data,
  confounders,
  confounders_tv,
  id,
  treatment = character(0)
) {
  # Sanity-check each confounder against the user's classification.
  # Two error modes:
  #   (a) a variable declared time-varying that doesn't actually vary
  #       within any individual (waste of a lag column; treat as baseline)
  #   (b) a variable declared baseline that DOES vary within some
  #       individuals (ICE will use only the baseline value, losing info)
  # We warn instead of abort because some designs legitimately have
  # variables that are constant in-sample but meaningful time-varying
  # in principle (e.g. a dummy for "first 6 months of follow-up" in a
  # dataset with exactly 6 months).
  #
  # Subtract the treatment columns from `baseline_vars`. Users
  # commonly write `confounders = ~ L + A:L` to encode an A x L
  # interaction on the outcome; `all.vars(~ L + A:L)` returns both
  # `L` and `A`, and without this subtraction `A` would trip the
  # within-individual variation check below. `A` is the treatment,
  # not a baseline confounder.
  baseline_vars <- setdiff(all.vars(confounders), treatment)
  tv_vars <- if (!is.null(confounders_tv)) {
    all.vars(confounders_tv)
  } else {
    character(0)
  }

  # Check tv_vars: for each, count distinct values per individual.
  # If every individual has exactly 1 unique value, the variable is
  # effectively baseline.
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
          ),
          class = "causatr_tv_confounder_constant"
        )
      }
    }
  }

  # Check baseline_vars: for each, warn if ANY individual has > 1
  # unique value. Unlike the tv check, a single violator is enough
  # to surface -- the baseline treatment would drop the within-person
  # variation entirely.
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
          ),
          class = "causatr_baseline_confounder_varying"
        )
      }
    }
  }
}
