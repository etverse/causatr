#' Check that a value is a single string
#'
#' @param x Value to check.
#' @param arg Argument name for error messages.
#' @param call Caller environment for error messages.
#' @return `NULL` invisibly; aborts if `x` is not a scalar character.
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

#' Check that a column exists in data
#'
#' @param data A data.frame or data.table.
#' @param col Character column name to look up.
#' @param arg Argument name for error messages.
#' @param call Caller environment for error messages.
#' @return `NULL` invisibly; aborts if `col` is not in `names(data)`.
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

#' Check that a value is a formula
#'
#' @param x Value to check.
#' @param arg Argument name for error messages.
#' @param call Caller environment for error messages.
#' @return `NULL` invisibly; aborts if `x` is not a formula.
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

#' Validate an interventions list
#'
#' @param x A named list of interventions (each `causatr_intervention`,
#'   `NULL`, or a named list of `causatr_intervention` for multivariate
#'   treatments).
#' @param call Caller environment for error messages.
#' @return `NULL` invisibly; aborts on invalid structure.
#' @noRd
check_intervention_list <- function(x, call = rlang::caller_env()) {
  # Top-level shape: must be a non-empty named list.
  if (!is.list(x) || length(x) == 0) {
    rlang::abort(
      "`interventions` must be a named list with at least one intervention.",
      call = call
    )
  }
  # Every element needs a name -- used as row labels in the results
  # table and as reference targets in the contrast step.
  if (is.null(names(x)) || any(names(x) == "")) {
    rlang::abort(
      "All elements of `interventions` must be named.",
      call = call
    )
  }
  # Names must be unique: duplicates otherwise slip through silently
  # and produce stale rows in the estimates / contrasts tables (the
  # second entry with the same name shadows the first in named-vector
  # indexing, but rbindlist / data.table still emits both).
  if (anyDuplicated(names(x))) {
    dups <- unique(names(x)[duplicated(names(x))])
    rlang::abort(
      paste0(
        "`interventions` has duplicated name(s): ",
        paste0("'", dups, "'", collapse = ", "),
        ". Each intervention must have a unique name."
      ),
      call = call
    )
  }
  # Per-element validation. Three valid shapes:
  #   (a) `NULL`           -- natural course (observed treatment as-is)
  #   (b) causatr_intervention -- bare intervention for scalar treatment
  #   (c) named list of causatr_intervention -- for multivariate treatment,
  #       one entry per treatment column (e.g. list(A1 = static(1), A2 = shift(-10)))
  for (nm in names(x)) {
    el <- x[[nm]]
    if (is.null(el)) {
      next
    }
    # Case (c) detection: plain list that isn't itself a
    # `causatr_intervention`. The class check is the discriminator --
    # a `causatr_intervention` is technically a list under the hood.
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

#' Check that an estimand override is compatible with the fit estimator
#'
#' @param estimand Character estimand requested in `contrast()`, or `NULL`.
#' @param fit_estimator Character causal estimator (`"gcomp"`, `"ipw"`, or
#'   `"matching"`).
#' @param fit_estimand Character estimand that was used at fitting time.
#' @param call Caller environment for error messages.
#' @return `NULL` invisibly; aborts if IPW/matching estimand is changed.
#' @noRd
check_estimand_compat <- function(
  estimand,
  fit_estimator,
  fit_estimand,
  call = rlang::caller_env()
) {
  # No override requested -- nothing to check.
  if (is.null(estimand)) {
    return(invisible(NULL))
  }

  # IPW and matching can't switch estimand at contrast time because
  # the weights / matched sets were estimated under the fit-time
  # estimand -- e.g. ATT weights upweight the control-over-treated
  # distribution, so re-averaging over "everyone" (ATE) under ATT
  # weights doesn't give you E[Y^a]. G-comp doesn't have this
  # problem: the outcome model is estimand-agnostic, and the estimand
  # only affects which rows we average predictions over.
  if (fit_estimator %in% c("ipw", "matching") && estimand != fit_estimand) {
    rlang::abort(
      paste0(
        "For estimator = '",
        fit_estimator,
        "', the estimand is fixed at fitting time because it determines the ",
        "weights. Refit with causat(estimand = '",
        estimand,
        "')."
      ),
      call = call
    )
  }
}

#' Check estimand-treatment compatibility
#'
#' ATT/ATC are only defined for binary point treatments.
#'
#' @param estimand Character estimand (`"ATE"`, `"ATT"`, or `"ATC"`).
#' @param treatment Character vector of treatment column name(s).
#' @param type `"point"` or `"longitudinal"`.
#' @param data Optional data.frame used to verify treatment is binary.
#' @param call Caller environment for error messages.
#' @return `NULL` invisibly; aborts on incompatible combinations.
#' @noRd
check_estimand_trt_compat <- function(
  estimand,
  treatment,
  type,
  data = NULL,
  call = rlang::caller_env()
) {
  # ATE is well-defined for any treatment (binary, continuous,
  # categorical, multivariate), so no check needed.
  if (estimand == "ATE") {
    return(invisible(NULL))
  }

  # ATT/ATC are defined as "average effect among the treated/controls",
  # which requires a natural binary 0/1 split of the population. For
  # continuous treatment there is no "treated group"; for longitudinal
  # there is no single baseline treatment to condition on. Multivariate
  # falls in the same bucket -- there's no unique "treated" level.
  # Factor or character-coded binary treatments are also rejected: the
  # downstream `contrast()` filter is `trt_sym == 1L`, so even a
  # two-level factor with levels `c("control", "treated")` silently
  # returns zero rows.
  msg <- paste0(
    "estimand = '",
    estimand,
    "' is only defined for binary point treatments coded as 0/1. ",
    "Use estimand = 'ATE' or subset = quote(...) for subgroup effects ",
    "(or recode the treatment as integer 0/1 if it already has two levels)."
  )

  if (type == "longitudinal") {
    rlang::abort(msg, call = call)
  }

  if (length(treatment) > 1L) {
    rlang::abort(msg, call = call)
  }

  # If we have data at this point, confirm the single treatment column
  # is actually coded 0/1 -- not just character values like "A"/"B",
  # which would fail silently later when contrast() tries to filter
  # on treatment == 1. The error message above already tells the user
  # to recode so they don't have to reverse-engineer the check.
  if (!is.null(data)) {
    trt_vals <- unique(stats::na.omit(data[[treatment]]))
    if (!all(trt_vals %in% c(0, 1))) {
      rlang::abort(msg, call = call)
    }
  }
}

#' Check compatibility between estimand and intervention types
#'
#' @description
#' Under the IPW / matching engines, the ATT and ATC estimands are
#' only well-defined when the intervention is static on a binary
#' treatment. The reasoning:
#'
#' - ATT = "the average treatment effect among those who actually
#'   received treatment". This requires (a) an unambiguous "treated"
#'   subpopulation (binary coding), and (b) a counterfactual
#'   treatment value that the subpopulation's observed treatment is
#'   being compared against -- i.e. a static target. For MTPs
#'   (`shift`, `scale_by`, `threshold`, `dynamic`, `ipsi`) the
#'   "among the treated" restriction is either undefined or exotic
#'   and the MTP literature does not use it.
#' - ATC is the symmetric statement for the untreated.
#'
#' Silently falling back to ATE weights under an ATT request would
#' return a pooled effect when the user asked for effect within a
#' subpopulation -- a silent estimand swap, exactly the kind of
#' mistake the package boundary checks exist to prevent. So we abort
#' at contrast time with `class = "causatr_bad_estimand_intervention"`
#' and point users to either `estimand = "ATE"` or
#' `estimator = "gcomp"` (which handles ATT/ATC under any
#' intervention natively via predict-then-average on the outcome
#' model).
#'
#' For `estimator = "gcomp"` the function is a no-op.
#'
#' Natural-course entries (`NULL`) are always allowed: the observed
#' marginal mean among the treated / controls is a well-defined
#' quantity for any estimator.
#'
#' @param estimand Character scalar: `"ATE"`, `"ATT"`, or `"ATC"`.
#'   Typically the *effective* estimand at contrast time -- i.e. the
#'   user's `estimand = ` override if present, otherwise
#'   `fit$estimand`.
#' @param interventions Named list of interventions from `contrast()`.
#'   Each element is a `causatr_intervention`, `NULL` (natural
#'   course), or a named list of `causatr_intervention` objects
#'   (multivariate treatment).
#' @param estimator Character scalar: `"gcomp"`, `"ipw"`, or
#'   `"matching"`.
#' @param call Caller environment for error messages.
#'
#' @return `NULL` invisibly; aborts with class
#'   `"causatr_bad_estimand_intervention"` on an invalid combination.
#'
#' @noRd
check_estimand_intervention_compat <- function(
  estimand,
  interventions,
  estimator,
  call = rlang::caller_env()
) {
  # G-comp handles every (estimand, intervention) combination
  # natively: the outcome model is estimand-agnostic, and the
  # estimand only affects which rows we average predictions over in
  # `compute_contrast()`. No gating needed.
  if (estimator == "gcomp") {
    return(invisible(NULL))
  }

  # ATE is well-defined under every intervention supported by the
  # IPW and matching engines.
  if (estimand == "ATE") {
    return(invisible(NULL))
  }

  # Collect non-static interventions (if any). Three shapes to
  # handle, mirroring `check_interventions_compat()`:
  #   - NULL                                -> natural course, skip
  #   - `causatr_intervention`              -> inspect `$type`
  #   - list of `causatr_intervention`      -> multivariate treatment;
  #                                           any non-static sub-intervention
  #                                           flags the whole regime
  bad <- character()
  for (nm in names(interventions)) {
    iv <- interventions[[nm]]
    if (is.null(iv)) {
      next
    }
    if (is.list(iv) && !inherits(iv, "causatr_intervention")) {
      for (sub_nm in names(iv)) {
        if (iv[[sub_nm]]$type != "static") {
          bad <- c(bad, paste0(nm, "$", sub_nm))
        }
      }
      next
    }
    if (iv$type != "static") {
      bad <- c(bad, nm)
    }
  }

  if (length(bad) == 0L) {
    return(invisible(NULL))
  }

  rlang::abort(
    c(
      paste0(
        "`estimand = '",
        estimand,
        "'` under `estimator = '",
        estimator,
        "'` only accepts static interventions."
      ),
      x = paste0(
        "Non-static intervention(s): ",
        paste0("'", bad, "'", collapse = ", "),
        "."
      ),
      i = "Use `estimand = 'ATE'` if you want the MTP / shift / IPSI effect on the full population.",
      i = "Use `estimator = 'gcomp'` if you need ATT/ATC under a non-static intervention -- gcomp handles this via predict-then-average on the outcome model, which works for any estimand x intervention combination."
    ),
    class = "causatr_bad_estimand_intervention",
    call = call
  )
}


#' Check for missing treatment values
#'
#' @param data A data.frame or data.table.
#' @param treatment Character vector of treatment column name(s).
#' @param censoring Character censoring column name, or `NULL`.
#' @param call Caller environment for error messages.
#' @return `NULL` invisibly; aborts if NAs found without censoring.
#' @noRd
check_treatment_nas <- function(
  data,
  treatment,
  censoring,
  call = rlang::caller_env()
) {
  # Missing treatment values must be handled explicitly -- silently
  # dropping them (as glm would by default via na.action) is wrong
  # because the dropped rows may be MAR, not MCAR, and biasing the
  # marginal mean. Two legitimate handling paths today:
  #   1. Censoring indicator -> IPCW via `censoring = "col"`
  #   2. Manual complete-case subset -> user removes rows before calling
  # (Multiple imputation via `causat_mice()` is not implemented; the
  # stub is unexported so it is not advertised in the abort hint.)
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
          i = "Or remove incomplete cases before calling `causat()`."
        ),
        call = call
      )
    }
  }
}

#' Validate all inputs to causat()
#'
#' @param data A data.frame or data.table.
#' @param outcome Character outcome column name.
#' @param treatment Character treatment column name(s).
#' @param confounders One-sided formula of baseline confounders.
#' @param confounders_tv One-sided formula of time-varying confounders, or
#'   `NULL`.
#' @param estimator Character causal estimator.
#' @param estimand Character estimand.
#' @param id Character ID column name, or `NULL`.
#' @param time Character time column name, or `NULL`.
#' @param history Positive integer or `Inf`.
#' @param call Caller environment for error messages.
#' @return `NULL` invisibly; aborts on any validation failure.
#' @noRd
check_causat_inputs <- function(
  data,
  outcome,
  treatment,
  confounders,
  confounders_tv,
  estimator,
  estimand,
  id,
  time,
  history,
  call = rlang::caller_env()
) {
  # `type` is inferred from the presence of id/time (both present =>
  # longitudinal, both absent => point). The cross-validation below
  # catches the "only one present" mistake.
  type <- if (!is.null(id) && !is.null(time)) "longitudinal" else "point"

  # Estimator x type compatibility: longitudinal data is ICE-only.
  # IPW/matching for longitudinal would need weightitMSM/matchit for
  # repeated measures, which is not supported.
  if (type == "longitudinal" && estimator %in% c("ipw", "matching")) {
    rlang::abort(
      paste0(
        "estimator = '",
        estimator,
        "' does not support longitudinal data. Use estimator = 'gcomp'."
      ),
      call = call
    )
  }

  # Multivariate treatments for IPW/matching are not supported; the
  # propensity / matched-design story is more involved than simple
  # per-component application.
  if (length(treatment) > 1L && estimator %in% c("ipw", "matching")) {
    rlang::abort(
      paste0(
        "Multivariate treatments are not supported for estimator = '",
        estimator,
        "'. Use estimator = 'gcomp' for joint interventions on multiple ",
        "treatments, or fit separate models for each treatment."
      ),
      call = call
    )
  }

  check_estimand_trt_compat(estimand, treatment, type, data = data, call = call)

  # Per-argument validation: each helper aborts on its own error
  # message, so we just call them in sequence.
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

  # Confounder columns must all exist. `all.vars()` extracts plain
  # variable names from the formula's LHS/RHS, stripping transforms
  # (`I(age^2)` -> `age`), which is what we want for column-name
  # checking.
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

  # Outcome and treatment must be distinct. Catches a common typo
  # where users accidentally pass the outcome column as the treatment.
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
  # id and time must come as a pair. `xor` catches the "only one
  # present" case; both NULL is fine (point treatment) and both
  # non-NULL is fine (longitudinal).
  if (xor(is.null(id), is.null(time))) {
    rlang::abort(
      "Both `id` and `time` must be provided together for longitudinal data.",
      call = call
    )
  }

  # `history` must be a positive integer or Inf. The two-step check
  # is ugly but necessary: `is_scalar_integer` rejects `5` (which R
  # treats as double literal), so we also accept scalar double that
  # equals its floor. `identical(., Inf)` is the clean way to test
  # for the special-case "all history".
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


#' Validate an external weight vector before handing it to a fitter
#'
#' @description
#' Up-front check on the `weights` argument passed to `causat()`.
#' Without this guard, non-finite or negative
#' weights silently fall through to the downstream GLM, which either
#' aborts with a cryptic message or (worse) produces NaN estimates.
#' Reject at the causatr boundary so users see a specific error with
#' the failing call site.
#'
#' Zero weights are allowed as a pass-through even though they carry
#' no information on their own -- users sometimes use zero weights to
#' implement a conditional subset, and the downstream fitter handles
#' it fine.
#'
#' @param weights Numeric vector or `NULL`.
#' @param n Expected length (number of rows in the data passed to
#'   `causat()`).
#' @param call Caller environment for error messages.
#'
#' @return `NULL` invisibly; aborts on invalid input.
#'
#' @noRd
check_weights <- function(weights, n, call = rlang::caller_env()) {
  if (is.null(weights)) {
    return(invisible(NULL))
  }
  if (!is.numeric(weights)) {
    rlang::abort(
      "`weights` must be numeric.",
      call = call
    )
  }
  if (length(weights) != n) {
    rlang::abort(
      paste0(
        "`weights` must have length equal to `nrow(data)` (",
        n,
        "), got ",
        length(weights),
        "."
      ),
      call = call
    )
  }
  if (anyNA(weights)) {
    rlang::abort(
      paste0(
        "`weights` contains ",
        sum(is.na(weights)),
        " missing value(s). ",
        "Drop those rows or impute before calling `causat()`."
      ),
      call = call
    )
  }
  if (any(!is.finite(weights))) {
    rlang::abort(
      "`weights` contains non-finite value(s) (Inf / NaN).",
      call = call
    )
  }
  if (any(weights < 0)) {
    rlang::abort(
      "`weights` must be non-negative.",
      call = call
    )
  }
  invisible(NULL)
}


#' Reject `na.action = na.exclude` forwarded through `...`
#'
#' @description
#' `stats::glm(..., na.action = na.exclude)` pads `residuals(model, "working")`
#' with `NA`s so its length equals the original data, *not* the post-omit
#' `model.matrix(model)` row count. The variance engine's
#' `prepare_model_if()` computes `r_score = residuals * working_weights`
#' and then relies on `length(r_score) == nrow(X_fit)` -- under `na.exclude`
#' this invariant breaks, recycling silently corrupts the IF, and R only
#' emits a "longer object length is not a multiple of shorter object length"
#' warning. Downstream SEs are meaningless.
#'
#' Rather than harden every place we touch residuals, we refuse `na.exclude`
#' at the `causat()` boundary. `na.omit` (the default) and `na.fail` are
#' the only sensible choices for a pipeline that builds its own
#' row-alignment bookkeeping from `fit_rows` / `model$na.action`.
#'
#' Verified via `/tmp/causatr_repro_issue7.R` on 2026-04-15: under
#' `na.action = na.exclude`, `prepare_model_if()` triggered the recycling
#' warning and returned a silently-wrong correction vector.
#'
#' @param ... Dots forwarded to `causat()`.
#' @param call Caller environment for error messages.
#'
#' @return `NULL` invisibly; aborts on `na.action = na.exclude`.
#'
#' @noRd
check_dots_na_action <- function(..., call = rlang::caller_env()) {
  dots <- list(...)
  if (!"na.action" %in% names(dots)) {
    return(invisible(NULL))
  }
  na_action <- dots$na.action
  # Accept `na.omit` (default) and `na.fail` (hard stop on NA).
  # Reject anything else -- notably `na.exclude`, the only other
  # base-R na.action, and any user-supplied function we can't reason
  # about. Match by identity to both the function and the name so
  # `na.action = na.exclude` and `na.action = "na.exclude"` both fail.
  ok <- FALSE
  if (is.function(na_action)) {
    ok <- identical(na_action, stats::na.omit) ||
      identical(na_action, stats::na.fail)
  } else if (is.character(na_action) && length(na_action) == 1L) {
    ok <- na_action %in% c("na.omit", "na.fail")
  }
  if (!ok) {
    rlang::abort(
      c(
        "`na.action` must be `na.omit` (default) or `na.fail`.",
        i = paste0(
          "causatr builds its own row-alignment bookkeeping from ",
          "`fit_rows` and the fitted model's `na.action` attribute. ",
          "`na.exclude` pads working residuals with NA and silently ",
          "corrupts the sandwich variance."
        ),
        i = "Drop NA rows before calling `causat()` or use `na.action = na.omit`."
      ),
      class = "causatr_bad_na_action",
      call = call
    )
  }
  invisible(NULL)
}
