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
  if (
    fit_estimator %in% c("ipw", "aipw", "matching") && estimand != fit_estimand
  ) {
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

  # ATT = E[Y^1 - Y^0 | A=1] and ATC = E[Y^1 - Y^0 | A=0] condition on
  # membership in a "treated" or "control" group, which presupposes a
  # binary 0/1 treatment. For continuous treatments there is no discrete
  # "treated group" to condition on; for longitudinal there is no single
  # baseline exposure defining membership; for multivariate there is no
  # unique "treated" level across the joint distribution.
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
      i = paste0(
        "Use `estimator = 'gcomp'` if you need ATT/ATC under a ",
        "non-static intervention -- gcomp handles this via ",
        "predict-then-average on the outcome model, which works ",
        "for any estimand x intervention combination."
      )
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
  target = NULL,
  call = rlang::caller_env()
) {
  # Missing treatment values must be handled explicitly -- silently
  # dropping them (as glm would by default via na.action) is wrong
  # because the dropped rows may be MAR, not MCAR, and biasing the
  # marginal mean. Three legitimate handling paths:
  #   1. Censoring indicator -> IPCW via `censoring = "col"`
  #   2. Multiple imputation -> impute with `mice::mice()`, pool with
  #      `causat_mice()` (imputes A from P(A | L) under MAR)
  #   3. Manual complete-case subset -> user removes rows before calling
  #
  # When transportability is active (`target` non-NULL), target rows
  # (S=0) have NA treatment by design -- they are not in the study and
  # have no observed A. Restrict the NA check to study rows (S=1).
  check_rows <- if (!is.null(target)) {
    data[[target]] == 1L
  } else {
    rep(TRUE, nrow(data))
  }

  trt_cols <- treatment
  for (col in trt_cols) {
    n_na <- sum(is.na(data[[col]][check_rows]))
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
            "Or impute upstream with `mice::mice()` and pool with ",
            "`causat_mice()` (multiple imputation for MAR treatment)."
          ),
          i = "Or remove incomplete cases before calling `causat()`."
        ),
        call = call
      )
    }
  }
}

#' Validate the `stratified` argument for stratified ICE
#'
#' @description
#' Stratified ICE fits a separate per-step outcome model for each level of
#' a baseline (time-invariant) grouping variable. This validator enforces
#' the three preconditions the feature requires:
#'
#' 1. The feature is only defined for longitudinal g-computation (ICE):
#'    point treatments and the IPW / AIPW / matching / SNM estimators have
#'    no per-step backward recursion to stratify. Rejected with class
#'    `causatr_stratified_not_ice`.
#' 2. The stratifying column must be **baseline** -- constant within every
#'    individual across time. A time-varying column would make stratum
#'    membership change across backward steps, which is not the supported
#'    estimand. Rejected with class `causatr_stratified_not_baseline`.
#' 3. The column must define a small number of discrete strata (factor,
#'    character, logical, or low-cardinality integer). A continuous column
#'    would yield one stratum per individual and leave every per-step
#'    model unidentified. Rejected with class `causatr_stratified_too_many`.
#'
#' Runs after `prepare_data()`, so `data` is a data.table and the
#' stratifying column has been retained in the keep set.
#'
#' @param data Prepared data.table.
#' @param stratified Character scalar column name, or `NULL` (no-op).
#' @param estimator Character estimator from `causat()`.
#' @param type `"point"` or `"longitudinal"`.
#' @param id Character ID column name, or `NULL`.
#' @param call Caller environment for error messages.
#'
#' @return `NULL` invisibly; aborts on any violation.
#'
#' @noRd
check_stratified <- function(
  data,
  stratified,
  estimator,
  type,
  id,
  call = rlang::caller_env()
) {
  if (is.null(stratified)) {
    return(invisible(NULL))
  }
  check_string(stratified, arg = "stratified", call = call)

  # (1) ICE-only. Reject before touching the column so the user gets the
  # most fundamental misuse first (e.g. `stratified` under IPW).
  if (!(estimator == "gcomp" && type == "longitudinal")) {
    rlang::abort(
      c(
        "`stratified` is only supported for longitudinal g-computation (ICE).",
        i = paste0(
          "Use `estimator = \"gcomp\"` with `id` and `time` for ",
          "stratified ICE."
        ),
        x = paste0(
          "Got estimator = \"",
          estimator,
          "\", type = \"",
          type,
          "\"."
        )
      ),
      class = "causatr_stratified_not_ice",
      call = call
    )
  }

  if (!stratified %in% names(data)) {
    rlang::abort(
      paste0("`stratified` column '", stratified, "' not found in `data`."),
      class = "causatr_stratified_not_found",
      call = call
    )
  }

  vals <- data[[stratified]]
  if (anyNA(vals)) {
    rlang::abort(
      c(
        paste0("`stratified` column '", stratified, "' contains NA values."),
        i = "The stratifying variable must be a complete baseline covariate."
      ),
      class = "causatr_stratified_na",
      call = call
    )
  }

  # (2) Time-invariance: every individual must carry a single stratum value.
  nu <- data[, data.table::uniqueN(get(stratified)), by = c(id)]$V1
  if (any(nu > 1L)) {
    rlang::abort(
      c(
        paste0(
          "`stratified` column '",
          stratified,
          "' varies within individuals."
        ),
        i = paste0(
          "Stratified ICE requires a baseline (time-invariant) grouping ",
          "variable."
        ),
        i = "Move time-varying structure into `confounders_tv` instead."
      ),
      class = "causatr_stratified_not_baseline",
      call = call
    )
  }

  # (3) Discrete, low-cardinality strata.
  n_levels <- data.table::uniqueN(vals)
  is_continuous <- is.double(vals) && any(vals != floor(vals))
  if (is_continuous || n_levels > 10L) {
    rlang::abort(
      c(
        paste0(
          "`stratified` column '",
          stratified,
          "' has ",
          n_levels,
          " distinct value(s)",
          if (is_continuous) " (continuous)" else "",
          "."
        ),
        i = "Stratified ICE needs a small number of discrete strata (<= 10).",
        i = "Discretise the variable (e.g. `cut()` into bins) before stratifying."
      ),
      class = "causatr_stratified_too_many",
      call = call
    )
  }

  invisible(NULL)
}

#' Validate the flexible-treatment ICE term (`treatment_form`)
#'
#' @description
#' `treatment_form` lets the treatment enter the per-step ICE outcome models
#' via a transformed term (`~ factor(A)`, `~ splines::ns(A, 3)`, ...) instead
#' of a bare numeric main effect, so a nonlinear or kinked counterfactual
#' response is not forced through a single linear slope. The intervention
#' still sets the numeric treatment column; only the model's design term
#' changes. This gate, called from [causat()] before any model is fitted,
#' enforces:
#'
#' 1. `treatment_form` is a one-sided formula. A two-sided formula or a
#'    non-formula is rejected with class `causatr_treatment_form_bad`.
#' 2. The flexible term is supported only for **longitudinal g-computation**
#'    (`estimator = "gcomp"` with `id` / `time`) -- point g-computation and
#'    the other estimators have no per-step ICE design matrix to retarget.
#'    Rejected with class `causatr_treatment_form_not_ice`.
#' 3. Every variable in the formula is a treatment column, and every treatment
#'    column appears. A covariate transform belongs in `confounders` /
#'    `confounders_tv`; an unreferenced treatment column would silently drop
#'    the exposure from the model. Either is rejected with class
#'    `causatr_treatment_form_bad`.
#'
#' @param treatment_form A one-sided formula, or `NULL` (the default, which is
#'   a no-op: the treatment enters as a bare numeric main effect).
#' @param treatment Character scalar or vector. Treatment column name(s).
#' @param estimator Character. The resolved estimator.
#' @param type Character. `"point"` or `"longitudinal"`.
#' @param call Caller environment for error reporting.
#'
#' @return Invisibly `NULL`. Called for its validation side effects.
#'
#' @seealso [check_stratified()], [ice_build_formula()]
#' @noRd
check_treatment_form <- function(
  treatment_form,
  treatment,
  estimator,
  type,
  call = rlang::caller_env()
) {
  if (is.null(treatment_form)) {
    return(invisible(NULL))
  }

  # (1) Must be a one-sided formula (`~ rhs`, length 2). A two-sided formula
  # (length 3) or a non-formula is a structural misuse -- reject before
  # inspecting terms so the message is about the shape, not the contents.
  if (!rlang::is_formula(treatment_form) || length(treatment_form) != 2L) {
    rlang::abort(
      c(
        "`treatment_form` must be a one-sided formula.",
        i = "For example `~ factor(A)` or `~ splines::ns(A, df = 3)`."
      ),
      class = "causatr_treatment_form_bad",
      call = call
    )
  }

  # (2) ICE-only. The flexible term retargets the per-step ICE design matrix;
  # point g-computation and the other estimators have no such per-step model.
  if (!(estimator == "gcomp" && type == "longitudinal")) {
    rlang::abort(
      c(
        paste0(
          "`treatment_form` is only supported for longitudinal ",
          "g-computation (ICE)."
        ),
        i = "Use `estimator = \"gcomp\"` with `id` and `time`.",
        x = paste0(
          "Got estimator = \"",
          estimator,
          "\", type = \"",
          type,
          "\"."
        )
      ),
      class = "causatr_treatment_form_not_ice",
      call = call
    )
  }

  # (3) Every term must reference a treatment column, and every treatment
  # column must appear. A term over a non-treatment variable belongs in
  # `confounders` / `confounders_tv`; an unreferenced treatment column would
  # silently drop the exposure from the per-step model.
  form_vars <- all.vars(treatment_form)
  extra <- setdiff(form_vars, treatment)
  if (length(extra) > 0L) {
    rlang::abort(
      c(
        "`treatment_form` may only reference the treatment column(s).",
        x = paste0(
          "Non-treatment variable(s): ",
          paste0("'", extra, "'", collapse = ", "),
          "."
        ),
        i = "Put covariate transforms in `confounders` / `confounders_tv`."
      ),
      class = "causatr_treatment_form_bad",
      call = call
    )
  }
  missing_trt <- setdiff(treatment, form_vars)
  if (length(missing_trt) > 0L) {
    rlang::abort(
      c(
        "Every treatment column must appear in `treatment_form`.",
        x = paste0(
          "Unreferenced treatment column(s): ",
          paste0("'", missing_trt, "'", collapse = ", "),
          "."
        )
      ),
      class = "causatr_treatment_form_bad",
      call = call
    )
  }

  invisible(NULL)
}

#' Validate natural-history MTP interventions against a fit
#'
#' @description
#' Natural-history modified treatment policies ([grace_period()] /
#' [carry_forward()]) are estimated by the augmented-data sequential regression
#' ([glmtp_iterate()]), which only exists for **longitudinal g-computation** on
#' a **discrete** treatment. This gate, called from [contrast()] when any
#' intervention is a `causatr_glmtp`, enforces:
#'
#' 1. `estimator = "gcomp"` and a longitudinal fit (`id` / `time`). IPW / AIPW /
#'    matching and point treatments have no augmented backward recursion.
#'    Rejected with class `causatr_glmtp_not_ice`.
#' 2. Transport (`target =`) is not supported for the augmented engine.
#'    Rejected with class `causatr_glmtp_not_ice`.
#' 3. A glmtp intervention is not mixed with a standard non-`NULL` intervention
#'    in the same contrast (a natural-course `NULL` reference is allowed and is
#'    routed through the identity policy). Rejected with class
#'    `causatr_glmtp_mixed`.
#' 4. The treatment is discrete with a tractable history enumeration -- delegated
#'    to [glmtp_support()] (class `causatr_glmtp_continuous_trt`) and
#'    [glmtp_check_tractable()] (class `causatr_glmtp_too_many`).
#'
#' A no-op when no intervention is a `causatr_glmtp`.
#'
#' @param fit A `causatr_fit`.
#' @param interventions Named list of interventions.
#' @param call Caller environment for error messages.
#'
#' @return `NULL` invisibly; aborts on any violation.
#'
#' @noRd
check_glmtp_compat <- function(fit, interventions, call = rlang::caller_env()) {
  is_g <- vapply(
    interventions,
    function(iv) inherits(iv, "causatr_glmtp"),
    logical(1)
  )
  if (!any(is_g)) {
    return(invisible(NULL))
  }

  # (1) ICE-only: longitudinal g-computation.
  if (!(fit$estimator == "gcomp" && fit$type == "longitudinal")) {
    rlang::abort(
      c(
        paste0(
          "Natural-history MTPs (`grace_period()` / `carry_forward()`) are ",
          "only supported for longitudinal g-computation."
        ),
        i = paste0(
          "Fit with `estimator = \"gcomp\"` and `id` / `time` columns ",
          "(longitudinal ICE)."
        ),
        x = paste0(
          "Got estimator = \"",
          fit$estimator,
          "\", type = \"",
          fit$type,
          "\"."
        )
      ),
      class = "causatr_glmtp_not_ice",
      call = call
    )
  }

  # (1b) Stratified ICE and the augmented engine are independent nuisance-model
  # choices that have not been composed; the augmented recursion would silently
  # ignore the stratifier. Reject rather than mislead.
  if (!is.null(fit$details$stratified)) {
    rlang::abort(
      c(
        paste0(
          "Natural-history MTPs are not supported together with stratified ",
          "ICE (`stratified =`)."
        ),
        i = "Refit without `stratified` to use `grace_period()` / `carry_forward()`."
      ),
      class = "causatr_glmtp_not_ice",
      call = call
    )
  }

  # (2) Transport is owned by a separate (pending) longitudinal-transport path.
  if (!is.null(fit$target)) {
    rlang::abort(
      c(
        "Natural-history MTPs do not support transportability (`target =`).",
        i = paste0(
          "The augmented sequential regression has no sampling-model channel; ",
          "drop `target` to estimate the study-population effect."
        )
      ),
      class = "causatr_glmtp_not_ice",
      call = call
    )
  }

  # (3) No mixing with standard non-NULL interventions. A NULL natural-course
  # reference is allowed (routed through the identity policy in glmtp_iterate()).
  is_std <- vapply(
    interventions,
    function(iv) !is.null(iv) && !inherits(iv, "causatr_glmtp"),
    logical(1)
  )
  if (any(is_std)) {
    rlang::abort(
      c(
        paste0(
          "Natural-history MTPs cannot be mixed with standard interventions ",
          "in one `contrast()` call."
        ),
        i = paste0(
          "Offending intervention(s): ",
          paste(shQuote(names(interventions)[is_std]), collapse = ", "),
          "."
        ),
        i = paste0(
          "Use only `grace_period()` / `carry_forward()` (a `NULL` ",
          "natural-course reference is allowed), or run a separate contrast."
        )
      ),
      class = "causatr_glmtp_mixed",
      call = call
    )
  }

  # (4) Discrete support + tractable enumeration. `glmtp_support()` rejects
  # continuous / factor / multivariate treatment; `glmtp_check_tractable()`
  # caps the worst-step blow-up. Use the budget from the first glmtp arm.
  budget <- interventions[[which(is_g)[1L]]]$budget %||% 1024L
  support <- glmtp_support(fit$data, fit$treatment, fit$censoring, call = call)
  glmtp_check_tractable(support, fit$details$n_times, budget, call = call)

  invisible(NULL)
}

#' Validate all inputs to causat()
#'
#' @param data A data.frame or data.table.
#' @param outcome Character outcome column name.
#' @param treatment Character treatment column name(s).
#' @param confounders One-sided formula of baseline confounders.
#' @param confounders_tv One-sided formula of time-varying confounders, or
#'   `NULL`.
#' @param confounders_outcome One-sided formula of outcome-model confounders,
#'   or `NULL`.
#' @param confounders_treatment One-sided formula of treatment-model
#'   confounders, or `NULL`.
#' @param confounders_censoring One-sided formula of censoring-model
#'   confounders, or `NULL`.
#' @param confounders_sampling One-sided formula of sampling-model confounders,
#'   or `NULL`.
#' @param confounders_tv_outcome One-sided formula of time-varying outcome-model
#'   confounders, or `NULL`.
#' @param confounders_tv_treatment One-sided formula of time-varying
#'   treatment-model confounders, or `NULL`.
#' @param estimator Character causal estimator.
#' @param estimand Character estimand.
#' @param id Character ID column name, or `NULL`.
#' @param time Character time column name, or `NULL`.
#' @param history Non-negative integer or `Inf`.
#' @param call Caller environment for error messages.
#' @return `NULL` invisibly; aborts on any validation failure.
#' @noRd
check_causat_inputs <- function(
  data,
  outcome,
  treatment,
  confounders,
  confounders_tv,
  confounders_outcome = NULL,
  confounders_treatment = NULL,
  confounders_censoring = NULL,
  confounders_sampling = NULL,
  confounders_tv_outcome = NULL,
  confounders_tv_treatment = NULL,
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

  # MatchIt operates on a single cross-section: it finds matched units
  # for each treated observation but has no notion of sequential treatment
  # periods or time-varying confounding. There is no `matchit` analogue
  # for the g-formula or product-of-ratios longitudinal weight chain.
  if (type == "longitudinal" && estimator == "matching") {
    rlang::abort(
      paste0(
        "estimator = '",
        estimator,
        "' does not support longitudinal data. Use estimator = 'gcomp' or 'ipw'."
      ),
      call = call
    )
  }

  # MatchIt expects a single binary or categorical treatment column with
  # no joint-treatment matching algorithm. IPW and AIPW handle
  # multivariate via sequential factorisation (`fit_treatment_models()`
  # + product density-ratio weights); per-component family / intervention
  # compatibility is checked downstream by the multivariate weight engine.
  if (length(treatment) > 1L && estimator == "matching") {
    rlang::abort(
      paste0(
        "Multivariate treatments are not supported for estimator = '",
        estimator,
        "'. Use estimator = 'gcomp' or estimator = 'ipw' for joint ",
        "interventions on multiple treatments."
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

  # `confounders` is NULL when the caller uses only per-component

  # formulas. Skip formula + column-existence checks when NULL; the
  # post-resolution estimator-requirement checks below will catch any
  # missing confounders.
  if (!is.null(confounders)) {
    check_formula(confounders, call = call)
  }
  check_col_exists(data, outcome, call = call)

  for (trt in treatment) {
    check_col_exists(data, trt, arg = "treatment", call = call)
  }

  # Confounder columns must all exist. `all.vars()` extracts plain
  # variable names from the formula's LHS/RHS, stripping transforms
  # (`I(age^2)` -> `age`), which is what we want for column-name
  # checking.
  if (!is.null(confounders)) {
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

  # Per-component confounder formulas: each is optional (NULL = fall back
  # to the unified `confounders` / `confounders_tv`). When supplied, must
  # be a formula whose variables all exist in `data`.
  per_component <- list(
    confounders_outcome = confounders_outcome,
    confounders_treatment = confounders_treatment,
    confounders_censoring = confounders_censoring,
    confounders_sampling = confounders_sampling,
    confounders_tv_outcome = confounders_tv_outcome,
    confounders_tv_treatment = confounders_tv_treatment
  )
  for (nm in names(per_component)) {
    val <- per_component[[nm]]
    if (!is.null(val)) {
      check_formula(val, arg = nm, call = call)
      comp_vars <- all.vars(val)
      missing_comp <- setdiff(comp_vars, names(data))
      if (length(missing_comp) > 0) {
        rlang::abort(
          paste0(
            nm,
            " variable(s) not found in `data`: ",
            paste(missing_comp, collapse = ", ")
          ),
          call = call
        )
      }
    }
  }

  # Post-resolution: verify that the estimator's required confounders are
  # available. The resolution logic mirrors causat.R: per-component
  # formula wins, then unified formula, then NULL. An estimator that
  # requires a component must have a non-NULL formula after resolution.
  # For longitudinal data (id + time supplied), time-varying
  # confounders alone are a valid source: per-period model formulas
  # are built from baseline_terms + tv_vars in the fit_* functions.
  conf_outcome <- confounders_outcome %||% confounders
  conf_treatment <- confounders_treatment %||% confounders
  is_longitudinal <- !is.null(id) && !is.null(time)

  has_outcome_conf <- !is.null(conf_outcome)
  if (!has_outcome_conf && is_longitudinal) {
    has_outcome_conf <- !is.null(confounders_tv_outcome %||% confounders_tv)
  }

  has_treatment_conf <- !is.null(conf_treatment)
  if (!has_treatment_conf && is_longitudinal) {
    has_treatment_conf <- !is.null(confounders_tv_treatment %||% confounders_tv)
  }

  needs_outcome <- estimator %in% c("gcomp", "aipw", "matching")
  needs_treatment <- estimator %in% c("ipw", "aipw", "matching", "snm")

  if (needs_outcome && !has_outcome_conf) {
    msg_suffix <- if (is_longitudinal) {
      "Supply `confounders_outcome`, `confounders`, or `confounders_tv`."
    } else {
      "Supply `confounders_outcome` or `confounders`."
    }
    rlang::abort(
      paste0(
        "`estimator = \"",
        estimator,
        "\"` requires outcome-model confounders. ",
        msg_suffix
      ),
      call = call
    )
  }

  if (needs_treatment && !has_treatment_conf) {
    msg_suffix <- if (is_longitudinal) {
      "Supply `confounders_treatment`, `confounders`, or `confounders_tv`."
    } else {
      "Supply `confounders_treatment` or `confounders`."
    }
    rlang::abort(
      paste0(
        "`estimator = \"",
        estimator,
        "\"` requires treatment-model confounders. ",
        msg_suffix
      ),
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
  # id and time must come as a pair. `xor` catches the "only one
  # present" case; both NULL is fine (point treatment) and both
  # non-NULL is fine (longitudinal).
  if (xor(is.null(id), is.null(time))) {
    rlang::abort(
      "Both `id` and `time` must be provided together for longitudinal data.",
      call = call
    )
  }

  # `history` controls the Markov lag order of the longitudinal treatment
  # model: history = 0 includes only current-period TV covariates (no
  # lags), history = k includes up to lag-k of TV covariates, and Inf
  # uses the full observed history. Must be a non-negative integer or Inf.
  # The two-step check is necessary: `is_scalar_integer` rejects `5` (which
  # R treats as a double literal), so we also accept a scalar double whose
  # value equals its floor. `identical(., Inf)` is the clean way to test
  # for the special-case "full history".
  if (!is.null(history)) {
    if (
      !rlang::is_scalar_double(history) &&
        !rlang::is_scalar_integer(history) &&
        !identical(history, Inf)
    ) {
      rlang::abort(
        "`history` must be a non-negative integer or `Inf`.",
        call = call
      )
    }
    if (!is.infinite(history) && (history < 0 || history != floor(history))) {
      rlang::abort(
        "`history` must be a non-negative integer or `Inf`.",
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


#' Validate `ipcw` and `censoring` parameter consistency
#'
#' @description
#' When `ipcw = TRUE`, the user wants causatr to fit an internal
#' censoring model. This requires a `censoring` column to be specified
#' and that column to be binary (0/1 with optional NA). Warns when the
#' censoring proportion is very high (>80%), signalling that IPCW
#' weights will be extreme.
#'
#' @param ipcw Logical. Whether built-in IPCW is requested.
#' @param censoring Character or NULL. Name of the censoring column.
#' @param censoring_col Optional. The actual censoring column values
#'   (for deeper validation). If NULL, only structural checks are run.
#' @param call Caller environment for error messages.
#'
#' @return `NULL` invisibly; aborts on invalid combinations.
#'
#' @noRd
check_ipcw_inputs <- function(
  ipcw,
  censoring,
  censoring_col = NULL,
  call = rlang::caller_env()
) {
  if (!isTRUE(ipcw)) {
    return(invisible(NULL))
  }

  if (is.null(censoring)) {
    rlang::abort(
      c(
        "`ipcw = TRUE` requires `censoring` to be specified.",
        i = paste0(
          "Provide the name of a censoring indicator ",
          "column (1 = censored, 0 = uncensored)."
        )
      ),
      class = "causatr_ipcw_no_censoring",
      call = call
    )
  }

  if (!is.null(censoring_col)) {
    vals <- stats::na.omit(as.integer(censoring_col))
    unique_vals <- sort(unique(vals))
    if (!all(unique_vals %in% c(0L, 1L))) {
      rlang::abort(
        c(
          paste0(
            "`",
            censoring,
            "` must be binary (0 = uncensored, 1 = censored). ",
            "Found values: ",
            paste(unique_vals, collapse = ", "),
            "."
          ),
          i = "Recode the censoring column to 0/1 before calling `causat()`."
        ),
        class = "causatr_ipcw_non_binary",
        call = call
      )
    }
    if (length(vals) > 0L) {
      p_cens <- mean(vals == 1L)
      if (p_cens > 0.80) {
        rlang::warn(
          c(
            paste0(
              round(100 * p_cens, 1),
              "% of observations are censored. ",
              "IPCW weights will be extreme and estimates may be unstable."
            ),
            i = "Consider whether the censoring model is correctly specified."
          ),
          class = "causatr_ipcw_high_censoring"
        )
      }
    }
  }

  invisible(NULL)
}


#' Validate transportability inputs
#'
#' @description
#' Checks that the `target` sampling indicator column is suitable for
#' fitting a sampling model \eqn{P(S = 1 \mid L)}.
#'
#' @param target Character. Column name of the sampling indicator.
#' @param target_col Numeric or integer vector. The sampling indicator values.
#' @param target_subset Character. `"target"` or `"all"`.
#' @param estimator Character. The estimator chosen in `causat()`.
#' @param call Caller environment for error messages.
#'
#' @noRd
check_transport_inputs <- function(
  target,
  target_col,
  target_subset = "target",
  estimator = "gcomp",
  estimand = "ATE",
  call = rlang::caller_env()
) {
  if (estimator == "matching") {
    rlang::abort(
      c(
        "`target` is not supported with `estimator = \"matching\"`.",
        i = paste0(
          "Transportability does not compose cleanly with matching. ",
          "Use `estimator = \"gcomp\"` or `\"ipw\"` instead."
        )
      ),
      class = "causatr_transport_matching",
      call = call
    )
  }

  if (estimand %in% c("ATT", "ATC")) {
    rlang::abort(
      c(
        paste0(
          "`estimand = \"",
          estimand,
          "\"` is not supported with `target` (transportability)."
        ),
        i = paste0(
          "Transport estimands are defined over the target population, ",
          "not the treated or control subgroup. Use `estimand = \"ATE\"`."
        )
      ),
      class = "causatr_transport_estimand",
      call = call
    )
  }

  if (anyNA(target_col)) {
    rlang::abort(
      c(
        paste0("Column `", target, "` contains NA values."),
        i = "The sampling indicator S must be complete (no NAs) for transportability."
      ),
      class = "causatr_transport_na_target",
      call = call
    )
  }

  vals <- sort(unique(as.integer(target_col)))
  if (!all(vals %in% c(0L, 1L))) {
    rlang::abort(
      c(
        paste0(
          "`",
          target,
          "` must be binary (0 = target, 1 = study). ",
          "Found values: ",
          paste(vals, collapse = ", "),
          "."
        ),
        i = "Recode the sampling indicator to 0/1 before calling `causat()`."
      ),
      class = "causatr_transport_non_binary",
      call = call
    )
  }

  if (!all(c(0L, 1L) %in% vals)) {
    rlang::abort(
      c(
        paste0(
          "`",
          target,
          "` must have both study (S = 1) and target (S = 0) rows. ",
          "Found only: ",
          paste(vals, collapse = ", "),
          "."
        ),
        i = "The sampling model requires rows from both populations."
      ),
      class = "causatr_transport_degenerate",
      call = call
    )
  }

  p_study <- mean(as.integer(target_col) == 1L)
  if (p_study > 0.95 || p_study < 0.05) {
    rlang::warn(
      c(
        paste0(
          round(100 * p_study, 1),
          "% of observations are study rows (S = 1). ",
          "Sampling weights may be extreme and estimates may be unstable."
        ),
        i = "Consider whether the sampling model is correctly specified."
      ),
      class = "causatr_transport_extreme_selection"
    )
  }

  invisible(NULL)
}


#' Resolve `cluster` argument to a vector aligned with the sandwich IF rows
#'
#' @description
#' Validates the user-supplied `cluster` argument to [contrast()] and
#' returns the cluster vector (or `NULL`) that the variance engine will
#' hand to `vcov_from_if(cluster = ...)`. Three shapes are accepted:
#'
#' - `NULL`: standard (ungrouped) sandwich aggregation.
#' - `"col"`: name of a column in `fit$data`; the vector is extracted
#'   from that column.
#' - numeric / character / factor of length `nrow(fit$data)`: used
#'   directly.
#'
#' NAs in the cluster vector are rejected upfront — `rowsum()` would
#' silently treat each NA as its own implicit group (via `factor()`
#' levels) and quietly produce an aggregated vcov that differs from a
#' user-specified cluster where the NA handling is explicit.
#'
#' Matching is intentionally rejected here rather than silently ignored,
#' because `variance_if_matching()` already clusters on `subclass` and a
#' user-supplied cluster would either conflict or double-count. Users
#' who genuinely need a design cluster outside the matched structure
#' should estimate on `gcomp` or `ipw`.
#'
#' @param cluster `NULL`, column name string, or length-n vector.
#' @param fit A `causatr_fit` object.
#' @param call Caller environment for error messages.
#'
#' @return `NULL` (no clustering) or a length-n vector suitable for
#'   `vcov_from_if(cluster = ...)`.
#'
#' @noRd
resolve_cluster <- function(cluster, fit, call = rlang::caller_env()) {
  if (is.null(cluster)) {
    return(NULL)
  }

  if (fit$estimator == "matching") {
    rlang::abort(
      c(
        "`cluster` is not supported for `estimator = \"matching\"`.",
        i = paste0(
          "Matching already aggregates IFs cluster-robustly on the ",
          "matched `subclass`. Adding a user-supplied cluster on top ",
          "would double-count within-subclass rows."
        ),
        i = "Use `estimator = \"gcomp\"` or `\"ipw\"` for a design cluster."
      ),
      class = "causatr_bad_cluster",
      call = call
    )
  }

  n <- nrow(fit$data)

  # Accept a column name or a length-n vector. Anything else (e.g. a
  # formula, a vector with the wrong length, a list) lands in the
  # final abort below.
  if (rlang::is_string(cluster)) {
    if (!cluster %in% names(fit$data)) {
      rlang::abort(
        paste0(
          "`cluster` column '",
          cluster,
          "' not found in `fit$data`."
        ),
        class = "causatr_bad_cluster",
        call = call
      )
    }
    vec <- fit$data[[cluster]]
  } else if (
    (is.numeric(cluster) || is.character(cluster) || is.factor(cluster)) &&
      length(cluster) == n
  ) {
    vec <- cluster
  } else {
    rlang::abort(
      c(
        "`cluster` must be NULL, a column name in `fit$data`, or a vector matching `nrow(fit$data)`.",
        i = paste0("Got length ", length(cluster), " vs n = ", n, ".")
      ),
      class = "causatr_bad_cluster",
      call = call
    )
  }

  if (anyNA(vec)) {
    rlang::abort(
      c(
        "`cluster` contains NA values.",
        i = paste0(
          "Clustering requires every row to belong to an identifiable ",
          "cluster. Drop NA rows or provide an explicit singleton id."
        )
      ),
      class = "causatr_bad_cluster",
      call = call
    )
  }

  vec
}
