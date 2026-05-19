#' Column names reserved for causatr-internal use
#'
#' @description
#' Used internally by ice() for the pseudo-outcome column. User data that
#' already carries any of these columns would silently collide with the
#' internal mutations, so we reject them up front via
#' `check_reserved_cols()`.
#'
#' @noRd
CAUSATR_RESERVED_COLS <- c(".pseudo_y", ".sampling_s")

#' Safely replay a fit function with stashed `...` arguments
#'
#' @description
#' Centralises the "combine base args + user `...` then `do.call()`"
#' pattern used by `ice_iterate()`, `refit_gcomp()`, `refit_ipw()`, and
#' `refit_matching()`. Hand-written `c(args, dots)` compositions repeated
#' across four sites had three known footguns (see the 2026-04-15
#' third-round dots audit):
#'
#'  - **R1** (T2): positional / unnamed entries in `dots` survive the
#'    `setdiff(names(dots), ...)` strip and collide with NSE in the
#'    target function's first argument (e.g. `stats::glm`'s `formula`).
#'  - **R2** (C5): duplicate named args silently take the first via
#'    `do.call()`, which is correct today but drifts if any of the four
#'    strip lists falls out of sync.
#'  - **R4**: the four strip lists are manually maintained and can
#'    diverge; adding a new reserved key to one site but forgetting
#'    another is a silent regression.
#'
#' This helper collapses all three into one place: unnamed entries are
#' dropped, any key already present in `base_args` (or explicitly listed
#' in `reserved`) is stripped from `dots`, and the rest flow through
#' `do.call()`.
#'
#' @param fn The function to call (e.g. `stats::glm`, `MatchIt::matchit`,
#'   `mgcv::gam`).
#' @param base_args Named list of arguments the caller is supplying
#'   explicitly. These always win over `dots` in case of name collision.
#'   Unnamed elements in `base_args` are preserved (they are meant to be
#'   positional, e.g. the first `formula` argument to `glm`).
#' @param dots User-supplied `...` list captured via `list(...)` at fit
#'   time and stashed on `fit$details$dots`. Can be `NULL`.
#' @param reserved Character vector of additional argument names to
#'   strip from `dots` -- useful for target-function parameters the
#'   caller sets implicitly (e.g. `MatchIt::matchit`'s `distance`, set
#'   by `refit_matching()` after the composition).
#'
#' @return The return value of `do.call(fn, ...)`.
#'
#' @noRd
replay_fit <- function(fn, base_args, dots = NULL, reserved = NULL) {
  dots <- dots %||% list()
  if (length(dots) > 0L) {
    # Drop unnamed entries (T2): `do.call(glm, list(formula, ..., ""=1))`
    # is unpredictable because glm uses NSE on its first (unnamed)
    # argument. Users should not be passing positional dots to
    # `causat()`, but defending against it is cheap.
    keep_named <- nzchar(names(dots))
    dots <- dots[keep_named]
    # Drop keys the caller set explicitly. `do.call()` takes the first
    # named duplicate, so today this is belt-and-braces -- but forgetting
    # the strip at a new call site would silently pick up the
    # user-supplied value and ignore the caller's.
    blocked <- c(names(base_args), reserved)
    if (length(blocked) > 0L) {
      dots <- dots[setdiff(names(dots), blocked)]
    }
  }
  do.call(fn, c(base_args, dots))
}

#' Assert that `data` does not carry any causatr-reserved column names
#'
#' @param data A data.frame / data.table.
#' @param which Optional subset of `CAUSATR_RESERVED_COLS` to check. Defaults
#'   to the full set.
#' @return `NULL` invisibly; aborts with a specific error listing the
#'   offending columns when one is present.
#' @noRd
check_reserved_cols <- function(data, which = CAUSATR_RESERVED_COLS) {
  bad <- intersect(which, names(data))
  if (length(bad) > 0L) {
    rlang::abort(
      paste0(
        "Column name(s) ",
        paste0("`", bad, "`", collapse = ", "),
        " are reserved by causatr internals. Rename the column(s) in your ",
        "input data before calling `causat()`."
      )
    )
  }
  invisible(NULL)
}


#' Check whether rows are uncensored
#'
#' A row is uncensored when the censoring column is `NA` or `0`.
#' Any other non-NA value (e.g. `1` for censored, `2` for a competing event)
#' is treated as censored.
#'
#' @param data A data.table.
#' @param censoring Character censoring column name, or `NULL`.
#' @return Logical vector of length `nrow(data)` (`TRUE` = uncensored).
#' @noRd
is_uncensored <- function(data, censoring) {
  # `censoring = NULL` is the "no censoring" shortcut -- every row is
  # treated as uncensored. This is the default for cross-sectional data.
  if (is.null(censoring)) {
    return(rep(TRUE, nrow(data)))
  }
  cens <- data[[censoring]]
  # NA censoring value is treated as uncensored: in typical
  # longitudinal datasets, C is only defined up to the moment of
  # dropout, and subsequent rows have NA. Treating NA as censored
  # would silently drop everyone from the fitting set.
  is.na(cens) | cens == 0L
}


#' Build a propensity score formula from confounders and treatment
#'
#' Strips any terms that involve the treatment variable before building
#' the two-sided formula `A ~ confounders`. When `confounders` includes
#' effect-modification terms like `A:sex`, those are excluded from the
#' propensity RHS (the treatment cannot appear on both sides). The
#' pure confounder terms (`L`, `sex`, `I(age^2)`, `L1:L2`) pass through.
#'
#' @param confounders One-sided formula of confounders.
#' @param treatment Character vector of treatment column name(s).
#' @return A two-sided formula: `treatment ~ confounder_terms`.
#' @noRd
build_ps_formula <- function(confounders, treatment) {
  em_info <- parse_effect_mod(confounders, treatment)
  # Use only the non-EM terms for the propensity model RHS. Including
  # `A:sex` on the RHS of a model for A would be circular (A appears on
  # both sides). Modifier main effects that appear as standalone terms
  # (e.g. `sex` in `~ L + sex + A:sex`) are already in `confounder_terms`
  # because `parse_effect_mod` only classifies the interaction itself as EM.
  ps_terms <- em_info$confounder_terms
  if (length(ps_terms) == 0L) {
    # Edge case: confounders formula had only EM terms (e.g. `~ A:sex`).
    # `reformulate("1", ...)` gives intercept-only, which is valid.
    ps_terms <- "1"
  }
  stats::reformulate(ps_terms, response = treatment)
}

#' Build an outcome model formula: response ~ treatment + confounders
#'
#' Extracts RHS term labels from `confounders` (preserving user-supplied
#' transforms such as `ns(age, 4)` or `L1*L2`), prepends the treatment
#' column name(s), and wraps the result in `stats::reformulate()`.
#'
#' @param response Character. LHS response column name.
#' @param treatment Character vector of treatment column name(s).
#' @param confounders One-sided formula.
#' @return Two-sided formula `response ~ treatment + confounder_terms`.
#' @noRd
build_outcome_formula <- function(response, treatment, confounders) {
  confounder_terms <- attr(stats::terms(confounders), "term.labels")
  stats::reformulate(c(treatment, confounder_terms), response = response)
}

#' Resolve per-component confounders from a causatr_fit object
#'
#' Each resolver returns the component-specific formula if set, falling
#' back to the deprecated \code{confounders} / \code{confounders_tv}
#' slot for backward compatibility.
#'
#' @param fit A \code{causatr_fit} object.
#' @return A one-sided formula or \code{NULL}.
#' @noRd
resolve_confounders_outcome <- function(fit) {
  fit$confounders_outcome %||% fit$confounders
}

#' @rdname resolve_confounders_outcome
#' @noRd
resolve_confounders_treatment <- function(fit) {
  fit$confounders_treatment %||% fit$confounders
}

#' @rdname resolve_confounders_outcome
#' @noRd
resolve_confounders_censoring <- function(fit) {
  fit$confounders_censoring %||% fit$confounders
}

#' @rdname resolve_confounders_outcome
#' @noRd
resolve_confounders_sampling <- function(fit) {
  fit$confounders_sampling %||% fit$confounders
}

#' @rdname resolve_confounders_outcome
#' @noRd
resolve_confounders_tv_outcome <- function(fit) {
  fit$confounders_tv_outcome %||% fit$confounders_tv
}

#' @rdname resolve_confounders_outcome
#' @noRd
resolve_confounders_tv_treatment <- function(fit) {
  fit$confounders_tv_treatment %||% fit$confounders_tv
}

#' Build lag-prefixed column names
#'
#' For every variable in `base_vars` and every lag order `k` in
#' `1:n_lags`, produces the string `"lag{k}_{var}"`.  Returns
#' `character(0)` when `n_lags == 0`.
#'
#' @param base_vars Character vector of base column names.
#' @param n_lags Non-negative integer number of lag orders.
#' @return Character vector of length `length(base_vars) * n_lags`.
#' @noRd
build_lag_terms <- function(base_vars, n_lags) {
  if (n_lags == 0L) {
    return(character(0L))
  }
  unlist(lapply(seq_len(n_lags), function(k) paste0("lag", k, "_", base_vars)))
}

#' Parse and validate confounders for treatment-touching terms
#'
#' Delegates to `parse_effect_mod()` to detect effect-modification terms,
#' then calls `check_em_compat()` to reject bare treatment terms (always
#' invalid). Returns the parsed `causatr_em_info` so callers can decide
#' whether to support or reject the EM terms for their method.
#'
#' @param confounders One-sided formula passed by the user.
#' @param treatment Character vector of treatment column name(s).
#' @param estimator Character. `"ipw"` or `"matching"`; used by downstream
#'   callers for method-specific validation.
#' @param data Optional data.frame/data.table for method-specific checks.
#' @return A `causatr_em_info` object (invisibly). Aborts if bare treatment
#'   terms are found.
#' @noRd
check_confounders_treatment <- function(
  confounders,
  treatment,
  estimator = NULL,
  data = NULL
) {
  em_info <- parse_effect_mod(confounders, treatment)
  # Bare treatment in confounders (e.g. `~ L + A`) is always invalid:
  # it puts A on both sides of the propensity model.
  check_em_compat(em_info, treatment, estimator, data)
  em_info
}

#' Get the logical vector of fitting rows
#'
#' A row is eligible for model fitting when it is uncensored AND has a
#' non-missing outcome.
#'
#' @param data A data.table.
#' @param outcome Character outcome column name.
#' @param censoring Character censoring column name, or `NULL`.
#' @return Logical vector of length `nrow(data)`.
#' @noRd
get_fit_rows <- function(data, outcome, censoring = NULL, target = NULL) {
  # Three conditions must hold for a row to enter the model fit:
  #
  # 1. Uncensored: a censored row carries an observed outcome that
  #    reflects censoring rather than the natural disease process;
  #    including it without IPCW weighting biases E[Y^a].
  #
  # 2. Non-missing outcome: rows with NA outcome would be silently
  #    dropped by glm's default `na.action = na.omit`, which shortens
  #    the fitted-value vector and misaligns `fit_rows`-indexed score
  #    matrices in the sandwich engine. Excluding them here keeps the
  #    bookkeeping correct and makes the exclusion explicit.
  #
  # 3. Study population (transport only): when `target` is non-NULL,
  #    restrict to S = 1 rows. Target rows (S = 0) have no observed Y
  #    or A and must not enter the outcome model fit.
  rows <- is_uncensored(data, censoring) & !is.na(data[[outcome]])
  if (!is.null(target)) {
    rows <- rows & (data[[target]] == 1L)
  }
  rows
}

#' Weighted or unweighted mean
#'
#' @param x Numeric vector.
#' @param w Numeric weight vector or `NULL` for unweighted mean.
#' @return Scalar mean.
#' @noRd
maybe_weighted_mean <- function(x, w = NULL) {
  if (!is.null(w)) {
    stats::weighted.mean(x, w, na.rm = TRUE)
  } else {
    mean(x, na.rm = TRUE)
  }
}

#' Check that an optional package is installed
#'
#' @param pkg Character package name.
#' @return `NULL` invisibly; aborts with installation instructions if absent.
#' @noRd
check_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    rlang::abort(
      paste0(
        "Package '",
        pkg,
        "' is required but not installed. ",
        "Install it with: install.packages('",
        pkg,
        "')"
      ),
      .call = FALSE
    )
  }
}


#' Check whether an intervention is stochastic
#'
#' @param iv A `causatr_intervention` object or `NULL`.
#' @return Logical scalar.
#' @noRd
is_stochastic_intervention <- function(iv) {
  !is.null(iv) &&
    inherits(iv, "causatr_intervention") &&
    iv$type == "stochastic"
}

#' Check whether an intervention (scalar or multivariate list) has any
#' stochastic component
#'
#' @param iv A `causatr_intervention`, a named list of them, or `NULL`.
#' @return Logical scalar.
#' @noRd
has_stochastic_component <- function(iv) {
  if (is.null(iv)) {
    return(FALSE)
  }
  if (inherits(iv, "causatr_intervention")) {
    return(iv$type == "stochastic")
  }
  if (is.list(iv)) {
    return(any(vapply(iv, is_stochastic_intervention, logical(1))))
  }
  FALSE
}

#' Extract n_mc from an intervention (scalar or multivariate list)
#'
#' For multivariate interventions, returns the maximum `n_mc` across
#' all stochastic components (all components use the same MC loop).
#'
#' @param iv A `causatr_intervention` or named list of them.
#' @return Integer `n_mc`, or `NULL` if no stochastic component.
#' @noRd
get_stochastic_n_mc <- function(iv) {
  if (is_stochastic_intervention(iv)) {
    return(iv$n_mc)
  }
  if (is.list(iv) && !inherits(iv, "causatr_intervention")) {
    nms <- vapply(
      iv,
      function(x) {
        if (is_stochastic_intervention(x)) x$n_mc else 0L
      },
      integer(1)
    )
    m <- max(nms)
    if (m > 0L) return(m)
  }
  NULL
}

#' Parallel-aware sapply for Monte Carlo loops
#'
#' Uses `future.apply::future_sapply()` when a non-sequential `future`
#' plan is registered; otherwise falls back to `vapply()`. The
#' `future.seed = TRUE` argument ensures reproducible RNG streams across
#' workers (L'Ecuyer-CMRG).
#'
#' @param X Integer sequence (e.g. `seq_len(n_mc)`).
#' @param FUN Function taking one element of `X` and returning a numeric
#'   vector of length `n_rows`.
#' @param n_rows Integer. Expected length of each FUN return value.
#' @return A numeric matrix of dimension `n_rows x length(X)`.
#' @noRd
mc_sapply <- function(X, FUN, n_rows) {
  use_future <- requireNamespace("future.apply", quietly = TRUE) &&
    requireNamespace("future", quietly = TRUE) &&
    !inherits(future::plan(), "sequential")
  if (use_future) {
    future.apply::future_sapply(
      X,
      FUN,
      future.seed = TRUE,
      simplify = TRUE
    )
  } else {
    vapply(X, FUN, numeric(n_rows))
  }
}
