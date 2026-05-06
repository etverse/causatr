#' Column names reserved for causatr-internal use
#'
#' @description
#' Used internally by ice() for the pseudo-outcome column. User data that
#' already carries any of these columns would silently collide with the
#' internal mutations, so we reject them up front via
#' `check_reserved_cols()`.
#'
#' @noRd
CAUSATR_RESERVED_COLS <- c(".pseudo_y")

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

#' Construct a `causatr_fit` object
#'
#' @param model Fitted model object (glm, gam, glm_weightit, etc.) or `NULL`
#'   (ICE defers fitting to `contrast()`).
#' @param data data.table of the full dataset.
#' @param treatment Character treatment column name(s).
#' @param outcome Character outcome column name.
#' @param confounders One-sided formula of baseline confounders.
#' @param confounders_tv One-sided formula of time-varying confounders or `NULL`.
#' @param family Character or family object describing the outcome distribution.
#' @param estimator Character causal estimator (`"gcomp"`, `"ipw"`, `"matching"`).
#' @param type `"point"` or `"longitudinal"`.
#' @param estimand Character estimand (`"ATE"`, `"ATT"`, `"ATC"`).
#' @param id Character ID column name or `NULL`.
#' @param time Character time column name or `NULL`.
#' @param censoring Character censoring column name or `NULL`.
#' @param history Integer Markov order for longitudinal.
#' @param numerator One-sided formula for stabilised weights or `NULL`.
#' @param weights_obj Legacy slot retained for `causatr_fit` backward
#'   compatibility; always `NULL` now that IPW is self-contained.
#' @param match_obj A `matchit` object (matching) or `NULL`.
#' @param call The original `causat()` call environment.
#' @param details Named list of estimator-specific metadata.
#' @return A list with class `"causatr_fit"`.
#' @noRd
new_causatr_fit <- function(
  model,
  data,
  treatment,
  outcome,
  confounders,
  confounders_tv,
  family,
  estimator,
  type,
  estimand,
  id,
  time,
  censoring,
  history,
  numerator,
  weights_obj,
  match_obj,
  call,
  details
) {
  structure(
    list(
      model = model,
      data = data,
      treatment = treatment,
      outcome = outcome,
      confounders = confounders,
      confounders_tv = confounders_tv,
      family = family,
      estimator = estimator,
      type = type,
      estimand = estimand,
      id = id,
      time = time,
      censoring = censoring,
      history = history,
      numerator = numerator,
      weights_obj = weights_obj,
      match_obj = match_obj,
      call = call,
      details = details
    ),
    class = "causatr_fit"
  )
}

#' Construct a `causatr_result` object
#'
#' @param estimates data.table of intervention-specific marginal means.
#' @param contrasts data.table of pairwise contrasts with SEs and CIs.
#' @param type Character contrast type (`"difference"`, `"ratio"`, `"or"`).
#' @param estimand Character estimand used.
#' @param ci_method Character CI method (`"sandwich"` or `"bootstrap"`).
#' @param reference Character name of the reference intervention or `NULL`.
#' @param interventions Named list of `causatr_intervention` objects.
#' @param n Integer sample size used for estimation.
#' @param estimator Character causal estimator.
#' @param vcov Variance-covariance matrix of marginal means.
#' @param call The original `contrast()` call environment.
#' @return A list with class `"causatr_result"`.
#' @noRd
new_causatr_result <- function(
  estimates,
  contrasts,
  type,
  estimand,
  ci_method,
  reference,
  interventions,
  n,
  estimator,
  family,
  fit_type,
  vcov,
  boot_t = NULL,
  boot_info = NULL,
  call
) {
  structure(
    list(
      estimates = estimates,
      contrasts = contrasts,
      type = type,
      estimand = estimand,
      ci_method = ci_method,
      reference = reference,
      interventions = interventions,
      n = n,
      estimator = estimator,
      family = family,
      fit_type = fit_type,
      vcov = vcov,
      boot_t = boot_t,
      # NULL when ci_method = "sandwich"; otherwise a 3-element list of
      # `n_requested`, `n_ok`, `n_fail` carried up from
      # `process_boot_results()` so `print.causatr_result()` and
      # downstream consumers can surface bootstrap failure rates without
      # re-deriving them from `boot_t`.
      boot_info = boot_info,
      call = call
    ),
    class = "causatr_result"
  )
}

#' Construct a `causatr_diag` object
#'
#' @description
#' Builds the nested per-intervention diagnostic container. Each entry in
#' `per_intervention` is a named list with `positivity`, `balance`, and
#' `weights` slots (any of which may be `NULL` when the underlying estimator
#' or treatment family does not support that diagnostic). The top-level
#' `positivity` / `balance` / `weights` slots are populated from the first
#' panel for backward compatibility with the flat shape that downstream
#' callers (`print.causatr_diag()`, `plot.causatr_diag()`, existing tests)
#' originally consumed.
#'
#' @param per_intervention Named list, one entry per intervention key, each
#'   itself a list with slots `positivity`, `balance`, `weights`.
#' @param match_quality data.table of match quality metrics or `NULL`. Lives
#'   at the top level rather than per-intervention because matching is done
#'   once at fit time and is intervention-agnostic.
#' @param estimator Character causal estimator.
#' @param fit_info Named list with summary metadata about the fit
#'   (`treatment_type`, `estimand`, `type`, `has_em`). Used by the print
#'   method to render the panel header.
#' @param fit The original `causatr_fit` (stored for `plot()` method).
#' @return A list with class `"causatr_diag"`.
#' @noRd
new_causatr_diag <- function(
  per_intervention,
  match_quality = NULL,
  estimator,
  fit_info = list(),
  fit = NULL
) {
  if (!is.list(per_intervention) || length(per_intervention) == 0L) {
    rlang::abort(
      "Internal error: `per_intervention` must be a non-empty named list."
    )
  }
  if (
    is.null(names(per_intervention)) || any(!nzchar(names(per_intervention)))
  ) {
    rlang::abort(
      "Internal error: every panel in `per_intervention` must be named."
    )
  }
  # First panel feeds the backward-compat top-level slots. The flat
  # access pattern (`diag$positivity`, `diag$balance`, `diag$weights`)
  # was the public API before the per-intervention redesign; preserving it lets every
  # existing test, print-tests, and downstream user expression keep
  # working unchanged when no `interventions =` argument was passed.
  first <- per_intervention[[1L]]
  structure(
    list(
      per_intervention = per_intervention,
      interventions = names(per_intervention),
      positivity = first$positivity,
      balance = first$balance,
      weights = first$weights,
      censoring = first$censoring,
      match_quality = match_quality,
      estimator = estimator,
      fit_info = fit_info,
      fit = fit
    ),
    class = "causatr_diag"
  )
}

#' Construct a `causatr_intervention` object
#'
#' @param type Character intervention type (e.g. `"static"`, `"shift"`).
#' @param params Named list of intervention parameters.
#' @return A list with class `"causatr_intervention"`.
#' @noRd
new_causatr_intervention <- function(type, params) {
  structure(
    c(list(type = type), params),
    class = "causatr_intervention"
  )
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

#' Resolve a family argument to a family object
#'
#' Accepts a character string (e.g. `"gaussian"`), a family function
#' (e.g. `stats::binomial`), or an already-evaluated family object
#' (e.g. `stats::binomial()`).
#'
#' @param family Character, function, or family object.
#' @return A family object (list with `$family`, `$link`, etc.).
#' @noRd
resolve_family <- function(family) {
  # Three-way dispatch to canonicalize every input to a family object.
  # Downstream code universally relies on `family$family` (character)
  # and `family$link` (character) being present -- normalizing here
  # means each call site doesn't need its own isinstance check.
  if (is.character(family)) {
    # Look up in stats:: namespace explicitly -- avoids picking up a
    # user-defined function with the same name in the caller's env.
    # Wrap in tryCatch because some stats:: functions share names with
    # family constructors (e.g. `stats::beta` is the beta distribution
    # function, not a GLM family; calling it returns a numeric value,
    # not a family object, and the `inherits(., "family")` guard below
    # catches that).
    fam_obj <- tryCatch(
      {
        fam_fn <- get(family, mode = "function", envir = asNamespace("stats"))
        fam_fn()
      },
      error = function(e) NULL
    )
    if (!is.null(fam_obj) && inherits(fam_obj, "family")) {
      return(fam_obj)
    }

    # Extended families from optional packages.
    if (family == "beta") {
      rlang::check_installed(
        "betareg",
        reason = "for beta regression outcomes (family = \"beta\")"
      )
      # `betareg::betareg()` doesn't accept a `family` argument; the beta
      # family object built here is used only by `is_binary_family()` (returns
      # FALSE so the non-binary estimation path is taken) and by variance
      # tier selection. The link slots are populated so downstream code that
      # calls `family$linkinv()` on predictions doesn't error.
      logit <- stats::make.link("logit")
      return(structure(
        list(
          family = "beta",
          link = "logit",
          linkfun = logit$linkfun,
          linkinv = logit$linkinv,
          mu.eta = logit$mu.eta,
          variance = function(mu) mu * (1 - mu)
        ),
        class = "family"
      ))
    }

    rlang::abort(paste0(
      "Unknown family: '",
      family,
      "'. ",
      "Supported: gaussian, binomial, poisson, quasibinomial, ",
      "quasipoisson, Gamma, inverse.gaussian, beta. ",
      "Or pass a family object directly."
    ))
  }
  if (is.function(family)) {
    return(family())
  }
  # Already a family object (list with $family, $link, etc.).
  family
}

#' Does a fitting function accept a `family` argument?
#'
#' Inspects the formal arguments of `fn` to determine whether it accepts
#' a named `family` parameter. Used to conditionally strip `family` from
#' the argument list when calling model functions that handle their own
#' family internally (e.g. `MASS::glm.nb`, `betareg::betareg`).
#'
#' @param fn A function (the user's `model_fn`).
#' @return Logical scalar.
#' @noRd
fn_accepts_family <- function(fn) {
  # `formals()` returns the declared parameter list without executing the
  # function; it's the correct way to inspect the signature of an arbitrary
  # function object. `MASS::glm.nb` and `betareg::betareg` do not declare
  # `family` in their formals, so this returns FALSE for both, causing the
  # fitter to omit the `family` argument rather than passing one those
  # functions would reject with a cryptic error.
  "family" %in% names(formals(fn))
}

#' Check whether a family describes a binary outcome
#'
#' @param family A character string, family object, or function.
#' @return Logical scalar.
#' @noRd
is_binary_family <- function(family) {
  if (is.null(family)) {
    return(FALSE)
  }
  fam <- tryCatch(resolve_family(family), error = function(e) NULL)
  if (!is.null(fam) && is.character(fam$family)) {
    return(fam$family %in% c("binomial", "quasibinomial"))
  }
  FALSE
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
get_fit_rows <- function(data, outcome, censoring = NULL) {
  # Both conditions must hold for a row to enter the model fit:
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
  is_uncensored(data, censoring) & !is.na(data[[outcome]])
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

#' Swap binomial for quasibinomial in weighted MSM fits
#'
#' @description
#' IPW outcome MSMs are fit with density-ratio weights, which are
#' non-integer. R's `glm.fit` warns "non-integer #successes in a
#' binomial glm!" for binomial families with non-integer weights.
#' quasibinomial produces identical coefficients, SEs, and
#' predictions but suppresses the warning.
#'
#' @param fam A family object (as returned by `resolve_family()`).
#' @return The same family, or `quasibinomial()` if the input was
#'   `binomial`.
#' @noRd
msm_family <- function(fam) {
  # IPW MSMs are fit with density-ratio weights, which are non-integer
  # fractions. `stats::glm` emits "non-integer #successes in a binomial
  # glm!" for every row when `family = binomial`. `quasibinomial` produces
  # identical coefficients, SEs, and predictions because it uses the same
  # IRLS update but skips the Pearson chi-squared check that triggers the
  # warning. The link is forwarded to preserve user-specified link choices
  # (e.g. `binomial("log")` -> `quasibinomial("log")`).
  if (is.list(fam) && identical(fam$family, "binomial")) {
    return(stats::quasibinomial(link = fam$link))
  }
  fam
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
