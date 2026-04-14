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
#' @param weights_obj A `weightit` object (IPW) or `NULL`.
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
#' @param balance Balance table (from cobalt or simple SMD computation).
#' @param positivity data.table of propensity score summaries.
#' @param weights data.table of weight distribution summaries (IPW) or `NULL`.
#' @param match_quality data.table of match quality metrics or `NULL`.
#' @param estimator Character causal estimator.
#' @param fit The original `causatr_fit` (stored for `plot()` method).
#' @return A list with class `"causatr_diag"`.
#' @noRd
new_causatr_diag <- function(
  balance,
  positivity,
  weights,
  match_quality,
  estimator,
  fit = NULL
) {
  structure(
    list(
      balance = balance,
      positivity = positivity,
      weights = weights,
      match_quality = match_quality,
      estimator = estimator,
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
  # `censoring = NULL` is the "no censoring" shortcut — every row is
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
  # This matches what glm/gam do internally so downstream code can
  # assume `family$family` and `family$link` are strings.
  if (is.character(family)) {
    # Look up in stats:: namespace explicitly — avoids picking up a
    # user-defined function with the same name in the caller's env.
    fam_fn <- tryCatch(
      get(family, mode = "function", envir = asNamespace("stats")),
      error = function(e) {
        rlang::abort(paste0("Unknown family: '", family, "'."))
      }
    )
    return(fam_fn())
  }
  if (is.function(family)) {
    return(family())
  }
  # Already a family object (list with $family, $link, etc.).
  family
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
#' @param confounders One-sided formula of confounders.
#' @param treatment Character treatment column name.
#' @return A two-sided formula: `treatment ~ confounder_terms`.
#' @noRd
build_ps_formula <- function(confounders, treatment) {
  # Turn `~ L1 + L2 + I(age^2)` (confounders) and `"A"` (treatment)
  # into `A ~ L1 + L2 + I(age^2)`. `term.labels` preserves user
  # transformations and interactions, which `reformulate()` then
  # reassembles verbatim into a two-sided formula.
  confounder_terms <- attr(stats::terms(confounders), "term.labels")
  stats::reformulate(confounder_terms, response = treatment)
}

#' Abort if `confounders` contains a term involving the treatment
#'
#' IPW and matching fitters wrap a hardcoded saturated MSM (`Y ~ A`) around a
#' `confounders`-driven propensity/match model. Any interaction term the user
#' writes as `A:modifier` (or `A * modifier`) cannot land anywhere useful:
#' it pollutes the PS formula with the treatment on both sides, and it is
#' silently dropped from the MSM, which `contrast()` averages over. That is
#' the core Phase 8 limitation (see `PHASE_8_INTERACTIONS.md`). Until Phase
#' 8 lands we abort early with a pointer to `estimator = "gcomp"` so users
#' do not silently get a homogeneous-effect estimate when they asked for a
#' heterogeneous one.
#'
#' Bare treatment terms in `confounders` (e.g. `~ L + A`) are also rejected
#' because `A` has no place in a propensity model of `A`.
#'
#' @param confounders One-sided formula passed by the user.
#' @param treatment Character vector of treatment column name(s).
#' @param estimator Character. `"ipw"` or `"matching"`; used in the error text.
#' @return `invisible(NULL)` on success; aborts otherwise.
#' @noRd
check_confounders_no_treatment <- function(confounders, treatment, estimator) {
  term_labels <- attr(stats::terms(confounders), "term.labels")
  if (length(term_labels) == 0L) {
    return(invisible(NULL))
  }

  # For each term, extract the variables it references via `all.vars()` on
  # the parsed expression. This catches `A:sex`, `sex:A`, `A*sex` (which
  # `terms()` has already expanded to main effects + interaction), and bare
  # `A`. It does not catch `I(A + 1)` on purpose: `I()` wraps are opaque to
  # `terms()` and users reaching for them are doing something advanced.
  offenders <- vapply(
    term_labels,
    function(tl) {
      vars <- all.vars(parse(text = tl)[[1L]])
      any(vars %in% treatment)
    },
    logical(1L)
  )

  if (any(offenders)) {
    bad <- term_labels[offenders]
    rlang::abort(
      c(
        paste0(
          "`confounders` contains term(s) involving the treatment, which ",
          "are not supported for `estimator = \"",
          estimator,
          "\"`."
        ),
        x = paste0("Offending term(s): ", paste(bad, collapse = ", "), "."),
        i = paste0(
          "IPW and matching wrap a saturated MSM `Y ~ A` around the ",
          "propensity/match model, so treatment-by-modifier interactions ",
          "cannot be estimated here."
        ),
        i = paste0(
          "Use `estimator = \"gcomp\"` for heterogeneous treatment effects, ",
          "or `by = \"modifier\"` in `contrast()` for stratum-specific ",
          "summaries of a homogeneous effect."
        ),
        i = "See `PHASE_8_INTERACTIONS.md` for the planned unified API."
      ),
      .call = FALSE
    )
  }
  invisible(NULL)
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
