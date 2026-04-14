#' Fit a matching-based causal model (point treatment)
#'
#' @description
#' Implements propensity-score matching (Hernán & Robins Ch. 15) by
#' delegating match construction to `MatchIt::matchit()` and fitting the
#' outcome model on the matched sample.
#'
#' ## Algorithm
#'
#' 1. Compute propensity scores and create matched sets via
#'    `MatchIt::matchit()`.
#' 2. Extract the matched data and match weights via
#'    `MatchIt::match.data()`.
#' 3. Fit the outcome model on the matched data: `Y ~ A` (or with
#'    covariates for doubly-robust estimation), weighted by the match
#'    weights.
#' 4. Inference uses cluster-robust SE via `sandwich::vcovCL()` with
#'    clustering on the matched-pair subclass.
#'
#' @param data data.table from `prepare_data()`.
#' @param outcome Character. Outcome column name.
#' @param treatment Character. Treatment column name (scalar only).
#' @param confounders One-sided formula of baseline confounders.
#' @param estimand Character. `"ATE"`, `"ATT"`, or `"ATC"`.
#' @param type Character. `"point"` or `"longitudinal"`.
#' @param weights Numeric vector or `NULL`. External observation weights.
#' @param call The original `causat()` call.
#' @param ... Passed to `MatchIt::matchit()` (e.g. `method = "nearest"`,
#'   `ratio = 1`, `caliper = 0.2`).
#'
#' @return A `causatr_fit` object with `match_obj` (the `matchit` object)
#'   and `model` (the `glm` fit on matched data).
#'
#' @noRd
fit_matching <- function(
  data,
  outcome,
  treatment,
  confounders,
  family,
  estimand,
  type,
  weights,
  call,
  ...
) {
  if (type == "longitudinal") {
    rlang::abort(
      "Matching is only supported for point treatments.",
      .call = FALSE
    )
  }

  # Build the treatment model formula: A ~ confounders.
  # MatchIt uses this to estimate propensity scores for matching.
  ps_formula <- build_ps_formula(confounders, treatment)

  # Fit rows: exclude missing outcomes before matching.
  fit_rows <- get_fit_rows(data, outcome)
  fit_data <- as.data.frame(data[fit_rows])

  # Step 1: Create matched sets via MatchIt.
  # For ATE, MatchIt::matchit() requires method = "full" (full matching) since
  # nearest-neighbour only supports ATT/ATC.  Allow user override via `...`.
  dots <- list(...)
  if (estimand == "ATE" && is.null(dots$method)) {
    check_pkg("optmatch")
    dots$method <- "full"
  }
  # Run MatchIt with the small-sample "Fewer control/treated units"
  # warning demoted to a message: that warning duplicates information
  # already surfaced by `diagnose()`'s match-quality report (matched
  # vs unmatched counts, weight distribution), and emitting it from
  # the fit step makes every analysis on a slightly imbalanced sample
  # look broken. Other MatchIt warnings still propagate to the user.
  m <- withCallingHandlers(
    do.call(
      MatchIt::matchit,
      c(list(ps_formula, data = fit_data, estimand = estimand), dots)
    ),
    warning = function(w) {
      if (grepl("Fewer (control|treated) units", conditionMessage(w))) {
        rlang::inform(conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    }
  )

  # Step 2: Extract matched data with match weights and subclass membership.
  matched_data <- MatchIt::match.data(m)

  # Step 3: Fit the outcome model on the matched sample.
  # Using the treatment-only formula Y ~ A (the matched design handles
  # confounding).  Match weights are applied via the `weights` argument.
  msm_formula <- stats::reformulate(treatment, response = outcome)

  matched_weights <- matched_data$weights
  if (!is.null(weights)) {
    # Align external weights to matched rows. `MatchIt::match.data()`
    # preserves rownames from `fit_data`, and our `fit_data` was
    # constructed via `as.data.frame(data[fit_rows])` which yields
    # rownames "1", ..., sum(fit_rows). So `matched_idx` is a
    # positional index into the length-`sum(fit_rows)` vector
    # `weights[fit_rows]`. Defensively guard against any of these
    # invariants being violated (e.g. a future refactor that uses
    # data.table directly and preserves the original rownames):
    # bail out with a clear error rather than silently producing
    # NA-tainted or misaligned weights.
    matched_idx <- as.integer(rownames(matched_data))
    fit_w <- weights[fit_rows]
    if (
      anyNA(matched_idx) ||
        any(matched_idx < 1L) ||
        any(matched_idx > length(fit_w))
    ) {
      rlang::abort(
        paste0(
          "Cannot align external weights to matched data: matched-row ",
          "indices fall outside `weights[fit_rows]` (length ",
          length(fit_w),
          "). This indicates a row-name invariant violation in `fit_data`."
        ),
        .call = FALSE
      )
    }
    matched_weights <- matched_weights * fit_w[matched_idx]
  }

  model <- stats::glm(
    msm_formula,
    data = matched_data,
    weights = matched_weights,
    family = resolve_family(family)
  )

  new_causatr_fit(
    model = model,
    data = data,
    treatment = treatment,
    outcome = outcome,
    confounders = confounders,
    confounders_tv = NULL,
    family = family,
    method = "matching",
    type = "point",
    estimand = estimand,
    id = NULL,
    time = NULL,
    censoring = NULL,
    history = 1L,
    numerator = NULL,
    weights_obj = NULL,
    match_obj = m,
    call = call,
    details = list(
      fit_rows = fit_rows,
      n_fit = nrow(matched_data),
      n_total = nrow(data),
      n_matched = nrow(matched_data),
      n_unmatched = sum(fit_rows) - nrow(matched_data),
      matched_data = data.table::as.data.table(matched_data),
      weights = weights
    )
  )
}
