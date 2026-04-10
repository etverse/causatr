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
  confounder_terms <- attr(stats::terms(confounders), "term.labels")
  ps_formula <- stats::reformulate(confounder_terms, response = treatment)

  # Fit rows: exclude missing outcomes before matching.
  fit_rows <- !is.na(data[[outcome]])
  fit_data <- as.data.frame(data[fit_rows])

  # Step 1: Create matched sets via MatchIt.
  # For ATE, MatchIt::matchit() requires method = "full" (full matching) since
  # nearest-neighbour only supports ATT/ATC.  Allow user override via `...`.
  dots <- list(...)
  if (estimand == "ATE" && is.null(dots$method)) {
    check_pkg("optmatch")
    dots$method <- "full"
  }
  m <- do.call(
    MatchIt::matchit,
    c(list(ps_formula, data = fit_data, estimand = estimand), dots)
  )

  # Step 2: Extract matched data with match weights and subclass membership.
  matched_data <- MatchIt::match.data(m)

  # Step 3: Fit the outcome model on the matched sample.
  # Using the treatment-only formula Y ~ A (the matched design handles
  # confounding).  Match weights are applied via the `weights` argument.
  msm_formula <- stats::reformulate(treatment, response = outcome)

  matched_weights <- matched_data$weights
  if (!is.null(weights)) {
    # Multiply external weights with match weights.
    # Need to align: matched_data rows come from fit_data[matched indices].
    matched_weights <- matched_weights *
      weights[fit_rows][as.integer(rownames(matched_data))]
  }

  msm_family <- if (is.character(family)) {
    get(family, mode = "function", envir = asNamespace("stats"))()
  } else if (is.function(family)) {
    family()
  } else {
    family
  }

  model <- stats::glm(
    msm_formula,
    data = matched_data,
    weights = matched_weights,
    family = msm_family
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
      matched_data = data.table::as.data.table(matched_data)
    )
  )
}
