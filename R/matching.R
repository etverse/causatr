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

  # MatchIt requires a binary treatment (`MatchIt:::check_treat`
  # aborts for 3+ levels). Intercept this upstream with a clearer
  # error so users don't hit the MatchIt internal message, and so
  # the rejection path is explicit in causatr's own namespace.
  trt_vals <- data[[treatment]]
  trt_levels <- if (is.factor(trt_vals)) {
    levels(droplevels(trt_vals))
  } else {
    unique(stats::na.omit(trt_vals))
  }
  if (length(trt_levels) > 2L) {
    rlang::abort(
      paste0(
        "Matching supports only binary treatments, but `",
        treatment,
        "` has ",
        length(trt_levels),
        " levels. Use `estimator = \"gcomp\"` or `estimator = \"ipw\"` ",
        "for categorical treatments."
      ),
      .call = FALSE
    )
  }

  # Parse effect-modification terms and reject bare treatment in
  # confounders (`~ L + A`). True EM terms (`A:sex`) are detected
  # and stored for downstream MSM expansion. The propensity formula
  # strips EM terms automatically via `build_ps_formula()`.
  em_info <- check_confounders_treatment(
    confounders,
    treatment,
    estimator = "matching"
  )

  # Matching effect modification is wired via MSM expansion from
  # `Y ~ A` to `Y ~ A + modifier + A:modifier`. For now, reject EM
  # terms until the MSM expansion is implemented.
  if (em_info$has_em) {
    em_labels <- vapply(em_info$em_terms, function(x) x$term, character(1L))
    rlang::abort(
      c(
        paste0(
          "Effect modification via `A:modifier` is not yet supported ",
          "for `estimator = \"matching\"`."
        ),
        x = paste0(
          "Offending term(s): ",
          paste(em_labels, collapse = ", "),
          "."
        ),
        i = paste0(
          "Use `estimator = \"gcomp\"` for heterogeneous treatment effects, ",
          "or `by = \"modifier\"` in `contrast()` for stratum-specific ",
          "summaries of a homogeneous effect."
        )
      ),
      class = "causatr_em_unsupported",
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

  matched_weights <- combine_match_and_external_weights(
    matched_data,
    external_weights = weights,
    fit_rows = fit_rows
  )

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
    estimator = "matching",
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
      weights = weights,
      # Capture MatchIt arguments (post ATE->full defaulting) so
      # `refit_matching()` replays the same matching specification
      # verbatim at bootstrap time. See B2 in the 2026-04-15 review:
      # dropping these produced bootstrap SEs that corresponded to
      # default nearest-neighbor matching even when the user asked
      # for caliper, ratio, distance, etc.
      dots = dots
    )
  )
}


#' Combine MatchIt match weights with optional external weights
#'
#' @description
#' `MatchIt::match.data()` returns a data frame whose rownames are the
#' positional indices of `fit_data` (which, in turn, was built via
#' `as.data.frame(data[fit_rows])`). To attach an external weight vector
#' (e.g. survey weights) to the matched sample we therefore need to:
#' (1) subset `external_weights` to the fit rows, and
#' (2) reindex by `as.integer(rownames(matched_data))` so each matched row
#'     picks up its original observation's weight.
#'
#' Used by both `fit_matching()` (at fit time) and `refit_matching()` (at
#' bootstrap time) so the row-name invariant is asserted in exactly one
#' place; a future change to `as.data.frame.data.table`'s rowname
#' behavior surfaces as a single loud abort rather than two silently
#' drifting call sites.
#'
#' @param matched_data Data frame returned by `MatchIt::match.data()`. Must
#'   carry a `weights` column (match weights) and rownames that parse as
#'   integer positional indices into `fit_rows`.
#' @param external_weights Numeric vector of length `nrow(data)` or `NULL`.
#'   The external (survey / IPCW) weights passed to `causat()`.
#' @param fit_rows Logical vector of length `nrow(data)` flagging which rows
#'   were passed to `MatchIt::matchit()`.
#'
#' @return Numeric vector of combined weights, length `nrow(matched_data)`.
#'
#' @noRd
combine_match_and_external_weights <- function(
  matched_data,
  external_weights,
  fit_rows
) {
  matched_weights <- matched_data$weights
  if (is.null(external_weights)) {
    return(matched_weights)
  }
  matched_idx <- as.integer(rownames(matched_data))
  fit_w <- external_weights[fit_rows]
  # Range + NA check: the invariant is that rownames(matched_data) are
  # positional indices in 1..sum(fit_rows). Any violation here means the
  # rowname assumption has been broken upstream (e.g. by a change in
  # `as.data.frame.data.table`), in which case the weight alignment
  # would silently corrupt the IF / bootstrap estimates.
  if (
    anyNA(matched_idx) ||
      any(matched_idx < 1L) ||
      any(matched_idx > length(fit_w))
  ) {
    rlang::abort(
      paste0(
        "Cannot align external weights to matched data: matched-row ",
        "indices fall outside `external_weights[fit_rows]` (length ",
        length(fit_w),
        "). This indicates a row-name invariant violation in `fit_data`."
      ),
      .call = FALSE
    )
  }
  matched_weights * fit_w[matched_idx]
}
