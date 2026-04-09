#' Fit a causal survival model
#'
#' @description
#' Convenience wrapper for causal survival analysis using pooled logistic
#' regression as a discrete-time hazard model (Hernán & Robins Ch. 17).
#'
#' ## Algorithm
#'
#' 1. Convert data to person-period format if not already long (using
#'    [to_person_period()]).
#' 2. Fit a pooled logistic regression for the discrete hazard:
#'    \eqn{logit Pr[D_{k+1} = 1 | survived to k, A, L]} =
#'    `time_formula + A + confounders`.
#' 3. For each intervention, predict individual-level hazards, compute
#'    survival as the cumulative product S_i(k) = prod(1 - h_i(m), m <= k),
#'    and average across individuals.
#' 4. Risk difference at time t = \eqn{(1 - S^{a1}(t)) - (1 - S^{a0}(t))}.
#'
#' When the per-interval hazard is small (< 0.1), the pooled logistic model
#' closely approximates a continuous-time Cox model (Technical Point 17.1).
#'
#' @param data A data frame or data.table.  Can be in wide format (one row
#'   per individual with a time-to-event column) or long person-period format
#'   (one row per person per time interval).  If wide, the data is
#'   auto-converted using [to_person_period()].
#' @param outcome Character. Name of the binary event indicator (1 = event
#'   occurred in this interval, 0 = survived / censored).
#' @param treatment Character. Name of the treatment variable.
#' @param confounders A one-sided formula specifying confounders.
#' @param id Character. Name of the individual ID variable.
#' @param time Character. Name of the time variable (interval index).
#' @param censoring Character or `NULL`. Name of the censoring indicator.
#'   If provided, rows where `censoring == 1` are excluded from fitting,
#'   and subsequent rows for that individual are also dropped.
#' @param competing Character or `NULL`. Name of a variable indicating the
#'   type of competing event (for competing risks analysis).
#' @param time_formula A one-sided formula specifying how time enters the
#'   hazard model. Default `~ splines::ns(time, 4)`. Use `~ factor(time)`
#'   for a fully saturated (non-parametric) baseline hazard.
#' @param weights Numeric vector or `NULL`. Pre-computed IPCW or survey
#'   weights.
#' @param ... Additional arguments passed to `glm()`.
#'
#' @return A `causatr_fit` object (with `type = "survival"`) suitable for
#'   use with [contrast()].
#'
#' @references
#' Hernán MA, Robins JM (2025). *Causal Inference: What If*. Chapman &
#' Hall/CRC. Chapter 17.
#'
#' @examples
#' \dontrun{
#' data("nhefs", package = "causatr")
#' fit_surv <- causat_survival(
#'   nhefs,
#'   outcome = "death",
#'   treatment = "qsmk",
#'   confounders = ~ sex + age + race + education +
#'     smokeintensity + smokeyrs + exercise + active + wt71,
#'   id = "seqn",
#'   time = "year"
#' )
#' result <- contrast(fit_surv,
#'   interventions = list(quit = static(1), continue = static(0)),
#'   type = "difference"
#' )
#' }
#'
#' @seealso [causat()], [contrast()], [to_person_period()]
#' @export
causat_survival <- function(
  data,
  outcome,
  treatment,
  confounders,
  id,
  time,
  censoring = NULL,
  competing = NULL,
  time_formula = ~ splines::ns(time, 4),
  weights = NULL,
  ...
) {
  if (!data.table::is.data.table(data)) {
    data <- data.table::as.data.table(data)
  } else {
    data <- data.table::copy(data)
  }

  # Validate required columns exist.
  for (col in c(outcome, treatment, id, time)) {
    if (!col %in% names(data)) {
      rlang::abort(paste0("Column '", col, "' not found in `data`."))
    }
  }

  # Ensure person-period format: data must have one row per id per time point.
  # Check by counting rows per id — if max == 1, the data is still in wide
  # format and needs conversion.
  rows_per_id <- data[, .N, by = c(id)]
  if (max(rows_per_id$N) == 1L) {
    rlang::abort(
      paste0(
        "Data appears to be in wide format (one row per individual). ",
        "Use `to_person_period()` to convert to long person-period format ",
        "before calling `causat_survival()`."
      )
    )
  }

  data.table::setkeyv(data, c(id, time))

  # Build the pooled logistic regression formula:
  # event ~ time_terms + treatment + confounders
  confounder_terms <- attr(stats::terms(confounders), "term.labels")
  time_terms <- attr(stats::terms(time_formula), "term.labels")
  rhs <- c(time_terms, treatment, confounder_terms)
  model_formula <- stats::reformulate(rhs, response = outcome)

  # Identify fitting rows: only rows where the individual has survived
  # up to this point (i.e. no prior event) and is not censored.
  # Create a "previously survived" indicator: for each individual, the
  # event column must be 0 at all prior time points.
  data[,
    prev_event := data.table::shift(
      cumsum(get(outcome)),
      n = 1L,
      fill = 0,
      type = "lag"
    ),
    by = c(id)
  ]

  fit_rows <- data[["prev_event"]] == 0
  if (!is.null(censoring)) {
    fit_rows <- fit_rows & (data[[censoring]] == 0 | is.na(data[[censoring]]))
  }
  fit_data <- data[fit_rows]

  model_weights <- if (!is.null(weights)) weights[fit_rows] else NULL

  # Fit the pooled logistic regression (discrete hazard model).
  model <- stats::glm(
    model_formula,
    data = fit_data,
    family = stats::binomial(),
    weights = model_weights,
    ...
  )

  # Store time points for prediction.
  time_points <- sort(unique(data[[time]]))

  new_causatr_fit(
    model = model,
    data = data,
    treatment = treatment,
    outcome = outcome,
    confounders = confounders,
    confounders_tv = NULL,
    family = "binomial",
    method = "gcomp",
    type = "survival",
    estimand = "ATE",
    id = id,
    time = time,
    censoring = censoring,
    history = 1L,
    numerator = NULL,
    weights_obj = NULL,
    match_obj = NULL,
    call = match.call(),
    details = list(
      fit_rows = fit_rows,
      n_fit = sum(fit_rows),
      n_total = nrow(data),
      time_points = time_points,
      time_formula = time_formula,
      competing = competing,
      model_fn = stats::glm
    )
  )
}
