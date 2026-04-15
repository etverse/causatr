#' Fit a causal survival model
#'
#' @description
#' Convenience wrapper for causal survival analysis using pooled logistic
#' regression as a discrete-time hazard model (Hernán & Robins Ch. 17).
#'
#' **Status: scaffolded.** The pooled logistic fit itself is implemented,
#' but the survival-curve contrast step in [contrast()] is not yet wired
#' up and aborts with an informative error. Competing-risks analysis (the
#' `competing` argument) is also not yet implemented and aborts at fit
#' time when supplied. Treat this function as experimental until both
#' paths land.
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
#'   type of competing event (for competing risks analysis). **Not yet
#'   implemented**: supplying a non-`NULL` value currently aborts. Reserved
#'   for a future release.
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

  # Reject user data that already carries a causatr-internal column
  # name (e.g. `.causatr_prev_event`). See R12 in the 2026-04-15 review.
  check_reserved_cols(data)

  # Same up-front weights validation as `causat()` — reject NA,
  # non-finite, negative, or mis-sized weight vectors so users see a
  # specific error instead of a cryptic GLM abort.
  check_weights(weights, nrow(data))

  # `competing` is reserved for Phase 6 (sub-distribution hazard /
  # Aalen-Johansen). The pooled-logistic path here does NOT apply
  # any competing-risks adjustment, so accepting a non-NULL value
  # would return a cause-deleted hazard model masquerading as a
  # competing-risks estimator. Abort until the implementation lands.
  if (!is.null(competing)) {
    rlang::abort(
      paste0(
        "Competing-risks survival analysis is not yet implemented ",
        "(planned for Phase 6). The `competing` argument is reserved ",
        "but currently has no effect; re-running with `competing = NULL` ",
        "would silently fit a plain cause-deleted hazard model, which ",
        "is biased in the presence of competing events. Track progress ",
        "in PHASE_6_SURVIVAL.md."
      ),
      .call = FALSE
    )
  }

  # Phase 6 is scaffolded: the pooled-logistic fit below runs to
  # completion, but `contrast()` on the resulting fit will abort with a
  # scaffold message. Warn at fit time as well so users learn about the
  # limitation before investing in a long fit rather than after.
  rlang::inform(
    c(
      "`causat_survival()` is scaffolded (Phase 6).",
      i = "The pooled-logistic hazard model will fit, but `contrast()` on survival fits is not yet implemented.",
      i = "Track progress in PHASE_6_SURVIVAL.md."
    ),
    .frequency = "regularly",
    .frequency_id = "causat_survival_scaffold"
  )

  # Sort once up front so every subsequent `by = c(id)` aggregation sees
  # a deterministic within-id row order. Previously the two aggregations
  # below ran on unsorted data and relied on data.table's first-seen
  # group ordering matching between the two passes. Any upstream reorder
  # would desync them. See R9 in the 2026-04-15 critical review.
  data.table::setkeyv(data, c(id, time))

  # Ensure person-period format: every id must appear at more than one
  # time point. Previous implementations used `max(rows_per_id$N) == 1L`
  # which only caught uniform wide format — a mixed frame with some
  # single-row ids and some multi-row ids would slip through and end up
  # fitting a degenerate risk set at the survival step. We now require
  # that every id has at least 2 rows (a minimally well-formed
  # person-period panel), and additionally that `time` actually varies
  # within each id so duplicated rows at the same time don't masquerade
  # as a long panel.
  rows_per_id <- data[, .N, by = c(id)]
  if (any(rows_per_id$N == 1L)) {
    n_wide <- sum(rows_per_id$N == 1L)
    rlang::abort(
      paste0(
        "Data does not appear to be in person-period format: ",
        n_wide,
        " of ",
        nrow(rows_per_id),
        " unique `",
        id,
        "` values have only a single row. ",
        "Use `to_person_period()` to convert to long format before ",
        "calling `causat_survival()`."
      )
    )
  }
  # `.SDcols = time` + `.SD[[1L]]` avoids a brittle `get(time)` call
  # in the j expression (data.table's NSE env binds column symbols
  # *and* the outer scope, and `get(time)` there can resolve to the
  # special-name binding instead of the user's variable).
  time_unique_per_id <- data[,
    list(nt = length(unique(.SD[[1L]]))),
    by = c(id),
    .SDcols = time
  ]
  if (any(time_unique_per_id$nt < rows_per_id$N)) {
    rlang::abort(
      paste0(
        "Some individuals have duplicated rows at the same `",
        time,
        "` value. ",
        "`causat_survival()` requires unique (",
        id,
        ", ",
        time,
        ") pairs."
      )
    )
  }

  # Build the pooled logistic regression formula for discrete-time
  # survival analysis (Hernán & Robins Ch. 17):
  #   event ~ time_terms + treatment + confounders
  # where `time_terms` is a flexible function of time (default
  # `~ ns(time, 4)`). Pooled logistic = run one logistic regression
  # over ALL (person, time) rows with a time-specific intercept
  # — equivalent to a discrete-time hazard model.
  confounder_terms <- attr(stats::terms(confounders), "term.labels")
  time_terms <- attr(stats::terms(time_formula), "term.labels")
  rhs <- c(time_terms, treatment, confounder_terms)
  model_formula <- stats::reformulate(rhs, response = outcome)

  # Fitting rows are the (person, time) cells where the individual
  # is still at risk: no prior event and not censored. The "at-risk"
  # convention is standard for hazard-based discrete-time models:
  # once an individual has the event, they drop out of the risk set
  # for all subsequent periods.
  #
  # `prev_event` is computed via within-id cumsum of the event
  # column, shifted by 1 lag so that the CURRENT row's indicator
  # reflects whether a prior event has occurred (not the current
  # one). Fill with 0 for the first period — everyone starts event-free.
  # Data is already keyed (id, time) from the up-front setkeyv, so
  # `shift` sees a deterministic within-id ordering.
  data[,
    .causatr_prev_event := data.table::shift(
      cumsum(get(outcome)),
      n = 1L,
      fill = 0,
      type = "lag"
    ),
    by = c(id)
  ]

  # Censoring rule: drop rows WHERE `censoring == 1` AND all subsequent
  # rows for that individual. The docstring promised this behavior but
  # the previous implementation only excluded the current row — a
  # subject censored at t=2 still contributed risk at t=3. The cumsum
  # + lag trick mirrors `prev_event` above: 0 at and before first
  # censor, >=1 after. Rows with `prev_cens > 0` OR `censoring == 1`
  # are both excluded. See B5 in the 2026-04-15 critical review.
  if (!is.null(censoring)) {
    data[,
      .causatr_prev_cens := data.table::shift(
        cumsum(get(censoring)),
        n = 1L,
        fill = 0,
        type = "lag"
      ),
      by = c(id)
    ]
    fit_rows <- data[[".causatr_prev_event"]] == 0 &
      data[[".causatr_prev_cens"]] == 0 &
      is_uncensored(data, censoring)
  } else {
    fit_rows <- data[[".causatr_prev_event"]] == 0
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
    estimator = "gcomp",
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
