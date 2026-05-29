#' Build the per-period propensity formula for longitudinal IPW
#'
#' @description
#' Counterpart to `ice_build_formula()` but for the **treatment density**
#' rather than the outcome. The k-th period's formula is
#' \deqn{A \sim \text{baseline} + L_k + \mathrm{lag}_1 A + \mathrm{lag}_1 L
#'       + \ldots + \mathrm{lag}_m A + \mathrm{lag}_m L,}
#' where `m = available_lags = min(time_idx, max_lag)`. Lag columns
#' that are entirely NA at the current period (which is what happens
#' at the first time step) are silently dropped.
#'
#' Effect-modification expansion is **not** applied here -- EM is rejected
#' upstream in `fit_longitudinal_ipw()`.
#'
#' @param treatment Character scalar. Treatment column name.
#' @param baseline_terms Character vector of baseline confounder term
#'   labels (already EM-stripped).
#' @param tv_vars Character vector of time-varying confounder column
#'   names.
#' @param available_lags Integer. Number of lag periods to include
#'   (`min(time_idx, max_lag)`).
#' @param data_at_time data.table of rows at the current period, used
#'   to drop all-NA lag columns.
#'
#' @return A two-sided formula `A ~ ...` ready for
#'   `fit_treatment_model()`.
#'
#' @noRd
build_longitudinal_ps_formula <- function(
  treatment,
  baseline_terms,
  tv_vars,
  available_lags,
  data_at_time
) {
  # Current-time TV confounders enter every period's formula. (The
  # treatment itself is the response, not a predictor.)
  rhs_dynamic <- tv_vars

  # Lag columns: `lag1_A`, `lag1_L`, ..., `lagm_A`, `lagm_L`.
  # The naming convention `lag{j}_{col}` matches what `create_lag_vars()`
  # in `prepare_data()` materialises; any mismatch here silently drops
  # the column in the all-NA guard below rather than causing an error.
  rhs_dynamic <- c(
    rhs_dynamic,
    build_lag_terms(c(treatment, tv_vars), available_lags)
  )

  # Drop columns that don't exist or are all-NA at this time. Mirrors
  # `ice_build_formula()`'s defensive guard. Most common case: at
  # `time_points[1]`, all `lag1_*` columns are NA by construction.
  valid <- vapply(
    rhs_dynamic,
    function(col) {
      col %in% names(data_at_time) && !all(is.na(data_at_time[[col]]))
    },
    logical(1)
  )
  rhs_dynamic <- rhs_dynamic[valid]

  # Baseline terms that overlap with TV vars are already represented
  # by the current-time TV entry; including them twice produces a
  # cosmetically wrong formula (R silently deduplicates, but the
  # formula object misleads anyone who inspects it).
  rhs_terms <- c(rhs_dynamic, setdiff(baseline_terms, rhs_dynamic))
  if (length(rhs_terms) == 0L) {
    # Intercept-only propensity at the first period when neither
    # baseline confounders nor TV confounders are supplied. Rare in
    # practice but the formula still needs an RHS.
    rhs_terms <- "1"
  }

  # For multivariate treatment, use the first component as a placeholder
  # response — callers immediately strip via remove_response() before
  # passing to fit_treatment_models().
  resp <- if (length(treatment) > 1L) treatment[1L] else treatment
  stats::reformulate(rhs_terms, response = resp)
}


#' Build the per-period numerator propensity formula for stabilized longitudinal IPW
#'
#' @description
#' Counterpart to `build_longitudinal_ps_formula()` for the per-period
#' numerator model `g_k(A_k | A_{1..k-1}, V)` under
#' `stabilize = "marginal"`. The chain rule keeps treatment lags in
#' the conditioning set; only the time-varying confounders are
#' dropped (the whole point of stabilization is to remove the
#' multiplicative L-dependence across periods that inflates the
#' cumulative product weight).
#'
#' Two paths:
#'
#' - `numerator = NULL` (default): condition on treatment lags only.
#'   The k-th formula is `A ~ lag1_A + ... + lag_m_A`.
#' - `numerator = ~ V`: user supplies an explicit conditioning set
#'   (typically baseline covariates). The k-th formula is
#'   `A ~ lag1_A + ... + lag_m_A + V`.
#'
#' Treatment lags are kept by default because dropping them would
#' break the chain-rule factorisation of the joint numerator density
#' (Robins, Hernán, Brumback 2000) -- the stabilizer must still
#' integrate to 1 over the joint support of the treatment trajectory.
#'
#' @param treatment Character scalar.
#' @param numerator One-sided formula or `NULL`.
#' @param available_lags Integer.
#' @param data_at_time data.table for the current period (used to drop
#'   all-NA lag columns at `time_points[1]`).
#' @return Two-sided formula `A ~ ...`.
#' @noRd
build_longitudinal_numerator_ps_formula <- function(
  treatment,
  numerator,
  available_lags,
  data_at_time
) {
  rhs_dynamic <- character(0L)
  if (available_lags > 0L) {
    for (lag_k in seq_len(available_lags)) {
      rhs_dynamic <- c(rhs_dynamic, paste0("lag", lag_k, "_", treatment))
    }
  }

  valid <- vapply(
    rhs_dynamic,
    function(col) {
      col %in% names(data_at_time) && !all(is.na(data_at_time[[col]]))
    },
    logical(1)
  )
  rhs_dynamic <- rhs_dynamic[valid]

  user_terms <- if (is.null(numerator)) {
    character(0L)
  } else {
    attr(stats::terms(numerator), "term.labels")
  }

  rhs_terms <- c(rhs_dynamic, user_terms)
  if (length(rhs_terms) == 0L) {
    rhs_terms <- "1"
  }

  resp <- if (length(treatment) > 1L) treatment[1L] else treatment
  stats::reformulate(rhs_terms, response = resp)
}


#' Strip the response from a two-sided formula
#'
#' @description
#' `fit_treatment_model()` takes a one-sided `confounders` formula, but
#' `build_longitudinal_ps_formula()` produces a two-sided `A ~ ...`
#' formula (so it can serve as the propensity formula directly). Strip
#' the LHS so we can reuse `fit_treatment_model()` without rewriting
#' it.
#'
#' @param ps_formula Two-sided formula.
#' @return One-sided formula with the same RHS.
#' @noRd
remove_response <- function(ps_formula) {
  stats::reformulate(
    attr(stats::terms(ps_formula), "term.labels"),
    response = NULL
  )
}
