#' Fit a longitudinal IPW model (per-period treatment density chain)
#'
#' @description
#' Longitudinal companion to `fit_ipw()` for time-varying treatments.
#' Fits the sequential factorisation
#' \deqn{f(A_0, A_1, \ldots, A_K \mid L) =
#'       \prod_{k=0}^{K} f(A_k \mid \bar A_{k-1}, \bar L_k)}
#' as `K + 1` independent per-period density models, one per
#' `time_points[k]`. Each per-period model is fit on the subset of
#' rows at `time == time_points[k]` with a formula that conditions on
#' baseline confounders, current-period time-varying confounders, and
#' lagged treatment / lagged time-varying confounders up to `history`.
#' The lag columns are already materialised by `prepare_data()` /
#' `create_lag_vars()`, so the per-period formula just references them
#' by name.
#'
#' Per-period density models are fit independently via
#' `fit_treatment_model()`, which makes the bread of the stacked
#' propensity system block-diagonal -- the variance engine sums `K`
#' single-model `apply_model_correction()` calls (same shape as the
#' multivariate IPW propensity correction).
#'
#' ## What is rejected at this fit boundary
#'
#' Longitudinal IPW currently supports:
#'
#' - univariate treatment (binary / continuous / count via
#'   `propensity_family`)
#' - `static()` / `shift()` / `scale_by()` / `dynamic()` interventions
#' - intercept-only Hajek MSM (no effect modification)
#' - unstabilized weights only (`stabilize = "none"`)
#'
#' Multivariate treatment, effect-modification (`A:modifier` in
#' `confounders`), `stabilize = "marginal"`, and `ipsi()` are rejected
#' upfront with classed errors. `threshold()` continues to be rejected
#' per period by `check_intervention_family_compat()`.
#'
#' @param data data.table from `prepare_data()` with lag columns and
#'   sorted by `(id, time)`.
#' @param outcome Character. Outcome column name (observed at the
#'   final time point only; earlier periods may carry NA).
#' @param treatment Character scalar. Treatment column name.
#' @param confounders One-sided formula of baseline (time-invariant)
#'   confounders.
#' @param confounders_tv One-sided formula of time-varying confounders
#'   or `NULL`.
#' @param family Character or family object for the final-period MSM.
#' @param estimand Character. Must be `"ATE"` for longitudinal data
#'   (`check_estimand_trt_compat()` already enforces this).
#' @param history Positive integer or `Inf`. Markov order for lag
#'   inclusion in the per-period propensity formula.
#' @param numerator One-sided formula or `NULL`. Reserved for chunk
#'   10b (stabilized weights); rejected here when non-`NULL`.
#' @param weights Numeric vector or `NULL`. External observation
#'   weights (e.g. survey weights, IPCW), aligned to `data` rows.
#' @param model_fn Function. Outcome MSM fitter; default `stats::glm`.
#' @param propensity_model_fn Function or `NULL`. Per-period
#'   propensity fitter. `NULL` reuses `model_fn`.
#' @param propensity_family Character or `NULL`. Explicit propensity
#'   family override (`"poisson"` / `"negbin"`); broadcast across all
#'   periods.
#' @param id Character. Individual ID column.
#' @param time Character. Time column.
#' @param call The original `causat()` call.
#' @param ... Forwarded to the per-period propensity fitter.
#'
#' @return A `causatr_fit` of `estimator = "ipw"`, `type =
#'   "longitudinal"` carrying:
#'   \describe{
#'     \item{`details$treatment_models_by_time`}{Length-`K+1` named
#'       list of `causatr_treatment_model` objects, named by string
#'       representation of `time_points[k]`.}
#'     \item{`details$time_points`}{Sorted unique time values.}
#'     \item{`details$n_times`}{Number of distinct time points.}
#'     \item{`details$max_lag`}{Resolved Markov order.}
#'     \item{`details$baseline_terms` / `tv_vars`}{Formula bookkeeping
#'       used to rebuild per-period propensity formulas at variance /
#'       bootstrap time.}
#'   }
#'
#' @noRd
fit_longitudinal_ipw <- function(
  data,
  outcome,
  treatment,
  confounders,
  confounders_tv,
  family,
  estimand,
  history,
  numerator,
  weights,
  model_fn,
  propensity_model_fn,
  propensity_family = NULL,
  stabilize = "none",
  id,
  time,
  call,
  ...
) {
  # Each abort is a classed error so callers can detect on `class`
  # rather than parsing error text.
  if (length(treatment) > 1L) {
    rlang::abort(
      c(
        "Multivariate longitudinal IPW is not supported.",
        i = "Use a single-treatment longitudinal IPW or `estimator = 'gcomp'` for joint interventions."
      ),
      class = "causatr_longitudinal_multivariate_pending"
    )
  }

  # `numerator =` is the user-facing knob for a custom numerator
  # formula. Under `stabilize = "none"` it has no meaning -- the
  # numerator is the same density as the denominator -- so reject the
  # combination upfront.
  if (!is.null(numerator)) {
    if (stabilize == "none") {
      rlang::abort(
        c(
          "`numerator =` requires `stabilize = 'marginal'`.",
          i = "Drop `numerator` or set `stabilize = 'marginal'`."
        ),
        class = "causatr_longitudinal_numerator_without_stabilize"
      )
    }
    check_formula(numerator)
  }

  # Effect modification under longitudinal IPW reuses the same parser
  # the point IPW path uses (`parse_effect_mod()`). EM terms appear in
  # `confounders` as `A:modifier`; the per-period propensity formulas
  # strip every treatment-touching term via `em_info$confounder_terms`
  # (so the propensity is `A ~ baseline + tv + lags` with no `A:`
  # interactions), and the final-period MSM expands from `Y ~ 1` to
  # `Y ~ 1 + modifier` via `build_ipw_msm_formula()`.
  #
  # Known limitation (Robins 2000): the modifier MUST be a
  # **baseline** covariate. A time-varying modifier in an MSM
  # conditions on a post-treatment variable, biasing the estimand
  # via mediator + collider paths. We don't enforce baseline-ness at
  # runtime because time-varying status isn't inferable from the data;
  # the constraint is doc-level.
  em_info <- check_confounders_treatment(
    confounders,
    treatment,
    estimator = "ipw"
  )

  # Sorted unique time points. Per-period propensity models index into
  # `time_points[k]`; lag columns (`lag1_A`, ...) created by
  # `prepare_data()` are aligned to this ordering because
  # `create_lag_vars()` keys on `(id, time)` before shifting.
  # Sorting here guarantees `treatment_models_by_time` is keyed in
  # chronological order even when the data arrives unsorted.
  time_points <- sort(unique(data[[time]]))
  n_times <- length(time_points)

  if (n_times < 2L) {
    rlang::abort(
      c(
        paste0(
          "`type = 'longitudinal'` requires at least 2 distinct time points (got ",
          n_times,
          ")."
        ),
        i = "Use `type = 'point'` if the data really has only one period per individual."
      ),
      class = "causatr_longitudinal_too_few_times"
    )
  }

  # Resolve the Markov order. Mirrors `fit_ice()`'s convention.
  max_lag <- if (is.infinite(history)) {
    n_times - 1L
  } else {
    as.integer(history)
  }

  # Baseline confounders enter every per-period propensity. Strip the
  # treatment-touching (EM) terms for parity with point IPW
  # `build_ps_formula()` -- here EM is rejected upfront, but using
  # `em_info$confounder_terms` keeps the call self-consistent if the
  # rejection is later relaxed.
  baseline_terms <- em_info$confounder_terms
  if (length(baseline_terms) == 0L) {
    baseline_terms <- character(0L)
  }

  # Time-varying confounder names (plain identifiers). The per-period
  # formula references the current-time `L` and the lag columns
  # `lag1_L`, `lag2_L`, ... by name.
  tv_vars <- if (!is.null(confounders_tv)) {
    all.vars(confounders_tv)
  } else {
    character(0)
  }

  # Capture user `...` so the bootstrap can replay the propensity
  # fitter configuration on every replicate via `refit_ipw()`.
  dots <- list(...)

  # Resolve the propensity fitter. For univariate longitudinal IPW the
  # default fitter is `model_fn` (typically `stats::glm`); a user can
  # override with `propensity_model_fn` (e.g. `mgcv::gam`). NB / count
  # families opt in via `propensity_family`.
  trt_family_first <- detect_treatment_family(data[[treatment]])
  prop_model_fn <- if (!is.null(propensity_model_fn)) {
    propensity_model_fn
  } else if (trt_family_first == "categorical") {
    nnet::multinom
  } else if (identical(propensity_family, "negbin")) {
    check_pkg("MASS")
    MASS::glm.nb
  } else {
    model_fn
  }

  # Per-period propensity fits. At each `time_points[k]`:
  #   1. Subset rows where `time == time_points[k]` (one row per id).
  #   2. Build the conditioning formula with available lags (`min(k-1,
  #      max_lag)`).
  #   3. Fit `f(A_k | A_{1..k-1}, L_{1..k}, baseline)` via
  #      `fit_treatment_model()`.
  #
  # The k-th model's `fit_rows` slot is local to its period subset --
  # so `compute_density_ratio_weights()` aligns its weight vector to
  # the period subset, not the full person-period data. Downstream we
  # carry the per-period subsets explicitly so the variance engine can
  # rebuild propensity-side design matrices.
  #
  # `treatment_models_by_time` is named by the string representation
  # of each time point (e.g. "0", "1", "2") so callers can retrieve
  # the k-th period model by name without positional indexing.
  treatment_models_by_time <- vector("list", n_times)
  fit_data_by_time <- vector("list", n_times)
  per_period_formula <- vector("list", n_times)

  for (k in seq_len(n_times)) {
    rows_k <- data[[time]] == time_points[k]
    data_k <- data[rows_k]

    # Available lag depth at time_idx (k - 1) following the ICE
    # convention -- `time_idx = 0` at the first period (no lags),
    # `time_idx = K` at the last (capped at max_lag).
    available_lags <- min(k - 1L, max_lag)

    ps_formula <- build_longitudinal_ps_formula(
      treatment = treatment,
      baseline_terms = baseline_terms,
      tv_vars = tv_vars,
      available_lags = available_lags,
      data_at_time = data_k
    )

    # `fit_treatment_model()` handles its own NA dropping and
    # `propensity_family` dispatch. The returned object exposes
    # `fit_rows` relative to `data_k`, `alpha_hat`, `X_prop`, etc.
    # Lag columns used here (e.g. `lag1_A`, `lag1_L`) are pre-computed
    # by `prepare_data()` / `create_lag_vars()`; referencing them by
    # name in the formula is safe because `build_longitudinal_ps_formula()`
    # drops any columns that are all-NA at this period.
    tm_args <- list(
      data = data_k,
      treatment = treatment,
      confounders = remove_response(ps_formula),
      model_fn = prop_model_fn,
      propensity_family = propensity_family
    )
    if (!is.null(weights)) {
      tm_args$weights <- weights[rows_k]
    }
    tm_k <- do.call(
      fit_treatment_model,
      c(tm_args, dots)
    )

    treatment_models_by_time[[k]] <- tm_k
    fit_data_by_time[[k]] <- data_k
    per_period_formula[[k]] <- ps_formula
  }
  names(treatment_models_by_time) <- as.character(time_points)
  names(fit_data_by_time) <- as.character(time_points)
  names(per_period_formula) <- as.character(time_points)

  # Stabilized weights: per-period numerator models g_k(A_k | A_{k-1},
  # ..., A_0) with L stripped from the conditioning. The full
  # density-ratio weight per period becomes
  #   w_k = g_k(d(A_k, L_k) | A_{1..k-1}^obs) * |Jac d^{-1}|
  #         / f_k(A_k | A_{1..k-1}^obs, L_{1..k}^obs).
  # Numerator parameters gamma are held fixed under the variance
  # engine's `numDeriv` perturbation of the denominator alpha (same
  # nuisance-fixed convention as multivariate IPW and
  # sigma / theta for Gaussian / negbin densities). Bootstrap refits
  # both gamma and alpha, capturing full variance.
  #
  # Custom `numerator` formula overrides the default "drop L" rule:
  # `numerator = ~ baseline` keeps only baseline confounders;
  # `numerator = ~ 1` reduces to a marginal-by-period density.
  numerator_models_by_time <- NULL
  if (stabilize == "marginal") {
    numerator_models_by_time <- vector("list", n_times)
    for (k in seq_len(n_times)) {
      data_k <- fit_data_by_time[[k]]
      available_lags <- min(k - 1L, max_lag)

      # Per-period numerator formula. Default behaviour (no `numerator`
      # argument): drop L (baseline + tv + lag_L); keep treatment lags
      # only. With a user-supplied `numerator` formula we honour it
      # verbatim as the conditioning RHS, plus the available treatment
      # lags so the per-period factorisation g_k(A_k | A_{1..k-1}, V)
      # still respects the chain rule. `V` is the user-supplied set.
      num_ps_formula <- build_longitudinal_numerator_ps_formula(
        treatment = treatment,
        numerator = numerator,
        available_lags = available_lags,
        data_at_time = data_k
      )

      num_args <- list(
        data = data_k,
        treatment = treatment,
        confounders = remove_response(num_ps_formula),
        model_fn = prop_model_fn,
        propensity_family = propensity_family
      )
      if (!is.null(weights)) {
        num_args$weights <- weights[data[[time]] == time_points[k]]
      }
      numerator_models_by_time[[k]] <- do.call(
        fit_treatment_model,
        c(num_args, dots)
      )
    }
    names(numerator_models_by_time) <- as.character(time_points)
  }

  # Final-period rows define the MSM scope and the marginal-mean
  # average. The Hajek intercept of `Y ~ 1` weighted by the cumulative
  # product weight reads `E[Y^d]` directly off the (sole) coefficient.
  # Only the final period carries a non-missing outcome in a typical
  # longitudinal study; earlier periods have `NA` in the outcome column.
  last_time <- time_points[n_times]
  rows_final <- data[[time]] == last_time
  fit_rows_final <- rows_final & !is.na(data[[outcome]])

  if (sum(fit_rows_final) == 0L) {
    rlang::abort(
      c(
        paste0(
          "No final-period (`time == ",
          last_time,
          "`) rows have non-missing outcome -- the longitudinal IPW MSM has nothing to fit."
        ),
        i = "Check that the outcome column is populated at the last time point."
      ),
      class = "causatr_longitudinal_no_final_outcome"
    )
  }

  # Placeholder MSM model. Same role as `fit_ipw()`'s placeholder
  # `Y ~ A` -- gives `print()` / `summary()` something to display
  # before any contrast is requested. The real per-intervention MSM
  # is refit inside `compute_ipw_contrast_longitudinal()`.
  # `Y ~ 1` (unweighted) rather than `Y ~ A` because in longitudinal
  # data the treatment varies by period, so `Y ~ A` would be ambiguous.
  fam_obj <- resolve_family(family)
  placeholder_args <- list(
    formula = stats::reformulate("1", response = outcome),
    data = data[fit_rows_final],
    family = fam_obj
  )
  if (!is.null(weights)) {
    placeholder_args$weights <- weights[fit_rows_final]
  }
  placeholder_model <- do.call(model_fn, placeholder_args)

  new_causatr_fit(
    model = placeholder_model,
    data = data,
    treatment = treatment,
    outcome = outcome,
    confounders = confounders,
    confounders_tv = confounders_tv,
    family = family,
    estimator = "ipw",
    type = "longitudinal",
    estimand = estimand,
    id = id,
    time = time,
    censoring = NULL,
    history = history,
    numerator = numerator,
    weights_obj = NULL,
    match_obj = NULL,
    call = call,
    details = list(
      time_points = time_points,
      n_times = n_times,
      max_lag = max_lag,
      baseline_terms = baseline_terms,
      tv_vars = tv_vars,
      em_info = em_info,
      treatment_models_by_time = treatment_models_by_time,
      numerator_models_by_time = numerator_models_by_time,
      fit_data_by_time = fit_data_by_time,
      per_period_formula = per_period_formula,
      fit_rows_final = fit_rows_final,
      n_fit_final = sum(fit_rows_final),
      n_total = nrow(data),
      weights = weights,
      dots = dots,
      model_fn = model_fn,
      propensity_model_fn = prop_model_fn,
      propensity_family = propensity_family,
      stabilize = stabilize,
      is_multivariate = FALSE
    )
  )
}


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

  rhs_terms <- c(rhs_dynamic, baseline_terms)
  if (length(rhs_terms) == 0L) {
    # Intercept-only propensity at the first period when neither
    # baseline confounders nor TV confounders are supplied. Rare in
    # practice but the formula still needs an RHS.
    rhs_terms <- "1"
  }

  stats::reformulate(rhs_terms, response = treatment)
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

  stats::reformulate(rhs_terms, response = treatment)
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


#' Compute per-intervention longitudinal IPW bundles
#'
#' @description
#' Longitudinal companion to `compute_ipw_contrast_point()`. For each
#' intervention:
#'
#' 1. Apply the intervention at every period (the same intervention
#'    is broadcast across `time_points[1:K]`).
#' 2. Build the cumulative density-ratio weight per individual via
#'    `compute_longitudinal_weights()`.
#' 3. Multiply external (survey / IPCW) weights from
#'    `fit$details$weights` evaluated on the final-period rows.
#' 4. Refit `Y ~ 1` weighted on the final-period rows; read
#'    \eqn{\hat\mu_a} off the intercept (Hajek).
#'
#' Returns the per-intervention bundles + the scalars
#' `compute_contrast()` needs to dispatch variance / bootstrap.
#'
#' @param fit A `causatr_fit` from `fit_longitudinal_ipw()`.
#' @param interventions Named list of `causatr_intervention` objects
#'   (or `NULL` for natural course).
#' @param target_within_first Logical vector flagging the target
#'   population at the first time point (mirrors ICE's convention).
#'
#' @return A list with `bundles`, `mu_hat`, `id_first` (character ids
#'   in target order), and `n_total` (`nrow(fit$data)`).
#'
#' @noRd
compute_ipw_contrast_longitudinal <- function(
  fit,
  interventions,
  target_within_first
) {
  data <- fit$data
  treatment <- fit$treatment
  outcome <- fit$outcome
  details <- fit$details
  treatment_models_by_time <- details$treatment_models_by_time
  numerator_models_by_time <- details$numerator_models_by_time
  fit_data_by_time <- details$fit_data_by_time
  time_points <- details$time_points
  n_times <- details$n_times
  id_col <- fit$id
  time_col <- fit$time
  fit_rows_final <- details$fit_rows_final
  fam_obj <- resolve_family(fit$family)
  model_fn <- details$model_fn
  ext_w <- details$weights

  # First-time-point rows define the per-individual weight ordering.
  # We average / aggregate per id over the final-period rows; the
  # natural id order comes from sorting the first-period subset.
  first_time <- time_points[1]
  rows_first <- data[[time_col]] == first_time
  ids_first <- as.character(data[rows_first][[id_col]])
  n_id <- length(ids_first)

  # Final-period rows are where we evaluate the MSM. Each id has one
  # final-period row; the cumulative weight per id is multiplied with
  # the external weight on that row before fitting the weighted Hajek
  # intercept.
  last_time <- time_points[n_times]
  data_final <- data[fit_rows_final]
  ids_final <- as.character(data_final[[id_col]])
  n_final <- nrow(data_final)

  # Map id -> position in the first-time ordering. Used to project the
  # cumulative weight (length n_id, in first-time order) onto the
  # final-period rows. Reverse lookup: given an id from `ids_final`,
  # `id_to_first_idx[id]` gives its row in `w_id` from
  # `compute_longitudinal_weights()`.
  id_to_first_idx <- stats::setNames(seq_len(n_id), ids_first)

  ext_w_final <- if (is.null(ext_w)) NULL else ext_w[fit_rows_final]

  # MSM formula. Default `Y ~ 1` (intercept-only Hájek). With
  # baseline-modifier effect modification (`A:modifier` in
  # `confounders`), expands to `Y ~ 1 + modifier` via the existing
  # `build_ipw_msm_formula()` so `predict()` returns
  # stratum-specific counterfactual means.
  em_info <- details$em_info
  msm_formula <- build_ipw_msm_formula(outcome, em_info)

  int_names <- names(interventions)

  # `ipsi()` is not supported under longitudinal IPW. Surface a clear
  # error rather than silently producing an untested result.
  for (nm in int_names) {
    iv <- interventions[[nm]]
    if (
      !is.null(iv) && inherits(iv, "causatr_intervention") && iv$type == "ipsi"
    ) {
      rlang::abort(
        c(
          "`ipsi()` is not supported under longitudinal IPW.",
          i = paste0("Intervention `", nm, "` is `ipsi()`."),
          i = "Use `static()`, `shift()`, `scale_by()`, or `dynamic()` instead, or switch to point IPW."
        ),
        class = "causatr_longitudinal_ipsi_pending"
      )
    }
  }

  bundles <- lapply(int_names, function(nm) {
    iv <- interventions[[nm]]

    # Per-id cumulative weight under the intervention. Internally the
    # routine loops over `time_points` and multiplies per-period
    # density ratios, returning a length-n_id vector ordered by
    # `ids_first`. NULL intervention -> all-ones (natural course).
    w_id <- compute_longitudinal_weights(
      treatment_models_by_time = treatment_models_by_time,
      fit_data_by_time = fit_data_by_time,
      ids_first = ids_first,
      id_col = id_col,
      intervention = iv,
      numerator_models_by_time = numerator_models_by_time
    )

    # Project per-id weights onto the final-period row order.
    w_final <- w_id[id_to_first_idx[ids_final]]

    # Compose with external weights. Same multiplicative convention as
    # `compute_ipw_contrast_point()`: density-ratio weights and survey
    # / IPCW weights enter as independent reweightings of the target
    # population.
    w_combined <- if (is.null(ext_w_final)) w_final else w_final * ext_w_final

    # Final-period weighted MSM. Intercept-only Hajek under
    # `Y ~ 1`; with EM the MSM is `Y ~ 1 + modifier` so `predict()`
    # returns stratum-specific counterfactual means. The treatment is
    # absorbed by the cumulative density-ratio weight.
    msm_args <- list(
      formula = msm_formula,
      data = data_final,
      family = msm_family(fam_obj),
      weights = w_combined
    )
    msm_model <- do.call(model_fn, msm_args)

    # Counterfactual marginal mean over the target subset of
    # individuals. `target_within_first` is length n_id; subset to the
    # final-period rows that correspond to target ids. Under EM the
    # MSM contains modifier columns, so `predict()` correctly returns
    # the per-row stratum mean -- aggregation by modifier value is
    # handled by the user via `by =` in `contrast()`.
    target_ids_final <- target_within_first[id_to_first_idx[ids_final]]
    preds_final <- stats::predict(msm_model, type = "response")
    valid_final <- target_ids_final & !is.na(preds_final)
    w_target_final <- if (!is.null(ext_w_final)) {
      ext_w_final[valid_final]
    } else {
      NULL
    }
    mu_hat_iv <- maybe_weighted_mean(preds_final[valid_final], w_target_final)

    list(
      intervention = iv,
      msm_model = msm_model,
      weights_id = w_id,
      weights_final = w_combined,
      preds_final = preds_final,
      valid_final = valid_final,
      target_ids_final = target_ids_final,
      mu_hat = mu_hat_iv
    )
  })
  names(bundles) <- int_names

  mu_hat <- vapply(bundles, function(b) b$mu_hat, numeric(1))
  names(mu_hat) <- int_names

  list(
    bundles = bundles,
    mu_hat = mu_hat,
    ids_first = ids_first,
    id_to_first_idx = id_to_first_idx,
    ids_final = ids_final,
    n_id = n_id,
    n_final = n_final,
    n_total = nrow(data)
  )
}
