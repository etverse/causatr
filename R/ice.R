#' Prepare a longitudinal ICE g-computation fit
#'
#' @description
#' Validates the longitudinal data structure, resolves the model families
#' for the final outcome and pseudo-outcome models, and stores all metadata
#' needed by `ice_iterate()` and `compute_contrast_ice()` to run the
#' backward iteration.
#'
#' No models are fitted at this stage -- the sequential models are
#' intervention-dependent and are fitted inside `ice_iterate()` during
#' [contrast()].
#'
#' @param data data.table with lag columns already created by
#'   `prepare_data()`.
#' @param outcome Character. Outcome column name (observed at final time).
#' @param treatment Character scalar or vector. Treatment column name(s).
#' @param confounders One-sided formula of baseline (time-invariant)
#'   confounders.
#' @param confounders_tv One-sided formula of time-varying confounders or
#'   `NULL`.
#' @param family Character or family object for the final outcome model.
#' @param estimand Character. Must be `"ATE"` for longitudinal data.
#' @param history Positive integer or `Inf`. Markov order for lagged
#'   predictors.
#' @param censoring Character or `NULL`. Name of the censoring indicator
#'   column (1 = censored, 0 = uncensored).
#' @param weights Numeric vector or `NULL`. External observation weights
#'   (e.g. IPCW).
#' @param model_fn Function with signature `(formula, data, family,
#'   weights, ...)`. Default `stats::glm`.
#' @param id Character. Name of the individual ID column.
#' @param time Character. Name of the time column.
#' @param call The original `causat()` call (for error messages).
#' @param ... Passed to `model_fn`.
#'
#' @return A `causatr_fit` object with `model = NULL` and all needed
#'   metadata in `details`.
#'
#' @noRd
fit_ice <- function(
  data,
  outcome,
  treatment,
  confounders,
  confounders_tv,
  family,
  estimand,
  history,
  censoring,
  weights,
  model_fn,
  id,
  time,
  call,
  ...
) {
  # Guard against a silent collision with any causatr-reserved column.
  # `ice_iterate()` writes predicted pseudo-outcomes into `.pseudo_y`,
  # and `causat_survival()` uses `.causatr_prev_event` / `_prev_cens`
  # for risk-set bookkeeping. If the user's data already carries any
  # of these, the in-place mutation would silently clobber it. The
  # single source of truth is `CAUSATR_RESERVED_COLS` in `R/utils.R`.
  check_reserved_cols(data)

  # Sorted unique time points define the backward iteration grid. The
  # ICE recursion walks from `time_points[n_times]` down to
  # `time_points[1]`, fitting one outcome (or pseudo-outcome) model per
  # step. Sorting here means downstream code can index by step number
  # without re-sorting and the backward pass is just `seq(n_times, 1)`.
  time_points <- sort(unique(data[[time]]))
  n_times <- length(time_points)

  # Extract baseline confounders as formula term *labels* (not variable
  # names). This preserves transformations like `I(age^2)` and
  # `factor(education)` -- we want the RHS expressions as written in the
  # user's `confounders = ~ ...` formula, not the base variable names.
  baseline_terms <- attr(stats::terms(confounders), "term.labels")

  # Parse effect-modification terms so `ice_build_formula()` can expand
  # `A:modifier` to include `lag1_A:modifier`, `lag2_A:modifier`, etc.
  # at each backward step. Without this, effect modification is only
  # captured at the current period, collapsing heterogeneity from
  # earlier treatment effects.
  em_info <- parse_effect_mod(confounders, treatment)

  # Time-varying confounders, in contrast, are tracked as plain variable
  # names. `ice_build_formula()` uses these to build lag column names
  # (`lag1_L`, `lag2_L`, ...) and to reference the current-time columns;
  # both require the bare identifier.
  tv_vars <- if (!is.null(confounders_tv)) {
    all.vars(confounders_tv)
  } else {
    character(0)
  }

  # Resolve the Markov order. `history = Inf` means "full history" --
  # cap at `n_times - 1` since there are no lags beyond the first time
  # point. An integer `history` lets the user explicitly restrict the
  # lag structure (e.g. `history = 1` for a first-order Markov model).
  max_lag <- if (is.infinite(history)) {
    n_times - 1L
  } else {
    as.integer(history)
  }

  # Resolve the family object. `causat()` accepts any of:
  #   (a) a string  -- e.g. `"gaussian"`  -> get() + call with no args
  #   (b) a closure -- e.g. `gaussian`    -> call with no args
  #   (c) a family  -- e.g. `gaussian()`  -> already resolved
  # Canonicalizing to case (c) here means downstream code can assume
  # `family_obj$family` and `family_obj$link` are both strings.
  family_obj <- if (is.character(family)) {
    get(family, mode = "function")()
  } else if (is.function(family)) {
    family()
  } else {
    family
  }

  # Binary outcome handling needs two different families: `binomial`
  # for the final-time model (where Y is actually 0/1), and
  # `quasibinomial` for every earlier step (where the "outcome" is a
  # predicted probability in (0, 1), not a Bernoulli draw). Fitting a
  # proper binomial on fractional data triggers the "non-integer
  # #successes" warning and is statistically incoherent; quasibinomial
  # uses the same logit link and identical score equations but drops
  # the integer check and estimates dispersion freely. This is the
  # fractional-logistic-regression approach of Papke & Wooldridge (1996),
  # applied to g-computation by Zivich et al. (2024, Section 3.2).
  family_pseudo <- if (family_obj$family == "binomial") {
    stats::quasibinomial(link = family_obj$link)
  } else {
    family_obj
  }

  # Capture user `...` so ICE bootstrap can replay the same model_fn
  # extras (e.g. mgcv::gam `method`, `gamma`) per-step. See B2.
  dots <- list(...)

  new_causatr_fit(
    model = NULL,
    data = data,
    treatment = treatment,
    outcome = outcome,
    confounders = confounders,
    confounders_tv = confounders_tv,
    family = family,
    estimator = "gcomp",
    type = "longitudinal",
    estimand = estimand,
    id = id,
    time = time,
    censoring = censoring,
    history = history,
    numerator = NULL,
    weights_obj = NULL,
    match_obj = NULL,
    call = call,
    details = list(
      time_points = time_points,
      n_times = n_times,
      baseline_terms = baseline_terms,
      tv_vars = tv_vars,
      max_lag = max_lag,
      em_info = em_info,
      model_fn = model_fn,
      family_outcome = family_obj,
      family_pseudo = family_pseudo,
      weights = weights,
      dots = dots
    )
  )
}


#' Build the outcome model formula for ICE at a given time step
#'
#' @description
#' Constructs the RHS formula for the ICE model at a given time index.
#' At time index k (0-based), available lags are `min(k, max_lag)`.
#' Treatment and time-varying confounder columns that are entirely `NA` at
#' the current step are dropped (handles the case where a time-varying
#' confounder is not measured at the earliest times).
#'
#' Baseline confounders are always included (they are time-invariant and
#' should never be `NA`).
#'
#' When effect-modification terms are present (e.g. `A:sex` in the
#' confounders formula), the function auto-expands them across available
#' treatment lags: at time_idx = 2 with max_lag = 2, `A:sex` produces
#' `A:sex`, `lag1_A:sex`, `lag2_A:sex`. This ensures the ICE outcome
#' model captures heterogeneous treatment effects from all prior periods,
#' not just the current one. Without this expansion, roughly half the
#' intended heterogeneity is collapsed.
#'
#' @param response Character. LHS variable name (`outcome` at final step,
#'   `".pseudo_y"` at earlier steps).
#' @param treatment Character scalar or vector. Treatment column name(s).
#' @param baseline_terms Character vector. RHS formula terms for baseline
#'   confounders (from `attr(terms(), "term.labels")`).
#' @param tv_vars Character vector. Plain column names of time-varying
#'   confounders.
#' @param time_idx Integer, 0-based. The current time index in the
#'   backward iteration (0 = first time, K = last time).
#' @param max_lag Integer. Maximum lag order (from `history`).
#' @param data_at_time data.table. Rows at the current time step,
#'   used to check for all-`NA` columns.
#' @param em_info A `causatr_em_info` object from `parse_effect_mod()`,
#'   or `NULL` if no EM terms are present.
#'
#' @return A formula object.
#'
#' @noRd
ice_build_formula <- function(
  response,
  treatment,
  baseline_terms,
  tv_vars,
  time_idx,
  max_lag,
  data_at_time,
  em_info = NULL
) {
  # The number of available lags at time index `k` is `min(k, max_lag)`:
  # at t = 0 there are zero lags regardless of max_lag; at t = 1 there
  # is one lag; and so on, capped by the user's Markov-order choice.
  available_lags <- min(time_idx, max_lag)
  lag_vars <- c(treatment, tv_vars)

  # RHS starts with the *current* time values of treatment and TV
  # confounders. These are always included (they are the "now" variables
  # the ICE model conditions on at this step).
  rhs_dynamic <- c(treatment, tv_vars)

  # Append lag columns in order: lag1_A, lag1_L, lag2_A, lag2_L, ...
  # The `prepare_data()` step materializes these as actual columns, so
  # we only need to reference them by name here.
  if (available_lags > 0L) {
    for (lag_k in seq_len(available_lags)) {
      for (v in lag_vars) {
        rhs_dynamic <- c(rhs_dynamic, paste0("lag", lag_k, "_", v))
      }
    }
  }

  # Guard against including columns that don't exist or are entirely
  # NA at this time step. Two cases this protects against:
  #   (a) A time-varying confounder first measured at t = 1 (not t = 0)
  #       -- at t = 0 the column is all NA, so we drop it from the formula.
  #   (b) `lag1_A` at t = 0 -- no past exists, column is all NA by
  #       construction, same drop.
  # Without this guard the glm would fail with a cryptic "contrasts can
  # be applied only to factors with 2 or more levels" or similar.
  valid <- vapply(
    rhs_dynamic,
    function(col) {
      col %in% names(data_at_time) && !all(is.na(data_at_time[[col]]))
    },
    logical(1)
  )
  # Inform (once per fit / per column) when a user-supplied column is
  # dropped. `lag*_` columns are expected to drop at t = 0, so we
  # filter them from the notification -- they are auto-generated by
  # prepare_data() and their drop is structural, not a user concern.
  # The expected (structural) drops we suppress:
  #   - lag columns at t = 0 (no prior period)
  #   - the current-time treatment when it's literally never measured
  #     (covered by the rhs_dynamic test below -- these are rare)
  dropped <- rhs_dynamic[!valid]
  user_dropped <- dropped[!grepl("^lag[0-9]+_", dropped)]
  if (length(user_dropped) > 0L) {
    rlang::inform(
      paste0(
        "ICE step (time_idx = ",
        time_idx,
        "): dropping all-NA column(s) from the outcome model formula: ",
        paste0("'", user_dropped, "'", collapse = ", "),
        ". Check your `confounders_tv` spec if this is unexpected."
      ),
      .frequency = "regularly",
      .frequency_id = paste0(
        "causatr_ice_dropped_",
        paste(user_dropped, collapse = "_")
      )
    )
  }
  rhs_dynamic <- rhs_dynamic[valid]

  # Baseline confounders are time-invariant and always available, so
  # they go in unconditionally. `baseline_terms` comes in as term labels
  # (e.g. `I(age^2)`) rather than variable names, which means
  # `reformulate()` reuses the user's original transformations verbatim.
  #
  # When EM terms are present (e.g. `A:sex`), auto-expand them to
  # include lag versions (`lag1_A:sex`, `lag2_A:sex`, ...) so the
  # outcome model captures heterogeneous treatment effects from all
  # available prior periods. The expansion terms go through the same
  # all-NA column validity check as the dynamic columns above.
  em_lag_terms <- character(0L)
  if (!is.null(em_info) && em_info$has_em) {
    for (em_term in em_info$em_terms) {
      lag_terms <- expand_em_lag_terms(em_term, available_lags)
      em_lag_terms <- c(em_lag_terms, lag_terms)
    }
    # Validity check: drop lag interaction terms that reference columns
    # which are all-NA at this time step (same guard as for rhs_dynamic).
    # This naturally drops `lag1_A:sex` at t = 0 where `lag1_A` doesn't
    # exist or is all-NA.
    if (length(em_lag_terms) > 0L) {
      em_lag_valid <- vapply(
        em_lag_terms,
        function(term) {
          # Extract the variable names from the term (e.g. "lag1_A:sex"
          # -> c("lag1_A", "sex")) and check all exist and are non-NA.
          term_vars <- all.vars(parse(text = term)[[1L]])
          all(vapply(
            term_vars,
            function(v) {
              v %in% names(data_at_time) && !all(is.na(data_at_time[[v]]))
            },
            logical(1L)
          ))
        },
        logical(1L)
      )
      em_lag_terms <- em_lag_terms[em_lag_valid]
    }
  }

  rhs_terms <- c(rhs_dynamic, baseline_terms, em_lag_terms)

  stats::reformulate(rhs_terms, response = response)
}


#' Apply intervention and update lag columns for longitudinal data
#'
#' @description
#' Applies the intervention to the **current-time** treatment column via
#' `apply_intervention()` and leaves all lag columns untouched. The
#' Robins iterated conditional expectation algorithm (and
#' `lmtp::lmtp_tmle()`'s sequential-regression shortcut) substitutes
#' only the current treatment at each backward step, so treatment-lag
#' columns (`lag1_A`, `lag2_A`, ...) must continue to hold the
#' **observed** A_{k-1}, A_{k-2}, ... in the prediction frame at step k.
#' Recomputing lag columns from the intervened treatment would
#' double-count any non-static intervention and return (K+1)*delta
#' instead of K*delta for `shift(delta)` on K time points (verified
#' against `lmtp::lmtp_tmle` in `test-simulation.R`).
#'
#' Time-varying confounder lags are **not** modified either -- ICE
#' conditions on the observed covariate history.
#'
#' @param data data.table. Full person-period data with lag columns.
#' @param treatment Character. Treatment column name(s).
#' @param intervention A `causatr_intervention` object or `NULL`.
#' @param id_col Character. Name of the individual ID column.
#' @param time_col Character. Name of the time column.
#'
#' @return A modified **copy** of `data` (original is never mutated).
#'
#' @noRd
ice_apply_intervention_long <- function(
  data,
  treatment,
  intervention,
  id_col,
  time_col
) {
  # The Robins iterated conditional expectation algorithm
  # (and `lmtp::lmtp_tmle()`'s sequential-regression shortcut)
  # substitutes only the **current-time** treatment at each
  # backward step. Lag columns must therefore hold the OBSERVED
  # A_{k-1}, A_{k-2}, ... in the prediction frame at step k, so
  # the chain of pseudo-outcome regressions back-substitutes
  # through the conditional-expectation tower correctly. We
  # apply the intervention to the current treatment column only
  # and leave both treatment-lag and covariate-lag columns
  # untouched. (Recomputing lag1_A from the intervened treatment
  # would double-count any non-static intervention via the
  # lag-column path AND the current-A prediction path, returning
  # (K+1)*delta instead of K*delta for `shift(delta)` on K time
  # points -- see test-simulation.R for the lmtp-validated truth
  # check.)
  data_iv <- apply_intervention(data, treatment, intervention)
  data_iv <- data.table::copy(data_iv)
  data.table::setkeyv(data_iv, c(id_col, time_col))
  data_iv
}


#' Run ICE backward iteration for a single intervention
#'
#' @description
#' Implements the ICE g-computation algorithm (Zivich et al. 2024):
#'
#' 1. Fit outcome model at the **final** time point:
#'    `E[Y | A_K, lags, L_K, baseline]` among uncensored with observed Y.
#' 2. Predict under intervention at final time -> pseudo-outcomes.
#' 3. At each **earlier** time step (backward):
#'    - Use pseudo-outcomes from the next step as the response.
#'    - Fit `E[pseudo | A_k, lags, L_k, baseline]` among uncensored with
#'      valid pseudo-outcomes (uses `quasibinomial` for binary outcomes).
#'    - Predict under intervention -> new pseudo-outcomes.
#' 4. The pseudo-outcomes at the **first** time step are the individual-level
#'    counterfactual expectations \eqn{\hat{Y}^*_0}.
#'
#' @param fit A `causatr_fit` of type `"longitudinal"` (from
#'   `fit_ice()`).
#' @param intervention A `causatr_intervention` object (or `NULL` for the
#'   natural course).
#'
#' @return A list with:
#'   \describe{
#'     \item{`pseudo_final`}{Numeric vector of \eqn{\hat{Y}^*_0}, one per individual
#'       at the first time point (ordered as they appear in the data).}
#'     \item{`models`}{Named list of fitted model objects, one per time
#'       point (keyed by the time-point value as character).}
#'     \item{`data_iv`}{data.table -- the intervention-modified data.}
#'     \item{`fit_ids`}{Named list of character vectors -- individual IDs
#'       in each model's fitting set (needed for sandwich variance).}
#'   }
#'
#' @noRd
ice_iterate <- function(fit, intervention) {
  data <- fit$data
  details <- fit$details
  outcome <- fit$outcome
  treatment <- fit$treatment
  id_col <- fit$id
  time_col <- fit$time
  censoring <- fit$censoring

  time_points <- details$time_points
  n_times <- details$n_times
  baseline_terms <- details$baseline_terms
  tv_vars <- details$tv_vars
  max_lag <- details$max_lag
  model_fn <- details$model_fn
  family_outcome <- details$family_outcome
  family_pseudo <- details$family_pseudo
  external_weights <- details$weights
  # User's stashed `...` for `model_fn`. The per-step `replay_fit()`
  # calls below take care of duplicate-key stripping in one place
  # (R/utils.R).
  model_fn_dots <- details$dots

  # Build the intervention-modified person-period frame once. `data_iv`
  # has the counterfactual treatment column and matching lags, but the
  # original covariate trajectory -- this is what we condition on below.
  data_iv <- ice_apply_intervention_long(
    data,
    treatment,
    intervention,
    id_col,
    time_col
  )

  # `pseudo` is the rolling vector of individual-level pseudo-outcomes.
  # Initialized to NA for every unique id; at the final step it gets
  # filled with predictions under the fitted outcome model, and each
  # backward step overwrites it with predictions from the pseudo-outcome
  # model at that step. After the full loop, `pseudo[first_ids]` is
  # \eqn{\hat Y^*_0} -- the individual counterfactual expectation at baseline.
  all_ids <- unique(data[[id_col]])
  id_chr <- as.character(all_ids)
  pseudo <- stats::setNames(rep(NA_real_, length(all_ids)), id_chr)

  # Storage indexed by time step. `fit_ids[[k]]` is the set of
  # individuals who contributed to fitting the k-th model; the sandwich
  # variance engine uses these to align per-step score residuals back
  # to the length-n IF vector.
  models <- vector("list", n_times)
  names(models) <- as.character(time_points)
  fit_ids <- vector("list", n_times)
  names(fit_ids) <- as.character(time_points)

  # -- Step 1: fit the outcome model at the FINAL time point.
  # This is the only step that uses the real observed outcome Y; every
  # earlier step uses a pseudo-outcome constructed by prior predictions.

  final_time <- time_points[n_times]
  final_idx <- n_times - 1L # 0-based time index for formula construction

  # Fitting mask at the final time: must be at the last time point AND
  # uncensored AND have an observed (non-NA) outcome. Each of these
  # conditions excludes a real failure mode -- censored individuals
  # never reach the final step; individuals with missing Y can't
  # contribute to a likelihood fit.
  mask_final <- data[[time_col]] == final_time
  uncens <- is_uncensored(data, censoring)
  fit_mask <- mask_final & uncens & !is.na(data[[outcome]])

  fit_data <- data[fit_mask]
  fit_ids[[n_times]] <- as.character(fit_data[[id_col]])

  em_info <- details$em_info

  # Formula construction is shared between all time steps; at the
  # final step the response is the real outcome, at earlier steps
  # it's `.pseudo_y`.
  formula_k <- ice_build_formula(
    outcome,
    treatment,
    baseline_terms,
    tv_vars,
    final_idx,
    max_lag,
    fit_data,
    em_info
  )

  # Build args for `model_fn` (default `stats::glm`) via do.call. We
  # go through do.call rather than a direct call because `glm()` uses
  # non-standard evaluation for `weights =`: passing `weights = NULL`
  # explicitly is NOT the same as omitting it, and do.call lets us
  # include the argument conditionally.
  model_args <- list(
    formula = formula_k,
    data = fit_data,
    family = family_outcome
  )
  if (!is.null(external_weights)) {
    model_args$weights <- external_weights[fit_mask]
  }
  models[[n_times]] <- replay_fit(model_fn, model_args, model_fn_dots)

  # Predict under the intervention for ALL uncensored individuals at
  # the final time (not just those in the fitting set). This is the
  # g-formula pattern: fit on a subset with observed outcomes, predict
  # for the whole target population under the counterfactual.
  pred_mask <- mask_final & uncens
  pred_data <- data_iv[pred_mask]
  preds <- stats::predict(
    models[[n_times]],
    newdata = pred_data,
    type = "response"
  )
  pred_ids <- as.character(data[pred_mask][[id_col]])
  pseudo[pred_ids] <- preds

  # -- Steps 2+: backward iteration (time K-1 down to time 0).
  # At each step we fit a "pseudo-outcome" model regressing the
  # already-filled `pseudo[i]` on the current time's treatment,
  # covariates, and lags, then overwrite `pseudo[i]` with predictions
  # under the intervention. This is the ICE algorithm -- the
  # conditional-mean integration collapses into a chain of regressions.

  for (step_i in seq(n_times - 1L, 1L, by = -1L)) {
    current_time <- time_points[step_i]
    time_idx <- step_i - 1L # 0-based

    mask_current <- data[[time_col]] == current_time
    mask_uncens <- mask_current & uncens

    # Pull pseudo-outcomes for uncensored individuals at this time.
    # Some may still be NA if they were censored later and never
    # received a prediction at the previous step -- `has_pseudo` filters
    # those out.
    current_ids <- as.character(data[mask_uncens][[id_col]])
    pseudo_y <- pseudo[current_ids]

    has_pseudo <- !is.na(pseudo_y)
    if (sum(has_pseudo) == 0L) {
      # A degenerate case: no valid pseudo-outcomes at this step. This
      # shouldn't happen on well-formed longitudinal data but abort
      # with a useful message if it does.
      rlang::abort(
        paste0(
          "No valid pseudo-outcomes at time ",
          current_time,
          " for ICE backward iteration."
        ),
        .call = FALSE
      )
    }

    # Materialize a copy of the fitting rows and attach the pseudo
    # response as `.pseudo_y`. `data.table::copy` + in-place `:=` is
    # necessary: without `copy` the mutation would leak into `data`.
    fit_data <- data.table::copy(data[mask_uncens][has_pseudo])
    fit_data[, .pseudo_y := pseudo_y[has_pseudo]]
    fit_ids[[step_i]] <- as.character(fit_data[[id_col]])

    # Build and fit the pseudo-outcome model. The `family_pseudo`
    # resolution in `fit_ice()` swapped binomial -> quasibinomial so
    # this call is happy with fractional responses.
    formula_k <- ice_build_formula(
      ".pseudo_y",
      treatment,
      baseline_terms,
      tv_vars,
      time_idx,
      max_lag,
      fit_data,
      em_info
    )

    # `do.call(model_fn, model_args_k)` rather than a direct call so
    # the `weights =` argument can be conditionally included: glm()
    # uses NSE on `weights`, and `weights = NULL` is not the same as
    # omitting it. External weights (IPCW, survey weights) MUST flow
    # through every backward pseudo-outcome model, otherwise the
    # stacked-EE sandwich for the longitudinal target is biased.
    # `mask_uncens` is length n_total; we subset it to `has_pseudo`
    # so the weight vector aligns with `fit_data` row-for-row.
    model_args_k <- list(
      formula = formula_k,
      data = fit_data,
      family = family_pseudo
    )
    if (!is.null(external_weights)) {
      model_args_k$weights <- external_weights[mask_uncens][has_pseudo]
    }
    models[[step_i]] <- replay_fit(model_fn, model_args_k, model_fn_dots)

    # Predict under intervention for ALL individuals at the current
    # time point (not just the fitting subset). This keeps the
    # population covariate trajectory intact as the recursion walks
    # backward toward baseline.
    pred_all <- data_iv[mask_current]
    preds <- stats::predict(
      models[[step_i]],
      newdata = pred_all,
      type = "response"
    )
    pred_ids_all <- as.character(data[mask_current][[id_col]])
    pseudo[pred_ids_all] <- preds
  }

  # After the loop, `pseudo` at the first time point holds the
  # individual-level counterfactual expectations \hat Y^*_{0,i}. The
  # marginal mean \hat\mu is just mean(pseudo_final) over the target
  # population -- computed in compute_contrast_ice(), not here.
  first_time <- time_points[1]
  rows_first <- data[[time_col]] == first_time
  first_ids <- as.character(data[rows_first][[id_col]])
  pseudo_final <- unname(pseudo[first_ids])

  list(
    pseudo_final = pseudo_final,
    models = models,
    data_iv = data_iv,
    fit_ids = fit_ids
  )
}
