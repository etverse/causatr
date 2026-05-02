#' Diagnostics for a fitted causal model
#'
#' @description
#' Computes diagnostics appropriate to the causal estimator and
#' treatment type:
#' - **Binary IPW**: propensity-score positivity (tail violations),
#'   covariate balance (SMDs via `cobalt`), weight distribution
#'   (treated / control / overall with ESS).
#' - **Continuous IPW**: density-range positivity (low-density tail
#'   observations), covariate balance (correlations), weight
#'   distribution (overall ESS).
#' - **Categorical IPW**: per-level probability positivity, covariate
#'   balance (pairwise SMDs via `cobalt`), weight distribution.
#' - **Count IPW** (Poisson / NB): density-range positivity,
#'   covariate balance, weight distribution.
#' - **Multivariate IPW**: per-component positivity, covariate
#'   balance on the first component, combined product-weight
#'   distribution.
#' - **Longitudinal IPW**: per-period positivity and weight
#'   distribution from `treatment_models_by_time`, per-period
#'   covariate balance.
#' - **Longitudinal gcomp (ICE)**: per-period covariate balance
#'   (no weights or positivity since ICE has no treatment model).
#' - `"matching"`: covariate balance before and after matching (via
#'   `cobalt`), match quality summary.
#' - `"gcomp"`: unadjusted covariate imbalance between treatment
#'   groups.
#'
#' **Per-intervention dispatch.** When `interventions =` is supplied,
#' each intervention spawns its own diagnostic panel: the positivity
#' summary is shared across panels (it depends only on the fitted
#' treatment density model, not the intervention), but the weight
#' distribution is computed per intervention. The default call
#' `diagnose(fit)` produces a single `obs` panel.
#'
#' @param fit A `causatr_fit` object returned by [causat()].
#' @param interventions Named list of `causatr_intervention` objects (or
#'   `NULL` entries for natural course). When `NULL` (the default), a single
#'   panel keyed `obs` is produced using the standard observed-treatment
#'   diagnostic for the estimator. When non-`NULL`, each named entry spawns
#'   its own per-intervention panel.
#' @param by Character scalar naming a baseline modifier for stratified
#'   balance reporting. Currently rejected with
#'   `causatr_diag_em_pending` -- the hook is reserved so future chunks can
#'   lift the rejection without a signature break.
#' @param stats Character vector. Balance statistics to compute. Passed to
#'   `cobalt::bal.tab()`. For binary treatments, valid options include `"m"`
#'   (standardised mean differences), `"v"` (variance ratios), and `"ks"`
#'   (Kolmogorov-Smirnov). Default `c("m", "v")`.
#' @param thresholds Named numeric vector. Balance thresholds for flagging
#'   imbalance, e.g. `c(m = 0.1, v = 2)`. Default `c(m = 0.1)`.
#' @param ps_bounds Numeric vector of length 2. Lower and upper bounds for
#'   flagging positivity violations. Default `c(0.025, 0.975)`.
#'
#' @return A `causatr_diag` object with slots:
#'   \describe{
#'     \item{`per_intervention`}{Named list of per-intervention panels. Each
#'       panel is itself a list with `positivity`, `balance`, `weights`
#'       slots (any of which may be `NULL`).}
#'     \item{`interventions`}{Character vector of intervention keys, in the
#'       order they appear in `per_intervention`.}
#'     \item{`positivity`, `balance`, `weights`}{Top-level shortcuts that
#'       point to the corresponding slots of the first panel; preserved for
#'       backward compatibility with the flat shape used pre-chunk 11a.}
#'     \item{`match_quality`}{data.table or `NULL`: match quality summary
#'       (matching only). Lives at the top level because matching is done
#'       once at fit time and is intervention-agnostic.}
#'     \item{`estimator`}{Character: the causal estimator.}
#'     \item{`fit_info`}{Named list with summary metadata
#'       (`treatment_type`, `estimand`, `type`, `has_em`).}
#'     \item{`fit`}{The original `causatr_fit` (stored for `plot()`).}
#'   }
#'
#' @details
#' ## Positivity
#' For binary treatment, fits a logistic regression of the treatment on the
#' confounders and flags individuals whose estimated propensity score falls
#' outside `ps_bounds`. The returned `positivity` table summarises the
#' propensity score distribution and the number/fraction of near-violations.
#' The propensity score is intervention-independent, so each panel carries
#' an identical positivity table.
#'
#' ## Balance (IPW and matching)
#' If the `cobalt` package is installed, balance is computed via
#' `cobalt::bal.tab()` on the propensity formula or matchit object. This
#' provides standardised mean differences (SMD), variance ratios, and KS
#' statistics. If `cobalt` is not installed, a simpler data.table-based SMD
#' comparison is returned. The chunk 11a balance view is the unadjusted
#' SMD across treatment groups -- post-weighting balance under specific
#' interventions / estimands lands in a later chunk.
#'
#' ## Weight distribution (IPW only)
#' For the default `obs` panel: summarises the observed-treatment
#' Horvitz-Thompson weights (`1/p` on treated, `1/(1-p)` on controls), which
#' is the standard Hernán & Robins Ch. 12 IPW weight diagnostic. For each
#' user-supplied intervention: summarises the per-intervention density-ratio
#' weight from `compute_density_ratio_weights()`. Both views report
#' mean / SD / min / max plus the effective sample size (Kish 1965)
#' `ESS = (sum w)^2 / sum(w^2)` for treated / control / overall groups.
#'
#' ## Match quality (matching only)
#' Reports the number matched, number discarded, and the fraction of the
#' original sample retained.
#'
#' @references
#' Greifer N (2024). cobalt: Covariate Balance Tables and Plots.
#' \url{https://ngreifer.github.io/cobalt/}
#'
#' Hernán MA, Robins JM (2020). *Causal Inference: What If*. CRC.
#'
#' @examples
#' \dontrun{
#' data("nhefs", package = "causatr")
#' fit <- causat(nhefs, outcome = "wt82_71", treatment = "qsmk",
#'               confounders = ~ sex + age + wt71,
#'               estimator = "ipw")
#' # Default: single observed-treatment panel.
#' diag <- diagnose(fit)
#' print(diag)
#' plot(diag)
#'
#' # Per-intervention dispatch:
#' diag2 <- diagnose(fit, interventions = list(a1 = static(1), a0 = static(0)))
#' print(diag2)
#' }
#'
#' @seealso [causat()], [plot.causatr_diag()]
#' @export
diagnose <- function(
  fit,
  interventions = NULL,
  by = NULL,
  stats = c("m", "v"),
  thresholds = c(m = 0.1),
  ps_bounds = c(0.025, 0.975)
) {
  if (!inherits(fit, "causatr_fit")) {
    rlang::abort("`fit` must be a `causatr_fit` object returned by `causat()`.")
  }

  check_diag_pending_gates(fit, interventions, by, call = rlang::current_env())

  # Resolve the panel set. `default_view = TRUE` flags the auto-injected
  # `obs` panel so the IPW handler knows to use the observed-treatment
  # HT view rather than the natural-course density-ratio view (which
  # collapses to all-1 weights under `compute_density_ratio_weights()`).
  default_view <- is.null(interventions)
  if (default_view) {
    interventions <- list(obs = NULL)
  } else {
    if (is.null(names(interventions)) || any(!nzchar(names(interventions)))) {
      rlang::abort(
        c(
          "Every entry in `interventions` must be named.",
          i = paste0(
            "Use e.g. `interventions = list(a1 = static(1), a0 = static(0))`."
          )
        ),
        .call = FALSE
      )
    }
  }

  # Longitudinal fits dispatch to a specialised path that builds
  # per-period positivity / balance / weight sub-tables.
  if (identical(fit$type, "longitudinal")) {
    return(diagnose_longitudinal(
      fit,
      interventions,
      default_view,
      stats,
      thresholds,
      ps_bounds
    ))
  }

  # Cache the propensity-score positivity table once: it depends only
  # on the fitted treatment density model, not on the intervention, so
  # recomputing it per panel would be wasteful. Same logic applies to
  # the unadjusted balance view (chunk 11a). Both are shared across
  # every panel below.
  positivity_shared <- compute_positivity(fit, ps_bounds)
  balance_shared <- compute_balance(fit, stats, thresholds)

  per_intervention <- lapply(seq_along(interventions), function(i) {
    nm <- names(interventions)[i]
    iv <- interventions[[i]]
    is_default_obs <- default_view && i == 1L
    weights_panel <- compute_panel_weights(fit, iv, is_default_obs)
    list(
      positivity = positivity_shared,
      balance = balance_shared,
      weights = weights_panel
    )
  })
  names(per_intervention) <- names(interventions)

  fit_info <- list(
    treatment_type = detect_diag_treatment_type(fit),
    estimand = fit$estimand,
    type = fit$type,
    has_em = !is.null(fit$details$em_info) &&
      length(fit$details$em_info$em_terms) > 0L
  )

  new_causatr_diag(
    per_intervention = per_intervention,
    match_quality = compute_match_quality(fit),
    estimator = fit$estimator,
    fit_info = fit_info,
    fit = fit
  )
}

#' Reject combinations not yet supported by `diagnose()`
#'
#' @description
#' Centralises the up-front rejection gates so each unsupported
#' combination aborts with a stable `causatr_diag_*_pending` classed
#' error. Two gates remain: IPW + ATT/ATC estimand (chunk 11d) and
#' `by =` effect-modification stratification (chunk 11d).
#'
#' @param fit A `causatr_fit` object.
#' @param interventions Named list (or `NULL`).
#' @param by Character scalar or `NULL`.
#' @param call Caller environment for error attribution (so the
#'   user-visible message names `diagnose()` rather than the helper).
#' @return `NULL` invisibly; aborts with the appropriate classed error.
#' @noRd
check_diag_pending_gates <- function(fit, interventions, by, call) {
  if (!is.null(by)) {
    rlang::abort(
      c(
        "`by =` stratified diagnostics are not yet implemented.",
        i = paste0(
          "Effect-modification-aware balance lands in a later chunk; ",
          "the `by =` argument is reserved so future chunks can lift this ",
          "rejection without a signature break."
        )
      ),
      class = "causatr_diag_em_pending",
      call = call
    )
  }

  if (fit$estimator != "ipw") {
    # gcomp + matching under ATT / ATC are unchanged from ATE in
    # chunk 11a: gcomp's balance is the unadjusted SMD (estimand-
    # agnostic), and matching's diagnostics flow from the matchit
    # object, which already reflects the requested estimand. Only
    # IPW needs the within-treated / within-control balance rewrite,
    # so the estimand gate fires only there.
    return(invisible(NULL))
  }

  if (fit$estimand %in% c("ATT", "ATC")) {
    rlang::abort(
      c(
        paste0(
          "`diagnose()` for IPW under `estimand = '",
          fit$estimand,
          "'` is not yet implemented."
        ),
        i = paste0(
          "Within-treated / within-control balance and the corresponding ",
          "weight summary land in a later chunk; the chunk 11a IPW balance ",
          "view assumes ATE."
        )
      ),
      class = "causatr_diag_estimand_pending",
      call = call
    )
  }

  invisible(NULL)
}

#' Build the per-intervention weight panel for a single intervention key
#'
#' @description
#' Three dispatch axes:
#'  1. Estimator: only `"ipw"` produces a non-NULL panel under chunk 11a;
#'     gcomp and matching panels return NULL (matching exposes its own
#'     metrics via the top-level `match_quality` slot, gcomp has no
#'     weight-distribution concept).
#'  2. Default-obs flag: when `is_default_obs = TRUE`, the panel uses the
#'     observed-treatment Horvitz-Thompson view (`1/p` and `1/(1-p)`),
#'     preserving the pre-chunk-11a default behaviour. Otherwise the panel
#'     uses `compute_density_ratio_weights()` for the per-intervention
#'     density-ratio weight.
#'  3. NULL intervention with `is_default_obs = FALSE`: a literal natural-
#'     course panel where every weight is 1. Useful when the user wants to
#'     compare the observed view (auto-injected `obs`) against the literal
#'     no-reweighting baseline.
#'
#' @param fit A `causatr_fit` object.
#' @param intervention A `causatr_intervention` object or `NULL`.
#' @param is_default_obs Logical. `TRUE` only for the auto-injected `obs`
#'   panel under `diagnose(fit)` (no `interventions =` argument).
#' @return A `data.table` with `group`, `n`, `mean`, `sd`, `min`, `max`,
#'   `ess` columns -- one row per group (`treated`, `control`, `overall`)
#'   -- or `NULL` for non-IPW estimators.
#' @noRd
compute_panel_weights <- function(fit, intervention, is_default_obs) {
  if (fit$estimator != "ipw") {
    return(NULL)
  }
  if (is_default_obs) {
    return(compute_weight_summary_observed(fit))
  }
  compute_weight_summary_intervention(fit, intervention)
}

#' Detect a coarse treatment-type label for the diagnostic header
#'
#' @param fit A `causatr_fit` object.
#' @return Character scalar describing the treatment type
#'   (`"binary"`, `"continuous"`, `"categorical"`, `"count"`,
#'   `"multivariate"`, or `"unknown"`).
#' @noRd
detect_diag_treatment_type <- function(fit) {
  if (length(fit$treatment) > 1L) {
    return("multivariate")
  }
  # Longitudinal IPW: read from the first period's treatment model.
  # Checked before the point-IPW branch because longitudinal fits
  # store `treatment_model` as a list-of-models, not a single model.
  if (
    fit$estimator == "ipw" &&
      !is.null(fit$details$treatment_models_by_time)
  ) {
    tm_first <- fit$details$treatment_models_by_time[[1]]
    fam <- tm_first$family
    return(switch(
      fam,
      bernoulli = "binary",
      gaussian = "continuous",
      categorical = "categorical",
      poisson = "count",
      negbin = "count",
      fam
    ))
  }
  if (
    fit$estimator == "ipw" &&
      !is.null(fit$details$treatment_model) &&
      inherits(fit$details$treatment_model, "causatr_treatment_model")
  ) {
    fam <- fit$details$treatment_model$family
    return(switch(
      fam,
      bernoulli = "binary",
      gaussian = "continuous",
      categorical = "categorical",
      poisson = "count",
      negbin = "count",
      fam
    ))
  }
  # gcomp / matching: read off the raw treatment column. Matching is
  # binary-only by construction; gcomp accepts whatever the outcome
  # model can handle, so the heuristic is the same one used inside
  # `compute_positivity()` (numeric 0/1 -> binary, factor / character
  # -> categorical, otherwise -> unknown).
  trt <- fit$data[[fit$treatment[1]]]
  if (is.factor(trt) || is.character(trt)) {
    return("categorical")
  }
  if (is.numeric(trt)) {
    uniq <- unique(stats::na.omit(trt))
    if (length(uniq) == 2L && all(uniq %in% c(0, 1))) {
      return("binary")
    }
    return("continuous")
  }
  "unknown"
}

#' Compute positivity diagnostics
#'
#' @description
#' Dispatches to the appropriate positivity helper based on treatment
#' type. Binary treatments get the classic propensity-score quantile
#' table with violation counts (Crump et al. 2009). Continuous and
#' count treatments get a density-range summary: quantiles of
#' \eqn{f(A_i \mid L_i)} and a count of low-density tail
#' observations. Categorical treatments get per-level
#' \eqn{P(A = k \mid L)} quantile summaries. Multivariate
#' treatments get a per-component positivity list.
#'
#' @param fit A `causatr_fit` object.
#' @param ps_bounds Numeric vector of length 2 defining violation
#'   thresholds for binary treatments, or the density tail quantile
#'   for non-binary treatments. Default `c(0.025, 0.975)`.
#' @return A data.table (or list of data.tables for multivariate),
#'   or `NULL` when no positivity diagnostic applies.
#' @noRd
compute_positivity <- function(fit, ps_bounds) {
  treatment <- fit$treatment

  # Multivariate: per-component positivity.
  if (length(treatment) > 1L) {
    return(compute_positivity_multivariate(fit, ps_bounds))
  }

  # IPW with a treatment density model: dispatch by family.
  if (fit$estimator == "ipw" && !is.null(fit$details$treatment_model)) {
    tm <- fit$details$treatment_model
    fam <- tm$family
    if (fam == "bernoulli") {
      return(compute_positivity_binary(fit, ps_bounds))
    }
    if (fam == "gaussian") {
      return(compute_positivity_density(fit, ps_bounds))
    }
    if (fam == "categorical") {
      return(compute_positivity_categorical(fit, ps_bounds))
    }
    if (fam %in% c("poisson", "negbin")) {
      return(compute_positivity_density(fit, ps_bounds))
    }
    return(NULL)
  }

  # Non-IPW estimators (gcomp, matching): binary PS check only.
  compute_positivity_binary(fit, ps_bounds)
}

#' Compute binary propensity-score positivity diagnostics
#'
#' @description
#' The classic propensity-score positivity table: quantiles of
#' \eqn{P(A = 1 \mid L)} plus violation counts. Only applicable
#' to binary 0/1 treatments.
#'
#' @param fit A `causatr_fit` object.
#' @param ps_bounds Numeric vector of length 2.
#' @return A data.table or `NULL`.
#' @noRd
compute_positivity_binary <- function(fit, ps_bounds) {
  treatment <- fit$treatment

  if (length(treatment) > 1L) {
    return(NULL)
  }

  data <- fit$data
  trt_vals <- unique(stats::na.omit(data[[treatment]]))
  if (!all(trt_vals %in% c(0, 1))) {
    return(NULL)
  }

  # Source propensity scores from the fitted model:
  #   IPW:      bernoulli density model in fit$details$propensity_model
  #   Matching: MatchIt distance vector (logistic PS by default)
  #   Gcomp:    fit a quick logistic regression for diagnostics only
  if (fit$estimator == "ipw") {
    tm <- fit$details$treatment_model
    if (is.null(tm) || tm$family != "bernoulli") {
      return(NULL)
    }
    ps <- as.numeric(stats::predict(
      fit$details$propensity_model,
      newdata = fit$data[fit$details$fit_rows],
      type = "response"
    ))
  } else if (!is.null(fit$match_obj)) {
    ps <- fit$match_obj$distance
    if (is.null(ps) || length(ps) == 0L) return(NULL)
  } else {
    ps_formula <- build_ps_formula(fit$confounders, treatment)
    fit_rows <- get_fit_rows(data, fit$outcome, fit$censoring)
    ps_model <- stats::glm(
      ps_formula,
      data = data[fit_rows],
      family = stats::binomial()
    )
    ps <- stats::fitted(ps_model)
  }

  # Crump et al. (2009) tail thresholds.
  n_low <- sum(ps < ps_bounds[1], na.rm = TRUE)
  n_high <- sum(ps > ps_bounds[2], na.rm = TRUE)
  n_total <- length(ps)

  data.table::data.table(
    statistic = c(
      "min",
      "q25",
      "median",
      "q75",
      "max",
      "n_below_lower",
      "n_above_upper",
      "n_violations",
      "pct_violations"
    ),
    value = c(
      min(ps, na.rm = TRUE),
      stats::quantile(ps, 0.25, na.rm = TRUE),
      stats::quantile(ps, 0.50, na.rm = TRUE),
      stats::quantile(ps, 0.75, na.rm = TRUE),
      max(ps, na.rm = TRUE),
      n_low,
      n_high,
      n_low + n_high,
      round(100 * (n_low + n_high) / n_total, 2)
    )
  )
}

#' Compute density-range positivity diagnostics
#'
#' @description
#' For continuous and count treatments, "positivity" means the fitted
#' conditional density \eqn{f(A_i \mid L_i)} evaluated at the
#' observed treatment value is bounded away from zero. A very small
#' density means the observed treatment is unlikely given the
#' covariates, yielding extreme density-ratio weights. This helper
#' reports quantiles of the fitted density plus a count of low-density
#' observations (those below the 1st percentile of the density
#' distribution).
#'
#' @param fit A `causatr_fit` with an IPW treatment model.
#' @param ps_bounds Numeric vector of length 2 (used to define the
#'   low-density fraction; observations with density below the
#'   `ps_bounds[1]`-th quantile are flagged).
#' @return A data.table with density quantiles and low-density counts.
#' @noRd
compute_positivity_density <- function(fit, ps_bounds) {
  tm <- fit$details$treatment_model
  fit_rows <- fit$details$fit_rows
  fit_data <- fit$data[fit_rows]
  a_obs <- fit_data[[fit$treatment[1]]]

  f_obs <- evaluate_density(tm, a_obs, fit_data)

  # Low-density threshold: observations below the 1st percentile of
  # the density distribution. These are the analogues of extreme PS
  # tails for binary treatment.
  low_thresh <- stats::quantile(f_obs, 0.01, na.rm = TRUE)
  n_low <- sum(f_obs < low_thresh, na.rm = TRUE)
  n_total <- length(f_obs)

  data.table::data.table(
    statistic = c(
      "min",
      "q25",
      "median",
      "q75",
      "max",
      "n_low_density",
      "pct_low_density"
    ),
    value = c(
      min(f_obs, na.rm = TRUE),
      stats::quantile(f_obs, 0.25, na.rm = TRUE),
      stats::quantile(f_obs, 0.50, na.rm = TRUE),
      stats::quantile(f_obs, 0.75, na.rm = TRUE),
      max(f_obs, na.rm = TRUE),
      n_low,
      round(100 * n_low / n_total, 2)
    )
  )
}

#' Compute per-level positivity for categorical treatments
#'
#' @description
#' For categorical treatments with \eqn{k \ge 2} levels, reports
#' quantiles of \eqn{P(A = k \mid L)} for each level, plus the
#' count of low-probability cells (predicted probability below 0.01).
#' A low predicted probability for a level means some individuals
#' are very unlikely to receive that treatment level given their
#' covariates -- a positivity concern.
#'
#' @param fit A `causatr_fit` with a categorical IPW treatment model.
#' @param ps_bounds Numeric vector of length 2 (unused, kept for
#'   interface consistency).
#' @return A data.table with per-level probability summaries.
#' @noRd
compute_positivity_categorical <- function(fit, ps_bounds) {
  tm <- fit$details$treatment_model
  fit_rows <- fit$details$fit_rows
  fit_data <- fit$data[fit_rows]
  a_obs <- fit_data[[fit$treatment[1]]]
  trt_levels <- tm$levels

  # Get the full predicted probability matrix from the multinomial
  # model. `predict(type = "probs")` returns n x K for K > 2 or a
  # vector for K = 2.
  prob_raw <- stats::predict(
    tm$model,
    newdata = fit_data,
    type = "probs"
  )
  if (is.null(dim(prob_raw))) {
    prob_mat <- cbind(1 - prob_raw, prob_raw)
    colnames(prob_mat) <- trt_levels
  } else {
    prob_mat <- prob_raw
  }

  # Per-level summary: quantiles + low-probability count.
  rows <- lapply(trt_levels, function(lev) {
    p_lev <- prob_mat[, lev]
    n_low <- sum(p_lev < 0.01, na.rm = TRUE)
    data.table::data.table(
      level = lev,
      min = min(p_lev, na.rm = TRUE),
      q25 = stats::quantile(p_lev, 0.25, na.rm = TRUE),
      median = stats::quantile(p_lev, 0.50, na.rm = TRUE),
      q75 = stats::quantile(p_lev, 0.75, na.rm = TRUE),
      max = max(p_lev, na.rm = TRUE),
      n_low_prob = n_low,
      pct_low_prob = round(100 * n_low / length(p_lev), 2)
    )
  })
  data.table::rbindlist(rows)
}

#' Compute per-component positivity for multivariate treatments
#'
#' @description
#' Loops over each component's `causatr_treatment_model` and calls
#' the appropriate single-component positivity helper. Returns a
#' named list of positivity tables, one per treatment component.
#'
#' @param fit A `causatr_fit` with multivariate IPW.
#' @param ps_bounds Numeric vector of length 2.
#' @return A named list of data.tables (one per component), or NULL
#'   for non-IPW fits.
#' @noRd
compute_positivity_multivariate <- function(fit, ps_bounds) {
  if (fit$estimator != "ipw") {
    return(NULL)
  }
  tms <- fit$details$treatment_models
  if (is.null(tms)) {
    return(NULL)
  }

  # Build a fake single-component fit for each component, then
  # dispatch to the family-appropriate helper.
  result <- lapply(names(tms), function(trt_name) {
    tm_k <- tms[[trt_name]]
    fake_fit <- list(
      treatment = trt_name,
      estimator = "ipw",
      data = fit$data,
      confounders = fit$confounders,
      outcome = fit$outcome,
      details = list(
        treatment_model = tm_k,
        propensity_model = tm_k$model,
        fit_rows = fit$details$fit_rows
      )
    )
    fam <- tm_k$family
    if (fam == "bernoulli") {
      compute_positivity_binary(fake_fit, ps_bounds)
    } else if (fam == "categorical") {
      compute_positivity_categorical(fake_fit, ps_bounds)
    } else {
      compute_positivity_density(fake_fit, ps_bounds)
    }
  })
  names(result) <- names(tms)
  result
}

#' Compute covariate balance (via cobalt or simple SMD fallback)
#'
#' @param fit A `causatr_fit` object.
#' @param stats Character vector of balance statistics for cobalt.
#' @param thresholds Named list of thresholds for cobalt.
#' @return A cobalt `bal.tab` object or a data.table of SMDs.
#' @noRd
compute_balance <- function(fit, stats, thresholds) {
  if (!rlang::is_installed("cobalt")) {
    return(compute_balance_simple(fit))
  }

  if (
    fit$estimator == "ipw" &&
      (!is.null(fit$details$treatment_model) ||
        !is.null(fit$details$treatment_models))
  ) {
    # The self-contained density-ratio engine has no `weightit` object
    # to feed cobalt directly. Call the formula interface on the
    # observed treatment -- this produces the "unadjusted" balance view.
    # For multivariate, use the first treatment component only --
    # `build_ps_formula()` requires a scalar response, and cobalt's
    # formula interface handles one treatment at a time. Per-component
    # balance across all K components would need K separate bal.tab
    # calls (reserved for a future chunk).
    trt <- if (length(fit$treatment) > 1L) {
      fit$treatment[1]
    } else {
      fit$treatment
    }
    ps_formula <- build_ps_formula(fit$confounders, trt)
    fit_rows <- fit$details$fit_rows
    cobalt::bal.tab(
      ps_formula,
      data = as.data.frame(fit$data[fit_rows]),
      stats = stats,
      thresholds = thresholds,
      binary = "std"
    )
  } else if (fit$estimator == "matching" && !is.null(fit$match_obj)) {
    cobalt::bal.tab(
      fit$match_obj,
      un = TRUE,
      stats = stats,
      thresholds = thresholds,
      binary = "std"
    )
  } else {
    ps_formula <- build_ps_formula(fit$confounders, fit$treatment)
    fit_rows <- get_fit_rows(fit$data, fit$outcome, fit$censoring)
    cobalt::bal.tab(
      ps_formula,
      data = as.data.frame(fit$data[fit_rows]),
      stats = stats,
      thresholds = thresholds,
      binary = "std"
    )
  }
}

#' Compute simple balance table without cobalt
#'
#' @description
#' Minimal fallback when cobalt isn't installed. For binary
#' treatments, computes unadjusted SMDs per confounder. For
#' continuous treatments, computes Pearson correlations between
#' treatment and each confounder. For categorical / multivariate,
#' returns NULL (no simple fallback -- cobalt is needed).
#'
#' @param fit A `causatr_fit` object.
#' @return A data.table of balance metrics, or `NULL`.
#' @noRd
compute_balance_simple <- function(fit) {
  data <- fit$data
  treatment <- fit$treatment
  outcome <- fit$outcome

  if (length(treatment) > 1L) {
    return(NULL)
  }

  fit_rows <- get_fit_rows(data, outcome, fit$censoring)
  d <- data[fit_rows]
  confounder_vars <- all.vars(fit$confounders)
  confounder_vars <- intersect(confounder_vars, names(d))

  trt_vals <- unique(stats::na.omit(d[[treatment]]))
  is_binary <- all(trt_vals %in% c(0, 1))

  if (is_binary) {
    return(compute_balance_simple_binary(d, treatment, confounder_vars))
  }

  # Continuous / count: Pearson correlation between treatment and
  # each numeric confounder.
  if (is.numeric(d[[treatment]])) {
    return(compute_balance_simple_corr(d, treatment, confounder_vars))
  }

  # Categorical: no simple fallback.
  NULL
}

#' SMD balance table for binary treatments (cobalt fallback)
#'
#' @param d data.table of analysis-sample rows.
#' @param treatment Character treatment column name.
#' @param confounder_vars Character vector of confounder names.
#' @return A data.table with `variable`, `mean_treated`,
#'   `mean_control`, `smd` columns.
#' @noRd
compute_balance_simple_binary <- function(d, treatment, confounder_vars) {
  rows_1 <- d[[treatment]] == 1
  rows_0 <- d[[treatment]] == 0

  results <- lapply(confounder_vars, function(v) {
    x <- d[[v]]
    if (!is.numeric(x)) {
      return(NULL)
    }
    m1 <- mean(x[rows_1], na.rm = TRUE)
    m0 <- mean(x[rows_0], na.rm = TRUE)
    # Rosenbaum & Rubin (1985) pooled SD.
    s_pooled <- sqrt(
      (stats::var(x[rows_1], na.rm = TRUE) +
        stats::var(x[rows_0], na.rm = TRUE)) /
        2
    )
    smd <- if (s_pooled > 0) (m1 - m0) / s_pooled else NA_real_
    data.table::data.table(
      variable = v,
      mean_treated = m1,
      mean_control = m0,
      smd = smd
    )
  })

  data.table::rbindlist(results[!vapply(results, is.null, logical(1))])
}

#' Correlation balance table for continuous treatments (cobalt fallback)
#'
#' @param d data.table of analysis-sample rows.
#' @param treatment Character treatment column name.
#' @param confounder_vars Character vector of confounder names.
#' @return A data.table with `variable` and `correlation` columns.
#' @noRd
compute_balance_simple_corr <- function(d, treatment, confounder_vars) {
  a <- d[[treatment]]
  results <- lapply(confounder_vars, function(v) {
    x <- d[[v]]
    if (!is.numeric(x)) {
      return(NULL)
    }
    r <- stats::cor(a, x, use = "pairwise.complete.obs")
    data.table::data.table(variable = v, correlation = r)
  })
  data.table::rbindlist(results[!vapply(results, is.null, logical(1))])
}

#' Observed-treatment IPW weight distribution summary
#'
#' @description
#' For binary treatments, reconstructs the Horvitz-Thompson weights
#' `1/p` on treated rows and `1/(1-p)` on controls and summarises by
#' arm. For non-binary treatments (continuous, categorical, count),
#' the observed-treatment density-ratio weight is identically 1
#' (natural-course view), so the summary reports unit weights. For
#' multivariate treatments, computes the product of per-component
#' natural-course weights (also unit). Used for the auto-injected
#' `obs` panel under `diagnose(fit)`.
#'
#' @param fit A `causatr_fit` with `estimator = "ipw"`.
#' @return A `data.table` with columns `group`, `n`, `mean`, `sd`,
#'   `min`, `max`, `ess`.
#' @noRd
compute_weight_summary_observed <- function(fit) {
  # Multivariate: natural-course product weight = 1 for all.
  if (isTRUE(fit$details$is_multivariate)) {
    fit_rows <- fit$details$fit_rows
    n <- sum(fit_rows)
    return(summarise_weights_overall(rep(1, n)))
  }

  tm <- fit$details$treatment_model
  fit_rows <- fit$details$fit_rows
  fit_data <- fit$data[fit_rows]

  # Binary: Horvitz-Thompson observed-arm view.
  if (tm$family == "bernoulli") {
    a_obs <- fit_data[[fit$treatment[1]]]
    p <- as.numeric(stats::predict(
      fit$details$propensity_model,
      newdata = fit_data,
      type = "response"
    ))
    w <- ifelse(a_obs == 1, 1 / p, 1 / (1 - p))
    return(summarise_weights_by_arm(w, a_obs))
  }

  # Non-binary (continuous / categorical / count): the "observed"
  # density-ratio weight f(A|L)/f(A|L) = 1. Report unit weights so
  # the panel is still meaningful (ESS = n, no instability).
  n <- sum(fit_rows)
  summarise_weights_overall(rep(1, n))
}

#' Per-intervention IPW weight distribution summary
#'
#' @description
#' Computes the density-ratio weight for a specific intervention and
#' summarises it. Binary treatments are split by arm (treated /
#' control / overall); non-binary treatments report overall only.
#' Multivariate treatments compute the combined product weight
#' across all components.
#'
#' Under static-arm interventions on a binary treatment, the
#' Horvitz-Thompson indicator zeroes the off-arm rows. The "control"
#' row for `static(1)` reports `mean = 0`, `ess = 0` -- this is the
#' correct M-estimation view.
#'
#' @param fit A `causatr_fit` with `estimator = "ipw"`.
#' @param intervention A `causatr_intervention` object (or `NULL` for
#'   natural course). For multivariate, a named list of per-component
#'   interventions.
#' @return A `data.table` with weight summary columns.
#' @noRd
compute_weight_summary_intervention <- function(fit, intervention) {
  # Multivariate: compute combined product weight via the multivariate
  # weight engine.
  if (isTRUE(fit$details$is_multivariate)) {
    tms <- fit$details$treatment_models
    fit_rows <- fit$details$fit_rows
    fit_data <- fit$data[fit_rows]
    w <- compute_density_ratio_weights_mv(
      tms,
      fit_data,
      intervention,
      fit$estimand
    )
    return(summarise_weights_overall(w))
  }

  tm <- fit$details$treatment_model
  fit_rows <- fit$details$fit_rows
  fit_data <- fit$data[fit_rows]

  w <- compute_density_ratio_weights(
    tm,
    fit_data,
    intervention,
    fit$estimand
  )

  # Binary: split by observed treatment arm.
  if (tm$family == "bernoulli") {
    a_obs <- fit_data[[fit$treatment[1]]]
    if (!isTRUE(all(tm$fit_rows))) {
      a_obs <- a_obs[tm$fit_rows]
    }
    return(summarise_weights_by_arm(w, a_obs))
  }

  # Non-binary: overall summary only (no arm split).
  summarise_weights_overall(w)
}

#' Aggregate a weight vector into the (treated / control / overall) summary
#' table used by `compute_weight_summary_*` helpers.
#'
#' @param w Numeric weight vector aligned with `a_obs`.
#' @param a_obs Observed binary treatment vector (length matches `w`).
#' @return A `data.table` with columns `group`, `n`, `mean`, `sd`, `min`,
#'   `max`, `ess`.
#' @noRd
summarise_weights_by_arm <- function(w, a_obs) {
  # Effective sample size (Kish 1965): a weighted sample with highly
  # variable weights has fewer "effective" observations than its
  # nominal n. The formula below is the ratio (sum w)^2 / sum w^2 --
  # equals n when all weights are equal, less otherwise. The 0/0
  # guard handles the all-zero off-arm case under static-arm
  # density-ratio weights (e.g. controls under `static(1)`); ESS = 0
  # is the right reading there because no row contributes.
  ess <- function(wts) {
    s2 <- sum(wts^2)
    if (s2 == 0) {
      return(0)
    }
    sum(wts)^2 / s2
  }
  masks <- list(a_obs == 1, a_obs == 0, rep(TRUE, length(w)))
  labels <- c("treated", "control", "overall")
  rows <- lapply(seq_along(labels), function(i) {
    w_sub <- w[masks[[i]]]
    data.table::data.table(
      group = labels[i],
      n = length(w_sub),
      mean = if (length(w_sub) == 0L) NA_real_ else mean(w_sub),
      sd = if (length(w_sub) < 2L) NA_real_ else stats::sd(w_sub),
      min = if (length(w_sub) == 0L) NA_real_ else min(w_sub),
      max = if (length(w_sub) == 0L) NA_real_ else max(w_sub),
      ess = if (length(w_sub) == 0L) 0 else ess(w_sub)
    )
  })
  data.table::rbindlist(rows)
}

#' Aggregate a weight vector into a single overall row
#'
#' @description
#' Used by non-binary treatment types where splitting by "treated /
#' control" is not meaningful. Reports a single `overall` row with
#' the same column schema as `summarise_weights_by_arm()` so
#' downstream print / plot code can handle both shapes uniformly.
#'
#' @param w Numeric weight vector.
#' @return A `data.table` with a single `overall` row.
#' @noRd
summarise_weights_overall <- function(w) {
  ess_val <- if (sum(w^2) == 0) 0 else sum(w)^2 / sum(w^2)
  data.table::data.table(
    group = "overall",
    n = length(w),
    mean = mean(w),
    sd = if (length(w) < 2L) NA_real_ else stats::sd(w),
    min = min(w),
    max = max(w),
    ess = ess_val
  )
}

#' Compute matching quality metrics
#'
#' @param fit A `causatr_fit` object (matching estimator only).
#' @return A data.table with total, matched, discarded counts and retention
#'   percentage, or `NULL` for non-matching fits.
#' @noRd
compute_match_quality <- function(fit) {
  if (fit$estimator != "matching" || is.null(fit$match_obj)) {
    return(NULL)
  }

  m <- fit$match_obj
  n_total <- length(m$weights)
  n_matched <- sum(m$weights > 0)
  n_discarded <- n_total - n_matched

  data.table::data.table(
    statistic = c(
      "n_total",
      "n_matched",
      "n_discarded",
      "pct_retained"
    ),
    value = c(
      n_total,
      n_matched,
      n_discarded,
      round(100 * n_matched / n_total, 1)
    )
  )
}


#' Diagnose longitudinal (ICE / longitudinal-IPW) fits
#'
#' @description
#' Dispatches per-period positivity, balance, and weight diagnostics
#' for longitudinal fits. Each panel stores `positivity`, `balance`,
#' and `weights` as named lists keyed by time-point strings
#' (e.g. `"0"`, `"1"`). For longitudinal IPW the cumulative weight
#' summary is appended under the `"cumulative"` key.
#'
#' ICE (longitudinal g-comp) fits have no treatment model, so
#' positivity and weights are NULL; only per-period balance is
#' reported.
#'
#' @param fit A longitudinal `causatr_fit`.
#' @param interventions Named list of interventions (already resolved).
#' @param default_view Logical: TRUE when no explicit interventions
#'   were requested.
#' @param stats Character: balance statistics for cobalt.
#' @param thresholds Named numeric: balance thresholds.
#' @param ps_bounds Numeric(2): positivity bounds.
#' @return A `causatr_diag` object.
#' @noRd
diagnose_longitudinal <- function(
  fit,
  interventions,
  default_view,
  stats,
  thresholds,
  ps_bounds
) {
  is_ipw <- fit$estimator == "ipw"
  time_points <- fit$details$time_points
  tp_labels <- as.character(time_points)

  # Positivity: per-period treatment density diagnostics (IPW only).
  positivity_shared <- if (is_ipw) {
    compute_positivity_longitudinal(fit, ps_bounds)
  }

  # Balance: per-period covariate balance.
  balance_shared <- compute_balance_longitudinal(
    fit,
    stats,
    thresholds
  )

  per_intervention <- lapply(seq_along(interventions), function(i) {
    iv <- interventions[[i]]
    is_default_obs <- default_view && i == 1L
    weights_panel <- if (is_ipw) {
      compute_weights_longitudinal(fit, iv, is_default_obs)
    }
    list(
      positivity = positivity_shared,
      balance = balance_shared,
      weights = weights_panel
    )
  })
  names(per_intervention) <- names(interventions)

  fit_info <- list(
    treatment_type = detect_diag_treatment_type(fit),
    estimand = fit$estimand,
    type = fit$type,
    has_em = !is.null(fit$details$em_info) &&
      length(fit$details$em_info$em_terms) > 0L
  )

  new_causatr_diag(
    per_intervention = per_intervention,
    match_quality = NULL,
    estimator = fit$estimator,
    fit_info = fit_info,
    fit = fit
  )
}


#' Per-period positivity for longitudinal IPW
#'
#' @description
#' Loops over `treatment_models_by_time` and calls the appropriate
#' single-period positivity helper on each. Returns a named list
#' keyed by time-point string, each element a positivity data.table.
#'
#' @param fit A longitudinal IPW `causatr_fit`.
#' @param ps_bounds Numeric(2).
#' @return Named list of data.tables.
#' @noRd
compute_positivity_longitudinal <- function(fit, ps_bounds) {
  tms_by_time <- fit$details$treatment_models_by_time
  fit_data_by_time <- fit$details$fit_data_by_time

  result <- lapply(names(tms_by_time), function(tp) {
    tm_k <- tms_by_time[[tp]]
    data_k <- fit_data_by_time[[tp]]
    # Build a fake single-period fit to reuse existing helpers.
    fake_fit <- list(
      treatment = fit$treatment,
      estimator = "ipw",
      data = data_k,
      confounders = fit$confounders,
      outcome = fit$outcome,
      details = list(
        treatment_model = tm_k,
        propensity_model = tm_k$model,
        fit_rows = tm_k$fit_rows
      )
    )
    fam <- tm_k$family
    if (fam == "bernoulli") {
      compute_positivity_binary(fake_fit, ps_bounds)
    } else if (fam == "categorical") {
      compute_positivity_categorical(fake_fit, ps_bounds)
    } else {
      compute_positivity_density(fake_fit, ps_bounds)
    }
  })
  names(result) <- names(tms_by_time)
  result
}


#' Per-period covariate balance for longitudinal fits
#'
#' @description
#' For longitudinal IPW: computes unadjusted balance at each period
#' using the per-period data subset and the treatment formula for
#' that period. For ICE: computes balance at each period using the
#' time-varying confounders on the per-period data subset.
#'
#' @param fit A longitudinal `causatr_fit`.
#' @param stats Character: balance statistics for cobalt.
#' @param thresholds Named numeric: balance thresholds.
#' @return Named list of cobalt `bal.tab` objects or data.tables
#'   (one per time point), or NULL.
#' @noRd
compute_balance_longitudinal <- function(fit, stats, thresholds) {
  time_points <- fit$details$time_points
  time_col <- fit$time
  treatment <- fit$treatment
  data <- fit$data

  # Confounders for each period: baseline + time-varying.
  conf_baseline <- fit$confounders
  conf_tv <- fit$confounders_tv

  result <- lapply(as.character(time_points), function(tp) {
    # Subset to this period's rows.
    rows_k <- data[[time_col]] == as.numeric(tp)
    d_k <- data[rows_k]

    # Build the confounder set: baseline terms + time-varying terms
    # that exist as columns in the per-period data.
    baseline_vars <- all.vars(conf_baseline)
    tv_vars <- if (!is.null(conf_tv)) all.vars(conf_tv) else character(0)
    all_vars <- unique(c(baseline_vars, tv_vars))
    available <- intersect(all_vars, names(d_k))

    if (length(available) == 0L) {
      return(NULL)
    }

    # Build a formula for cobalt: treatment ~ available confounders.
    ps_formula <- stats::reformulate(available, response = treatment)

    if (rlang::is_installed("cobalt")) {
      tryCatch(
        cobalt::bal.tab(
          ps_formula,
          data = as.data.frame(d_k),
          stats = stats,
          thresholds = thresholds,
          binary = "std"
        ),
        error = function(e) NULL
      )
    } else {
      # Simple fallback.
      trt_vals <- unique(stats::na.omit(d_k[[treatment]]))
      is_binary <- all(trt_vals %in% c(0, 1))
      if (is_binary) {
        compute_balance_simple_binary(
          d_k,
          treatment,
          available
        )
      } else if (is.numeric(d_k[[treatment]])) {
        compute_balance_simple_corr(
          d_k,
          treatment,
          available
        )
      }
    }
  })
  names(result) <- as.character(time_points)
  # Drop NULL entries (periods with no confounders).
  result <- result[!vapply(result, is.null, logical(1))]
  if (length(result) == 0L) NULL else result
}


#' Per-period weight diagnostics for longitudinal IPW
#'
#' @description
#' Computes per-period density-ratio weight summaries and the
#' cumulative (product) weight summary. For the default `obs` panel,
#' each period's "observed" weight is the Horvitz-Thompson weight
#' for binary treatments or unit weight for non-binary. For
#' user-supplied interventions, per-period density-ratio weights are
#' computed and then multiplied into the cumulative weight.
#'
#' @param fit A longitudinal IPW `causatr_fit`.
#' @param intervention A `causatr_intervention` or NULL.
#' @param is_default_obs Logical.
#' @return Named list of data.tables (one per time point +
#'   `"cumulative"`).
#' @noRd
compute_weights_longitudinal <- function(
  fit,
  intervention,
  is_default_obs
) {
  tms_by_time <- fit$details$treatment_models_by_time
  fit_data_by_time <- fit$details$fit_data_by_time
  time_points <- fit$details$time_points
  tp_labels <- as.character(time_points)
  id_col <- fit$id
  stabilize <- !is.null(fit$details$numerator_models_by_time)

  # First-period ids define the canonical ordering.
  first_tp <- tp_labels[1]
  ids_first <- as.character(
    fit_data_by_time[[first_tp]][[id_col]]
  )
  n_id <- length(ids_first)

  # Per-period weight matrix (n_id x K).
  W_per_period <- matrix(1, nrow = n_id, ncol = length(tp_labels))
  period_summaries <- vector("list", length(tp_labels))
  names(period_summaries) <- tp_labels

  for (k in seq_along(tp_labels)) {
    tp <- tp_labels[k]
    tm_k <- tms_by_time[[tp]]
    data_k <- fit_data_by_time[[tp]]
    ids_k <- as.character(data_k[[id_col]])

    if (is_default_obs) {
      # Default "observed" view: HT weights for binary, unit for rest.
      if (tm_k$family == "bernoulli") {
        a_obs <- data_k[[fit$treatment[1]]]
        p <- as.numeric(stats::predict(
          tm_k$model,
          newdata = data_k,
          type = "response"
        ))
        w_k <- ifelse(a_obs == 1, 1 / p, 1 / (1 - p))
      } else {
        w_k <- rep(1, nrow(data_k))
      }
    } else if (is.null(intervention)) {
      w_k <- rep(1, nrow(data_k))
    } else {
      if (stabilize) {
        tm_num_k <- fit$details$numerator_models_by_time[[tp]]
        w_k <- compute_stabilized_period_weight(
          tm_denom = tm_k,
          tm_num = tm_num_k,
          data = data_k,
          intervention = intervention
        )
      } else {
        w_k <- compute_density_ratio_weights(
          treatment_model = tm_k,
          data = data_k,
          intervention = intervention,
          estimand = fit$estimand
        )
      }
    }

    # Align to canonical id ordering.
    period_ids <- ids_k[tm_k$fit_rows]
    w_aligned <- rep(1, n_id)
    pos <- match(period_ids, ids_first)
    w_aligned[pos] <- w_k
    W_per_period[, k] <- w_aligned

    # Per-period summary.
    if (
      tm_k$family == "bernoulli" &&
        (is_default_obs || is.null(intervention))
    ) {
      a_obs <- data_k[[fit$treatment[1]]]
      period_summaries[[tp]] <- summarise_weights_by_arm(
        w_k,
        a_obs[tm_k$fit_rows]
      )
    } else {
      period_summaries[[tp]] <- summarise_weights_overall(w_k)
    }
  }

  # Cumulative (product) weight across all periods.
  w_cumulative <- apply(W_per_period, 1, prod)
  period_summaries[["cumulative"]] <- summarise_weights_overall(
    w_cumulative
  )

  period_summaries
}
