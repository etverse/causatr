#' Diagnostics for a fitted causal model
#'
#' @description
#' Computes diagnostics appropriate to the causal estimator:
#' - **All estimators**: positivity checks (flags covariate strata where the
#'   probability of treatment is near 0 or 1).
#' - `"ipw"` with binary treatment: covariate balance via `cobalt` on the
#'   propensity formula, plus a weight-distribution summary (one panel per
#'   intervention, or the observed-treatment Horvitz-Thompson summary when
#'   `interventions = NULL`).
#' - `"matching"`: covariate balance before and after matching (via `cobalt`),
#'   match quality summary (% matched, caliper info).
#' - `"gcomp"`: unadjusted covariate imbalance between treatment groups.
#'
#' **Per-intervention dispatch.** When `interventions =` is supplied, each
#' intervention spawns its own diagnostic panel: the propensity-score summary
#' is shared across panels (it depends only on the fitted treatment density
#' model, not the intervention), but the weight-distribution summary is
#' computed via `compute_density_ratio_weights()` for that specific
#' intervention. The default call `diagnose(fit)` produces a single `obs`
#' panel using the standard Horvitz-Thompson observed-treatment view
#' (`1/p` on treated rows, `1/(1-p)` on controls) for IPW.
#'
#' **Point-treatment and binary IPW only in this chunk.** The full rewrite
#' (continuous / categorical / count / multivariate IPW; longitudinal ICE;
#' ATT / ATC balance; effect-modification stratification; per-intervention
#' propensity histograms; weight-distribution plots; `vignettes/diagnostics.qmd`)
#' lands in subsequent chunks. Each currently-unsupported combination aborts
#' with a `causatr_diag_*_pending` classed error.
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

  # Order matters: longitudinal is checked first because the
  # downstream `treatment_model` family / multivariate gates only
  # apply to point-treatment IPW. ICE longitudinal fits do not have
  # a flat `treatment_model` slot to inspect. Caller env is threaded
  # through so the abort messages identify `diagnose()` as the
  # offending call rather than the internal helper.
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

#' Reject combinations not yet supported by the chunk 11a `diagnose()` rewrite
#'
#' @description
#' Centralises the up-front rejection gates so each unsupported combination
#' aborts with a stable `causatr_diag_*_pending` classed error. Order is
#' load-bearing: longitudinal is checked first because subsequent gates
#' inspect `fit$details$treatment_model`, which only exists for point-
#' treatment IPW. Inside IPW the multivariate gate beats the family gate
#' (multivariate fits store `treatment_models`, not a singular
#' `treatment_model`).
#'
#' @param fit A `causatr_fit` object.
#' @param interventions Named list (or `NULL`).
#' @param by Character scalar or `NULL`.
#' @param call Caller environment for error attribution (so the
#'   user-visible message names `diagnose()` rather than the helper).
#' @return `NULL` invisibly; aborts with the appropriate classed error.
#' @noRd
check_diag_pending_gates <- function(fit, interventions, by, call) {
  if (identical(fit$type, "longitudinal")) {
    rlang::abort(
      c(
        "`diagnose()` for longitudinal fits is not yet implemented.",
        i = paste0(
          "Per-period positivity / balance / weight diagnostics for ICE ",
          "fits land in a later chunk of the `diagnose()` rewrite."
        )
      ),
      class = "causatr_diag_longitudinal_pending",
      call = call
    )
  }

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

  if (isTRUE(fit$details$is_multivariate)) {
    rlang::abort(
      c(
        "`diagnose()` for multivariate IPW is not yet implemented.",
        i = paste0(
          "Per-component propensity / weight diagnostics for multivariate ",
          "IPW land in a later chunk."
        )
      ),
      class = "causatr_diag_multivariate_pending",
      call = call
    )
  }

  tm <- fit$details$treatment_model
  if (is.null(tm)) {
    rlang::abort(
      "Internal error: IPW point-treatment fit is missing `details$treatment_model`."
    )
  }
  fam <- tm$family
  if (identical(fam, "gaussian")) {
    rlang::abort(
      c(
        "`diagnose()` for continuous-treatment IPW is not yet implemented.",
        i = paste0(
          "Density-range positivity and weight-distribution diagnostics for ",
          "gaussian treatment models land in a later chunk."
        )
      ),
      class = "causatr_diag_continuous_pending",
      call = call
    )
  }
  if (identical(fam, "categorical")) {
    rlang::abort(
      c(
        "`diagnose()` for categorical-treatment IPW is not yet implemented.",
        i = paste0(
          "Per-level positivity and weight-distribution diagnostics for ",
          "multinomial treatment models land in a later chunk."
        )
      ),
      class = "causatr_diag_categorical_pending",
      call = call
    )
  }
  if (fam %in% c("poisson", "negbin")) {
    rlang::abort(
      c(
        "`diagnose()` for count-treatment IPW is not yet implemented.",
        i = paste0(
          "Discrete-density positivity and weight-distribution diagnostics ",
          "for Poisson / negative-binomial treatment models land in a later ",
          "chunk."
        )
      ),
      class = "causatr_diag_count_pending",
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
  if (fit$estimator == "ipw" && !is.null(fit$details$treatment_model)) {
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

#' Compute propensity score positivity diagnostics
#'
#' @param fit A `causatr_fit` object (binary treatment only).
#' @param ps_bounds Numeric vector of length 2 defining violation thresholds.
#' @return A data.table with PS quantiles and violation counts, or `NULL` for
#'   non-binary treatments.
#' @noRd
compute_positivity <- function(fit, ps_bounds) {
  treatment <- fit$treatment

  # Positivity is only meaningful for a single binary treatment --
  # that's where "probability of treatment" is a scalar we can
  # threshold. Multivariate / categorical / continuous treatments
  # return NULL (diagnose() skips positivity for those).
  if (length(treatment) > 1L) {
    return(NULL)
  }

  data <- fit$data
  trt_vals <- unique(stats::na.omit(data[[treatment]]))
  if (!all(trt_vals %in% c(0, 1))) {
    return(NULL)
  }

  # Source the propensity scores from whatever's already been fit:
  #   - IPW:      pull from the self-contained bernoulli density model
  #               stashed in `fit$details$propensity_model`.
  #   - Matching: reuse the MatchIt distance vector (PS when the
  #               distance method is logistic, which is the default).
  #   - G-comp:   no PS was fit, so run a quick logistic regression
  #               of treatment on confounders to get one. This is
  #               purely for diagnostics -- it doesn't affect estimation.
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

  # Count violations on both tails. `ps_bounds = c(0.025, 0.975)` is
  # the Crump et al. (2009) default for "extreme overlap zones";
  # individuals outside this range have near-zero probability of
  # being in one arm and contribute unstable weights (or can't be
  # matched at all).
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

  if (fit$estimator == "ipw" && !is.null(fit$details$treatment_model)) {
    # The self-contained density-ratio engine has no `weightit` object
    # to feed cobalt directly. Call the formula interface on the
    # observed treatment -- this produces the "unadjusted" standardised
    # mean differences, which is the most universal balance view the
    # density-ratio engine can surface without committing to one
    # specific intervention's post-weighting balance. Per-intervention
    # post-weighting balance lands in a later chunk.
    ps_formula <- build_ps_formula(fit$confounders, fit$treatment)
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

#' Compute simple SMD balance table without cobalt
#'
#' @param fit A `causatr_fit` object with binary treatment.
#' @return A data.table with columns `variable`, `mean_treated`,
#'   `mean_control`, and `smd`, or `NULL` for non-binary treatments.
#' @noRd
compute_balance_simple <- function(fit) {
  # Minimal fallback when cobalt isn't installed: compute unadjusted
  # standardised mean differences (SMDs) for each confounder, one
  # row per confounder. Not as rich as cobalt (no weighted SMDs, no
  # variance ratios, no factor-level breakdown), but enough to
  # surface gross imbalance in a no-dependency environment.
  data <- fit$data
  treatment <- fit$treatment
  outcome <- fit$outcome

  # Same narrowing as compute_positivity: binary scalar treatment only.
  if (length(treatment) > 1L) {
    return(NULL)
  }

  trt_vals <- unique(stats::na.omit(data[[treatment]]))
  if (!all(trt_vals %in% c(0, 1))) {
    return(NULL)
  }

  # Drop rows with missing outcome and censored rows -- same
  # fitting-row definition as the main pipeline, so balance is
  # computed on the analysis sample (not the pre-filter raw data).
  fit_rows <- get_fit_rows(data, outcome, fit$censoring)
  d <- data[fit_rows]
  confounder_vars <- all.vars(fit$confounders)
  confounder_vars <- intersect(confounder_vars, names(d))

  rows_1 <- d[[treatment]] == 1
  rows_0 <- d[[treatment]] == 0

  # Per-confounder SMD. Factors/characters are skipped (return NULL
  # and get filtered by `!vapply(., is.null, .)` below) -- a proper
  # balance table for categoricals would split into levels, which is
  # what cobalt::bal.tab() does.
  results <- lapply(confounder_vars, function(v) {
    x <- d[[v]]
    if (!is.numeric(x)) {
      return(NULL)
    }
    m1 <- mean(x[rows_1], na.rm = TRUE)
    m0 <- mean(x[rows_0], na.rm = TRUE)
    # Pooled SD denominator per Rosenbaum & Rubin (1985): sqrt of
    # the average of per-group variances. Threshold convention is
    # |SMD| < 0.1 for good balance.
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

#' Observed-treatment IPW weight distribution summary (binary)
#'
#' @description
#' Reconstructs the observed-treatment Horvitz-Thompson weights `1/p` on
#' treated rows and `1/(1-p)` on controls from the stashed bernoulli
#' density model, and summarises them by arm plus an overall row. This is
#' the canonical Hernan & Robins Ch. 12 IPW weight diagnostic, decoupled
#' from any specific intervention. Used for the auto-injected `obs` panel
#' under `diagnose(fit)` (no `interventions =` argument).
#'
#' @param fit A `causatr_fit` with `estimator = "ipw"` and a bernoulli
#'   treatment model.
#' @return A `data.table` with columns `group`, `n`, `mean`, `sd`, `min`,
#'   `max`, `ess` -- one row per group (`treated`, `control`, `overall`).
#'
#' @noRd
compute_weight_summary_observed <- function(fit) {
  fit_rows <- fit$details$fit_rows
  fit_data <- fit$data[fit_rows]
  a_obs <- fit_data[[fit$treatment[1]]]

  p <- as.numeric(stats::predict(
    fit$details$propensity_model,
    newdata = fit_data,
    type = "response"
  ))
  w <- ifelse(a_obs == 1, 1 / p, 1 / (1 - p))
  summarise_weights_by_arm(w, a_obs)
}

#' Per-intervention IPW weight distribution summary (binary)
#'
#' @description
#' Reconstructs the per-intervention density-ratio weight via
#' `compute_density_ratio_weights()` and summarises it by arm plus an
#' overall row. The summary structure (treated / control / overall with
#' ESS) mirrors `compute_weight_summary_observed()` so downstream
#' consumers can render either panel uniformly.
#'
#' Under static-arm interventions on a binary treatment, the
#' Horvitz-Thompson indicator zeroes the off-arm rows. The "control" row
#' for `static(1)` therefore reports `mean = 0`, `ess = 0` -- this is the
#' correct M-estimation view of the ATE-arm functional and is intended,
#' not a bug.
#'
#' @param fit A `causatr_fit` with `estimator = "ipw"` and a bernoulli
#'   treatment model.
#' @param intervention A `causatr_intervention` object (or `NULL` for
#'   literal natural course).
#' @return A `data.table` with the same columns as
#'   `compute_weight_summary_observed()`.
#'
#' @noRd
compute_weight_summary_intervention <- function(fit, intervention) {
  tm <- fit$details$treatment_model
  fit_rows <- fit$details$fit_rows
  fit_data <- fit$data[fit_rows]
  a_obs <- fit_data[[fit$treatment[1]]]

  # `compute_density_ratio_weights()` indexes into `data` via
  # `tm$fit_rows` (relative to `fit_data`, i.e. the outcome-filtered
  # subset, per the 2026-04-16 outcome-NA fix). When the propensity
  # mask is a strict subset of the outcome mask (e.g. confounder NAs)
  # the returned weight vector is shorter than `a_obs`; align by
  # subsetting `a_obs` to the same propensity mask before grouping.
  w <- compute_density_ratio_weights(
    tm,
    fit_data,
    intervention,
    fit$estimand
  )
  if (!isTRUE(all(tm$fit_rows))) {
    a_obs <- a_obs[tm$fit_rows]
  }
  summarise_weights_by_arm(w, a_obs)
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
