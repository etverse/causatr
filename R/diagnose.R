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
#' @param by Character scalar naming a baseline variable for stratified
#'   balance reporting. When non-`NULL`, covariate balance is computed
#'   within each stratum of `by` via `cobalt::bal.tab(..., cluster = by)`.
#'   The variable must be present in the data.
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
#'       backward compatibility with the pre-restructure flat shape.}
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
#' comparison is returned. Balance is the unadjusted SMD across treatment
#' groups; post-weighting balance under specific interventions or estimands
#' is not computed.
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

  check_diag_by_arg(fit, by, call = rlang::current_env())

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
  # the unadjusted balance view. Both are shared across
  # every panel below.
  positivity_shared <- compute_positivity(fit, ps_bounds)
  balance_shared <- compute_balance(fit, stats, thresholds, by = by)
  censoring_shared <- compute_censoring_diagnostics(fit)

  per_intervention <- lapply(seq_along(interventions), function(i) {
    nm <- names(interventions)[i]
    iv <- interventions[[i]]
    is_default_obs <- default_view && i == 1L
    weights_panel <- compute_panel_weights(fit, iv, is_default_obs)
    pct_panel <- if (is_default_obs) NULL else compute_pct_intervened(fit, iv)
    list(
      positivity = positivity_shared,
      balance = balance_shared,
      weights = weights_panel,
      censoring = censoring_shared,
      pct_intervened = pct_panel
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

#' Validate `by =` argument for stratified diagnostics
#'
#' @description
#' Checks that the `by` variable exists in the data and is a
#' baseline variable. Aborts with an informative error if the
#' variable is missing or not present in the data.
#'
#' @param fit A `causatr_fit` object.
#' @param by Character scalar or `NULL`.
#' @param call Caller environment for error attribution.
#' @return `NULL` invisibly.
#' @noRd
check_diag_by_arg <- function(fit, by, call) {
  if (is.null(by)) {
    return(invisible(NULL))
  }

  if (!is.character(by) || length(by) != 1L) {
    rlang::abort(
      "`by` must be a single character string naming a variable.",
      call = call
    )
  }

  if (!by %in% names(fit$data)) {
    rlang::abort(
      c(
        paste0("`by = \"", by, "\"` not found in the data."),
        i = "The `by` variable must be a column in the data passed to `causat()`."
      ),
      call = call
    )
  }

  invisible(NULL)
}

#' Build the per-intervention weight panel for a single intervention key
#'
#' @description
#' Three dispatch axes:
#'  1. Estimator: only `"ipw"` produces a non-NULL panel;
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
