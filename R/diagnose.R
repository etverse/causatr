#' Diagnostics for a fitted causal model
#'
#' @description
#' Computes diagnostics appropriate to the causal estimator:
#' - **All estimators**: positivity checks (flags covariate strata where the
#'   probability of treatment is near 0 or 1).
#' - `"ipw"` with binary treatment: covariate balance via `cobalt` on the
#'   propensity formula, plus an observed-treatment weight distribution
#'   summary (`1/p` on treated rows, `1/(1-p)` on controls) with mean, SD,
#'   min, max, and effective sample size.
#' - `"matching"`: covariate balance before and after matching (via `cobalt`),
#'   match quality summary (% matched, caliper info).
#' - `"gcomp"`: unadjusted covariate imbalance between treatment groups.
#'
#' **Point-treatment fits only.** `diagnose()` aborts on a longitudinal
#' (ICE) `fit`. For longitudinal data, inspect per-period propensity /
#' balance tables manually.
#'
#' **Binary treatment only for IPW weight summaries.** Continuous and
#' categorical IPW fits abort with `causatr_diag_unsupported_family`: the
#' density-ratio weight under those families depends on the intervention,
#' and `diagnose()` is an intervention-free diagnostic.
#'
#' @param fit A `causatr_fit` object returned by [causat()].
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
#'     \item{`balance`}{`cobalt::bal.tab` object (if cobalt installed) or a
#'       data.table of SMDs. Covariate balance summary.}
#'     \item{`positivity`}{data.table: propensity score summary and count of
#'       near-violations.}
#'     \item{`weights`}{data.table or `NULL`: weight distribution summary
#'       (IPW only).}
#'     \item{`match_quality`}{data.table or `NULL`: match quality summary
#'       (matching only).}
#'     \item{`estimator`}{Character: the causal estimator.}
#'     \item{`fit`}{The original `causatr_fit` (stored for `plot()`).}
#'   }
#'
#' @details
#' ## Positivity
#' For binary treatment, fits a logistic regression of the treatment on the
#' confounders and flags individuals whose estimated propensity score falls
#' outside `ps_bounds`. The returned `positivity` table summarises the
#' propensity score distribution and the number/fraction of near-violations.
#'
#' ## Balance (IPW and matching)
#' If the `cobalt` package is installed, balance is computed via
#' `cobalt::bal.tab()` on the internal `weightit` or `matchit` object. This
#' provides standardised mean differences (SMD), variance ratios, and KS
#' statistics before and after adjustment. If `cobalt` is not installed, a
#' simpler data.table-based SMD comparison is returned.
#'
#' ## Weight distribution (IPW only)
#' Summarises the IPW weights: mean, SD, min, max, and the effective sample
#' size (ESS) for the treated and control groups.
#'
#' ## Match quality (matching only)
#' Reports the number matched, number discarded, and the fraction of the
#' original sample retained.
#'
#' @references
#' Greifer N (2024). cobalt: Covariate Balance Tables and Plots.
#' \url{https://ngreifer.github.io/cobalt/}
#'
#' @examples
#' \dontrun{
#' data("nhefs", package = "causatr")
#' fit <- causat(nhefs, outcome = "wt82_71", treatment = "qsmk",
#'               confounders = ~ sex + age + wt71,
#'               estimator = "ipw")
#' diag <- diagnose(fit)
#' print(diag)
#' plot(diag)
#' }
#'
#' @seealso [causat()], [plot.causatr_diag()]
#' @export
diagnose <- function(
  fit,
  stats = c("m", "v"),
  thresholds = c(m = 0.1),
  ps_bounds = c(0.025, 0.975)
) {
  if (!inherits(fit, "causatr_fit")) {
    rlang::abort("`fit` must be a `causatr_fit` object returned by `causat()`.")
  }

  # Reject longitudinal fits up front rather than running the
  # point-treatment path on them, which silently produces a degenerate
  # positivity table and crashes in cobalt's printer.
  if (identical(fit$type, "longitudinal")) {
    rlang::abort(
      c(
        "`diagnose()` is not supported for longitudinal fits.",
        i = paste0(
          "Run `diagnose()` on a point-treatment subset of the data ",
          "(e.g. baseline with `time == min(time)`)."
        )
      ),
      .call = FALSE
    )
  }

  # IPW gets a dedicated binary-static dispatch because the observed-
  # treatment weight summary `1/p` / `1/(1-p)` is only well-defined for a
  # binary-treatment bernoulli density model. Continuous and categorical
  # density families are gated out here rather than silently returning
  # half-populated tables from `compute_weight_summary()`.
  if (fit$estimator == "ipw") {
    return(diagnose_ipw_point(fit, stats, thresholds, ps_bounds))
  }

  positivity <- compute_positivity(fit, ps_bounds)
  balance <- compute_balance(fit, stats, thresholds)
  match_quality <- compute_match_quality(fit)

  new_causatr_diag(
    balance = balance,
    positivity = positivity,
    weights = NULL,
    match_quality = match_quality,
    estimator = fit$estimator,
    fit = fit
  )
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

  # Positivity is only meaningful for a single binary treatment —
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
  #               purely for diagnostics — it doesn't affect estimation.
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
    # observed treatment — this produces the "unadjusted" standardised
    # mean differences, which is the most universal balance view the
    # density-ratio engine can surface without committing to one
    # specific intervention's post-weighting balance.
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

  # Drop rows with missing outcome and censored rows — same
  # fitting-row definition as the main pipeline, so balance is
  # computed on the analysis sample (not the pre-filter raw data).
  fit_rows <- get_fit_rows(data, outcome, fit$censoring)
  d <- data[fit_rows]
  confounder_vars <- all.vars(fit$confounders)
  confounder_vars <- intersect(confounder_vars, names(d))

  rows_1 <- d[[treatment]] == 1
  rows_0 <- d[[treatment]] == 0

  # Per-confounder SMD. Factors/characters are skipped (return NULL
  # and get filtered by `!vapply(., is.null, .)` below) — a proper
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
#' Reconstructs the observed-treatment Horvitz–Thompson weights `1/p` on
#' treated rows and `1/(1-p)` on controls from the stashed bernoulli
#' density model, and summarises them by arm plus an overall row. This is
#' the canonical Hernán & Robins Ch. 12 IPW weight diagnostic, decoupled
#' from any specific intervention.
#'
#' @param fit A `causatr_fit` with `estimator = "ipw"` and a bernoulli
#'   treatment model.
#' @return A `data.table` with columns `group`, `n`, `mean`, `sd`, `min`,
#'   `max`, `ess` — one row per group (`treated`, `control`, `overall`).
#'
#' @noRd
compute_weight_summary <- function(fit) {
  fit_rows <- fit$details$fit_rows
  fit_data <- fit$data[fit_rows]
  a_obs <- fit_data[[fit$treatment[1]]]

  # Effective sample size (Kish 1965): a weighted sample with highly
  # variable weights has fewer "effective" observations than its
  # nominal n. The formula below is the ratio (sum w)^2 / sum w^2 —
  # equals n when all weights are equal, less otherwise. Used as a
  # quick heuristic for "did IPW destroy my power?".
  ess <- function(wts) sum(wts)^2 / sum(wts^2)

  p <- as.numeric(stats::predict(
    fit$details$propensity_model,
    newdata = fit_data,
    type = "response"
  ))
  w <- ifelse(a_obs == 1, 1 / p, 1 / (1 - p))
  masks <- list(a_obs == 1, a_obs == 0, rep(TRUE, length(w)))
  labels <- c("treated", "control", "overall")

  rows <- lapply(seq_along(labels), function(i) {
    w_sub <- w[masks[[i]]]
    data.table::data.table(
      group = labels[i],
      n = length(w_sub),
      mean = mean(w_sub),
      sd = stats::sd(w_sub),
      min = min(w_sub),
      max = max(w_sub),
      ess = ess(w_sub)
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

#' Diagnose a self-contained IPW point-treatment fit (binary)
#'
#' @description
#' Assembles positivity, balance, and weight-distribution diagnostics for
#' an `estimator = "ipw"` fit whose treatment density model is bernoulli.
#' Continuous (gaussian) and categorical treatment families are rejected
#' with `causatr_diag_unsupported_family` because the density-ratio
#' weight under those families depends on the intervention, and
#' `diagnose()` is an intervention-free diagnostic.
#'
#' @param fit A `causatr_fit` with `estimator = "ipw"`.
#' @param stats,thresholds,ps_bounds Forwarded from `diagnose()`.
#' @return A `causatr_diag` object with the usual slots
#'   (`balance`, `positivity`, `weights`, `match_quality = NULL`).
#'
#' @noRd
diagnose_ipw_point <- function(fit, stats, thresholds, ps_bounds) {
  tm <- fit$details$treatment_model
  if (is.null(tm) || tm$family != "bernoulli") {
    fam <- if (is.null(tm)) "<unknown>" else tm$family
    rlang::abort(
      paste0(
        "Weight diagnostics are supported for binary treatment under ",
        "`estimator = \"ipw\"`; this fit uses treatment family '",
        fam,
        "'."
      ),
      class = "causatr_diag_unsupported_family",
      .call = FALSE
    )
  }

  # Binary bernoulli IPW: positivity uses the stashed PS, balance uses
  # the formula interface on the fit rows, and the weight summary uses
  # the observed-treatment Horvitz–Thompson weights. All three helpers
  # already know how to read from `fit$details`, so this wrapper is a
  # thin dispatcher that forwards their outputs into `new_causatr_diag()`.
  positivity <- compute_positivity(fit, ps_bounds)
  balance <- compute_balance(fit, stats, thresholds)
  weights_summary <- compute_weight_summary(fit)

  new_causatr_diag(
    balance = balance,
    positivity = positivity,
    weights = weights_summary,
    match_quality = NULL,
    estimator = fit$estimator,
    fit = fit
  )
}
