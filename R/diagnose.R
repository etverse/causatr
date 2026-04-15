#' Diagnostics for a fitted causal model
#'
#' @description
#' Computes diagnostics appropriate to the causal estimator:
#' - **All estimators**: positivity checks (flags covariate strata where the
#'   probability of treatment is near 0 or 1).
#' - `"ipw"`: covariate balance before and after weighting (via `cobalt`),
#'   weight distribution summary (mean, SD, max, effective sample size).
#' - `"matching"`: covariate balance before and after matching (via `cobalt`),
#'   match quality summary (% matched, caliper info).
#' - `"gcomp"`: unadjusted covariate imbalance between treatment groups.
#'
#' **Point-treatment fits only.** `diagnose()` aborts on a longitudinal
#' (ICE) `fit`. For longitudinal data, inspect per-period propensity /
#' balance tables manually.
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

  positivity <- compute_positivity(fit, ps_bounds)
  balance <- compute_balance(fit, stats, thresholds)
  weights_summary <- compute_weight_summary(fit)
  match_quality <- compute_match_quality(fit)

  new_causatr_diag(
    balance = balance,
    positivity = positivity,
    weights = weights_summary,
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
  #   - IPW:      reuse the WeightIt-computed PS directly.
  #   - Matching: reuse the MatchIt distance vector (PS when the
  #               distance method is logistic, which is the default).
  #   - G-comp:   no PS was fit, so run a quick logistic regression
  #               of treatment on confounders to get one. This is
  #               purely for diagnostics — it doesn't affect estimation.
  if (!is.null(fit$weights_obj)) {
    ps <- fit$weights_obj$ps
  } else if (fit$estimator == "ipw" && !is.null(fit$details$propensity_model)) {
    # Self-contained IPW: propensity scores come from the stashed
    # treatment density model. For binary treatment the response-
    # scale prediction is the propensity directly.
    tm <- fit$details$treatment_model
    if (!is.null(tm) && tm$family == "bernoulli") {
      ps <- as.numeric(stats::predict(
        fit$details$propensity_model,
        newdata = fit$data[fit$details$fit_rows],
        type = "response"
      ))
    } else {
      return(NULL)
    }
  } else if (!is.null(fit$match_obj)) {
    ps <- fit$match_obj$distance
    if (is.null(ps) || length(ps) == 0L) return(NULL)
  } else {
    ps_formula <- build_ps_formula(fit$confounders, treatment)
    fit_rows <- get_fit_rows(data, fit$outcome)
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

  if (fit$estimator == "ipw" && !is.null(fit$weights_obj)) {
    cobalt::bal.tab(
      fit$weights_obj,
      un = TRUE,
      stats = stats,
      thresholds = thresholds,
      binary = "std"
    )
  } else if (fit$estimator == "ipw" && !is.null(fit$details$treatment_model)) {
    # Self-contained IPW: no WeightIt weightit object to feed cobalt.
    # Fall back to the formula-based `bal.tab()` call on the observed
    # treatment — this produces the "unadjusted" balance table
    # (standardised mean differences before weighting), which is the
    # most universal diagnostic the density-ratio engine can hand
    # cobalt without committing to one specific intervention's
    # post-weighting balance.
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
    fit_rows <- get_fit_rows(fit$data, fit$outcome)
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

  # Drop rows with missing outcome — same fitting-row definition as
  # the main pipeline, so balance is computed on the analysis sample
  # (not the pre-filter raw data).
  fit_rows <- !is.na(data[[outcome]])
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

#' Compute IPW weight distribution summary
#'
#' @param fit A `causatr_fit` object (IPW estimator only).
#' @return A data.table with mean, SD, min, max, and ESS by treatment group,
#'   or `NULL` for non-IPW fits.
#' @noRd
compute_weight_summary <- function(fit) {
  if (fit$estimator != "ipw") {
    return(NULL)
  }

  if (is.null(fit$weights_obj)) {
    # Self-contained IPW: reconstruct the observed-treatment IPW
    # weights from the stashed density model. Binary natural-course
    # ATE weights are `1/p` for treated rows and `1/(1-p)` for
    # controls; continuous / categorical treatment families are not
    # summarised here because "the IPW weight" depends on the
    # intervention, which `diagnose()` is not scoped to.
    tm <- fit$details$treatment_model
    if (is.null(tm)) {
      return(NULL)
    }
    fit_rows <- fit$details$fit_rows
    fit_data <- fit$data[fit_rows]
    a_obs <- fit_data[[fit$treatment[1]]]
    ess <- function(wts) sum(wts)^2 / sum(wts^2)

    if (tm$family == "bernoulli") {
      # Binary ATE weights: `1/p` for treated, `1/(1-p)` for
      # controls. This is the canonical "observed-treatment" IPW
      # weight summary from Hernán & Robins Ch. 12.
      p <- as.numeric(stats::predict(
        fit$details$propensity_model,
        newdata = fit_data,
        type = "response"
      ))
      w <- ifelse(a_obs == 1, 1 / p, 1 / (1 - p))
      masks <- list(a_obs == 1, a_obs == 0, rep(TRUE, length(w)))
      labels <- c("treated", "control", "overall")
    } else if (tm$family == "gaussian") {
      # Continuous treatment: inverse-density weight `1/f(A|L)` at
      # the observed treatment. The density is a Gaussian pdf with
      # conditional mean from the fitted linear model and residual
      # SD from `tm$sigma`. Only an "overall" summary row is
      # returned — there are no discrete groups to stratify on.
      f_obs <- evaluate_density(tm, a_obs, fit_data)
      w <- 1 / f_obs
      masks <- list(rep(TRUE, length(w)))
      labels <- "overall"
    } else {
      return(NULL)
    }

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
    return(data.table::rbindlist(rows))
  }

  w <- fit$weights_obj$weights
  trt <- fit$data[[fit$treatment[1]]]
  fit_rows <- !is.na(fit$data[[fit$outcome]])
  trt <- trt[fit_rows]

  # Effective sample size (Kish 1965): a weighted sample with highly
  # variable weights has fewer "effective" observations than its
  # nominal n. The formula below is the ratio (sum w)^2 / sum w^2 —
  # equals n when all weights are equal, less otherwise. Used as a
  # quick heuristic for "did IPW destroy my power?".
  ess <- function(wts) sum(wts)^2 / sum(wts^2)

  # WeightIt's `$treat` attribute carries a `treat.type` tag that
  # tells us whether the treatment is binary, multinomial, or
  # continuous. The weight summary uses it to pick group labels:
  # by arm (binary / multinomial) or overall only (continuous, no
  # discrete groups). A missing tag indicates a non-standard or
  # serialized WeightIt fit and would silently mislabel binary /
  # multinomial weight summaries, so abort instead of defaulting.
  treat_type <- attr(fit$weights_obj$treat, "treat.type")
  if (is.null(treat_type)) {
    rlang::abort(
      paste0(
        "WeightIt object is missing the `treat.type` attribute. This ",
        "indicates a non-standard or serialized WeightIt fit. Refit ",
        "the model with `causat(..., estimator = 'ipw')` so causatr can ",
        "label weight summaries correctly."
      ),
      .call = FALSE
    )
  }

  if (treat_type == "binary") {
    # WeightIt encodes the treatment in `fit$weights_obj$treat` and
    # exposes its levels via `levels(treat)` (factor) or via
    # `unique(treat)` (numeric). Previous versions of this function
    # hard-coded `trt == 1` / `trt == 0`, which silently produced empty
    # "treated" / "control" masks for any non-0/1 binary coding (e.g.
    # factor levels `c("control", "treated")` or integer codes
    # `c(1, 2)`). We now pull the two levels from WeightIt directly and
    # pick the second one as "treated" (matching `glm()`'s convention
    # that the second factor level is the one being modelled).
    trt_col <- fit$weights_obj$treat
    lev <- if (is.factor(trt_col)) {
      levels(trt_col)
    } else {
      sort(unique(stats::na.omit(as.vector(trt_col))))
    }
    if (length(lev) != 2L) {
      rlang::abort(
        "Binary WeightIt treatment must have exactly 2 levels."
      )
    }
    treated_lev <- lev[2]
    control_lev <- lev[1]
    # Labels stay as the canonical "treated" / "control" / "overall"
    # triple (snapshot-tested and print-friendly); the numeric / factor
    # level value enters the mask only, not the label. Downstream
    # `print.causatr_diag()` can always display both by joining on
    # this summary with the user's own metadata if needed.
    group_labels <- c("treated", "control", "overall")
    group_masks <- list(
      trt == treated_lev,
      trt == control_lev,
      rep(TRUE, length(w))
    )
  } else if (treat_type == "multinomial") {
    trt_levels <- sort(unique(stats::na.omit(trt)))
    group_labels <- c(as.character(trt_levels), "overall")
    group_masks <- c(
      lapply(trt_levels, function(lev) trt == lev),
      list(rep(TRUE, length(w)))
    )
  } else {
    group_labels <- "overall"
    group_masks <- list(rep(TRUE, length(w)))
  }

  rows <- lapply(seq_along(group_labels), function(i) {
    w_sub <- w[group_masks[[i]]]
    data.table::data.table(
      group = group_labels[i],
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
