#' Diagnostics for a fitted causal model
#'
#' @description
#' Computes diagnostics appropriate to the estimation method:
#' - **All methods**: positivity checks (flags covariate strata where the
#'   probability of treatment is near 0 or 1).
#' - `"ipw"`: covariate balance before and after weighting (via `cobalt`),
#'   weight distribution summary (mean, SD, max, effective sample size).
#' - `"matching"`: covariate balance before and after matching (via `cobalt`),
#'   match quality summary (% matched, caliper info).
#' - `"gcomp"`: unadjusted covariate imbalance between treatment groups.
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
#'     \item{`method`}{Character: the estimation method.}
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
#'               method = "ipw")
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

  positivity <- compute_positivity(fit, ps_bounds)
  balance <- compute_balance(fit, stats, thresholds)
  weights_summary <- compute_weight_summary(fit)
  match_quality <- compute_match_quality(fit)

  new_causatr_diag(
    balance = balance,
    positivity = positivity,
    weights = weights_summary,
    match_quality = match_quality,
    method = fit$method,
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

  if (length(treatment) > 1L) {
    return(NULL)
  }

  data <- fit$data
  trt_vals <- unique(stats::na.omit(data[[treatment]]))
  if (!all(trt_vals %in% c(0, 1))) {
    return(NULL)
  }

  if (!is.null(fit$weights_obj)) {
    ps <- fit$weights_obj$ps
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

  if (fit$method == "ipw" && !is.null(fit$weights_obj)) {
    cobalt::bal.tab(
      fit$weights_obj,
      un = TRUE,
      stats = stats,
      thresholds = thresholds,
      binary = "std"
    )
  } else if (fit$method == "matching" && !is.null(fit$match_obj)) {
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
  data <- fit$data
  treatment <- fit$treatment
  outcome <- fit$outcome

  if (length(treatment) > 1L) {
    return(NULL)
  }

  trt_vals <- unique(stats::na.omit(data[[treatment]]))
  if (!all(trt_vals %in% c(0, 1))) {
    return(NULL)
  }

  fit_rows <- !is.na(data[[outcome]])
  d <- data[fit_rows]
  confounder_vars <- all.vars(fit$confounders)
  confounder_vars <- intersect(confounder_vars, names(d))

  rows_1 <- d[[treatment]] == 1
  rows_0 <- d[[treatment]] == 0

  results <- lapply(confounder_vars, function(v) {
    x <- d[[v]]
    if (!is.numeric(x)) {
      return(NULL)
    }
    m1 <- mean(x[rows_1], na.rm = TRUE)
    m0 <- mean(x[rows_0], na.rm = TRUE)
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
#' @param fit A `causatr_fit` object (IPW method only).
#' @return A data.table with mean, SD, min, max, and ESS by treatment group,
#'   or `NULL` for non-IPW fits.
#' @noRd
compute_weight_summary <- function(fit) {
  if (fit$method != "ipw" || is.null(fit$weights_obj)) {
    return(NULL)
  }

  w <- fit$weights_obj$weights
  trt <- fit$data[[fit$treatment[1]]]
  fit_rows <- !is.na(fit$data[[fit$outcome]])
  trt <- trt[fit_rows]

  ess <- function(wts) sum(wts)^2 / sum(wts^2)

  treat_type <- attr(fit$weights_obj$treat, "treat.type") %||% "continuous"

  if (treat_type == "binary") {
    group_labels <- c("treated", "control", "overall")
    group_masks <- list(trt == 1, trt == 0, rep(TRUE, length(w)))
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
#' @param fit A `causatr_fit` object (matching method only).
#' @return A data.table with total, matched, discarded counts and retention
#'   percentage, or `NULL` for non-matching fits.
#' @noRd
compute_match_quality <- function(fit) {
  if (fit$method != "matching" || is.null(fit$match_obj)) {
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
