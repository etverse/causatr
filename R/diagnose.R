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
    confounder_terms <- attr(stats::terms(fit$confounders), "term.labels")
    ps_formula <- stats::reformulate(confounder_terms, response = treatment)
    fit_rows <- !is.na(data[[fit$outcome]])
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
    confounder_terms <- attr(stats::terms(fit$confounders), "term.labels")
    ps_formula <- stats::reformulate(
      confounder_terms,
      response = fit$treatment
    )
    fit_rows <- !is.na(fit$data[[fit$outcome]])
    cobalt::bal.tab(
      ps_formula,
      data = as.data.frame(fit$data[fit_rows]),
      stats = stats,
      thresholds = thresholds,
      binary = "std"
    )
  }
}

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

  data.table::data.table(
    group = c("treated", "control", "overall"),
    n = c(sum(trt == 1), sum(trt == 0), length(w)),
    mean = c(
      mean(w[trt == 1]),
      mean(w[trt == 0]),
      mean(w)
    ),
    sd = c(
      stats::sd(w[trt == 1]),
      stats::sd(w[trt == 0]),
      stats::sd(w)
    ),
    min = c(
      min(w[trt == 1]),
      min(w[trt == 0]),
      min(w)
    ),
    max = c(
      max(w[trt == 1]),
      max(w[trt == 0]),
      max(w)
    ),
    ess = c(
      ess(w[trt == 1]),
      ess(w[trt == 0]),
      ess(w)
    )
  )
}

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
