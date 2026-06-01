#' Tidy a causatr result
#'
#' @description
#' Returns a tidy data frame of intervention means and/or pairwise contrasts
#' from a `causatr_result` object, compatible with the
#' [broom](https://broom.tidymodels.org/) ecosystem.
#'
#' @param x A `causatr_result` object.
#' @param which Character. What to tidy: `"contrasts"` (default),
#'   `"means"`, or `"all"` (both).
#' @param conf.int Logical. Include confidence interval columns? Default
#'   `TRUE`.
#' @param conf.level Numeric. Confidence level for intervals. Default `0.95`.
#' @param ... Currently unused.
#'
#' @return A data.frame with columns `term`, `estimate`, `std.error`,
#'   `conf.low`, `conf.high`, and `type` (either `"mean"` or `"contrast"`).
#'
#' @examples
#' \dontrun{
#' result <- contrast(fit, interventions = list(a1 = static(1), a0 = static(0)))
#' tidy(result)
#' tidy(result, which = "means")
#' tidy(result, conf.level = 0.99)
#' }
#'
#' @seealso [coef.causatr_result()], [confint.causatr_result()]
#' @export
tidy.causatr_result <- function(
  x,
  which = c("contrasts", "means", "all"),
  conf.int = TRUE,
  conf.level = 0.95,
  ...
) {
  which <- rlang::arg_match(which)
  # Two-sided normal critical value for means CIs. Contrast CIs are
  # read from the stored `ci_lower`/`ci_upper` (already computed with
  # log-scale delta for ratios/ORs), but means use Wald symmetric
  # intervals recomputed here in case `conf.level` differs from what
  # `contrast()` originally used.
  z <- stats::qnorm((1 + conf.level) / 2)

  # Multiple-imputation results use per-row t critical values (Barnard-Rubin
  # or Boot MI degrees of freedom) instead of the normal quantile, so the
  # tidy CIs match what the pooling engines stored. `z_means`/`z_contr` fall
  # back to the scalar normal `z` for ordinary (non-MI) results.
  alpha_t <- (1 - conf.level) / 2
  mi <- attr(x, "mi_details")
  is_mi <- !is.null(mi) && x$ci_method %in% c("rubin", "boot_mi")
  z_means <- if (is_mi) stats::qt(1 - alpha_t, df = mi$estimates$df) else z
  z_contr <- if (is_mi && !is.null(mi$contrasts)) {
    stats::qt(1 - alpha_t, df = mi$contrasts$df)
  } else {
    z
  }

  # SNM results have a "parameter" column; other estimators use
  # "intervention". Pick whichever is present.
  term_col <- if ("parameter" %in% names(x$estimates)) {
    x$estimates$parameter
  } else {
    x$estimates$intervention
  }
  row_type <- if (identical(x$estimator, "snm")) "parameter" else "mean"

  # Mean (or blip parameter) row per entry. `std.error` is the
  # broom-convention name for SE; our internal table uses `se`.
  # `type` labels each row so a consumer can filter `means` vs
  # `contrasts` vs `parameter` after the `all` rbind below.
  means_df <- data.frame(
    term = term_col,
    estimate = x$estimates$estimate,
    std.error = x$estimates$se,
    type = row_type,
    stringsAsFactors = FALSE
  )
  if (conf.int) {
    means_df$conf.low <- means_df$estimate - z_means * means_df$std.error
    means_df$conf.high <- means_df$estimate + z_means * means_df$std.error
  }

  # Contrast CIs: for difference contrasts, recompute at the user's
  # conf.level (symmetric Wald). For ratio/OR, the stored CIs came
  # through the log-scale delta method in contrast(), which symmetric
  # Wald cannot reproduce — keep the original CIs.
  n_contrasts <- nrow(x$contrasts)
  contrasts_df <- data.frame(
    term = x$contrasts$comparison,
    estimate = x$contrasts$estimate,
    std.error = x$contrasts$se,
    type = rep("contrast", n_contrasts),
    stringsAsFactors = FALSE
  )
  if (conf.int) {
    if (x$type == "difference") {
      contrasts_df$conf.low <- x$contrasts$estimate - z_contr * x$contrasts$se
      contrasts_df$conf.high <- x$contrasts$estimate + z_contr * x$contrasts$se
    } else {
      # Ratio/OR: log-scale delta method CI.
      # se_log = se_linear / estimate (derivable from stored fields).
      se_log <- x$contrasts$se / x$contrasts$estimate
      contrasts_df$conf.low <- exp(log(x$contrasts$estimate) - z_contr * se_log)
      contrasts_df$conf.high <- exp(
        log(x$contrasts$estimate) + z_contr * se_log
      )
    }
  }

  # Carry the `by` subgroup column through if it's present in the
  # source. This means `tidy()` + `dplyr::filter(by == "...")` Just
  # Works for effect-modification results.
  has_by <- "by" %in% names(x$estimates)
  if (has_by) {
    means_df$by <- x$estimates$by
    contrasts_df$by <- x$contrasts$by
  }

  # Switch on the requested shape. `"all"` rbinds means first then
  # contrasts, consistent with how print.causatr_result shows them.
  result <- switch(
    which,
    contrasts = contrasts_df,
    means = means_df,
    all = rbind(means_df, contrasts_df)
  )

  # `rbind()` preserves the original rownames (1,2,...) from each
  # component, producing discontiguous indices (1,2,1,2) for `which =
  # "all"`. Reset to sequential so downstream code doesn't trip on the
  # duplicates (e.g. broom's `bind_rows()` checks rownames uniqueness).
  rownames(result) <- NULL
  result
}

#' Glance at a causatr result
#'
#' @description
#' Returns a one-row data frame of model-level summaries from a
#' `causatr_result` object, compatible with the
#' [broom](https://broom.tidymodels.org/) ecosystem.
#'
#' @param x A `causatr_result` object.
#' @param ... Currently unused.
#'
#' @return A one-row data.frame with columns `estimator`, `estimand`,
#'   `contrast_type`, `ci_method`, `n`, and `n_interventions`.
#'
#' @examples
#' \dontrun{
#' result <- contrast(fit, interventions = list(a1 = static(1), a0 = static(0)))
#' glance(result)
#' }
#'
#' @seealso [coef.causatr_result()]
#' @export
glance.causatr_result <- function(x, ...) {
  data.frame(
    estimator = x$estimator,
    estimand = x$estimand,
    contrast_type = x$type,
    ci_method = x$ci_method,
    n = x$n,
    n_interventions = nrow(x$estimates),
    stringsAsFactors = FALSE
  )
}
