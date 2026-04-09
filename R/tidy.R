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
  z <- stats::qnorm((1 + conf.level) / 2)

  means_df <- data.frame(
    term = x$estimates$intervention,
    estimate = x$estimates$estimate,
    std.error = x$estimates$se,
    type = "mean",
    stringsAsFactors = FALSE
  )
  if (conf.int) {
    means_df$conf.low <- means_df$estimate - z * means_df$std.error
    means_df$conf.high <- means_df$estimate + z * means_df$std.error
  }

  contrasts_df <- data.frame(
    term = x$contrasts$comparison,
    estimate = x$contrasts$estimate,
    std.error = x$contrasts$se,
    type = "contrast",
    stringsAsFactors = FALSE
  )
  if (conf.int) {
    contrasts_df$conf.low <- x$contrasts$ci_lower
    contrasts_df$conf.high <- x$contrasts$ci_upper
  }

  has_by <- "by" %in% names(x$estimates)
  if (has_by) {
    means_df$by <- x$estimates$by
    contrasts_df$by <- x$contrasts$by
  }

  result <- switch(
    which,
    contrasts = contrasts_df,
    means = means_df,
    all = rbind(means_df, contrasts_df)
  )

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
#' @return A one-row data.frame with columns `method`, `estimand`,
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
    method = x$method,
    estimand = x$estimand,
    contrast_type = x$type,
    ci_method = x$ci_method,
    n = x$n,
    n_interventions = nrow(x$estimates),
    stringsAsFactors = FALSE
  )
}
