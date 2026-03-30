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
#'
#' @return A `causatr_diag` object with slots:
#'   \describe{
#'     \item{`balance`}{data.table or `cobalt::bal.tab` object: covariate
#'       balance summary.}
#'     \item{`positivity`}{data.table: treatment probability by covariate
#'       strata, flagging near-violations.}
#'     \item{`weights`}{list or `NULL`: weight distribution summary
#'       (IPW only).}
#'     \item{`match_quality`}{list or `NULL`: match quality summary
#'       (matching only).}
#'     \item{`method`}{Character: the estimation method.}
#'   }
#'
#' @details
#' ## Positivity
#' For binary treatment, flags strata where the estimated propensity score
#' is below 0.025 or above 0.975 (configurable). For continuous treatment,
#' flags regions of the covariate space with very low treatment density.
#'
#' ## Balance (IPW and matching)
#' If the `cobalt` package is installed, balance is computed via
#' `cobalt::bal.tab()` on the internal `weightit` or `matchit` object. This
#' provides standardised mean differences (SMD) and KS statistics before and
#' after adjustment. If `cobalt` is not installed, a simpler summary is
#' returned.
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
diagnose <- function(fit) {
  if (!inherits(fit, "causatr_fit")) {
    rlang::abort("`fit` must be a `causatr_fit` object returned by `causat()`.")
  }
  rlang::abort("diagnose() is not yet implemented.", .call = FALSE)
}
