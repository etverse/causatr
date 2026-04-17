#' Fit causal models across multiply-imputed datasets (not implemented)
#'
#' @description
#' **Not implemented.** This function is an internal placeholder for a
#' multiple-imputation workflow. Calling it always errors. Handle
#' missing data with a complete-case analysis or impute upstream and
#' call [causat()] on a single completed dataset.
#'
#' @section Design sketch:
#' Apply [causat()] and [contrast()] across all imputed datasets in a
#' `mids` object (from `mice::mice()`), then pool point estimates and
#' variances using Rubin's rules:
#' - **Pooled estimate**: mean of per-imputation estimates.
#' - **Total variance**: within-imputation variance + between-imputation
#'   variance + between-imputation correction.
#' - **CIs**: normal approximation on the pooled estimate.
#'
#' @param imp A `mids` object returned by `mice::mice()`.
#' @param outcome Character. Passed to [causat()].
#' @param treatment Character. Passed to [causat()].
#' @param confounders A one-sided formula. Passed to [causat()].
#' @param interventions A named list of interventions. Passed to [contrast()].
#' @param estimator Character. Passed to [causat()]. Default `"gcomp"`.
#' @param family Character or family object. Passed to [causat()].
#' @param estimand Character. Passed to [causat()].
#' @param type Character. Contrast scale. Passed to [contrast()].
#' @param ci_method Character. Variance method. Passed to [contrast()].
#' @param conf_level Numeric. Confidence level. Default `0.95`.
#' @param ... Additional arguments passed to [causat()].
#'
#' @return A `causatr_result` with pooled estimates and variances following
#'   Rubin's rules.
#'
#' @details
#' ## Rubin's rules
#' Let \eqn{\hat{Q}_i} be the estimate from imputation \eqn{i} and
#' \eqn{U_i} its variance. The pooled estimate is
#' \eqn{\bar{Q} = (1/m) \sum \hat{Q}_i}. The total variance is:
#'
#' \deqn{T = \bar{U} + B + B/m}
#'
#' where \eqn{\bar{U} = (1/m) \sum U_i} (within-imputation variance) and
#' \eqn{B = (1/(m-1)) \sum (\hat{Q}_i - \bar{Q})^2} (between-imputation
#' variance).
#'
#' Degrees of freedom follow the standard Barnard-Rubin approximation.
#'
#' @references
#' Rubin DB (1987). *Multiple Imputation for Nonresponse in Surveys*.
#' Wiley.
#'
#' van Buuren S, Groothuis-Oudshoorn K (2011). mice: Multivariate Imputation
#' by Chained Equations in R. *Journal of Statistical Software* 45(3):1-67.
#'
#' @examples
#' \dontrun{
#' library(mice)
#'
#' # Step 1: impute
#' imp <- mice(data, m = 20, method = "pmm")
#'
#' # Step 2: fit + contrast across imputations
#' result <- causat_mice(
#'   imp,
#'   outcome = "Y",
#'   treatment = "A",
#'   confounders = ~ L1 + L2,
#'   interventions = list(treat_all = static(1), treat_none = static(0)),
#'   type = "difference"
#' )
#' }
#'
#' @seealso [causat()], [contrast()]
#' @keywords internal
#' @noRd
causat_mice <- function(
  imp,
  outcome,
  treatment,
  confounders,
  interventions,
  estimator = "gcomp",
  family = "gaussian",
  estimand = "ATE",
  type = "difference",
  ci_method = "sandwich",
  conf_level = 0.95,
  ...
) {
  check_pkg("mice")

  if (!inherits(imp, "mids")) {
    rlang::abort(
      paste0(
        "`imp` must be a `mids` object returned by `mice::mice()`. ",
        "Got an object of class: ",
        paste(class(imp), collapse = ", ")
      ),
      .call = FALSE
    )
  }

  rlang::abort(
    "causat_mice() is not implemented.",
    .call = FALSE
  )
}
