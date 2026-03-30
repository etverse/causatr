#' Compute causal contrasts from a fitted model
#'
#' @description
#' Standardises outcome model predictions under each named intervention and
#' reports pairwise causal contrasts (differences, ratios, or odds ratios)
#' with uncertainty estimates.
#'
#' For all three estimation methods, `contrast()` computes E\[Y^a\] by setting
#' each individual's treatment to the intervened value and averaging the fitted
#' outcome model predictions over the target population. The methods differ in
#' how the outcome model was fitted and how the variance is estimated:
#' - `"gcomp"`: standard `glm`/`gam` on the full data; sandwich SE via stacked
#'   estimating equations.
#' - `"ipw"`: `glm_weightit()` fit weighted for the target estimand (ATE/ATT);
#'   M-estimation sandwich SE that accounts for weight estimation uncertainty.
#' - `"matching"`: `glm()` on the matched sample with match weights; SE via
#'   cluster-robust sandwich (`sandwich::vcovCL(vcov = ~subclass)`) to account
#'   for pair membership.
#'
#' @param fit A `causatr_fit` object returned by [causat()].
#' @param interventions A named list of `causatr_intervention` objects created
#'   by [static()], [shift()], [dynamic()], [scale()], [threshold()], or
#'   [ipsi()]. At least one element is required.
#'
#'   **Note:** Non-static interventions (`shift()`, `scale()`, `threshold()`,
#'   `dynamic()`, `ipsi()`) are only supported for `method = "gcomp"`.
#'   For IPW and matching, the weights/matched sets were estimated under the
#'   observed treatment distribution; applying a different regime requires
#'   re-calling [causat()] with updated data. This restriction applies
#'   regardless of treatment type (binary, categorical, or continuous).
#' @param type Character. The contrast scale: `"difference"` (default),
#'   `"ratio"`, or `"or"` (odds ratio). All pairwise contrasts are reported.
#' @param ci_method Character. The variance/CI method: `"sandwich"` (default,
#'   valid for GLM-based fits), `"bootstrap"`, or `"delta"`.
#' @param n_boot Integer. Number of bootstrap replications when
#'   `ci_method = "bootstrap"`. Default `500`.
#' @param conf_level Numeric. Confidence level for intervals. Default `0.95`.
#' @param by Character or `NULL`. Name of a variable to stratify estimates by
#'   (effect modification). If provided, E\[Y^a\] is computed within each
#'   level of `by`.
#'
#' @return A `causatr_result` object with slots:
#'   \describe{
#'     \item{`estimates`}{data.table with one row per intervention:
#'       `intervention`, `estimate`, `se`, `ci_lower`, `ci_upper`.}
#'     \item{`contrasts`}{data.table with one row per pairwise comparison:
#'       `comparison`, `estimate`, `se`, `ci_lower`, `ci_upper`.}
#'     \item{`type`}{Contrast scale.}
#'     \item{`ci_method`}{Inference method.}
#'     \item{`interventions`}{The intervention list.}
#'     \item{`n`}{Sample size.}
#'     \item{`method`}{Estimation method from the fit.}
#'     \item{`vcov`}{Full variance-covariance matrix for all E\[Y^a\].}
#'     \item{`call`}{The original call.}
#'   }
#'
#' @details
#' ## Standardisation
#' For each intervention `a`, `contrast()` evaluates the intervention function
#' on the observed data to obtain the intervened treatment vector `a(Lᵢ)`,
#' then computes:
#' ```
#' E[Y^a] = (1/n) Σᵢ Ê[Y | A = a(Lᵢ), Lᵢ]
#' ```
#' This operation is the same for all three methods — the difference is in
#' the outcome model stored in `fit$model` (standard `glm`/`gam` for gcomp;
#' `glm_weightit()` for IPW; `glm()` on the matched sample for matching)
#' and how uncertainty propagates to the variance estimate.
#'
#' ## Supported treatment and intervention types by method
#' | Method | Treatment types | Intervention types |
#' |---|---|---|
#' | `"gcomp"` | binary, categorical, continuous | all |
#' | `"ipw"` | binary, categorical, continuous (WeightIt) | `static()` only |
#' | `"matching"` | binary, categorical, continuous | `static()` only |
#'
#' ## Variance estimation
#' - `"sandwich"`: Stacked estimating equations propagate outcome model
#'   uncertainty through to the marginal mean. For IPW fits, M-estimation
#'   includes the propensity score estimating equations so weight uncertainty
#'   is accounted for. For matching fits, cluster-robust standard errors
#'   cluster on matched pair (`subclass`).
#' - `"bootstrap"`: Resamples individuals (entire pipeline refitted `n_boot`
#'   times). Respects cluster structure for clustered data. Parallelisable
#'   via a `future` plan.
#' - `"delta"`: Explicit delta method for contrast transformations (applied
#'   internally for ratio/OR contrasts when `ci_method = "sandwich"`).
#'
#' @references
#' Hernán MA, Robins JM (2025). *Causal Inference: What If*. Chapman &
#' Hall/CRC. Chapters 12–13.
#'
#' Greifer N (2024). WeightIt: Weighting Methods for Covariate Balancing.
#' \url{https://ngreifer.github.io/WeightIt/}
#'
#' Imai K, King G, Stuart EA (2011). Misunderstandings between experimentalists
#' and observationalists about causal inference. *Journal of the Royal
#' Statistical Society* Series A 171:481–502.
#'
#' Zivich PN, Ross RK, Shook-Sa BE, Cole SR, Edwards JK (2024). Empirical
#' sandwich variance estimator for iterated conditional expectation
#' g-computation. *Statistics in Medicine* 43:5562–5572.
#'
#' @examples
#' \dontrun{
#' data("nhefs", package = "causatr")
#' fit <- causat(nhefs, outcome = "wt82_71", treatment = "qsmk",
#'               confounders = ~ sex + age + wt71)
#'
#' # Mean difference with sandwich SE
#' result <- contrast(fit,
#'   interventions = list(quit = static(1), continue = static(0)),
#'   type = "difference"
#' )
#' print(result)
#'
#' # Risk ratio with bootstrap CI
#' contrast(fit,
#'   interventions = list(quit = static(1), continue = static(0)),
#'   type = "ratio",
#'   ci_method = "bootstrap",
#'   n_boot = 500
#' )
#' }
#'
#' @seealso [causat()], [static()], [shift()], [dynamic()],
#'   [coef.causatr_result()], [confint.causatr_result()]
#' @export
contrast <- function(
  fit,
  interventions,
  type = c("difference", "ratio", "or"),
  ci_method = c("sandwich", "bootstrap", "delta"),
  n_boot = 500L,
  conf_level = 0.95,
  by = NULL
) {
  call <- match.call()

  if (!inherits(fit, "causatr_fit")) {
    rlang::abort("`fit` must be a `causatr_fit` object returned by `causat()`.")
  }

  type <- rlang::arg_match(type)
  ci_method <- rlang::arg_match(ci_method)
  check_intervention_list(interventions)
  check_interventions_compat(fit$method, interventions)

  compute_contrast(
    fit,
    interventions,
    type,
    ci_method,
    n_boot,
    conf_level,
    by,
    call
  )
}

#' @noRd
check_interventions_compat <- function(
  method,
  interventions,
  call = rlang::caller_env()
) {
  if (method %in% c("ipw", "matching")) {
    non_static <- vapply(
      interventions,
      function(iv) iv$type != "static",
      logical(1)
    )
    if (any(non_static)) {
      rlang::abort(
        paste0(
          "Non-static interventions (shift, dynamic, scale, threshold, ipsi) ",
          "are not supported for method = '",
          method,
          "'. ",
          "The weights/matched sets were estimated under the original ",
          "treatment regime and are not valid under a different intervention. ",
          "Use method = 'gcomp' instead."
        ),
        call = call
      )
    }
  }
}

#' @noRd
compute_contrast <- function(
  fit,
  interventions,
  type,
  ci_method,
  n_boot,
  conf_level,
  by,
  call
) {
  rlang::abort("contrast() is not yet implemented.", .call = FALSE)
}
