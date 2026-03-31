#' Compute causal contrasts from a fitted model
#'
#' @description
#' Standardises outcome model predictions under each named intervention and
#' reports pairwise causal contrasts (differences, ratios, or odds ratios)
#' with uncertainty estimates.
#'
#' For all three estimation methods, `contrast()` computes E\[Y^a\] by setting
#' each individual's treatment to the intervened value and averaging the fitted
#' outcome model predictions over the target population. The target population
#' is controlled by `estimand` (or `subset` for subgroup effects). The methods
#' differ in how the outcome model was fitted and how the variance is estimated:
#' - `"gcomp"`: standard `glm`/`gam` on the full data; sandwich SE via stacked
#'   estimating equations.
#' - `"ipw"`: `glm_weightit()` fit weighted for the target estimand (ATE/ATT);
#'   M-estimation sandwich SE that accounts for weight estimation uncertainty.
#' - `"matching"`: `glm()` on the matched sample with match weights; SE via
#'   cluster-robust sandwich (`sandwich::vcovCL(vcov = ~subclass)`) to account
#'   for pair membership.
#'
#' @param fit A `causatr_fit` object returned by [causat()].
#' @param interventions A named list of interventions. Each element must be
#'   one of:
#'   - A `causatr_intervention` object created by [static()], [shift()],
#'     [dynamic()], [scale()], [threshold()], or [ipsi()].
#'   - `NULL`, meaning the natural course (observed treatment values are used
#'     as-is). This is the reference for modified treatment policies. Note: for
#'     binary treatments, the natural-course mean is the observed mean, not a
#'     counterfactual — a warning is issued.
#'   - A named list of `causatr_intervention` objects, one per treatment
#'     variable, for multivariate (joint) treatments. Each sub-list must name
#'     every treatment variable specified in `causat()`.
#'
#'   **Note:** Non-static interventions (`shift()`, `scale()`, `threshold()`,
#'   `dynamic()`, `ipsi()`) are only supported for `method = "gcomp"`.
#'   For IPW and matching, the weights/matched sets were estimated under the
#'   observed treatment distribution; applying a different regime requires
#'   re-calling [causat()] with updated data.
#' @param type Character. The contrast scale: `"difference"` (default),
#'   `"ratio"`, or `"or"` (odds ratio). All pairwise contrasts are reported.
#' @param estimand Character or `NULL`. The target estimand: `"ATE"`,
#'   `"ATT"`, or `"ATC"`. For `method = "gcomp"`, overrides the estimand
#'   stored in `fit` (allowing one fit to produce multiple estimands). For
#'   `method = "ipw"` or `"matching"`, must match the estimand used at fitting
#'   time — changing it aborts with an informative message. If `NULL`,
#'   defaults to `fit$estimand`. Mutually exclusive with `subset`.
#' @param subset A quoted expression (`quote(...)`) defining a subgroup to
#'   average over instead of an estimand. Evaluated in the context of the
#'   fitted data. Works for any treatment type. For example,
#'   `subset = quote(age > 50)` or `subset = quote(A == 1)` (the latter
#'   is equivalent to `estimand = "ATT"` for binary treatments). Mutually
#'   exclusive with `estimand`.
#' @param reference Character or `NULL`. Name of the reference intervention
#'   (the denominator/subtracted value for pairwise contrasts). Default: the
#'   first intervention in the list. Only relevant when `type = "difference"`
#'   or `"ratio"` and there are more than two interventions. For categorical
#'   treatments, use this to specify the reference level.
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
#'     \item{`estimand`}{`"ATE"`, `"ATT"`, `"ATC"`, or `"subset"`.}
#'     \item{`ci_method`}{Inference method.}
#'     \item{`reference`}{Name of the reference intervention.}
#'     \item{`interventions`}{The intervention list.}
#'     \item{`n`}{Number of individuals averaged over.}
#'     \item{`method`}{Estimation method from the fit.}
#'     \item{`vcov`}{Full variance-covariance matrix for all E\[Y^a\].}
#'     \item{`call`}{The original call.}
#'   }
#'
#' @details
#' ## Estimands and standardisation
#' For each intervention `a`, `contrast()` evaluates the intervention function
#' on the target rows to obtain the intervened treatment vector `a(Lᵢ)`, then
#' computes:
#' ```
#' E[Y^a] = (1/|S|) Σᵢ∈S Ê[Y | A = a(Lᵢ), Lᵢ]
#' ```
#' where `S` is the set of rows determined by the estimand:
#' - `"ATE"`: all rows (full population)
#' - `"ATT"`: rows where `A == 1` (observed treated)
#' - `"ATC"`: rows where `A == 0` (observed controls)
#' - `subset`: rows satisfying the quoted expression
#'
#' For `"gcomp"`, one fit supports multiple estimands because the outcome
#' model is the same — only the rows averaged over change. For `"ipw"` and
#' `"matching"`, the estimand is baked into the weights/matching and cannot
#' be changed after fitting.
#'
#' ## Estimand applicability
#' | Treatment type | ATE | ATT/ATC | subset |
#' |---|---|---|---|
#' | Binary point | Yes | Yes | Yes |
#' | Continuous point | Yes | No (abort) | Yes |
#' | Categorical point | Yes | No (abort) | Yes |
#' | Multivariate point | Yes | No (abort) | Yes |
#' | Longitudinal | Yes | No (abort) | Yes |
#'
#' ## Treatment types and intervention support
#' | Method | Treatment types | Intervention types |
#' |---|---|---|
#' | `"gcomp"` | binary, categorical, continuous, multivariate | all |
#' | `"ipw"` | binary, categorical, continuous (WeightIt) | `static()` only |
#' | `"matching"` | binary, categorical, continuous | `static()` only |
#'
#' ## Variance estimation
#' - `"sandwich"`: Stacked estimating equations propagate outcome model
#'   uncertainty through to the marginal mean.
#' - `"bootstrap"`: Resamples individuals (entire pipeline refitted `n_boot`
#'   times). Respects cluster structure for longitudinal data.
#' - `"delta"`: Explicit delta method for ratio/OR contrasts (applied
#'   internally when `ci_method = "sandwich"`).
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
#' # ATE: mean difference with sandwich SE
#' result <- contrast(fit,
#'   interventions = list(quit = static(1), continue = static(0)),
#'   type = "difference"
#' )
#'
#' # ATT: override estimand in contrast() (gcomp only)
#' result_att <- contrast(fit,
#'   interventions = list(quit = static(1), continue = static(0)),
#'   estimand = "ATT"
#' )
#'
#' # Subgroup effect: age > 50
#' result_sub <- contrast(fit,
#'   interventions = list(quit = static(1), continue = static(0)),
#'   subset = quote(age > 50)
#' )
#'
#' # Continuous treatment: shift with NULL reference
#' fit_cont <- causat(nhefs, outcome = "wt82_71",
#'                    treatment = "smokeintensity",
#'                    confounders = ~ sex + age + wt71)
#' contrast(fit_cont,
#'   interventions = list(reduce10 = shift(-10), observed = NULL),
#'   type = "difference"
#' )
#'
#' # Categorical treatment: three arms, reference = "radio"
#' contrast(fit_cat,
#'   interventions = list(chemo = static("A"), radio = static("B"),
#'                        combo = static("C")),
#'   type = "difference",
#'   reference = "radio"
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
  estimand = NULL,
  subset = NULL,
  reference = NULL,
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

  if (!is.null(estimand) && !is.null(subset)) {
    rlang::abort("Specify either 'estimand' or 'subset', not both.")
  }

  if (!is.null(estimand)) {
    estimand <- rlang::arg_match(estimand, c("ATE", "ATT", "ATC"))
    check_estimand_compat(estimand, fit$method, fit$estimand)
    check_estimand_treatment_compat(
      estimand,
      fit$treatment,
      fit$type
    )
  }

  check_intervention_list(interventions)
  check_interventions_compat(fit$method, interventions)

  if (!is.null(reference) && !reference %in% names(interventions)) {
    rlang::abort(
      paste0(
        "`reference` ('",
        reference,
        "') must be the name of one of the interventions."
      )
    )
  }

  compute_contrast(
    fit,
    interventions,
    type,
    estimand,
    subset,
    reference,
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
      function(iv) {
        if (is.null(iv)) return(FALSE)
        if (is.list(iv) && !inherits(iv, "causatr_intervention")) {
          return(any(vapply(iv, function(sub) sub$type != "static", logical(1))))
        }
        iv$type != "static"
      },
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
  estimand,
  subset,
  reference,
  ci_method,
  n_boot,
  conf_level,
  by,
  call
) {
  rlang::abort("contrast() is not yet implemented.", .call = FALSE)
}
