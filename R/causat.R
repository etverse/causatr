#' Fit a causal model
#'
#' @description
#' Prepares the causal estimation pipeline for a given method. For
#' `"gcomp"`, fits the conditional outcome model E\[Y | A, L\] that will be
#' used by [contrast()] for standardisation. For `"ipw"`, estimates
#' propensity-score-based weights (via `WeightIt::weightit()`) that will be
#' used for weighted estimation in [contrast()]. For `"matching"`, creates
#' matched sets (via `MatchIt::matchit()`) that will be used for matched
#' estimation in [contrast()].
#'
#' For longitudinal data (`id` and `time` provided), `"gcomp"` uses ICE
#' g-computation (Zivich et al., 2024): outcome models are fitted at each
#' time point via backward iteration.
#'
#' @param data A data frame or data.table.
#' @param outcome Character. Name of the outcome variable.
#' @param treatment Character. Name of the treatment variable.
#' @param confounders A one-sided formula specifying confounders,
#'   e.g. `~ L1 + L2`. Interactions and transformations are allowed,
#'   e.g. `~ L1 * L2 + splines::ns(age, 4)`.
#' @param method Character. Estimation method: `"gcomp"` (default), `"ipw"`,
#'   or `"matching"`. IPW requires the `WeightIt` package; matching requires
#'   the `MatchIt` package.
#' @param family Character or family object. The outcome model family for
#'   `"gcomp"` (e.g. `"gaussian"`, `"binomial"`). Passed to `glm()` or
#'   `mgcv::gam()`. Ignored for `"ipw"` and `"matching"`.
#' @param id Character or `NULL`. Name of the individual ID variable. Must be
#'   provided together with `time` for longitudinal data.
#' @param time Character or `NULL`. Name of the time variable. Must be
#'   provided together with `id` for longitudinal data.
#' @param censoring Character or `NULL`. Name of the censoring indicator
#'   variable (1 = censored, 0 = uncensored).
#' @param weights Numeric vector or `NULL`. Pre-computed observation weights
#'   (e.g. survey weights or externally computed IPCW). For `"gcomp"`,
#'   passed to `glm()`. For `"ipw"`, multiplied with the estimated propensity
#'   weights.
#' @param ... Additional arguments passed to the underlying estimation
#'   function: `WeightIt::weightit()` for `method = "ipw"` (e.g.
#'   `estimand = "ATE"`, `method = "glm"`); `MatchIt::matchit()` for
#'   `method = "matching"` (e.g. `method = "nearest"`, `ratio = 1`).
#'
#' @return A `causatr_fit` object with slots:
#'   \describe{
#'     \item{`model`}{Fitted model object(s): `glm`/`gam` for `"gcomp"`;
#'       `NULL` for `"ipw"` and `"matching"` (weights/matched data are
#'       stored in `weights_obj` / `match_obj` instead).}
#'     \item{`data`}{data.table used for fitting.}
#'     \item{`treatment`, `outcome`, `confounders`, `family`}{Model spec.}
#'     \item{`method`}{`"gcomp"`, `"ipw"`, or `"matching"`.}
#'     \item{`type`}{`"point"` or `"longitudinal"`.}
#'     \item{`id`, `time`, `censoring`}{Longitudinal identifiers.}
#'     \item{`weights_obj`}{`weightit` object (IPW only).}
#'     \item{`match_obj`}{`matchit` object (matching only).}
#'     \item{`call`}{The original call.}
#'     \item{`details`}{Internal diagnostics list.}
#'   }
#'
#' @details
#' ## G-computation (`method = "gcomp"`)
#' Fits E\[Y | A, L\] using `glm()` (or `mgcv::gam()` if the formula
#' contains `s()` or `te()` terms). [contrast()] standardises over the
#' observed covariate distribution to obtain E\[Y^a\] = (1/n) Σᵢ Ê\[Y | A=a,
#' Lᵢ\] under each intervention. Supports any intervention type: static,
#' dynamic, shift, threshold, scale, IPSI.
#'
#' For longitudinal data, uses ICE g-computation (Zivich et al., 2024):
#' outcome models are fitted one per time point via backward iteration.
#' The resulting pseudo-outcomes carry the intervention-specific causal
#' information forward.
#'
#' ## IPW (`method = "ipw"`)
#' Calls `WeightIt::weightit()` to estimate propensity-score-based weights
#' for the desired `estimand` (ATE by default; pass `estimand = "ATT"` via
#' `...`), then fits a weighted outcome model via `WeightIt::glm_weightit()`.
#' [contrast()] standardises predictions from that model under each
#' intervention. Variance via M-estimation accounts for weight estimation
#' uncertainty. Supports binary, categorical, and continuous treatments
#' (anything `WeightIt` supports). Only `static()` interventions are
#' supported, because the weights were estimated under the observed treatment
#' distribution; re-estimating under a different regime (shift, dynamic, etc.)
#' requires re-calling `causat()`.
#'
#' ## Matching (`method = "matching"`)
#' Calls `MatchIt::matchit()` to create matched sets (binary, categorical,
#' and continuous treatments are supported — see `MatchIt` documentation),
#' then fits the outcome model on the matched sample. [contrast()] standardises
#' predictions from that model. Variance via cluster-robust sandwich SE on
#' matched pairs. Only `static()` interventions are supported for the same
#' reason as IPW.
#'
#' ## Identifiability assumptions
#' All methods assume: (1) exchangeability (no unmeasured confounding given
#' L), (2) positivity (every individual has positive probability of each
#' treatment value given L), (3) consistency (the observed outcome under the
#' observed treatment equals the potential outcome). Positivity is checked
#' automatically and a warning is issued if near-violations are detected.
#'
#' @references
#' Hernán MA, Robins JM (2025). *Causal Inference: What If*. Chapman &
#' Hall/CRC. Chapters 12 (IPW), 13 (g-formula), 15 (matching), 21 (ICE).
#'
#' Zivich PN, Ross RK, Shook-Sa BE, Cole SR, Edwards JK (2024). Empirical
#' sandwich variance estimator for iterated conditional expectation
#' g-computation. *Statistics in Medicine* 43:5562–5572.
#'
#' @examples
#' \dontrun{
#' data("nhefs", package = "causatr")
#'
#' # Point treatment, g-computation (default)
#' fit <- causat(
#'   nhefs,
#'   outcome = "wt82_71",
#'   treatment = "qsmk",
#'   confounders = ~ sex + age + race + education +
#'     smokeintensity + smokeyrs + exercise + active + wt71
#' )
#'
#' # IPW (requires WeightIt)
#' fit_ipw <- causat(
#'   nhefs,
#'   outcome = "wt82_71",
#'   treatment = "qsmk",
#'   confounders = ~ sex + age + race + education +
#'     smokeintensity + smokeyrs + exercise + active + wt71,
#'   method = "ipw",
#'   estimand = "ATE"
#' )
#'
#' # Matching (requires MatchIt)
#' fit_match <- causat(
#'   nhefs,
#'   outcome = "wt82_71",
#'   treatment = "qsmk",
#'   confounders = ~ sex + age + race + education +
#'     smokeintensity + smokeyrs + exercise + active + wt71,
#'   method = "matching",
#'   match_method = "nearest"
#' )
#' }
#'
#' @seealso [contrast()], [diagnose()], [causat_survival()],
#'   [static()], [shift()], [dynamic()]
#' @export
causat <- function(
  data,
  outcome,
  treatment,
  confounders,
  method = c("gcomp", "ipw", "matching"),
  family = "gaussian",
  id = NULL,
  time = NULL,
  censoring = NULL,
  weights = NULL,
  ...
) {
  call <- match.call()
  method <- rlang::arg_match(method)

  check_causat_inputs(
    data,
    outcome = outcome,
    treatment = treatment,
    confounders = confounders,
    method = method,
    id = id,
    time = time
  )

  data <- prepare_data(
    data,
    outcome = outcome,
    treatment = treatment,
    confounders = confounders,
    id = id,
    time = time,
    censoring = censoring
  )

  type <- if (!is.null(id) && !is.null(time)) "longitudinal" else "point"

  switch(
    method,
    gcomp = fit_gcomp(
      data,
      outcome,
      treatment,
      confounders,
      family,
      type,
      weights,
      call,
      ...
    ),
    ipw = fit_ipw(
      data,
      outcome,
      treatment,
      confounders,
      type,
      weights,
      call,
      ...
    ),
    matching = fit_matching(
      data,
      outcome,
      treatment,
      confounders,
      type,
      weights,
      call,
      ...
    )
  )
}
