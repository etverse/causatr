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
#' @param treatment Character scalar or character vector. Name(s) of the
#'   treatment variable(s). Pass a character vector for multivariate
#'   (joint) treatments, e.g. `treatment = c("A1", "A2")`.
#' @param confounders A one-sided formula specifying baseline (time-invariant)
#'   confounders, e.g. `~ L1 + L2`. Interactions and transformations are
#'   allowed, e.g. `~ L1 * L2 + splines::ns(age, 4)`. For longitudinal data,
#'   these confounders are constant within each individual (measured at
#'   baseline) and enter every time-step outcome model.
#' @param confounders_tv A one-sided formula or `NULL`. Time-varying
#'   confounders for longitudinal data (e.g. `~ CD4 + viral_load`). These
#'   change over time within individuals and enter the outcome model at each
#'   ICE step alongside their lagged values (controlled by `history`). Ignored
#'   for point treatments. If `NULL`, no time-varying confounders are used.
#' @param method Character. Estimation method: `"gcomp"` (default), `"ipw"`,
#'   or `"matching"`. IPW requires the `WeightIt` package; matching requires
#'   the `MatchIt` package.
#' @param family Character or family object. The outcome model family for
#'   `"gcomp"` (e.g. `"gaussian"`, `"binomial"`). Passed to `glm()` or
#'   `mgcv::gam()`. Ignored for `"ipw"` and `"matching"`.
#' @param estimand Character. The target estimand: `"ATE"` (default),
#'   `"ATT"`, or `"ATC"`. `"ATT"` and `"ATC"` are only defined for binary
#'   point treatments. For `"gcomp"`, the estimand stored here is used as the
#'   default in [contrast()] but can be overridden there. For `"ipw"` and
#'   `"matching"`, the estimand is fixed at fitting time because it determines
#'   the weights or match direction; it cannot be changed in [contrast()].
#' @param id Character or `NULL`. Name of the individual ID variable. Must be
#'   provided together with `time` for longitudinal data.
#' @param time Character or `NULL`. Name of the time variable. Must be
#'   provided together with `id` for longitudinal data.
#' @param censoring Character or `NULL`. Name of the censoring indicator
#'   variable (1 = censored, 0 = uncensored). For longitudinal data,
#'   censoring is time-varying: `C_k = 1` means the individual dropped out
#'   at time `k`.
#' @param history Positive integer or `Inf`. Markov order for longitudinal
#'   models: how many lagged time points of treatment and time-varying
#'   confounders to include in each ICE outcome model. Default `1` (one lag).
#'   `Inf` includes the full history. Ignored for point treatments.
#' @param numerator A one-sided formula or `NULL`. Numerator formula for
#'   stabilized IPW weights in longitudinal models. Defaults to baseline
#'   confounders only (no time-varying confounders), which gives the standard
#'   stabilized weights. Only relevant for `method = "ipw"` with longitudinal
#'   data.
#' @param weights Numeric vector or `NULL`. Pre-computed observation weights
#'   (e.g. survey weights or externally computed IPCW). For `"gcomp"`,
#'   passed to `glm()`. For `"ipw"`, multiplied with the estimated propensity
#'   weights.
#' @param model_fn Function. The fitting function for the outcome model in
#'   g-computation. Must accept `(formula, data, family, weights, ...)`.
#'   Default `stats::glm`; pass `mgcv::gam` for GAMs, `MASS::glm.nb` for
#'   negative-binomial, etc. Ignored for `"ipw"` and `"matching"`.
#' @param ... Additional arguments passed to the underlying estimation
#'   function: `WeightIt::weightit()` for `method = "ipw"` (e.g.
#'   `method = "glm"`); `MatchIt::matchit()` for `method = "matching"` (e.g.
#'   `method = "nearest"`, `ratio = 1`).
#'
#' @return A `causatr_fit` object with slots:
#'   \describe{
#'     \item{`model`}{Fitted model object(s): `glm`/`gam` for `"gcomp"`;
#'       `NULL` for `"ipw"` and `"matching"` (weights/matched data are
#'       stored in `weights_obj` / `match_obj` instead).}
#'     \item{`data`}{data.table used for fitting.}
#'     \item{`treatment`, `outcome`, `confounders`, `confounders_tv`,
#'       `family`}{Model spec.}
#'     \item{`method`}{`"gcomp"`, `"ipw"`, or `"matching"`.}
#'     \item{`type`}{`"point"` or `"longitudinal"`.}
#'     \item{`estimand`}{`"ATE"`, `"ATT"`, or `"ATC"`.}
#'     \item{`id`, `time`, `censoring`}{Longitudinal identifiers.}
#'     \item{`history`}{Markov order for longitudinal ICE models.}
#'     \item{`numerator`}{Numerator formula for longitudinal IPW.}
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
#' target population (controlled by `estimand`) to obtain E\[Y^a\] under each
#' intervention. For `"gcomp"`, the estimand can be overridden in
#' [contrast()] — fit once, contrast with multiple estimands.
#'
#' For longitudinal data, uses ICE g-computation (Zivich et al., 2024):
#' outcome models are fitted one per time point via backward iteration,
#' conditioning on baseline confounders (`confounders`), time-varying
#' confounders (`confounders_tv`), and their lags up to `history` steps.
#'
#' ## IPW (`method = "ipw"`)
#' Calls `WeightIt::weightit()` to estimate propensity-score-based weights
#' for the desired `estimand`, then fits a weighted outcome model via
#' `WeightIt::glm_weightit()`. The estimand is fixed at fitting time.
#' Supports binary, categorical, and continuous treatments (anything
#' `WeightIt` supports). Only `static()` interventions are supported
#' in [contrast()], because the weights were estimated under the observed
#' treatment distribution.
#'
#' For longitudinal data, calls `WeightIt::weightitMSM()`. The denominator
#' model at each time `k` includes baseline confounders, concurrent
#' time-varying confounders, and lagged treatment. The numerator model
#' includes only baseline confounders and lagged treatment (standard
#' stabilized weights), unless overridden via `numerator`.
#'
#' ## Matching (`method = "matching"`)
#' Calls `MatchIt::matchit()` to create matched sets. The estimand is
#' fixed at fitting time. Only `static()` interventions are supported
#' in [contrast()].
#'
#' ## Estimands
#' | Estimand | Population averaged over | Applicability |
#' |---|---|---|
#' | `"ATE"` | All individuals | Always |
#' | `"ATT"` | Observed treated (`A = 1`) | Binary point treatment only |
#' | `"ATC"` | Observed controls (`A = 0`) | Binary point treatment only |
#'
#' For continuous, categorical, multivariate, or longitudinal treatments,
#' use `estimand = "ATE"` and pass a `subset` expression to [contrast()]
#' for subgroup effects.
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
#' # ATT via g-computation (override estimand in contrast())
#' fit_att <- causat(
#'   nhefs,
#'   outcome = "wt82_71",
#'   treatment = "qsmk",
#'   confounders = ~ sex + age + race + education +
#'     smokeintensity + smokeyrs + exercise + active + wt71,
#'   estimand = "ATT"
#' )
#'
#' # IPW (requires WeightIt) — estimand fixed at fit time
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
#'   estimand = "ATT"
#' )
#'
#' # Longitudinal with time-varying confounders
#' fit_long <- causat(
#'   data_long,
#'   outcome = "Y",
#'   treatment = "A",
#'   confounders = ~ sex + race + baseline_age,
#'   confounders_tv = ~ CD4 + viral_load,
#'   id = "id",
#'   time = "time",
#'   history = 1L
#' )
#'
#' # Multivariate treatment
#' fit_multi <- causat(
#'   data,
#'   outcome = "Y",
#'   treatment = c("A1", "A2"),
#'   confounders = ~ L1 + L2
#' )
#' }
#'
#' @seealso [contrast()], [diagnose()], [causat_survival()], [causat_mice()],
#'   [static()], [shift()], [dynamic()]
#' @export
causat <- function(
  data,
  outcome,
  treatment,
  confounders,
  confounders_tv = NULL,
  method = c("gcomp", "ipw", "matching"),
  family = "gaussian",
  estimand = c("ATE", "ATT", "ATC"),
  id = NULL,
  time = NULL,
  censoring = NULL,
  history = 1L,
  numerator = NULL,
  weights = NULL,
  model_fn = stats::glm,
  ...
) {
  call <- match.call()
  method <- rlang::arg_match(method)
  estimand <- rlang::arg_match(estimand)

  type <- if (!is.null(id) && !is.null(time)) "longitudinal" else "point"

  check_causat_inputs(
    # nolint: object_usage_linter
    data,
    outcome = outcome,
    treatment = treatment,
    confounders = confounders,
    confounders_tv = confounders_tv,
    method = method,
    estimand = estimand,
    id = id,
    time = time,
    history = history
  )

  data <- prepare_data(
    # nolint: object_usage_linter
    data,
    outcome = outcome,
    treatment = treatment,
    confounders = confounders,
    confounders_tv = confounders_tv,
    id = id,
    time = time,
    censoring = censoring,
    history = history
  )

  check_treatment_nas(data, treatment, censoring)

  switch(
    method,
    gcomp = fit_gcomp(
      data,
      outcome,
      treatment,
      confounders,
      confounders_tv,
      family,
      estimand,
      type,
      history,
      censoring,
      weights,
      model_fn,
      id,
      time,
      call,
      ...
    ),
    ipw = fit_ipw(
      data,
      outcome,
      treatment,
      confounders,
      confounders_tv,
      estimand,
      type,
      history,
      numerator,
      weights,
      call,
      ...
    ),
    matching = fit_matching(
      data,
      outcome,
      treatment,
      confounders,
      estimand,
      type,
      weights,
      call,
      ...
    )
  )
}
