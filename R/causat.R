#' Fit a causal model
#'
#' @description
#' Prepares the causal estimation pipeline for a given estimator. For
#' `"gcomp"`, fits the conditional outcome model E\[Y | A, L\] that will be
#' used by [contrast()] for standardisation. For `"ipw"`, fits the
#' conditional treatment density \eqn{f(A \mid L)} that the self-contained
#' density-ratio engine uses to build per-intervention Hájek weights for
#' the MSM refit inside [contrast()]. For `"matching"`, creates matched
#' sets (via `MatchIt::matchit()`) that will be used for matched
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
#' @param estimator Character. Causal estimator: `"gcomp"` (default), `"ipw"`,
#'   or `"matching"`. IPW uses a self-contained density-ratio engine
#'   (no runtime dependency on `WeightIt`); matching requires the
#'   `MatchIt` package. Note: `"matching"` is restricted to **binary
#'   point treatments** (MatchIt does not support multi-category or
#'   continuous treatments); use `"gcomp"` or `"ipw"` for those cases.
#' @param family Character or family object. The outcome model family
#'   (e.g. `"gaussian"`, `"binomial"`, `stats::quasibinomial()`). Used by
#'   all methods: passed to the outcome model for `"gcomp"`, to the
#'   per-intervention weighted MSM refit for `"ipw"`, and to the outcome
#'   model on matched data for `"matching"`.
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
#'   variable (1 = censored, 0 = uncensored). **This is a row filter, not
#'   IPCW:** rows where `censoring != 0` are excluded from model fitting,
#'   but no censoring model is fit and no inverse-probability weights are
#'   computed. For g-computation with a correctly specified outcome model,
#'   this is sufficient under MAR censoring (the regression surface is
#'   unchanged). For IPW under MAR censoring, supply pre-computed IPCW
#'   weights via `weights =`. For longitudinal data, censoring is
#'   time-varying: `C_k = 1` means the individual dropped out at time `k`;
#'   ICE backward iteration restricts each step's model to uncensored
#'   individuals (Zivich et al., 2024).
#' @param history Positive integer or `Inf`. Markov order for longitudinal
#'   models: how many lagged time points of treatment and time-varying
#'   confounders to include in each ICE outcome model. Default `1` (one lag).
#'   `Inf` includes the full history. Ignored for point treatments.
#' @param numerator A one-sided formula or `NULL`. Numerator formula for
#'   stabilized IPW weights in longitudinal models. Defaults to baseline
#'   confounders only (no time-varying confounders), which gives the standard
#'   stabilized weights. Only relevant for `estimator = "ipw"` with longitudinal
#'   data.
#' @param weights Numeric vector or `NULL`. Pre-computed observation weights
#'   (e.g. survey weights or externally computed IPCW). For `"gcomp"`,
#'   passed to `glm()`. For `"ipw"`, multiplied with the estimated propensity
#'   weights.
#' @param type Character or `NULL`. Whether the data are `"point"` (single
#'   time-point per individual) or `"longitudinal"` (repeated measures).
#'   Default `NULL` auto-detects: `"longitudinal"` when both `id` and `time`
#'   are provided, `"point"` otherwise. Passing `type = "longitudinal"`
#'   explicitly requires `id` and `time`.
#' @param model_fn Function. The fitting function for the outcome model in
#'   g-computation. Must accept `(formula, data, family, weights, ...)`.
#'   Default `stats::glm`; pass `mgcv::gam` for GAMs, `MASS::glm.nb` for
#'   negative-binomial, etc. For `"ipw"`, `model_fn` also fits the
#'   placeholder `Y ~ A` display model and is the default propensity
#'   fitter when `propensity_model_fn` is `NULL`. Ignored for
#'   `"matching"`.
#' @param propensity_model_fn Function or `NULL`. IPW only. The fitting
#'   function for the conditional treatment density \eqn{f(A \mid L)}
#'   used to build density-ratio weights. Must accept the same
#'   `(formula, data, family, weights, ...)` signature. `NULL`
#'   (default) reuses `model_fn`. Pass `mgcv::gam` for a flexible
#'   propensity.
#' @param propensity_family Character or `NULL`. IPW only. Explicit
#'   treatment density family: `"poisson"` or `"negbin"` for count
#'   treatments. `NULL` (default) auto-detects from the treatment
#'   values (bernoulli / gaussian / categorical). Auto-detection never
#'   infers count — use this parameter to opt in. For `"negbin"`,
#'   `MASS::glm.nb` is auto-selected as the propensity fitter unless
#'   `propensity_model_fn` is explicitly provided. Ignored for
#'   `"gcomp"` and `"matching"`.
#' @param ... Additional arguments passed to the underlying estimation
#'   function. For `estimator = "ipw"`, dots are forwarded into the
#'   user's `propensity_model_fn` via `fit_treatment_model()` (e.g.
#'   smoothing arguments for `mgcv::gam`). For `estimator = "matching"`,
#'   dots are forwarded into `MatchIt::matchit()` (e.g.
#'   `method = "nearest"`, `ratio = 1`, `caliper = 0.2`).
#'
#' @return A `causatr_fit` object with slots:
#'   \describe{
#'     \item{`model`}{Fitted model object(s): `glm`/`gam` for `"gcomp"`;
#'       a placeholder `Y ~ A` weighted MSM for `"ipw"` (the density
#'       model lives in `details$propensity_model`); the matched-data
#'       outcome model for `"matching"`.}
#'     \item{`data`}{data.table used for fitting.}
#'     \item{`treatment`, `outcome`, `confounders`, `confounders_tv`,
#'       `family`}{Model spec.}
#'     \item{`estimator`}{`"gcomp"`, `"ipw"`, or `"matching"`.}
#'     \item{`type`}{`"point"` or `"longitudinal"`.}
#'     \item{`estimand`}{`"ATE"`, `"ATT"`, or `"ATC"`.}
#'     \item{`id`, `time`, `censoring`}{Longitudinal identifiers.}
#'     \item{`history`}{Markov order for longitudinal ICE models.}
#'     \item{`numerator`}{Numerator formula for longitudinal IPW.}
#'     \item{`weights_obj`}{Legacy slot, always `NULL`.}
#'     \item{`match_obj`}{`matchit` object (matching only).}
#'     \item{`call`}{The original call.}
#'     \item{`details`}{Internal diagnostics list.}
#'   }
#'
#' @details
#' ## G-computation (`estimator = "gcomp"`)
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
#' ## IPW (`estimator = "ipw"`)
#' Fits a conditional treatment density \eqn{f(A \mid L)} via
#' `propensity_model_fn` (defaults to `model_fn`). Each intervention
#' requested in [contrast()] builds its own density-ratio weight
#' vector and refits a weighted marginal structural model. Supports
#' binary (via HT or IPSI), continuous (via smooth shift / scale_by),
#' and dynamic-on-binary interventions. `threshold()` and continuous
#' `dynamic()` rules are routed to `estimator = "gcomp"` with a
#' pointer. Categorical treatments are not supported under IPW;
#' `estimator = "gcomp"` handles them via predict-then-average.
#'
#' **Longitudinal IPW is not supported.** `estimator = "ipw"` with `id`
#' and `time` aborts with an informative error. Use
#' `estimator = "gcomp"` (which uses ICE g-computation for longitudinal
#' data) for time-varying treatments.
#'
#' ## Matching (`estimator = "matching"`)
#' Calls `MatchIt::matchit()` to create matched sets. The estimand is
#' fixed at fitting time. Only `static()` interventions are supported
#' in [contrast()]. **Matching is binary-only**: `MatchIt` does not
#' support categorical (k > 2 levels) or continuous treatments, and
#' causatr aborts with a clear error in both cases. Use
#' `estimator = "gcomp"` or `estimator = "ipw"` for categorical /
#' continuous treatments. Longitudinal matching is also unsupported.
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
#' ## Missing data handling
#'
#' **Outcome (Y) missing.** Rows with `NA` outcome are excluded from model
#' fitting via `get_fit_rows()`. If a `censoring` column is provided, rows
#' with `censoring != 0` are also excluded. [contrast()] predicts on ALL
#' rows (including those with missing Y) and standardizes over the target
#' population. Under MCAR, this complete-case analysis is unbiased. Under
#' MAR (censoring depends on A and/or L), g-computation with a correctly
#' specified outcome model is still consistent because the regression
#' surface E\[Y | A, L\] is unchanged by the censoring mechanism (Hernán &
#' Robins, Ch. 13). For IPW under MAR, IPCW weights are needed — supply
#' them via `weights =` (built-in IPCW is planned for Phase 11). The
#' `censoring =` parameter is a **row filter**, not IPCW: no censoring
#' model is fit internally.
#'
#' **Treatment (A) missing.** `causat()` aborts if any treatment values
#' are `NA` and no `censoring` column is provided. You cannot define a
#' counterfactual without knowing the treatment value. Remove rows with
#' missing A before calling `causat()` (unbiased under MCAR), or use
#' multiple imputation (planned via `causat_mice()`).
#'
#' **Covariates (L) missing.** For `estimator = "gcomp"`, the outcome
#' model's `na.action = na.omit` drops rows with NA covariates at fit
#' time. `predict()` returns NA for rows with missing L, which
#' [contrast()] excludes from the target population. Under MCAR this is
#' unbiased. For `estimator = "ipw"`, the propensity model also drops
#' NA-covariate rows; if this creates a row-count mismatch with the
#' outcome mask, `causat()` aborts with a clear error. For
#' `estimator = "matching"`, MatchIt handles covariate NAs internally.
#' Under MAR, multiple imputation is needed (planned via
#' `causat_mice()`).
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
#' # IPW — self-contained density-ratio engine; estimand fixed at
#' # fit time. A non-default `propensity_model_fn` (e.g. `mgcv::gam`)
#' # swaps in a smooth propensity model without touching the rest of
#' # the pipeline.
#' fit_ipw <- causat(
#'   nhefs,
#'   outcome = "wt82_71",
#'   treatment = "qsmk",
#'   confounders = ~ sex + age + race + education +
#'     smokeintensity + smokeyrs + exercise + active + wt71,
#'   estimator = "ipw",
#'   estimand = "ATE"
#' )
#'
#' # Matching (requires MatchIt). `method = "nearest"` is MatchIt's own arg.
#' fit_match <- causat(
#'   nhefs,
#'   outcome = "wt82_71",
#'   treatment = "qsmk",
#'   confounders = ~ sex + age + race + education +
#'     smokeintensity + smokeyrs + exercise + active + wt71,
#'   estimator = "matching",
#'   estimand = "ATT",
#'   method = "nearest"
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
#' @seealso [contrast()], [diagnose()], [causat_survival()],
#'   [static()], [shift()], [dynamic()]
#' @export
causat <- function(
  data,
  outcome,
  treatment,
  confounders,
  confounders_tv = NULL,
  estimator = c("gcomp", "ipw", "matching"),
  family = "gaussian",
  estimand = c("ATE", "ATT", "ATC"),
  id = NULL,
  time = NULL,
  censoring = NULL,
  history = 1L,
  numerator = NULL,
  weights = NULL,
  type = NULL,
  model_fn = stats::glm,
  propensity_model_fn = NULL,
  propensity_family = NULL,
  ...
) {
  # Capture the call for later display in print/summary of the result.
  call <- match.call()
  estimator <- rlang::arg_match(estimator)
  estimand <- rlang::arg_match(estimand)

  # Auto-detect point vs longitudinal from the presence of id/time.
  # `type` lets the user force one or the other — useful for tests
  # and for the rare case of longitudinal-shaped data the user wants
  # to analyze cross-sectionally at a single time point.
  has_long <- !is.null(id) && !is.null(time)

  if (is.null(type)) {
    type <- if (has_long) "longitudinal" else "point"
  } else {
    type <- rlang::arg_match(type, values = c("point", "longitudinal"))
    if (type == "longitudinal" && !has_long) {
      rlang::abort(
        '`type = "longitudinal"` requires both `id` and `time` to be specified.',
        call = call
      )
    }
  }

  # All structural validation happens here. By the time check_causat_inputs()
  # returns, the arguments are guaranteed consistent and any missing columns
  # are surfaced with clear error messages — the downstream fit_* functions
  # don't need to re-validate.
  check_causat_inputs(
    # nolint: object_usage_linter
    data,
    outcome = outcome,
    treatment = treatment,
    confounders = confounders,
    confounders_tv = confounders_tv,
    estimator = estimator,
    estimand = estimand,
    id = id,
    time = time,
    history = history
  )

  # prepare_data() does three things:
  #   1. Coerces `data` to data.table (for fast subsetting + `:=`).
  #   2. For longitudinal data, sorts by (id, time) and materializes
  #      lag columns (`lag1_A`, `lag2_A`, ...) up to `history`.
  #   3. Validates that the person-period structure is rectangular
  #      (every id observed at every time, or consistently censored).
  # All downstream fit_* functions assume the prepared shape.
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

  # NA check on treatment values: if any are missing, user must either
  # provide a censoring column (IPCW), use mice imputation, or remove
  # incomplete cases manually. We do this AFTER prepare_data() because
  # lag materialization is what actually exposes the NAs at baseline.
  check_treatment_nas(data, treatment, censoring)

  # Validate external weights up front. Earlier versions silently
  # passed non-finite / negative weights through to the fit step,
  # where GLM's own check sometimes aborted with a cryptic message and
  # sometimes silently produced NaN estimates. Reject at the causatr
  # boundary with a specific error so users know which call site is
  # the problem.
  check_weights(weights, nrow(data))

  # Refuse `na.action = na.exclude` forwarded through `...`. The
  # variance engine assumes `length(residuals(m, "working")) ==
  # nrow(model.matrix(m))`, which `na.exclude` violates by padding with
  # NAs — recycling then silently corrupts the IF and sandwich SEs.
  # See check_dots_na_action() for the full rationale.
  check_dots_na_action(..., call = call)

  # Dispatch to the estimator-specific fitter. Each returns a
  # `causatr_fit` with the same S3 class and slot structure, which
  # contrast() and diagnose() then consume uniformly.
  switch(
    estimator,
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
      family,
      estimand,
      type,
      history,
      numerator,
      weights,
      model_fn,
      propensity_model_fn,
      propensity_family,
      call,
      ...
    ),
    matching = fit_matching(
      data,
      outcome,
      treatment,
      confounders,
      family,
      estimand,
      type,
      weights,
      call,
      ...
    ),
    # Defensive default: unreachable under normal use because
    # rlang::arg_match(estimator) above restricts `estimator` to the
    # allowed set. Kept so a future refactor that loosens arg_match
    # cannot silently return NULL from causat().
    rlang::abort(c(
      paste0("Unknown `estimator` '", estimator, "'."),
      i = "Must be one of: 'gcomp', 'ipw', 'matching'."
    ))
  )
}
