#' Fit a causal model
#'
#' @description
#' Prepares the causal estimation pipeline for a given estimator. For
#' `"gcomp"`, fits the conditional outcome model E\[Y | A, L\] that will be
#' used by [contrast()] for standardisation. For `"ipw"`, fits the
#' conditional treatment density \eqn{f(A \mid L)} that the self-contained
#' density-ratio engine uses to build per-intervention Hajek weights for
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
#'   baseline) and enter every time-step outcome model. **Soft-deprecated:**
#'   prefer the per-component formulas (`confounders_outcome`,
#'   `confounders_treatment`, etc.) for finer control. When `NULL` (the new
#'   default), at least the per-component formulas required by the chosen
#'   estimator must be supplied.
#' @param confounders_tv A one-sided formula or `NULL`. Time-varying
#'   confounders for longitudinal data (e.g. `~ CD4 + viral_load`). These
#'   change over time within individuals and enter the outcome model at each
#'   ICE step alongside their lagged values (controlled by `history`). Ignored
#'   for point treatments. If `NULL`, no time-varying confounders are used.
#'   **Soft-deprecated:** prefer `confounders_tv_outcome` and
#'   `confounders_tv_treatment` for finer control.
#' @param confounders_outcome One-sided formula or `NULL`. Confounders
#'   for the outcome model only. When non-`NULL`, overrides `confounders`
#'   for the outcome model (g-comp, AIPW outcome side, ICE backward
#'   iteration, matching MSM). Falls back to `confounders` when `NULL`.
#' @param confounders_treatment One-sided formula or `NULL`. Confounders
#'   for the treatment (propensity) model only. When non-`NULL`, overrides
#'   `confounders` for the propensity model (IPW, AIPW treatment side,
#'   matching distance formula). Falls back to `confounders` when `NULL`.
#' @param confounders_censoring One-sided formula or `NULL`. Confounders
#'   for the censoring model when `ipcw = TRUE`. Falls back to
#'   `confounders` when `NULL`.
#' @param confounders_sampling One-sided formula or `NULL`. Confounders
#'   for the sampling model \eqn{P(S = 1 \mid L)} when `target` is
#'   non-`NULL`. Falls back to `confounders` when `NULL`.
#' @param confounders_tv_outcome One-sided formula or `NULL`. Time-varying
#'   confounders for the outcome model only (longitudinal data). Falls back
#'   to `confounders_tv` when `NULL`.
#' @param confounders_tv_treatment One-sided formula or `NULL`. Time-varying
#'   confounders for the treatment model only (longitudinal data). Falls
#'   back to `confounders_tv` when `NULL`.
#' @param estimator Character. Causal estimator: `"gcomp"` (default), `"ipw"`,
#'   `"aipw"`, `"matching"`, or `"snm"`. IPW uses a self-contained
#'   density-ratio engine (no runtime dependency on `WeightIt`); AIPW is
#'   doubly-robust (outcome model + propensity weights); matching requires
#'   the `MatchIt` package; SNM fits a structural nested mean model via
#'   g-estimation.
#'   Note: `"matching"` is restricted to **binary point treatments**
#'   (MatchIt does not support multi-category or continuous treatments);
#'   use `"gcomp"` or `"ipw"` for those cases.
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
#'   variable (1 = censored, 0 = uncensored). Without `ipcw = TRUE`,
#'   this is a **row filter only:** rows where `censoring != 0` are
#'   excluded from model fitting but no censoring model is fit. For
#'   g-computation with a correctly specified outcome model, row
#'   filtering is sufficient under MAR censoring. For IPW under MAR
#'   censoring, set `ipcw = TRUE` or supply pre-computed IPCW weights
#'   via `weights =`. For longitudinal data, censoring is time-varying:
#'   `C_k = 1` means the individual dropped out at time `k`.
#' @param ipcw Logical. If `TRUE`, fit an internal censoring model
#'   \eqn{P(C = 0 \mid A, L)} and compute stabilized IPCW weights that
#'   are composed multiplicatively with any external `weights`. Requires
#'   `censoring` to be specified. Default `FALSE` preserves the legacy
#'   row-filter-only behaviour. For IPW under MAR censoring this is
#'   essential; for g-comp it provides doubly-robust protection.
#'
#'   **Variance note:** For point treatments the sandwich SE fully accounts
#'   for censoring-model estimation uncertainty (Channel 2 correction). For
#'   longitudinal treatments (ICE, longitudinal IPW) the Channel 2 correction
#'   for per-period censoring models is omitted, making the sandwich SE
#'   conservative. Use `ci_method = "bootstrap"` for tighter CIs in
#'   longitudinal IPCW settings.
#' @param censoring_model_fn Function or `NULL`. The fitting function
#'   for the censoring model when `ipcw = TRUE`. Must accept
#'   `(formula, data, family, weights, ...)`. Default `NULL` uses
#'   `stats::glm`. Ignored when `ipcw = FALSE`.
#' @param history Non-negative integer or `Inf`. Markov lag order for
#'   longitudinal models: how many lagged time points of treatment and
#'   time-varying confounders to include in each per-period model.
#'   Default `Inf` (full history). `0` means no lags (current-period
#'   TV covariates only). Ignored for point treatments.
#' @param numerator A one-sided formula or `NULL`. Numerator formula for
#'   stabilized IPW weights in longitudinal models. Defaults to baseline
#'   confounders only (no time-varying confounders), which gives the standard
#'   stabilized weights. Only relevant for `estimator = "ipw"` with longitudinal
#'   data.
#' @param cluster Character or `NULL`. Name of a column in `data`
#'   identifying cluster membership (e.g. site, household, PSU id).
#'   Stored in the fit and preserved through `prepare_data()` so that
#'   `contrast()` can default its own `cluster =` argument to this
#'   column, producing a cluster-robust sandwich without the user
#'   having to repeat the column name. Users can still override at
#'   contrast time by passing `cluster = ` explicitly. Forwarded to
#'   matching is not allowed (matching clusters on its own `subclass`).
#' @param weights Numeric vector, a `survey::svydesign` object, or
#'   `NULL`. Pre-computed observation weights (e.g. survey weights or
#'   externally computed IPCW). When a `survey.design` object is
#'   supplied, `stats::weights()` is used to extract the sampling
#'   weights and the first-stage cluster (PSU) is auto-propagated
#'   into the fit's `cluster` slot so the contrast-time sandwich is
#'   cluster-robust by default; override by passing an explicit
#'   `cluster =` alongside. For `"gcomp"`,
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
#'   infers count -- use this parameter to opt in. For `"negbin"`,
#'   `MASS::glm.nb` is auto-selected as the propensity fitter unless
#'   `propensity_model_fn` is explicitly provided. Ignored for
#'   `"gcomp"` and `"matching"`.
#' @param stabilize Character. Multivariate IPW only. Controls the
#'   numerator density in each per-component weight factor.
#'   \itemize{
#'     \item `"none"` (default): numerator conditions on the same
#'       covariates as the denominator. Per-component weight is
#'       \eqn{f_k(d_k^{-1}(A_k) \mid A_{1..k-1}^{\mathrm{obs}}, L) \cdot |\mathrm{Jac}| / f_k(A_k \mid A_{1..k-1}^{\mathrm{obs}}, L)}.
#'     \item `"marginal"`: numerator drops \eqn{L}, conditioning only
#'       on prior treatments (Robins, Hernán, Brumback 2000). Fits a
#'       second per-component density model \eqn{g_k(A_k \mid A_{1..k-1})}
#'       and uses it in the numerator; the denominator stays at the
#'       full-L density. This dampens the multiplicative L-dependence
#'       across \eqn{K} factors and typically reduces weight variance
#'       substantially on datasets where multiple components have
#'       well-predicted propensities.
#'   }
#'   Ignored for univariate IPW, `"gcomp"`, and `"matching"`.
#' @param target Character or `NULL`. Column name of a binary 0/1 sampling
#'   indicator S, where S = 1 identifies study-population rows (observed A
#'   and Y) and S = 0 identifies target-population rows (only L observed).
#'   When non-`NULL`, enables transportability/generalizability estimation:
#'   the sampling model \eqn{P(S = 1 \mid L)} is fit on all rows and stored
#'   on the fit object for use in [contrast()]. Default `NULL` produces a
#'   study-population estimand with no transportability adjustment.
#'   **Not supported with `estimator = "matching"`**.
#' @param sampling_model_fn Function or `NULL`. Fitting function for the
#'   sampling model \eqn{P(S = 1 \mid L)} when `target` is non-`NULL`.
#'   Must accept `(formula, data, family, ...)`. Default `NULL` uses
#'   `stats::glm` with `family = binomial()`. Ignored when `target = NULL`.
#' @param target_subset Character. Which rows define the target population
#'   for transportability/generalizability. `"target"` (default) restricts
#'   the target to S = 0 rows only (transportability: the study is external
#'   to the target population). `"all"` uses all rows S = 0 and S = 1
#'   (generalizability: the study is a biased subsample of the target).
#'   Ignored when `target = NULL`.
#' @param treatment_free One-sided formula or `NULL`. SNM only. Specifies
#'   the treatment-free outcome model \eqn{E[Y \mid L]}, a nuisance model
#'   whose predictions are subtracted from Y before g-estimation. This
#'   projects out the \eqn{L \to Y} association and reduces the variance
#'   of \eqn{\hat\psi} without changing the point estimate (both approaches
#'   are consistent). DTRreg calls this the `tf.mod` argument. The formula
#'   should contain only baseline confounders (no treatment or modifiers).
#'   Default `NULL` uses the standard residualized-treatment moment condition
#'   without a treatment-free model. Ignored for non-SNM estimators.
#' @param stratified Character or `NULL`. Name of a baseline
#'   (time-invariant) column to stratify the ICE outcome models on.
#'   **Longitudinal g-computation (ICE) only.** When set, a separate
#'   outcome / pseudo-outcome model is fitted for each stratum at every
#'   backward step instead of a single pooled model -- the right choice
#'   when the outcome--treatment relationship differs structurally across
#'   baseline subgroups (e.g. different functional form by sex). The column
#'   must be discrete, complete, and constant within each individual.
#'   Both variance paths are available in [contrast()]: the ID-cluster
#'   bootstrap (`ci_method = "bootstrap"`) and the analytic per-stratum
#'   stacked-EE sandwich (`ci_method = "sandwich"`). Default `NULL` fits
#'   pooled models. Ignored / rejected for point treatments and non-gcomp
#'   estimators.
#' @param treatment_form One-sided formula or `NULL`. Controls how the
#'   treatment enters the per-step ICE outcome models, e.g.
#'   `~ factor(A)` (categorical / ordinal dose) or
#'   `~ splines::ns(A, df = 3)` (smooth dose-response). **Longitudinal
#'   g-computation (ICE) only**, including natural-history MTPs
#'   ([grace_period()] / [carry_forward()]). The intervention always sets
#'   the *numeric* treatment column; only the model's design term changes,
#'   so a nonlinear or kinked counterfactual response (e.g. a capped dose)
#'   is no longer forced through a single linear slope. Lag terms are
#'   expanded automatically (`factor(A)` -> `factor(lag1_A)`). Every term
#'   must reference a treatment column. Under `factor(A)` an intervention
#'   must keep the treatment within observed levels. Default `NULL` enters
#'   the treatment as a bare numeric main effect (the historical
#'   behaviour). Rejected for point treatments and non-gcomp estimators.
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
#' [contrast()] -- fit once, contrast with multiple estimands.
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
#' **Longitudinal IPW** is supported via sequential density-ratio weights
#' (cumulative product across time points) with optional stabilization.
#' Supply `id` and `time` to activate.
#'
#' ## AIPW (`estimator = "aipw"`)
#' Doubly-robust estimation combining the outcome model (as in g-comp)
#' with propensity-score weights (as in IPW). The AIPW estimator is
#' consistent if either the outcome model or the treatment model is
#' correctly specified. Supports binary, continuous, categorical,
#' count, and multivariate treatments (point); longitudinal AIPW
#' is also supported via sequential outcome models and cumulative
#' density-ratio weights.
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
#' surface E\[Y | A, L\] is unchanged by the censoring mechanism (Hernan &
#' Robins, Ch. 13). For IPW under MAR, IPCW weights are needed -- supply
#' them via `weights =`. The
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
#' Hernan MA, Robins JM (2025). *Causal Inference: What If*. Chapman &
#' Hall/CRC. Chapters 12 (IPW), 13 (g-formula), 15 (matching), 21 (ICE).
#'
#' Zivich PN, Ross RK, Shook-Sa BE, Cole SR, Edwards JK (2024). Empirical
#' sandwich variance estimator for iterated conditional expectation
#' g-computation. *Statistics in Medicine* 43:5562-5572.
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
#'     smokeintensity + smokeyrs + exercise + active + wt71,
#'   model_fn = stats::glm
#' )
#'
#' # ATT via g-computation (override estimand in contrast())
#' fit_att <- causat(
#'   nhefs,
#'   outcome = "wt82_71",
#'   treatment = "qsmk",
#'   confounders = ~ sex + age + race + education +
#'     smokeintensity + smokeyrs + exercise + active + wt71,
#'   estimand = "ATT",
#'   model_fn = stats::glm
#' )
#'
#' # IPW -- self-contained density-ratio engine; estimand fixed at
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
#'   estimand = "ATE",
#'   model_fn = stats::glm,
#'   propensity_model_fn = stats::glm
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
#'   history = Inf,
#'   model_fn = stats::glm
#' )
#'
#' # Multivariate treatment
#' fit_multi <- causat(
#'   data,
#'   outcome = "Y",
#'   treatment = c("A1", "A2"),
#'   confounders = ~ L1 + L2,
#'   model_fn = stats::glm
#' )
#'
#' # Negative binomial outcome (count data with overdispersion)
#' fit_nb <- causat(
#'   data,
#'   outcome = "Y",
#'   treatment = "A",
#'   confounders = ~ L1 + L2,
#'   model_fn = MASS::glm.nb
#' )
#'
#' # Beta regression outcome (proportions/rates in (0, 1))
#' fit_beta <- causat(
#'   data,
#'   outcome = "Y",
#'   treatment = "A",
#'   confounders = ~ L1 + L2,
#'   model_fn = betareg::betareg,
#'   family = "beta"
#' )
#' }
#'
#' @seealso [contrast()], [diagnose()], [static()], [shift()], [dynamic()]
#' @export
causat <- function(
  data,
  outcome,
  treatment,
  confounders = NULL,
  confounders_tv = NULL,
  confounders_outcome = NULL,
  confounders_treatment = NULL,
  confounders_censoring = NULL,
  confounders_sampling = NULL,
  confounders_tv_outcome = NULL,
  confounders_tv_treatment = NULL,
  estimator = c("gcomp", "ipw", "aipw", "matching", "snm"),
  family = "gaussian",
  estimand = c("ATE", "ATT", "ATC"),
  id = NULL,
  time = NULL,
  censoring = NULL,
  ipcw = FALSE,
  censoring_model_fn = NULL,
  history = Inf,
  numerator = NULL,
  weights = NULL,
  cluster = NULL,
  type = NULL,
  model_fn = stats::glm,
  propensity_model_fn = NULL,
  propensity_family = NULL,
  stabilize = c("none", "marginal"),
  target = NULL,
  sampling_model_fn = NULL,
  target_subset = c("target", "all"),
  treatment_free = NULL,
  stratified = NULL,
  treatment_form = NULL,
  ...
) {
  stabilize <- rlang::arg_match(stabilize)
  target_subset <- rlang::arg_match(target_subset)
  # Capture the call for later display in print/summary of the result.
  call <- match.call()
  # `estimator`, not `method`: avoids shadowing `MatchIt::matchit(method = ...)`,
  # which is forwarded verbatim through `...` to MatchIt. Using the same
  # argument name would create an ambiguous named match in do.call().
  estimator <- rlang::arg_match(estimator)
  estimand <- rlang::arg_match(estimand)

  # Resolve per-component confounder formulas. Each per-component
  # formula falls back to the unified formula when NULL, so callers
  # that only supply `confounders` / `confounders_tv` get identical
  # behaviour to before this feature was added.
  conf_outcome <- confounders_outcome %||% confounders
  conf_treatment <- confounders_treatment %||% confounders
  conf_censoring <- confounders_censoring %||% confounders
  conf_sampling <- confounders_sampling %||% confounders
  conf_tv_outcome <- confounders_tv_outcome %||% confounders_tv
  conf_tv_treatment <- confounders_tv_treatment %||% confounders_tv

  # Warn when the user relies on silent defaults for model fitters.
  # `missing()` detects whether the argument was supplied at the call site,
  # regardless of whether the default value happens to match what the user
  # would have chosen. The option gate lets tests suppress globally.
  suppress <- isTRUE(getOption("causatr.suppress_default_warnings"))

  if (!suppress && missing(model_fn) && !estimator %in% c("matching", "snm")) {
    rlang::warn(
      c(
        "`model_fn` not specified; defaulting to `stats::glm`.",
        i = paste0(
          "Set `model_fn` explicitly ",
          "(e.g. `model_fn = stats::glm` or `model_fn = mgcv::gam`)."
        )
      ),
      class = "causatr_model_fn_default"
    )
  }

  if (
    !suppress &&
      missing(propensity_model_fn) &&
      estimator %in% c("ipw", "aipw", "snm")
  ) {
    rlang::warn(
      c(
        "`propensity_model_fn` not specified; defaulting to `stats::glm`.",
        i = paste0(
          "Set `propensity_model_fn` explicitly if you need a ",
          "different fitter (e.g. `mgcv::gam`)."
        )
      ),
      class = "causatr_propensity_fn_default"
    )
  }

  # treatment_free is SNM-only; reject for other estimators.
  if (!is.null(treatment_free) && estimator != "snm") {
    rlang::abort(
      c(
        "`treatment_free` is only supported for `estimator = \"snm\"`.",
        i = paste0(
          "The treatment-free outcome model is an efficiency device ",
          "specific to g-estimation. For other estimators, variance ",
          "reduction comes from model specification or doubly-robust ",
          "augmentation (AIPW)."
        )
      ),
      class = "causatr_treatment_free_not_snm"
    )
  }

  # Separate fit (causat) from inference (contrast): the fit object carries
  # the nuisance models and data, but not the estimand-specific marginal
  # means. This lets one g-comp fit produce ATE, ATT, ATC, and subgroup
  # effects in subsequent contrast() calls without refitting the outcome model.
  # IPW and matching bake the estimand into the weights/matching at fit time,
  # so their fit is single-estimand -- contrast() enforces this.

  # Auto-detect point vs longitudinal from the presence of id/time.
  # `type` lets the user force one or the other -- useful for tests
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
  # are surfaced with clear error messages -- the downstream fit_* functions
  # don't need to re-validate.
  check_causat_inputs(
    # nolint: object_usage_linter
    data,
    outcome = outcome,
    treatment = treatment,
    confounders = confounders,
    confounders_tv = confounders_tv,
    confounders_outcome = confounders_outcome,
    confounders_treatment = confounders_treatment,
    confounders_censoring = confounders_censoring,
    confounders_sampling = confounders_sampling,
    confounders_tv_outcome = confounders_tv_outcome,
    confounders_tv_treatment = confounders_tv_treatment,
    estimator = estimator,
    estimand = estimand,
    id = id,
    time = time,
    history = history
  )

  # `weights = svydesign_obj` path: unpack the sampling weights into a
  # numeric vector and adopt the design's first-stage cluster (PSU) as
  # the default contrast-time cluster unless the user passed their own
  # `cluster =`. This keeps the convenience of
  # `causat(..., weights = svydesign(...))` while letting users opt
  # out of the cluster propagation if they want a non-clustered
  # sandwich on survey-weighted data.
  if (inherits(weights, "survey.design")) {
    unpacked <- unpack_svydesign(
      weights,
      data = data,
      user_cluster = cluster,
      call = call
    )
    weights <- unpacked$weights
    cluster <- unpacked$cluster
  }

  # prepare_data() does three things:
  #   1. Coerces `data` to data.table (for fast subsetting + `:=`).
  #   2. For longitudinal data, sorts by (id, time) and materializes
  #      lag columns (`lag1_A`, `lag2_A`, ...) up to `history`.
  #   3. Validates that the person-period structure is rectangular
  #      (every id observed at every time, or consistently censored).
  # All downstream fit_* functions assume the prepared shape.
  # `cluster` (if provided) must exist in `data` and is passed to
  # `prepare_data()` so the column survives the column-stripping step
  # and is available to `contrast()` at variance time.
  if (!is.null(cluster)) {
    check_string(cluster)
    if (!cluster %in% names(data)) {
      rlang::abort(
        paste0(
          "`cluster` column '",
          cluster,
          "' not found in `data`."
        ),
        class = "causatr_bad_cluster",
        call = call
      )
    }
    if (estimator == "matching") {
      rlang::abort(
        c(
          "`cluster` is not supported for `estimator = \"matching\"`.",
          i = paste0(
            "Matching already aggregates IFs cluster-robustly on the ",
            "matched `subclass`. Use `estimator = \"gcomp\"` or ",
            "`\"ipw\"` if you need a design cluster."
          )
        ),
        class = "causatr_bad_cluster",
        call = call
      )
    }
  }

  data <- prepare_data(
    # nolint: object_usage_linter
    data,
    outcome = outcome,
    treatment = treatment,
    confounders = confounders,
    confounders_tv = confounders_tv,
    confounders_outcome = confounders_outcome,
    confounders_treatment = confounders_treatment,
    confounders_censoring = confounders_censoring,
    confounders_sampling = confounders_sampling,
    confounders_tv_outcome = confounders_tv_outcome,
    confounders_tv_treatment = confounders_tv_treatment,
    id = id,
    time = time,
    censoring = censoring,
    history = history,
    cluster = cluster,
    target = target,
    stratified = stratified
  )

  # Stratified ICE: validate the stratifying column on the prepared data
  # (data.table; the column survived `prepare_data()`'s keep set). Rejects
  # non-ICE estimators, time-varying / continuous / NA stratifiers. Runs
  # here so a misuse is caught before any model is fitted.
  check_stratified(data, stratified, estimator, type, id, call = call)

  # Flexible-treatment ICE term: validate `treatment_form` (ICE-only,
  # one-sided formula referencing the treatment column(s)) before any model
  # is fitted. The intervention still sets the numeric treatment column; only
  # the per-step design term changes.
  check_treatment_form(treatment_form, treatment, estimator, type, call = call)

  # NA check on treatment values: if any are missing, user must either
  # provide a censoring column (IPCW), use mice imputation, or remove
  # incomplete cases manually. We do this AFTER prepare_data() because
  # lag materialization is what actually exposes the NAs at baseline.
  check_treatment_nas(data, treatment, censoring, target = target)

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
  # NAs -- recycling then silently corrupts the IF and sandwich SEs.
  # See check_dots_na_action() for the full rationale.
  check_dots_na_action(..., call = call)

  # Transportability: validate S column and fit P(S=1|L) before dispatch.
  # The sampling model is stored on the fit object but not composed into
  # `weights` -- the sampling-weight role differs by estimator (gcomp:
  # standardize over target rows; IPW: multiply with treatment weights).
  sampling_details <- NULL
  if (!is.null(target)) {
    check_transport_inputs(
      target = target,
      target_col = data[[target]],
      target_subset = target_subset,
      estimator = estimator,
      estimand = estimand,
      call = call
    )
    samp_fn <- sampling_model_fn %||% stats::glm
    # For longitudinal data, S is constant across time for each individual.
    # Fit the sampling model on first-period rows only to avoid K-fold
    # inflation of the likelihood and the sandwich score.
    if (type == "longitudinal") {
      first_t <- min(data[[time]])
      rows_first_samp <- data[[time]] == first_t
      samp_data <- data[rows_first_samp]
      samp_weights <- if (is.null(weights)) NULL else weights[rows_first_samp]
    } else {
      samp_data <- data
      samp_weights <- weights
    }
    samp_model <- fit_sampling_model(
      data = samp_data,
      target = target,
      confounders = conf_sampling,
      treatment = treatment,
      model_fn = samp_fn,
      weights = samp_weights
    )
    sampling_details <- list(
      transport = TRUE,
      sampling_model = samp_model,
      sampling_model_fn = samp_fn,
      target_subset = target_subset
    )
  }

  # Built-in IPCW: fit a censoring model and compose stabilized
  # weights with any external weights BEFORE dispatching to the
  # estimator-specific fitter. Point treatments use a single
  # censoring model; longitudinal treatments use per-period
  # censoring models via fit_censoring_models_longitudinal().
  ipcw_details <- NULL
  if (ipcw && type == "point") {
    check_ipcw_inputs(
      ipcw,
      censoring,
      censoring_col = data[[censoring]],
      call = call
    )
    cens_model_fn <- censoring_model_fn %||% stats::glm
    cens_model <- fit_censoring_model(
      data = data,
      censoring = censoring,
      treatment = treatment,
      confounders = conf_censoring,
      model_fn = cens_model_fn,
      weights = weights
    )
    ipcw_w <- compute_ipcw_weights(
      cens_model,
      n_total = nrow(data),
      censoring_col = as.integer(data[[censoring]]),
      stabilize = TRUE
    )
    ipcw_details <- list(
      ipcw = TRUE,
      censoring_model = cens_model,
      ipcw_weights = ipcw_w,
      censoring_model_fn = cens_model_fn,
      weights_pre_ipcw = weights
    )
    weights <- if (is.null(weights)) ipcw_w else weights * ipcw_w
    check_weights(weights, nrow(data))
  } else if (ipcw && type == "longitudinal") {
    check_ipcw_inputs(
      ipcw,
      censoring,
      censoring_col = data[[censoring]],
      call = call
    )
    cens_model_fn <- censoring_model_fn %||% stats::glm
    time_points_ipcw <- sort(unique(data[[time]]))
    cens_result <- fit_censoring_models_longitudinal(
      data = data,
      censoring = censoring,
      treatment = treatment,
      confounders = conf_censoring,
      confounders_tv = conf_tv_outcome,
      model_fn = cens_model_fn,
      id = id,
      time = time,
      time_points = time_points_ipcw,
      history = history,
      weights = weights
    )
    ipcw_w <- cens_result$cumulative_weights
    ipcw_details <- list(
      ipcw = TRUE,
      censoring_models = cens_result$models,
      ipcw_weights = ipcw_w,
      censoring_model_fn = cens_model_fn,
      weights_pre_ipcw = weights,
      per_period_weights = cens_result$per_period_weights
    )
    weights <- if (is.null(weights)) ipcw_w else weights * ipcw_w
    check_weights(weights, nrow(data))
  }

  # Dispatch to the estimator-specific fitter. Each returns a
  # `causatr_fit` with the same S3 class and slot structure, which
  # contrast() and diagnose() then consume uniformly.
  #
  # ICE (gcomp longitudinal) defers fitting the K period outcome models to
  # contrast() rather than doing it here: the pseudo-outcome sequence at each
  # time step depends on the intervention being applied (which regime to set at
  # period k determines the backward-recursion starting values). Because two
  # different interventions yield two different sets of period models, there is
  # no single "the outcome model" to store on the fit; contrast() runs
  # ice_iterate() per intervention instead. The fit$model slot holds NULL for
  # ICE; fit$data and fit$details carry everything ice_iterate() needs.
  fit <- switch(
    estimator,
    gcomp = fit_gcomp(
      data,
      outcome,
      treatment,
      conf_outcome,
      conf_tv_outcome,
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
      target = target,
      confounders_outcome = confounders_outcome,
      confounders_treatment = confounders_treatment,
      confounders_censoring = confounders_censoring,
      confounders_sampling = confounders_sampling,
      confounders_tv_outcome = confounders_tv_outcome,
      confounders_tv_treatment = confounders_tv_treatment,
      stratified = stratified,
      treatment_form = treatment_form,
      ...
    ),
    ipw = fit_ipw(
      data,
      outcome,
      treatment,
      conf_treatment,
      conf_tv_treatment,
      family,
      estimand,
      type,
      history,
      numerator,
      weights,
      model_fn,
      propensity_model_fn,
      propensity_family,
      stabilize = stabilize,
      id = id,
      time = time,
      call = call,
      target = target,
      confounders_outcome = conf_outcome,
      confounders_outcome_raw = confounders_outcome,
      confounders_treatment_raw = confounders_treatment,
      confounders_censoring_raw = confounders_censoring,
      confounders_sampling_raw = confounders_sampling,
      confounders_tv_outcome_raw = confounders_tv_outcome,
      confounders_tv_treatment_raw = confounders_tv_treatment,
      ...
    ),
    aipw = fit_aipw(
      data,
      outcome,
      treatment,
      conf_outcome,
      conf_tv_outcome,
      family,
      estimand,
      type,
      history,
      censoring,
      weights,
      model_fn,
      propensity_model_fn,
      propensity_family,
      stabilize = stabilize,
      id = id,
      time = time,
      call = call,
      target = target,
      confounders_treatment = conf_treatment,
      confounders_tv_treatment = conf_tv_treatment,
      confounders_outcome_raw = confounders_outcome,
      confounders_treatment_raw = confounders_treatment,
      confounders_censoring_raw = confounders_censoring,
      confounders_sampling_raw = confounders_sampling,
      confounders_tv_outcome_raw = confounders_tv_outcome,
      confounders_tv_treatment_raw = confounders_tv_treatment,
      ...
    ),
    matching = fit_matching(
      data,
      outcome,
      treatment,
      conf_outcome,
      family,
      estimand,
      type,
      weights,
      call,
      confounders_treatment = conf_treatment,
      confounders_outcome_raw = confounders_outcome,
      confounders_treatment_raw = confounders_treatment,
      confounders_censoring_raw = confounders_censoring,
      confounders_sampling_raw = confounders_sampling,
      confounders_tv_outcome_raw = confounders_tv_outcome,
      confounders_tv_treatment_raw = confounders_tv_treatment,
      ...
    ),
    snm = fit_snm(
      data,
      outcome,
      treatment,
      conf_treatment,
      conf_tv_treatment,
      family,
      estimand,
      type,
      history,
      weights,
      propensity_model_fn,
      propensity_family,
      id = id,
      time = time,
      call = call,
      target = target,
      confounders_outcome = conf_outcome,
      confounders_outcome_raw = confounders_outcome,
      confounders_treatment_raw = confounders_treatment,
      confounders_censoring_raw = confounders_censoring,
      confounders_sampling_raw = confounders_sampling,
      confounders_tv_outcome_raw = confounders_tv_outcome,
      confounders_tv_treatment_raw = confounders_tv_treatment,
      treatment_free = treatment_free,
      model_fn = model_fn,
      ...
    ),
    # Defensive default: unreachable under normal use because
    # rlang::arg_match(estimator) above restricts `estimator` to the
    # allowed set. Kept so a future refactor that loosens arg_match
    # cannot silently return NULL from causat().
    rlang::abort(c(
      paste0("Unknown `estimator` '", estimator, "'."),
      i = "Must be one of: 'gcomp', 'ipw', 'aipw', 'matching', 'snm'."
    ))
  )

  # Stash the cluster column name on the fit so `contrast()` can
  # default its own `cluster =` argument without the user having to
  # repeat it. The cluster is not used during fitting -- it only enters
  # the sandwich variance aggregation in contrast(). Threading it through
  # every `fit_*()` signature would touch four files for one assignment;
  # we do it once here instead.
  if (!is.null(cluster)) {
    fit$details$cluster <- cluster
  }

  if (!is.null(ipcw_details)) {
    for (nm in names(ipcw_details)) {
      fit$details[[nm]] <- ipcw_details[[nm]]
    }
    # Ensure the censoring column name is on the fit object even for
    # estimators that don't receive it directly (e.g. fit_ipw sets
    # censoring = NULL). The variance engine needs it to look up the
    # censoring column in fit$data.
    if (is.null(fit$censoring)) {
      fit$censoring <- censoring
    }
  }

  if (!is.null(sampling_details)) {
    for (nm in names(sampling_details)) {
      fit$details[[nm]] <- sampling_details[[nm]]
    }
    fit$target <- target
    fit$target_subset <- target_subset
  }

  fit
}


#' Unpack a `survey::svydesign` object into weights + cluster
#'
#' @description
#' Translates a `survey.design` object into the `weights` numeric vector
#' and `cluster` column name that `causat()` internally consumes. The
#' design's inverse sampling probabilities are extracted via
#' `stats::weights()` (which dispatches to `survey:::weights.survey.design`),
#' and the first-stage cluster (PSU) is adopted
#' as the default contrast-time cluster unless the user has supplied
#' their own `cluster =` at the `causat()` boundary (in which case the
#' user's value is kept and the PSU is ignored).
#'
#' Validates that the design's `nrow(design$variables) == nrow(data)` so
#' the weight vector aligns row-for-row with `data`. Requires the
#' `survey` package to be installed (Suggests-only); aborts with a
#' classed error otherwise.
#'
#' @param design A `survey::svydesign` / `survey.design` object.
#' @param data The `data` frame passed to `causat()`.
#' @param user_cluster The user's explicit `cluster =` argument (or
#'   `NULL`); the design's PSU is only adopted when this is `NULL`.
#' @param call Caller environment for error messages.
#'
#' @return A list with two components: `weights` (numeric vector
#'   aligned to `data`) and `cluster` (either the adopted PSU column
#'   name, the user's explicit override, or `NULL` if neither is
#'   available).
#'
#' @noRd
unpack_svydesign <- function(
  design,
  data,
  user_cluster = NULL,
  call = rlang::caller_env()
) {
  check_pkg("survey")

  n_design <- nrow(design$variables)
  n_data <- nrow(data)
  if (n_design != n_data) {
    rlang::abort(
      c(
        paste0(
          "`svydesign` row count (",
          n_design,
          ") does not match `data` row count (",
          n_data,
          ")."
        ),
        i = paste0(
          "Pass the same `data` to `causat()` and `svydesign()`, with ",
          "rows in identical order. The weight vector must align ",
          "row-for-row."
        )
      ),
      class = "causatr_bad_svydesign",
      call = call
    )
  }

  # `weights()` is a generic; `survey` ships
  # `weights.survey.design()` which returns the design's sampling
  # weights. We go through `stats::weights()` (the generic) rather
  # than referencing the method directly so this works whether the
  # dispatch table finds the method via `survey:::weights.survey.design`
  # or via the base `weights.default`.
  w <- stats::weights(design)
  if (length(w) != n_data) {
    rlang::abort(
      paste0(
        "`survey::weights(design)` returned ",
        length(w),
        " values but `data` has ",
        n_data,
        " rows."
      ),
      class = "causatr_bad_svydesign",
      call = call
    )
  }

  # First-stage PSU. `design$cluster` is a data.frame with one column
  # per cluster stage (column 1 = first stage / primary sampling unit).
  # We adopt column 1 only when (a) the user hasn't given an explicit
  # override, (b) the column name survives into `data` unchanged
  # (which is always true when the design was built on the same
  # `data`), and (c) the cluster is nontrivial (> 1 unique level -- a
  # single-PSU design has no within-cluster correlation to model).
  adopted_cluster <- user_cluster
  if (is.null(user_cluster) && !is.null(design$cluster)) {
    psu_name <- names(design$cluster)[1]
    if (
      !is.null(psu_name) &&
        psu_name %in% names(data) &&
        length(unique(design$cluster[[1]])) > 1L
    ) {
      adopted_cluster <- psu_name
    }
  }

  list(weights = as.numeric(w), cluster = adopted_cluster)
}
