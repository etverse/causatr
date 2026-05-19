#' Compute causal contrasts from a fitted model
#'
#' @description
#' Standardises outcome model predictions under each named intervention and
#' reports pairwise causal contrasts (differences, ratios, or odds ratios)
#' with uncertainty estimates.
#'
#' For all three causal estimators, `contrast()` computes E\[Y^a\] by setting
#' each individual's treatment to the intervened value and averaging the fitted
#' outcome model predictions over the target population. The target population
#' is controlled by `estimand` (or `subset` for subgroup effects). The
#' estimators differ only in how the outcome model was fitted; variance is
#' computed by the
#' **unified influence-function engine** `variance_if()` (`R/variance_if.R`),
#' which handles all four cases (g-comp, IPW, matching, ICE) via one entry
#' point. See `vignettes/variance-theory.qmd` for the derivation and the
#' `R/variance_if.R` roxygen header for the architecture.
#'
#' - `"gcomp"`: standard `glm`/`gam` on the full data; one-model IF
#'   correction via `prepare_model_if()` + `apply_model_correction()`.
#' - `"ipw"`: density-ratio weighted MSM refit per intervention; IF
#'   Channel 2 is the stacked-sandwich MSM-plus-propensity correction
#'   assembled by `compute_ipw_if_self_contained_one()`, with the
#'   cross-derivative \eqn{A_{\beta\alpha}} computed numerically by
#'   `numDeriv::jacobian()`.
#' - `"matching"`: `glm()` on the matched sample with match weights; IF is
#'   computed on the matched sample and aggregated cluster-robustly on
#'   `matched$subclass` via `vcov_from_if(cluster = ...)`.
#'
#' @param fit A `causatr_fit` object returned by [causat()].
#' @param interventions A named list of interventions. Each element must be
#'   one of:
#'   - A `causatr_intervention` object created by [static()], [shift()],
#'     [dynamic()], [scale_by()], [threshold()], [ipsi()], or
#'     [stochastic()].
#'   - `NULL`, meaning the natural course (observed treatment values are used
#'     as-is). The natural course is the standard reference for modified
#'     treatment policies on continuous treatments (e.g. `shift(-10)` vs
#'     `NULL`; Diaz et al. 2023) and for longitudinal dynamic regimes
#'     (Hernan & Robins Ch. 21). For binary treatments, the
#'     natural-course marginal mean equals E\[Y\] under the observed
#'     treatment mechanism (by consistency); use `static(1)` vs
#'     `static(0)` for the conventional ATE. Supported for all estimators.
#'   - A named list of `causatr_intervention` objects, one per treatment
#'     variable, for multivariate (joint) treatments. Each sub-list must name
#'     every treatment variable specified in `causat()`.
#'
#'   **Note:** Non-static interventions (`shift()`, `scale_by()`, `threshold()`,
#'   `dynamic()`, `ipsi()`, `stochastic()`) are only supported for
#'   `estimator = "gcomp"`.
#'   For IPW and matching, the weights/matched sets were estimated under the
#'   observed treatment distribution; applying a different regime requires
#'   re-calling [causat()] with updated data. `ipsi()` additionally requires
#'   a fitted propensity model and currently aborts for all estimators.
#' @param type Character. The contrast scale: `"difference"` (default),
#'   `"ratio"`, or `"or"` (odds ratio). All pairwise contrasts are reported.
#' @param estimand Character or `NULL`. The target estimand: `"ATE"`,
#'   `"ATT"`, or `"ATC"`. For `estimator = "gcomp"`, overrides the estimand
#'   stored in `fit` (allowing one fit to produce multiple estimands). For
#'   `estimator = "ipw"` or `"matching"`, must match the estimand used at
#'   fitting time -- changing it aborts with an informative message. If
#'   `NULL`, defaults to `fit$estimand`. Mutually exclusive with `subset`.
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
#' @param ci_method Character. The variance/CI approach: `"sandwich"` (default)
#'   or `"bootstrap"`.
#' @param n_boot Integer. Number of bootstrap replications when
#'   `ci_method = "bootstrap"`. Default `500`.
#' @param conf_level Numeric. Confidence level for intervals. Default `0.95`.
#' @param by Character or `NULL`. Name of a variable to stratify estimates by
#'   (effect modification). If provided, E\[Y^a\] is computed within each
#'   level of `by`.
#' @param parallel Character. Parallelisation backend for bootstrap
#'   (`ci_method = "bootstrap"` only). `"no"` (default) runs sequentially;
#'   `"multicore"` uses forked processes via [parallel::mclapply()] (Unix
#'   only); `"snow"` uses socket clusters via [parallel::parLapply()]
#'   (cross-platform); `"future"` dispatches replicates through
#'   [future.apply::future_lapply()] so any active [future::plan()] is
#'   honoured (local multisession, remote cluster, `future.batchtools`,
#'   etc.). The `"no" / "multicore" / "snow"` paths go through
#'   [boot::boot()]; `"future"` bypasses it to let `future` own the
#'   scheduling, which means the `ncpus` argument is ignored for
#'   `"future"` (the plan's worker count governs). Ignored when
#'   `ci_method = "sandwich"`.
#' @param ncpus Integer. Number of CPU cores for parallel bootstrap. Default
#'   `getOption("boot.ncpus", 1L)`. Passed directly to [boot::boot()].
#' @param cluster Character or `NULL`. Name of a column in `fit$data`
#'   identifying cluster membership (e.g. site, household, PSU id). When
#'   non-`NULL`, the sandwich variance is aggregated cluster-robustly:
#'   per-individual influence functions are **summed within each cluster
#'   before squaring** (Liang & Zeger 1986), and the resulting
#'   \eqn{(1/n^2)\sum_c (\sum_{i \in c} \mathrm{IF}_i)^2} estimator is
#'   numerically equivalent to the stacked-EE cluster-robust sandwich and
#'   to [sandwich::vcovCL()] applied on the final predict-then-average
#'   step. Works for `estimator = "gcomp"`, `"ipw"`, and `"ice"`; aborts
#'   under `"matching"` (matching already clusters on its own
#'   `subclass`, and layering a user cluster on top would double-count).
#'   For ICE with a longitudinal fit, the cluster vector is read from the
#'   first-time-point rows. Cluster is unused when
#'   `ci_method = "bootstrap"` for point fits; ICE bootstrap already
#'   resamples entire individual trajectories.
#'
#' @return A `causatr_result` object with slots:
#'   \describe{
#'     \item{`estimates`}{data.table with one row per intervention:
#'       `intervention`, `estimate`, `se`, `ci_lower`, `ci_upper`.}
#'     \item{`contrasts`}{data.table with one row per pairwise comparison:
#'       `comparison`, `estimate`, `se`, `ci_lower`, `ci_upper`.}
#'     \item{`type`}{Contrast scale.}
#'     \item{`estimand`}{`"ATE"`, `"ATT"`, `"ATC"`, or `"subset"`.}
#'     \item{`ci_method`}{Inference approach.}
#'     \item{`reference`}{Name of the reference intervention.}
#'     \item{`interventions`}{The intervention list.}
#'     \item{`n`}{Number of individuals averaged over.}
#'     \item{`estimator`}{Causal estimator from the fit.}
#'     \item{`vcov`}{Full variance-covariance matrix for all E\[Y^a\].}
#'     \item{`call`}{The original call.}
#'   }
#'
#' @details
#' ## Estimands and standardisation
#' For each intervention `a`, `contrast()` evaluates the intervention function
#' on the target rows to obtain the intervened treatment vector
#' \eqn{a(L_i)}, then computes:
#'
#' \deqn{E[Y^a] = (1/|S|) \sum_{i \in S} \hat{E}[Y | A = a(L_i), L_i]}
#' where `S` is the set of rows determined by the estimand:
#' - `"ATE"`: all rows (full population)
#' - `"ATT"`: rows where `A == 1` (observed treated)
#' - `"ATC"`: rows where `A == 0` (observed controls)
#' - `subset`: rows satisfying the quoted expression
#'
#' For `"gcomp"`, one fit supports multiple estimands because the outcome
#' model is the same -- only the rows averaged over change. For `"ipw"` and
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
#' | `"ipw"` | binary, continuous | all deterministic |
#' | `"matching"` | binary | `static()` only |
#'
#' ## Variance estimation
#' - `"sandwich"`: The unified influence-function engine `variance_if()`.
#'   Per-individual IF = Channel 1 (direct covariate-sampling term) +
#'   Channel 2 (one correction per nuisance model). Aggregated via
#'   `vcov_from_if()`, with cluster-robust aggregation for matching.
#'   Mathematically equivalent to the stacked M-estimation sandwich; see
#'   `vignettes/variance-theory.qmd` Sections 3-4.
#' - `"bootstrap"`: Resamples individuals (entire pipeline refitted `n_boot`
#'   times). Respects cluster structure for longitudinal data.
#' - The delta method is applied internally for ratio / odds-ratio
#'   contrasts on top of the `"sandwich"` vcov; it is not a separate
#'   `ci_method` option.
#'
#' @references
#' Hernan MA, Robins JM (2025). *Causal Inference: What If*. Chapman &
#' Hall/CRC. Chapters 12-13.
#'
#' Imai K, King G, Stuart EA (2008). Misunderstandings between experimentalists
#' and observationalists about causal inference. *Journal of the Royal
#' Statistical Society* Series A 171:481-502.
#'
#' Zivich PN, Ross RK, Shook-Sa BE, Cole SR, Edwards JK (2024). Empirical
#' sandwich variance estimator for iterated conditional expectation
#' g-computation. *Statistics in Medicine* 43:5562-5572.
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
#'
#' # Stochastic intervention: random treatment assignment
#' sampler <- function(data, trt) {
#'   rbinom(nrow(data), 1, plogis(0.3 * data$age))
#' }
#' contrast(fit,
#'   interventions = list(
#'     stoch = stochastic(sampler, n_mc = 100L),
#'     always = static(1)
#'   ),
#'   type = "difference"
#' )
#' }
#'
#' @seealso [causat()], [static()], [shift()], [dynamic()],
#'   [stochastic()], [coef.causatr_result()],
#'   [confint.causatr_result()]
#' @export
contrast <- function(
  fit,
  interventions,
  type = c("difference", "ratio", "or"),
  estimand = NULL,
  subset = NULL,
  reference = NULL,
  ci_method = c("sandwich", "bootstrap"),
  n_boot = 500L,
  conf_level = 0.95,
  by = NULL,
  parallel = getOption("boot.parallel", "no"),
  ncpus = getOption("boot.ncpus", 1L),
  cluster = NULL
) {
  # `parallel` accepts the same three values as `boot::boot()` and
  # defaults to `getOption("boot.parallel", "no")`, so a session-wide
  # `options(boot.parallel = "multicore")` flips bootstrap parallelism
  # on for every `contrast()` call without per-call plumbing.
  parallel <- match.arg(parallel, c("no", "multicore", "snow", "future"))
  # Capture the caller's frame ONCE, at the top level, so a `subset`
  # expression referring to a session variable (e.g. `quote(age >
  # cutoff)`) can still resolve `cutoff` deep inside
  # the compute_contrast, variance_bootstrap, ice_variance_bootstrap paths
  # where `parent.frame()` no longer points at the user. Threading
  # `subset_env` explicitly instead of re-capturing downstream is the
  # only correct way -- see B1 in the 2026-04-15 critical review.
  subset_env <- parent.frame()
  # Capture the original call so that the returned `causatr_result`
  # can echo it in its `print` and `summary` methods. `match.call()`
  # here -- not at compute_contrast() -- so the recorded call reflects
  # what the *user* typed, not the internal dispatch.
  call <- match.call()

  if (!inherits(fit, "causatr_fit")) {
    rlang::abort("`fit` must be a `causatr_fit` object returned by `causat()`.")
  }

  # Canonicalize character choices. `rlang::arg_match()` both validates
  # and returns the single chosen element, so `type` and `ci_method`
  # are guaranteed scalar after these lines.
  type <- rlang::arg_match(type)
  ci_method <- rlang::arg_match(ci_method)

  # `subset` must be a *quoted* expression, not a logical vector. We
  # accept language objects so the expression can be evaluated against
  # the data inside `compute_contrast()` -- this is how Hernan/Robins-
  # style subgroup effects like `quote(age > 50)` work.
  if (!is.null(subset) && !is.language(subset)) {
    rlang::abort(
      paste0(
        "`subset` must be a quoted expression (e.g. `quote(age > 50)`), ",
        "not a ",
        class(subset)[1],
        "."
      )
    )
  }

  # `estimand` and `subset` are mutually exclusive: both would fight
  # over defining the target population. An explicit `estimand` picks
  # a pre-specified group (e.g. ATT = "treated at baseline"); a
  # `subset` picks a user-defined expression. Pick one.
  if (!is.null(estimand) && !is.null(subset)) {
    rlang::abort("Specify either 'estimand' or 'subset', not both.")
  }

  # Estimand compatibility checks: only g-comp can override the fitted
  # estimand (one outcome model, re-average over a different target).
  # IPW and matching bake the estimand into the weights/matching, so
  # overriding here would silently break the inverse-weighting identity.
  if (!is.null(estimand)) {
    estimand <- rlang::arg_match(estimand, c("ATE", "ATT", "ATC"))
    check_estimand_compat(estimand, fit$estimator, fit$estimand)
    check_estimand_trt_compat(
      estimand,
      fit$treatment,
      fit$type,
      data = fit$data
    )
  }

  # Validate the interventions list -- names, types, and estimator
  # compatibility (IPW/matching only accept static(), ICE accepts all).
  check_intervention_list(interventions)
  check_interventions_compat(fit$estimator, interventions)

  # Estimand x intervention gate.
  # ATT / ATC are only well-defined under static interventions for the IPW
  # and matching engines; silently falling back to ATE weights under an
  # ATT request would return a pooled effect when the user asked for a
  # subpopulation effect. Abort upfront with the
  # `causatr_bad_estimand_intervention` class. We pass the **effective**
  # estimand -- the user's `estimand = ` override if present, otherwise
  # `fit$estimand` -- so a fit made under ATT that is contrast-time
  # overridden to ATE (gcomp only) still goes through the right branch.
  effective_estimand <- if (is.null(estimand)) fit$estimand else estimand
  check_estimand_intervention_compat(
    effective_estimand,
    interventions,
    fit$estimator
  )

  # `reference` names the intervention used as the contrast denominator
  # (pairwise a_j vs a_ref). Must exist in the named list.
  if (!is.null(reference) && !reference %in% names(interventions)) {
    rlang::abort(
      paste0(
        "`reference` ('",
        reference,
        "') must be the name of one of the interventions."
      )
    )
  }

  # `by` stratifies results by levels of a data column -- e.g.
  # `by = "sex"` returns separate estimates per sex category. The
  # compute_contrast() loop handles the actual stratification.
  if (!is.null(by)) {
    check_string(by)
    if (!by %in% names(fit$data)) {
      rlang::abort(
        paste0("`by` variable '", by, "' not found in fitted data.")
      )
    }
  }

  # `parallel` is already normalised at the top of `contrast()` via
  # an explicit `match.arg()` against the full choice list, so we
  # cannot use `rlang::arg_match()` here (it derives choices from
  # the formal's default, which is a single resolved option string).

  # Cluster resolution: validate and extract the cluster vector (or
  # return NULL) before dispatching, so downstream variance branches
  # can always trust a validated vector or a clean NULL. `cluster`
  # passed at the fit stage is read back from `fit$details$cluster`
  # and treated as the default when the user doesn't override it --
  # this is what lets `causat(weights = svydesign(...))` auto-
  # propagate its PSU into the sandwich without a second argument.
  if (is.null(cluster)) {
    cluster <- fit$details$cluster
  }
  cluster_vec <- resolve_cluster(cluster, fit)

  # Hand off to the internal engine. Everything above is argument
  # validation; the actual math lives in `compute_contrast()` and
  # its estimator-specific delegates.
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
    parallel,
    ncpus,
    call,
    subset_env,
    cluster_vec = cluster_vec
  )
}

#' Check that interventions are compatible with the causal estimator
#'
#' IPW and matching only support static interventions; non-static
#' interventions require g-computation.
#'
#' @param estimator Character causal estimator.
#' @param interventions Named list of `causatr_intervention` objects.
#' @param call Caller environment for error messages.
#' @return `NULL` invisibly; aborts if incompatible interventions are found.
#' @noRd
check_interventions_compat <- function(
  estimator,
  interventions,
  call = rlang::caller_env()
) {
  # Matching is the only estimator gated here. IPW handles non-static
  # interventions via the self-contained density-ratio engine
  # (`fit_treatment_model()` + `make_weight_fn()` + per-intervention
  # weighted MSM refit inside `compute_contrast()`); the per-intervention
  # compatibility check lives in `check_intervention_family_compat()`
  # in `R/ipw_weights.R` and fires inside `compute_density_ratio_weights()`.
  # G-comp supports every intervention type because it re-predicts from
  # an outcome model that doesn't care about the treatment-assignment
  # mechanism.
  if (estimator == "matching") {
    # For each element in the interventions list, determine whether it
    # represents a non-static regime. Three shapes to handle:
    #   - NULL entry         -> natural course (always allowed)
    #   - bare intervention  -> check `$type` directly
    #   - list of interventions (multivariate treatment) -> check every sub
    non_static <- vapply(
      interventions,
      function(iv) {
        if (is.null(iv)) {
          return(FALSE)
        }
        # Multivariate treatment: `iv` is a plain list of sub-interventions
        # (one per treatment component), not a `causatr_intervention`.
        # Any non-static sub flags the whole regime as non-static.
        if (is.list(iv) && !inherits(iv, "causatr_intervention")) {
          return(any(vapply(
            iv,
            function(sub) sub$type != "static",
            logical(1)
          )))
        }
        iv$type != "static"
      },
      logical(1)
    )
    if (any(non_static)) {
      # Hard abort with a pointer to the fix. The semantic problem is
      # that IPW weights are computed once at fit time under the
      # observed treatment distribution; matching picks controls once
      # against the observed treated group. Neither is valid under a
      # shift / MTP / dynamic rule, because the density ratio or
      # balancing property no longer holds. G-comp's outcome model is
      # intervention-agnostic, so it handles all regimes uniformly.
      rlang::abort(
        paste0(
          "Non-static interventions ",
          "(shift, dynamic, scale, threshold, ipsi, stochastic) ",
          "are not supported for estimator = '",
          estimator,
          "'. ",
          "The weights/matched sets were estimated under the original ",
          "treatment regime and are not valid under a different intervention. ",
          "Use estimator = 'gcomp' instead."
        ),
        call = call
      )
    }
  }
}

#' Predict under an intervention (deterministic or stochastic MC)
#'
#' @description
#' For deterministic interventions, returns a single prediction vector.
#' For stochastic interventions, runs `n_mc` Monte Carlo draws through the
#' sampler, predicts under each, and returns the row-wise MC average.
#'
#' @param model Fitted model object (GLM, GAM, or custom).
#' @param data data.table. Full dataset.
#' @param treatment Character. Treatment column name(s).
#' @param iv A `causatr_intervention` or `NULL`.
#'
#' @return A list with:
#'   - `preds`: length-n numeric prediction vector (MC-averaged for stochastic).
#'   - `data_a`: one counterfactual data.table (last MC draw for stochastic).
#'
#' @noRd
predict_under_intervention <- function(model, data, treatment, iv) {
  if (!has_stochastic_component(iv)) {
    data_a <- apply_intervention(data, treatment, iv)
    preds <- stats::predict(model, newdata = data_a, type = "response")
    return(list(preds = preds, data_a = data_a))
  }

  n_mc <- get_stochastic_n_mc(iv)
  one_draw <- function(m) {
    data_a_m <- apply_intervention(data, treatment, iv)
    stats::predict(model, newdata = data_a_m, type = "response")
  }
  preds_mat <- mc_sapply(seq_len(n_mc), one_draw, n_rows = nrow(data))
  data_a_last <- apply_intervention(data, treatment, iv)
  preds_avg <- rowMeans(preds_mat, na.rm = TRUE)
  list(preds = preds_avg, data_a = data_a_last)
}


#' Fit a treatment model for MC marginalization under transport
#'
#' When an MTP intervention (shift, scale_by, threshold, dynamic) is
#' combined with transportability, target-population rows have no
#' observed treatment. This helper fits a simple GLM for
#' \eqn{f(A \mid L, S = 1)} on study rows so that treatment values
#' can be drawn (or enumerated for binary) to marginalize the
#' outcome-model prediction over the study treatment distribution.
#'
#' @param data data.table. Full dataset (study + target).
#' @param treatment Character. Treatment column name.
#' @param confounders One-sided formula.
#' @param fit_rows Logical vector. Study rows (S = 1) with observed outcome.
#' @return A fitted GLM object.
#' @noRd
fit_mc_treatment_model <- function(data, treatment, confounders, fit_rows) {
  study_data <- data[fit_rows]
  trt_vals <- study_data[[treatment]]
  is_binary <- length(unique(stats::na.omit(trt_vals))) <= 2L &&
    all(trt_vals %in% c(0L, 1L, 0, 1), na.rm = TRUE)
  fam <- if (is_binary) stats::binomial() else stats::gaussian()
  rhs <- deparse(confounders[[2]], width.cutoff = 500L)
  fml <- stats::as.formula(paste(treatment, "~", rhs))
  stats::glm(fml, data = study_data, family = fam)
}


#' MC-marginalize outcome predictions for target rows
#'
#' For target rows where treatment \eqn{A} is unobserved, computes
#' \deqn{E_{A \mid L, S=1}[\hat{m}(d(A, L), L) \mid L_i]}
#' by either exact enumeration (binary treatment) or Monte Carlo
#' integration (continuous treatment). Used by gcomp and AIPW Term 1
#' when MTP interventions are combined with transportability.
#'
#' @param outcome_model Fitted outcome model.
#' @param mc_data data.table. Target rows needing marginalization.
#' @param treatment Character. Treatment column name.
#' @param iv A `causatr_intervention` or NULL.
#' @param treatment_model Fitted GLM from `fit_mc_treatment_model()`.
#' @param n_mc Integer. MC draws for continuous treatment (ignored for binary).
#' @return Numeric vector of marginalized predictions (length = nrow(mc_data)).
#' @noRd
mc_marginalize_preds <- function(
  outcome_model,
  mc_data,
  treatment,
  iv,
  treatment_model,
  n_mc = 50L
) {
  is_binary <- treatment_model$family$family == "binomial"
  if (is_binary) {
    mc_marginalize_binary(
      outcome_model,
      mc_data,
      treatment,
      iv,
      treatment_model
    )
  } else {
    mc_marginalize_continuous(
      outcome_model,
      mc_data,
      treatment,
      iv,
      treatment_model,
      n_mc
    )
  }
}


#' Exact marginalization for binary treatment
#'
#' Enumerates both treatment values weighted by the propensity:
#' \deqn{\hat\mu_i = \hat p_i \cdot \hat m(d(1, L_i), L_i) +
#'   (1 - \hat p_i) \cdot \hat m(d(0, L_i), L_i)}
#'
#' @param outcome_model Fitted outcome model (from g-comp).
#' @param mc_data data.table of target rows to marginalize over.
#' @param treatment Character. Treatment column name.
#' @param iv A `causatr_intervention` object or `NULL` (natural course).
#' @param treatment_model Fitted binomial GLM for \eqn{P(A=1 \mid L, S=1)}.
#'
#' @return Numeric vector of marginalized predictions (one per target row).
#'
#' @noRd
mc_marginalize_binary <- function(
  outcome_model,
  mc_data,
  treatment,
  iv,
  treatment_model
) {
  p1 <- stats::predict(treatment_model, newdata = mc_data, type = "response")

  # A = 0 branch
  d0 <- data.table::copy(mc_data)
  d0[, (treatment) := 0L]
  if (!is.null(iv)) {
    d0 <- apply_single_intervention(d0, treatment, iv)
  }
  preds_0 <- stats::predict(outcome_model, newdata = d0, type = "response")

  # A = 1 branch
  d1 <- data.table::copy(mc_data)
  d1[, (treatment) := 1L]
  if (!is.null(iv)) {
    d1 <- apply_single_intervention(d1, treatment, iv)
  }
  preds_1 <- stats::predict(outcome_model, newdata = d1, type = "response")

  (1 - p1) * preds_0 + p1 * preds_1
}


#' MC marginalization for continuous treatment
#'
#' Draws \eqn{A_1, \ldots, A_M \sim f(A \mid L, S=1)} from the fitted
#' treatment model, applies the intervention, predicts, and averages.
#'
#' @param outcome_model Fitted outcome model (from g-comp).
#' @param mc_data data.table of target rows to marginalize over.
#' @param treatment Character. Treatment column name.
#' @param iv A `causatr_intervention` object or `NULL` (natural course).
#' @param treatment_model Fitted gaussian GLM for \eqn{E[A \mid L, S=1]}.
#' @param n_mc Integer. Number of Monte Carlo draws per target row.
#'
#' @return Numeric vector of marginalized predictions (one per target row).
#'
#' @noRd
mc_marginalize_continuous <- function(
  outcome_model,
  mc_data,
  treatment,
  iv,
  treatment_model,
  n_mc
) {
  mu_a <- stats::predict(treatment_model, newdata = mc_data, type = "response")
  sigma_a <- sqrt(
    sum(stats::residuals(treatment_model, type = "response")^2) /
      treatment_model$df.residual
  )

  n_rows <- nrow(mc_data)
  preds_sum <- rep(0, n_rows)

  draw_data <- data.table::copy(mc_data)
  for (m in seq_len(n_mc)) {
    a_draw <- stats::rnorm(n_rows, mean = mu_a, sd = sigma_a)
    data.table::set(draw_data, j = treatment, value = a_draw)
    if (!is.null(iv)) {
      draw_data <- apply_single_intervention(draw_data, treatment, iv)
    }
    preds_m <- stats::predict(
      outcome_model,
      newdata = draw_data,
      type = "response"
    )
    preds_sum <- preds_sum + preds_m
  }

  preds_sum / n_mc
}


#' Core standardisation engine for causal contrasts
#'
#' @description
#' Implements the g-formula standardisation algorithm (Hernan & Robins Ch. 13):
#' for each named intervention, sets each individual's treatment to the
#' intervened value (via `apply_intervention()`), predicts outcomes from the
#' fitted model, averages over the target population, then computes pairwise
#' contrasts with uncertainty estimates.
#'
#' For each intervention a, computes \eqn{E[Y^a]} by averaging model
#' predictions over the target population rows.
#'
#' @param fit A `causatr_fit` object.
#' @param interventions Named list of interventions (or `NULL`).
#' @param type Contrast scale: `"difference"`, `"ratio"`, `"or"`.
#' @param estimand `"ATE"`, `"ATT"`, `"ATC"`, or `NULL`.
#' @param subset Quoted expression or `NULL`.
#' @param reference Name of the reference intervention or `NULL`.
#' @param ci_method `"sandwich"` (IF variance) or `"bootstrap"`.
#' @param n_boot Bootstrap replications (for `"bootstrap"`).
#' @param conf_level Confidence level (e.g. 0.95).
#' @param by Stratification variable or `NULL`.
#' @param call The original `contrast()` call.
#'
#' @return A `causatr_result` object.
#'
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
  parallel,
  ncpus,
  call,
  subset_env = parent.frame(),
  cluster_vec = NULL
) {
  data <- fit$data
  int_names <- names(interventions)

  # Resolve the effective estimand: if the caller passed `estimand =`
  # to `contrast()`, use that override (g-comp only -- the check above
  # already blocked IPW/matching). Otherwise fall back to whatever was
  # recorded at fitting time. This is what lets a single g-comp fit
  # produce ATE, ATT, ATC, and subgroup effects from one model.
  est <- if (!is.null(estimand)) estimand else fit$estimand

  # -- `by` branch: effect modification.
  # When `by = "sex"` is given, we recursively call compute_contrast()
  # once per level and stitch the results into one combined table.
  # This isn't just a user-facing convenience -- it's the only way to
  # get the variance engine to re-define the target population per
  # level, since vcov is computed conditional on the target.
  if (!is.null(by)) {
    # B8 (2026-04-15 review): enumerate levels from the *fit rows*
    # (rows the outcome / propensity / match model was actually fit on),
    # not the full dataset. Levels that appear only in censored / NA-
    # outcome rows produced empty target populations downstream, and
    # the inner `compute_contrast()` then aborted the whole call with
    # "target population is empty". With ATT/ATC layered on top, a
    # level with no treated units also killed every other stratum.
    #
    # Fall back to the full data for longitudinal (ICE) fits because
    # `fit_rows` there refers to intermediate backward-iteration rows,
    # not a single outcome-model fit -- the full person-period frame is
    # the safer source for level enumeration.
    fit_rows_for_by <- fit$details$fit_rows
    if (
      !is.null(fit_rows_for_by) &&
        is.logical(fit_rows_for_by) &&
        length(fit_rows_for_by) == nrow(data)
    ) {
      by_source <- data[fit_rows_for_by][[by]]
    } else {
      by_source <- data[[by]]
    }
    by_vals <- sort(unique(stats::na.omit(by_source)))
    by_sym <- as.name(by)

    # If the caller *also* asked for ATT/ATC, that adds a
    # treatment-value restriction on top of the by-level restriction.
    # We build an `est_subset` expression (e.g. `A == 1` for ATT) and
    # AND it into every level's combined subset below.
    est_subset <- NULL
    if (est %in% c("ATT", "ATC")) {
      trt_sym <- as.name(fit$treatment[1])
      trt_val <- if (est == "ATT") 1L else 0L
      est_subset <- bquote(.(trt_sym) == .(trt_val))
    }

    # Per-level compute_contrast() call. We pass `estimand = NULL` and
    # a *combined* subset expression -- this collapses the ATT/ATC
    # selection and the by-level selection into a single quoted
    # predicate, so the inner call treats it as a subgroup request.
    #
    # B8: wrap the inner call in tryCatch and skip levels whose target
    # population turns out to be empty (e.g. a `by` level with no
    # treated units under ATT). Without this, one unbalanced stratum
    # killed every other level's estimate. Emit a single warning
    # listing the skipped levels once at the end.
    skipped <- character()
    results_list <- lapply(by_vals, function(lev) {
      by_subset <- bquote(.(by_sym) == .(lev))
      combined_subset <- if (!is.null(subset)) {
        bquote(.(subset) & .(by_subset))
      } else {
        by_subset
      }
      if (!is.null(est_subset)) {
        combined_subset <- bquote(.(combined_subset) & .(est_subset))
      }
      tryCatch(
        compute_contrast(
          fit,
          interventions,
          type,
          estimand = NULL,
          subset = combined_subset,
          reference,
          ci_method,
          n_boot,
          conf_level,
          by = NULL,
          parallel,
          ncpus,
          call,
          subset_env = subset_env,
          cluster_vec = cluster_vec
        ),
        error = function(e) {
          # Match on the classed abort from build_point_channel_pieces()
          # (C3 second-round review). Class-based matching survives
          # any future wording drift in the error message.
          if (inherits(e, "causatr_empty_target")) {
            skipped <<- c(skipped, as.character(lev))
            return(NULL)
          }
          rlang::abort(conditionMessage(e), parent = e)
        }
      )
    })
    names(results_list) <- as.character(by_vals)
    # Drop skipped levels from the stitched output.
    keep <- !vapply(results_list, is.null, logical(1))
    results_list <- results_list[keep]
    if (length(skipped) > 0L) {
      rlang::warn(
        paste0(
          "Skipped `by` level(s) with empty target population: ",
          paste(skipped, collapse = ", "),
          ". This typically means the level has no rows satisfying the ",
          "estimand (e.g. no treated units under ATT)."
        )
      )
    }
    if (length(results_list) == 0L) {
      rlang::abort(
        "All `by` levels have empty target populations; no effect estimates."
      )
    }

    # Stitch per-level tables: each level contributes its `estimates`
    # and `contrasts` data.tables, augmented with the by-level label
    # and its target-population size. `data.table::copy()` is required
    # because we're mutating each table in place.
    est_list <- lapply(names(results_list), function(lev) {
      dt <- data.table::copy(results_list[[lev]]$estimates)
      dt[, by := lev]
      dt[, n_by := results_list[[lev]]$n]
      dt
    })
    con_list <- lapply(names(results_list), function(lev) {
      dt <- data.table::copy(results_list[[lev]]$contrasts)
      dt[, by := lev]
      dt[, n_by := results_list[[lev]]$n]
      dt
    })

    combined_est <- data.table::rbindlist(est_list)
    combined_con <- data.table::rbindlist(con_list)

    return(
      new_causatr_result(
        estimates = combined_est,
        contrasts = combined_con,
        type = type,
        estimand = if (!is.null(subset)) "subset" else est,
        ci_method = ci_method,
        reference = if (!is.null(reference)) reference else int_names[1],
        interventions = interventions,
        n = sum(vapply(
          results_list,
          function(r) r$n,
          integer(1)
        )),
        estimator = fit$estimator,
        family = fit$family,
        fit_type = fit$type,
        vcov = lapply(results_list, function(r) r$vcov),
        boot_t = lapply(results_list, function(r) r$boot_t),
        boot_info = lapply(results_list, function(r) r$boot_info),
        call = call
      )
    )
  }

  # Compute marginal means + variance.
  #
  # Both point-treatment and longitudinal (ICE) paths produce the same three
  # outputs: mu_hat (named vector of marginal means), vcov_mat (k x k vcov
  # matrix), and n_target (number of individuals averaged over).
  #
  # Everything downstream (SEs, estimates table, pairwise contrasts, result
  # assembly) is shared between point and longitudinal g-computation.
  #
  # Bundles are computed per intervention and then diffed (rather than
  # predicting a single contrast directly) because each intervention induces
  # its own density-ratio weights (IPW), pseudo-outcome sequence (ICE), or
  # augmentation term (AIPW). The variance engine must see the full k-vector
  # of mu_hats to compute the joint vcov -- which is required for correct
  # contrast SEs that account for the covariance between mu_a1 and mu_a0.

  if (fit$type == "longitudinal") {
    # -- Longitudinal data: dispatch on estimator.
    # The target population for longitudinal data is defined at the
    # FIRST time point: every unique individual shows up once at
    # baseline, so that's where we enumerate the target. Subsetting is
    # evaluated against `data` but then collapsed to a baseline-only
    # logical vector (`target_within_first`) because that's what the
    # downstream variance engine (ICE or longitudinal IPW) expects.
    time_col <- fit$time
    first_time <- fit$details$time_points[1]
    rows_first <- data[[time_col]] == first_time

    if (!is.null(subset)) {
      target_baseline <- rows_first &
        as.logical(eval(subset, envir = data, enclos = subset_env))
    } else {
      target_baseline <- rows_first
    }
    target_within_first <- target_baseline[rows_first]
    n_target <- sum(target_within_first)

    boot_t <- NULL
    boot_info <- NULL

    if (fit$estimator == "ipw") {
      # Longitudinal IPW. Per-intervention bundles refit
      # the final-period weighted Hajek MSM under the cumulative
      # density-ratio weight. Variance dispatches through
      # `variance_if_ipw_longitudinal()` (sandwich) or
      # `ipw_longitudinal_variance_bootstrap()` (bootstrap, id-
      # clustered).
      ipw_long <- compute_ipw_contrast_longitudinal(
        fit,
        interventions,
        target_within_first
      )
      mu_hat <- ipw_long$mu_hat

      if (ci_method == "sandwich") {
        vcov_mat <- variance_if(
          fit,
          ipw_bundles = ipw_long$bundles,
          target_within_first = target_within_first,
          cluster_vec = cluster_vec
        )
      } else {
        boot_res <- ipw_longitudinal_variance_bootstrap(
          fit,
          interventions,
          n_boot,
          target_within_first,
          est,
          subset,
          parallel,
          ncpus,
          subset_env = subset_env
        )
        vcov_mat <- boot_res$vcov
        boot_t <- boot_res$boot_t
        boot_info <- boot_res$boot_info
      }
    } else if (fit$estimator == "aipw") {
      # Longitudinal AIPW (ICE-AIPW, Bang & Robins 2005). Composes
      # ICE g-computation with density-ratio weighting, augmenting
      # inside the backward loop at each time step.
      aipw_long <- compute_aipw_contrast_longitudinal(
        fit,
        interventions,
        target_within_first
      )
      ice_aipw_results <- aipw_long$ice_aipw_results
      mu_hat <- aipw_long$mu_hat

      if (ci_method == "sandwich") {
        vcov_mat <- variance_if(
          fit,
          ice_aipw_results = ice_aipw_results,
          target_within_first = target_within_first,
          cluster_vec = cluster_vec
        )
      } else {
        boot_res <- aipw_longitudinal_variance_bootstrap(
          fit,
          interventions,
          n_boot,
          target_within_first,
          est,
          subset,
          parallel,
          ncpus,
          subset_env = subset_env
        )
        vcov_mat <- boot_res$vcov
        boot_t <- boot_res$boot_t
        boot_info <- boot_res$boot_info
      }
    } else {
      # ICE g-computation (default longitudinal path). Each call
      # returns per-individual pseudo-outcomes at baseline along with
      # the fitted chain of models the sandwich variance engine
      # consumes.

      # For transport, override target_within_first with the transport
      # target population (S=0 for transportability, all for
      # generalizability). Point IPW/AIPW use sampling weights instead;
      # gcomp standardises the pseudo-outcome mean over the target.
      if (!is.null(fit$target)) {
        transport_rows_ice <- transport_target_idx(
          data,
          fit$target,
          fit$target_subset
        )
        target_within_first <- transport_rows_ice[rows_first]
        n_target <- sum(target_within_first)
      }

      ice_results <- stats::setNames(
        lapply(interventions, function(iv) ice_iterate(fit, iv)),
        int_names
      )

      ext_w <- fit$details$weights
      if (!is.null(ext_w)) {
        w_first <- ext_w[rows_first]
        w_target_ice <- w_first[target_within_first]
      } else {
        w_target_ice <- NULL
      }
      mu_hat <- vapply(
        ice_results,
        function(res) {
          maybe_weighted_mean(
            res$pseudo_final[target_within_first],
            w_target_ice
          )
        },
        numeric(1)
      )
      names(mu_hat) <- int_names

      if (ci_method == "sandwich") {
        vcov_mat <- variance_if(
          fit,
          ice_results = ice_results,
          target_within_first = target_within_first,
          cluster_vec = cluster_vec
        )
      } else {
        boot_res <- ice_variance_bootstrap(
          fit,
          interventions,
          n_boot,
          target_within_first,
          est,
          subset,
          parallel,
          ncpus,
          subset_env = subset_env
        )
        vcov_mat <- boot_res$vcov
        boot_t <- boot_res$boot_t
        boot_info <- boot_res$boot_info
      }
    }
  } else if (fit$estimator == "ipw") {
    # Point-treatment IPW. The density-ratio weights depend on the
    # intervention, so `compute_ipw_contrast_point()` refits a weighted
    # MSM per intervention and returns one bundle per intervention.
    # Variance dispatch then goes through `variance_if()` (sandwich)
    # or `variance_bootstrap()` (bootstrap), uniformly with the other
    # estimators.
    #
    # Hajek estimand: the MSM is Y ~ 1 (intercept-only), not Y ~ A.
    # The treatment effect is fully absorbed into the density-ratio weight
    # w_i = f(d(A_i) | L_i) / f(A_i | L_i). The intercept of the weighted
    # regression then directly estimates E[Y^a] = (sum_i w_i Y_i) /
    # (sum_i w_i) (Hajek 1971). Using Y ~ A instead would double-count
    # the treatment by conditioning on both the weight and the A column.
    target_idx <- get_target_idx(data, fit$treatment, est, subset, subset_env)
    target_idx <- transport_target_idx(data, fit$target, fit$target_subset) %||%
      target_idx
    n_target <- sum(target_idx)
    if (n_target == 0L) {
      rlang::abort(
        "compute_contrast(): target population is empty.",
        class = "causatr_empty_target"
      )
    }

    ipw_point <- compute_ipw_contrast_point(fit, interventions, target_idx)
    mu_hat <- ipw_point$mu_hat

    boot_t <- NULL
    boot_info <- NULL
    if (ci_method == "sandwich") {
      vcov_mat <- variance_if(
        fit,
        target_idx = target_idx,
        mu_hat = mu_hat,
        ipw_bundles = ipw_point$bundles,
        ipw_fit_idx = ipw_point$fit_idx,
        ipw_n_total = ipw_point$n_total,
        cluster_vec = cluster_vec
      )
    } else {
      boot_res <- variance_bootstrap(
        fit,
        interventions,
        n_boot,
        target_idx,
        est,
        subset,
        parallel,
        ncpus,
        subset_env = subset_env
      )
      vcov_mat <- boot_res$vcov
      boot_t <- boot_res$boot_t
      boot_info <- boot_res$boot_info
    }
  } else if (fit$estimator == "aipw") {
    # Point-treatment AIPW. Combines the outcome-model predictions
    # (g-comp piece) with density-ratio weighted residuals (IPW
    # augmentation). Both nuisance models were fit at causat() time;
    # the contrast assembles the doubly-robust functional per
    # intervention.
    target_idx <- get_target_idx(data, fit$treatment, est, subset, subset_env)
    target_idx <- transport_target_idx(data, fit$target, fit$target_subset) %||%
      target_idx
    n_target <- sum(target_idx)
    if (n_target == 0L) {
      rlang::abort(
        "compute_contrast(): target population is empty.",
        class = "causatr_empty_target"
      )
    }

    aipw_point <- compute_aipw_contrast_point(fit, interventions, target_idx)
    mu_hat <- aipw_point$mu_hat

    boot_t <- NULL
    boot_info <- NULL

    # Reject sandwich for MTP + transport: the MC marginalization
    # introduces treatment-model dependence not captured by the IF.
    is_transport_aipw <- isTRUE(fit$details$transport)
    any_mc_aipw <- is_transport_aipw &&
      any(vapply(
        interventions,
        needs_observed_treatment,
        logical(1)
      ))
    if (ci_method == "sandwich" && any_mc_aipw) {
      rlang::abort(
        c(
          paste0(
            "Sandwich variance is not supported for MTP interventions ",
            "(shift, scale_by, threshold, dynamic) combined with ",
            "transportability under AIPW."
          ),
          i = paste0(
            "The MC marginalization over the study treatment distribution ",
            "introduces dependence on the treatment model whose influence ",
            "function is not yet propagated through Term 1."
          ),
          i = "Use `ci_method = \"bootstrap\"` for valid inference."
        ),
        class = "causatr_mtp_transport_sandwich"
      )
    }

    if (ci_method == "sandwich") {
      vcov_mat <- variance_if(
        fit,
        target_idx = target_idx,
        mu_hat = mu_hat,
        aipw_bundles = aipw_point$bundles,
        aipw_fit_idx = aipw_point$fit_idx,
        aipw_n_total = aipw_point$n_total,
        cluster_vec = cluster_vec
      )
    } else {
      boot_res <- variance_bootstrap(
        fit,
        interventions,
        n_boot,
        target_idx,
        est,
        subset,
        parallel,
        ncpus,
        subset_env = subset_env
      )
      vcov_mat <- boot_res$vcov
      boot_t <- boot_res$boot_t
      boot_info <- boot_res$boot_info
    }
  } else {
    # -- Point-treatment g-computation / matching.
    # Single outcome model, predict once per intervention, average
    # over the target population. Matching reuses the same
    # predict-then-average path because its matched-sample outcome
    # model is already the marginal structural model, so marginal-
    # mean predictions read off the MSM directly.
    model <- fit$model

    # Logical vector (length n) flagging the target population.
    # ATE: all n rows (rep(TRUE, n)).
    # ATT: rows where A == 1 at baseline (standardise over the treated).
    # ATC: rows where A == 0 at baseline (standardise over the controls).
    # subset: eval(subset_expr, data) -- overrides the estimand keyword.
    # This vector is passed to the variance engine so the IF aggregation
    # averages over the same n_target rows that mu_hat was computed on.
    target_idx <- get_target_idx(data, fit$treatment, est, subset, subset_env)
    target_idx <- transport_target_idx(data, fit$target, fit$target_subset) %||%
      target_idx

    # Predict E[Y | A = a(L_i), L_i] under each intervention. For
    # deterministic interventions this is a single apply + predict; for
    # stochastic interventions, MC-average over n_mc draws.
    pred_results <- lapply(interventions, function(iv) {
      predict_under_intervention(model, data, fit$treatment, iv)
    })
    preds_list <- lapply(pred_results, `[[`, "preds")
    data_a_list <- lapply(pred_results, `[[`, "data_a")

    # MTP + transport: MC-marginalize predictions for target rows
    # where treatment is unobserved. For each MTP intervention,
    # replace the NA predictions on target rows with
    # E_{A|L,S=1}[m(d(A,L), L)] computed via exact enumeration
    # (binary) or Monte Carlo integration (continuous).
    is_transport <- isTRUE(fit$details$transport)
    any_needs_mc <- is_transport &&
      any(vapply(
        interventions,
        needs_observed_treatment,
        logical(1)
      ))
    if (any_needs_mc && length(fit$treatment) == 1L) {
      mc_rows <- is.na(data[[fit$treatment]])
      if (any(mc_rows)) {
        mc_tm <- fit_mc_treatment_model(
          data,
          fit$treatment,
          resolve_confounders_treatment(fit),
          fit$details$fit_rows
        )
        mc_data <- data[mc_rows]
        for (i in seq_along(interventions)) {
          if (needs_observed_treatment(interventions[[i]])) {
            preds_list[[i]][mc_rows] <- mc_marginalize_preds(
              model,
              mc_data,
              fit$treatment,
              interventions[[i]],
              mc_tm
            )
          }
        }
      }
    }

    # Handle NA predictions (e.g. rows with missing confounders).
    # Intersect a "valid-across-all-interventions" mask with the target
    # to avoid losing rows whose prediction is NA under *any* regime
    # we care about. Inform (not warn) if the drop is nontrivial:
    # for the canonical NHEFS workflow this fires every time because
    # 117 rows have missing education, and the book accepts the
    # exclusion as a data-hygiene fact rather than a problem. Keeping
    # it as a `warn()` made every NHEFS-using test surface a noisy
    # WARN line. `inform()` is the right semantic level -- visible to
    # users but not flagged by testthat as a real warning.
    valid_preds <- Reduce(
      `&`,
      lapply(preds_list, function(p) !is.na(p)),
      init = rep(TRUE, nrow(data))
    )
    n_dropped <- sum(!valid_preds & target_idx)
    if (n_dropped > 0L) {
      rlang::inform(
        paste0(
          n_dropped,
          " row(s) with NA predictions excluded from the ",
          "target population."
        )
      )
    }
    target_idx <- target_idx & valid_preds
    n_target <- sum(target_idx)

    # Marginal mean: weighted or plain average over target rows.
    ext_w <- fit$details$weights
    w_target <- if (!is.null(ext_w)) ext_w[target_idx] else NULL
    mu_hat <- vapply(
      preds_list,
      function(p) maybe_weighted_mean(p[target_idx], w_target),
      numeric(1)
    )
    names(mu_hat) <- int_names

    # Variance: sandwich via the IF engine (shared across gcomp / ipw
    # / matching -- the engine's dispatcher picks the right branch from
    # `fit$estimator`) or nonparametric bootstrap.
    boot_t <- NULL
    boot_info <- NULL
    if (ci_method == "sandwich" && any_needs_mc) {
      rlang::abort(
        c(
          paste0(
            "Sandwich variance is not supported for MTP interventions ",
            "(shift, scale_by, threshold, dynamic) combined with ",
            "transportability under g-computation."
          ),
          i = paste0(
            "The MC marginalization over the study treatment distribution ",
            "introduces dependence on a treatment model whose influence ",
            "function is not yet accounted for in the sandwich."
          ),
          i = "Use `ci_method = \"bootstrap\"` for valid inference."
        ),
        class = "causatr_mtp_transport_sandwich"
      )
    }
    if (ci_method == "sandwich") {
      vcov_mat <- variance_if(
        fit,
        interventions = interventions,
        data_a_list = data_a_list,
        preds_list = preds_list,
        mu_hat = mu_hat,
        target_idx = target_idx,
        cluster_vec = cluster_vec
      )
    } else {
      boot_res <- variance_bootstrap(
        fit,
        interventions,
        n_boot,
        target_idx,
        est,
        subset,
        parallel,
        ncpus,
        subset_env = subset_env
      )
      vcov_mat <- boot_res$vcov
      boot_t <- boot_res$boot_t
      boot_info <- boot_res$boot_info
    }
  }

  rownames(vcov_mat) <- int_names
  colnames(vcov_mat) <- int_names

  # Marginal-mean SEs come from the diagonal of vcov. `pmax(., 0)`
  # guards against tiny negative values from floating-point roundoff
  # in the IF aggregation -- they'd otherwise become NaN under sqrt.
  se_means <- sqrt(pmax(diag(vcov_mat), 0))
  # Two-sided normal critical value. Using qnorm(0.5 + level/2) gives
  # the right half-width for any conf_level without hand-coding 1.96.
  z <- stats::qnorm((1 + conf_level) / 2)

  # First output: per-intervention marginal-mean table.
  estimates_dt <- data.table::data.table(
    intervention = int_names,
    estimate = mu_hat,
    se = se_means,
    ci_lower = mu_hat - z * se_means,
    ci_upper = mu_hat + z * se_means
  )

  # Reference for pairwise contrasts. If the user didn't name one,
  # default to the first intervention in list order -- matches how
  # users conventionally write `list(treat, control)` with the
  # control as the second element.
  ref_name <- if (!is.null(reference)) reference else int_names[1]
  non_ref <- setdiff(int_names, ref_name)
  mu_ref <- mu_hat[ref_name]
  idx_ref <- which(int_names == ref_name)

  # Tolerance for boundary checks on marginal means. A predicted mean
  # this close to 0 or 1 makes ratios / odds ratios numerically
  # unstable (the log-scale delta method divides by mu_ref and
  # 1 - mu_ref), so we refuse to compute rather than return Inf / NaN.
  tol_edge <- sqrt(.Machine$double.eps)

  if (type == "ratio" && abs(mu_ref) < tol_edge) {
    rlang::abort(paste0(
      "Reference intervention '",
      ref_name,
      "' has a marginal mean of ",
      mu_ref,
      ". The risk/mean ratio is undefined."
    ))
  }
  if (type == "or" && (abs(mu_ref) < tol_edge || abs(1 - mu_ref) < tol_edge)) {
    rlang::abort(paste0(
      "Reference intervention '",
      ref_name,
      "' has a marginal mean of ",
      mu_ref,
      ". The odds ratio is undefined when the probability is 0 or 1."
    ))
  }
  # Ratios and odds ratios use a log-scale delta method CI (see below).
  # log() of a negative or zero ratio is NaN / -Inf, so we must reject
  # any comparison where the reference and an alternative have opposite
  # signs -- legal for Gaussian outcomes where mu can be negative, fatal
  # for the exp(log(R) +/- z se_log) CI construction. Odds ratios inherit
  # the same concern via the OR formula's mu/(1-mu) term, which can
  # cross zero for Gaussian outcomes. Abort with a clear message
  # pointing to `type = "difference"` or a binomial/poisson/gamma
  # refit so users are not handed silent NaN CIs.
  if (type %in% c("ratio", "or")) {
    neg_mu <- mu_hat[!is.na(mu_hat)] <= 0
    if (any(neg_mu)) {
      bad <- names(mu_hat)[which(mu_hat <= 0)]
      rlang::abort(paste0(
        "type = '",
        type,
        "' requires all intervention-specific marginal means to be ",
        "strictly positive, but '",
        paste(bad, collapse = "', '"),
        "' ",
        if (length(bad) == 1L) "has" else "have",
        " a non-positive estimate (mu_hat = ",
        paste(round(mu_hat[bad], 4), collapse = ", "),
        "). The log-scale delta method CI is undefined for non-positive ",
        "ratios. Use `type = \"difference\"` for Gaussian outcomes, or ",
        "refit with a binomial / poisson / gamma family for ratios."
      ))
    }
    if (type == "or") {
      # R6 (2026-04-15 review): mirror the NA filter from the `<= 0`
      # branch above. Previously `any(mu_hat >= 1)` returned NA when any
      # mu_hat was NA, and `if (NA)` aborted with "missing value where
      # TRUE/FALSE needed" -- surfacing a bogus OR-validation error on
      # a NA-prediction that the target population had already excluded.
      gt1_mu <- mu_hat[!is.na(mu_hat)] >= 1
      if (any(gt1_mu)) {
        bad <- names(mu_hat)[which(mu_hat >= 1)]
        rlang::abort(paste0(
          "type = 'or' requires all intervention-specific marginal means ",
          "to lie strictly in (0, 1) (they are probabilities), but '",
          paste(bad, collapse = "', '"),
          "' ",
          if (length(bad) == 1L) "has" else "have",
          " a mean >= 1 (mu_hat = ",
          paste(round(mu_hat[bad], 4), collapse = ", "),
          "). Refit with `family = \"binomial\"` / `\"quasibinomial\"` ",
          "or use `type = \"difference\"`."
        ))
      }
    }
  }

  # Pairwise contrasts a_j vs a_ref via the delta method on the vcov.
  # The vcov from `variance_if()` is on the (mu_1, mu_2, ..., mu_k)
  # scale, so for differences we can read the variance straight off;
  # for ratios / ORs we project through the appropriate gradient.
  contrasts_list <- lapply(non_ref, function(nm) {
    mu_a <- mu_hat[nm]
    idx_a <- which(int_names == nm)

    if (type == "difference") {
      # Difference: delta = mu_a - mu_ref.
      # Var(delta) = Var(mu_a) + Var(mu_ref) - 2 Cov(mu_a, mu_ref).
      # The cross term is the one dropped when people incorrectly add
      # per-intervention SEs in quadrature -- our IF engine keeps the
      # full covariance so we use the proper formula here.
      est_c <- mu_a - mu_ref
      var_c <- vcov_mat[idx_a, idx_a] +
        vcov_mat[idx_ref, idx_ref] -
        2 * vcov_mat[idx_a, idx_ref]
      se_c <- sqrt(max(var_c, 0))
      ci_lo <- est_c - z * se_c
      ci_hi <- est_c + z * se_c
    } else if (type == "ratio") {
      if (abs(mu_a) < tol_edge) {
        rlang::abort(paste0(
          "Intervention '",
          nm,
          "' has a marginal mean of ",
          mu_a,
          ". The risk/mean ratio is undefined (log-scale CI requires log(0))."
        ))
      }
      # Delta method for R := mu_a / mu_ref.
      # Gradient w.r.t. mu_a is 1/mu_ref;
      # w.r.t. mu_ref is -mu_a / mu_ref^2.
      # Linear-scale SE from grad^T V grad, then convert to log-scale
      # CI: log(R) +/- z * se_log, se_log = se / R. Log-scale CIs respect
      # the (0, Inf) support of a ratio and have better coverage than
      # Wald CIs, which can produce negative lower bounds.
      est_c <- mu_a / mu_ref
      grad <- c(1 / mu_ref, -mu_a / mu_ref^2)
      sub_v <- vcov_mat[c(idx_a, idx_ref), c(idx_a, idx_ref)]
      se_c <- sqrt(max(as.numeric(t(grad) %*% sub_v %*% grad), 0))
      se_log <- se_c / est_c
      ci_lo <- exp(log(est_c) - z * se_log)
      ci_hi <- exp(log(est_c) + z * se_log)
    } else {
      if (abs(mu_a) < tol_edge || abs(1 - mu_a) < tol_edge) {
        rlang::abort(paste0(
          "Intervention '",
          nm,
          "' has a marginal mean of ",
          mu_a,
          ". The odds ratio is undefined when the probability is 0 or 1."
        ))
      }
      # Odds ratio:
      #   OR = [mu_a/(1 - mu_a)] / [mu_ref/(1 - mu_ref)]
      # Delta method gradients:
      # Gradient w.r.t. mu_a is OR / (mu_a * (1 - mu_a));
      # w.r.t. mu_ref is -OR / (mu_ref * (1 - mu_ref)).
      # These come from differentiating log(OR) w.r.t. each mu and
      # multiplying by OR. Same log-scale CI pattern as ratios.
      est_c <- (mu_a / (1 - mu_a)) / (mu_ref / (1 - mu_ref))
      grad <- c(
        est_c / (mu_a * (1 - mu_a)),
        -est_c / (mu_ref * (1 - mu_ref))
      )
      sub_v <- vcov_mat[c(idx_a, idx_ref), c(idx_a, idx_ref)]
      se_c <- sqrt(max(as.numeric(t(grad) %*% sub_v %*% grad), 0))
      se_log <- se_c / est_c
      ci_lo <- exp(log(est_c) - z * se_log)
      ci_hi <- exp(log(est_c) + z * se_log)
    }

    data.table::data.table(
      comparison = paste0(nm, " vs ", ref_name),
      estimate = est_c,
      se = se_c,
      ci_lower = ci_lo,
      ci_upper = ci_hi
    )
  })

  contrasts_dt <- if (length(contrasts_list) > 0) {
    data.table::rbindlist(contrasts_list)
  } else {
    data.table::data.table(
      comparison = character(0),
      estimate = numeric(0),
      se = numeric(0),
      ci_lower = numeric(0),
      ci_upper = numeric(0)
    )
  }

  new_causatr_result(
    estimates = estimates_dt,
    contrasts = contrasts_dt,
    type = type,
    estimand = if (!is.null(subset)) "subset" else est,
    ci_method = ci_method,
    reference = ref_name,
    interventions = interventions,
    n = n_target,
    estimator = fit$estimator,
    family = fit$family,
    fit_type = fit$type,
    vcov = vcov_mat,
    boot_t = boot_t,
    boot_info = boot_info,
    call = call
  )
}

#' Determine which rows of `data` belong to the target population
#'
#' @description
#' Returns a logical vector of length `nrow(data)` indicating which rows
#' should be included when averaging predictions to estimate \eqn{E[Y^a]}.
#'
#' @param data A data.table.
#' @param treatment Character scalar or vector. Treatment column name(s).
#' @param estimand Character. `"ATE"` (all rows), `"ATT"` (treated rows where
#'   `A == 1`), or `"ATC"` (control rows where `A == 0`).
#' @param subset A quoted expression or `NULL`. When provided, overrides
#'   `estimand` and selects rows satisfying the expression evaluated in the
#'   context of `data`.
#' @param subset_env Environment in which to resolve free variables referenced
#'   by `subset` (e.g. a session-scoped `cutoff`). Must be captured at
#'   `contrast()`'s top level -- `parent.frame()` at this call site would be
#'   the internal dispatch frame, not the user's.
#'
#' @return Logical vector of length `nrow(data)`.
#'
#' @noRd
get_target_idx <- function(
  data,
  treatment,
  estimand,
  subset,
  subset_env = parent.frame()
) {
  # A quoted subset expression always takes priority over the estimand
  # keyword. `data` is a data.table (list-like), so `eval()` can resolve
  # column names directly; `subset_env` is the user's caller frame
  # threaded from `contrast()`, needed for `quote(age > cutoff)` to
  # resolve `cutoff` from the session.
  if (!is.null(subset)) {
    sel <- as.logical(eval(subset, envir = data, enclos = subset_env))
    if (length(sel) != nrow(data)) {
      rlang::abort(
        paste0(
          "`subset` expression must evaluate to a length-",
          nrow(data),
          " logical vector; got length ",
          length(sel),
          "."
        )
      )
    }
    return(sel)
  }
  if (estimand == "ATE") {
    # Average over all individuals in the dataset.
    return(rep(TRUE, nrow(data)))
  }
  # ATT and ATC are defined on the first (or only) treatment variable.
  trt_vals <- data[[treatment[1]]]
  if (estimand == "ATT") {
    return(!is.na(trt_vals) & trt_vals == 1)
  }
  !is.na(trt_vals) & trt_vals == 0
}


#' Apply transport target override to a target_idx logical vector
#'
#' @description
#' For transportability/generalizability, the target population is defined
#' by the sampling indicator S (not by estimand or treatment status).
#' Call this after `get_target_idx()` to override when `target` is non-NULL.
#'
#' This is a separate function (not folded into `get_target_idx()`) so that
#' `get_target_idx()` keeps its original signature — future-backend workers
#' load the installed package and would fail if signature changes.
#'
#' @param data A data.table.
#' @param target Column name of the sampling indicator, or `NULL`.
#' @param target_subset `"target"` (standardize over S=0) or `"all"` (all rows).
#'
#' @return A logical vector of length `nrow(data)`, or `NULL` if `target` is
#'   `NULL` (signals no transport override — caller keeps existing `target_idx`).
#' @noRd
transport_target_idx <- function(data, target, target_subset) {
  if (is.null(target)) {
    return(NULL)
  }
  if (identical(target_subset, "target")) {
    return(data[[target]] == 0L)
  }
  # target_subset == "all": generalizability — all rows are the target
  rep(TRUE, nrow(data))
}
