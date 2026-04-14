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
#' differ only in how the outcome model was fitted; variance is computed by the
#' **unified influence-function engine** `variance_if()` (`R/variance_if.R`),
#' which handles all four cases (g-comp, IPW, matching, ICE) via one entry
#' point. See `vignettes/variance-theory.qmd` for the derivation and the
#' `R/variance_if.R` roxygen header for the architecture.
#'
#' - `"gcomp"`: standard `glm`/`gam` on the full data; one-model IF
#'   correction via `prepare_model_if()` + `apply_model_correction()`.
#' - `"ipw"`: `glm_weightit()` fit weighted for the target estimand (ATE/ATT);
#'   IF Channel 2 delegated to `prepare_propensity_if()` (WeightIt shortcut via
#'   `sandwich::estfun(asympt = TRUE)`), which already accounts for weight
#'   estimation uncertainty.
#' - `"matching"`: `glm()` on the matched sample with match weights; IF is
#'   computed on the matched sample and aggregated cluster-robustly on
#'   `matched$subclass` via `vcov_from_if(cluster = ...)`.
#'
#' @param fit A `causatr_fit` object returned by [causat()].
#' @param interventions A named list of interventions. Each element must be
#'   one of:
#'   - A `causatr_intervention` object created by [static()], [shift()],
#'     [dynamic()], [scale_by()], [threshold()], or [ipsi()].
#'   - `NULL`, meaning the natural course (observed treatment values are used
#'     as-is). The natural course is the standard reference for modified
#'     treatment policies on continuous treatments (e.g. `shift(-10)` vs
#'     `NULL`; Díaz et al. 2023) and for longitudinal dynamic regimes
#'     (Hernán & Robins Ch. 21). For binary treatments, the
#'     natural-course marginal mean equals E\[Y\] under the observed
#'     treatment mechanism (by consistency); use `static(1)` vs
#'     `static(0)` for the conventional ATE. Supported for all methods.
#'   - A named list of `causatr_intervention` objects, one per treatment
#'     variable, for multivariate (joint) treatments. Each sub-list must name
#'     every treatment variable specified in `causat()`.
#'
#'   **Note:** Non-static interventions (`shift()`, `scale_by()`, `threshold()`,
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
#' @param ci_method Character. The variance/CI method: `"sandwich"` (default)
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
#'   (cross-platform). Passed directly to [boot::boot()]. Ignored when
#'   `ci_method = "sandwich"`.
#' @param ncpus Integer. Number of CPU cores for parallel bootstrap. Default
#'   `getOption("boot.ncpus", 1L)`. Passed directly to [boot::boot()].
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
#' E\[Y^a\] = (1/|S|) Σᵢ∈S Ê\[Y | A = a(Lᵢ), Lᵢ\]
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
#' - `"sandwich"`: The unified influence-function engine `variance_if()`.
#'   Per-individual IF = Channel 1 (direct covariate-sampling term) +
#'   Channel 2 (one correction per nuisance model). Aggregated via
#'   `vcov_from_if()`, with cluster-robust aggregation for matching.
#'   Mathematically equivalent to the stacked M-estimation sandwich; see
#'   `vignettes/variance-theory.qmd` Sections 3–4.
#' - `"bootstrap"`: Resamples individuals (entire pipeline refitted `n_boot`
#'   times). Respects cluster structure for longitudinal data.
#' - The delta method is applied internally for ratio / odds-ratio
#'   contrasts on top of the `"sandwich"` vcov; it is not a separate
#'   `ci_method` option.
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
  ci_method = c("sandwich", "bootstrap"),
  n_boot = 500L,
  conf_level = 0.95,
  by = NULL,
  parallel = c("no", "multicore", "snow"),
  ncpus = getOption("boot.ncpus", 1L)
) {
  # Capture the original call so that the returned `causatr_result`
  # can echo it in its `print` and `summary` methods. `match.call()`
  # here — not at compute_contrast() — so the recorded call reflects
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
  # the data inside `compute_contrast()` — this is how Hernán/Robins-
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
    check_estimand_compat(estimand, fit$method, fit$estimand)
    check_estimand_trt_compat(
      estimand,
      fit$treatment,
      fit$type,
      data = fit$data
    )
  }

  # Validate the interventions list — names, types, and method
  # compatibility (IPW/matching only accept static(), ICE accepts all).
  check_intervention_list(interventions)
  check_interventions_compat(fit$method, interventions)

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

  # `by` stratifies results by levels of a data column — e.g.
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

  parallel <- rlang::arg_match(parallel)

  # Hand off to the internal engine. Everything above is argument
  # validation; the actual math lives in `compute_contrast()` and
  # its method-specific delegates.
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
    call
  )
}

#' Check that interventions are compatible with the estimation method
#'
#' IPW and matching only support static interventions; non-static
#' interventions require g-computation.
#'
#' @param method Character estimation method.
#' @param interventions Named list of `causatr_intervention` objects.
#' @param call Caller environment for error messages.
#' @return `NULL` invisibly; aborts if incompatible interventions are found.
#' @noRd
check_interventions_compat <- function(
  method,
  interventions,
  call = rlang::caller_env()
) {
  # Only IPW and matching are restricted. G-comp supports every
  # intervention type because it re-predicts from an outcome model
  # that doesn't care about the treatment-assignment mechanism.
  if (method %in% c("ipw", "matching")) {
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
      # Phase 4 plan: self-contained IPW for non-static interventions
      # (see PHASE_4_INTERVENTIONS_SELF_IPW.md).
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

#' Core standardisation engine for causal contrasts
#'
#' @description
#' Implements the g-formula standardisation algorithm (Hernán & Robins Ch. 13):
#' for each named intervention, sets each individual's treatment to the
#' intervened value (via `apply_intervention()`), predicts outcomes from the
#' fitted model, averages over the target population, then computes pairwise
#' contrasts with uncertainty estimates.
#'
#' For each intervention a, computes \eqn{E[Y^a]} by averaging model
#' predictions over the target population rows.
#'
#' @param fit A `causatr_fit` object.
#' @param interventions Named list of `causatr_intervention` objects (or `NULL`).
#' @param type Contrast scale: `"difference"`, `"ratio"`, or `"or"`.
#' @param estimand Character or `NULL`. `"ATE"`, `"ATT"`, or `"ATC"`.
#' @param subset Quoted expression or `NULL`. Rows to average over.
#' @param reference Character or `NULL`. Name of the reference intervention.
#' @param ci_method Character. `"sandwich"` (IF variance via `variance_if()`)
#'   or `"bootstrap"`.
#' @param n_boot Integer. Bootstrap replications (for `ci_method = "bootstrap"`).
#' @param conf_level Numeric. Confidence level (e.g. 0.95).
#' @param by Character or `NULL`. Stratification variable for effect modification.
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
  call
) {
  data <- fit$data
  int_names <- names(interventions)

  # Resolve the effective estimand: if the caller passed `estimand =`
  # to `contrast()`, use that override (g-comp only — the check above
  # already blocked IPW/matching). Otherwise fall back to whatever was
  # recorded at fitting time. This is what lets a single g-comp fit
  # produce ATE, ATT, ATC, and subgroup effects from one model.
  est <- if (!is.null(estimand)) estimand else fit$estimand

  # ── `by` branch: effect modification.
  # When `by = "sex"` is given, we recursively call compute_contrast()
  # once per level and stitch the results into one combined table.
  # This isn't just a user-facing convenience — it's the only way to
  # get the variance engine to re-define the target population per
  # level, since vcov is computed conditional on the target.
  if (!is.null(by)) {
    by_vals <- sort(unique(stats::na.omit(data[[by]])))
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
    # a *combined* subset expression — this collapses the ATT/ATC
    # selection and the by-level selection into a single quoted
    # predicate, so the inner call treats it as a subgroup request.
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
        call
      )
    })
    names(results_list) <- as.character(by_vals)

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
        method = fit$method,
        family = fit$family,
        fit_type = fit$type,
        vcov = lapply(results_list, function(r) r$vcov),
        boot_t = lapply(results_list, function(r) r$boot_t),
        call = call
      )
    )
  }

  # Compute marginal means + variance.
  #
  # Both point-treatment and longitudinal (ICE) paths produce the same three
  # outputs: mu_hat (named vector of marginal means), vcov_mat (k × k vcov
  # matrix), and n_target (number of individuals averaged over).
  #
  # Everything downstream (SEs, estimates table, pairwise contrasts, result
  # assembly) is shared between point and longitudinal g-computation.

  if (fit$type == "survival") {
    rlang::abort(
      paste0(
        "Survival curve estimation via contrast() is not yet implemented. ",
        "causat_survival() currently fits a pooled logistic model; survival ",
        "curve contrasts are planned for a future release."
      ),
      .call = FALSE
    )
  }

  if (fit$type == "longitudinal") {
    # ── Longitudinal ICE g-computation.
    # The target population for longitudinal data is defined at the
    # FIRST time point: every unique individual shows up once at
    # baseline, so that's where we enumerate the target. Subsetting is
    # evaluated against `data` but then collapsed to a baseline-only
    # logical vector (`target_within_first`) because that's what the
    # downstream ICE variance engine expects.
    time_col <- fit$time
    first_time <- fit$details$time_points[1]
    rows_first <- data[[time_col]] == first_time

    if (!is.null(subset)) {
      # Evaluate the user's subset expression in the environment of
      # the full data.table, then AND with the baseline-row mask. This
      # lets users write things like `quote(age > 50)` and have it
      # interpreted over the first time point's covariates.
      target_baseline <- rows_first &
        as.logical(eval(subset, envir = as.list(data)))
    } else {
      target_baseline <- rows_first
    }
    # Collapse to a length-n_first logical — the per-individual target
    # mask used by `variance_if_ice_one()` and the mu_hat average.
    target_within_first <- target_baseline[rows_first]
    n_target <- sum(target_within_first)

    # Run ICE backward iteration once per intervention. Each call
    # returns the vector of individual pseudo-outcomes at baseline
    # (\hat Y^*_{0,i}) along with the fitted chain of models, which
    # the sandwich variance engine consumes directly.
    ice_results <- stats::setNames(
      lapply(interventions, function(iv) ice_iterate(fit, iv)),
      int_names
    )

    # Marginal mean: weighted average of pseudo-outcomes over the
    # target population. `maybe_weighted_mean()` is NA-safe and
    # handles the NULL-weights fall-through.
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

    # Variance: sandwich via the IF engine (default, fast) or
    # bootstrap via `ice_variance_bootstrap()` (clustered on `id`).
    boot_t <- NULL
    if (ci_method == "sandwich") {
      vcov_mat <- variance_if(
        fit,
        ice_results = ice_results,
        target_within_first = target_within_first
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
        ncpus
      )
      vcov_mat <- boot_res$vcov
      boot_t <- boot_res$boot_t
    }
  } else {
    # ── Point-treatment g-computation / IPW / matching.
    # Single outcome model, predict once per intervention, average
    # over the target population. This is the standard parametric
    # g-formula; IPW and matching reuse the same predict-then-average
    # path because their weighted / matched-sample outcome model is
    # already the marginal structural model, so marginal-mean
    # predictions read off the MSM directly.
    model <- fit$model

    # Logical vector (length n) flagging the target population.
    # Determined by `est` (ATE -> everyone; ATT -> treated at baseline;
    # ATC -> controls at baseline) and the optional `subset`.
    target_idx <- get_target_idx(data, fit$treatment, est, subset)

    # Build one counterfactual data.table per intervention. Each is a
    # copy of `data` with the treatment column(s) overwritten according
    # to the intervention rule (static, shift, dynamic, …).
    data_a_list <- lapply(interventions, function(iv) {
      apply_intervention(data, fit$treatment, iv)
    })

    # Predict E[Y | A = a(L_i), L_i] under each intervention. This is
    # the g-computation step: the outcome model was fit on observed A,
    # but we plug in the counterfactual A value row-by-row.
    preds_list <- lapply(data_a_list, function(da) {
      predict(model, newdata = da, type = "response")
    })

    # Handle NA predictions (e.g. rows with missing confounders).
    # Intersect a "valid-across-all-interventions" mask with the target
    # to avoid losing rows whose prediction is NA under *any* regime
    # we care about. Inform (not warn) if the drop is nontrivial:
    # for the canonical NHEFS workflow this fires every time because
    # 117 rows have missing education, and the book accepts the
    # exclusion as a data-hygiene fact rather than a problem. Keeping
    # it as a `warn()` made every NHEFS-using test surface a noisy
    # WARN line. `inform()` is the right semantic level — visible to
    # users but not flagged by testthat as a real warning.
    valid_preds <- Reduce(`&`, lapply(preds_list, function(p) !is.na(p)))
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
    # / matching — the engine's dispatcher picks the right branch from
    # `fit$method`) or nonparametric bootstrap.
    boot_t <- NULL
    if (ci_method == "sandwich") {
      vcov_mat <- variance_if(
        fit,
        interventions = interventions,
        data_a_list = data_a_list,
        preds_list = preds_list,
        mu_hat = mu_hat,
        target_idx = target_idx
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
        ncpus
      )
      vcov_mat <- boot_res$vcov
      boot_t <- boot_res$boot_t
    }
  }

  rownames(vcov_mat) <- int_names
  colnames(vcov_mat) <- int_names

  # Marginal-mean SEs come from the diagonal of vcov. `pmax(., 0)`
  # guards against tiny negative values from floating-point roundoff
  # in the IF aggregation — they'd otherwise become NaN under sqrt.
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
  # default to the first intervention in list order — matches how
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

  # Pairwise contrasts a_j vs a_ref via the delta method on the vcov.
  # The vcov from `variance_if()` is on the (mu_1, mu_2, ..., mu_k)
  # scale, so for differences we can read the variance straight off;
  # for ratios / ORs we project through the appropriate gradient.
  contrasts_list <- lapply(non_ref, function(nm) {
    mu_a <- mu_hat[nm]
    idx_a <- which(int_names == nm)

    if (type == "difference") {
      # Var(mu_a - mu_ref) = Var(mu_a) + Var(mu_ref) - 2 Cov(mu_a, mu_ref).
      # The cross term is the one dropped when people incorrectly add
      # per-intervention SEs in quadrature — our IF engine keeps the
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
      # Delta method for R = mu_a / mu_ref:
      #   dR/dmu_a    = 1/mu_ref
      #   dR/dmu_ref  = -mu_a/mu_ref^2
      # Linear-scale SE from grad^T V grad, then convert to log-scale
      # CI: log(R) ± z * se_log, se_log = se / R. Log-scale CIs respect
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
      #   dOR/dmu_a   = OR / (mu_a * (1 - mu_a))
      #   dOR/dmu_ref = -OR / (mu_ref * (1 - mu_ref))
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
    method = fit$method,
    family = fit$family,
    fit_type = fit$type,
    vcov = vcov_mat,
    boot_t = boot_t,
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
#'
#' @return Logical vector of length `nrow(data)`.
#'
#' @noRd
get_target_idx <- function(data, treatment, estimand, subset) {
  # A quoted subset expression always takes priority over the estimand keyword.
  if (!is.null(subset)) {
    return(as.logical(eval(subset, envir = as.list(data))))
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
