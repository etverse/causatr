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
  ci_method = c("sandwich", "bootstrap"),
  n_boot = 500L,
  conf_level = 0.95,
  by = NULL,
  parallel = c("no", "multicore", "snow"),
  ncpus = getOption("boot.ncpus", 1L)
) {
  call <- match.call()

  if (!inherits(fit, "causatr_fit")) {
    rlang::abort("`fit` must be a `causatr_fit` object returned by `causat()`.")
  }

  type <- rlang::arg_match(type)
  ci_method <- rlang::arg_match(ci_method)

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

  if (!is.null(estimand) && !is.null(subset)) {
    rlang::abort("Specify either 'estimand' or 'subset', not both.")
  }

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

  if (!is.null(by)) {
    check_string(by)
    if (!by %in% names(fit$data)) {
      rlang::abort(
        paste0("`by` variable '", by, "' not found in fitted data.")
      )
    }
  }

  parallel <- rlang::arg_match(parallel)

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
  if (method %in% c("ipw", "matching")) {
    non_static <- vapply(
      interventions,
      function(iv) {
        if (is.null(iv)) {
          return(FALSE)
        }
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
#' @param ci_method Character. `"sandwich"`, `"bootstrap"`, or `"delta"`.
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

  # Resolve estimand: use the caller's override when provided, otherwise fall
  # back to the one stored at fitting time.
  est <- if (!is.null(estimand)) estimand else fit$estimand

  if (!is.null(by)) {
    by_vals <- sort(unique(stats::na.omit(data[[by]])))
    by_sym <- as.name(by)

    est_subset <- NULL
    if (est %in% c("ATT", "ATC")) {
      trt_sym <- as.name(fit$treatment[1])
      trt_val <- if (est == "ATT") 1L else 0L
      est_subset <- bquote(.(trt_sym) == .(trt_val))
    }

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
    # Longitudinal ICE g-computation: run backward iteration per intervention.
    # Run the full backward iteration per intervention.  Target population
    # is defined over individuals at the first time point.
    time_col <- fit$time
    first_time <- fit$details$time_points[1]
    rows_first <- data[[time_col]] == first_time

    if (!is.null(subset)) {
      target_baseline <- rows_first &
        as.logical(eval(subset, envir = as.list(data)))
    } else {
      target_baseline <- rows_first
    }
    # target_within_first: logical vector over first-time rows only.
    target_within_first <- target_baseline[rows_first]
    n_target <- sum(target_within_first)

    # Run ICE backward iteration for each intervention.
    ice_results <- stats::setNames(
      lapply(interventions, function(iv) ice_iterate(fit, iv)),
      int_names
    )

    # Marginal mean = (weighted) average pseudo-outcomes at first time.
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
          res$pseudo_final[target_within_first], w_target_ice
        )
      },
      numeric(1)
    )
    names(mu_hat) <- int_names

    boot_t <- NULL
    if (ci_method == "sandwich") {
      vcov_mat <- ice_variance_sandwich(
        fit,
        ice_results,
        target_within_first
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
    # Point-treatment g-computation: single model, predict-then-average.
    # Single outcome model: predict under each intervention, average over
    # the target population (standard parametric g-formula).
    model <- fit$model

    # Logical vector (length n) that flags the target population rows.
    target_idx <- get_target_idx(data, fit$treatment, est, subset)

    # Create one counterfactual dataset per intervention:
    # copy `data` and set each individual's treatment to the intervened value.
    data_a_list <- lapply(interventions, function(iv) {
      apply_intervention(data, fit$treatment, iv)
    })

    # Predict E[Y | A = a(L_i), L_i] under each intervention.
    preds_list <- lapply(data_a_list, function(da) {
      predict(model, newdata = da, type = "response")
    })

    valid_preds <- Reduce(`&`, lapply(preds_list, function(p) !is.na(p)))
    n_dropped <- sum(!valid_preds & target_idx)
    if (n_dropped > 0L) {
      rlang::warn(
        paste0(
          n_dropped,
          " row(s) with NA predictions excluded from the ",
          "target population."
        )
      )
    }
    target_idx <- target_idx & valid_preds
    n_target <- sum(target_idx)

    # Marginal mean = (weighted) average predictions over target rows.
    ext_w <- fit$details$weights
    w_target <- if (!is.null(ext_w)) ext_w[target_idx] else NULL
    mu_hat <- vapply(
      preds_list,
      function(p) maybe_weighted_mean(p[target_idx], w_target),
      numeric(1)
    )
    names(mu_hat) <- int_names

    boot_t <- NULL
    if (ci_method == "sandwich") {
      vcov_mat <- variance_sandwich(
        fit,
        data_a_list,
        preds_list,
        target_idx
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

  # SE for each marginal mean from the diagonal of the vcov matrix.
  se_means <- sqrt(pmax(diag(vcov_mat), 0))
  z <- stats::qnorm((1 + conf_level) / 2)

  estimates_dt <- data.table::data.table(
    intervention = int_names,
    estimate = mu_hat,
    se = se_means,
    ci_lower = mu_hat - z * se_means,
    ci_upper = mu_hat + z * se_means
  )

  # Reference intervention defaults to the first one in the list.
  ref_name <- if (!is.null(reference)) reference else int_names[1]
  non_ref <- setdiff(int_names, ref_name)
  mu_ref <- mu_hat[ref_name]
  idx_ref <- which(int_names == ref_name)

  # Tolerance for boundary checks on marginal means. A predicted mean
  # this close to 0 or 1 makes ratios / odds ratios numerically unstable.
  tol_edge <- sqrt(.Machine$double.eps)

  if (type == "ratio" && abs(mu_ref) < tol_edge) {
    rlang::abort(paste0(
      "Reference intervention '", ref_name, "' has a marginal mean of ",
      mu_ref, ". The risk/mean ratio is undefined."
    ))
  }
  if (type == "or" && (abs(mu_ref) < tol_edge || abs(1 - mu_ref) < tol_edge)) {
    rlang::abort(paste0(
      "Reference intervention '", ref_name, "' has a marginal mean of ",
      mu_ref, ". The odds ratio is undefined when the probability is 0 or 1."
    ))
  }

  # Pairwise contrasts vs. the reference, with SE via delta method on the vcov.
  contrasts_list <- lapply(non_ref, function(nm) {
    mu_a <- mu_hat[nm]
    idx_a <- which(int_names == nm)

    if (type == "difference") {
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
          "Intervention '", nm, "' has a marginal mean of ", mu_a,
          ". The risk/mean ratio is undefined (log-scale CI requires log(0))."
        ))
      }
      # Delta method SE on the linear scale, then log-scale CI.
      # Log-scale CIs respect the (0, Inf) support and have better
      # coverage than Wald CIs, which can produce negative lower bounds.
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
          "Intervention '", nm, "' has a marginal mean of ", mu_a,
          ". The odds ratio is undefined when the probability is 0 or 1."
        ))
      }
      # OR = [mu_a/(1-mu_a)] / [mu_ref/(1-mu_ref)]. Log-scale CI.
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
