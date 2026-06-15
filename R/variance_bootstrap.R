#' Run a bootstrap via `boot::boot()` or via the `future` backend
#'
#' @description
#' Thin dispatcher between `boot::boot()`'s built-in parallelism
#' (`parallel = "no" / "multicore" / "snow"`) and the `future` backend
#' (`parallel = "future"`). The future path bypasses `boot::boot()`
#' entirely, draws `R` resample index vectors, evaluates `statistic(d,
#' idx)` via [future.apply::future_lapply()] so any active
#' [future::plan()] is respected, and assembles a minimal
#' `boot::boot`-compatible result (with slot `$t`) so
#' `process_boot_results()` consumes it uniformly.
#'
#' Using `future_lapply` rather than e.g. `boot::boot(parallel = "snow",
#' cl = future_cluster)` keeps us compatible with remote / cluster
#' / cloud futures the user has configured (e.g. `plan(cluster,
#' workers = <hostname_list>)` or `future.batchtools`). The tradeoff is
#' that we lose `boot::boot`'s own seed-aware PRNG, so the future path
#' seeds each worker via `future.seed = TRUE` to guarantee
#' reproducibility under `set.seed()` at the dispatch site.
#'
#' The `stype = "i"` in the stub is what tells `boot.ci()` and the
#' process function that `$t` holds per-replicate statistics keyed by
#' integer resample indices, matching `boot::boot(sim = "ordinary")`.
#'
#' @param data Data passed as the first argument to `statistic`.
#' @param statistic A `function(d, indices)` that returns a numeric
#'   vector of length `k` per replicate.
#' @param R Integer. Number of bootstrap replicates.
#' @param parallel Character. `"no"`, `"multicore"`, `"snow"`, or
#'   `"future"`.
#' @param ncpus Integer. Forwarded to `boot::boot()` for the non-future
#'   backends; ignored for `"future"` (the plan determines worker
#'   count).
#'
#' @return A `boot`-compatible list with at minimum `$t` (an `R x k`
#'   matrix of replicate statistics). Sufficient for
#'   `process_boot_results()`.
#'
#' @noRd
dispatch_boot <- function(data, statistic, R, parallel, ncpus) {
  if (parallel != "future") {
    return(boot::boot(
      data = data,
      statistic = statistic,
      R = R,
      parallel = parallel,
      ncpus = ncpus
    ))
  }

  check_pkg("future.apply")

  # Build per-replicate resample indices up front. For numeric-index
  # `boot_fn` (the shared point-treatment path) this is the full index
  # vector `1..n`; for the ICE cluster bootstrap the "data" is already
  # a vector of unique ids and `boot_fn` indexes into it. Either way,
  # `boot::boot` draws `length(data)` indices per replicate from
  # `1..length(data)` with replacement under `sim = "ordinary"`, so we
  # replicate that contract locally.
  n <- if (is.data.frame(data) || data.table::is.data.table(data)) {
    nrow(data)
  } else {
    length(data)
  }

  indices_list <- replicate(
    R,
    sample.int(n, n, replace = TRUE),
    simplify = FALSE
  )

  # `future.seed = TRUE` gives each worker a deterministic L'Ecuyer
  # stream so the `set.seed()` call at the dispatch site governs
  # every replicate -- matching `boot::boot`'s per-call determinism.
  results <- future.apply::future_lapply(
    indices_list,
    function(idx) statistic(data, idx),
    future.seed = TRUE
  )

  # Assemble an (R x k) matrix matching `boot::boot()$t`. Failed
  # replicates (NA vectors from the caller's tryCatch) slot in as NA
  # rows; `process_boot_results()` already filters those via
  # `complete.cases()`.
  t_mat <- do.call(rbind, results)

  list(t = t_mat, R = R, sim = "ordinary", stype = "i")
}


#' Process bootstrap results into vcov matrix and boot_t
#'
#' Shared post-processing for `variance_bootstrap()` and
#' `ice_variance_bootstrap()`: extracts complete replicates, warns on
#' failures, and computes the variance-covariance matrix.
#'
#' @param boot_res A `boot::boot` result object.
#' @param int_names Character vector of intervention names.
#' @param n_boot Integer. Total number of requested replicates.
#' @return A list with `vcov` (k x k matrix), `boot_t` (matrix of
#'   successful replicates), and `boot_info` (a 3-element list of
#'   `n_requested`, `n_ok`, `n_fail`) so downstream consumers
#'   (`contrast()`, `print.causatr_result()`) can surface the bootstrap
#'   failure rate without re-deriving it from `boot_t`.
#' @noRd
process_boot_results <- function(boot_res, int_names, n_boot) {
  # `boot_res$t` is an (R x k) matrix of statistics. Failed replicates
  # are flagged by NA rows (each per-rep function wraps its body in
  # tryCatch and returns `rep(NA_real_, k)` on error), so
  # `complete.cases()` identifies the usable ones.
  #
  # We also track per-intervention failure counts in the returned
  # `boot_info$n_fail_by_int`. This is the only way to see whether
  # failures cluster on one intervention (e.g. a rare static value
  # that triggers separation in every resample) versus are spread
  # across the whole replicate (e.g. a bad factor-level sample): the
  # whole-row drop biases the vcov more in the former case. Downstream
  # `print.causatr_result()` surfaces the vector so users spot the
  # pattern.
  t_mat <- boot_res$t
  colnames(t_mat) <- int_names
  complete_rows <- stats::complete.cases(t_mat)
  n_ok <- sum(complete_rows)
  n_fail <- n_boot - n_ok
  # Per-intervention failure counts: for each column, count NA rows.
  # This is strictly `>= n_fail / k` (a failed replicate contributes
  # an NA to every column), but is a better diagnostic because a
  # heterogeneous failure pattern reveals which intervention caused
  # the failure.
  n_fail_by_int <- vapply(
    seq_along(int_names),
    function(j) sum(is.na(t_mat[, j])),
    integer(1)
  )
  names(n_fail_by_int) <- int_names

  if (n_ok < 2L) {
    # With fewer than 2 successful replicates the sample variance is
    # undefined. Return a NA vcov so downstream `sqrt(diag(.))` yields
    # NA SE rather than a misleading zero.
    rlang::warn(
      "Bootstrap produced fewer than 2 non-NA replicates; SE estimates will be NA."
    )
    vcov_mat <- matrix(
      NA_real_,
      nrow = length(int_names),
      ncol = length(int_names)
    )
  } else {
    if (n_fail > 0L) {
      pct <- round(100 * n_fail / n_boot, 1)
      # Two different warnings for the same problem -- the 20%
      # threshold is an ad-hoc "this is getting bad" line. High
      # failure rates usually mean a fragile pipeline (e.g. factor
      # levels absent from some resamples), and the resulting SEs
      # are biased because the successful replicates aren't a
      # random sample of the sampling distribution.
      if (pct > 20) {
        rlang::warn(paste0(
          n_fail,
          " of ",
          n_boot,
          " bootstrap replicates (",
          pct,
          "%) failed. High failure rate may indicate model ",
          "instability; variance estimates may be unreliable."
        ))
      } else {
        rlang::warn(paste0(
          n_fail,
          " of ",
          n_boot,
          " bootstrap replicates (",
          pct,
          "%) failed and were discarded. ",
          "Variance is estimated from the ",
          n_ok,
          " successful replicates."
        ))
      }
    }
    # The bootstrap vcov is just the sample covariance of the
    # successful replicates (Davison & Hinkley 1997, Sec.2.5.3).
    vcov_mat <- stats::var(t_mat[complete_rows, , drop = FALSE])
  }

  rownames(vcov_mat) <- int_names
  colnames(vcov_mat) <- int_names

  boot_t <- t_mat[complete_rows, , drop = FALSE]
  colnames(boot_t) <- int_names

  list(
    vcov = vcov_mat,
    boot_t = boot_t,
    boot_info = list(
      n_requested = n_boot,
      n_ok = n_ok,
      n_fail = n_fail,
      n_fail_by_int = n_fail_by_int
    )
  )
}

#' Confidence-interval bounds from a bootstrap replicate vector
#'
#' @description
#' Computes a lower/upper confidence bound from a vector of bootstrap
#' replicate statistics, by one of two standard methods (Efron & Tibshirani
#' 1993, ch. 13; Davison & Hinkley 1997, sec. 5.2):
#'
#' \itemize{
#'   \item `"percentile"`: the empirical \eqn{\alpha/2} and \eqn{1-\alpha/2}
#'     quantiles of the replicates. Transformation-respecting and bounded by
#'     the estimand's support, so it is appropriate for probabilities, ratios,
#'     and odds ratios.
#'   \item `"normal"`: \eqn{\hat\theta \pm z_{1-\alpha/2}\,\widehat{sd}}, where
#'     \eqn{\widehat{sd}} is the standard deviation of the replicates. Symmetric;
#'     equals the Wald interval built from the bootstrap standard error.
#' }
#'
#' Both are computed from the same replicate vector, so offering both costs no
#' extra resampling.
#'
#' @param reps Numeric vector of bootstrap replicate statistics (one per
#'   successful replicate). Non-finite entries are dropped.
#' @param estimate Numeric scalar. The point estimate, used to centre the
#'   `"normal"` interval (ignored for `"percentile"`).
#' @param level Numeric scalar in (0, 1). Confidence level.
#' @param method `"percentile"` or `"normal"`.
#' @returns A length-2 numeric vector named `c("lower", "upper")`. Returns
#'   `c(NA, NA)` when fewer than two finite replicates are available.
#' @noRd
boot_ci_block <- function(reps, estimate, level, method) {
  alpha <- (1 - level) / 2
  reps <- reps[is.finite(reps)]
  if (length(reps) < 2L) {
    return(c(lower = NA_real_, upper = NA_real_))
  }
  bounds <- if (method == "percentile") {
    stats::quantile(reps, c(alpha, 1 - alpha), names = FALSE)
  } else {
    # Normal approximation: centre on the point estimate (not the replicate
    # mean) so the interval is the Wald CI from the bootstrap SE -- this keeps
    # `boot_ci = "normal"` identical to the legacy vcov-based bounds.
    z <- stats::qnorm(1 - alpha)
    estimate + c(-1, 1) * z * stats::sd(reps)
  }
  c(lower = bounds[1], upper = bounds[2])
}

#' Bootstrap variance-covariance matrix for marginal means
#'
#' @description
#' Estimates \eqn{Var(\hat{\mu})} by resampling the entire estimation pipeline `n_boot`
#' times (Hernan & Robins, Technical Point 13.1).  Works for all causal
#' estimators (g-comp, IPW, matching).
#'
#' ## Algorithm
#'
#' For each bootstrap replicate b = 1, ..., B:
#' 1. Draw n individuals with replacement from the full dataset.
#' 2. Refit the full estimation pipeline on the bootstrap sample:
#'    - **g-comp**: refit the outcome model on uncensored rows.
#'    - **IPW**: re-estimate weights + refit the weighted MSM.
#'    - **matching**: re-match + refit on the matched sample.
#' 3. For each intervention, apply it and compute the marginal mean.
#' 4. Collect the k-vector of marginal means.
#'
#' The \eqn{k \times k} bootstrap vcov is `var(boot_replicates)`.
#'
#' @param fit A `causatr_fit` object.
#' @param interventions Named list of `causatr_intervention` objects.
#' @param n_boot Positive integer. Number of bootstrap replicates.
#' @param target_idx Logical vector (length n) flagging target-population rows.
#' @param est Character. Estimand string (`"ATE"`, `"ATT"`, or `"ATC"`).
#' @param subset Quoted expression or `NULL`.
#'
#' @return A named \eqn{k \times k} variance-covariance matrix.
#'
#' @noRd
variance_bootstrap <- function(
  fit,
  interventions,
  n_boot,
  target_idx,
  est,
  subset,
  parallel = "no",
  ncpus = 1L,
  subset_env = parent.frame(),
  trim = 1
) {
  data <- fit$data
  int_names <- names(interventions)
  estimator <- fit$estimator
  outcome <- fit$outcome
  treatment <- fit$treatment
  censoring <- fit$censoring

  # boot_fn: called by boot::boot() for each replicate.
  # The bootstrap is correct here -- not just a computational convenience --
  # because it resamples individuals and refits the ENTIRE pipeline
  # (treatment model + outcome model + standardisation) on each resample.
  # This reproduces the full sampling distribution including the
  # uncertainty from nuisance estimation, which a naive vcov of predictions
  # under a single fit would miss. The IF sandwich is an analytic
  # approximation to the same quantity; bootstrap is the non-parametric
  # alternative that avoids assuming parametric models for the IF terms.
  #
  # The entire body is wrapped in tryCatch because bootstrap samples may
  # cause downstream failures (e.g. factor levels absent from uncensored
  # rows but present in prediction data). Failed replicates return NA and
  # are excluded from the variance calculation -- this is the standard
  # bootstrap approach (Davison & Hinkley, 1997, Sec.2.5.3).
  boot_fn <- function(d, indices) {
    tryCatch(
      {
        d_b <- d[indices]
        # Resample weights alongside data rows. When IPCW is active,
        # use the pre-IPCW weights — each refit_*() function refits
        # the censoring model and recomposes IPCW weights internally.
        # Resampling the pre-IPCW weights rather than the final composed
        # weights ensures the censoring model is re-estimated on the
        # bootstrap sample's censoring pattern, not the original's.
        orig_w <- if (isTRUE(fit$details$ipcw)) {
          fit$details$weights_pre_ipcw
        } else {
          fit$details$weights
        }
        w_b <- if (!is.null(orig_w)) orig_w[indices] else NULL

        # Returns the transport/generalizability override for target_idx, or
        # NULL when no transport is active (NULL-safe for %||% chaining).
        transport_override <- function(data) {
          if (is.null(fit$target)) {
            return(NULL)
          }
          if (identical(fit$target_subset, "target")) {
            return(data[[fit$target]] == 0L)
          }
          rep(TRUE, nrow(data))
        }

        # IPW replicates refit the density model on `d_b` via
        # `refit_ipw()` and then reuse the analytic
        # `compute_ipw_contrast_point()` path: one weighted MSM per
        # intervention, reading off the marginal mean from each. This
        # keeps the bootstrap and sandwich paths structurally
        # identical -- both funnel through the same per-intervention
        # MSM builder -- so a plumbing regression in one surfaces in
        # the other.
        if (estimator == "ipw") {
          # Silence the GLM "non-integer #successes in a binomial
          # glm" warning (fires whenever external continuous-valued
          # weights enter a weighted logistic propensity fit) and
          # the `causatr_singular_bread` warning on thin resamples.
          # Both are expected noise under bootstrap and would
          # otherwise flood the console or be escalated to errors
          # under `testthat` edition 3's warning-as-error rule.
          repl <- withCallingHandlers(
            {
              fit_b <- refit_ipw(fit, d_b, weights = w_b)
              target_idx_b <- get_target_idx(
                d_b,
                treatment,
                est,
                subset,
                subset_env = subset_env
              )
              target_idx_b <- transport_override(d_b) %||% target_idx_b
              ipw_point_b <- compute_ipw_contrast_point(
                fit_b,
                interventions,
                target_idx_b,
                trim = trim
              )
              ipw_point_b$mu_hat
            },
            warning = function(w) {
              if (inherits(w, "causatr_singular_bread")) {
                invokeRestart("muffleWarning")
              }
              msg <- conditionMessage(w)
              if (
                grepl(
                  "non-integer #successes in a binomial glm",
                  msg,
                  fixed = TRUE
                )
              ) {
                invokeRestart("muffleWarning")
              }
            }
          )
          return(repl)
        }

        if (estimator == "aipw") {
          repl <- withCallingHandlers(
            {
              fit_b <- refit_aipw(fit, d_b, weights = w_b)
              target_idx_b <- get_target_idx(
                d_b,
                treatment,
                est,
                subset,
                subset_env = subset_env
              )
              target_idx_b <- transport_override(d_b) %||% target_idx_b
              aipw_point_b <- compute_aipw_contrast_point(
                fit_b,
                interventions,
                target_idx_b,
                trim = trim
              )
              aipw_point_b$mu_hat
            },
            warning = function(w) {
              if (inherits(w, "causatr_singular_bread")) {
                invokeRestart("muffleWarning")
              }
              msg <- conditionMessage(w)
              if (
                grepl(
                  "non-integer #successes in a binomial glm",
                  msg,
                  fixed = TRUE
                )
              ) {
                invokeRestart("muffleWarning")
              }
            }
          )
          return(repl)
        }

        # Demote only the specific warnings we expect to emit on
        # nearly every replicate and that would otherwise flood the
        # console: GLM "fitted probabilities numerically 0 or 1",
        # `X'WX` near-singular on thin resamples, and MatchIt's
        # "Fewer control/treated units". Everything else -- non-
        # convergence, family mismatches, factor-level surprises -- is
        # allowed to surface so users can spot pipeline instability
        # instead of seeing a silent NA column in `boot_t`.
        model_b <- withCallingHandlers(
          refit_model(fit, d_b, weights = w_b),
          warning = function(w) {
            # T7 (2026-04-15 third-round review): match the singular
            # bread warning by its `causatr_singular_bread` class,
            # not by substring -- refactor-safe.
            if (inherits(w, "causatr_singular_bread")) {
              invokeRestart("muffleWarning")
            }
            msg <- conditionMessage(w)
            if (
              grepl(
                "fitted probabilities numerically 0 or 1",
                msg,
                fixed = TRUE
              ) ||
                grepl("Fewer (control|treated) units", msg)
            ) {
              invokeRestart("muffleWarning")
            }
          }
        )

        target_idx_b <- get_target_idx(
          d_b,
          treatment,
          est,
          subset,
          subset_env = subset_env
        )
        target_idx_b <- transport_override(d_b) %||% target_idx_b

        # MTP + transport: fit a treatment model on the bootstrap
        # sample's study rows and MC-marginalize predictions for
        # target rows where A is unobserved.
        is_transport_b <- !is.null(fit$target)
        mc_tm_b <- NULL
        mc_rows_b <- NULL
        any_mc_b <- FALSE
        if (is_transport_b && length(treatment) == 1L) {
          mc_rows_b <- is.na(d_b[[treatment]])
          if (
            any(mc_rows_b) &&
              any(vapply(
                interventions,
                needs_observed_treatment,
                logical(1)
              ))
          ) {
            any_mc_b <- TRUE
            fit_rows_b <- get_fit_rows(
              d_b,
              fit$outcome,
              fit$censoring,
              target = fit$target
            )
            mc_tm_b <- fit_mc_treatment_model(
              d_b,
              treatment,
              resolve_confounders_treatment(fit),
              fit_rows_b
            )
          }
        }

        vapply(
          interventions,
          function(iv) {
            data_a_b <- apply_intervention(d_b, treatment, iv)
            pred_a_b <- predict(
              model_b,
              newdata = data_a_b,
              type = "response"
            )
            # MC-marginalize target-row predictions for MTP interventions
            if (any_mc_b && needs_observed_treatment(iv)) {
              pred_a_b[mc_rows_b] <- mc_marginalize_preds(
                model_b,
                d_b[mc_rows_b],
                treatment,
                iv,
                mc_tm_b
              )
            }
            valid <- target_idx_b & !is.na(pred_a_b)
            maybe_weighted_mean(
              pred_a_b[valid],
              if (!is.null(w_b)) w_b[valid]
            )
          },
          numeric(1)
        )
      },
      error = function(e) rep(NA_real_, length(int_names))
    )
  }

  # Point-treatment bootstrap resamples rows (individuals) with replacement.
  # For longitudinal data, the equivalent functions in
  # variance_bootstrap_longitudinal.R resample ENTIRE individual trajectories
  # (all person-period rows for a given id move together). This "cluster
  # bootstrap" preserves within-individual correlation across time -- sampling
  # rows independently would break the temporal dependence structure and
  # underestimate uncertainty by treating repeated observations as independent.
  boot_res <- dispatch_boot(
    data = data,
    statistic = boot_fn,
    R = n_boot,
    parallel = parallel,
    ncpus = ncpus
  )

  process_boot_results(boot_res, int_names, n_boot)
}

#' Bootstrap variance for a multinomial-outcome g-computation
#'
#' @description
#' The multinomial analogue of `variance_bootstrap()`. Each replicate refits
#' the `nnet::multinom` outcome model and, for every intervention, averages
#' the predicted class probabilities over the target rows to obtain a
#' K-vector \eqn{P(Y = k \mid do(A = a))}. The replicate statistic is the
#' full \eqn{k_{int} \times K} surface flattened in class-major order (class
#' outer, intervention inner), so the shared `process_boot_results()` engine
#' computes one big covariance from which the per-class \eqn{k_{int} \times
#' k_{int}} blocks are sliced. Per-class contrasts only need the within-class
#' block, so the cross-class covariance is not retained.
#'
#' @param fit A `causatr_fit` with a multinomial outcome model.
#' @param interventions Named list of `causatr_intervention` objects.
#' @param n_boot Positive integer. Number of replicates.
#' @param target_idx Logical vector flagging target-population rows.
#' @param est Estimand string.
#' @param subset Quoted subset expression or `NULL`.
#' @param parallel,ncpus Parallelism controls forwarded to `dispatch_boot()`.
#' @param class_labels Character vector of the K outcome class labels.
#' @param int_names Character vector of intervention names.
#' @param subset_env Environment for resolving `subset`.
#' @return A list with `vcov` (named list of K per-class \eqn{k_{int} \times
#'   k_{int}} matrices), `boot_t` (named list of K replicate matrices, one per
#'   class), and `boot_info`.
#' @noRd
variance_bootstrap_multinom <- function(
  fit,
  interventions,
  n_boot,
  target_idx,
  est,
  subset,
  parallel = "no",
  ncpus = 1L,
  class_labels,
  int_names,
  subset_env = parent.frame()
) {
  data <- fit$data
  treatment <- fit$treatment
  k_int <- length(int_names)
  k_class <- length(class_labels)

  boot_fn <- function(d, indices) {
    tryCatch(
      {
        d_b <- d[indices]
        orig_w <- fit$details$weights
        w_b <- if (!is.null(orig_w)) orig_w[indices] else NULL

        # Muffle the same near-singular / separation warnings the scalar
        # bootstrap demotes; everything else surfaces so users can spot
        # pipeline instability instead of a silent NA column.
        model_b <- withCallingHandlers(
          refit_model(fit, d_b, weights = w_b),
          warning = function(w) {
            if (inherits(w, "causatr_singular_bread")) {
              invokeRestart("muffleWarning")
            }
            msg <- conditionMessage(w)
            if (
              grepl(
                "fitted probabilities numerically 0 or 1",
                msg,
                fixed = TRUE
              )
            ) {
              invokeRestart("muffleWarning")
            }
          }
        )

        target_idx_b <- get_target_idx(
          d_b,
          treatment,
          est,
          subset,
          subset_env = subset_env
        )

        # Per intervention: average each class column over valid target rows.
        per_int <- lapply(interventions, function(iv) {
          data_a_b <- apply_intervention(d_b, treatment, iv)
          probs_b <- predict_outcome(model_b, data_a_b)
          valid <- target_idx_b & !apply(is.na(probs_b), 1L, any)
          apply(probs_b[valid, , drop = FALSE], 2L, function(col) {
            maybe_weighted_mean(col, if (!is.null(w_b)) w_b[valid])
          })
        })

        # Flatten class-major (class outer, intervention inner) so the
        # class-c block occupies positions (c-1)*k_int + (1:k_int).
        unlist(lapply(seq_len(k_class), function(ci) {
          vapply(per_int, `[`, numeric(1), ci)
        }))
      },
      error = function(e) rep(NA_real_, k_int * k_class)
    )
  }

  boot_res <- dispatch_boot(
    data = data,
    statistic = boot_fn,
    R = n_boot,
    parallel = parallel,
    ncpus = ncpus
  )

  # Class-major flat labels matching boot_fn's ordering.
  flat_names <- paste(
    rep(int_names, times = k_class),
    rep(class_labels, each = k_int),
    sep = ":"
  )
  processed <- process_boot_results(boot_res, flat_names, n_boot)

  # Slice the big covariance into per-class k_int-by-k_int blocks keyed by
  # class label (per-class contrasts only need the within-class covariance),
  # restoring intervention names on the rows/cols. The flat replicate matrix
  # `boot_t` is kept whole (class-major columns) so `confint()` can read
  # percentile intervals that line up row-for-row with the class-major
  # `estimates` table.
  vcov_list <- vector("list", k_class)
  names(vcov_list) <- class_labels
  for (ci in seq_len(k_class)) {
    idx_c <- (ci - 1L) * k_int + seq_len(k_int)
    block <- processed$vcov[idx_c, idx_c, drop = FALSE]
    dimnames(block) <- list(int_names, int_names)
    vcov_list[[ci]] <- block
  }

  list(
    vcov = vcov_list,
    boot_t = processed$boot_t,
    boot_info = processed$boot_info
  )
}

#' Refit the full estimation pipeline on a bootstrap sample
#'
#' @description
#' Dispatches to the appropriate refitting logic based on the causal
#' estimator stored in `fit$estimator`.
#'
#' @param fit A `causatr_fit` object (original fit, for extracting formulas etc.).
#' @param d_b A data.table bootstrap sample.
#'
#' @return A fitted model object (glm, glm_weightit, etc.) or NULL on failure.
#'
#' @noRd
refit_model <- function(fit, d_b, weights = NULL) {
  if (fit$estimator == "gcomp") {
    refit_gcomp(fit, d_b, weights = weights)
  } else if (fit$estimator == "ipw") {
    refit_ipw(fit, d_b, weights = weights)
  } else if (fit$estimator == "aipw") {
    refit_aipw(fit, d_b, weights = weights)
  } else if (fit$estimator == "matching") {
    refit_matching(fit, d_b, weights = weights)
  } else {
    rlang::abort(
      paste0(
        "Bootstrap is not supported for estimator = '",
        fit$estimator,
        "'."
      ),
      class = "causatr_unknown_estimator",
      .call = FALSE
    )
  }
}

#' Refit g-comp outcome model on a bootstrap sample
#'
#' @param fit A `causatr_fit` object.
#' @param d_b A data.table bootstrap sample.
#' @return A fitted model object.
#' @noRd
refit_gcomp <- function(fit, d_b, weights = NULL) {
  model_formula <- stats::formula(fit$model)
  family <- fit$model$family
  model_fn <- fit$details$model_fn
  censoring <- fit$censoring
  outcome <- fit$outcome

  # IPCW: refit censoring model and compose weights
  if (isTRUE(fit$details$ipcw)) {
    ipcw_w_b <- refit_censoring_weights(fit, d_b)
    weights <- if (is.null(weights)) ipcw_w_b else weights * ipcw_w_b
  }

  fit_rows_b <- get_fit_rows(d_b, outcome, censoring, target = fit$target)

  args <- list(model_formula, data = d_b[fit_rows_b])
  if (fn_accepts_family(model_fn)) {
    args$family <- family
  }
  if (!is.null(weights)) {
    args$weights <- weights[fit_rows_b]
  }
  # Replay the user's original `...` via the central `replay_fit()`
  # helper so duplicate-key / positional-dot handling is identical to
  # every other refit site. See the 2026-04-15 third-round dots audit.
  replay_fit(model_fn, args, fit$details$dots)
}

#' Refit the IPW pipeline on a bootstrap sample
#'
#' @description
#' Replays `fit_ipw()` on the resampled data, reusing every
#' ingredient of the original fit: the user-supplied
#' `propensity_model_fn`, outcome `model_fn`, `confounders`,
#' `estimand`, stashed `dots`, and -- crucially -- the resampled slice
#' of the external weight vector. Returns a fresh `causatr_fit`
#' that `boot_fn()` feeds straight into
#' `compute_ipw_contrast_point()` to read off the replicate's
#' marginal-mean vector. Bootstrap and sandwich thus share the same
#' per-intervention MSM builder, so plumbing regressions surface
#' symmetrically.
#'
#' @param fit The original `causatr_fit` from `fit_ipw()`.
#' @param d_b Resampled `data.table`.
#' @param weights Resampled external weight vector, or `NULL`.
#'
#' @return A fresh `causatr_fit` object for the bootstrap replicate.
#'
#' @noRd
refit_ipw <- function(fit, d_b, weights = NULL) {
  # IPCW: refit censoring model and compose weights before the
  # full IPW pipeline replay
  if (isTRUE(fit$details$ipcw)) {
    ipcw_w_b <- refit_censoring_weights(fit, d_b)
    weights <- if (is.null(weights)) ipcw_w_b else weights * ipcw_w_b
  }

  # `do.call()` with the stashed `dots` replays the exact propensity
  # fitter configuration (smoother arguments, family overrides, etc.)
  # the user passed to the original `causat()` call. Every non-dots
  # argument is rebuilt from the original fit's slots.
  #
  # Note on `call`: `do.call()` evaluates its argument list in the
  # calling frame before dispatching, so putting `call = fit$call`
  # inside `args` would try to re-evaluate the user's original
  # `causat(data = d_local, ...)` expression in whatever frame
  # `boot::boot` is running -- which is not the user's session -- and
  # explode on a missing-object lookup. We pass `call = NULL` into
  # `fit_ipw()` and patch the original back in after the refit.
  args <- list(
    data = d_b,
    outcome = fit$outcome,
    treatment = fit$treatment,
    confounders = resolve_confounders_treatment(fit),
    confounders_tv = resolve_confounders_tv_treatment(fit),
    confounders_outcome = fit$details$confounders_outcome,
    confounders_outcome_raw = fit$details$confounders_outcome,
    confounders_treatment_raw = fit$details$confounders_treatment,
    confounders_tv_outcome_raw = fit$details$confounders_tv_outcome,
    confounders_tv_treatment_raw = fit$details$confounders_tv_treatment,
    family = fit$family,
    estimand = fit$estimand,
    type = fit$type,
    history = fit$history,
    numerator = fit$numerator,
    weights = weights,
    model_fn = fit$details$model_fn,
    propensity_model_fn = fit$details$propensity_model_fn,
    propensity_family = fit$details$propensity_family,
    stabilize = fit$details$stabilize %||% "none",
    call = NULL,
    target = fit$target
  )
  # Longitudinal IPW needs the id / time slots threaded through to
  # `fit_longitudinal_ipw()`; point IPW ignores them.
  if (fit$type == "longitudinal") {
    args$id <- fit$id
    args$time <- fit$time
  }
  fit_b <- do.call(fit_ipw, c(args, fit$details$dots))
  fit_b$call <- fit$call

  # Transport: refit sampling model on bootstrap replicate. For
  # longitudinal fits, the sampling model must use first-period rows only
  # (same as the original fit) to avoid K-fold inflation.
  if (!is.null(fit$target)) {
    if (fit$type == "longitudinal") {
      first_t_b <- fit$details$time_points[1]
      rows_first_b <- d_b[[fit$time]] == first_t_b
      d_samp_b <- d_b[rows_first_b]
      w_samp_b <- if (is.null(weights)) NULL else weights[rows_first_b]
    } else {
      d_samp_b <- d_b
      w_samp_b <- weights
    }
    samp_model_b <- fit_sampling_model(
      d_samp_b,
      fit$target,
      resolve_confounders_sampling(fit),
      fit$treatment,
      model_fn = fit$details$sampling_model_fn,
      weights = w_samp_b
    )
    fit_b$details$transport <- TRUE
    fit_b$details$sampling_model <- samp_model_b
    fit_b$details$sampling_model_fn <- fit$details$sampling_model_fn
    fit_b$details$target_subset <- fit$target_subset
    fit_b$target <- fit$target
    fit_b$target_subset <- fit$target_subset
  }

  fit_b
}

#' Refit both AIPW nuisance models on a bootstrap sample
#'
#' @description
#' Replays `fit_aipw()` on the resampled data, refitting both the
#' outcome model \eqn{E[Y \mid A, L]} and the treatment density
#' \eqn{f(A \mid L)}. Returns a fresh `causatr_fit` that
#' `compute_aipw_contrast_point()` consumes.
#'
#' @param fit The original `causatr_fit` from `fit_aipw()`.
#' @param d_b Resampled `data.table`.
#' @param weights Resampled external weight vector, or `NULL`.
#'
#' @return A fresh `causatr_fit` object for the bootstrap replicate.
#'
#' @noRd
refit_aipw <- function(fit, d_b, weights = NULL) {
  if (isTRUE(fit$details$ipcw)) {
    ipcw_w_b <- refit_censoring_weights(fit, d_b)
    weights <- if (is.null(weights)) ipcw_w_b else weights * ipcw_w_b
  }

  if (fit$type == "longitudinal") {
    args <- list(
      data = d_b,
      outcome = fit$outcome,
      treatment = fit$treatment,
      confounders = resolve_confounders_outcome(fit),
      confounders_tv = resolve_confounders_tv_outcome(fit),
      family = fit$family,
      estimand = fit$estimand,
      history = fit$history,
      censoring = fit$censoring,
      weights = weights,
      model_fn = fit$details$model_fn,
      propensity_model_fn = fit$details$propensity_model_fn,
      propensity_family = fit$details$propensity_family,
      id = fit$id,
      time = fit$time,
      call = NULL,
      confounders_treatment = resolve_confounders_treatment(fit),
      confounders_tv_treatment = resolve_confounders_tv_treatment(fit)
    )
    fit_b <- do.call(
      fit_aipw_longitudinal,
      c(args, fit$details$dots)
    )
    fit_b$call <- fit$call
    return(fit_b)
  }

  args <- list(
    data = d_b,
    outcome = fit$outcome,
    treatment = fit$treatment,
    confounders = resolve_confounders_outcome(fit),
    family = fit$family,
    estimand = fit$estimand,
    censoring = fit$censoring,
    weights = weights,
    model_fn = fit$details$model_fn,
    propensity_model_fn = fit$details$propensity_model_fn,
    propensity_family = fit$details$propensity_family,
    stabilize = fit$details$stabilize %||% "none",
    call = NULL,
    target = fit$target,
    confounders_treatment = resolve_confounders_treatment(fit)
  )
  fit_b <- do.call(fit_aipw_point, c(args, fit$details$dots))
  fit_b$call <- fit$call

  # Transport: refit sampling model on bootstrap replicate
  if (!is.null(fit$target)) {
    samp_model_b <- fit_sampling_model(
      d_b,
      fit$target,
      resolve_confounders_sampling(fit),
      fit$treatment,
      model_fn = fit$details$sampling_model_fn,
      weights = weights
    )
    fit_b$details$transport <- TRUE
    fit_b$details$sampling_model <- samp_model_b
    fit_b$details$sampling_model_fn <- fit$details$sampling_model_fn
    fit_b$details$target_subset <- fit$target_subset
    fit_b$target <- fit$target
    fit_b$target_subset <- fit$target_subset
  }

  fit_b
}

#' Re-match and refit outcome model on a bootstrap sample
#'
#' Re-matching on each bootstrap sample is essential: the matched pairs are
#' determined by the propensity scores estimated on the original data. If we
#' kept the original matched set and only refit the outcome model, the bootstrap
#' distribution would not capture uncertainty from the matching step itself,
#' producing anti-conservative SEs. Full re-matching is computationally heavier
#' but the only valid approach.
#'
#' @param fit A `causatr_fit` object.
#' @param d_b A data.table bootstrap sample.
#' @return A `glm` model fit on the matched bootstrap data.
#' @noRd
refit_matching <- function(fit, d_b, weights = NULL) {
  # IPCW: refit censoring model and compose weights before matching
  if (isTRUE(fit$details$ipcw)) {
    ipcw_w_b <- refit_censoring_weights(fit, d_b)
    weights <- if (is.null(weights)) {
      ipcw_w_b
    } else {
      weights * ipcw_w_b
    }
  }

  ps_formula <- build_ps_formula(
    resolve_confounders_treatment(fit),
    fit$treatment
  )

  fit_rows_b <- get_fit_rows(d_b, fit$outcome)
  fit_data_b <- as.data.frame(d_b[fit_rows_b])

  # Replay the original MatchIt arguments stashed by `fit_matching()`
  # (caliper, ratio, distance, method overrides, ...). `fit_matching()`
  # captures `dots` AFTER the ATE->"full" defaulting, so the stashed
  # dots already carry `method =` when required. See B2 / audit E.
  base_args <- list(ps_formula, data = fit_data_b, estimand = fit$estimand)
  if (fit$estimand == "ATE") {
    check_pkg("optmatch")
  }
  m_b <- replay_fit(MatchIt::matchit, base_args, fit$details$dots)
  matched_b <- MatchIt::match.data(m_b)

  # Combine match weights with external weights via the shared helper
  # (defined in R/matching.R). Using the same code path as `fit_matching()`
  # means any row-name invariant violation fails loudly at the bootstrap
  # boundary rather than silently producing NA-tainted or misaligned
  # weights on a subset of replicates.
  matched_weights <- combine_match_and_external_weights(
    matched_b,
    external_weights = weights,
    fit_rows = fit_rows_b
  )

  # Replay the same MSM formula used at fit time: `Y ~ A` without EM,
  # `Y ~ A + modifier + A:modifier` with EM. `em_info` was stashed in
  # `fit$details` by `fit_matching()`.
  em_info <- fit$details$em_info
  msm_formula <- if (!is.null(em_info)) {
    build_matching_msm_formula(fit$outcome, fit$treatment, em_info)
  } else {
    stats::reformulate(fit$treatment, response = fit$outcome)
  }
  # `glm()` evaluates `weights` in the formula's environment. Reset
  # to the current frame so `matched_weights` is found.
  environment(msm_formula) <- environment()
  stats::glm(
    msm_formula,
    data = matched_b,
    weights = matched_weights,
    family = fit$model$family
  )
}
