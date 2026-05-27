#' Compute SNM contrast: blip parameters or averaged blip effect
#'
#' @description
#' Three paths depending on whether `treatment_values` and `by` are supplied:
#'
#' **Path A** (`treatment_values = NULL`): Returns the blip parameter
#' table directly — each \eqn{\psi_j} is its own estimand with SE and CI.
#'
#' **Path B** (`treatment_values = c(a0, a1)`, `by = NULL`): Computes the
#' population-averaged blip effect
#' \deqn{\Delta = (a_1 - a_0)\bigl(\hat\psi_0 + \sum_j \hat\psi_j
#'   \bar{m}_j\bigr)}
#' with delta-method variance from the blip vcov.
#'
#' **Path C** (`treatment_values = c(a0, a1)`, `by = "col"`): Computes
#' per-stratum averaged blip effects. \eqn{\hat\psi} is estimated globally
#' on all fit rows; per-stratum means \eqn{\bar{m}_s = \mathrm{mean}(M_s)}
#' are averaged over the stratum-specific covariate distribution. This
#' implements the stratum-conditional averaged blip without re-estimating
#' \eqn{\hat\psi} per stratum — the delta-method SE uses the global vcov
#' with stratum-specific gradient \eqn{g_s = (a_1 - a_0)\,\bar{m}_s}.
#'
#' @references
#' Robins JM (1994). Correcting for non-compliance in randomized trials
#' using structural nested mean models. *Comm Stat Theory Methods*,
#' 23(8), 2379–2412.
#'
#' Vansteelandt S, Joffe M (2014). Structural nested models and G-estimation:
#' the partially realized promise. *Statistical Science*, 29(4), 707–731.
#'
#' @param fit A `causatr_fit` with `estimator = "snm"`.
#' @param treatment_values Numeric vector of length 2 or `NULL`.
#' @param ci_method `"sandwich"` or `"bootstrap"`.
#' @param conf_level Numeric confidence level.
#' @param n_boot Positive integer. Number of bootstrap replicates.
#' @param parallel Character. `"no"`, `"multicore"`, `"snow"`, or
#'   `"future"`.
#' @param ncpus Integer. Number of CPU cores for parallel backends.
#' @param cluster_vec Character, integer, or factor vector of cluster IDs
#'   (length = `nrow(fit$data)`) or `NULL` for i.i.d. aggregation. When
#'   non-`NULL`, the Liang-Zeger cluster-robust sandwich is used.
#' @param by Character column name for per-stratum averaged blip, or
#'   `NULL`. Requires `treatment_values`. Not supported for longitudinal
#'   or categorical-treatment SNM fits.
#' @param call The original `contrast()` call.
#' @return A `causatr_result` object.
#' @noRd
compute_snm_contrast <- function(
  fit,
  treatment_values,
  ci_method,
  conf_level,
  n_boot = 500L,
  parallel = "no",
  ncpus = 1L,
  cluster_vec = NULL,
  by = NULL,
  call
) {
  # `by` requires treatment_values: the blip parameter table (Path A) is
  # already a complete per-parameter breakdown; stratifying it is undefined.
  if (!is.null(by) && is.null(treatment_values)) {
    rlang::abort(
      c(
        "`by` requires `treatment_values` for SNM fits.",
        i = paste0(
          "Specify `treatment_values = c(a0, a1)` to average the blip ",
          "over the stratum-specific covariate distribution."
        )
      ),
      class = "causatr_snm_by_needs_treatment_values"
    )
  }

  # Validate by: must be a scalar string naming a column in fit$data.
  # Do this before the expensive blip computation.
  if (!is.null(by)) {
    check_string(by)
    if (!by %in% names(fit$data)) {
      rlang::abort(
        paste0("`by` column '", by, "' not found in fitted data."),
        class = "causatr_snm_by_not_found"
      )
    }
  }

  # Compute point estimates regardless of ci_method — sandwich and
  # bootstrap share the same psi_hat from the g-estimating equation.
  if (fit$type == "longitudinal") {
    snm_result <- compute_snm_blip_longitudinal(fit)
  } else {
    snm_result <- compute_snm_blip_point(fit)
  }
  psi_hat <- snm_result$psi_hat
  p_psi <- length(psi_hat)
  z <- stats::qnorm((1 + conf_level) / 2)
  n_target <- snm_result$n_obs

  # Variance: sandwich or bootstrap
  if (ci_method == "sandwich") {
    if (fit$type == "longitudinal") {
      # Longitudinal variance uses the ID structure for clustering internally.
      vcov_psi <- variance_if_snm_longitudinal(fit, snm_result)
    } else {
      vcov_psi <- variance_if_snm(fit, snm_result, cluster_vec = cluster_vec)
    }
  } else {
    # Bootstrap: resample and re-estimate blip parameters
    boot_res <- snm_variance_bootstrap(
      fit,
      treatment_values = treatment_values,
      n_boot = n_boot,
      parallel = parallel,
      ncpus = ncpus
    )
    vcov_psi <- boot_res$vcov
  }

  if (is.null(treatment_values)) {
    # Path A: blip parameter table — one row per parameter.
    # For longitudinal fits, psi_hat contains per-stage parameters
    # (e.g. stage0_psi_intercept, stage1_psi_intercept, ...).
    se_psi <- sqrt(pmax(diag(vcov_psi), 0))
    estimates_dt <- data.table::data.table(
      parameter = names(psi_hat),
      estimate = as.numeric(psi_hat),
      se = se_psi,
      ci_lower = psi_hat - z * se_psi,
      ci_upper = psi_hat + z * se_psi
    )
    contrasts_dt <- data.table::data.table(
      comparison = character(0),
      estimate = numeric(0),
      se = numeric(0),
      ci_lower = numeric(0),
      ci_upper = numeric(0)
    )
    vcov_out <- vcov_psi
  } else {
    # Path B: population-averaged blip effect.
    # For longitudinal fits, treatment_values is not meaningful (the
    # per-stage blip table IS the estimand). Reject with guidance.
    if (fit$type == "longitudinal") {
      rlang::abort(
        c(
          paste0(
            "`treatment_values` is not supported for longitudinal SNM fits."
          ),
          i = paste0(
            "Longitudinal SNMs return per-stage blip parameters directly. ",
            "Use `contrast(fit)` without `treatment_values`."
          )
        ),
        class = "causatr_snm_long_no_treatment_values"
      )
    }

    # Categorical treatments have per-level blip parameters — the
    # population-averaged blip effect via treatment_values only applies
    # to scalar (binary/continuous/count) treatments.
    if (snm_result$td$is_categorical) {
      rlang::abort(
        c(
          paste0(
            "`treatment_values` is not supported for categorical ",
            "treatment SNM fits."
          ),
          i = paste0(
            "Categorical SNMs estimate per-level blip parameters ",
            "(one set per non-reference level). Use `contrast(fit)` ",
            "without `treatment_values` to see them."
          )
        ),
        class = "causatr_snm_cat_no_treatment_values"
      )
    }

    if (length(treatment_values) != 2L) {
      rlang::abort(
        "`treatment_values` must be a numeric vector of length 2 (e.g. c(0, 1)).",
        class = "causatr_snm_bad_treatment_values"
      )
    }
    a0 <- treatment_values[1]
    a1 <- treatment_values[2]
    delta_a <- a1 - a0

    M <- snm_result$M

    # Helper to compute the averaged blip and its SE for one covariate mean.
    # g = (a1 - a0) * m_bar — the gradient of the averaged blip w.r.t. psi.
    # Delta-method variance: g' vcov g. Bootstrap: scalar vcov[1,1].
    compute_avg_blip <- function(m_bar) {
      avg <- delta_a * sum(psi_hat * m_bar)
      if (ci_method == "bootstrap") {
        se <- sqrt(max(vcov_psi[1, 1], 0))
      } else {
        grad <- delta_a * m_bar
        se <- sqrt(max(as.numeric(t(grad) %*% vcov_psi %*% grad), 0))
      }
      list(avg = avg, se = se)
    }

    comparison_label <- paste0("a=", a1, " vs a=", a0)

    if (is.null(by)) {
      # Path B: pooled averaged blip over all fit rows.
      # gamma_avg = (a1 - a0) * (psi_0 + sum psi_j * mean(m_j))
      m_bar <- colMeans(M)
      out <- compute_avg_blip(m_bar)
      avg_effect <- out$avg
      se_effect <- out$se

      estimates_dt <- data.table::data.table(
        parameter = "avg_blip_effect",
        estimate = avg_effect,
        se = se_effect,
        ci_lower = avg_effect - z * se_effect,
        ci_upper = avg_effect + z * se_effect
      )
      contrasts_dt <- data.table::data.table(
        comparison = comparison_label,
        estimate = avg_effect,
        se = se_effect,
        ci_lower = avg_effect - z * se_effect,
        ci_upper = avg_effect + z * se_effect
      )
      vcov_out <- vcov_psi
    } else {
      # Path C: by-stratified averaged blip.
      # psi_hat is estimated globally on all fit rows; per-stratum gradient
      # g_s = (a1 - a0) * colMeans(M[in stratum s, ]) captures the
      # stratum-specific covariate distribution. Delta-method SE uses the
      # global vcov — no re-estimation per stratum.
      fit_rows_lv <- snm_result$fit_rows
      by_vec <- fit$data[[by]][fit_rows_lv]
      by_levels <- sort(unique(stats::na.omit(by_vec)))

      est_rows <- lapply(by_levels, function(lv) {
        idx <- which(by_vec == lv)
        m_bar_s <- colMeans(M[idx, , drop = FALSE])
        out_s <- compute_avg_blip(m_bar_s)
        n_s <- length(idx)
        data.table::data.table(
          parameter = paste0("avg_blip_", by, "_", lv),
          estimate = out_s$avg,
          se = out_s$se,
          ci_lower = out_s$avg - z * out_s$se,
          ci_upper = out_s$avg + z * out_s$se,
          by = as.character(lv),
          n_by = n_s
        )
      })
      estimates_dt <- data.table::rbindlist(est_rows)

      con_rows <- lapply(by_levels, function(lv) {
        idx <- which(by_vec == lv)
        m_bar_s <- colMeans(M[idx, , drop = FALSE])
        out_s <- compute_avg_blip(m_bar_s)
        n_s <- length(idx)
        data.table::data.table(
          comparison = comparison_label,
          estimate = out_s$avg,
          se = out_s$se,
          ci_lower = out_s$avg - z * out_s$se,
          ci_upper = out_s$avg + z * out_s$se,
          by = as.character(lv),
          n_by = n_s
        )
      })
      contrasts_dt <- data.table::rbindlist(con_rows)

      # vcov_out is the full psi vcov (shared across strata)
      vcov_out <- vcov_psi
    }
  }

  new_causatr_result(
    estimates = estimates_dt,
    contrasts = contrasts_dt,
    type = "difference",
    estimand = "ATE",
    ci_method = ci_method,
    reference = NULL,
    interventions = NULL,
    n = n_target,
    estimator = "snm",
    family = fit$family,
    fit_type = fit$type,
    vcov = vcov_out,
    call = call
  )
}
