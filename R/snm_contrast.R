#' Compute SNM contrast: blip parameters or averaged blip effect
#'
#' @description
#' Two paths depending on whether `treatment_values` is supplied:
#'
#' **Path A** (`treatment_values = NULL`): Returns the blip parameter
#' table directly — each \eqn{\psi_j} is its own estimand with SE and CI.
#'
#' **Path B** (`treatment_values = c(a0, a1)`): Computes the
#' population-averaged blip effect
#' \deqn{\Delta = \frac{1}{n}\sum_i \bigl[\gamma(a_1, L_i; \hat\psi)
#'   - \gamma(a_0, L_i; \hat\psi)\bigr]
#'   = (a_1 - a_0)\bigl(\hat\psi_0 + \sum_j \hat\psi_j \bar{m}_j\bigr)}
#' with delta-method variance from the blip vcov.
#'
#' @references
#' Robins JM (1994). Correcting for non-compliance in randomized trials
#' using structural nested mean models. *Comm Stat Theory Methods*,
#' 23(8), 2379–2412.
#'
#' @param fit A `causatr_fit` with `estimator = "snm"`.
#' @param treatment_values Numeric vector of length 2 or `NULL`.
#' @param ci_method `"sandwich"` or `"bootstrap"`.
#' @param conf_level Numeric confidence level.
#' @param n_boot Positive integer. Number of bootstrap replicates.
#' @param parallel Character. `"no"`, `"multicore"`, `"snow"`, or
#'   `"future"`.
#' @param ncpus Integer. Number of CPU cores for parallel backends.
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
  call
) {
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
      vcov_psi <- variance_if_snm_longitudinal(fit, snm_result)
    } else {
      vcov_psi <- variance_if_snm(fit, snm_result)
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

    # For linear blip: gamma(a, l; psi) = a * (psi_0 + sum psi_j m_j)
    # Average effect = (a1 - a0) * (psi_0 + sum psi_j * mean(m_j))
    M <- snm_result$M
    m_bar <- colMeans(M)
    avg_effect <- delta_a * sum(psi_hat * m_bar)

    if (ci_method == "bootstrap") {
      # Bootstrap vcov is already for the averaged blip effect (scalar)
      se_effect <- sqrt(max(vcov_psi[1, 1], 0))
    } else {
      # Delta method: grad_j = delta_a * mean(m_j)
      grad <- delta_a * m_bar
      var_effect <- as.numeric(t(grad) %*% vcov_psi %*% grad)
      se_effect <- sqrt(max(var_effect, 0))
    }

    comparison_label <- paste0("a=", a1, " vs a=", a0)
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
