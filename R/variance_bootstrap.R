#' Bootstrap variance–covariance matrix for marginal means
#'
#' @description
#' Estimates Var(μ̂) by resampling the entire estimation pipeline `n_boot`
#' times (Hernán & Robins, Technical Point 13.1).  Works for all estimation
#' methods (g-comp, IPW, matching).
#'
#' ## Algorithm
#'
#' For each bootstrap replicate b = 1, …, B:
#' 1. Draw n individuals with replacement from the full dataset.
#' 2. Refit the full estimation pipeline on the bootstrap sample:
#'    - **g-comp**: refit the outcome model on uncensored rows.
#'    - **IPW**: re-estimate weights + refit the weighted MSM.
#'    - **matching**: re-match + refit on the matched sample.
#' 3. For each intervention, apply it and compute the marginal mean.
#' 4. Collect the k-vector of marginal means.
#'
#' The k × k bootstrap vcov is `var(boot_replicates)`.
#'
#' @param fit A `causatr_fit` object.
#' @param interventions Named list of `causatr_intervention` objects.
#' @param n_boot Positive integer. Number of bootstrap replicates.
#' @param target_idx Logical vector (length n) flagging target-population rows.
#' @param est Character. Estimand string (`"ATE"`, `"ATT"`, or `"ATC"`).
#' @param subset Quoted expression or `NULL`.
#'
#' @return A named k × k variance–covariance matrix.
#'
#' @noRd
variance_bootstrap <- function(
  fit,
  interventions,
  n_boot,
  target_idx,
  est,
  subset
) {
  data <- fit$data
  int_names <- names(interventions)
  method <- fit$method
  outcome <- fit$outcome
  treatment <- fit$treatment
  censoring <- fit$censoring

  # boot_fn: called by boot::boot() for each replicate.
  boot_fn <- function(d, indices) {
    d_b <- d[indices]

    # Refit the full pipeline on the bootstrap sample, dispatching by method.
    # Suppress warnings from WeightIt/MatchIt about balance/weights in
    # edge resamples (these are expected and do not affect the bootstrap).
    model_b <- tryCatch(
      suppressWarnings(refit_model(fit, d_b)),
      error = function(e) NULL
    )

    if (is.null(model_b)) {
      return(rep(NA_real_, length(int_names)))
    }

    # Recompute target-population rows on the bootstrap sample.
    target_idx_b <- get_target_idx(d_b, treatment, est, subset)

    # For each intervention, apply it and compute the marginal mean.
    vapply(
      interventions,
      function(iv) {
        data_a_b <- apply_intervention(d_b, treatment, iv)
        pred_a_b <- predict(model_b, newdata = data_a_b, type = "response")
        valid <- target_idx_b & !is.na(pred_a_b)
        mean(pred_a_b[valid])
      },
      numeric(1)
    )
  }

  boot_res <- boot::boot(data = data, statistic = boot_fn, R = n_boot)

  t_mat <- boot_res$t
  complete_rows <- stats::complete.cases(t_mat)

  if (sum(complete_rows) < 2L) {
    rlang::warn(
      "Bootstrap produced fewer than 2 non-NA replicates; SE estimates will be NA."
    )
    vcov_mat <- matrix(
      NA_real_,
      nrow = length(int_names),
      ncol = length(int_names)
    )
  } else {
    vcov_mat <- stats::var(t_mat[complete_rows, , drop = FALSE])
  }

  rownames(vcov_mat) <- int_names
  colnames(vcov_mat) <- int_names
  vcov_mat
}

#' Refit the full estimation pipeline on a bootstrap sample
#'
#' @description
#' Dispatches to the appropriate refitting logic based on the estimation
#' method stored in `fit$method`.
#'
#' @param fit A `causatr_fit` object (original fit, for extracting formulas etc.).
#' @param d_b A data.table bootstrap sample.
#'
#' @return A fitted model object (glm, glm_weightit, etc.) or NULL on failure.
#'
#' @noRd
refit_model <- function(fit, d_b) {
  if (fit$method == "gcomp") {
    refit_gcomp(fit, d_b)
  } else if (fit$method == "ipw") {
    refit_ipw(fit, d_b)
  } else if (fit$method == "matching") {
    refit_matching(fit, d_b)
  } else {
    rlang::abort(
      paste0("Bootstrap not yet supported for method = '", fit$method, "'."),
      .call = FALSE
    )
  }
}

#' @noRd
refit_gcomp <- function(fit, d_b) {
  model_formula <- stats::formula(fit$model)
  family <- fit$model$family
  model_fn <- fit$details$model_fn
  censoring <- fit$censoring
  outcome <- fit$outcome

  fit_rows_b <- rep(TRUE, nrow(d_b))
  if (!is.null(censoring)) {
    cens_col <- d_b[[censoring]]
    fit_rows_b <- fit_rows_b & !is.na(cens_col) & (cens_col == 0)
  }
  fit_rows_b <- fit_rows_b & !is.na(d_b[[outcome]])

  model_fn(model_formula, data = d_b[fit_rows_b], family = family)
}

#' @noRd
refit_ipw <- function(fit, d_b) {
  # Re-estimate propensity-score weights.
  confounder_terms <- attr(stats::terms(fit$confounders), "term.labels")
  ps_formula <- stats::reformulate(confounder_terms, response = fit$treatment)

  fit_rows_b <- !is.na(d_b[[fit$outcome]])
  fit_data_b <- d_b[fit_rows_b]

  w_b <- WeightIt::weightit(
    ps_formula,
    data = fit_data_b,
    estimand = fit$estimand
  )

  # Refit the weighted MSM.
  msm_formula <- stats::reformulate(fit$treatment, response = fit$outcome)
  WeightIt::glm_weightit(
    msm_formula,
    data = fit_data_b,
    weightit = w_b,
    family = stats::gaussian()
  )
}

#' @noRd
refit_matching <- function(fit, d_b) {
  confounder_terms <- attr(stats::terms(fit$confounders), "term.labels")
  ps_formula <- stats::reformulate(confounder_terms, response = fit$treatment)

  fit_rows_b <- !is.na(d_b[[fit$outcome]])
  fit_data_b <- as.data.frame(d_b[fit_rows_b])

  m_b <- MatchIt::matchit(
    ps_formula,
    data = fit_data_b,
    estimand = fit$estimand
  )
  matched_b <- MatchIt::match.data(m_b)

  msm_formula <- stats::reformulate(fit$treatment, response = fit$outcome)
  stats::glm(
    msm_formula,
    data = matched_b,
    weights = matched_b$weights,
    family = stats::gaussian()
  )
}
