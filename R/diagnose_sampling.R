#' Sampling-model diagnostics for transportability / generalizability
#'
#' @description
#' Computes a summary table for the sampling model \eqn{P(S = 1 \mid L)}
#' and the derived sampling weights. Returns `NULL` when the fit does
#' not use transportability (i.e. `target = NULL`).
#'
#' The returned table includes:
#' - Population sizes (study vs target)
#' - Sampling-score distribution (quantiles of \eqn{\hat P(S=1 \mid L)})
#' - Per-group mean sampling score (overlap diagnostic)
#' - Sampling-weight distribution on study rows (mean, sd, min, max, ESS)
#' - Count of extreme sampling weights (positivity violation flag)
#'
#' @param fit A `causatr_fit` object.
#' @return A `data.table` with `statistic` and `value` columns, or
#'   `NULL` if transportability is not active.
#' @noRd
compute_sampling_diagnostics <- function(fit) {
  if (is.null(fit$target) || is.null(fit$details$sampling_model)) {
    return(NULL)
  }

  data <- fit$data
  target_col <- fit$target
  target_subset <- fit$target_subset
  sm <- fit$details$sampling_model

  s_vals <- data[[target_col]]
  study_idx <- which(s_vals == 1L)
  target_idx <- which(s_vals == 0L)
  n_study <- length(study_idx)
  n_target <- length(target_idx)

  # P(S=1|L) on all rows via the sampling model's design matrix.
  pred_terms <- stats::delete.response(stats::terms(sm$model))
  X_all <- stats::model.matrix(
    pred_terms,
    data = data,
    xlev = sm$model$xlevels
  )
  p_all <- stats::plogis(as.numeric(X_all %*% sm$gamma_hat))

  # Sampling weights on study rows.
  w_S <- compute_sampling_weights(sm, data, target_col, target_subset)
  w_study <- w_S[study_idx]

  # ESS on study rows (Kish 1965).
  ess <- if (sum(w_study^2) == 0) 0 else sum(w_study)^2 / sum(w_study^2)

  # Extreme weight count: study-row weights exceeding 4x the 99th
  # percentile signal severe positivity violations in the sampling model
  # (analogous to extreme propensity-score flags).
  q99 <- stats::quantile(w_study, 0.99)
  n_extreme <- sum(w_study > 4 * q99)

  data.table::data.table(
    statistic = c(
      "n_study",
      "n_target",
      "pct_study",
      "target_subset",
      "p_study_min",
      "p_study_q05",
      "p_study_median",
      "p_study_q95",
      "p_study_max",
      "p_study_mean_s1",
      "p_study_mean_s0",
      "sw_mean",
      "sw_sd",
      "sw_min",
      "sw_max",
      "sw_ess",
      "sw_n_extreme"
    ),
    value = c(
      n_study,
      n_target,
      round(100 * n_study / (n_study + n_target), 2),
      NA_real_,
      round(min(p_all), 4),
      round(stats::quantile(p_all, 0.05), 4),
      round(stats::median(p_all), 4),
      round(stats::quantile(p_all, 0.95), 4),
      round(max(p_all), 4),
      round(mean(p_all[study_idx]), 4),
      round(mean(p_all[target_idx]), 4),
      round(mean(w_study), 4),
      round(if (length(w_study) < 2L) NA_real_ else stats::sd(w_study), 4),
      round(min(w_study), 4),
      round(max(w_study), 4),
      round(ess, 2),
      n_extreme
    )
  )
}
