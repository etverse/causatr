#' Censoring diagnostics for IPCW fits
#'
#' @description
#' Computes a summary table for the censoring model and IPCW weights.
#' Returns `NULL` when the fit does not use built-in IPCW.
#'
#' The returned table includes:
#' - IPCW weight distribution (mean, sd, min, max, ESS among uncensored)
#' - Censoring model positivity (quantiles of P(C=0|A,L))
#' - Censoring prevalence
#'
#' For longitudinal fits, returns a named list of per-period tables
#' plus a `"cumulative"` entry for the full cumulative weights.
#'
#' @param fit A `causatr_fit` object.
#' @return A `data.table` (point) or named list of `data.table`s
#'   (longitudinal), or `NULL` if IPCW is not active.
#' @noRd
compute_censoring_diagnostics <- function(fit) {
  if (!isTRUE(fit$details$ipcw)) {
    return(NULL)
  }

  data <- fit$data
  censoring <- fit$censoring

  if (fit$type == "longitudinal") {
    return(compute_censoring_diagnostics_long(fit))
  }

  ipcw_w <- fit$details$ipcw_weights
  cens_col <- as.integer(data[[censoring]])
  uncens <- cens_col == 0L & !is.na(cens_col)

  w_uncens <- ipcw_w[uncens]
  n_total <- nrow(data)
  n_cens <- sum(cens_col == 1L, na.rm = TRUE)

  cm <- fit$details$censoring_model
  p_uncens <- cm$p_uncensored

  ess <- sum(w_uncens)^2 / sum(w_uncens^2)

  data.table::data.table(
    statistic = c(
      "n_total",
      "n_censored",
      "pct_censored",
      "ipcw_mean",
      "ipcw_sd",
      "ipcw_min",
      "ipcw_max",
      "ipcw_ess",
      "p_uncens_min",
      "p_uncens_q05",
      "p_uncens_median",
      "p_uncens_q95",
      "p_uncens_max"
    ),
    value = round(
      c(
        n_total,
        n_cens,
        100 * n_cens / n_total,
        mean(w_uncens),
        stats::sd(w_uncens),
        min(w_uncens),
        max(w_uncens),
        ess,
        min(p_uncens),
        stats::quantile(p_uncens, 0.05),
        stats::median(p_uncens),
        stats::quantile(p_uncens, 0.95),
        max(p_uncens)
      ),
      4
    )
  )
}


#' Per-period censoring diagnostics for longitudinal IPCW
#'
#' @param fit A longitudinal `causatr_fit` with `details$ipcw`.
#' @return Named list of `data.table`s (one per period + "cumulative").
#' @noRd
compute_censoring_diagnostics_long <- function(fit) {
  data <- fit$data
  censoring <- fit$censoring
  time_col <- fit$time
  time_points <- fit$details$time_points
  cens_models <- fit$details$censoring_models
  ppw <- fit$details$per_period_weights
  ipcw_w <- fit$details$ipcw_weights
  cens_col <- as.integer(data[[censoring]])

  result <- list()

  for (k in seq_along(time_points)) {
    rows_k <- data[[time_col]] == time_points[k]
    cens_k <- cens_col[rows_k]
    n_k <- sum(rows_k)
    n_cens_k <- sum(cens_k == 1L, na.rm = TRUE)

    w_k <- ppw[[k]]
    uncens_k <- cens_k == 0L & !is.na(cens_k)
    w_unc_k <- w_k[uncens_k]

    cm_k <- cens_models[[k]]

    if (n_cens_k == 0L || is.null(cm_k)) {
      result[[as.character(time_points[k])]] <-
        data.table::data.table(
          statistic = c(
            "n_total",
            "n_censored",
            "pct_censored"
          ),
          value = c(n_k, 0, 0)
        )
      next
    }

    p_unc_k <- cm_k$p_uncensored
    ess_k <- sum(w_unc_k)^2 / sum(w_unc_k^2)

    result[[as.character(time_points[k])]] <-
      data.table::data.table(
        statistic = c(
          "n_total",
          "n_censored",
          "pct_censored",
          "ipcw_mean",
          "ipcw_sd",
          "ipcw_min",
          "ipcw_max",
          "ipcw_ess",
          "p_uncens_min",
          "p_uncens_q05",
          "p_uncens_median",
          "p_uncens_q95",
          "p_uncens_max"
        ),
        value = round(
          c(
            n_k,
            n_cens_k,
            100 * n_cens_k / n_k,
            mean(w_unc_k),
            stats::sd(w_unc_k),
            min(w_unc_k),
            max(w_unc_k),
            ess_k,
            min(p_unc_k),
            stats::quantile(p_unc_k, 0.05),
            stats::median(p_unc_k),
            stats::quantile(p_unc_k, 0.95),
            max(p_unc_k)
          ),
          4
        )
      )
  }

  # Cumulative weights across all periods
  uncens_all <- cens_col == 0L & !is.na(cens_col)
  w_all <- ipcw_w[uncens_all]
  n_cens_all <- sum(cens_col == 1L, na.rm = TRUE)
  ess_cum <- sum(w_all)^2 / sum(w_all^2)

  result[["cumulative"]] <- data.table::data.table(
    statistic = c(
      "n_total",
      "n_censored",
      "pct_censored",
      "ipcw_mean",
      "ipcw_sd",
      "ipcw_min",
      "ipcw_max",
      "ipcw_ess"
    ),
    value = round(
      c(
        nrow(data),
        n_cens_all,
        100 * n_cens_all / nrow(data),
        mean(w_all),
        stats::sd(w_all),
        min(w_all),
        max(w_all),
        ess_cum
      ),
      4
    )
  )

  result
}
