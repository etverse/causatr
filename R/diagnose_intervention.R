#' Compute the percentage of observations affected by an intervention
#'
#' For each observation, compares the observed treatment value to the
#' value under intervention. Returns a one-row `data.table` with
#' `n_total`, `n_changed`, and `pct_changed`. Returns `NULL` for
#' interventions where the counterfactual treatment is undefined
#' (`ipsi()`) or when the intervention is `NULL` (natural course).
#'
#' @param fit A `causatr_fit` object.
#' @param intervention A `causatr_intervention` or `NULL`.
#' @return A `data.table` or `NULL`.
#' @noRd
compute_pct_intervened <- function(fit, intervention) {
  if (is.null(intervention)) {
    return(NULL)
  }
  if (
    inherits(intervention, "causatr_intervention") &&
      identical(intervention$type, "ipsi")
  ) {
    return(NULL)
  }

  treatment <- fit$treatment
  data <- fit$data

  # For multivariate treatments, apply per-component and aggregate
  if (length(treatment) > 1L) {
    return(compute_pct_intervened_multivariate(
      data,
      treatment,
      intervention
    ))
  }

  observed <- data[[treatment]]
  data_iv <- apply_intervention(
    data.table::copy(data),
    treatment,
    intervention
  )
  intervened <- data_iv[[treatment]]

  n_total <- length(observed)
  n_changed <- sum(observed != intervened, na.rm = TRUE)

  data.table::data.table(
    n_total = n_total,
    n_changed = n_changed,
    pct_changed = round(100 * n_changed / n_total, 1)
  )
}


#' Feasibility for multivariate treatments (per-component)
#'
#' @param data data.table.
#' @param treatment Character vector of treatment column names.
#' @param intervention Named list of per-component interventions.
#' @return A `data.table` with one row per treatment component.
#' @noRd
compute_pct_intervened_multivariate <- function(
  data,
  treatment,
  intervention
) {
  data_iv <- apply_intervention(
    data.table::copy(data),
    treatment,
    intervention
  )

  rows <- lapply(treatment, function(trt) {
    observed <- data[[trt]]
    intervened <- data_iv[[trt]]
    n_total <- length(observed)
    n_changed <- sum(observed != intervened, na.rm = TRUE)
    data.table::data.table(
      component = trt,
      n_total = n_total,
      n_changed = n_changed,
      pct_changed = round(100 * n_changed / n_total, 1)
    )
  })

  data.table::rbindlist(rows)
}


#' Feasibility for longitudinal fits (per-period)
#'
#' Computes percent intervened per time period and overall.
#'
#' @param fit A longitudinal `causatr_fit`.
#' @param intervention A `causatr_intervention` or `NULL`.
#' @return A named list of `data.table`s keyed by time period, or `NULL`.
#' @noRd
compute_pct_intervened_longitudinal <- function(fit, intervention) {
  if (is.null(intervention)) {
    return(NULL)
  }
  if (
    inherits(intervention, "causatr_intervention") &&
      identical(intervention$type, "ipsi")
  ) {
    return(NULL)
  }

  treatment <- fit$treatment
  data <- fit$data
  time_col <- fit$time
  time_points <- sort(unique(data[[time_col]]))

  data_iv <- apply_intervention(
    data.table::copy(data),
    treatment,
    intervention
  )

  result <- list()
  n_changed_total <- 0L
  n_total_total <- 0L

  for (tp in time_points) {
    mask <- data[[time_col]] == tp
    obs_k <- data[[treatment]][mask]
    iv_k <- data_iv[[treatment]][mask]
    n_k <- length(obs_k)
    changed_k <- sum(obs_k != iv_k, na.rm = TRUE)

    result[[as.character(tp)]] <- data.table::data.table(
      n_total = n_k,
      n_changed = changed_k,
      pct_changed = round(100 * changed_k / n_k, 1)
    )

    n_changed_total <- n_changed_total + changed_k
    n_total_total <- n_total_total + n_k
  }

  result[["overall"]] <- data.table::data.table(
    n_total = n_total_total,
    n_changed = n_changed_total,
    pct_changed = round(100 * n_changed_total / n_total_total, 1)
  )

  result
}
