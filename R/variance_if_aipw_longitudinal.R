#' Longitudinal AIPW sandwich variance (ICE-AIPW)
#'
#' @description
#' Analytic sandwich for the ICE-AIPW estimator (Bang & Robins 2005),
#' assembled as a full stacked M-estimation sandwich
#' \eqn{V = B^{-1} M B^{-\top}/n}. The stacked estimating equation
#' \eqn{\psi(\theta)} concatenates the per-period propensity scores, the
#' per-step ICE outcome / pseudo-outcome scores, and the marginal-mean
#' equation; the bread is the numerical Jacobian of the summed score
#' (see `build_aipw_long_psi()`). Because the bread captures every
#' block-triangular cross-term -- including the baseline-pseudo-regression
#' block -- the sandwich is consistent on **both balanced and unbalanced**
#' panels (monotone dropout / censoring row-filter), unlike a
#' forward-cascade approximation of the block-triangular bread inverse.
#'
#' Derives from composition of Zivich et al. (2024, *Stat. Med.* 43:5562-5572,
#' arXiv:2306.10976) for the ICE outcome-model chain and Shook-Sa et al.
#' (2025, *Biometrics* 81(2):ujaf054, arXiv:2404.16166) for the
#' point-treatment AIPW propensity correction.
#'
#' **DR caveat:** This sandwich SE is consistent only when **both**
#' nuisance models (outcome + propensity) are correctly specified.
#' Under one-model misspecification the AIPW point estimate remains
#' consistent (DR property), but the sandwich SE is not. For DR-valid
#' variance under misspecification, use `ci_method = "bootstrap"`
#' (Rotnitzky et al. 2012, *Biometrika*).
#'
#' @param fit `causatr_fit` with `estimator = "aipw"`,
#'   `type = "longitudinal"`.
#' @param ice_aipw_results Named list of `ice_aipw_iterate()` outputs.
#' @param target_within_first Logical vector over first-time-point
#'   individuals flagging the target population.
#' @param cluster_vec Optional cluster vector for clustered SE.
#' @param trim Numeric upper-quantile weight truncation (1 = none).
#'
#' @return A named `k x k` variance-covariance matrix.
#'
#' @noRd
variance_if_aipw_longitudinal <- function(
  fit,
  ice_aipw_results,
  target_within_first,
  cluster_vec = NULL,
  trim = 1
) {
  data <- fit$data
  details <- fit$details
  int_names <- names(ice_aipw_results)
  id_col <- fit$id
  time_col <- fit$time
  ext_w <- details$weights
  time_points <- details$time_points

  first_time <- time_points[1]
  rows_first <- data[[time_col]] == first_time
  all_ids <- as.character(data[rows_first][[id_col]])
  n <- length(all_ids)
  id_to_idx <- stats::setNames(seq_len(n), all_ids)

  target <- target_within_first
  if (is.null(ext_w)) {
    w_t <- rep(1, sum(target))
  } else {
    w_first <- ext_w[rows_first]
    w_t <- w_first[target]
  }
  sum_w_target <- sum(w_t)
  if (sum_w_target <= 0) {
    rlang::abort(
      "variance_if_aipw_longitudinal(): target-population weights sum to 0.",
      class = "causatr_empty_target"
    )
  }

  cluster_first <- if (is.null(cluster_vec)) {
    NULL
  } else {
    cluster_vec[rows_first]
  }

  # Per-intervention id-level influence functions from the stacked
  # estimating equation. Propensity coefficients are shared across
  # interventions, so the cross-intervention covariance flows through the
  # shared propensity score rows when `vcov_from_if()` crossproducts the
  # IF columns.
  IF_list <- lapply(int_names, function(nm) {
    aipw_long_stacked_if(
      fit,
      ice_aipw_results[[nm]],
      all_ids,
      n,
      id_to_idx,
      target,
      w_t,
      sum_w_target,
      trim = trim
    )
  })
  names(IF_list) <- int_names

  vcov_from_if(IF_list, n, int_names, cluster = cluster_first)
}
