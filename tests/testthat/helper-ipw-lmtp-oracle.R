# lmtp + manual IPSI contrast-level oracle helpers for non-static IPW.
#
# These helpers power the T-oracle5 and T-oracle6 tests in
# test-ipw-lmtp-oracle.R. lmtp is a point-estimate oracle only — its SE
# uses the EIF, not the M-estimation sandwich, so only point estimates
# are compared. The manual IPSI oracle computes Kennedy (2019) closed-
# form weights from the same propensity model causatr fits internally.
# See PHASE_4_INTERVENTIONS_SELF_IPW.md §9 for the design rationale.

#' Run lmtp::lmtp_sdr() for a point-treatment shift and extract E[Y(A+delta)]
#'
#' @description
#' Wraps `lmtp::lmtp_sdr()` (sequentially doubly robust) with
#' `learners_trt = "SL.glm"`, `learners_outcome = "SL.glm"`, and
#' `folds = 1` (no cross-fitting) to match causatr's parametric GLM
#' setup on the same data. Returns the marginal mean point estimate
#' under the shifted treatment.
#'
#' `lmtp_ipw()` was removed in lmtp >= 1.5.0 (IPW requires correctly
#' pre-specified parametric models). Under correct specification of
#' both treatment and outcome GLMs (guaranteed by the linear DGP),
#' SDR is consistent for the same target parameter as causatr's
#' density-ratio IPW, making it a valid point-estimate oracle.
#'
#' @param data Data frame with outcome, treatment, and covariate columns.
#' @param outcome Character scalar, outcome column name.
#' @param treatment Character scalar, treatment column name (continuous).
#' @param baseline Character vector of baseline covariate names.
#' @param delta Numeric scalar, the shift magnitude (A* = A + delta).
#'
#' @return A list with:
#'   - `theta` — point estimate E[Y(A + delta)].
#'   - `lmtp_obj` — raw lmtp fit for debugging on failure.
lmtp_oracle_shift <- function(data, outcome, treatment, baseline, delta) {
  shift_fn <- function(data, trt) {
    data[[trt]] + delta
  }

  fit <- lmtp::lmtp_sdr(
    data = data,
    trt = treatment,
    outcome = outcome,
    baseline = baseline,
    shift = shift_fn,
    outcome_type = "continuous",
    learners_trt = "SL.glm",
    learners_outcome = "SL.glm",
    folds = 1
  )

  # lmtp >= 1.5 stores the point estimate in an S4 `influence_func_estimate`
  # object at `fit$estimate@x`.
  list(
    theta = fit$estimate@x,
    lmtp_obj = fit
  )
}


#' Compute the IPSI Hájek weighted mean by hand (Kennedy 2019)
#'
#' @description
#' Fits a logistic propensity model `A ~ confounders` and computes the
#' Kennedy (2019) closed-form IPSI weights
#'   w_i = (delta * A_i + (1 - A_i)) / (delta * p_i + (1 - p_i))
#' then returns the Hájek weighted mean `sum(w * Y) / sum(w)` as the
#' oracle estimate of E[Y(ipsi(delta))].
#'
#' This is a self-contained first-principles computation with no
#' dependency on causatr internals, so agreement between this oracle
#' and `contrast(fit, ipsi(delta))` proves the IPSI weight path and
#' the intercept-only MSM are wired correctly.
#'
#' @param data Data frame with outcome, treatment, and covariate columns.
#' @param outcome Character scalar, outcome column name.
#' @param treatment Character scalar, treatment column name (binary 0/1).
#' @param ps_formula One-sided formula for the propensity-score model
#'   (e.g. `~ L`).
#' @param delta Numeric scalar > 0, the IPSI odds-ratio multiplier.
#'
#' @return A list with:
#'   - `theta` — Hájek weighted mean under IPSI weights.
#'   - `weights` — the raw IPSI weight vector.
#'   - `ps_model` — the fitted propensity GLM.
manual_ipsi_oracle <- function(data, outcome, treatment, ps_formula, delta) {
  A <- data[[treatment]]
  Y <- data[[outcome]]

  # Fit the same logistic propensity model causatr uses internally
  # (default propensity_model_fn = stats::glm, family = binomial).
  ps_rhs <- paste(deparse(ps_formula[[length(ps_formula)]]), collapse = " ")
  ps_full <- stats::as.formula(paste(treatment, "~", ps_rhs))
  ps_model <- stats::glm(ps_full, data = data, family = stats::binomial())
  p <- unname(stats::predict(ps_model, type = "response"))

  # Kennedy (2019) closed-form IPSI weight. The intervention shifts the
  # propensity odds by delta: p*(L) = delta*p(L) / (1 - p(L) + delta*p(L)).
  # The resulting density-ratio weight under the Bernoulli model is:
  #   w_i = (delta * A_i + (1 - A_i)) / (delta * p_i + (1 - p_i))
  w <- (delta * A + (1 - A)) / (delta * p + (1 - p))

  # Hájek weighted mean: E[Y(ipsi(delta))] = sum(w*Y) / sum(w).
  theta <- sum(w * Y) / sum(w)

  list(
    theta = theta,
    weights = w,
    ps_model = ps_model
  )
}


#' Compute the longitudinal IPSI Hájek mean by hand (per-period Kennedy product)
#'
#' @description
#' Longitudinal companion to `manual_ipsi_oracle()`. Fits a separate
#' logistic propensity model per period, computes the Kennedy (2019)
#' closed-form IPSI weight at each period, takes the per-individual
#' product across periods, and returns the Hájek weighted mean of the
#' final-period outcome:
#'   W_i = prod_t (delta * A_{t,i} + (1 - A_{t,i})) /
#'                (delta * p_{t,i} + (1 - p_{t,i}))
#'   theta = sum(W * Y_final) / sum(W).
#'
#' This is a self-contained first-principles computation with no
#' dependency on causatr internals, so agreement between this oracle and
#' `contrast(fit, ipsi(delta))` on a longitudinal IPW fit proves the
#' per-period IPSI weight product, the cumulative weight alignment, and
#' the intercept-only final-period MSM are wired correctly.
#'
#' Callers supply one data frame per period (one row per individual,
#' carrying every column the corresponding per-period propensity formula
#' references, including any `lag*_*` columns) and a matching list of
#' two-sided propensity formulas. The order of `data_by_period` and
#' `formulas` must agree, and every per-period data frame must share the
#' same individual ordering as `y_final`.
#'
#' @param data_by_period List of per-period data frames (one row per
#'   individual at that period), ordered earliest to latest.
#' @param formulas List of two-sided propensity formulas (`A ~ ...`),
#'   one per period, aligned with `data_by_period`.
#' @param treatment Character scalar, treatment column name (binary 0/1).
#' @param y_final Numeric vector of final-period outcomes, ordered to
#'   match the rows of each `data_by_period` element.
#' @param delta Numeric scalar > 0, the IPSI odds-ratio multiplier.
#'
#' @return A list with:
#'   - `theta` — Hájek weighted mean under the per-period IPSI product.
#'   - `weights` — the cumulative per-individual IPSI weight vector.
#'   - `ps_models` — list of the fitted per-period propensity GLMs.
manual_long_ipsi_oracle <- function(
  data_by_period,
  formulas,
  treatment,
  y_final,
  delta
) {
  n_periods <- length(data_by_period)
  n_id <- nrow(data_by_period[[1L]])

  # Per-period Kennedy weight, accumulated into a per-individual product.
  W <- rep(1, n_id)
  ps_models <- vector("list", n_periods)
  for (t in seq_len(n_periods)) {
    dat_t <- data_by_period[[t]]
    a_t <- dat_t[[treatment]]
    ps_models[[t]] <- stats::glm(
      formulas[[t]],
      data = dat_t,
      family = stats::binomial()
    )
    p_t <- unname(stats::predict(ps_models[[t]], type = "response"))
    w_t <- (delta * a_t + (1 - a_t)) / (delta * p_t + (1 - p_t))
    W <- W * w_t
  }

  theta <- sum(W * y_final) / sum(W)

  list(
    theta = theta,
    weights = W,
    ps_models = ps_models
  )
}
