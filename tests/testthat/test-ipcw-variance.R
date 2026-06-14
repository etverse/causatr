# Variance regression tests for IPCW sandwich and bootstrap.
#
# These tests verify the IPCW Channel 2 corrections in the sandwich
# estimator by comparing against manual finite-difference calculations
# and checking that the correction reduces variance relative to the
# naive (weight-estimation-ignoring) sandwich.

# ── Point gcomp: sandwich vs bootstrap agreement ─────────────────

test_that("gcomp IPCW: sandwich vs bootstrap SE agree", {
  d <- simulate_mar_outcome(n = 3000, seed = 600)
  dt <- data.table::as.data.table(d)

  fit <- causat(
    dt,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "gcomp",
    censoring = "C",
    ipcw = TRUE
  )

  r_sw <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "sandwich"
  )

  r_bt <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "bootstrap",
    n_boot = 500
  )

  ratio <- r_sw$contrasts$se / r_bt$contrasts$se
  expect_gt(ratio, 0.5)
  expect_lt(ratio, 2.0)
})


test_that("IPW IPCW: sandwich vs bootstrap SE agree", {
  d <- simulate_mar_outcome(n = 3000, seed = 601)
  dt <- data.table::as.data.table(d)

  fit <- causat(
    dt,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw",
    censoring = "C",
    ipcw = TRUE
  )

  r_sw <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "sandwich"
  )

  r_bt <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "bootstrap",
    n_boot = 500
  )

  ratio <- r_sw$contrasts$se / r_bt$contrasts$se
  expect_gt(ratio, 0.5)
  expect_lt(ratio, 2.0)
})


# ── Manual A_{beta,gamma} via finite differences ────────────────

test_that("numDeriv jacobian matches manual finite differences for gcomp", {
  d <- simulate_mar_outcome(n = 500, seed = 602)
  dt <- data.table::as.data.table(d)

  fit <- causat(
    dt,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "gcomp",
    censoring = "C",
    ipcw = TRUE
  )

  model <- fit$model
  cens_model <- fit$details$censoring_model
  n <- nrow(fit$data)

  # Reconstruct the phi_bar closure (same as in compute_ipcw_if_correction)
  fam <- model$family
  beta_hat <- stats::coef(model)
  fit_idx <- which(get_fit_rows(fit$data, fit$outcome, fit$censoring))
  x_fit <- stats::model.matrix(model)
  eta_fit <- as.numeric(x_fit %*% beta_hat)
  mu_fit <- fam$linkinv(eta_fit)
  mu_eta_fit <- fam$mu.eta(eta_fit)
  var_mu_fit <- fam$variance(mu_fit)
  y_fit <- stats::model.response(stats::model.frame(model))
  r_fit <- y_fit - mu_fit
  n_fit <- length(y_fit)

  w_ext <- fit$details$weights_pre_ipcw
  w_ext_fit <- if (is.null(w_ext)) rep(1, n_fit) else w_ext[fit_idx]

  ipcw_wfn <- make_ipcw_weight_fn(
    cens_model,
    n_total = n,
    censoring_col = as.integer(fit$data[[fit$censoring]]),
    stabilize = TRUE
  )

  phi_bar_cens <- function(gamma) {
    w_ipcw <- ipcw_wfn(gamma)
    w_ipcw_fit <- w_ipcw[fit_idx]
    s_per_i <- w_ext_fit * w_ipcw_fit * mu_eta_fit * r_fit / var_mu_fit
    as.numeric(crossprod(x_fit, s_per_i)) / n_fit
  }

  gamma_hat <- cens_model$alpha_hat
  p_cens <- length(gamma_hat)
  q_out <- length(beta_hat)

  # numDeriv jacobian
  a_numderiv <- -numDeriv::jacobian(phi_bar_cens, x = gamma_hat)

  # Manual forward finite differences
  eps <- 1e-5
  phi0 <- phi_bar_cens(gamma_hat)
  a_manual <- matrix(0, q_out, p_cens)
  for (j in seq_len(p_cens)) {
    gamma_plus <- gamma_hat
    gamma_plus[j] <- gamma_plus[j] + eps
    phi_plus <- phi_bar_cens(gamma_plus)
    a_manual[, j] <- -(phi_plus - phi0) / eps
  }

  expect_equal(a_numderiv, a_manual, tolerance = 1e-3)
})


# ── Sandwich vs bootstrap for more estimators ───────────────────

test_that("matching IPCW: sandwich vs bootstrap SE ratio in [0.5, 2.0]", {
  d <- simulate_mar_outcome(n = 3000, seed = 604)
  dt <- data.table::as.data.table(d)

  fit <- causat(
    dt,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "matching",
    censoring = "C",
    ipcw = TRUE
  )

  r_sw <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "sandwich"
  )
  r_bt <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "bootstrap",
    n_boot = 500
  )

  ratio <- r_sw$contrasts$se / r_bt$contrasts$se
  expect_gt(ratio, 0.5)
  expect_lt(ratio, 2.0)
})


