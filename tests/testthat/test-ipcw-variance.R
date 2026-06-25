# Variance regression tests for IPCW sandwich and bootstrap.
#
# These tests verify the IPCW Channel 2 corrections in the sandwich
# estimator by comparing against manual finite-difference calculations
# and checking that the correction reduces variance relative to the
# naive (weight-estimation-ignoring) sandwich.

# ── Point gcomp: sandwich vs bootstrap agreement ─────────────────

test_that("gcomp IPCW: sandwich vs bootstrap SE agree", {
  skip_if_fast()
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


# Tight oracle: the IPCW gcomp sandwich must equal the full estimated-gamma
# M-estimation sandwich, including the censoring model's contribution to the
# marginal-mean variance through BOTH the outcome model (indirect) and the
# IPCW-weighted average (direct). We build that sandwich here independently of
# causatr's channel-wise IF assembly -- stacking the logistic censoring score,
# the IPCW-weighted Gaussian outcome score, and the IPCW-weighted marginal-mean
# equations into theta = (gamma, beta, mu_a1, mu_a0), then reading the sandwich
# off the EE Jacobian (numerical bread) and the empirical meat. Reusing
# causatr's fitted nuisance coefficients (they solve the EE) isolates the
# variance computation: a different recipe (bread^-1 meat bread^-T) must land on
# the same SEs. The direct term is what makes the per-class MEAN SEs agree; the
# difference SE agrees with or without it (the term cancels in contrasts).
test_that("gcomp IPCW sandwich matches the estimated-gamma M-estimation stack", {
  d <- simulate_mar_outcome(n = 2500, seed = 909)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "gcomp",
    censoring = "C",
    ipcw = TRUE
  )
  sw <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "sandwich"
  )

  # --- independent estimated-gamma M-estimation sandwich --------------------
  n <- nrow(d)
  unc <- as.numeric(d$C == 0)
  Xc <- cbind(1, d$A, d$L) # censoring design [1, A, L]
  Xb <- cbind(1, d$A, d$L) # Gaussian outcome design [1, A, L]
  Xa1 <- cbind(1, 1, d$L)
  Xa0 <- cbind(1, 0, d$L)
  y0 <- ifelse(unc == 1, d$Y, 0) # zeroed on censored rows (masked by w)
  gamma_hat <- fit$details$censoring_model$alpha_hat
  beta_hat <- stats::coef(fit$model)

  expit <- function(x) 1 / (1 + exp(-x))
  stacked <- function(th) {
    g <- th[1:3]
    b <- th[4:6]
    mu1 <- th[7]
    mu0 <- th[8]
    w <- (1 / expit(Xc %*% g)) * unc # unstabilized; p_marg cancels
    cr <- unc - expit(Xc %*% g)
    rb <- y0 - Xb %*% b # Gaussian residual (identity link)
    cbind(
      Xc * as.numeric(cr),
      Xb * as.numeric(w * rb),
      w * (mu1 - Xa1 %*% b),
      w * (mu0 - Xa0 %*% b)
    )
  }
  # mu at the fitted coefficients (IPCW-weighted g-formula average).
  w_hat <- (1 / expit(Xc %*% gamma_hat)) * unc
  sw_tot <- sum(w_hat)
  mu1 <- sum(w_hat * (Xa1 %*% beta_hat)) / sw_tot
  mu0 <- sum(w_hat * (Xa0 %*% beta_hat)) / sw_tot
  theta <- c(gamma_hat, beta_hat, mu1, mu0)

  q <- length(theta)
  eps <- 1e-6
  J <- matrix(0, q, q)
  ee <- function(th) colSums(stacked(th))
  for (j in seq_len(q)) {
    tp <- theta
    tm <- theta
    tp[j] <- tp[j] + eps
    tm[j] <- tm[j] - eps
    J[, j] <- (ee(tp) - ee(tm)) / (2 * eps)
  }
  bread_inv <- solve(-J / n)
  psi <- stacked(theta)
  meat <- crossprod(psi) / n
  V <- bread_inv %*% meat %*% t(bread_inv) / n
  ref_se <- c(a1 = sqrt(V[7, 7]), a0 = sqrt(V[8, 8]))
  ref_diff_se <- sqrt(V[7, 7] + V[8, 8] - 2 * V[7, 8])

  # causatr's per-class mean SEs and the difference SE must match the stack.
  expect_equal(sw$estimates$estimate, c(mu1, mu0), tolerance = 1e-6)
  expect_equal(unname(sw$estimates$se), unname(ref_se), tolerance = 1e-4)
  expect_equal(sw$contrasts$se, ref_diff_se, tolerance = 1e-4)
})


# The direct dmu/dgamma term must be active and materially reduce the
# marginal-mean SE relative to the indirect-only (known-gamma) sandwich -- a
# guard that the efficiency-gain channel is wired, not silently zero. Flipping
# the `ipcw` flag off suppresses the whole Channel-3 correction; the point
# estimates are invariant, the per-class mean SEs shrink under the full sandwich,
# and the difference SE is unchanged (the term cancels in contrasts).
test_that("gcomp IPCW direct-path term is active and cancels in contrasts", {
  d <- simulate_mar_outcome(n = 2500, seed = 911)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "gcomp",
    censoring = "C",
    ipcw = TRUE
  )
  ivs <- list(a1 = static(1), a0 = static(0))
  full <- contrast(
    fit,
    interventions = ivs,
    type = "difference",
    ci_method = "sandwich"
  )

  fit_known <- fit
  fit_known$details$ipcw <- FALSE
  known <- contrast(
    fit_known,
    interventions = ivs,
    type = "difference",
    ci_method = "sandwich"
  )

  expect_equal(
    full$estimates$estimate,
    known$estimates$estimate,
    tolerance = 1e-10
  )
  # Estimating the censoring model reduces the marginal-mean SE.
  ratio <- full$estimates$se / known$estimates$se
  expect_true(all(ratio <= 1 + 1e-8))
  expect_true(min(ratio) < 0.99)
  # The reduction cancels in the difference contrast.
  expect_equal(full$contrasts$se, known$contrasts$se, tolerance = 1e-3)
})


# Survey/external weights compose with IPCW: the censoring model is survey-
# weighted, the outcome model and the marginal average carry survey x IPCW
# weights, and the direct-path dw/dgamma keeps the full (survey x IPCW) weight.
# The independent estimated-gamma M-estimation stack carries the same survey
# weight in every block, so the weighted IPCW sandwich must match it.
test_that("weighted gcomp IPCW sandwich matches the weighted M-estimation stack", {
  d <- simulate_mar_outcome(n = 2500, seed = 912)
  set.seed(99)
  d$sw <- stats::runif(nrow(d), 0.5, 2)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "gcomp",
    censoring = "C",
    ipcw = TRUE,
    weights = d$sw
  )
  sw <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "sandwich"
  )

  n <- nrow(d)
  unc <- as.numeric(d$C == 0)
  svy <- d$sw
  Xc <- cbind(1, d$A, d$L)
  Xb <- cbind(1, d$A, d$L)
  Xa1 <- cbind(1, 1, d$L)
  Xa0 <- cbind(1, 0, d$L)
  y0 <- ifelse(unc == 1, d$Y, 0)
  gamma_hat <- fit$details$censoring_model$alpha_hat
  beta_hat <- stats::coef(fit$model)
  expit <- function(x) 1 / (1 + exp(-x))

  stacked <- function(th) {
    g <- th[1:3]
    b <- th[4:6]
    mu1 <- th[7]
    mu0 <- th[8]
    w <- svy * (1 / expit(Xc %*% g)) * unc # total = survey x IPCW
    cr <- unc - expit(Xc %*% g)
    rb <- y0 - Xb %*% b
    cbind(
      Xc * as.numeric(svy * cr), # survey-weighted censoring score
      Xb * as.numeric(w * rb),
      w * (mu1 - Xa1 %*% b),
      w * (mu0 - Xa0 %*% b)
    )
  }
  w_hat <- svy * (1 / expit(Xc %*% gamma_hat)) * unc
  sw_tot <- sum(w_hat)
  mu1 <- sum(w_hat * (Xa1 %*% beta_hat)) / sw_tot
  mu0 <- sum(w_hat * (Xa0 %*% beta_hat)) / sw_tot
  theta <- c(gamma_hat, beta_hat, mu1, mu0)

  q <- length(theta)
  eps <- 1e-6
  J <- matrix(0, q, q)
  ee <- function(th) colSums(stacked(th))
  for (j in seq_len(q)) {
    tp <- theta
    tm <- theta
    tp[j] <- tp[j] + eps
    tm[j] <- tm[j] - eps
    J[, j] <- (ee(tp) - ee(tm)) / (2 * eps)
  }
  bread_inv <- solve(-J / n)
  meat <- crossprod(stacked(theta)) / n
  V <- bread_inv %*% meat %*% t(bread_inv) / n
  ref_se <- c(sqrt(V[7, 7]), sqrt(V[8, 8]))
  ref_diff_se <- sqrt(V[7, 7] + V[8, 8] - 2 * V[7, 8])

  expect_equal(sw$estimates$estimate, c(mu1, mu0), tolerance = 1e-6)
  expect_equal(unname(sw$estimates$se), ref_se, tolerance = 1e-4)
  expect_equal(sw$contrasts$se, ref_diff_se, tolerance = 1e-4)
})


# Regression guard for the "no direct-term gap" finding on point IPW: the Hajek
# MSM mean IS the model parameter, so its gamma-dependence is already in the
# indirect cross-term and the sandwich is well-calibrated WITHOUT a separate
# direct term. Monte-Carlo: the mean sandwich SE must track the empirical
# sampling SD (unlike the plug-in g-formula, which needs the direct term).
test_that("IPW IPCW marginal-mean sandwich tracks the empirical SD (no gap)", {
  skip_if_fast()
  gen <- function(n) {
    L <- stats::rnorm(n)
    A <- stats::rbinom(n, 1, stats::plogis(0.4 * L))
    Yb <- stats::rbinom(n, 1, stats::plogis(-0.3 + 0.9 * A + 0.6 * L))
    C <- stats::rbinom(n, 1, stats::plogis(-0.6 + 0.5 * A + 0.7 * L))
    Yo <- Yb
    Yo[C == 1] <- NA
    data.frame(L = L, A = A, Y = Yo, C = C)
  }
  ivs <- list(a1 = static(1), a0 = static(0))
  R <- 300L
  nn <- 2500L
  est <- matrix(NA_real_, R, 2L)
  se <- matrix(NA_real_, R, 2L)
  set.seed(5)
  for (r in seq_len(R)) {
    dat <- gen(nn)
    fit <- tryCatch(
      causat(
        dat,
        "Y",
        "A",
        ~L,
        estimator = "ipw",
        censoring = "C",
        ipcw = TRUE
      ),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      next
    }
    sw <- tryCatch(
      contrast(
        fit,
        interventions = ivs,
        type = "difference",
        ci_method = "sandwich"
      ),
      error = function(e) NULL
    )
    if (is.null(sw)) {
      next
    }
    est[r, ] <- sw$estimates$estimate
    se[r, ] <- sw$estimates$se
  }
  ok <- stats::complete.cases(est)
  emp_sd <- apply(est[ok, ], 2L, stats::sd)
  mean_se <- colMeans(se[ok, ])
  # Two-sided: the IPW sandwich SE neither over- nor under-states the truth.
  expect_equal(mean_se / emp_sd, c(1, 1), tolerance = 0.06)
})


test_that("IPW IPCW: sandwich vs bootstrap SE agree", {
  skip_if_fast()
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
  skip_if_fast()
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
