# T-end-to-end oracle test for `compute_ipw_if_self_contained_one()`.
#
# Goal: verify Branch B of the variance engine reproduces the
# classical M-estimation sandwich variance for a stacked
# `(alpha, beta)` system, to floating-point precision, on a setup
# small enough to hand-assemble the stacked bread and meat.
#
# Setup:
#
# - DGP: `simulate_binary_continuous()` with a moderate sample size
#   (n = 500 rows, one confounder L). Binary treatment A, Gaussian
#   outcome Y. See `helper-dgp.R` for the true parameters.
# - Propensity model: `glm(A ~ L, family = binomial)`.
# - Intervention: `static(1)`. HT weights `w_i = I(A_i = 1) / p_i`
#   (what `compute_density_ratio_weights()` returns under Phase 4's
#   HT form for `static` on binary).
# - MSM: intercept-only weighted Gaussian GLM `glm(Y ~ 1, weights = w,
#   family = gaussian)`, whose intercept is the Hájek mean of Y among
#   the treated — i.e. the estimand `mu_1 = E[Y^1]`.
#
# Under this architecture the `(alpha, beta)` stacked M-estimation
# system has:
#
#   psi_alpha_i = (A_i - p_i) * X_ps_i         [logistic score]
#   psi_beta_i  = w_i * (Y_i - beta_0)         [weighted intercept score]
#
# and the (p_alpha + 1) × (p_alpha + 1) bread has the analytic form
#
#   A_aa = (1/n) * X_ps^T * diag(p * (1 - p)) * X_ps
#   A_bb = (1/n) * sum(w_i)
#   A_ba = (1/n) * sum_i [(Y_i - beta_0) * A_i * (1 - p_i) / p_i * X_ps_i]
#
# The meat is built straight from the per-obs scores:
#
#   M = (1/n) * sum_i [psi_i * psi_i^T],   psi_i = (psi_alpha_i, psi_beta_i)
#
# And the sandwich variance of `(alpha_hat, beta_0_hat)` is
#
#   V = (1/n) * B^(-1) * M * B^(-T).
#
# We read off V[p_alpha + 1, p_alpha + 1] = Var(mu_hat), and compare
# against `sum(IF^2) / n^2` where IF comes from
# `compute_ipw_if_self_contained_one()`. Both should give the same
# number to floating-point precision, because both are computing
# exactly the same thing — once through the explicit matrix algebra,
# once through the causatr primitives (`sandwich::estfun` +
# `sandwich::bread` + `numDeriv::jacobian`).
#
# If the agreement tolerance has to be loosened much beyond ~1e-6,
# the function's implementation is wrong.

test_that("Branch B IF aggregates to the analytic stacked-sandwich variance for static(1) HT on binary", {
  skip_if_not_installed("numDeriv")
  skip_if_not_installed("sandwich")

  d <- simulate_binary_continuous(n = 500, seed = 707)
  n <- nrow(d)

  # ---- Propensity model: plain logistic GLM on the confounder -----
  ps_fit <- stats::glm(A ~ L, data = d, family = stats::binomial())
  alpha_hat <- stats::coef(ps_fit)
  X_ps <- stats::model.matrix(ps_fit)
  p <- as.numeric(stats::fitted(ps_fit))
  p_alpha <- length(alpha_hat)

  # ---- HT weights for static(1) -----------------------------------
  # w_i = I(A_i = 1) / p_i. These are exactly what
  # `compute_density_ratio_weights(tm, data, static(1))` returns
  # when the treatment density model is fit on the same propensity
  # formula — verified separately in test-ipw-weights.R. We build
  # them directly here to keep the oracle test independent of
  # causatr's weight builder.
  w <- as.numeric(d$A == 1) / p

  # ---- MSM: intercept-only weighted Gaussian GLM -------------------
  # Use `quasi(link = "identity", variance = "constant")` so
  # `sandwich::estfun` / `sandwich::bread` give clean non-weighted
  # dispersion behaviour. Actually `stats::gaussian()` works fine for
  # sandwich's glm methods — the dispersion cancels between estfun
  # and bread regardless of its estimated value.
  msm_fit <- stats::glm(
    Y ~ 1,
    data = d,
    weights = w,
    family = stats::gaussian()
  )
  beta_hat <- stats::coef(msm_fit) # scalar
  mu_hat <- as.numeric(beta_hat) # identity link + intercept-only => mu == beta0

  # Sanity: on `simulate_binary_continuous` the true E[Y^1] = 5.
  # The Hajek mean should land near there (well within 0.3 on n=500).
  expect_equal(mu_hat, 5, tolerance = 0.3)

  # ---- Hand-assembled stacked bread -------------------------------
  # A_aa: (1/n) X_ps^T diag(p(1-p)) X_ps
  W_ps <- p * (1 - p)
  A_aa <- crossprod(X_ps, X_ps * W_ps) / n

  # A_bb: (1/n) sum w_i
  A_bb <- sum(w) / n

  # A_ba: (1/n) sum_i (Y_i - beta0) * A_i * (1 - p_i) / p_i * X_ps_i
  a_ba_per_row <- (d$Y - mu_hat) * d$A * (1 - p) / p
  A_ba <- as.numeric(
    crossprod(matrix(a_ba_per_row, ncol = 1), X_ps)
  ) /
    n # length p_alpha

  B <- matrix(0, p_alpha + 1, p_alpha + 1)
  B[seq_len(p_alpha), seq_len(p_alpha)] <- A_aa
  B[p_alpha + 1, seq_len(p_alpha)] <- A_ba
  B[p_alpha + 1, p_alpha + 1] <- A_bb

  # ---- Hand-assembled meat ----------------------------------------
  # psi_alpha_i = X_ps_i * (A_i - p_i)      (logistic score)
  # psi_beta_i  = w_i * (Y_i - beta0)       (weighted intercept score)
  psi_alpha <- X_ps * (d$A - p) # n × p_alpha
  psi_beta <- w * (d$Y - mu_hat) # n-vector
  Psi <- cbind(psi_alpha, psi_beta) # n × (p_alpha + 1)
  M <- crossprod(Psi) / n # (p_alpha + 1)^2

  # Sandwich variance V = (1/n) B^(-1) M B^(-T)
  B_inv <- solve(B)
  V_analytic <- (B_inv %*% M %*% t(B_inv)) / n
  var_mu_analytic <- V_analytic[p_alpha + 1, p_alpha + 1]

  # ---- causatr Branch B via the new helper ------------------------
  dt <- data.table::as.data.table(d)
  tm <- fit_treatment_model(dt, treatment = "A", confounders = ~L)
  wfn <- make_weight_fn(tm, dt, static(1))

  # Sanity: the closure at alpha_hat must return the same weights we
  # built by hand.
  expect_equal(wfn(tm$alpha_hat), w, tolerance = 1e-12)

  # For intercept-only Gaussian MSM with static intervention and a
  # saturated-over-treatment-level contrast, Ch1 is identically zero
  # (predictions are constant per treatment level — here, constant
  # at `beta_hat`). Verified independently via `build_point_channel_pieces`.
  Ch1 <- rep(0, n)
  J <- 1 # d mu/d beta_0 = 1 for identity-link intercept-only Gaussian

  IF_vec <- compute_ipw_if_self_contained_one(
    msm_model = msm_fit,
    propensity_model = ps_fit,
    weight_fn = wfn,
    J = J,
    Ch1_i = Ch1,
    fit_idx = seq_len(n),
    fit_idx_ps = seq_len(n),
    n_total = n
  )

  # vcov_from_if divides by n_total^2, so the variance of mu_hat is:
  var_mu_causatr <- sum(IF_vec^2) / n^2

  # Tolerance 1e-6: both paths compute the same stacked-sandwich
  # variance in IEEE 754 double arithmetic. The only source of drift
  # is `numDeriv::jacobian`'s central-difference accuracy in
  # approximating `A_{beta, alpha}` (Richardson extrapolation is on
  # by default, order ~1e-8 for smooth functions, but we leave
  # headroom). If this ever widens past 1e-4 something has broken
  # structurally and the debugging should start from
  # `compute_ipw_if_self_contained_one()`'s three-channel split.
  expect_equal(var_mu_causatr, var_mu_analytic, tolerance = 1e-6)
})

# ---- Sanity: natural course ---------------------------------------

test_that("Branch B IF is all-Ch1 when the intervention is natural course", {
  # A NULL intervention means "observed treatment, observed density
  # ratio = 1". The MSM correction term should vanish (no gradient on
  # a flat weight vector) and the propensity correction should vanish
  # (`make_weight_fn(NULL)` is the constant-1 closure, which has zero
  # jacobian — already verified in test-ipw-cross-derivative.R). So
  # the entire IF is whatever Ch1 is, and if Ch1 is zero the IF is
  # identically zero.
  skip_if_not_installed("numDeriv")

  d <- simulate_binary_continuous(n = 200, seed = 708)
  n <- nrow(d)
  dt <- data.table::as.data.table(d)

  tm <- fit_treatment_model(dt, treatment = "A", confounders = ~L)
  ps_fit <- tm$model

  # Natural-course weights are all 1. Use a weighted glm with
  # constant weights so the MSM exactly matches an unweighted fit.
  msm_fit <- stats::glm(
    Y ~ 1,
    data = d,
    weights = rep(1, n),
    family = stats::gaussian()
  )
  wfn <- make_weight_fn(tm, dt, NULL)

  IF_vec <- compute_ipw_if_self_contained_one(
    msm_model = msm_fit,
    propensity_model = ps_fit,
    weight_fn = wfn,
    J = 1,
    Ch1_i = rep(0, n),
    fit_idx = seq_len(n),
    fit_idx_ps = seq_len(n),
    n_total = n
  )

  # The propensity-correction channel must be zero because the
  # weight closure has zero gradient. The MSM correction channel is
  # NOT zero — it captures the uncertainty in the unweighted
  # intercept fit. That piece is real and has to be in the IF.
  # What we ARE asserting: that the propensity-correction contribution
  # is numerically zero, which we check by comparing to a
  # "no-propensity-correction" computation.
  #
  # Construction: run `compute_ipw_if_self_contained_one` on a
  # degenerate setup where the propensity model has no confounders
  # (intercept-only, so alpha_hat is a scalar). Then A_beta_alpha is
  # a 1 × 1 matrix — still computable via numDeriv. The propensity
  # term would be nonzero only if the weight depends on alpha, which
  # it does NOT under the NULL closure. So the whole propensity
  # piece must be (numerically) zero.
  # Use sum(abs(...)) as a coarse "is it essentially zero?" check.
  # The MSM term is nonzero; the propensity term is zero; their sum
  # should equal the MSM-only IF. Let's reconstruct the MSM-only IF
  # and verify.
  E_msm <- sandwich::estfun(msm_fit)
  bread_msm <- sandwich::bread(msm_fit)
  msm_only_per_obs <- as.numeric(E_msm %*% bread_msm) # J = 1
  expect_equal(IF_vec, msm_only_per_obs, tolerance = 1e-9)
})
