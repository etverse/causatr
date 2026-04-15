# Hand-derived cross-derivative tests for the Phase 4 IPW weight closures.
#
# The variance engine's Branch B body will compute the cross-derivative
#
#   A_{beta, alpha} = d/d alpha (average MSM score at fitted beta)
#
# numerically via `numDeriv::jacobian()` applied to a function whose only
# argument is `alpha`. The function in question factors as "weighted MSM
# score given weights w(alpha)", and all the `alpha` dependence flows
# through `w(alpha)`. So the correctness of the numerical cross-derivative
# reduces to the correctness of `numDeriv::jacobian(weight_fn, alpha_hat)`:
# if that matches the hand-derived analytic gradient `dw/dalpha`, the
# downstream MSM composition with `apply_model_correction()` is just matrix
# algebra on the existing primitives.
#
# This file pins the four cases where the analytic derivative is tractable:
#
#   1. static(1) on binary (HT indicator weight on the treated)
#   2. static(0) on binary (HT indicator weight on the controls)
#   3. shift(delta) on continuous Gaussian (pushforward density ratio)
#   4. ipsi(delta) on binary (Kennedy closed form)
#
# Each test block derives dw/dalpha symbolically in a comment, encodes the
# derivation in R, runs `numDeriv::jacobian(weight_fn, alpha_hat)`, and
# asserts elementwise agreement at ~1e-6. The tolerance is limited by the
# central-difference precision of `numDeriv`, not by the analytic gradient.
#
# If any of these ever drift, it means either (a) `make_weight_fn()` has a
# sign / Jacobian bug, or (b) the closure is capturing the wrong variable
# from the environment — both of which would silently corrupt the
# propensity-uncertainty correction in the sandwich. Either way, the tests
# here are what catches it.

bc_setup <- function(n = 400, seed = 701) {
  d <- simulate_binary_continuous(n = n, seed = seed)
  dt <- data.table::as.data.table(d)
  tm <- fit_treatment_model(dt, treatment = "A", confounders = ~L)
  list(data = dt, tm = tm)
}

cc_setup <- function(n = 400, seed = 702) {
  d <- simulate_continuous_continuous(n = n, seed = seed)
  dt <- data.table::as.data.table(d)
  tm <- fit_treatment_model(dt, treatment = "A", confounders = ~L)
  list(data = dt, tm = tm)
}

# ---- 1. static(1) on binary ---------------------------------------

test_that("hand-derived gradient matches numDeriv::jacobian for static(1)", {
  skip_if_not_installed("numDeriv")

  # Derivation. For a logistic propensity
  #
  #   p_i = plogis(X_i^T alpha),    d p_i / d alpha = p_i * (1 - p_i) * X_i,
  #
  # the static(1) weight is
  #
  #   w_i(alpha) = I(A_obs_i == 1) / p_i.
  #
  # Differentiating:
  #
  #   d w_i / d alpha_k
  #     = I(A = 1) * (-1 / p_i^2) * (p_i * (1 - p_i)) * X_{i,k}
  #     = -I(A = 1) * ((1 - p_i) / p_i) * X_{i,k}.
  #
  # Treated rows contribute -((1-p)/p) * X; control rows contribute 0.

  s <- bc_setup()
  tm <- s$tm
  dat <- s$data[tm$fit_rows]
  X <- tm$X_prop
  a_obs <- dat$A
  p <- stats::plogis(as.numeric(X %*% tm$alpha_hat))

  # Analytic gradient: n x p_alpha matrix.
  # Row i: -I(A_i = 1) * ((1 - p_i) / p_i) * X_i^T.
  factor_per_row <- -(a_obs == 1) * (1 - p) / p
  J_analytic <- factor_per_row * X

  wfn <- make_weight_fn(tm, s$data, static(1))
  J_numeric <- numDeriv::jacobian(wfn, x = tm$alpha_hat)

  expect_equal(dim(J_numeric), dim(J_analytic))
  # `J_analytic` inherits `assign` / `dimnames` from `tm$X_prop`
  # (which is a `model.matrix()` result); `J_numeric` is a plain
  # matrix. Drop attributes before comparing so the numeric content
  # is the only thing under test.
  expect_equal(
    J_numeric,
    J_analytic,
    tolerance = 1e-6,
    ignore_attr = TRUE
  )

  # Control rows must be exactly zero in both (the weight is
  # identically 0 on the control arm under static(1)).
  ctrl <- a_obs == 0
  expect_true(all(abs(J_analytic[ctrl, ]) < 1e-12))
  expect_true(all(abs(J_numeric[ctrl, ]) < 1e-8))
})

# ---- 2. static(0) on binary ---------------------------------------

test_that("hand-derived gradient matches numDeriv::jacobian for static(0)", {
  skip_if_not_installed("numDeriv")

  # Derivation. For static(0),
  #
  #   w_i(alpha) = I(A_obs_i == 0) / (1 - p_i).
  #
  # Differentiating through q_i = 1 - p_i:
  #
  #   d q_i / d alpha_k = -p_i * (1 - p_i) * X_{i,k},
  #   d w_i / d alpha_k = I(A = 0) * (-1 / q_i^2) * (-p_i * (1 - p_i)) * X_{i,k}
  #                     = I(A = 0) * (p_i / (1 - p_i)) * X_{i,k}.
  #
  # Treated rows contribute 0; control rows contribute (p/(1-p)) * X.

  s <- bc_setup()
  tm <- s$tm
  dat <- s$data[tm$fit_rows]
  X <- tm$X_prop
  a_obs <- dat$A
  p <- stats::plogis(as.numeric(X %*% tm$alpha_hat))

  factor_per_row <- (a_obs == 0) * p / (1 - p)
  J_analytic <- factor_per_row * X

  wfn <- make_weight_fn(tm, s$data, static(0))
  J_numeric <- numDeriv::jacobian(wfn, x = tm$alpha_hat)

  # `J_analytic` inherits `assign` / `dimnames` from `tm$X_prop`
  # (which is a `model.matrix()` result); `J_numeric` is a plain
  # matrix. Drop attributes before comparing so the numeric content
  # is the only thing under test.
  expect_equal(
    J_numeric,
    J_analytic,
    tolerance = 1e-6,
    ignore_attr = TRUE
  )

  trt <- a_obs == 1
  expect_true(all(abs(J_analytic[trt, ]) < 1e-12))
  expect_true(all(abs(J_numeric[trt, ]) < 1e-8))
})

# ---- 3. shift(delta) on continuous --------------------------------

test_that("hand-derived gradient matches numDeriv::jacobian for shift(delta)", {
  skip_if_not_installed("numDeriv")

  # Derivation. Under a Gaussian treatment model
  #
  #   A_i | L_i ~ N(mu_i = X_i^T alpha, sigma^2),
  #
  # the pushforward density ratio for shift(delta) is
  #
  #   w_i(alpha)
  #     = dnorm(A_i - delta, mu_i, sigma) / dnorm(A_i, mu_i, sigma)
  #     = exp( (1 / (2 sigma^2)) * [(A_i - mu_i)^2 - (A_i - delta - mu_i)^2] )
  #     = exp( delta * (2 (A_i - mu_i) - delta) / (2 sigma^2) ).
  #
  # Differentiating (sigma is fixed):
  #
  #   d log w_i / d alpha_k
  #     = delta / (2 sigma^2) * d/d alpha_k [2 (A_i - mu_i) - delta]
  #     = delta / (2 sigma^2) * (-2 X_{i,k})
  #     = -delta * X_{i,k} / sigma^2
  #
  # so
  #
  #   d w_i / d alpha_k = -w_i * delta * X_{i,k} / sigma^2.
  #
  # The sign reflects the usual MTP interpretation: a larger
  # propensity-regression intercept slides the fitted mean right, which
  # (for delta < 0) shrinks the weight at A_obs < mu and stretches it at
  # A_obs > mu.

  s <- cc_setup()
  tm <- s$tm
  dat <- s$data[tm$fit_rows]
  X <- tm$X_prop
  a_obs <- dat$A
  mu <- as.numeric(X %*% tm$alpha_hat)
  sigma <- tm$sigma
  delta <- -0.5

  # w_i evaluated at alpha_hat
  w <- stats::dnorm(a_obs - delta, mu, sigma) /
    stats::dnorm(a_obs, mu, sigma)

  # Row i: -w_i * delta / sigma^2 * X_i^T
  factor_per_row <- -w * delta / (sigma^2)
  J_analytic <- factor_per_row * X

  wfn <- make_weight_fn(tm, s$data, shift(delta))
  J_numeric <- numDeriv::jacobian(wfn, x = tm$alpha_hat)

  # `J_analytic` inherits `assign` / `dimnames` from `tm$X_prop`
  # (which is a `model.matrix()` result); `J_numeric` is a plain
  # matrix. Drop attributes before comparing so the numeric content
  # is the only thing under test.
  expect_equal(
    J_numeric,
    J_analytic,
    tolerance = 1e-6,
    ignore_attr = TRUE
  )

  # Sanity: a shift of 0 should give a zero jacobian (the weight is
  # identically 1, independent of alpha).
  wfn0 <- make_weight_fn(tm, s$data, shift(0))
  J0 <- numDeriv::jacobian(wfn0, x = tm$alpha_hat)
  expect_true(all(abs(J0) < 1e-9))
})

# ---- 4. ipsi(delta) on binary -------------------------------------

test_that("hand-derived gradient matches numDeriv::jacobian for ipsi(delta)", {
  skip_if_not_installed("numDeriv")

  # Derivation. Kennedy's (2019) closed form:
  #
  #   w_i(alpha) = N_i / D_i(alpha),
  #
  # where
  #
  #   N_i      = delta * A_obs_i + (1 - A_obs_i)  (constant in alpha)
  #   D_i      = delta * p_i + (1 - p_i)
  #   p_i      = plogis(X_i^T alpha),
  #   dp/dalpha_k = p_i (1 - p_i) X_{i,k}.
  #
  # Differentiating:
  #
  #   dD_i / d alpha_k = (delta - 1) * p_i (1 - p_i) * X_{i,k},
  #   d w_i / d alpha_k
  #     = N_i * (-1 / D_i^2) * dD_i / d alpha_k
  #     = - N_i * (delta - 1) * p_i (1 - p_i) * X_{i,k} / D_i^2.
  #
  # At delta = 1: D_i = 1, w_i = 1, and the (delta - 1) factor zeroes the
  # whole gradient — correct, because IPSI with delta = 1 is the natural
  # course.

  s <- bc_setup()
  tm <- s$tm
  dat <- s$data[tm$fit_rows]
  X <- tm$X_prop
  a_obs <- dat$A
  p <- stats::plogis(as.numeric(X %*% tm$alpha_hat))

  delta <- 2
  N <- delta * a_obs + (1 - a_obs)
  D <- delta * p + (1 - p)
  factor_per_row <- -N * (delta - 1) * p * (1 - p) / (D^2)
  J_analytic <- factor_per_row * X

  wfn <- make_weight_fn(tm, s$data, ipsi(delta))
  J_numeric <- numDeriv::jacobian(wfn, x = tm$alpha_hat)

  # `J_analytic` inherits `assign` / `dimnames` from `tm$X_prop`
  # (which is a `model.matrix()` result); `J_numeric` is a plain
  # matrix. Drop attributes before comparing so the numeric content
  # is the only thing under test.
  expect_equal(
    J_numeric,
    J_analytic,
    tolerance = 1e-6,
    ignore_attr = TRUE
  )

  # delta = 1 should give an identically zero jacobian.
  wfn1 <- make_weight_fn(tm, s$data, ipsi(1))
  J1 <- numDeriv::jacobian(wfn1, x = tm$alpha_hat)
  expect_true(all(abs(J1) < 1e-9))
})

# ---- 5. NULL (natural course) -------------------------------------

test_that("natural-course weight closure has identically zero gradient", {
  skip_if_not_installed("numDeriv")

  # The NULL intervention corresponds to w_i(alpha) = 1 for all i and
  # all alpha. Its gradient is identically zero — which means the
  # propensity-uncertainty correction for a natural-course contrast
  # vanishes, as it should (there's no intervention to invert).
  s <- bc_setup()
  wfn <- make_weight_fn(s$tm, s$data, NULL)
  J <- numDeriv::jacobian(wfn, x = s$tm$alpha_hat)
  expect_true(all(abs(J) < 1e-12))
})
