# Unit tests for R/ipw_weights.R (Phase 4 foundation layer).
#
# The weight builder has three branches, all of which must be pinned:
#
#   1. Horvitz-Thompson indicator:
#      static() / dynamic() on discrete treatments.
#      w_i = I(A_obs_i == target_i) / f(A_obs_i | L_i).
#
#   2. Smooth density ratio:
#      shift() / scale_by() / threshold() on continuous treatments.
#      w_i = f(d(A_i, L_i) | L_i) / f(A_i | L_i), Gaussian-based.
#
#   3. IPSI closed form:
#      w_i = (delta * A_i + (1 - A_i)) / (delta * p_i + (1 - p_i)).
#
# Plus the weight_fn closure must agree with compute_density_ratio_weights()
# at `alpha = alpha_hat` (critical because the variance engine's Branch B
# uses the closure via `numDeriv::jacobian`).
#
# Finally, a truth test: the Hájek mean built from static(1) / static(0)
# weights on `simulate_binary_continuous()` must approximate the analytic
# E[Y^1] = 5 and E[Y^0] = 2 within Monte Carlo tolerance. This is the
# smallest end-to-end sanity check that "these are actually IPW weights."

# ---- helpers --------------------------------------------------------

bc_tm <- function(n = 2000, seed = 11) {
  d <- simulate_binary_continuous(n = n, seed = seed)
  dt <- data.table::as.data.table(d)
  tm <- fit_treatment_model(dt, treatment = "A", confounders = ~L)
  list(data = dt, tm = tm, df = d)
}

cc_tm <- function(n = 2000, seed = 12) {
  d <- simulate_continuous_continuous(n = n, seed = seed)
  dt <- data.table::as.data.table(d)
  tm <- fit_treatment_model(dt, treatment = "A", confounders = ~L)
  list(data = dt, tm = tm, df = d)
}

hajek <- function(w, y) sum(w * y) / sum(w)

# ---- compute_density_ratio_weights() -------------------------------

test_that("NULL intervention returns ones", {
  s <- bc_tm(n = 200)
  w <- compute_density_ratio_weights(s$tm, s$data, NULL)
  expect_equal(w, rep(1, sum(s$tm$fit_rows)))
})

test_that("static(1) on binary returns I(A=1)/p (HT form)", {
  s <- bc_tm(n = 500)
  w <- compute_density_ratio_weights(s$tm, s$data, static(1))

  fit_data <- s$data[s$tm$fit_rows]
  p <- stats::predict(s$tm$model, newdata = fit_data, type = "response")
  a_obs <- fit_data$A
  expected <- as.numeric(a_obs == 1) / ifelse(a_obs == 1, p, 1 - p)
  expect_equal(w, unname(expected), tolerance = 1e-12)

  # Controls contribute zero weight.
  expect_true(all(w[a_obs == 0] == 0))
  # Treated contribute 1/p, which is finite and positive.
  expect_true(all(is.finite(w[a_obs == 1])) && all(w[a_obs == 1] > 0))
})

test_that("static(0) on binary returns I(A=0)/(1-p) (HT form)", {
  s <- bc_tm(n = 500)
  w <- compute_density_ratio_weights(s$tm, s$data, static(0))

  fit_data <- s$data[s$tm$fit_rows]
  p <- stats::predict(s$tm$model, newdata = fit_data, type = "response")
  a_obs <- fit_data$A
  expected <- as.numeric(a_obs == 0) / ifelse(a_obs == 1, p, 1 - p)
  expect_equal(w, unname(expected), tolerance = 1e-12)

  expect_true(all(w[a_obs == 1] == 0))
})

test_that("static(v) on continuous treatment aborts with a helpful error", {
  s <- cc_tm(n = 200)
  expect_error(
    compute_density_ratio_weights(s$tm, s$data, static(1)),
    "degenerate"
  )
})

test_that("shift(delta) on continuous returns the pushforward density ratio", {
  s <- cc_tm(n = 500)
  delta <- -0.5
  w <- compute_density_ratio_weights(s$tm, s$data, shift(delta))

  fit_data <- s$data[s$tm$fit_rows]
  mu <- stats::predict(s$tm$model, newdata = fit_data, type = "response")
  sigma <- s$tm$sigma
  # Pushforward of f under d(a) = a + delta: f_d(y|l) = f(y - delta|l).
  # Weight at A_obs: f(A_obs - delta|l) / f(A_obs|l).
  # NOTE: this is NOT `f(A_obs + delta|l) / f(A_obs|l)` — that naive
  # form (the one our first-pass code carried) gives E[Y^{shift(-delta)}]
  # instead of E[Y^{shift(delta)}]. Verified analytically during the
  # Phase 4 foundation-chunk review.
  expected <- stats::dnorm(fit_data$A - delta, mu, sigma) /
    stats::dnorm(fit_data$A, mu, sigma)
  expect_equal(w, unname(expected), tolerance = 1e-12)

  # A zero-shift intervention should return ones exactly.
  w0 <- compute_density_ratio_weights(s$tm, s$data, shift(0))
  expect_equal(w0, rep(1, length(w0)), tolerance = 1e-12)
})

test_that("scale_by(factor) on continuous returns the pushforward density ratio with Jacobian", {
  s <- cc_tm(n = 500)
  fct <- 0.75
  w <- compute_density_ratio_weights(s$tm, s$data, scale_by(fct))

  fit_data <- s$data[s$tm$fit_rows]
  mu <- stats::predict(s$tm$model, newdata = fit_data, type = "response")
  sigma <- s$tm$sigma
  # Pushforward of f under d(a) = c*a, d^-1(y) = y/c, |Jac|=1/|c|.
  # f_d(y|l) = f(y/c|l) / |c|, so weight = f(A/c|l) / (|c| * f(A|l)).
  expected <- stats::dnorm(fit_data$A / fct, mu, sigma) /
    (abs(fct) * stats::dnorm(fit_data$A, mu, sigma))
  expect_equal(w, unname(expected), tolerance = 1e-12)

  # scale_by(1) is the identity → weights all 1.
  w1 <- compute_density_ratio_weights(s$tm, s$data, scale_by(1))
  expect_equal(w1, rep(1, length(w1)), tolerance = 1e-12)
})

test_that("threshold(lo, hi) under IPW is rejected with a gcomp pointer", {
  # The pushforward of a continuous density under a boundary clamp
  # has point masses at `lo` and `hi`, so the density ratio w.r.t.
  # the fitted Gaussian is not well-defined. The IPW engine rejects
  # the combination and points users to `estimator = 'gcomp'`.
  s <- cc_tm(n = 200)
  expect_error(
    compute_density_ratio_weights(s$tm, s$data, threshold(-1, 1)),
    "not supported"
  )
  expect_error(
    make_weight_fn(s$tm, s$data, threshold(-1, 1)),
    "not supported"
  )
})

test_that("scale_by(0) is rejected as a degenerate MTP", {
  s <- cc_tm(n = 200)
  expect_error(
    compute_density_ratio_weights(s$tm, s$data, scale_by(0)),
    "collapses"
  )
})

test_that("shift() on binary treatment aborts with a compat error", {
  s <- bc_tm(n = 200)
  expect_error(
    compute_density_ratio_weights(s$tm, s$data, shift(0.5)),
    "numeric continuous"
  )
})

test_that("dynamic(rule) on binary returns the HT indicator weight", {
  s <- bc_tm(n = 500)
  # Rule: treat iff L > 0.
  iv <- dynamic(function(d, a) as.numeric(d$L > 0))
  w <- compute_density_ratio_weights(s$tm, s$data, iv)

  fit_data <- s$data[s$tm$fit_rows]
  p <- stats::predict(s$tm$model, newdata = fit_data, type = "response")
  a_obs <- fit_data$A
  rule_vals <- as.numeric(fit_data$L > 0)
  expected <- as.numeric(a_obs == rule_vals) /
    ifelse(a_obs == 1, p, 1 - p)
  expect_equal(w, unname(expected), tolerance = 1e-12)

  # Rows where the observed treatment disagrees with the rule contribute zero.
  disagree <- a_obs != rule_vals
  expect_true(all(w[disagree] == 0))
})

test_that("dynamic() on continuous treatment aborts", {
  s <- cc_tm(n = 200)
  iv <- dynamic(function(d, a) d$L)
  expect_error(
    compute_density_ratio_weights(s$tm, s$data, iv),
    "degenerate"
  )
})

test_that("ipsi(delta) on binary returns Kennedy closed-form", {
  s <- bc_tm(n = 500)
  delta <- 2
  w <- compute_density_ratio_weights(s$tm, s$data, ipsi(delta))

  fit_data <- s$data[s$tm$fit_rows]
  p <- unname(stats::predict(s$tm$model, newdata = fit_data, type = "response"))
  a_obs <- fit_data$A
  expected <- (delta * a_obs + (1 - a_obs)) / (delta * p + (1 - p))
  expect_equal(w, expected, tolerance = 1e-12)

  # Sanity: delta = 1 recovers the natural course (ratio is 1 for
  # everyone because the propensity is unchanged).
  w1 <- compute_density_ratio_weights(s$tm, s$data, ipsi(1))
  expect_equal(w1, rep(1, length(w1)), tolerance = 1e-12)
})

test_that("ipsi() on continuous treatment aborts", {
  s <- cc_tm(n = 200)
  expect_error(
    compute_density_ratio_weights(s$tm, s$data, ipsi(2)),
    "binary"
  )
})

# ---- truth test -----------------------------------------------------

test_that("Hajek means from static(1) / static(0) recover analytic E[Y^a]", {
  # simulate_binary_continuous: Y = 2 + 3*A + 1.5*L + noise
  # => E[Y^1] = 5, E[Y^0] = 2, ATE = 3.
  # Under a CORRECTLY SPECIFIED logistic propensity (A ~ L, which is the
  # true DGP), the Hajek-IPW means converge to the truth, so a moderately
  # large n should hit the target to well within Monte Carlo tolerance.
  s <- bc_tm(n = 5000, seed = 101)
  fit_data <- s$data[s$tm$fit_rows]
  y <- fit_data$Y

  w1 <- compute_density_ratio_weights(s$tm, s$data, static(1))
  w0 <- compute_density_ratio_weights(s$tm, s$data, static(0))

  mu1 <- hajek(w1, y)
  mu0 <- hajek(w0, y)

  expect_equal(mu1, 5, tolerance = 0.15)
  expect_equal(mu0, 2, tolerance = 0.15)
  expect_equal(mu1 - mu0, 3, tolerance = 0.2)
})

test_that("Hajek mean from shift(delta) on continuous matches analytic target", {
  # simulate_continuous_continuous: Y = 1 + 2*A + L + noise
  # counterfactual under shift(delta): A^d = A + delta, so
  #   Y^d = 1 + 2*(A + delta) + L + noise
  #   E[Y^{shift(delta)}] = 1 + 2*(E[A] + delta) + E[L]
  #                      = 1 + 2*(1 + delta) + 0
  #                      = 3 + 2 * delta     (since E[A] = 1, E[L] = 0)
  # shift(-1) => 3 - 2 = 1
  # shift(0)  => 3
  # Contrast shift(-1) - shift(0) = -2, matching the helper-dgp comment.
  #
  # The correct IPW weight for shift(delta) is the pushforward density
  # ratio: w_i = f(A_obs_i - delta | L_i) / f(A_obs_i | L_i). Getting
  # the sign wrong here was the bug that this test first uncovered.
  s <- cc_tm(n = 5000, seed = 102)
  fit_data <- s$data[s$tm$fit_rows]
  y <- fit_data$Y

  w_shift <- compute_density_ratio_weights(s$tm, s$data, shift(-1))
  w_none <- compute_density_ratio_weights(s$tm, s$data, shift(0))

  mu_shift <- hajek(w_shift, y)
  mu_none <- hajek(w_none, y)

  expect_equal(mu_shift, 1, tolerance = 0.25)
  expect_equal(mu_none, 3, tolerance = 0.15)
  expect_equal(mu_shift - mu_none, -2, tolerance = 0.3)
})

test_that("Hajek mean from scale_by on continuous matches analytic target", {
  # Y = 1 + 2*A + L, counterfactual under scale_by(c): A^d = c*A.
  # E[Y^{scale(c)}] = 1 + 2*c*E[A] + 0 = 1 + 2*c.
  # For c = 0.5 => 2. For c = 1 => 3 (natural course).
  s <- cc_tm(n = 5000, seed = 103)
  fit_data <- s$data[s$tm$fit_rows]
  y <- fit_data$Y

  w_half <- compute_density_ratio_weights(s$tm, s$data, scale_by(0.5))
  mu_half <- hajek(w_half, y)
  expect_equal(mu_half, 2, tolerance = 0.4)
})

# ---- make_weight_fn() closure consistency ---------------------------

test_that("make_weight_fn() at alpha_hat matches compute_density_ratio_weights()", {
  # Pin each intervention × family combo where the closure is wired.
  # The variance engine calls `weight_fn(alpha_hat)` and expects the
  # result to be indistinguishable from the fit-time weight vector —
  # any drift here silently corrupts numDeriv::jacobian's base point.
  s_b <- bc_tm(n = 300, seed = 21)
  s_c <- cc_tm(n = 300, seed = 22)

  combos <- list(
    list(s = s_b, iv = NULL),
    list(s = s_b, iv = static(1)),
    list(s = s_b, iv = static(0)),
    list(
      s = s_b,
      iv = dynamic(function(d, a) as.numeric(d$L > 0))
    ),
    list(s = s_b, iv = ipsi(1.8)),
    list(s = s_c, iv = NULL),
    list(s = s_c, iv = shift(-0.5)),
    list(s = s_c, iv = scale_by(0.8))
  )

  for (combo in combos) {
    w_direct <- compute_density_ratio_weights(
      combo$s$tm,
      combo$s$data,
      combo$iv
    )
    wfn <- make_weight_fn(combo$s$tm, combo$s$data, combo$iv)
    w_closure <- wfn(combo$s$tm$alpha_hat)
    expect_equal(w_closure, w_direct, tolerance = 1e-10)
  }
})

test_that("make_weight_fn(NULL) closure returns ones regardless of alpha", {
  s <- bc_tm(n = 200)
  wfn <- make_weight_fn(s$tm, s$data, NULL)
  expect_equal(wfn(s$tm$alpha_hat), rep(1, sum(s$tm$fit_rows)))
  expect_equal(wfn(rep(0, length(s$tm$alpha_hat))), rep(1, sum(s$tm$fit_rows)))
  expect_equal(
    wfn(rep(999, length(s$tm$alpha_hat))),
    rep(1, sum(s$tm$fit_rows))
  )
})

test_that("numDeriv::jacobian() on make_weight_fn() returns a non-trivial matrix for static(1)", {
  skip_if_not_installed("numDeriv")
  s <- bc_tm(n = 200, seed = 31)
  wfn <- make_weight_fn(s$tm, s$data, static(1))

  J <- numDeriv::jacobian(wfn, x = s$tm$alpha_hat)

  expect_equal(dim(J), c(sum(s$tm$fit_rows), length(s$tm$alpha_hat)))
  # For static(1) the weight is either 0 (controls) or 1/p (treated), so
  # d/d_alpha w != 0 on the treated rows — the jacobian cannot be all-zero.
  expect_true(any(abs(J) > 1e-6))
  # Control rows (a_obs == 0) have weight identically 0 under static(1),
  # so their row of J must be all zero.
  fit_data <- s$data[s$tm$fit_rows]
  zero_rows <- fit_data$A == 0
  expect_true(all(abs(J[zero_rows, ]) < 1e-10))
})

test_that("numDeriv::jacobian() on make_weight_fn() returns a non-trivial matrix for shift(-0.5)", {
  skip_if_not_installed("numDeriv")
  s <- cc_tm(n = 200, seed = 32)
  wfn <- make_weight_fn(s$tm, s$data, shift(-0.5))

  J <- numDeriv::jacobian(wfn, x = s$tm$alpha_hat)
  expect_equal(dim(J), c(sum(s$tm$fit_rows), length(s$tm$alpha_hat)))
  expect_true(any(abs(J) > 1e-6))
})

test_that("numDeriv::jacobian() on ipsi(2) closure is non-zero and finite", {
  skip_if_not_installed("numDeriv")
  s <- bc_tm(n = 200, seed = 33)
  wfn <- make_weight_fn(s$tm, s$data, ipsi(2))

  J <- numDeriv::jacobian(wfn, x = s$tm$alpha_hat)
  expect_true(all(is.finite(J)))
  expect_true(any(abs(J) > 1e-6))
})

# ---- positivity guard -----------------------------------------------

test_that("check_density_positivity() aborts on zero / non-finite values", {
  # Direct unit test on the guard helper. Going through
  # `compute_density_ratio_weights()` with a perfectly-separable
  # logistic fit does NOT reliably trigger an exact zero —
  # `glm.fit` converges to finite coefficients and predictions
  # saturate at ~0.9999 rather than exactly 1 — so we exercise the
  # guard directly on synthetic density vectors.
  expect_error(check_density_positivity(c(0.5, 0, 0.7), "test"), "positivity")
  expect_error(
    check_density_positivity(c(0.5, NaN, 0.7), "test"),
    "positivity"
  )
  expect_error(
    check_density_positivity(c(0.5, Inf, 0.7), "test"),
    "positivity"
  )
  expect_error(
    check_density_positivity(c(0.5, -0.1, 0.7), "test"),
    "positivity"
  )
  expect_silent(check_density_positivity(c(0.5, 0.6, 0.7), "test"))
})
