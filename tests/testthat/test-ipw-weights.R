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


# Fifth-round critical review, S2: warn when the intervened density is
# near-zero for most observations, indicating the intervention pushes
# treatment outside the fitted distribution's support.
# Repro script: /tmp/causatr_repro_s2_f_int_positivity.R
test_that("warn_intervened_density_near_zero fires when > 80% near-zero", {
  # All near-zero: should warn
  f_bad <- rep(1e-15, 100)
  expect_warning(
    warn_intervened_density_near_zero(f_bad, "test"),
    class = "causatr_near_zero_intervened_density"
  )
})

test_that("warn_intervened_density_near_zero is silent for healthy densities", {
  f_good <- rep(0.5, 100)
  expect_silent(warn_intervened_density_near_zero(f_good, "test"))
})

test_that("shift with extreme delta warns about near-zero intervened density", {
  set.seed(42)
  n <- 200
  L <- rnorm(n)
  A <- 5 + 0.5 * L + rnorm(n, sd = 0.5)
  Y <- 2 + 0.5 * A + L + rnorm(n)
  d <- data.table::data.table(Y = Y, A = A, L = L)

  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw"
  )

  # shift(-6) pushes treatment well outside the support
  expect_warning(
    contrast(
      fit,
      interventions = list(shifted = shift(-6), baseline = shift(0)),
      type = "difference"
    ),
    class = "causatr_near_zero_intervened_density"
  )
})

# ---- make_weight_fn() defensive guards ------------------------------------
#
# The aborts in make_weight_fn() are nominally unreachable through the
# public API because `check_intervention_family_compat()` rejects
# every (intervention, family) combination they would catch. The tests
# below construct evil inputs (wrong-class objects, hand-rolled
# intervention objects with skipped validation) to exercise the guards
# directly. The motivation is the user's standing rule: every function
# in R/ gets a test, including defensive paths -- a regression that
# silently bypasses CIFC would otherwise produce nonsense weights.

test_that("make_weight_fn() rejects a non-treatment_model `treatment_model`", {
  expect_error(
    make_weight_fn(list(family = "bernoulli"), data.frame(), static(1)),
    "causatr_treatment_model"
  )
})

test_that("make_weight_fn() aborts on scale_by(0) for a Gaussian treatment", {
  s <- cc_tm(n = 200)
  # `compute_density_ratio_weights()` has its own scale_by(0) guard at
  # line ~250; this test exercises the duplicate guard in
  # `make_weight_fn()` itself, which protects the closure path used by
  # `numDeriv::jacobian()` inside the variance engine.
  expect_error(
    make_weight_fn(s$tm, s$data, scale_by(0)),
    "collapses"
  )
})

test_that("make_weight_fn() aborts on scale_by(0) for a count treatment", {
  set.seed(31)
  d <- data.table::as.data.table(data.frame(
    A = rpois(200, lambda = 3),
    L = rnorm(200)
  ))
  tm <- fit_treatment_model(
    d,
    treatment = "A",
    confounders = ~L,
    propensity_family = "poisson"
  )
  # The `scale_by(0)` guard for count treatments fires in
  # `check_intervention_family_compat()` first (the integer-preservation
  # check trips on `0 / 0 = NaN`), so to exercise the Poisson branch's
  # internal `if (fct == 0)` guard at line ~588 we sidestep CIFC by
  # mocking it to a no-op. This mirrors the gaussian test above and
  # protects the closure path against a future regression that bypasses
  # the upstream check.
  testthat::local_mocked_bindings(
    check_intervention_family_compat = function(...) invisible(NULL)
  )
  expect_error(
    make_weight_fn(tm, d, scale_by(0)),
    "collapses"
  )
})

test_that("make_weight_fn() defensive: gaussian + threshold internal-error guard", {
  s <- cc_tm(n = 200)
  # `check_intervention_family_compat()` rejects threshold/gaussian
  # upstream; mock it to a no-op to reach the internal-error abort
  # inside `make_weight_fn()`'s gaussian branch.
  testthat::local_mocked_bindings(
    check_intervention_family_compat = function(...) invisible(NULL)
  )
  expect_error(
    make_weight_fn(s$tm, s$data, threshold(-1, 1)),
    "Internal error.*gaussian branch"
  )
})

test_that("make_weight_fn() defensive: count + threshold internal-error guard", {
  set.seed(32)
  d <- data.table::as.data.table(data.frame(
    A = rpois(200, lambda = 3),
    L = rnorm(200)
  ))
  tm <- fit_treatment_model(
    d,
    treatment = "A",
    confounders = ~L,
    propensity_family = "poisson"
  )
  testthat::local_mocked_bindings(
    check_intervention_family_compat = function(...) invisible(NULL)
  )
  expect_error(
    make_weight_fn(tm, d, threshold(0, 5)),
    "Internal error.*count branch"
  )
})

test_that("make_weight_fn() defensive: unknown family final-fallback abort", {
  # Forge a treatment_model with an unrecognised family tag. The final
  # `rlang::abort()` at the bottom of make_weight_fn() catches this
  # without going through any of the family-specific branches.
  s <- bc_tm(n = 100)
  evil_tm <- s$tm
  evil_tm$family <- "frobnitz"
  testthat::local_mocked_bindings(
    check_intervention_family_compat = function(...) invisible(NULL)
  )
  expect_error(
    make_weight_fn(evil_tm, s$data, static(1)),
    "does not handle"
  )
})

# ---- apply_intervention_to_values() ---------------------------------------
#
# Used internally by the HT branch of make_weight_fn() for static and
# dynamic on bernoulli/categorical. The shift / scale / threshold / ipsi
# branches are reached only via direct calls (production code skips
# this helper for those). Tests below cover the full switch and every
# defensive guard inside the dynamic branch.

test_that("apply_intervention_to_values() static returns rep(value)", {
  out <- apply_intervention_to_values(
    static(1),
    data = data.frame(),
    a_obs = c(0, 1, 0, 1)
  )
  expect_equal(out, c(1, 1, 1, 1))
})

test_that("apply_intervention_to_values() shift returns A + delta", {
  out <- apply_intervention_to_values(
    shift(0.5),
    data = data.frame(),
    a_obs = c(0, 1, 2)
  )
  expect_equal(out, c(0.5, 1.5, 2.5))
})

test_that("apply_intervention_to_values() scale returns A * factor", {
  out <- apply_intervention_to_values(
    scale_by(2),
    data = data.frame(),
    a_obs = c(0, 1, 2)
  )
  expect_equal(out, c(0, 2, 4))
})

test_that("apply_intervention_to_values() threshold clamps to [lo, hi]", {
  out <- apply_intervention_to_values(
    threshold(0, 1),
    data = data.frame(),
    a_obs = c(-2, 0.5, 3)
  )
  expect_equal(out, c(0, 0.5, 1))
})

test_that("apply_intervention_to_values() dynamic length-mismatch aborts", {
  bad_rule <- dynamic(function(d, a) c(0, 1)) # length 2, should be 3
  expect_error(
    apply_intervention_to_values(
      bad_rule,
      data = data.frame(L = c(1, 2, 3)),
      a_obs = c(0, 1, 0)
    ),
    "length 3"
  )
})

test_that("apply_intervention_to_values() dynamic non-numeric for numeric A aborts", {
  bad_rule <- dynamic(function(d, a) rep("x", length(a)))
  expect_error(
    apply_intervention_to_values(
      bad_rule,
      data = data.frame(L = c(1, 2)),
      a_obs = c(0.1, 0.2)
    ),
    "non-numeric"
  )
})

test_that("apply_intervention_to_values() dynamic char-with-unknown-level for factor A aborts", {
  a_factor <- factor(c("a", "b", "a"), levels = c("a", "b"))
  bad_rule <- dynamic(function(d, a) c("a", "z", "a"))
  expect_error(
    apply_intervention_to_values(
      bad_rule,
      data = data.frame(L = c(1, 2, 3)),
      a_obs = a_factor
    ),
    "level\\(s\\) not in"
  )
})

test_that("apply_intervention_to_values() dynamic char-with-known-level coerces back to factor", {
  a_factor <- factor(c("a", "b", "a"), levels = c("a", "b"))
  good_rule <- dynamic(function(d, a) c("b", "a", "b"))
  out <- apply_intervention_to_values(
    good_rule,
    data = data.frame(L = c(1, 2, 3)),
    a_obs = a_factor
  )
  expect_s3_class(out, "factor")
  expect_equal(levels(out), c("a", "b"))
  expect_equal(as.character(out), c("b", "a", "b"))
})

test_that("apply_intervention_to_values() dynamic factor with mismatched levels aborts", {
  a_factor <- factor(c("a", "b", "a"), levels = c("a", "b"))
  wrong_factor <- factor(c("a", "b", "a"), levels = c("a", "b", "c"))
  bad_rule <- dynamic(function(d, a) wrong_factor)
  expect_error(
    apply_intervention_to_values(
      bad_rule,
      data = data.frame(L = c(1, 2, 3)),
      a_obs = a_factor
    ),
    "mismatched levels"
  )
})

test_that("apply_intervention_to_values() dynamic non-factor non-character for factor A aborts", {
  a_factor <- factor(c("a", "b", "a"), levels = c("a", "b"))
  bad_rule <- dynamic(function(d, a) c(1.5, 2.5, 1.5))
  expect_error(
    apply_intervention_to_values(
      bad_rule,
      data = data.frame(L = c(1, 2, 3)),
      a_obs = a_factor
    ),
    "non-factor, non-character"
  )
})

test_that("apply_intervention_to_values() rejects ipsi (caller-bug guard)", {
  # IPSI never reaches this helper through the public API -- it's a
  # closed-form weight branch that bypasses the HT path entirely. The
  # guard exists so a future refactor that accidentally routes IPSI
  # through here aborts loudly instead of producing nonsense.
  expect_error(
    apply_intervention_to_values(
      ipsi(2),
      data = data.frame(L = c(1, 2)),
      a_obs = c(0, 1)
    ),
    "should not be called with an IPSI"
  )
})

# ---- check_intervention_family_compat() rejection branches ---------------
#
# CIFC is the upstream gate that the runtime weight builders rely on.
# Most happy-path branches are exercised transitively by the rest of
# this file; the tests below pin the rarer rejection messages so a
# regression that flips one of them surfaces directly.

test_that("check_intervention_family_compat() rejects threshold on a binary treatment", {
  s <- bc_tm(n = 100)
  # Binary fam = "bernoulli"; falls into the `iv_type == "threshold"
  # && fam != "gaussian" && !is_count` branch.
  expect_error(
    check_intervention_family_compat(threshold(0, 1), s$tm),
    "require a numeric continuous treatment"
  )
})

test_that("check_intervention_family_compat() rejects threshold on a categorical treatment", {
  set.seed(61)
  d <- data.table::as.data.table(data.frame(
    A = factor(sample(letters[1:3], 200, replace = TRUE)),
    L = rnorm(200)
  ))
  tm <- fit_treatment_model(
    d,
    treatment = "A",
    confounders = ~L,
    model_fn = nnet::multinom
  )
  expect_error(
    check_intervention_family_compat(threshold(0, 1), tm),
    "require a numeric continuous treatment"
  )
})

test_that("check_intervention_family_compat() falls back to model.response for scale on count", {
  # When the treatment_model$model lacks `$data` (some fitters drop it),
  # CIFC reconstructs the observed treatment via
  # `model.response(model.frame(.))`. Fit a count model, blank out
  # `$data`, and check that `scale_by()` validation still works on the
  # reconstructed vector. We use a non-integer-preserving factor so the
  # fallback path completes and aborts on the integer-support guard.
  set.seed(62)
  d <- data.table::as.data.table(data.frame(
    A = rpois(150, lambda = 4),
    L = rnorm(150)
  ))
  tm <- fit_treatment_model(
    d,
    treatment = "A",
    confounders = ~L,
    propensity_family = "poisson"
  )
  tm$model$data <- NULL # force the model.response fallback
  # `scale_by(2.5)` gives `a / 2.5`, which is rarely integer-valued
  # for integer `a` -- triggers the "not integer inverse" abort.
  expect_error(
    check_intervention_family_compat(scale_by(2.5), tm),
    "does not produce integer inverse"
  )
})

# ---- compute_density_ratio_weights() defensive guards --------------------

test_that("compute_density_ratio_weights() rejects a non-treatment_model `treatment_model`", {
  expect_error(
    compute_density_ratio_weights(
      list(family = "bernoulli"),
      data.frame(),
      static(1)
    ),
    "causatr_treatment_model"
  )
})

test_that("compute_density_ratio_weights() defensive: unreachable internal-error fallback", {
  # CIFC rejects threshold on gaussian upstream, so the final
  # internal-error abort in compute_density_ratio_weights() is only
  # reachable if a future refactor accidentally bypasses CIFC. Mock
  # CIFC to a no-op to drive that fallback path.
  s <- cc_tm(n = 100)
  testthat::local_mocked_bindings(
    check_intervention_family_compat = function(...) invisible(NULL)
  )
  expect_error(
    compute_density_ratio_weights(s$tm, s$data, threshold(-1, 1)),
    "Internal error.*has no branch"
  )
})

# ---- ht_bayes_numerator() per-estimand branches --------------------------
#
# `ht_bayes_numerator()` returns the f*_i = f(A* | L_i) factor used by
# the HT density-ratio weight under each estimand. ATE returns 1 (no
# propensity-score uncertainty), ATT returns p(L), ATC returns 1 - p(L).
# Tests below pin each branch + the two internal-error guards.

test_that("ht_bayes_numerator(ATE) returns the constant 1", {
  s <- bc_tm(n = 100)
  out <- ht_bayes_numerator("ATE", s$tm, s$data, "bernoulli")
  expect_identical(out, 1)
})

test_that("ht_bayes_numerator(ATT) returns predicted p(L)", {
  s <- bc_tm(n = 200)
  out <- ht_bayes_numerator("ATT", s$tm, s$data, "bernoulli")
  expected <- unname(stats::predict(
    s$tm$model,
    newdata = s$data,
    type = "response"
  ))
  expect_equal(out, expected)
})

test_that("ht_bayes_numerator(ATC) returns 1 - p(L)", {
  s <- bc_tm(n = 200)
  out <- ht_bayes_numerator("ATC", s$tm, s$data, "bernoulli")
  expected <- 1 -
    unname(stats::predict(
      s$tm$model,
      newdata = s$data,
      type = "response"
    ))
  expect_equal(out, expected)
})

test_that("ht_bayes_numerator() aborts when ATT/ATC requested for non-Bernoulli", {
  s <- cc_tm(n = 100)
  expect_error(
    ht_bayes_numerator("ATT", s$tm, s$data, "gaussian"),
    "non-Bernoulli treatment family"
  )
})

test_that("ht_bayes_numerator() aborts on unknown estimand", {
  s <- bc_tm(n = 100)
  expect_error(
    ht_bayes_numerator("FROBNITZ", s$tm, s$data, "bernoulli"),
    "unknown estimand"
  )
})

test_that("apply_intervention_to_values() rejects an unknown intervention type", {
  evil_iv <- structure(
    list(type = "frobnitz"),
    class = c("causatr_intervention", "list")
  )
  expect_error(
    apply_intervention_to_values(
      evil_iv,
      data = data.frame(),
      a_obs = c(0, 1)
    ),
    "Unknown intervention type"
  )
})
