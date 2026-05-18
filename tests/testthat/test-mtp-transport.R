# Tests for MTP (modified treatment policy) + transportability (Phase 17l).
#
# MTP interventions (shift, scale_by, threshold, dynamic) combined with
# transportability require MC marginalization: target rows (S=0) have no
# observed treatment, so we integrate E_{A|L,S=1}[m(d(A,L), L)] over
# the study treatment distribution.
#
# DGP: simulate_mtp_transport() in helper-dgp.R
#   L ~ N(0,1); P(S=1|L) = expit(-0.5 + L)
#   A|L,S=1 ~ N(0.5 + 0.3L, 1); Y|A,L,S=1 ~ N(2 + 3A + 1.5L + AL, 1)
#
# The A:L interaction must be included in confounders for gcomp to
# recover the heterogeneous treatment effect.
#
# Truth for shift(delta) vs natural:
#   E[Y^{A+delta}] - E[Y^A] = 3*delta + delta*E[L|target]

# Helper: extract contrast estimate from a causatr_result
get_contrast_est <- function(res) {
  res$contrasts$estimate[1L]
}

# -- Gcomp + shift + transportability: truth recovery --------------------

test_that("gcomp + shift + transport: recovers analytical truth", {
  d <- simulate_mtp_transport(n = 8000, seed = 101)
  delta <- 1

  fit <- causat(d, "Y", "A", ~ L + A:L, target = "S")

  res <- contrast(
    fit,
    interventions = list(shifted = shift(delta), natural = NULL),
    reference = "natural",
    type = "difference",
    ci_method = "bootstrap",
    n_boot = 200
  )

  truth <- 3 * delta + delta * mean(d$L[d$S == 0])
  expect_equal(get_contrast_est(res), truth, tolerance = 0.15)
})


# -- Gcomp + scale_by + transportability: truth recovery ------------------

test_that("gcomp + scale_by + transport: recovers truth", {
  d <- simulate_mtp_transport(n = 8000, seed = 102)
  # scale_by(2): A_d = 2*A. Y(a,L) = 2 + 3a + 1.5L + aL
  # Y(2a,L) - Y(a,L) = 3a + aL = a(3 + L)
  # E[Y^{2A} - Y^A | target] = E_L[E_A[A(3+L) | L, S=1] | target]
  #   = E_L[(0.5+0.3L)(3+L) | target]

  fit <- causat(d, "Y", "A", ~ L + A:L, target = "S")

  res <- contrast(
    fit,
    interventions = list(scaled = scale_by(2), natural = NULL),
    reference = "natural",
    type = "difference",
    ci_method = "bootstrap",
    n_boot = 200
  )

  L_target <- d$L[d$S == 0]
  mu_A_given_L <- 0.5 + 0.3 * L_target
  truth <- mean(mu_A_given_L * (3 + L_target))
  expect_equal(get_contrast_est(res), truth, tolerance = 0.20)
})


# -- Gcomp + threshold + transportability: smoke test ---------------------

test_that("gcomp + threshold + transport: runs without error", {
  d <- simulate_mtp_transport(n = 4000, seed = 103)

  fit <- causat(d, "Y", "A", ~ L + A:L, target = "S")
  res <- contrast(
    fit,
    interventions = list(clamped = threshold(-1, 1), natural = NULL),
    reference = "natural",
    ci_method = "bootstrap",
    n_boot = 50
  )

  expect_s3_class(res, "causatr_result")
  expect_true(is.finite(get_contrast_est(res)))
})


# -- Gcomp + dynamic + transportability: truth recovery -------------------

test_that("gcomp + dynamic + transport: recovers truth", {
  # dynamic() that doubles treatment: d(data, trt) = 2 * trt
  # Same truth as scale_by(2).
  d <- simulate_mtp_transport(n = 8000, seed = 104)

  fit <- causat(d, "Y", "A", ~ L + A:L, target = "S")

  double_trt <- dynamic(function(data, trt) 2 * trt)
  res <- contrast(
    fit,
    interventions = list(doubled = double_trt, natural = NULL),
    reference = "natural",
    type = "difference",
    ci_method = "bootstrap",
    n_boot = 200
  )

  L_target <- d$L[d$S == 0]
  mu_A_given_L <- 0.5 + 0.3 * L_target
  truth <- mean(mu_A_given_L * (3 + L_target))
  expect_equal(get_contrast_est(res), truth, tolerance = 0.20)
})


# -- Gcomp + shift + generalizability: truth recovery ---------------------

test_that("gcomp + shift + generalizability: recovers truth", {
  d <- simulate_mtp_transport(n = 8000, seed = 105)
  delta <- 1

  fit <- causat(d, "Y", "A", ~ L + A:L, target = "S", target_subset = "all")

  res <- contrast(
    fit,
    interventions = list(shifted = shift(delta), natural = NULL),
    reference = "natural",
    type = "difference",
    ci_method = "bootstrap",
    n_boot = 200
  )

  # Under generalizability, target = all rows
  truth <- 3 * delta + delta * mean(d$L)
  expect_equal(get_contrast_est(res), truth, tolerance = 0.15)
})


# -- IPW + shift + transportability: truth recovery -----------------------

test_that("IPW + shift + transport: recovers truth", {
  d <- simulate_mtp_transport(n = 8000, seed = 106)
  delta <- 1

  fit <- causat(
    d,
    "Y",
    "A",
    ~L,
    estimator = "ipw",
    target = "S"
  )

  res <- contrast(
    fit,
    interventions = list(shifted = shift(delta), natural = NULL),
    reference = "natural",
    type = "difference",
    ci_method = "sandwich"
  )

  truth <- 3 * delta + delta * mean(d$L[d$S == 0])
  expect_equal(get_contrast_est(res), truth, tolerance = 0.20)
})


# -- IPW + ipsi + transportability: smoke test ----------------------------
# ipsi() is only for binary treatment, so use binary DGP

test_that("IPW + ipsi + transport: runs and produces finite estimates", {
  d <- simulate_transport(n = 6000, seed = 107)

  fit <- causat(
    d,
    "Y",
    "A",
    ~L,
    estimator = "ipw",
    target = "S"
  )

  res <- contrast(
    fit,
    interventions = list(a1 = static(1), ipsi_d = ipsi(0.5)),
    reference = "a1",
    type = "difference",
    ci_method = "sandwich"
  )

  expect_s3_class(res, "causatr_result")
  expect_true(is.finite(get_contrast_est(res)))
})


# -- AIPW + shift + transportability: truth recovery ----------------------

test_that("AIPW + shift + transport: recovers truth", {
  d <- simulate_mtp_transport(n = 8000, seed = 108)
  delta <- 1

  fit <- causat(
    d,
    "Y",
    "A",
    ~ L + A:L,
    estimator = "aipw",
    target = "S"
  )

  res <- contrast(
    fit,
    interventions = list(shifted = shift(delta), natural = NULL),
    reference = "natural",
    type = "difference",
    ci_method = "bootstrap",
    n_boot = 200
  )

  truth <- 3 * delta + delta * mean(d$L[d$S == 0])
  expect_equal(get_contrast_est(res), truth, tolerance = 0.20)
})


# -- AIPW + shift DR: wrong outcome, still consistent --------------------

test_that("AIPW + shift + transport: DR under wrong outcome model", {
  set.seed(109)
  n <- 8000
  L <- rnorm(n)
  S <- rbinom(n, 1, plogis(-0.5 + L))
  A <- ifelse(S == 1L, rnorm(n, 0.5 + 0.3 * L, 1), NA_real_)
  Y <- ifelse(
    S == 1L,
    2 + 3 * A + 1.5 * L + A * L + 0.5 * L^2 + rnorm(n),
    NA_real_
  )
  d <- data.frame(Y = Y, A = A, L = L, S = S)

  # Outcome model is misspecified (linear A:L, missing L^2 term).
  # Treatment model and sampling model are correct.
  # IPW augmentation term should correct the outcome-model bias.
  fit <- causat(
    d,
    "Y",
    "A",
    ~ L + A:L,
    estimator = "aipw",
    target = "S"
  )

  delta <- 1
  res <- contrast(
    fit,
    interventions = list(shifted = shift(delta), natural = NULL),
    reference = "natural",
    type = "difference",
    ci_method = "bootstrap",
    n_boot = 200
  )

  truth <- 3 * delta + delta * mean(d$L[d$S == 0])
  expect_equal(get_contrast_est(res), truth, tolerance = 0.25)
})


# -- AIPW + shift + transportability: smoke test --------------------------

test_that("AIPW + shift + transport: smoke test runs without error", {
  d <- simulate_mtp_transport(n = 4000, seed = 110)

  fit <- causat(
    d,
    "Y",
    "A",
    ~ L + A:L,
    estimator = "aipw",
    target = "S"
  )

  res <- contrast(
    fit,
    interventions = list(shifted = shift(0.5), natural = NULL),
    reference = "natural",
    ci_method = "bootstrap",
    n_boot = 50
  )

  expect_s3_class(res, "causatr_result")
  expect_true(is.finite(get_contrast_est(res)))
})


# -- Cross-estimator agreement: gcomp shift ~ IPW shift -------------------

test_that("gcomp shift ~ IPW shift under correct specification", {
  d <- simulate_mtp_transport(n = 8000, seed = 111)
  delta <- 1

  fit_gc <- causat(d, "Y", "A", ~ L + A:L, target = "S")
  fit_ipw <- causat(d, "Y", "A", ~L, estimator = "ipw", target = "S")

  res_gc <- contrast(
    fit_gc,
    interventions = list(shifted = shift(delta), natural = NULL),
    reference = "natural",
    ci_method = "bootstrap",
    n_boot = 100
  )
  res_ipw <- contrast(
    fit_ipw,
    interventions = list(shifted = shift(delta), natural = NULL),
    reference = "natural",
    ci_method = "sandwich"
  )

  expect_equal(
    get_contrast_est(res_gc),
    get_contrast_est(res_ipw),
    tolerance = 0.25
  )
})


# -- Bootstrap SE is finite and reasonable --------------------------------

test_that("gcomp + shift + transport: bootstrap SE is finite and reasonable", {
  d <- simulate_mtp_transport(n = 6000, seed = 112)

  fit <- causat(d, "Y", "A", ~ L + A:L, target = "S")

  res <- contrast(
    fit,
    interventions = list(shifted = shift(1), natural = NULL),
    reference = "natural",
    ci_method = "bootstrap",
    n_boot = 200
  )

  se <- res$contrasts$std_error
  expect_true(all(is.finite(se)))
  expect_true(all(se > 0))
  expect_true(all(se < 2))
})


# -- Sandwich rejection for gcomp + MTP + transport -----------------------

test_that("gcomp + shift + transport: sandwich is rejected", {
  d <- simulate_mtp_transport(n = 2000, seed = 113)

  fit <- causat(d, "Y", "A", ~ L + A:L, target = "S")

  expect_error(
    contrast(
      fit,
      interventions = list(shifted = shift(1), natural = NULL),
      reference = "natural",
      ci_method = "sandwich"
    ),
    class = "causatr_mtp_transport_sandwich"
  )
})


# -- Sandwich rejection for AIPW + MTP + transport ------------------------

test_that("AIPW + shift + transport: sandwich is rejected", {
  d <- simulate_mtp_transport(n = 2000, seed = 114)

  fit <- causat(d, "Y", "A", ~ L + A:L, estimator = "aipw", target = "S")

  expect_error(
    contrast(
      fit,
      interventions = list(shifted = shift(1), natural = NULL),
      reference = "natural",
      ci_method = "sandwich"
    ),
    class = "causatr_mtp_transport_sandwich"
  )
})


# -- IPW + shift + transport: sandwich IS supported -----------------------

test_that("IPW + shift + transport: sandwich works (no MC marginalization)", {
  d <- simulate_mtp_transport(n = 4000, seed = 115)

  fit <- causat(d, "Y", "A", ~L, estimator = "ipw", target = "S")

  res <- contrast(
    fit,
    interventions = list(shifted = shift(1), natural = NULL),
    reference = "natural",
    ci_method = "sandwich"
  )

  expect_s3_class(res, "causatr_result")
  se <- res$contrasts$std_error
  expect_true(all(is.finite(se)))
  expect_true(all(se > 0))
})


# -- Static intervention: no MC marginalization needed --------------------

test_that("gcomp + static + transport: sandwich still works (no MC needed)", {
  d <- simulate_mtp_transport(n = 4000, seed = 116)

  fit <- causat(d, "Y", "A", ~ L + A:L, target = "S")

  # static() sets A to a fixed value — no observed A needed on target rows
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich"
  )

  expect_s3_class(res, "causatr_result")
  expect_true(is.finite(get_contrast_est(res)))
})


# -- Binary treatment + dynamic + transport: exact marginalization --------

test_that("gcomp + dynamic + binary treatment + transport: exact marginalization", {
  set.seed(117)
  n <- 6000
  L <- rnorm(n)
  S <- rbinom(n, 1, plogis(-0.5 + L))
  A <- ifelse(S == 1L, rbinom(n, 1, plogis(0.3 * L)), NA_integer_)
  Y <- ifelse(S == 1L, 2 + 3 * A + 1.5 * L + A * L + rnorm(n), NA_real_)
  d <- data.frame(Y = Y, A = as.integer(A), L = L, S = S)

  fit <- causat(d, "Y", "A", ~ L + A:L, target = "S")

  # dynamic rule: "always treat" — d(data, trt) = 1
  # Same as static(1), so truth = 3 + E[L|S=0]
  always_treat <- dynamic(function(data, trt) rep(1L, nrow(data)))
  res <- contrast(
    fit,
    interventions = list(always = always_treat, a0 = static(0)),
    reference = "a0",
    ci_method = "bootstrap",
    n_boot = 100
  )

  truth <- 3 + mean(d$L[d$S == 0])
  expect_equal(get_contrast_est(res), truth, tolerance = 0.15)
})
