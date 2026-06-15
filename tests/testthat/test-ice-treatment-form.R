# Tests for the flexible-treatment ICE term (Phase 22b-5, `treatment_form =`).
# By default the treatment enters every per-step ICE outcome model as a bare
# numeric main effect; a kinked or non-monotone counterfactual dose-response is
# then misspecified. `treatment_form = ~ factor(A)` / `~ splines::ns(A, df)` lets
# the treatment enter flexibly while the intervention still sets the *numeric*
# treatment column -- only the model's design term changes. The truth oracles
# below are exact analytic contrasts on data-generating processes with NO
# treatment -> covariate feedback, so the marginal counterfactual mean has a
# closed form and the contrast cancels the (model-independent) covariate term.

# ---------------------------------------------------------------------------
# DGPs with exogenous time-varying covariates (A does not affect future L), so
# E[Y^{static(a)}] = b0 + tau * f(a) + 0.5 * E[L0] is analytic and the contrast
# between two regimes drops the covariate term entirely.
# ---------------------------------------------------------------------------

# Categorical treatment {0,1,2} with a NON-MONOTONE level effect f = (0, 3, 1):
# the bare-numeric model fits a single slope through a non-monotone response and
# is badly misspecified; factor(A) carries one coefficient per level.
make_cat_ice <- function(n, seed) {
  set.seed(seed)
  tau <- 2L
  f_lvl <- c(0, 3, 1) # f(0), f(1), f(2)
  L0 <- stats::rnorm(n)
  A <- matrix(0L, n, tau)
  L <- matrix(0, n, tau)
  for (t in seq_len(tau)) {
    Lt <- 0.5 * L0 + stats::rnorm(n) # exogenous: independent of A
    At <- pmin(2L, stats::rpois(n, exp(-0.3 + 0.5 * Lt))) # confounded by L
    A[, t] <- At
    L[, t] <- Lt
  }
  Y <- 1 + f_lvl[A[, 1] + 1] + f_lvl[A[, 2] + 1] + 0.5 * L0 + stats::rnorm(n)
  d <- data.frame(
    id = rep(seq_len(n), each = tau),
    t = rep(seq_len(tau), n),
    L0 = rep(L0, each = tau),
    A = as.vector(t(A)),
    L = as.vector(t(L)),
    Y = NA_real_
  )
  d$Y[d$t == tau] <- Y
  d
}

# Continuous treatment with a CURVED (quadratic) dose-response (A - 1)^2: a
# linear-in-A model cannot represent the curvature, so predicting under a shift
# is biased; splines::ns(A, df) recovers it.
make_ns_ice <- function(n, seed) {
  set.seed(seed)
  tau <- 2L
  L0 <- stats::rnorm(n)
  A <- matrix(0, n, tau)
  L <- matrix(0, n, tau)
  for (t in seq_len(tau)) {
    Lt <- 0.5 * L0 + stats::rnorm(n) # exogenous
    At <- 0.5 * Lt + stats::rnorm(n) # continuous, confounded by L
    A[, t] <- At
    L[, t] <- Lt
  }
  Y <- 1 + (A[, 1] - 1)^2 + (A[, 2] - 1)^2 + 0.5 * L0 + stats::rnorm(n)
  d <- data.frame(
    id = rep(seq_len(n), each = tau),
    t = rep(seq_len(tau), n),
    L0 = rep(L0, each = tau),
    A = as.vector(t(A)),
    L = as.vector(t(L)),
    Y = NA_real_
  )
  d$Y[d$t == tau] <- Y
  d
}

fit_tf <- function(d, treatment_form = NULL) {
  causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    id = "id",
    time = "t",
    estimator = "gcomp",
    history = Inf,
    treatment_form = treatment_form
  )
}

est_of <- function(res, name) {
  res$estimates$estimate[res$estimates$intervention == name]
}

# ---------------------------------------------------------------------------
# Backward compatibility + formula construction
# ---------------------------------------------------------------------------

test_that("treatment_form = NULL and ~ A produce an identical estimate", {
  # Explicitly naming the bare term must reproduce the default to numerical
  # precision -- the substitution path collapses `A` -> `lag1_A` exactly like
  # the historical string concatenation.
  d <- make_cat_ice(2000L, 1L)
  fit_default <- fit_tf(d, NULL)
  fit_bare <- fit_tf(d, ~A)
  mu_default <- mean(ice_iterate(fit_default, static(1))$pseudo_final)
  mu_bare <- mean(ice_iterate(fit_bare, static(1))$pseudo_final)
  expect_equal(mu_default, mu_bare, tolerance = 1e-12)
})

test_that("ice_build_formula expands a flexible treatment term across lags", {
  # At time index 2 with max_lag 2, `factor(A)` must appear at the current
  # period and at both lags (`factor(lag1_A)`, `factor(lag2_A)`), with the TV
  # confounder and baseline terms intact.
  dat <- data.table::data.table(
    A = c(0, 1, 2),
    lag1_A = c(0, 0, 1),
    lag2_A = c(0, 0, 0),
    L = c(0.1, 0.2, 0.3),
    lag1_L = c(0, 0.1, 0.2),
    lag2_L = c(0, 0, 0.1),
    L0 = c(1, 1, 1)
  )
  f <- ice_build_formula(
    response = "Y",
    treatment = "A",
    baseline_terms = "L0",
    tv_vars = "L",
    tv_terms = "L",
    time_idx = 2L,
    max_lag = 2L,
    data_at_time = dat,
    em_info = NULL,
    treatment_terms = "factor(A)"
  )
  labels <- attr(stats::terms(f), "term.labels")
  expect_true(all(
    c("factor(A)", "factor(lag1_A)", "factor(lag2_A)", "L0") %in% labels
  ))
  # The bare-numeric default must NOT introduce a factor term.
  f_bare <- ice_build_formula(
    response = "Y",
    treatment = "A",
    baseline_terms = "L0",
    tv_vars = "L",
    tv_terms = "L",
    time_idx = 1L,
    max_lag = 2L,
    data_at_time = dat,
    em_info = NULL,
    treatment_terms = NULL
  )
  expect_true(all(
    c("A", "lag1_A") %in% attr(stats::terms(f_bare), "term.labels")
  ))
  expect_false(any(grepl(
    "factor",
    attr(stats::terms(f_bare), "term.labels")
  )))
})

# ---------------------------------------------------------------------------
# Truth-based: factor(A) and ns(A) recover the analytic truth; the bare model
# does not. The contrast cancels the model-independent baseline/covariate term.
# ---------------------------------------------------------------------------

test_that("factor(A) recovers a non-monotone categorical dose-response; bare-numeric does not", {
  # Truth: static(1) - static(0) = 2 * (f(1) - f(0)) = 2 * 3 = 6, independent of
  # the L distribution (the baseline term cancels in the contrast).
  d <- make_cat_ice(8000L, 1L)
  ivs <- list(a1 = static(1), a0 = static(0))

  fit_factor <- fit_tf(d, ~ factor(A))
  rf <- contrast(fit_factor, interventions = ivs, type = "difference")
  diff_factor <- est_of(rf, "a1") - est_of(rf, "a0")
  # Correctly specified factor model: the finite-sample gap to truth is ~0.013
  # at this n; assert it stays within a few sampling SE.
  expect_lt(abs(diff_factor - 6), 0.04)

  fit_bare <- fit_tf(d, NULL)
  rb <- contrast(fit_bare, interventions = ivs, type = "difference")
  diff_bare <- est_of(rb, "a1") - est_of(rb, "a0")
  # The bare-numeric slope through a non-monotone response is ~4.3 off truth and
  # strictly worse than factor(A) (engine-necessity for the flexible term).
  expect_gt(abs(diff_bare - 6), 2)
  expect_lt(abs(diff_factor - 6), abs(diff_bare - 6))
})

test_that("splines::ns(A) recovers a curved continuous dose-response; linear does not", {
  # Truth: shift(0.5) - shift(0) = sum_t E[2*delta*(A_t - 1) + delta^2]
  #      = tau * (-2*delta + delta^2) = 2 * (-1 + 0.25) = -1.5, independent of
  # Var(A) (the contrast cancels the variance and baseline terms).
  d <- make_ns_ice(8000L, 1L)
  ivs <- list(s5 = shift(0.5), s0 = shift(0))

  fit_ns <- fit_tf(d, ~ splines::ns(A, df = 4))
  rn <- contrast(fit_ns, interventions = ivs, type = "difference")
  diff_ns <- est_of(rn, "s5") - est_of(rn, "s0")
  # The df=4 spline captures the curvature; gap to truth is ~0.006 at this n.
  expect_lt(abs(diff_ns - (-1.5)), 0.03)

  fit_lin <- fit_tf(d, NULL)
  rl <- contrast(fit_lin, interventions = ivs, type = "difference")
  diff_lin <- est_of(rl, "s5") - est_of(rl, "s0")
  # A linear-in-A model misses the curvature and is ~0.5 off truth.
  expect_gt(abs(diff_lin - (-1.5)), 0.3)
  expect_lt(abs(diff_ns - (-1.5)), abs(diff_lin - (-1.5)))
})

# ---------------------------------------------------------------------------
# Variance: the analytic ICE sandwich handles a flexible design matrix (factor
# terms are just extra columns), so it agrees with the ID-cluster bootstrap.
# ---------------------------------------------------------------------------

test_that("plain-ICE factor(A) sandwich SE agrees with the bootstrap SE", {
  skip_if_fast()
  d <- make_cat_ice(3000L, 5L)
  fit <- fit_tf(d, ~ factor(A))
  ivs <- list(a1 = static(1), a0 = static(0))
  rs <- contrast(
    fit,
    interventions = ivs,
    type = "difference",
    ci_method = "sandwich"
  )
  set.seed(1)
  rb <- contrast(
    fit,
    interventions = ivs,
    type = "difference",
    ci_method = "bootstrap",
    n_boot = 150L
  )
  se_s <- rs$estimates$se[rs$estimates$intervention == "a1"]
  se_b <- rb$estimates$se[rb$estimates$intervention == "a1"]
  expect_true(is.finite(se_s) && se_s > 0)
  # Relative agreement; the bootstrap SE carries Monte-Carlo noise at n_boot=150
  # (~6%), so a 15% band is a couple of MC-SE around the exact sandwich value.
  expect_equal(se_s, se_b, tolerance = 0.15)
})

# ---------------------------------------------------------------------------
# Rejections
# ---------------------------------------------------------------------------

test_that("treatment_form is rejected outside longitudinal g-computation", {
  dpt <- data.frame(Y = rnorm(80), A = rbinom(80, 1, 0.5), L = rnorm(80))
  expect_error(
    causat(
      dpt,
      outcome = "Y",
      treatment = "A",
      confounders = ~L,
      estimator = "gcomp",
      treatment_form = ~ factor(A)
    ),
    class = "causatr_treatment_form_not_ice"
  )

  d <- make_cat_ice(400L, 2L)
  expect_error(
    causat(
      d,
      outcome = "Y",
      treatment = "A",
      confounders = ~L0,
      confounders_tv = ~L,
      id = "id",
      time = "t",
      estimator = "ipw",
      treatment_form = ~ factor(A)
    ),
    class = "causatr_treatment_form_not_ice"
  )
})

test_that("treatment_form rejects non-treatment variables and non-formula input", {
  d <- make_cat_ice(400L, 3L)
  # A covariate transform belongs in `confounders` / `confounders_tv`.
  expect_error(
    fit_tf(d, ~ factor(L)),
    class = "causatr_treatment_form_bad"
  )
  # A two-sided formula is a structural misuse.
  expect_error(
    fit_tf(d, Y ~ factor(A)),
    class = "causatr_treatment_form_bad"
  )
  # A non-formula is rejected on shape.
  expect_error(
    fit_tf(d, "factor(A)"),
    class = "causatr_treatment_form_bad"
  )
})

test_that("factor(A) errors when an intervention leaves the observed support", {
  # `factor(A)` stores the observed levels; a static regime to an unseen level
  # makes predict() fail (documented constraint). ns(A) extrapolates instead.
  d <- make_cat_ice(500L, 4L)
  fit <- fit_tf(d, ~ factor(A))
  expect_error(
    contrast(fit, interventions = list(a = static(3)), type = "difference"),
    regexp = "new level"
  )
})
