# Complex DGP tests: nonlinear confounding, multiple confounders, GAMs, splines.
# These test that the estimation methods work correctly when the true
# data-generating mechanism is more realistic than the simple linear DGPs
# in helper-dgp.R.

# DGP: Nonlinear confounding -----------------------------------------
# True ATE = 3 (constant treatment effect, nonlinear confounding).
# GLM Y ~ A + L is misspecified; correct model needs splines or GAM.
simulate_nonlinear <- function(n = 3000, seed = 42) {
  set.seed(seed)
  L <- rnorm(n)
  ps <- plogis(sin(L) + 0.3 * L^2)
  A <- rbinom(n, 1, ps)
  Y <- 1 + 3 * A + sin(2 * L) + L^2 + rnorm(n)
  data.frame(Y = Y, A = A, L = L)
}

# DGP: Multiple confounders + interactions ----------------------------
# True ATE = 3. Requires all confounders + L1:L2 interaction.
simulate_multi_confounder <- function(n = 3000, seed = 42) {
  set.seed(seed)
  L1 <- rnorm(n)
  L2 <- rnorm(n)
  L3 <- rnorm(n)
  L4 <- rbinom(n, 1, 0.5)
  ps <- plogis(0.5 * L1 + 0.3 * L2 - 0.2 * L3 + 0.4 * L4)
  A <- rbinom(n, 1, ps)
  Y <- 2 + 3 * A + L1 + 0.5 * L2 + L1 * L2 + 0.8 * L4 + rnorm(n)
  data.frame(Y = Y, A = A, L1 = L1, L2 = L2, L3 = L3, L4 = L4)
}


# ============================================================
# NONLINEAR CONFOUNDING × GCOMP
# ============================================================

test_that("gcomp × nonlinear confounding × GLM with splines recovers ATE", {
  d <- simulate_nonlinear(n = 3000, seed = 42)
  fit <- causat(d, outcome = "Y", treatment = "A",
    confounders = ~ splines::ns(L, 5)
  )
  res <- contrast(fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0"
  )
  expect_equal(res$contrasts$estimate, 3, tolerance = 0.15)
})

test_that("gcomp × nonlinear confounding × misspecified GLM is biased", {
  d <- simulate_nonlinear(n = 3000, seed = 42)
  fit <- causat(d, outcome = "Y", treatment = "A", confounders = ~L)
  res <- contrast(fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0"
  )
  expect_true(abs(res$contrasts$estimate - 3) > 0.3)
})

test_that("gcomp × nonlinear confounding × GAM recovers ATE", {
  skip_if_not_installed("mgcv")
  d <- simulate_nonlinear(n = 3000, seed = 42)
  fit <- causat(d, outcome = "Y", treatment = "A",
    confounders = ~ s(L),
    model_fn = mgcv::gam
  )
  res <- contrast(fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0"
  )
  expect_equal(res$contrasts$estimate, 3, tolerance = 0.15)
})


# ============================================================
# GAM × SANDWICH vs BOOTSTRAP
# ============================================================

test_that("gcomp × GAM × sandwich vs bootstrap SE agreement", {
  skip_if_not_installed("mgcv")
  d <- simulate_nonlinear(n = 2000, seed = 42)
  fit <- causat(d, outcome = "Y", treatment = "A",
    confounders = ~ s(L),
    model_fn = mgcv::gam
  )
  res_sw <- contrast(fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich"
  )
  res_bt <- suppressWarnings(contrast(fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "bootstrap", n_boot = 200
  ))
  ratio <- res_sw$contrasts$se / res_bt$contrasts$se
  expect_true(ratio > 0.5 && ratio < 2.0)
})

test_that("gcomp × GLM with splines × sandwich vs bootstrap SE agreement", {
  d <- simulate_nonlinear(n = 2000, seed = 42)
  fit <- causat(d, outcome = "Y", treatment = "A",
    confounders = ~ splines::ns(L, 4)
  )
  res_sw <- contrast(fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich"
  )
  res_bt <- suppressWarnings(contrast(fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "bootstrap", n_boot = 200
  ))
  ratio <- res_sw$contrasts$se / res_bt$contrasts$se
  expect_true(ratio > 0.5 && ratio < 2.0)
})


# ============================================================
# MULTIPLE CONFOUNDERS × ALL METHODS
# ============================================================

test_that("gcomp × multiple confounders + interactions recovers ATE", {
  d <- simulate_multi_confounder(n = 3000, seed = 42)
  fit <- causat(d, outcome = "Y", treatment = "A",
    confounders = ~ L1 + L2 + L3 + L4 + L1:L2
  )
  res <- contrast(fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0"
  )
  expect_equal(res$contrasts$estimate, 3, tolerance = 0.15)
})

test_that("ipw × multiple confounders recovers ATE", {
  d <- simulate_multi_confounder(n = 3000, seed = 42)
  fit <- causat(d, outcome = "Y", treatment = "A",
    confounders = ~ L1 + L2 + L3 + L4,
    method = "ipw"
  )
  res <- contrast(fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0"
  )
  expect_equal(res$contrasts$estimate, 3, tolerance = 0.15)
})


# ============================================================
# TRIANGULATION × NONLINEAR CONFOUNDING
# ============================================================

test_that("triangulation × nonlinear confounding × flexible models", {
  skip_if_not_installed("mgcv")
  d <- simulate_nonlinear(n = 3000, seed = 42)

  fit_gc <- causat(d, outcome = "Y", treatment = "A",
    confounders = ~ splines::ns(L, 5)
  )
  fit_ipw <- causat(d, outcome = "Y", treatment = "A",
    confounders = ~ splines::ns(L, 5),
    method = "ipw"
  )
  fit_gam <- causat(d, outcome = "Y", treatment = "A",
    confounders = ~ s(L),
    model_fn = mgcv::gam
  )

  ate_gc <- contrast(fit_gc,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0"
  )$contrasts$estimate
  ate_ipw <- contrast(fit_ipw,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0"
  )$contrasts$estimate
  ate_gam <- contrast(fit_gam,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0"
  )$contrasts$estimate

  expect_equal(ate_gc, 3, tolerance = 0.15)
  expect_equal(ate_ipw, 3, tolerance = 0.15)
  expect_equal(ate_gam, 3, tolerance = 0.15)
})
