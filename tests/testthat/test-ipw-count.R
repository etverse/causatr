# Tests for count treatment (Poisson / negative binomial) IPW engine.
# Chunk 3j of Phase 4.
#
# DGP (simulate_count_treatment):
#   L ~ N(0, 1)
#   A | L ~ Poisson(exp(0.5 * L))
#   Y | A, L ~ N(2 + 1.5*A + L, sd = 1)
#
# True shift(delta) contrast = 1.5 * delta (linear in A).
# True shift(1) vs shift(0) = 1.5.

# ---- Treatment model unit tests ----

test_that("fit_treatment_model() with propensity_family = 'poisson' returns Poisson fit", {
  d <- simulate_count_treatment(n = 500, seed = 1)
  dt <- data.table::as.data.table(d)

  tm <- fit_treatment_model(
    dt,
    treatment = "A",
    confounders = ~L,
    propensity_family = "poisson"
  )

  expect_s3_class(tm, "causatr_treatment_model")
  expect_identical(tm$family, "poisson")
  expect_null(tm$sigma)
  expect_null(tm$theta)
  expect_equal(sum(tm$fit_rows), nrow(d))
  expect_equal(length(tm$alpha_hat), ncol(tm$X_prop))

  # Coefficients match a plain Poisson GLM on the same data
  ref_fit <- stats::glm(A ~ L, data = d, family = stats::poisson())
  expect_equal(
    unname(tm$alpha_hat),
    unname(stats::coef(ref_fit)),
    tolerance = 1e-10
  )
})

test_that("fit_treatment_model() with propensity_family = 'negbin' returns NB fit", {
  skip_if_not_installed("MASS")
  d <- simulate_count_treatment(n = 500, seed = 2)
  dt <- data.table::as.data.table(d)

  tm <- fit_treatment_model(
    dt,
    treatment = "A",
    confounders = ~L,
    model_fn = MASS::glm.nb,
    propensity_family = "negbin"
  )

  expect_s3_class(tm, "causatr_treatment_model")
  expect_identical(tm$family, "negbin")
  expect_null(tm$sigma)
  expect_true(is.numeric(tm$theta) && tm$theta > 0)
  expect_equal(sum(tm$fit_rows), nrow(d))
  expect_equal(length(tm$alpha_hat), ncol(tm$X_prop))

  # Coefficients match a plain MASS::glm.nb on the same data
  ref_fit <- MASS::glm.nb(A ~ L, data = d)
  expect_equal(
    unname(tm$alpha_hat),
    unname(stats::coef(ref_fit)),
    tolerance = 1e-10
  )
  expect_equal(tm$theta, ref_fit$theta, tolerance = 1e-6)
})

test_that("evaluate_density() Poisson returns dpois per row", {
  d <- simulate_count_treatment(n = 300, seed = 3)
  dt <- data.table::as.data.table(d)
  tm <- fit_treatment_model(
    dt,
    treatment = "A",
    confounders = ~L,
    propensity_family = "poisson"
  )

  fit_data <- dt[tm$fit_rows]
  a_obs <- fit_data$A
  lambda <- stats::predict(tm$model, newdata = fit_data, type = "response")

  f_obs <- evaluate_density(tm, a_obs, fit_data)
  expect_equal(f_obs, stats::dpois(a_obs, lambda), tolerance = 1e-12)

  # Shifted treatment values
  f_shift <- evaluate_density(tm, a_obs - 1L, fit_data)
  expect_equal(f_shift, stats::dpois(a_obs - 1L, lambda), tolerance = 1e-12)
})

test_that("evaluate_density() negbin returns dnbinom per row", {
  skip_if_not_installed("MASS")
  d <- simulate_count_treatment(n = 300, seed = 4)
  dt <- data.table::as.data.table(d)
  tm <- fit_treatment_model(
    dt,
    treatment = "A",
    confounders = ~L,
    model_fn = MASS::glm.nb,
    propensity_family = "negbin"
  )

  fit_data <- dt[tm$fit_rows]
  a_obs <- fit_data$A
  lambda <- stats::predict(tm$model, newdata = fit_data, type = "response")

  f_obs <- evaluate_density(tm, a_obs, fit_data)
  expect_equal(
    f_obs,
    stats::dnbinom(a_obs, mu = lambda, size = tm$theta),
    tolerance = 1e-12
  )
})

# ---- Truth-based Poisson test (T-count-poisson) ----

test_that("ipw × count trt × Poisson × shift(1) recovers ATE = 1.5", {
  # True shift(1) - shift(0) = 1.5 (linear in A, coefficient = 1.5).
  d <- simulate_count_treatment(n = 5000, seed = 10)

  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw",
    propensity_family = "poisson"
  )

  res <- contrast(
    fit,
    interventions = list(s1 = shift(1), s0 = NULL),
    reference = "s0",
    ci_method = "sandwich"
  )

  # Point estimate within 0.2 of truth
  expect_equal(res$contrasts$estimate[1], 1.5, tolerance = 0.2)
  # SEs are finite and positive
  expect_true(all(is.finite(res$estimates$se) & res$estimates$se > 0))
  # Contrast SE is finite and positive
  expect_true(is.finite(res$contrasts$se[1]) && res$contrasts$se[1] > 0)
})

# ---- NB parity test (T-count-negbin) ----

test_that("ipw × count trt × negbin agrees with Poisson on Poisson DGP", {
  # When the true DGP is Poisson, a NB fit (which nests Poisson at

  # theta -> Inf) should produce the same point estimate.
  skip_if_not_installed("MASS")
  d <- simulate_count_treatment(n = 5000, seed = 10)

  fit_pois <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw",
    propensity_family = "poisson"
  )
  res_pois <- contrast(
    fit_pois,
    interventions = list(s1 = shift(1), s0 = NULL),
    reference = "s0",
    ci_method = "sandwich"
  )

  fit_nb <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw",
    propensity_family = "negbin"
  )
  res_nb <- contrast(
    fit_nb,
    interventions = list(s1 = shift(1), s0 = NULL),
    reference = "s0",
    ci_method = "sandwich"
  )

  # Point estimates agree within 0.15 (NB nests Poisson; finite-sample
  # differences from estimating the extra theta parameter).
  expect_equal(
    res_nb$contrasts$estimate[1],
    res_pois$contrasts$estimate[1],
    tolerance = 0.15
  )
  # Both recover the truth
  expect_equal(res_nb$contrasts$estimate[1], 1.5, tolerance = 0.25)
})

# ---- Rejection tests (T-count-reject) ----

test_that("shift(0.5) on count treatment aborts", {
  d <- simulate_count_treatment(n = 200, seed = 20)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw",
    propensity_family = "poisson"
  )

  expect_snapshot(
    contrast(
      fit,
      interventions = list(s = shift(0.5)),
      ci_method = "sandwich"
    ),
    error = TRUE
  )
})

test_that("scale_by(2) on count treatment aborts (inverse not integer)", {
  # A / 2 is non-integer for odd A values.
  d <- simulate_count_treatment(n = 200, seed = 21)

  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw",
    propensity_family = "poisson"
  )

  expect_snapshot(
    contrast(
      fit,
      interventions = list(s = scale_by(2)),
      ci_method = "sandwich"
    ),
    error = TRUE
  )
})

test_that("static() on count treatment aborts", {
  d <- simulate_count_treatment(n = 200, seed = 22)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw",
    propensity_family = "poisson"
  )

  expect_snapshot(
    contrast(
      fit,
      interventions = list(s = static(2)),
      ci_method = "sandwich"
    ),
    error = TRUE
  )
})

test_that("dynamic() on count treatment aborts", {
  d <- simulate_count_treatment(n = 200, seed = 23)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw",
    propensity_family = "poisson"
  )

  expect_snapshot(
    contrast(
      fit,
      interventions = list(
        s = dynamic(function(data, a) pmax(a - 1L, 0L))
      ),
      ci_method = "sandwich"
    ),
    error = TRUE
  )
})

test_that("threshold() on count treatment aborts", {
  d <- simulate_count_treatment(n = 200, seed = 24)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw",
    propensity_family = "poisson"
  )

  expect_snapshot(
    contrast(
      fit,
      interventions = list(s = threshold(0, 5)),
      ci_method = "sandwich"
    ),
    error = TRUE
  )
})

test_that("ipsi() on count treatment aborts", {
  d <- simulate_count_treatment(n = 200, seed = 25)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw",
    propensity_family = "poisson"
  )

  expect_snapshot(
    contrast(
      fit,
      interventions = list(s = ipsi(2)),
      ci_method = "sandwich"
    ),
    error = TRUE
  )
})

# ---- Scale_by with integer-preserving factor succeeds ----

test_that("ipw × count trt × Poisson × scale_by(1) is a no-op", {
  # scale_by(1) is the identity intervention: a_eval = a_obs / 1 = a_obs
  # (trivially integer), so this must succeed.
  d <- simulate_count_treatment(n = 2000, seed = 30)

  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw",
    propensity_family = "poisson"
  )

  res <- contrast(
    fit,
    interventions = list(s1 = scale_by(1), s0 = NULL),
    reference = "s0",
    ci_method = "sandwich"
  )

  # scale_by(1) should return a contrast of ~0 vs natural course
  expect_equal(res$contrasts$estimate[1], 0, tolerance = 0.05)
})
