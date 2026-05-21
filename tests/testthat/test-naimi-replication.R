# Naimi, Cole & Kennedy (2017)-inspired longitudinal DGP.
#
# Two time points, binary treatment, binary time-varying confounder,
# continuous outcome (CD4 count analogue). Treatment-confounder feedback
# at t=1 makes naive regression biased; ICE g-computation, longitudinal
# IPW, and longitudinal AIPW should all recover the MC truth.

# -- ICE g-comp recovers MC truth -----------------------------------------

test_that("Naimi DGP: ICE recovers ATE(always vs never) ≈ MC truth", {
  sim <- simulate_naimi_longitudinal(n = 5000, seed = 42)
  d <- sim$data
  truth <- sim$truth_ate

  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~age,
    confounders_tv = ~Ltv,
    id = "id",
    time = "time"
  )

  res <- contrast(
    fit,
    interventions = list(always = static(1), never = static(0)),
    reference = "never",
    ci_method = "sandwich"
  )

  ate <- res$contrasts$estimate[1]
  expect_equal(ate, truth, tolerance = 0.02)
  expect_true(is.finite(res$contrasts$se[1]) && res$contrasts$se[1] > 0)
})


# -- Longitudinal IPW recovers MC truth -----------------------------------

test_that("Naimi DGP: longitudinal IPW recovers ATE ≈ MC truth", {
  sim <- simulate_naimi_longitudinal(n = 5000, seed = 42)
  d <- sim$data
  truth <- sim$truth_ate

  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~age,
    confounders_tv = ~Ltv,
    id = "id",
    time = "time",
    estimator = "ipw"
  )

  res <- contrast(
    fit,
    interventions = list(always = static(1), never = static(0)),
    reference = "never",
    ci_method = "sandwich"
  )

  ate <- res$contrasts$estimate[1]
  expect_equal(ate, truth, tolerance = 0.02)
  expect_true(is.finite(res$contrasts$se[1]) && res$contrasts$se[1] > 0)
})


# -- Longitudinal AIPW recovers MC truth ----------------------------------

test_that("Naimi DGP: longitudinal AIPW recovers ATE ≈ MC truth", {
  sim <- simulate_naimi_longitudinal(n = 5000, seed = 42)
  d <- sim$data
  truth <- sim$truth_ate

  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~age,
    confounders_tv = ~Ltv,
    id = "id",
    time = "time",
    estimator = "aipw",
    propensity_model_fn = stats::glm
  )

  res <- contrast(
    fit,
    interventions = list(always = static(1), never = static(0)),
    reference = "never",
    ci_method = "sandwich"
  )

  ate <- res$contrasts$estimate[1]
  expect_equal(ate, truth, tolerance = 0.02)
  expect_true(is.finite(res$contrasts$se[1]) && res$contrasts$se[1] > 0)
})


# -- Cross-method agreement on Naimi DGP ----------------------------------

test_that("Naimi DGP: ICE, IPW, AIPW agree within 2", {
  sim <- simulate_naimi_longitudinal(n = 5000, seed = 42)
  d <- sim$data

  ivs <- list(always = static(1), never = static(0))

  fit_ice <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~age,
    confounders_tv = ~Ltv,
    id = "id",
    time = "time"
  )
  fit_ipw <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~age,
    confounders_tv = ~Ltv,
    id = "id",
    time = "time",
    estimator = "ipw"
  )
  fit_aipw <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~age,
    confounders_tv = ~Ltv,
    id = "id",
    time = "time",
    estimator = "aipw",
    propensity_model_fn = stats::glm
  )

  res_ice <- contrast(fit_ice, interventions = ivs, reference = "never")
  res_ipw <- contrast(fit_ipw, interventions = ivs, reference = "never")
  res_aipw <- contrast(fit_aipw, interventions = ivs, reference = "never")

  ates <- c(
    ice = res_ice$contrasts$estimate[1],
    ipw = res_ipw$contrasts$estimate[1],
    aipw = res_aipw$contrasts$estimate[1]
  )

  pairs <- combn(ates, 2)
  for (j in seq_len(ncol(pairs))) {
    expect_lt(abs(pairs[1, j] - pairs[2, j]), 2)
  }
})


# -- Naive regression should be biased ------------------------------------

test_that("Naimi DGP: naive cross-sectional regression is biased", {
  sim <- simulate_naimi_longitudinal(n = 5000, seed = 42)
  truth <- sim$truth_ate

  d_final <- sim$data[sim$data$time == 1, ]

  naive_fit <- stats::lm(Y ~ A + Ltv + age, data = d_final)
  naive_ate <- stats::coef(naive_fit)["A"]

  # Naive estimate should be substantially different from truth.
  expect_gt(abs(naive_ate - truth), 100)
})
