test_that("causat(estimator = 'ipw') fits on simple data", {
  df <- data.frame(
    Y = c(1, 2, 3, 4, 5, 6, 7, 8),
    A = c(0, 0, 0, 0, 1, 1, 1, 1),
    L = c(1, 2, 3, 4, 1, 2, 3, 4)
  )
  fit <- suppressWarnings(causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw"
  ))
  expect_s3_class(fit, "causatr_fit")
  expect_equal(fit$estimator, "ipw")
  # Self-contained IPW: the treatment density model lives in `$details`;
  # `$weights_obj` is NULL and `$model` is only a display placeholder.
  expect_true(!is.null(fit$details$treatment_model))
  expect_true(!is.null(fit$details$propensity_model))
  expect_null(fit$weights_obj)
})

test_that("IPW with external weights stores in details", {
  d <- simulate_binary_continuous(n = 500, seed = 42)
  w <- runif(nrow(d), 0.5, 2)
  # `suppressWarnings()` here silences GLM's "non-integer #successes
  # in a binomial glm" warning, which fires whenever external
  # continuous-valued weights enter a logistic propensity fit. It is
  # harmless (the GLM fits exactly the same quasibinomial-style
  # weighted score equations) and appears for every weighted IPW
  # fit, so we keep the signal in the test clean.
  fit <- suppressWarnings(causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw",
    weights = w
  ))
  expect_true(!is.null(fit$details$weights))
  expect_equal(length(fit$details$weights), nrow(d))
  expect_false(".causatr_w" %in% names(fit$data))
})

test_that("IPW bootstrap with external weights gives finite SE", {
  d <- simulate_binary_continuous(n = 500, seed = 42)
  w <- runif(nrow(d), 0.5, 2)
  fit <- suppressWarnings(causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw",
    weights = w
  ))
  res <- suppressWarnings(contrast(
    fit,
    interventions = list(a0 = static(0), a1 = static(1)),
    ci_method = "bootstrap",
    n_boot = 50
  ))
  expect_true(all(is.finite(res$contrasts$se)))
  expect_true(res$contrasts$se > 0)
})

test_that("IPW continuous treatment fits and contrasts", {
  d <- simulate_continuous_continuous(n = 2000, seed = 42)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw"
  )
  expect_s3_class(fit, "causatr_fit")
  expect_equal(fit$estimator, "ipw")
  expect_equal(fit$details$treatment_model$family, "gaussian")
})

test_that("causat(estimator = 'ipw') rejects longitudinal data", {
  df <- data.frame(
    Y = c(0, 1, 0, 1),
    A = c(0, 1, 0, 1),
    L = c(1, 1, 2, 2),
    id = c(1, 1, 2, 2),
    time = c(0, 1, 0, 1)
  )
  expect_error(
    causat(
      df,
      outcome = "Y",
      treatment = "A",
      confounders = ~L,
      estimator = "ipw",
      id = "id",
      time = "time"
    ),
    "longitudinal"
  )
})

test_that("IPW with propensity_model_fn = mgcv::gam fits", {
  skip_if_not_installed("mgcv")
  d <- simulate_binary_continuous(n = 400, seed = 11)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~ s(L),
    estimator = "ipw",
    propensity_model_fn = mgcv::gam
  )
  expect_s3_class(fit, "causatr_fit")
  expect_true(inherits(fit$details$propensity_model, "gam"))
})

test_that("IPW rejects static(1) on a continuous treatment", {
  d <- simulate_continuous_continuous(n = 400, seed = 1)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw"
  )
  expect_snapshot(
    error = TRUE,
    contrast(
      fit,
      interventions = list(a = static(1), b = static(0)),
      ci_method = "sandwich"
    )
  )
})

test_that("IPW accepts A:modifier interaction terms in confounders", {
  # EM terms expand the per-intervention MSM from `Y ~ 1` to
  # `Y ~ 1 + modifier` (chunk 6b). The propensity formula strips the
  # EM term automatically via `build_ps_formula()`.
  d <- simulate_binary_continuous(n = 200, seed = 1)
  d$sex <- rbinom(nrow(d), 1, 0.5)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L + sex + A:sex,
    estimator = "ipw"
  )
  expect_s3_class(fit, "causatr_fit")
  expect_true(fit$details$em_info$has_em)
})

test_that("IPW rejects bare treatment in confounders", {
  d <- simulate_binary_continuous(n = 200, seed = 2)
  expect_error(
    causat(
      d,
      outcome = "Y",
      treatment = "A",
      confounders = ~ L + A,
      estimator = "ipw"
    ),
    class = "causatr_bare_treatment_in_confounders"
  )
})

# `fit_ipw()` rejects scientifically out-of-scope combinations up
# front: longitudinal IPW (Phase 10 work) and multivariate treatment
# (Phase 8 work). The aborts ship today as actionable error messages
# pointing users to the supported alternatives.
test_that("IPW rejects type = 'longitudinal' upstream and via fit_ipw()", {
  # Two paths: `check_causat_inputs()` rejects it at the top of
  # `causat()` (the user-facing path), and `fit_ipw()` has its own
  # defensive abort for the same condition (the dispatcher path used
  # by the bootstrap refitter). The test exercises both.
  d <- simulate_binary_continuous(n = 100, seed = 411)
  d$id <- seq_len(nrow(d))
  d$t <- 0L
  expect_error(
    causat(
      d,
      outcome = "Y",
      treatment = "A",
      confounders = ~L,
      estimator = "ipw",
      id = "id",
      time = "t",
      type = "longitudinal"
    ),
    "does not support longitudinal data"
  )
  # Direct call to fit_ipw with type = "longitudinal" hits the
  # internal abort branch.
  expect_error(
    fit_ipw(
      data = data.table::as.data.table(d),
      outcome = "Y",
      treatment = "A",
      confounders = ~L,
      family = "gaussian",
      estimand = "ATE",
      type = "longitudinal",
      call = NULL
    ),
    "Longitudinal IPW is not supported"
  )
})

test_that("IPW accepts multivariate treatment via sequential factorisation", {
  set.seed(412)
  d <- data.frame(
    Y = rnorm(100),
    A1 = rbinom(100, 1, 0.5),
    A2 = rbinom(100, 1, 0.5),
    L = rnorm(100)
  )
  fit <- causat(
    d,
    outcome = "Y",
    treatment = c("A1", "A2"),
    confounders = ~L,
    estimator = "ipw"
  )
  expect_true(isTRUE(fit$details$is_multivariate))
  expect_length(fit$details$treatment_models, 2L)
  expect_true(inherits(
    fit$details$treatment_models,
    "causatr_treatment_models"
  ))
})

# Defensive guard: `compute_ipw_contrast_point()` aborts with an
# "Internal error" if a caller hands it an IPW fit whose
# `$details$treatment_model` is missing. `fit_ipw()` always populates it,
# so this branch only fires if a downstream consumer mutates the fit
# object. The test strips the slot to exercise the guard.
test_that("compute_ipw_contrast_point aborts on missing treatment_model", {
  d <- simulate_binary_continuous(n = 100, seed = 3)
  fit <- suppressWarnings(causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw"
  ))
  fit$details$treatment_model <- NULL
  expect_error(
    compute_ipw_contrast_point(
      fit,
      interventions = list(a1 = static(1), a0 = static(0)),
      target_idx = rep(TRUE, nrow(d))
    ),
    "treatment_model"
  )
})
