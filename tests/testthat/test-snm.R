# --- Chunk 18a: routing, validation, and rejection tests --------------------

test_that("causat() accepts estimator = 'snm' and returns causatr_fit", {
  d <- simulate_binary_continuous(n = 500, seed = 42)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "snm",
    propensity_model_fn = stats::glm
  )
  expect_s3_class(fit, "causatr_fit")
  expect_equal(fit$estimator, "snm")
  expect_equal(fit$type, "point")
})

test_that("SNM fit stores treatment model and blip specification", {
  d <- simulate_binary_continuous(n = 500, seed = 42)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "snm",
    propensity_model_fn = stats::glm
  )
  expect_true(!is.null(fit$details$treatment_model))
  expect_s3_class(fit$details$treatment_model, "causatr_treatment_model")
  expect_true(!is.null(fit$details$blip_spec))
  expect_equal(fit$details$blip_type, "linear")
  expect_equal(fit$details$blip_spec$param_names, "psi_intercept")
  expect_equal(fit$details$blip_spec$n_params, 1L)
})

test_that("SNM fit with EM terms stores modifier columns in blip spec", {
  d <- simulate_effect_mod(n = 500, seed = 42)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L + sex + A:sex,
    estimator = "snm",
    propensity_model_fn = stats::glm
  )
  expect_equal(
    fit$details$blip_spec$modifier_cols,
    "sex"
  )
  expect_equal(
    fit$details$blip_spec$param_names,
    c("psi_intercept", "psi_sex")
  )
  expect_equal(fit$details$blip_spec$n_params, 2L)
  expect_true(fit$details$em_info$has_em)
})

test_that("SNM fit with per-component confounders works", {
  d <- simulate_effect_mod(n = 500, seed = 42)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders_outcome = ~ L + sex + A:sex,
    confounders_treatment = ~ L + sex,
    estimator = "snm",
    propensity_model_fn = stats::glm
  )
  expect_s3_class(fit, "causatr_fit")
  expect_equal(fit$details$blip_spec$modifier_cols, "sex")
  expect_equal(fit$details$blip_spec$n_params, 2L)
})

test_that("SNM without EM terms informs user (blip = constant ATE)", {
  d <- simulate_binary_continuous(n = 500, seed = 42)
  expect_message(
    causat(
      d,
      outcome = "Y",
      treatment = "A",
      confounders = ~L,
      estimator = "snm",
      propensity_model_fn = stats::glm
    ),
    class = "causatr_snm_no_em"
  )
})

test_that("SNM with EM terms does not emit no-EM message", {
  d <- simulate_effect_mod(n = 500, seed = 42)
  expect_no_message(
    causat(
      d,
      outcome = "Y",
      treatment = "A",
      confounders = ~ L + sex + A:sex,
      estimator = "snm",
      propensity_model_fn = stats::glm
    ),
    class = "causatr_snm_no_em"
  )
})

test_that("SNM supports continuous treatment", {
  d <- simulate_continuous_continuous(n = 500, seed = 42)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "snm",
    propensity_model_fn = stats::glm
  )
  expect_s3_class(fit, "causatr_fit")
  expect_equal(fit$details$treatment_model$family, "gaussian")
})

test_that("SNM fit stores data and outcome correctly", {
  d <- simulate_binary_continuous(n = 500, seed = 42)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "snm",
    propensity_model_fn = stats::glm
  )
  expect_equal(fit$outcome, "Y")
  expect_equal(fit$treatment, "A")
  expect_equal(nrow(fit$data), 500L)
  expect_null(fit$model)
})

# --- Rejection paths --------------------------------------------------------

test_that("SNM rejects multivariate treatment", {
  d <- data.frame(
    Y = rnorm(100),
    A1 = rbinom(100, 1, 0.5),
    A2 = rbinom(100, 1, 0.5),
    L = rnorm(100)
  )
  expect_error(
    causat(
      d,
      outcome = "Y",
      treatment = c("A1", "A2"),
      confounders = ~L,
      estimator = "snm",
      propensity_model_fn = stats::glm
    ),
    class = "causatr_snm_multivariate"
  )
})

test_that("SNM rejects longitudinal data (pending 18d)", {
  d <- data.frame(
    id = rep(1:50, each = 2),
    time = rep(0:1, 50),
    Y = rnorm(100),
    A = rbinom(100, 1, 0.5),
    L = rnorm(100)
  )
  suppressWarnings(
    expect_error(
      causat(
        d,
        outcome = "Y",
        treatment = "A",
        confounders = ~L,
        estimator = "snm",
        id = "id",
        time = "time",
        propensity_model_fn = stats::glm
      ),
      class = "causatr_snm_longitudinal_pending"
    )
  )
})

test_that("contrast() rejects interventions for SNM fit", {
  d <- simulate_binary_continuous(n = 500, seed = 42)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "snm",
    propensity_model_fn = stats::glm
  )
  expect_error(
    contrast(
      fit,
      interventions = list(a1 = static(1), a0 = static(0))
    ),
    class = "causatr_snm_no_interventions"
  )
})

test_that("SNM requires treatment confounders", {
  d <- simulate_binary_continuous(n = 500, seed = 42)
  expect_error(
    causat(
      d,
      outcome = "Y",
      treatment = "A",
      estimator = "snm",
      propensity_model_fn = stats::glm
    ),
    "requires treatment-model confounders"
  )
})

test_that("SNM does not emit model_fn default warning", {
  d <- simulate_binary_continuous(n = 500, seed = 42)
  withr::local_options(causatr.suppress_default_warnings = FALSE)
  expect_no_warning(
    causat(
      d,
      outcome = "Y",
      treatment = "A",
      confounders = ~L,
      estimator = "snm",
      propensity_model_fn = stats::glm,
      model_fn = stats::glm
    ),
    class = "causatr_model_fn_default"
  )
})

test_that("SNM emits propensity_model_fn default warning when omitted", {
  d <- simulate_binary_continuous(n = 500, seed = 42)
  withr::local_options(causatr.suppress_default_warnings = FALSE)
  expect_warning(
    causat(
      d,
      outcome = "Y",
      treatment = "A",
      confounders = ~L,
      estimator = "snm"
    ),
    class = "causatr_propensity_fn_default"
  )
})

# --- build_blip_spec --------------------------------------------------------

test_that("build_blip_spec handles no-EM case", {
  em_info <- parse_effect_mod(~ L + M, "A")
  spec <- build_blip_spec(em_info, "A")
  expect_equal(spec$modifier_cols, character(0L))
  expect_equal(spec$param_names, "psi_intercept")
  expect_equal(spec$n_params, 1L)
})

test_that("build_blip_spec handles single modifier", {
  em_info <- parse_effect_mod(~ L + A:sex, "A")
  spec <- build_blip_spec(em_info, "A")
  expect_equal(spec$modifier_cols, "sex")
  expect_equal(spec$param_names, c("psi_intercept", "psi_sex"))
  expect_equal(spec$n_params, 2L)
})

test_that("build_blip_spec handles multiple modifiers", {
  em_info <- parse_effect_mod(~ L + A:sex + A:age, "A")
  spec <- build_blip_spec(em_info, "A")
  expect_equal(sort(spec$modifier_cols), c("age", "sex"))
  expect_equal(spec$n_params, 3L)
})
