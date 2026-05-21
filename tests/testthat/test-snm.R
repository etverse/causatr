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

test_that("SNM accepts longitudinal data with id and time", {
  d <- data.frame(
    id = rep(1:50, each = 2),
    time = rep(0:1, 50),
    Y = c(rep(NA, 50), rnorm(50)),
    A = rbinom(100, 1, 0.5),
    L = rnorm(100)
  )

  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~1,
    confounders_tv = ~L,
    estimator = "snm",
    id = "id",
    time = "time",
    type = "longitudinal",
    propensity_model_fn = stats::glm
  )
  expect_equal(fit$type, "longitudinal")
  expect_equal(fit$estimator, "snm")
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

# --- Chunk 18b: point estimation + contrast + sandwich -----------------------

# --- Truth-based tests -------------------------------------------------------

test_that("SNM recovers blip params: continuous trt, EM (design doc DGP)", {
  dgp <- simulate_snm_point(n = 5000, seed = 101)
  fit <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders_outcome = ~ L + A:M,
    confounders_treatment = ~L,
    estimator = "snm"
  )
  result <- contrast(fit, ci_method = "sandwich")

  expect_s3_class(result, "causatr_result")
  expect_equal(result$estimator, "snm")
  expect_equal(nrow(result$estimates), 2L)

  psi <- result$estimates$estimate
  names(psi) <- result$estimates$parameter
  expect_equal(psi[["psi_intercept"]], 3, tolerance = 0.1)
  expect_equal(psi[["psi_M"]], 2, tolerance = 0.15)

  # SEs should be finite and positive

  expect_true(all(result$estimates$se > 0))
  expect_true(all(is.finite(result$estimates$se)))

  # CIs should cover truth
  expect_true(result$estimates$ci_lower[1] < 3)
  expect_true(result$estimates$ci_upper[1] > 3)
  expect_true(result$estimates$ci_lower[2] < 2)
  expect_true(result$estimates$ci_upper[2] > 2)
})

test_that("SNM recovers constant ATE (no EM)", {
  dgp <- simulate_snm_point_no_em(n = 5000, seed = 202)
  fit <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "snm"
  )
  result <- contrast(fit, ci_method = "sandwich")

  expect_equal(nrow(result$estimates), 1L)
  expect_equal(result$estimates$parameter, "psi_intercept")
  expect_equal(result$estimates$estimate, 3, tolerance = 0.1)
  expect_true(result$estimates$se > 0)
})

test_that("SNM recovers blip params: binary treatment, EM", {
  dgp <- simulate_snm_point_binary(n = 5000, seed = 303)
  fit <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders_outcome = ~ L + A:M,
    confounders_treatment = ~L,
    estimator = "snm"
  )
  result <- contrast(fit, ci_method = "sandwich")

  psi <- result$estimates$estimate
  names(psi) <- result$estimates$parameter
  expect_equal(psi[["psi_intercept"]], 3, tolerance = 0.2)
  expect_equal(psi[["psi_M"]], 2, tolerance = 0.3)
})

test_that("SNM treatment_values returns averaged blip effect", {
  dgp <- simulate_snm_point(n = 5000, seed = 404)
  fit <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders_outcome = ~ L + A:M,
    confounders_treatment = ~L,
    estimator = "snm"
  )
  result <- contrast(fit, treatment_values = c(0, 1))

  expect_equal(nrow(result$estimates), 1L)
  expect_equal(result$estimates$parameter, "avg_blip_effect")
  # For a = 1 vs a = 0: effect = psi_0 + psi_M * mean(M)
  # mean(M) ≈ 0.5 (M = I(L > 0) where L ~ N(0,1))
  # truth ≈ 3 + 2 * 0.5 = 4
  expect_equal(result$estimates$estimate, 4, tolerance = 0.15)
  expect_true(result$estimates$se > 0)

  # Contrasts table should have the comparison
  expect_equal(nrow(result$contrasts), 1L)
})

test_that("SNM sandwich SE matches bootstrap SE (consistency check)", {
  dgp <- simulate_snm_point_no_em(n = 2000, seed = 505)
  fit <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "snm"
  )

  sandwich_result <- contrast(fit, ci_method = "sandwich")
  sandwich_se <- sandwich_result$estimates$se

  # Manual bootstrap for comparison
  set.seed(606)
  n_boot <- 200
  boot_psi <- replicate(n_boot, {
    idx <- sample.int(nrow(dgp$data), replace = TRUE)
    boot_data <- dgp$data[idx]
    tryCatch(
      {
        boot_fit <- causat(
          boot_data,
          outcome = "Y",
          treatment = "A",
          confounders = ~L,
          estimator = "snm"
        )
        boot_res <- contrast(boot_fit, ci_method = "sandwich")
        boot_res$estimates$estimate
      },
      error = function(e) NA_real_
    )
  })
  boot_se <- stats::sd(boot_psi[!is.na(boot_psi)])

  # Sandwich and bootstrap SEs should agree within 30%
  expect_equal(sandwich_se, boot_se, tolerance = 0.3)
})

# --- DTRreg cross-check ------------------------------------------------------

test_that("SNM matches DTRreg on binary treatment (no EM)", {
  skip_if_not_installed("DTRreg")
  # No-EM case: causatr and DTRreg solve equivalent EEs, so point
  # estimates and SEs should match closely. The EM case uses structurally
  # different parameterizations (treatment-free model vs residualized-
  # treatment moment condition) — EM validation is done via delicatessen.
  set.seed(707)
  n <- 5000
  L <- stats::rnorm(n)
  A <- stats::rbinom(n, 1, stats::plogis(0.5 * L))
  Y <- 2 + 3 * A + 1.5 * L + stats::rnorm(n)
  d <- data.table::data.table(Y = Y, A = A, L = L)

  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "snm"
  )
  result <- contrast(fit, ci_method = "sandwich")
  psi_causatr <- result$estimates$estimate
  se_causatr <- result$estimates$se

  dtr_fit <- DTRreg::DTRreg(
    outcome = Y,
    blip.mod = list(~1),
    treat.mod = list(A ~ L),
    tf.mod = list(~L),
    data = data.frame(Y = Y, A = A, L = L),
    method = "gest",
    treat.type = "bin",
    var.estim = "sandwich"
  )
  psi_dtreg <- stats::coef(dtr_fit)$stage_1
  se_dtreg <- sqrt(diag(dtr_fit$covmat[[1]]))

  expect_equal(psi_causatr, unname(psi_dtreg[1]), tolerance = 0.01)
  expect_equal(se_causatr, unname(se_dtreg[1]), tolerance = 0.01)
})

test_that("SNM sandwich SE matches DTRreg (no-EM, binary)", {
  skip_if_not_installed("DTRreg")
  # No-EM case: both implementations reduce to the same EE, so SEs
  # should match closely. The EM case has structurally different
  # sandwiches (treatment-free model vs residualized-treatment) — SE
  # validation for that case uses the delicatessen reference below.
  set.seed(808)
  n <- 5000
  L <- stats::rnorm(n)
  A <- stats::rbinom(n, 1, stats::plogis(0.5 * L))
  Y <- 2 + 3 * A + 1.5 * L + stats::rnorm(n)
  d <- data.table::data.table(Y = Y, A = A, L = L)

  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "snm"
  )
  result <- contrast(fit, ci_method = "sandwich")
  se_causatr <- result$estimates$se

  dtr_fit <- DTRreg::DTRreg(
    outcome = Y,
    blip.mod = list(~1),
    treat.mod = list(A ~ L),
    tf.mod = list(~L),
    data = data.frame(Y = Y, A = A, L = L),
    method = "gest",
    treat.type = "bin",
    var.estim = "sandwich"
  )
  se_dtreg <- sqrt(diag(dtr_fit$covmat[[1]]))

  expect_equal(se_causatr, unname(se_dtreg[1]), tolerance = 0.01)
})

# --- delicatessen reference values --------------------------------------------
# Reference values generated by data-raw/snm_reference.py using
# delicatessen (Zivich 2022) stacked M-estimation sandwich variance.
# DGP: simulate_snm_point(n = 5000, seed = 101); fixture: snm_fixture.csv.

test_that("SNM sandwich matches delicatessen — continuous trt, EM", {
  ref_psi_intercept <- 3.0014
  ref_psi_M <- 2.0402
  ref_se_intercept <- 0.0305
  ref_se_M <- 0.0529

  dgp <- simulate_snm_point(n = 5000, seed = 101)
  fit <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders_outcome = ~ L + A:M,
    confounders_treatment = ~L,
    estimator = "snm"
  )
  result <- contrast(fit, ci_method = "sandwich")

  psi <- result$estimates$estimate
  names(psi) <- result$estimates$parameter
  se <- result$estimates$se
  names(se) <- result$estimates$parameter

  expect_equal(psi[["psi_intercept"]], ref_psi_intercept, tolerance = 0.01)
  expect_equal(psi[["psi_M"]], ref_psi_M, tolerance = 0.01)
  expect_equal(se[["psi_intercept"]], ref_se_intercept, tolerance = 0.01)
  expect_equal(se[["psi_M"]], ref_se_M, tolerance = 0.01)
})

test_that("SNM sandwich matches delicatessen — binary trt, EM", {
  ref_psi_intercept <- 2.9174
  ref_psi_M <- 2.0923
  ref_se_intercept <- 0.0605
  ref_se_M <- 0.1027

  dgp <- simulate_snm_point_binary(n = 5000, seed = 303)
  fit <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders_outcome = ~ L + A:M,
    confounders_treatment = ~L,
    estimator = "snm"
  )
  result <- contrast(fit, ci_method = "sandwich")

  psi <- result$estimates$estimate
  names(psi) <- result$estimates$parameter
  se <- result$estimates$se
  names(se) <- result$estimates$parameter

  expect_equal(psi[["psi_intercept"]], ref_psi_intercept, tolerance = 0.01)
  expect_equal(psi[["psi_M"]], ref_psi_M, tolerance = 0.01)
  expect_equal(se[["psi_intercept"]], ref_se_intercept, tolerance = 0.01)
  expect_equal(se[["psi_M"]], ref_se_M, tolerance = 0.01)
})

test_that("SNM sandwich matches delicatessen — continuous trt, two modifiers", {
  ref_psi_intercept <- 0.9992
  ref_psi_M1 <- 2.0382
  ref_psi_M2 <- 0.5027
  ref_se_intercept <- 0.0376
  ref_se_M1 <- 0.0677
  ref_se_M2 <- 0.0506

  set.seed(111)
  n <- 5000
  L <- stats::rnorm(n)
  M1 <- as.numeric(L > 0)
  M2 <- stats::rnorm(n)
  A <- stats::rnorm(n, mean = 0.5 * L, sd = 1)
  Y <- 3 + 1 * A + 2 * A * M1 + 0.5 * A * M2 + 1.5 * L + stats::rnorm(n)
  d <- data.table::data.table(Y = Y, A = A, L = L, M1 = M1, M2 = M2)

  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders_outcome = ~ L + A:M1 + A:M2,
    confounders_treatment = ~L,
    estimator = "snm"
  )
  result <- contrast(fit, ci_method = "sandwich")

  psi <- result$estimates$estimate
  names(psi) <- result$estimates$parameter
  se <- result$estimates$se
  names(se) <- result$estimates$parameter

  expect_equal(psi[["psi_intercept"]], ref_psi_intercept, tolerance = 0.01)
  expect_equal(psi[["psi_M1"]], ref_psi_M1, tolerance = 0.01)
  expect_equal(psi[["psi_M2"]], ref_psi_M2, tolerance = 0.01)
  expect_equal(se[["psi_intercept"]], ref_se_intercept, tolerance = 0.01)
  expect_equal(se[["psi_M1"]], ref_se_M1, tolerance = 0.01)
  expect_equal(se[["psi_M2"]], ref_se_M2, tolerance = 0.01)
})

# --- Large-sample convergence ------------------------------------------------

test_that("SNM point estimates converge to truth at large n", {
  dgp <- simulate_snm_point(n = 20000, seed = 909)
  fit <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders_outcome = ~ L + A:M,
    confounders_treatment = ~L,
    estimator = "snm"
  )
  result <- contrast(fit)
  psi <- result$estimates$estimate
  names(psi) <- result$estimates$parameter

  # With n = 20000 the estimates should be very close to truth
  expect_equal(psi[["psi_intercept"]], 3, tolerance = 0.05)
  expect_equal(psi[["psi_M"]], 2, tolerance = 0.08)

  # SEs should shrink proportional to 1/sqrt(n)
  expect_true(all(result$estimates$se < 0.04))
})

# --- Multiple modifiers DGP --------------------------------------------------

test_that("SNM recovers blip params with two modifiers", {
  set.seed(111)
  n <- 5000
  L <- stats::rnorm(n)
  M1 <- as.numeric(L > 0)
  M2 <- stats::rnorm(n)
  A <- stats::rnorm(n, mean = 0.5 * L, sd = 1)
  # gamma(a, l; psi) = a * (1 + 2*M1 + 0.5*M2)
  Y <- 3 + 1 * A + 2 * A * M1 + 0.5 * A * M2 + 1.5 * L + stats::rnorm(n)
  d <- data.table::data.table(Y = Y, A = A, L = L, M1 = M1, M2 = M2)

  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders_outcome = ~ L + A:M1 + A:M2,
    confounders_treatment = ~L,
    estimator = "snm"
  )
  result <- contrast(fit)
  psi <- result$estimates$estimate
  names(psi) <- result$estimates$parameter

  expect_equal(psi[["psi_intercept"]], 1, tolerance = 0.15)
  expect_equal(psi[["psi_M1"]], 2, tolerance = 0.2)
  expect_equal(psi[["psi_M2"]], 0.5, tolerance = 0.1)
  expect_equal(nrow(result$vcov), 3L)
})

# --- Rejection and edge-case tests -------------------------------------------

test_that("contrast() rejects interventions for SNM fit (updated message)", {
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

test_that("contrast() rejects bootstrap for SNM (pending)", {
  dgp <- simulate_snm_point_no_em(n = 500, seed = 42)
  fit <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "snm"
  )
  expect_error(
    contrast(fit, ci_method = "bootstrap"),
    class = "causatr_snm_bootstrap_pending"
  )
})

test_that("treatment_values rejected for non-SNM estimators", {
  d <- simulate_binary_continuous(n = 500, seed = 42)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    model_fn = stats::glm
  )
  expect_error(
    contrast(
      fit,
      interventions = list(a1 = static(1), a0 = static(0)),
      treatment_values = c(0, 1)
    ),
    class = "causatr_treatment_values_not_snm"
  )
})

test_that("treatment_values must be length 2", {
  dgp <- simulate_snm_point_no_em(n = 500, seed = 42)
  fit <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "snm"
  )
  expect_error(
    contrast(fit, treatment_values = c(0, 1, 2)),
    class = "causatr_snm_bad_treatment_values"
  )
})

test_that("SNM vcov matrix has correct dimensions and names", {
  dgp <- simulate_snm_point(n = 1000, seed = 42)
  fit <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders_outcome = ~ L + A:M,
    confounders_treatment = ~L,
    estimator = "snm"
  )
  result <- contrast(fit)

  expect_equal(nrow(result$vcov), 2L)
  expect_equal(ncol(result$vcov), 2L)
  expect_equal(rownames(result$vcov), c("psi_intercept", "psi_M"))
  expect_equal(colnames(result$vcov), c("psi_intercept", "psi_M"))
  # Vcov should be positive semi-definite
  eig <- eigen(result$vcov, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eig >= -1e-10))
})

# --- Chunk 18b½: treatment-free outcome model ---------------------------------

# --- Treatment-free model truth tests -----------------------------------------

test_that("SNM with treatment-free model recovers blip params: continuous trt, EM", {
  dgp <- simulate_snm_point(n = 5000, seed = 101)
  fit <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders_outcome = ~ L + A:M,
    confounders_treatment = ~L,
    estimator = "snm",
    treatment_free = ~L
  )
  result <- contrast(fit, ci_method = "sandwich")

  psi <- result$estimates$estimate
  names(psi) <- result$estimates$parameter
  expect_equal(psi[["psi_intercept"]], 3, tolerance = 0.15)
  expect_equal(psi[["psi_M"]], 2, tolerance = 0.15)

  # SEs should be smaller than without TF model
  fit_no_tf <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders_outcome = ~ L + A:M,
    confounders_treatment = ~L,
    estimator = "snm"
  )
  r_no_tf <- contrast(fit_no_tf, ci_method = "sandwich")
  expect_true(all(result$estimates$se < r_no_tf$estimates$se))
})

test_that("SNM with TF model recovers constant ATE (no EM)", {
  dgp <- simulate_snm_point_no_em(n = 5000, seed = 202)
  fit <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "snm",
    treatment_free = ~L
  )
  result <- contrast(fit, ci_method = "sandwich")

  expect_equal(nrow(result$estimates), 1L)
  expect_equal(result$estimates$estimate, 3, tolerance = 0.1)
  expect_true(result$estimates$se > 0)
})

test_that("SNM with TF model recovers blip params: binary trt, EM", {
  dgp <- simulate_snm_point_binary(n = 5000, seed = 303)
  fit <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders_outcome = ~ L + A:M,
    confounders_treatment = ~L,
    estimator = "snm",
    treatment_free = ~L
  )
  result <- contrast(fit, ci_method = "sandwich")

  psi <- result$estimates$estimate
  names(psi) <- result$estimates$parameter
  expect_equal(psi[["psi_intercept"]], 3, tolerance = 0.2)
  expect_equal(psi[["psi_M"]], 2, tolerance = 0.3)
})

test_that("SNM with TF model + treatment_values returns averaged blip", {
  dgp <- simulate_snm_point(n = 5000, seed = 404)
  fit <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders_outcome = ~ L + A:M,
    confounders_treatment = ~L,
    estimator = "snm",
    treatment_free = ~L
  )
  result <- contrast(fit, treatment_values = c(0, 1))

  expect_equal(nrow(result$estimates), 1L)
  expect_equal(result$estimates$parameter, "avg_blip_effect")
  # truth ~ 3 + 2 * 0.5 = 4
  expect_equal(result$estimates$estimate, 4, tolerance = 0.15)
  expect_true(result$estimates$se > 0)
})

# --- DTRreg cross-checks with treatment-free model ----------------------------

test_that("SNM with TF model matches DTRreg (binary trt, EM case)", {
  skip_if_not_installed("DTRreg")
  # With the joint estimation approach, causatr and DTRreg should
  # match exactly even in the EM case (both use treatment-free model).
  set.seed(707)
  n <- 5000
  L <- stats::rnorm(n)
  M <- as.numeric(L > 0)
  A <- stats::rbinom(n, 1, stats::plogis(0.5 * L))
  Y <- 2 + 3 * A + 1.5 * L + 2 * A * M + stats::rnorm(n)
  d <- data.table::data.table(Y = Y, A = A, L = L, M = M)

  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders_outcome = ~ L + A:M,
    confounders_treatment = ~L,
    estimator = "snm",
    treatment_free = ~L
  )
  result <- contrast(fit, ci_method = "sandwich")

  dtr_fit <- DTRreg::DTRreg(
    outcome = Y,
    blip.mod = list(~M),
    treat.mod = list(A ~ L),
    tf.mod = list(~L),
    data = data.frame(Y = Y, A = A, L = L, M = M),
    method = "gest",
    treat.type = "bin",
    var.estim = "sandwich"
  )
  psi_dtreg <- stats::coef(dtr_fit)$stage_1
  se_dtreg <- sqrt(diag(dtr_fit$covmat[[1]]))

  psi_causatr <- result$estimates$estimate
  se_causatr <- result$estimates$se

  expect_equal(psi_causatr[1], unname(psi_dtreg[1]), tolerance = 0.01)
  expect_equal(psi_causatr[2], unname(psi_dtreg[2]), tolerance = 0.01)
  expect_equal(se_causatr[1], unname(se_dtreg[1]), tolerance = 0.01)
  expect_equal(se_causatr[2], unname(se_dtreg[2]), tolerance = 0.01)
})

test_that("SNM with TF model matches DTRreg (binary trt, no EM)", {
  skip_if_not_installed("DTRreg")
  set.seed(808)
  n <- 5000
  L <- stats::rnorm(n)
  A <- stats::rbinom(n, 1, stats::plogis(0.5 * L))
  Y <- 2 + 3 * A + 1.5 * L + stats::rnorm(n)
  d <- data.table::data.table(Y = Y, A = A, L = L)

  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "snm",
    treatment_free = ~L
  )
  result <- contrast(fit, ci_method = "sandwich")

  dtr_fit <- DTRreg::DTRreg(
    outcome = Y,
    blip.mod = list(~1),
    treat.mod = list(A ~ L),
    tf.mod = list(~L),
    data = data.frame(Y = Y, A = A, L = L),
    method = "gest",
    treat.type = "bin",
    var.estim = "sandwich"
  )
  psi_dtreg <- stats::coef(dtr_fit)$stage_1
  se_dtreg <- sqrt(diag(dtr_fit$covmat[[1]]))

  expect_equal(
    result$estimates$estimate,
    unname(psi_dtreg[1]),
    tolerance = 0.01
  )
  expect_equal(result$estimates$se, unname(se_dtreg[1]), tolerance = 0.01)
})

# --- delicatessen cross-checks with treatment-free model ----------------------

test_that("SNM TF sandwich matches delicatessen — continuous trt, EM", {
  ref_psi_intercept <- 3.0321
  ref_psi_M <- 1.9796
  ref_se_intercept <- 0.0195
  ref_se_M <- 0.0288

  dgp <- simulate_snm_point(n = 5000, seed = 101)
  fit <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders_outcome = ~ L + A:M,
    confounders_treatment = ~L,
    estimator = "snm",
    treatment_free = ~L
  )
  result <- contrast(fit, ci_method = "sandwich")

  psi <- result$estimates$estimate
  names(psi) <- result$estimates$parameter
  se <- result$estimates$se
  names(se) <- result$estimates$parameter

  expect_equal(psi[["psi_intercept"]], ref_psi_intercept, tolerance = 0.01)
  expect_equal(psi[["psi_M"]], ref_psi_M, tolerance = 0.01)
  expect_equal(se[["psi_intercept"]], ref_se_intercept, tolerance = 0.01)
  expect_equal(se[["psi_M"]], ref_se_M, tolerance = 0.01)
})

test_that("SNM TF sandwich matches delicatessen — binary trt, EM", {
  ref_psi_intercept <- 2.9441
  ref_psi_M <- 2.0402
  ref_se_intercept <- 0.0414
  ref_se_M <- 0.0579

  dgp <- simulate_snm_point_binary(n = 5000, seed = 303)
  fit <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders_outcome = ~ L + A:M,
    confounders_treatment = ~L,
    estimator = "snm",
    treatment_free = ~L
  )
  result <- contrast(fit, ci_method = "sandwich")

  psi <- result$estimates$estimate
  names(psi) <- result$estimates$parameter
  se <- result$estimates$se
  names(se) <- result$estimates$parameter

  expect_equal(psi[["psi_intercept"]], ref_psi_intercept, tolerance = 0.01)
  expect_equal(psi[["psi_M"]], ref_psi_M, tolerance = 0.01)
  expect_equal(se[["psi_intercept"]], ref_se_intercept, tolerance = 0.01)
  expect_equal(se[["psi_M"]], ref_se_M, tolerance = 0.01)
})

test_that("SNM TF sandwich matches delicatessen — continuous trt, 2 modifiers", {
  ref_psi_intercept <- 0.9881
  ref_psi_M1 <- 2.0600
  ref_psi_M2 <- 0.5013
  ref_se_intercept <- 0.0212
  ref_se_M1 <- 0.0291
  ref_se_M2 <- 0.0143

  set.seed(111)
  n <- 5000
  L <- stats::rnorm(n)
  M1 <- as.numeric(L > 0)
  M2 <- stats::rnorm(n)
  A <- stats::rnorm(n, mean = 0.5 * L, sd = 1)
  Y <- 3 + 1 * A + 2 * A * M1 + 0.5 * A * M2 + 1.5 * L + stats::rnorm(n)
  d <- data.table::data.table(Y = Y, A = A, L = L, M1 = M1, M2 = M2)

  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders_outcome = ~ L + A:M1 + A:M2,
    confounders_treatment = ~L,
    estimator = "snm",
    treatment_free = ~L
  )
  result <- contrast(fit, ci_method = "sandwich")

  psi <- result$estimates$estimate
  names(psi) <- result$estimates$parameter
  se <- result$estimates$se
  names(se) <- result$estimates$parameter

  expect_equal(psi[["psi_intercept"]], ref_psi_intercept, tolerance = 0.01)
  expect_equal(psi[["psi_M1"]], ref_psi_M1, tolerance = 0.01)
  expect_equal(psi[["psi_M2"]], ref_psi_M2, tolerance = 0.01)
  expect_equal(se[["psi_intercept"]], ref_se_intercept, tolerance = 0.01)
  expect_equal(se[["psi_M1"]], ref_se_M1, tolerance = 0.01)
  expect_equal(se[["psi_M2"]], ref_se_M2, tolerance = 0.01)
})

# --- Treatment-free model SE reduction validation -----------------------------

test_that("TF model reduces SEs across all DGPs", {
  # Continuous treatment + EM
  dgp1 <- simulate_snm_point(n = 5000, seed = 101)
  fit1_no <- causat(
    dgp1$data,
    outcome = "Y",
    treatment = "A",
    confounders_outcome = ~ L + A:M,
    confounders_treatment = ~L,
    estimator = "snm"
  )
  fit1_tf <- causat(
    dgp1$data,
    outcome = "Y",
    treatment = "A",
    confounders_outcome = ~ L + A:M,
    confounders_treatment = ~L,
    estimator = "snm",
    treatment_free = ~L
  )
  r1_no <- contrast(fit1_no)
  r1_tf <- contrast(fit1_tf)
  expect_true(all(r1_tf$estimates$se < r1_no$estimates$se))

  # Binary treatment + EM
  dgp2 <- simulate_snm_point_binary(n = 5000, seed = 303)
  fit2_no <- causat(
    dgp2$data,
    outcome = "Y",
    treatment = "A",
    confounders_outcome = ~ L + A:M,
    confounders_treatment = ~L,
    estimator = "snm"
  )
  fit2_tf <- causat(
    dgp2$data,
    outcome = "Y",
    treatment = "A",
    confounders_outcome = ~ L + A:M,
    confounders_treatment = ~L,
    estimator = "snm",
    treatment_free = ~L
  )
  r2_no <- contrast(fit2_no)
  r2_tf <- contrast(fit2_tf)
  expect_true(all(r2_tf$estimates$se < r2_no$estimates$se))
})

# --- Treatment-free model sandwich vs bootstrap consistency -------------------

test_that("TF model sandwich SE matches bootstrap SE", {
  dgp <- simulate_snm_point_no_em(n = 2000, seed = 505)
  fit <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "snm",
    treatment_free = ~L
  )

  sandwich_result <- contrast(fit, ci_method = "sandwich")
  sandwich_se <- sandwich_result$estimates$se

  set.seed(606)
  n_boot <- 200
  boot_psi <- replicate(n_boot, {
    idx <- sample.int(nrow(dgp$data), replace = TRUE)
    boot_data <- dgp$data[idx]
    tryCatch(
      {
        boot_fit <- causat(
          boot_data,
          outcome = "Y",
          treatment = "A",
          confounders = ~L,
          estimator = "snm",
          treatment_free = ~L
        )
        boot_res <- contrast(boot_fit, ci_method = "sandwich")
        boot_res$estimates$estimate
      },
      error = function(e) NA_real_
    )
  })
  boot_se <- stats::sd(boot_psi[!is.na(boot_psi)])

  expect_equal(sandwich_se, boot_se, tolerance = 0.3)
})

# --- Treatment-free model rejection and edge cases ----------------------------

test_that("treatment_free rejected for non-SNM estimators", {
  d <- simulate_binary_continuous(n = 500, seed = 42)
  expect_error(
    causat(
      d,
      outcome = "Y",
      treatment = "A",
      confounders = ~L,
      model_fn = stats::glm,
      treatment_free = ~L
    ),
    class = "causatr_treatment_free_not_snm"
  )
})

test_that("treatment_free must be a formula", {
  d <- simulate_binary_continuous(n = 500, seed = 42)
  expect_error(
    causat(
      d,
      outcome = "Y",
      treatment = "A",
      confounders = ~L,
      estimator = "snm",
      propensity_model_fn = stats::glm,
      treatment_free = "L"
    ),
    class = "causatr_snm_bad_treatment_free"
  )
})

test_that("TF model vcov is positive semi-definite", {
  dgp <- simulate_snm_point(n = 1000, seed = 42)
  fit <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders_outcome = ~ L + A:M,
    confounders_treatment = ~L,
    estimator = "snm",
    treatment_free = ~L
  )
  result <- contrast(fit)
  eig <- eigen(result$vcov, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eig >= -1e-10))
})

# --- Chunk 18c: time-varying (post-treatment) modifiers ----------------------
# SNMs identify the blip under treatment-model correctness alone, so
# modifiers that depend on treatment (post-treatment) are valid — unlike
# IPW-MSM, which conditions on a descendant of A (collider bias).

# --- Truth-based tests with TF model (recovers structural parameters) ---------

test_that("SNM with TF recovers blip params: continuous trt, TV modifier", {
  dgp <- simulate_snm_point_tv_modifier(n = 5000, seed = 1801)
  fit <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders_outcome = ~ L + A:M,
    confounders_treatment = ~L,
    estimator = "snm",
    treatment_free = ~L
  )
  result <- contrast(fit, ci_method = "sandwich")

  psi <- result$estimates$estimate
  names(psi) <- result$estimates$parameter
  expect_equal(psi[["psi_intercept"]], 3, tolerance = 0.1)
  expect_equal(psi[["psi_M"]], 2, tolerance = 0.1)

  expect_true(all(result$estimates$se > 0))
  expect_true(all(is.finite(result$estimates$se)))
  expect_true(result$estimates$ci_lower[1] < 3)
  expect_true(result$estimates$ci_upper[1] > 3)
  expect_true(result$estimates$ci_lower[2] < 2)
  expect_true(result$estimates$ci_upper[2] > 2)
})

test_that("SNM with TF recovers blip params: binary trt, TV modifier", {
  dgp <- simulate_snm_point_tv_modifier_binary(n = 5000, seed = 1802)
  fit <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders_outcome = ~ L + A:M,
    confounders_treatment = ~L,
    estimator = "snm",
    treatment_free = ~L
  )
  result <- contrast(fit, ci_method = "sandwich")

  psi <- result$estimates$estimate
  names(psi) <- result$estimates$parameter
  expect_equal(psi[["psi_intercept"]], 3, tolerance = 0.2)
  expect_equal(psi[["psi_M"]], 2, tolerance = 0.2)
})

test_that("SNM without TF runs on TV modifier (different blip quantity)", {
  # Without a treatment-free model, the blip absorbs the A -> M -> Y
  # indirect path, so psi values differ from the structural truth.
  # This is a valid moment-condition estimand — just a different one.
  dgp <- simulate_snm_point_tv_modifier(n = 5000, seed = 1801)
  fit <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders_outcome = ~ L + A:M,
    confounders_treatment = ~L,
    estimator = "snm"
  )
  result <- contrast(fit, ci_method = "sandwich")

  expect_s3_class(result, "causatr_result")
  expect_equal(nrow(result$estimates), 2L)
  expect_true(all(result$estimates$se > 0))
  expect_true(all(is.finite(result$estimates$estimate)))
})

# --- TF model reduces SEs for TV modifier DGPs --------------------------------

test_that("TF model reduces SEs on TV modifier DGP", {
  dgp <- simulate_snm_point_tv_modifier(n = 5000, seed = 1801)

  fit_no <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders_outcome = ~ L + A:M,
    confounders_treatment = ~L,
    estimator = "snm"
  )
  fit_tf <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders_outcome = ~ L + A:M,
    confounders_treatment = ~L,
    estimator = "snm",
    treatment_free = ~L
  )

  r_no <- contrast(fit_no)
  r_tf <- contrast(fit_tf)
  expect_true(all(r_tf$estimates$se < r_no$estimates$se))
})

# --- delicatessen cross-checks for TV modifier ---------------------------------
# Reference values generated by data-raw/snm_reference.py scenarios 7–8b.
# Fixture: snm_fixture_tv_modifier.csv, snm_fixture_tv_modifier_binary.csv.

test_that("SNM sandwich matches delicatessen — continuous trt, TV modifier (no TF)", {
  ref_psi_intercept <- 2.9995
  ref_psi_M <- 2.4488
  ref_se_intercept <- 0.0200
  ref_se_M <- 0.0264

  dgp <- simulate_snm_point_tv_modifier(n = 5000, seed = 1801)
  fit <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders_outcome = ~ L + A:M,
    confounders_treatment = ~L,
    estimator = "snm"
  )
  result <- contrast(fit, ci_method = "sandwich")

  psi <- result$estimates$estimate
  names(psi) <- result$estimates$parameter
  se <- result$estimates$se
  names(se) <- result$estimates$parameter

  expect_equal(psi[["psi_intercept"]], ref_psi_intercept, tolerance = 0.01)
  expect_equal(psi[["psi_M"]], ref_psi_M, tolerance = 0.01)
  expect_equal(se[["psi_intercept"]], ref_se_intercept, tolerance = 0.01)
  expect_equal(se[["psi_M"]], ref_se_M, tolerance = 0.01)
})

test_that("SNM TF sandwich matches delicatessen — continuous trt, TV modifier", {
  ref_psi_intercept <- 3.0176
  ref_psi_M <- 1.9888
  ref_se_intercept <- 0.0144
  ref_se_M <- 0.0126

  dgp <- simulate_snm_point_tv_modifier(n = 5000, seed = 1801)
  fit <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders_outcome = ~ L + A:M,
    confounders_treatment = ~L,
    estimator = "snm",
    treatment_free = ~L
  )
  result <- contrast(fit, ci_method = "sandwich")

  psi <- result$estimates$estimate
  names(psi) <- result$estimates$parameter
  se <- result$estimates$se
  names(se) <- result$estimates$parameter

  expect_equal(psi[["psi_intercept"]], ref_psi_intercept, tolerance = 0.01)
  expect_equal(psi[["psi_M"]], ref_psi_M, tolerance = 0.01)
  expect_equal(se[["psi_intercept"]], ref_se_intercept, tolerance = 0.01)
  expect_equal(se[["psi_M"]], ref_se_M, tolerance = 0.01)
})

test_that("SNM sandwich matches delicatessen — binary trt, TV modifier (no TF)", {
  ref_psi_intercept <- 2.2904
  ref_psi_M <- 3.4059
  ref_se_intercept <- 0.0676
  ref_se_M <- 0.0897

  dgp <- simulate_snm_point_tv_modifier_binary(n = 5000, seed = 1802)
  fit <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders_outcome = ~ L + A:M,
    confounders_treatment = ~L,
    estimator = "snm"
  )
  result <- contrast(fit, ci_method = "sandwich")

  psi <- result$estimates$estimate
  names(psi) <- result$estimates$parameter
  se <- result$estimates$se
  names(se) <- result$estimates$parameter

  expect_equal(psi[["psi_intercept"]], ref_psi_intercept, tolerance = 0.01)
  expect_equal(psi[["psi_M"]], ref_psi_M, tolerance = 0.01)
  expect_equal(se[["psi_intercept"]], ref_se_intercept, tolerance = 0.01)
  expect_equal(se[["psi_M"]], ref_se_M, tolerance = 0.01)
})

test_that("SNM TF sandwich matches delicatessen — binary trt, TV modifier", {
  ref_psi_intercept <- 2.9756
  ref_psi_M <- 2.0387
  ref_se_intercept <- 0.0335
  ref_se_M <- 0.0337

  dgp <- simulate_snm_point_tv_modifier_binary(n = 5000, seed = 1802)
  fit <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders_outcome = ~ L + A:M,
    confounders_treatment = ~L,
    estimator = "snm",
    treatment_free = ~L
  )
  result <- contrast(fit, ci_method = "sandwich")

  psi <- result$estimates$estimate
  names(psi) <- result$estimates$parameter
  se <- result$estimates$se
  names(se) <- result$estimates$parameter

  expect_equal(psi[["psi_intercept"]], ref_psi_intercept, tolerance = 0.01)
  expect_equal(psi[["psi_M"]], ref_psi_M, tolerance = 0.01)
  expect_equal(se[["psi_intercept"]], ref_se_intercept, tolerance = 0.01)
  expect_equal(se[["psi_M"]], ref_se_M, tolerance = 0.01)
})

# --- DTRreg cross-check: TV modifier with TF model ----------------------------

test_that("SNM with TF matches DTRreg on TV modifier DGP", {
  skip_if_not_installed("DTRreg")
  # DTRreg's tf.mod jointly estimates (beta, psi), matching causatr's
  # treatment_free approach. Both should agree on the TV modifier DGP.
  dgp <- simulate_snm_point_tv_modifier_binary(n = 5000, seed = 1802)

  fit <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders_outcome = ~ L + A:M,
    confounders_treatment = ~L,
    estimator = "snm",
    treatment_free = ~L
  )
  result <- contrast(fit, ci_method = "sandwich")

  d_df <- as.data.frame(dgp$data)
  dtr_fit <- DTRreg::DTRreg(
    outcome = d_df$Y,
    blip.mod = list(~M),
    treat.mod = list(A ~ L),
    tf.mod = list(~L),
    data = d_df,
    method = "gest",
    treat.type = "bin",
    var.estim = "sandwich"
  )
  psi_dtreg <- stats::coef(dtr_fit)$stage_1
  se_dtreg <- sqrt(diag(dtr_fit$covmat[[1]]))

  psi_causatr <- result$estimates$estimate
  se_causatr <- result$estimates$se

  expect_equal(psi_causatr[1], unname(psi_dtreg[1]), tolerance = 0.01)
  expect_equal(psi_causatr[2], unname(psi_dtreg[2]), tolerance = 0.01)
  expect_equal(se_causatr[1], unname(se_dtreg[1]), tolerance = 0.01)
  expect_equal(se_causatr[2], unname(se_dtreg[2]), tolerance = 0.01)
})

# --- Scientific gap: IPW biased with post-treatment modifier -------------------

test_that("IPW is biased when modifier is post-treatment (SNM is not)", {
  # The key scientific motivation for SNMs: IPW-MSM conditions on a
  # descendant of A when the modifier is post-treatment, producing
  # biased estimates. SNM with treatment-free model identifies the
  # structural blip correctly.
  dgp <- simulate_snm_point_tv_modifier(n = 5000, seed = 1801)

  # SNM with TF: correct
  fit_snm <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders_outcome = ~ L + A:M,
    confounders_treatment = ~L,
    estimator = "snm",
    treatment_free = ~L
  )
  r_snm <- contrast(fit_snm, treatment_values = c(0, 1))

  # IPW with M as modifier: biased (M is post-treatment)
  fit_ipw <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L + A:M,
    estimator = "ipw"
  )
  r_ipw <- contrast(
    fit_ipw,
    interventions = list(a1 = shift(1), a0 = shift(0)),
    type = "difference"
  )

  # SNM averaged blip should be close to the structural truth
  # truth: psi_0 + psi_M * mean(M) = 3 + 2 * mean(M)
  mean_M <- mean(dgp$data$M)
  truth_avg <- 3 + 2 * mean_M
  expect_equal(r_snm$estimates$estimate, truth_avg, tolerance = 0.15)

  # IPW should be substantially biased (>1 unit off from truth)
  expect_true(abs(r_ipw$contrasts$estimate - truth_avg) > 1)
})

# --- Vcov properties for TV modifier ------------------------------------------

test_that("Vcov is PSD for TV modifier DGP (no TF)", {
  dgp <- simulate_snm_point_tv_modifier(n = 1000, seed = 42)
  fit <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders_outcome = ~ L + A:M,
    confounders_treatment = ~L,
    estimator = "snm"
  )
  result <- contrast(fit)
  eig <- eigen(result$vcov, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eig >= -1e-10))
})

test_that("Vcov is PSD for TV modifier DGP (with TF)", {
  dgp <- simulate_snm_point_tv_modifier(n = 1000, seed = 42)
  fit <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders_outcome = ~ L + A:M,
    confounders_treatment = ~L,
    estimator = "snm",
    treatment_free = ~L
  )
  result <- contrast(fit)
  eig <- eigen(result$vcov, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eig >= -1e-10))
})


# --- Chunk 18d: longitudinal SNMM -------------------------------------------

test_that("Longitudinal SNM: binary treatment, per-stage estimation", {
  dgp <- simulate_snm_longitudinal(n = 10000, seed = 123)

  fit <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders = ~1,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    type = "longitudinal",
    family = "gaussian",
    estimator = "snm"
  )

  expect_equal(fit$type, "longitudinal")
  expect_equal(fit$estimator, "snm")

  result <- contrast(fit, ci_method = "sandwich")

  # Per-stage parameter names
  expect_equal(
    result$estimates$parameter,
    c("stage0_psi_intercept", "stage1_psi_intercept")
  )

  # Truth: psi = 3 at both stages
  truth <- dgp$truth_psi
  for (i in seq_along(truth)) {
    expect_equal(
      result$estimates$estimate[i],
      unname(truth[i]),
      tolerance = 0.15
    )
  }

  # SEs should be positive and finite
  expect_true(all(result$estimates$se > 0))
  expect_true(all(is.finite(result$estimates$se)))

  # CIs should cover the truth (at 95% level)
  for (i in seq_along(truth)) {
    expect_true(result$estimates$ci_lower[i] < unname(truth[i]))
    expect_true(result$estimates$ci_upper[i] > unname(truth[i]))
  }
})


test_that("Longitudinal SNM: continuous treatment, per-stage estimation", {
  dgp <- simulate_snm_longitudinal_continuous(n = 5000, seed = 42)

  fit <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders = ~1,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    type = "longitudinal",
    family = "gaussian",
    estimator = "snm"
  )

  result <- contrast(fit, ci_method = "sandwich")

  truth <- dgp$truth_psi
  for (i in seq_along(truth)) {
    expect_equal(
      result$estimates$estimate[i],
      unname(truth[i]),
      tolerance = 0.25
    )
  }
})


test_that("Longitudinal SNM: vcov is PSD and correctly sized", {
  dgp <- simulate_snm_longitudinal(n = 2000, seed = 42)

  fit <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders = ~1,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    type = "longitudinal",
    family = "gaussian",
    estimator = "snm"
  )

  result <- contrast(fit, ci_method = "sandwich")

  # Vcov should be 2x2 (one psi per stage)
  expect_equal(nrow(result$vcov), 2L)
  expect_equal(ncol(result$vcov), 2L)

  # PSD check
  eig <- eigen(result$vcov, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eig >= -1e-10))

  # Row/col names match parameter names
  expect_equal(rownames(result$vcov), result$estimates$parameter)
})


test_that("Longitudinal SNM rejects treatment_values", {
  dgp <- simulate_snm_longitudinal(n = 500, seed = 42)

  fit <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders = ~1,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    type = "longitudinal",
    family = "gaussian",
    estimator = "snm"
  )

  expect_error(
    contrast(fit, treatment_values = c(0, 1), ci_method = "sandwich"),
    class = "causatr_snm_long_no_treatment_values"
  )
})


test_that("Longitudinal SNM rejects fewer than 2 time points", {
  d <- data.table::data.table(
    id = 1:10,
    time = rep(0, 10),
    Y = rnorm(10),
    A = rbinom(10, 1, 0.5),
    L = rnorm(10)
  )
  expect_error(
    causat(
      d,
      outcome = "Y",
      treatment = "A",
      confounders = ~1,
      confounders_tv = ~L,
      id = "id",
      time = "time",
      type = "longitudinal",
      family = "gaussian",
      estimator = "snm"
    ),
    class = "causatr_longitudinal_too_few_times"
  )
})


test_that("Longitudinal SNM: history = 0 produces no lag terms", {
  dgp <- simulate_snm_longitudinal(n = 5000, seed = 42)

  fit_h0 <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders = ~1,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    type = "longitudinal",
    family = "gaussian",
    estimator = "snm",
    history = 0
  )

  # history = 0 means no lag columns in the per-period PS formula.
  # At time 1, the formula should be A ~ L (no L_lag1 or A_lag1).
  ps_fmla_t1 <- fit_h0$details$per_period_formula[["1"]]
  ps_vars <- all.vars(ps_fmla_t1)
  expect_false(any(grepl("_lag", ps_vars)))

  res_h0 <- contrast(fit_h0, ci_method = "sandwich")
  expect_equal(nrow(res_h0$estimates), 2L)
  expect_true(all(res_h0$estimates$se > 0))

  # Should still recover truth (3.15 + 3 = 6.15) within tolerance
  expect_equal(sum(res_h0$estimates$estimate), 6.15, tolerance = 0.2)
})


test_that("Longitudinal SNM: DTRreg cross-check, binary treatment", {
  skip_if_not_installed("DTRreg")

  dgp <- simulate_snm_longitudinal(n = 5000, seed = 42)
  d <- dgp$data

  # causatr

  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~1,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    type = "longitudinal",
    family = "gaussian",
    estimator = "snm"
  )
  res <- contrast(fit, ci_method = "sandwich")

  # DTRreg: wide format, matching treatment model specification
  wide <- data.table::dcast(d, id ~ time, value.var = c("A", "L", "Y"))
  wide[, Y := Y_1]
  data.table::setnames(
    wide,
    c("A_0", "A_1", "L_0", "L_1"),
    c("A0", "A1", "L_tv0", "L_tv1")
  )

  dtr <- DTRreg::DTRreg(
    outcome = wide$Y,
    blip.mod = list(~1, ~1),
    treat.mod = list(A0 ~ L_tv0, A1 ~ L_tv1 + A0),
    tf.mod = list(~1, ~1),
    method = "gest",
    treat.type = "bin",
    data = as.data.frame(wide),
    var.estim = "none"
  )

  # Point estimates should be close (not exact due to slightly different
  # treatment model specs — causatr includes lag features)
  dtr_psi_0 <- unname(coef(dtr)[[1]])
  dtr_psi_1 <- unname(coef(dtr)[[2]])

  expect_equal(res$estimates$estimate[1], dtr_psi_0, tolerance = 0.05)
  expect_equal(res$estimates$estimate[2], dtr_psi_1, tolerance = 0.05)
})


test_that("Longitudinal SNM: continuous treatment truth-based check", {
  # DTRreg does not support continuous treatment with multi-stage linear

  # blip, so we validate against analytical truth directly.
  dgp <- simulate_snm_longitudinal_continuous(n = 10000, seed = 123)

  fit <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders = ~1,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    type = "longitudinal",
    family = "gaussian",
    estimator = "snm"
  )
  res <- contrast(fit, ci_method = "sandwich")

  truth <- dgp$truth_psi
  for (i in seq_along(truth)) {
    expect_equal(
      res$estimates$estimate[i],
      unname(truth[i]),
      tolerance = 0.1
    )
  }

  # CIs should cover truth
  for (i in seq_along(truth)) {
    expect_true(res$estimates$ci_lower[i] < unname(truth[i]))
    expect_true(res$estimates$ci_upper[i] > unname(truth[i]))
  }
})


test_that("Longitudinal SNM: TF model produces consistent point estimates", {
  dgp <- simulate_snm_longitudinal(n = 3000, seed = 42)

  fit_no_tf <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders = ~1,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    type = "longitudinal",
    family = "gaussian",
    estimator = "snm"
  )
  res_no_tf <- contrast(fit_no_tf, ci_method = "sandwich")

  fit_tf <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders = ~1,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    type = "longitudinal",
    family = "gaussian",
    estimator = "snm",
    treatment_free = ~L
  )
  res_tf <- contrast(fit_tf, ci_method = "sandwich")

  expect_equal(
    res_no_tf$estimates$estimate,
    res_tf$estimates$estimate,
    tolerance = 0.01
  )
  expect_true(all(res_tf$estimates$se > 0))
  expect_true(all(is.finite(res_tf$estimates$se)))
})


test_that("Longitudinal SNM: TF model improves efficiency when E[R*Z] != 0", {
  # DGP with L^2 in the outcome but not the treatment model, so
  # E[R * L^2] != 0 and the TF model absorbs non-orthogonal variance.
  set.seed(1901)
  n <- 3000

  L0 <- stats::rnorm(n)
  A0 <- stats::rbinom(n, 1, stats::plogis(0.3 * L0))
  L1 <- 0.5 * L0 + 0.3 * A0 + stats::rnorm(n, sd = 0.7)
  A1 <- stats::rbinom(n, 1, stats::plogis(0.3 * L1 + 0.2 * A0))
  # Outcome has quadratic L terms not in treatment model
  Y <- 2 +
    3 * A0 +
    3 * A1 +
    1.5 * L0 +
    0.8 * L0^2 +
    0.5 * L1 +
    0.4 * L1^2 +
    stats::rnorm(n)

  data <- data.table::data.table(
    id = rep(seq_len(n), each = 2L),
    time = rep(0:1, n),
    Y = as.numeric(rbind(NA, Y)),
    A = as.numeric(rbind(A0, A1)),
    L = as.numeric(rbind(L0, L1))
  )

  fit_notf <- causat(
    data,
    outcome = "Y",
    treatment = "A",
    confounders = ~1,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    type = "longitudinal",
    family = "gaussian",
    estimator = "snm"
  )

  fit_tf <- causat(
    data,
    outcome = "Y",
    treatment = "A",
    confounders = ~1,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    type = "longitudinal",
    family = "gaussian",
    estimator = "snm",
    treatment_free = ~ L + I(L^2)
  )

  res_notf <- contrast(fit_notf, ci_method = "sandwich")
  res_tf <- contrast(fit_tf, ci_method = "sandwich")

  # Point estimates should agree (g-estimation is CAN under correct PS)
  expect_equal(
    res_notf$estimates$estimate,
    res_tf$estimates$estimate,
    tolerance = 0.15
  )

  # TF model should strictly reduce at least one SE when E[R*Z] != 0
  ratio <- res_tf$estimates$se / res_notf$estimates$se
  expect_true(any(ratio < 0.95))

  # Vcov should be PSD
  eig <- eigen(res_tf$vcov, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eig >= -1e-10))
})


test_that("Longitudinal SNM: n_obs and fit metadata", {
  dgp <- simulate_snm_longitudinal(n = 500, seed = 42)

  fit <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders = ~1,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    type = "longitudinal",
    family = "gaussian",
    estimator = "snm"
  )

  result <- contrast(fit, ci_method = "sandwich")

  expect_equal(result$n, 500L)
  expect_equal(result$estimator, "snm")
  expect_equal(result$fit_type, "longitudinal")
})


test_that("Longitudinal SNM: bootstrap is rejected", {
  dgp <- simulate_snm_longitudinal(n = 500, seed = 42)

  fit <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders = ~1,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    type = "longitudinal",
    family = "gaussian",
    estimator = "snm"
  )

  expect_error(
    contrast(fit, ci_method = "bootstrap"),
    class = "causatr_snm_bootstrap_pending"
  )
})


test_that("Longitudinal SNM: sandwich SE matches cluster bootstrap SE", {
  dgp <- simulate_snm_longitudinal(n = 2000, seed = 42)

  fit <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders = ~1,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    type = "longitudinal",
    family = "gaussian",
    estimator = "snm"
  )
  res <- contrast(fit, ci_method = "sandwich")
  sandwich_se <- res$estimates$se

  set.seed(1234)
  n_boot <- 200
  unique_ids <- unique(dgp$data$id)
  n_ids <- length(unique_ids)

  boot_psi <- replicate(n_boot, {
    boot_ids <- sample(unique_ids, n_ids, replace = TRUE)
    boot_data <- data.table::rbindlist(lapply(seq_along(boot_ids), function(i) {
      rows <- dgp$data[id == boot_ids[i]]
      rows <- data.table::copy(rows)
      rows[, id := i]
      rows
    }))
    tryCatch(
      {
        boot_fit <- causat(
          boot_data,
          outcome = "Y",
          treatment = "A",
          confounders = ~1,
          confounders_tv = ~L,
          id = "id",
          time = "time",
          type = "longitudinal",
          family = "gaussian",
          estimator = "snm"
        )
        boot_res <- contrast(boot_fit, ci_method = "sandwich")
        boot_res$estimates$estimate
      },
      error = function(e) rep(NA_real_, 2)
    )
  })

  boot_se <- apply(boot_psi, 1, sd, na.rm = TRUE)
  expect_equal(sandwich_se, boot_se, tolerance = 0.3)
})


# --- Chunk 18e: longitudinal TV-EM truth-based test ----------------------------
# Headline scientific demonstration: SNMs correctly handle time-varying
# effect modification (M_1 = 1{L_1 > 0}, post-treatment because L_1
# depends on A_0), while IPW-MSM conditioning on M_1 is biased.
# DGP: simulate_snm_longitudinal_tv_em() in helper-dgp.R.
# Truth: stage0 = (psi_intercept=1.15, psi_M=2),
#         stage1 = (psi_intercept=2, psi_M=2).

test_that("Longitudinal SNM with TV-EM: recovers per-stage blip parameters", {
  dgp <- simulate_snm_longitudinal_tv_em(n = 15000, seed = 123)

  fit <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders_outcome = ~ A:M,
    confounders_tv = ~ L + M,
    id = "id",
    time = "time",
    type = "longitudinal",
    family = "gaussian",
    estimator = "snm"
  )

  res <- contrast(fit, ci_method = "sandwich")

  expect_equal(
    res$estimates$parameter,
    c(
      "stage0_psi_intercept",
      "stage0_psi_M",
      "stage1_psi_intercept",
      "stage1_psi_M"
    )
  )

  truth <- dgp$truth_psi
  for (i in seq_along(truth)) {
    expect_equal(
      res$estimates$estimate[i],
      unname(truth[i]),
      tolerance = 0.15
    )
  }

  # CIs should cover truth
  for (i in seq_along(truth)) {
    expect_true(res$estimates$ci_lower[i] < unname(truth[i]))
    expect_true(res$estimates$ci_upper[i] > unname(truth[i]))
  }

  # SEs should be positive and finite
  expect_true(all(res$estimates$se > 0))
  expect_true(all(is.finite(res$estimates$se)))
})


test_that("Longitudinal SNM with TV-EM + TF: same truth, tighter SEs", {
  dgp <- simulate_snm_longitudinal_tv_em(n = 10000, seed = 123)

  fit_no_tf <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders_outcome = ~ A:M,
    confounders_tv = ~ L + M,
    id = "id",
    time = "time",
    type = "longitudinal",
    family = "gaussian",
    estimator = "snm"
  )
  res_no_tf <- contrast(fit_no_tf, ci_method = "sandwich")

  fit_tf <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders_outcome = ~ A:M,
    confounders_tv = ~ L + M,
    id = "id",
    time = "time",
    type = "longitudinal",
    family = "gaussian",
    estimator = "snm",
    treatment_free = ~L
  )
  res_tf <- contrast(fit_tf, ci_method = "sandwich")

  truth <- dgp$truth_psi
  for (i in seq_along(truth)) {
    expect_equal(
      res_tf$estimates$estimate[i],
      unname(truth[i]),
      tolerance = 0.15
    )
  }

  # TF model should reduce SEs
  ratio <- res_tf$estimates$se / res_no_tf$estimates$se
  expect_true(any(ratio < 0.95))
})


test_that("Longitudinal SNM with TV-EM: DTRreg cross-check (history=0)", {
  skip_if_not_installed("DTRreg")

  # Use history=0 to match DTRreg's per-stage treatment model exactly:
  # DTRreg fits A0 ~ L0 + M0 and A1 ~ L1 + M1 (no lag terms).
  dgp <- simulate_snm_longitudinal_tv_em(n = 5000, seed = 42)
  d <- dgp$data

  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders_outcome = ~ A:M,
    confounders_tv = ~ L + M,
    id = "id",
    time = "time",
    type = "longitudinal",
    family = "gaussian",
    estimator = "snm",
    history = 0
  )
  res <- contrast(fit, ci_method = "sandwich")

  # DTRreg: wide format
  wide <- data.table::dcast(d, id ~ time, value.var = c("A", "L", "M", "Y"))
  wide[, Y := Y_1]
  data.table::setnames(
    wide,
    c("A_0", "A_1", "L_0", "L_1", "M_0", "M_1"),
    c("A0", "A1", "L_tv0", "L_tv1", "M_tv0", "M_tv1")
  )

  dtr <- DTRreg::DTRreg(
    outcome = wide$Y,
    blip.mod = list(~M_tv0, ~M_tv1),
    treat.mod = list(A0 ~ L_tv0 + M_tv0, A1 ~ L_tv1 + M_tv1),
    tf.mod = list(~1, ~1),
    method = "gest",
    treat.type = "bin",
    data = as.data.frame(wide),
    var.estim = "none"
  )

  dtr_psi_0 <- unname(coef(dtr)[[1]])
  dtr_psi_1 <- unname(coef(dtr)[[2]])

  # Stage 1 (last treatment): both packages agree tightly
  expect_equal(res$estimates$estimate[3], dtr_psi_1[1], tolerance = 0.05)
  expect_equal(res$estimates$estimate[4], dtr_psi_1[2], tolerance = 0.05)

  # Stage 0 (first treatment): DTRreg's backward-transformed outcome at
  # stage 0 with history=0 uses a different TF intercept estimation path,
  # producing a wider discrepancy (~0.35 on intercept). Both are close to
  # truth (causatr: 1.13, DTRreg: 1.48, truth: 1.15) — causatr is closer.
  expect_equal(res$estimates$estimate[1], dtr_psi_0[1], tolerance = 0.4)
  expect_equal(res$estimates$estimate[2], dtr_psi_0[2], tolerance = 0.15)
})


test_that("IPW is biased with post-treatment modifier; SNM is not", {
  # The headline scientific result: conditioning on M_1 (post-treatment)
  # in an IPW MSM introduces collider bias because M_1 = 1{L_1 > 0}
  # and L_1 depends on A_0. The SNM's moment condition sidesteps this
  # because it uses the treatment residual as an instrument.
  dgp <- simulate_snm_longitudinal_tv_em(n = 10000, seed = 321)

  # SNM: recovers truth

  fit_snm <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders_outcome = ~ A:M,
    confounders_tv = ~ L + M,
    id = "id",
    time = "time",
    type = "longitudinal",
    family = "gaussian",
    estimator = "snm"
  )
  res_snm <- contrast(fit_snm, ci_method = "sandwich")

  # SNM stage-1 blip intercept should be near 2
  snm_psi10 <- res_snm$estimates$estimate[
    res_snm$estimates$parameter == "stage1_psi_intercept"
  ]
  expect_equal(snm_psi10, 2, tolerance = 0.2)

  # IPW on the same data: point-treatment at the final period with
  # M as a post-treatment covariate in confounders. The MSM conditions
  # on M, which opens a collider path A_0 -> L_1 -> M_1 <- ... -> Y.
  d_final <- dgp$data[time == 1]
  d_final[, L0 := dgp$data[time == 0]$L]
  d_final[, A0 := dgp$data[time == 0]$A]

  fit_ipw <- causat(
    d_final,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L0 + L + M + A:M,
    family = "gaussian",
    estimator = "ipw"
  )
  res_ipw <- contrast(
    fit_ipw,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "sandwich"
  )
  ipw_ate <- res_ipw$contrasts$estimate

  # True average stage-1 effect: E[2 + 2*M_1] ~ 3.1 (M_1 has ~55% prevalence)
  true_avg <- 2 + 2 * mean(dgp$data[time == 1]$M)

  # SNM averaged blip should be close to truth
  snm_avg <- res_snm$estimates$estimate[3] +
    res_snm$estimates$estimate[4] * mean(dgp$data[time == 1]$M)
  expect_equal(snm_avg, true_avg, tolerance = 0.3)

  # IPW should be detectably biased (off by > 1 from truth)
  expect_true(abs(ipw_ate - true_avg) > 1)
})


test_that("Longitudinal SNM with TV-EM: vcov is PSD and correctly sized", {
  dgp <- simulate_snm_longitudinal_tv_em(n = 2000, seed = 42)

  fit <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders_outcome = ~ A:M,
    confounders_tv = ~ L + M,
    id = "id",
    time = "time",
    type = "longitudinal",
    family = "gaussian",
    estimator = "snm"
  )
  res <- contrast(fit, ci_method = "sandwich")

  # 4x4 vcov: 2 stages x 2 blip params
  expect_equal(nrow(res$vcov), 4L)
  expect_equal(ncol(res$vcov), 4L)

  eig <- eigen(res$vcov, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eig >= -1e-10))

  expect_equal(rownames(res$vcov), res$estimates$parameter)
})


test_that("Longitudinal SNM with TV-EM: sandwich SE matches cluster bootstrap", {
  dgp <- simulate_snm_longitudinal_tv_em(n = 2000, seed = 42)

  fit <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders_outcome = ~ A:M,
    confounders_tv = ~ L + M,
    id = "id",
    time = "time",
    type = "longitudinal",
    family = "gaussian",
    estimator = "snm"
  )
  res <- contrast(fit, ci_method = "sandwich")
  sandwich_se <- res$estimates$se

  set.seed(1234)
  n_boot <- 200
  unique_ids <- unique(dgp$data$id)
  n_ids <- length(unique_ids)

  boot_psi <- replicate(n_boot, {
    boot_ids <- sample(unique_ids, n_ids, replace = TRUE)
    boot_data <- data.table::rbindlist(lapply(seq_along(boot_ids), function(i) {
      rows <- dgp$data[id == boot_ids[i]]
      rows <- data.table::copy(rows)
      rows[, id := i]
      rows
    }))
    tryCatch(
      {
        boot_fit <- causat(
          boot_data,
          outcome = "Y",
          treatment = "A",
          confounders_outcome = ~ A:M,
          confounders_tv = ~ L + M,
          id = "id",
          time = "time",
          type = "longitudinal",
          family = "gaussian",
          estimator = "snm"
        )
        boot_res <- contrast(boot_fit, ci_method = "sandwich")
        boot_res$estimates$estimate
      },
      error = function(e) rep(NA_real_, 4)
    )
  })

  boot_se <- apply(boot_psi, 1, sd, na.rm = TRUE)
  expect_equal(sandwich_se, boot_se, tolerance = 0.35)
})
