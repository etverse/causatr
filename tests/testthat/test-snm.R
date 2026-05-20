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
  fit <- suppressMessages(causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "snm"
  ))
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
  fit <- suppressMessages(causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "snm"
  ))

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
        boot_fit <- suppressMessages(causat(
          boot_data,
          outcome = "Y",
          treatment = "A",
          confounders = ~L,
          estimator = "snm"
        ))
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

  fit <- suppressMessages(causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "snm"
  ))
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

  fit <- suppressMessages(causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "snm"
  ))
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
  fit <- suppressMessages(causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "snm"
  ))
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
  fit <- suppressMessages(causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "snm"
  ))
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
  fit <- suppressMessages(causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "snm",
    treatment_free = ~L
  ))
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

  fit <- suppressMessages(causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "snm",
    treatment_free = ~L
  ))
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
  fit <- suppressMessages(causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "snm",
    treatment_free = ~L
  ))

  sandwich_result <- contrast(fit, ci_method = "sandwich")
  sandwich_se <- sandwich_result$estimates$se

  set.seed(606)
  n_boot <- 200
  boot_psi <- replicate(n_boot, {
    idx <- sample.int(nrow(dgp$data), replace = TRUE)
    boot_data <- dgp$data[idx]
    tryCatch(
      {
        boot_fit <- suppressMessages(causat(
          boot_data,
          outcome = "Y",
          treatment = "A",
          confounders = ~L,
          estimator = "snm",
          treatment_free = ~L
        ))
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
