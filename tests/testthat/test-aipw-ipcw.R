# Tests for AIPW + IPCW composition (Phase 16o).
#
# Uses DGP-M2b (non-linear outcome with interaction, differential
# censoring) where IPCW genuinely matters: differential censoring
# biases complete-case AIPW but the triply-weighted DR estimator
# corrects for it.
#
# Validates:
#   - AIPW + IPCW point estimate recovers truth (ATE = 3)
#   - Sandwich SE is finite and CI covers truth
#   - AIPW + IPCW reduces bias vs naive complete-case AIPW
#   - Sandwich vs bootstrap SE agreement
#   - DR property: wrong outcome + correct propensity + IPCW recovers truth
#   - AIPW + IPCW + external weights composes correctly
#   - Details stashing: ipcw, censoring_model, ipcw_weights present

test_that("AIPW + IPCW recovers ATE on DGP-M2b", {
  d <- simulate_mar_outcome_complex(n = 5000, seed = 300)
  dt <- data.table::as.data.table(d)

  fit <- causat(
    dt,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L1 + L2,
    estimator = "aipw",
    censoring = "C",
    ipcw = TRUE,
    propensity_model_fn = stats::glm
  )

  result <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "sandwich"
  )

  ate <- -result$contrasts$estimate
  expect_equal(ate, 3, tolerance = 0.15)
})

test_that("AIPW + IPCW: sandwich SE finite and CI covers truth", {
  d <- simulate_mar_outcome_complex(n = 3000, seed = 301)
  dt <- data.table::as.data.table(d)

  fit <- causat(
    dt,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L1 + L2,
    estimator = "aipw",
    censoring = "C",
    ipcw = TRUE,
    propensity_model_fn = stats::glm
  )

  result <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "sandwich"
  )

  se <- result$contrasts$se
  expect_true(is.finite(se))
  expect_true(se > 0)
  ci_lo <- result$contrasts$ci_lower
  ci_hi <- result$contrasts$ci_upper
  expect_true(ci_lo < -3 && ci_hi > -3)
})

test_that("AIPW + IPCW reduces bias vs naive complete-case AIPW", {
  d <- simulate_mar_outcome_complex(n = 5000, seed = 310)
  dt <- data.table::as.data.table(d)

  fit_ipcw <- causat(
    dt,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L1 + L2,
    estimator = "aipw",
    censoring = "C",
    ipcw = TRUE,
    propensity_model_fn = stats::glm
  )
  r_ipcw <- contrast(
    fit_ipcw,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "sandwich"
  )

  fit_naive <- causat(
    dt,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L1 + L2,
    estimator = "aipw",
    censoring = "C",
    ipcw = FALSE,
    propensity_model_fn = stats::glm
  )
  r_naive <- contrast(
    fit_naive,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "sandwich"
  )

  ate_ipcw <- -r_ipcw$contrasts$estimate
  ate_naive <- -r_naive$contrasts$estimate
  expect_lt(abs(ate_ipcw - 3), abs(ate_naive - 3))
})

test_that("AIPW + IPCW: sandwich vs bootstrap SE ratio in (0.5, 2.0)", {
  skip_if_fast()
  d <- simulate_mar_outcome_complex(n = 2000, seed = 320)
  dt <- data.table::as.data.table(d)

  fit <- causat(
    dt,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L1 + L2,
    estimator = "aipw",
    censoring = "C",
    ipcw = TRUE,
    propensity_model_fn = stats::glm
  )

  r_sand <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "sandwich"
  )

  r_boot <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "bootstrap",
    n_boot = 200
  )

  ratio <- r_sand$contrasts$se / r_boot$contrasts$se
  expect_gt(ratio, 0.5)
  expect_lt(ratio, 2.0)
})

test_that("AIPW + IPCW DR: wrong outcome model still recovers truth", {
  skip_if_fast()
  d <- simulate_mar_outcome_complex(n = 5000, seed = 330)
  dt <- data.table::as.data.table(d)
  # Misspecify outcome: omit the interaction and quadratic terms
  # that the DGP uses (A*L1 and L1^2). Propensity + censoring correct.
  fit <- causat(
    dt,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L1 + L2,
    estimator = "aipw",
    censoring = "C",
    ipcw = TRUE,
    propensity_model_fn = stats::glm
  )

  r_boot <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "bootstrap",
    n_boot = 300
  )

  ate <- -r_boot$contrasts$estimate
  expect_equal(ate, 3, tolerance = 0.3)
  ci_lo <- r_boot$contrasts$ci_lower
  ci_hi <- r_boot$contrasts$ci_upper
  expect_true(ci_lo < -3 && ci_hi > -3)
})

test_that("AIPW + IPCW + external weights composes correctly", {
  d <- simulate_mar_outcome_complex(n = 3000, seed = 340)
  dt <- data.table::as.data.table(d)
  ext_w <- runif(nrow(dt), 0.5, 1.5)

  fit <- suppressWarnings(causat(
    dt,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L1 + L2,
    estimator = "aipw",
    censoring = "C",
    ipcw = TRUE,
    weights = ext_w,
    propensity_model_fn = stats::glm
  ))

  result <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "sandwich"
  )

  se <- result$contrasts$se
  expect_true(is.finite(se))
  expect_true(se > 0)
  ate <- -result$contrasts$estimate
  expect_equal(ate, 3, tolerance = 0.3)
})

test_that("AIPW + IPCW stashes correct details", {
  d <- simulate_mar_outcome_complex(n = 1000, seed = 350)
  dt <- data.table::as.data.table(d)

  fit <- causat(
    dt,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L1 + L2,
    estimator = "aipw",
    censoring = "C",
    ipcw = TRUE,
    propensity_model_fn = stats::glm
  )

  expect_true(isTRUE(fit$details$ipcw))
  expect_s3_class(fit$details$censoring_model, "causatr_censoring_model")
  expect_true(is.numeric(fit$details$ipcw_weights))
  expect_equal(length(fit$details$ipcw_weights), nrow(dt))
  expect_true(all(fit$details$ipcw_weights[dt$C == 1] == 0))
  expect_true(all(fit$details$ipcw_weights[dt$C == 0] > 0))
})


# AIPW + IPCW + transport: sandwich variance must not error on
# `target_fit` lookup. Before the fix, `target_fit` was only defined
# in the non-transport branch but referenced by the IPCW censoring
# correction unconditionally.
# Critical review round 2026-05-17, Issue #1.
# Repro: /tmp/causatr_repro_ipcw_transport.R
test_that("AIPW + IPCW + transport: sandwich SE runs and agrees with bootstrap", {
  skip_if_fast()
  set.seed(500)
  n <- 2000
  L <- rnorm(n)
  ps_s <- plogis(-0.3 + 0.8 * L)
  S <- rbinom(n, 1, ps_s)

  A <- ifelse(S == 1L, rbinom(n, 1, plogis(0.2 + 0.3 * L)), NA_integer_)
  Y <- ifelse(S == 1L, 1 + 2 * A + 0.5 * L + rnorm(n), NA_real_)
  C <- ifelse(S == 1L, rbinom(n, 1, plogis(-2 + 0.3 * L)), NA_integer_)
  Y[C == 1L & !is.na(C)] <- NA_real_

  dt <- data.table::data.table(Y = Y, A = A, L = L, S = S, C = C)

  fit <- causat(
    dt,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "aipw",
    target = "S",
    censoring = "C",
    ipcw = TRUE,
    propensity_model_fn = stats::glm
  )

  result_sw <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "sandwich"
  )

  # Target ATE = 2 + E[0 | S=0] = 2 (no interaction), so contrast ~ -2
  ate <- result_sw$contrasts$estimate
  expect_equal(ate, -2, tolerance = 0.3)

  se_sw <- result_sw$contrasts$se
  expect_true(all(is.finite(se_sw)))
  expect_true(all(se_sw > 0))

  # Sandwich vs bootstrap agreement
  result_bt <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "bootstrap",
    n_boot = 200
  )
  se_bt <- result_bt$contrasts$se
  ratio <- se_sw / se_bt
  expect_true(ratio > 0.5 && ratio < 2.0)
})
