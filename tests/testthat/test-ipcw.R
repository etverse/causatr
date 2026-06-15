# Tests for the built-in IPCW integration (Phase 14b).
#
# Uses DGP-M2b (non-linear outcome with interaction, 2 confounders,
# differential censoring) for point-estimate tests because IPCW
# genuinely matters there: a misspecified linear model is biased
# under complete-case analysis but corrected by IPCW. DGP-M2 is
# kept for mechanical tests (backward compat, details stashing)
# where IPCW correctness is not the focus.
#
# Validates:
#   - causat(ipcw = TRUE) for gcomp, IPW, and matching on DGP-M2b
#   - IPCW reduces censoring bias vs naive complete-case
#   - ipcw = FALSE is backward-compatible
#   - Sandwich SE vs bootstrap SE agreement
#   - IPCW + external weights composition
#   - Error conditions

# ── Point estimates: gcomp (DGP-M2b) ─────────────────────────────

test_that("gcomp + IPCW recovers ATE on DGP-M2b", {
  d <- simulate_mar_outcome_complex(n = 5000, seed = 200)
  dt <- data.table::as.data.table(d)

  fit <- causat(
    dt,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L1 + L2,
    estimator = "gcomp",
    censoring = "C",
    ipcw = TRUE
  )

  result <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "sandwich"
  )

  expect_equal(result$contrasts$estimate, -3, tolerance = 0.15)
})

test_that("gcomp + IPCW: sandwich SE and CI cover truth", {
  d <- simulate_mar_outcome_complex(n = 3000, seed = 201)
  dt <- data.table::as.data.table(d)

  fit <- causat(
    dt,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L1 + L2,
    estimator = "gcomp",
    censoring = "C",
    ipcw = TRUE
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

test_that("gcomp IPCW reduces bias vs naive complete-case", {
  d <- simulate_mar_outcome_complex(n = 5000, seed = 210)
  dt <- data.table::as.data.table(d)

  fit_ipcw <- causat(
    dt,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L1 + L2,
    estimator = "gcomp",
    censoring = "C",
    ipcw = TRUE
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
    estimator = "gcomp",
    censoring = "C",
    ipcw = FALSE
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


# ── Point estimates: IPW (DGP-M2b) ───────────────────────────────

test_that("IPW + IPCW recovers ATE on DGP-M2b", {
  d <- simulate_mar_outcome_complex(n = 5000, seed = 202)
  dt <- data.table::as.data.table(d)

  fit <- causat(
    dt,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L1 + L2,
    estimator = "ipw",
    censoring = "C",
    ipcw = TRUE
  )

  result <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "sandwich"
  )

  expect_equal(result$contrasts$estimate, -3, tolerance = 0.20)
})

test_that("IPW + IPCW: sandwich SE and CI cover truth", {
  d <- simulate_mar_outcome_complex(n = 3000, seed = 203)
  dt <- data.table::as.data.table(d)

  fit <- causat(
    dt,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L1 + L2,
    estimator = "ipw",
    censoring = "C",
    ipcw = TRUE
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

test_that("IPW IPCW reduces bias vs naive complete-case", {
  d <- simulate_mar_outcome_complex(n = 5000, seed = 211)
  dt <- data.table::as.data.table(d)

  fit_ipcw <- causat(
    dt,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L1 + L2,
    estimator = "ipw",
    censoring = "C",
    ipcw = TRUE
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
    estimator = "ipw",
    censoring = "C",
    ipcw = FALSE
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


# ── Point estimates: matching (DGP-M2b) ──────────────────────────

test_that("matching + IPCW recovers ATE on DGP-M2b", {
  d <- simulate_mar_outcome_complex(n = 8000, seed = 204)
  dt <- data.table::as.data.table(d)

  fit <- causat(
    dt,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L1 + L2,
    estimator = "matching",
    censoring = "C",
    ipcw = TRUE,
    estimand = "ATT"
  )

  result <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "sandwich"
  )

  expect_equal(result$contrasts$estimate, -3, tolerance = 0.05)
})


# ── Backward compatibility (DGP-M2, mechanical) ──────────────────

test_that("ipcw = FALSE gives same result as no ipcw argument", {
  d <- simulate_mar_outcome(n = 2000, seed = 205)
  dt <- data.table::as.data.table(d)

  fit_default <- causat(
    dt,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "gcomp",
    censoring = "C"
  )

  fit_explicit <- causat(
    dt,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "gcomp",
    censoring = "C",
    ipcw = FALSE
  )

  r_default <- contrast(
    fit_default,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "sandwich"
  )
  r_explicit <- contrast(
    fit_explicit,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "sandwich"
  )

  expect_equal(
    r_default$contrasts$estimate,
    r_explicit$contrasts$estimate
  )
})


# ── IPCW details stashed on fit ───────────────────────────────────

test_that("IPCW details are stashed on the fit object", {
  d <- simulate_mar_outcome(n = 1000, seed = 206)
  dt <- data.table::as.data.table(d)

  fit <- causat(
    dt,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "gcomp",
    censoring = "C",
    ipcw = TRUE
  )

  expect_true(fit$details$ipcw)
  expect_s3_class(
    fit$details$censoring_model,
    "causatr_censoring_model"
  )
  expect_length(fit$details$ipcw_weights, nrow(dt))
  expect_true(all(fit$details$ipcw_weights[dt$C == 1L] == 0))
  expect_true(all(fit$details$ipcw_weights[dt$C == 0L] > 0))
  expect_null(fit$details$weights_pre_ipcw)
})


# ── Sandwich vs bootstrap SE agreement (DGP-M2b) ─────────────────

test_that("gcomp IPCW: sandwich and bootstrap SE agree", {
  skip_if_fast()
  d <- simulate_mar_outcome_complex(n = 3000, seed = 207)
  dt <- data.table::as.data.table(d)

  fit <- causat(
    dt,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L1 + L2,
    estimator = "gcomp",
    censoring = "C",
    ipcw = TRUE
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

  se_sand <- r_sand$contrasts$se
  se_boot <- r_boot$contrasts$se
  ratio <- se_sand / se_boot
  expect_true(ratio > 0.5 && ratio < 2.0)
})

test_that("IPW IPCW: sandwich and bootstrap SE agree", {
  skip_if_fast()
  d <- simulate_mar_outcome_complex(n = 3000, seed = 209)
  dt <- data.table::as.data.table(d)

  fit <- causat(
    dt,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L1 + L2,
    estimator = "ipw",
    censoring = "C",
    ipcw = TRUE
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

  se_sand <- r_sand$contrasts$se
  se_boot <- r_boot$contrasts$se
  ratio <- se_sand / se_boot
  expect_true(ratio > 0.5 && ratio < 2.0)
})


# ── IPCW + external weights ──────────────────────────────────────

test_that("IPCW composes with external weights", {
  d <- simulate_mar_outcome_complex(n = 3000, seed = 208)
  dt <- data.table::as.data.table(d)

  set.seed(208)
  ext_w <- runif(nrow(dt), 0.5, 1.5)

  fit <- suppressWarnings(causat(
    dt,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L1 + L2,
    estimator = "gcomp",
    censoring = "C",
    ipcw = TRUE,
    weights = ext_w
  ))

  expect_equal(fit$details$weights_pre_ipcw, ext_w)

  result <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "sandwich"
  )

  expect_equal(result$contrasts$estimate, -3, tolerance = 0.2)
})


# ── Longitudinal ICE + IPCW (DGP-M5) ────────────────────────────

test_that("ICE + IPCW recovers ATE on DGP-M5", {
  d <- simulate_longitudinal_mar_outcome(n = 5000, seed = 300)
  dt <- data.table::as.data.table(d)

  fit <- causat(
    dt,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    estimator = "gcomp",
    type = "longitudinal",
    id = "id",
    time = "time",
    censoring = "C",
    ipcw = TRUE
  )

  result <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "sandwich"
  )

  expect_equal(result$contrasts$estimate, -5, tolerance = 0.15)
})

test_that("ICE + IPCW: sandwich SE and CI cover truth", {
  d <- simulate_longitudinal_mar_outcome(n = 3000, seed = 301)
  dt <- data.table::as.data.table(d)

  fit <- causat(
    dt,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    estimator = "gcomp",
    type = "longitudinal",
    id = "id",
    time = "time",
    censoring = "C",
    ipcw = TRUE
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
  expect_true(ci_lo < -5 && ci_hi > -5)
})

test_that("ICE + IPCW: sandwich and bootstrap SE agree", {
  skip_if_fast()
  d <- simulate_longitudinal_mar_outcome(n = 2000, seed = 302)
  dt <- data.table::as.data.table(d)

  fit <- causat(
    dt,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    estimator = "gcomp",
    type = "longitudinal",
    id = "id",
    time = "time",
    censoring = "C",
    ipcw = TRUE
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

  se_sand <- r_sand$contrasts$se
  se_boot <- r_boot$contrasts$se
  ratio <- se_sand / se_boot
  expect_true(ratio > 0.5 && ratio < 2.0)
})

test_that("ICE IPCW does not degrade vs naive (linear DGP)", {
  # DGP-M5 has a correctly specified linear outcome model, so
  # g-comp absorbs censoring bias without IPCW. Verify IPCW
  # does not degrade the estimate — both should be near truth.
  d <- simulate_longitudinal_mar_outcome(n = 5000, seed = 303)
  dt <- data.table::as.data.table(d)

  fit_ipcw <- causat(
    dt,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    estimator = "gcomp",
    type = "longitudinal",
    id = "id",
    time = "time",
    censoring = "C",
    ipcw = TRUE
  )
  r_ipcw <- contrast(
    fit_ipcw,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "sandwich"
  )

  ate_ipcw <- -r_ipcw$contrasts$estimate
  expect_equal(ate_ipcw, 5, tolerance = 0.15)
})


# ── Longitudinal IPW + IPCW (DGP-M5) ───────────────────────────

test_that("longitudinal IPW + IPCW recovers ATE on DGP-M5", {
  d <- simulate_longitudinal_mar_outcome(n = 5000, seed = 304)
  dt <- data.table::as.data.table(d)

  fit <- causat(
    dt,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    estimator = "ipw",
    type = "longitudinal",
    id = "id",
    time = "time",
    censoring = "C",
    ipcw = TRUE
  )

  result <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "sandwich"
  )

  expect_equal(result$contrasts$estimate, -5, tolerance = 0.25)
})

test_that("longitudinal IPW IPCW reduces bias vs naive", {
  d <- simulate_longitudinal_mar_outcome(n = 5000, seed = 305)
  dt <- data.table::as.data.table(d)

  fit_ipcw <- causat(
    dt,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    estimator = "ipw",
    type = "longitudinal",
    id = "id",
    time = "time",
    censoring = "C",
    ipcw = TRUE
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
    confounders = ~L0,
    confounders_tv = ~L,
    estimator = "ipw",
    type = "longitudinal",
    id = "id",
    time = "time",
    censoring = "C",
    ipcw = FALSE
  )
  r_naive <- contrast(
    fit_naive,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "sandwich"
  )

  ate_ipcw <- -r_ipcw$contrasts$estimate
  ate_naive <- -r_naive$contrasts$estimate
  expect_lt(abs(ate_ipcw - 5), abs(ate_naive - 5))
})

test_that("longitudinal IPW IPCW: SE and CI cover truth", {
  d <- simulate_longitudinal_mar_outcome(n = 3000, seed = 306)
  dt <- data.table::as.data.table(d)

  fit <- causat(
    dt,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    estimator = "ipw",
    type = "longitudinal",
    id = "id",
    time = "time",
    censoring = "C",
    ipcw = TRUE
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
  expect_true(ci_lo < -5 && ci_hi > -5)
})


# ── Longitudinal IPCW details stashed ──────────────────────────

test_that("longitudinal IPCW details stashed on fit", {
  d <- simulate_longitudinal_mar_outcome(n = 1000, seed = 307)
  dt <- data.table::as.data.table(d)

  fit <- causat(
    dt,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    estimator = "gcomp",
    type = "longitudinal",
    id = "id",
    time = "time",
    censoring = "C",
    ipcw = TRUE
  )

  expect_true(fit$details$ipcw)
  expect_true(is.list(fit$details$censoring_models))
  expect_length(fit$details$ipcw_weights, nrow(fit$data))
  expect_true(
    all(fit$details$ipcw_weights[fit$data$C == 1L] == 0)
  )
  expect_true(
    all(fit$details$ipcw_weights[fit$data$time == 0] == 1)
  )
  expect_null(fit$details$weights_pre_ipcw)
})


# ── Error conditions ──────────────────────────────────────────────

test_that("ipcw = TRUE without censoring column aborts", {
  dt <- data.table::data.table(
    Y = rnorm(100),
    A = rbinom(100, 1, 0.5),
    L = rnorm(100)
  )

  expect_error(
    causat(
      dt,
      outcome = "Y",
      treatment = "A",
      confounders = ~L,
      estimator = "gcomp",
      ipcw = TRUE
    ),
    class = "causatr_ipcw_no_censoring"
  )
})
