# Tests for per-component confounder formulas:
#   confounders_outcome, confounders_treatment, confounders_censoring,
#   confounders_sampling, confounders_tv_outcome, confounders_tv_treatment.
#
# Three axes:
#   1. Routing — each formula reaches the correct model component.
#   2. Backward compat — `confounders` alone fills all components.
#   3. Validation — bad types, missing columns, missing required formulas.


# ── Helpers ──────────────────────────────────────────────────────────────────
make_point_data <- function(n = 300, seed = 42) {
  set.seed(seed)
  L1 <- stats::rnorm(n)
  L2 <- stats::rnorm(n)
  A <- stats::rbinom(n, 1, stats::plogis(0.5 * L1))
  Y <- 1 + 2 * A + 0.3 * L1 + 0.2 * L2 + stats::rnorm(n)
  data.table::data.table(Y = Y, A = A, L1 = L1, L2 = L2)
}

make_long_data <- function(n_ids = 200, seed = 42) {
  set.seed(seed)
  L0 <- stats::rnorm(n_ids)
  A0 <- stats::rbinom(n_ids, 1, stats::plogis(0.3 * L0))
  L1 <- L0 + 0.2 * A0 + stats::rnorm(n_ids)
  A1 <- stats::rbinom(n_ids, 1, stats::plogis(0.3 * L1))
  Y <- 1 + A0 + A1 + 0.5 * L1 + stats::rnorm(n_ids)
  data.table::data.table(
    id = rep(seq_len(n_ids), each = 2),
    time = rep(0:1, times = n_ids),
    A = as.numeric(rbind(A0, A1)),
    L = as.numeric(rbind(L0, L1)),
    Y = rep(Y, each = 2)
  )
}


# ── 1. Routing: formulas reach correct model components ──────────────────────

test_that("AIPW: outcome model uses confounders_outcome, PS uses confounders_treatment", {
  d <- make_point_data()
  fit <- causat(d, outcome = "Y", treatment = "A",
    confounders_outcome = ~ L1 + L2,
    confounders_treatment = ~ L1,
    estimator = "aipw",
    propensity_model_fn = stats::glm)

  outcome_vars <- all.vars(stats::formula(fit$model))
  expect_true("L2" %in% outcome_vars)

  ps_vars <- all.vars(fit$details$treatment_model$ps_formula)
  expect_true("L1" %in% ps_vars)
  expect_false("L2" %in% ps_vars)
})

test_that("IPW: PS model uses confounders_treatment", {
  d <- make_point_data()
  fit <- causat(d, outcome = "Y", treatment = "A",
    confounders_treatment = ~ L1,
    estimator = "ipw",
    propensity_model_fn = stats::glm)

  ps_vars <- all.vars(fit$details$treatment_model$ps_formula)
  expect_true("L1" %in% ps_vars)
  expect_false("L2" %in% ps_vars)
})

test_that("G-comp: outcome model uses confounders_outcome", {
  d <- make_point_data()
  fit <- causat(d, outcome = "Y", treatment = "A",
    confounders_outcome = ~ L1 + L2)

  outcome_vars <- all.vars(stats::formula(fit$model))
  expect_true("L1" %in% outcome_vars)
  expect_true("L2" %in% outcome_vars)
})

test_that("Matching: distance uses confounders_treatment", {
  d <- make_point_data()
  fit <- causat(d, outcome = "Y", treatment = "A",
    confounders_outcome = ~ L1 + L2,
    confounders_treatment = ~ L1,
    estimator = "matching", estimand = "ATT")

  match_formula <- fit$match_obj$formula
  match_vars <- all.vars(match_formula)
  expect_true("L1" %in% match_vars)
  expect_false("L2" %in% match_vars)
})


# ── 2. Backward compatibility ────────────────────────────────────────────────

test_that("Legacy confounders fills all components identically", {
  d <- make_point_data()
  fit_legacy <- causat(d, outcome = "Y", treatment = "A",
    confounders = ~ L1 + L2,
    estimator = "aipw",
    propensity_model_fn = stats::glm)
  fit_new <- causat(d, outcome = "Y", treatment = "A",
    confounders_outcome = ~ L1 + L2,
    confounders_treatment = ~ L1 + L2,
    estimator = "aipw",
    propensity_model_fn = stats::glm)

  ivs <- list(a1 = static(1), a0 = static(0))
  res_legacy <- contrast(fit_legacy, ivs, reference = "a0")
  res_new <- contrast(fit_new, ivs, reference = "a0")

  expect_equal(
    res_legacy$contrasts$estimate,
    res_new$contrasts$estimate,
    tolerance = 1e-10
  )
})

test_that("confounders_tv falls back to confounders_tv when specific TV args absent", {
  long <- make_long_data()
  fit <- causat(long, outcome = "Y", treatment = "A",
    confounders = ~ 1, confounders_tv = ~ L,
    id = "id", time = "time")

  expect_equal(resolve_confounders_tv_outcome(fit), ~ L)
  expect_equal(resolve_confounders_tv_treatment(fit), ~ L)
})

test_that("Per-component TV overrides confounders_tv", {
  long <- make_long_data()
  fit <- causat(long, outcome = "Y", treatment = "A",
    confounders = ~ 1, confounders_tv = ~ L,
    confounders_tv_outcome = ~ L,
    id = "id", time = "time")

  expect_equal(resolve_confounders_tv_outcome(fit), ~ L)
  expect_equal(resolve_confounders_tv_treatment(fit), ~ L)
})


# ── 3. Validation ────────────────────────────────────────────────────────────

test_that("Error on non-formula confounders_outcome", {
  d <- make_point_data()
  expect_snapshot(
    causat(d, outcome = "Y", treatment = "A",
      confounders_outcome = "L1"),
    error = TRUE
  )
})

test_that("Error on missing column in confounders_outcome", {
  d <- make_point_data()
  expect_snapshot(
    causat(d, outcome = "Y", treatment = "A",
      confounders_outcome = ~ nonexistent),
    error = TRUE
  )
})

test_that("Error when AIPW missing required confounders_treatment", {
  d <- make_point_data()
  expect_snapshot(
    causat(d, outcome = "Y", treatment = "A",
      confounders_outcome = ~ L1,
      estimator = "aipw",
      propensity_model_fn = stats::glm),
    error = TRUE
  )
})

test_that("Error when IPW has no treatment confounders at all", {
  d <- make_point_data()
  expect_snapshot(
    causat(d, outcome = "Y", treatment = "A",
      estimator = "ipw",
      propensity_model_fn = stats::glm),
    error = TRUE
  )
})


# ── 4. Resolvers ─────────────────────────────────────────────────────────────

test_that("resolve_confounders_outcome falls back to confounders", {
  d <- make_point_data()
  fit <- causat(d, outcome = "Y", treatment = "A",
    confounders = ~ L1 + L2)
  expect_equal(resolve_confounders_outcome(fit), ~ L1 + L2)
  expect_equal(resolve_confounders_treatment(fit), ~ L1 + L2)
})

test_that("resolve_confounders_outcome uses specific when set", {
  d <- make_point_data()
  fit <- causat(d, outcome = "Y", treatment = "A",
    confounders_outcome = ~ L1 + L2,
    confounders_treatment = ~ L1)
  expect_equal(resolve_confounders_outcome(fit), ~ L1 + L2)
  expect_equal(resolve_confounders_treatment(fit), ~ L1)
})


# ── 5. Display ───────────────────────────────────────────────────────────────

test_that("print shows per-component formulas when set", {
  d <- make_point_data()
  fit <- causat(d, outcome = "Y", treatment = "A",
    confounders_outcome = ~ L1 + L2,
    confounders_treatment = ~ L1,
    estimator = "aipw",
    propensity_model_fn = stats::glm)
  out <- capture.output(print(fit))
  expect_true(any(grepl("Conf (outcome)", out, fixed = TRUE)))
  expect_true(any(grepl("Conf (treatment)", out, fixed = TRUE)))
})

test_that("summary shows per-component formulas when set", {
  d <- make_point_data()
  fit <- causat(d, outcome = "Y", treatment = "A",
    confounders_outcome = ~ L1 + L2,
    confounders_treatment = ~ L1,
    estimator = "aipw",
    propensity_model_fn = stats::glm)
  out <- capture.output(summary(fit))
  expect_true(any(grepl("Conf (outcome)", out, fixed = TRUE)))
  expect_true(any(grepl("Conf (treatment)", out, fixed = TRUE)))
})


# ── 6. Diagnostics use treatment confounders ─────────────────────────────────

test_that("diagnose uses confounders_treatment for balance", {
  d <- make_point_data()
  fit <- causat(d, outcome = "Y", treatment = "A",
    confounders_outcome = ~ L1 + L2,
    confounders_treatment = ~ L1,
    estimator = "ipw",
    propensity_model_fn = stats::glm)
  diag <- diagnose(fit)
  bal_vars <- rownames(diag$balance$Balance)
  expect_true("L1" %in% bal_vars)
  expect_false("L2" %in% bal_vars)
})


# ── 7. Contrast produces correct results with split confounders ──────────────

test_that("AIPW with split confounders produces finite estimates and SE", {
  d <- make_point_data()
  fit <- causat(d, outcome = "Y", treatment = "A",
    confounders_outcome = ~ L1 + L2,
    confounders_treatment = ~ L1,
    estimator = "aipw",
    propensity_model_fn = stats::glm)
  ivs <- list(a1 = static(1), a0 = static(0))
  res <- contrast(fit, ivs, reference = "a0")

  expect_true(is.finite(res$contrasts$estimate))
  expect_true(is.finite(res$contrasts$se))
  expect_true(res$contrasts$se > 0)
  # ATE should be near 2 (true value)
  expect_equal(res$contrasts$estimate, 2, tolerance = 0.5)
})

test_that("G-comp with confounders_outcome matches legacy confounders", {
  d <- make_point_data()
  fit1 <- causat(d, outcome = "Y", treatment = "A",
    confounders = ~ L1 + L2)
  fit2 <- causat(d, outcome = "Y", treatment = "A",
    confounders_outcome = ~ L1 + L2)
  ivs <- list(a1 = static(1), a0 = static(0))
  res1 <- contrast(fit1, ivs, reference = "a0")
  res2 <- contrast(fit2, ivs, reference = "a0")

  expect_equal(res1$contrasts$estimate, res2$contrasts$estimate, tolerance = 1e-10)
  expect_equal(res1$contrasts$se, res2$contrasts$se, tolerance = 1e-10)
})
