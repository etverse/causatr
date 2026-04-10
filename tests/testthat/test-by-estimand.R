# Tests for by + estimand interaction across methods, estimands, and outcome types.
# DGP functions are defined in helper-dgp.R (auto-loaded by testthat).
#
# These tests use DGPs with heterogeneous treatment effects (τ depends on L)
# so that ATT ≠ ATE ≠ ATC within each sex stratum, ensuring we can distinguish
# correct from incorrect estimand targeting.

tol <- 0.6


# ============================================================
# G-COMP × BY × CONTINUOUS OUTCOME
# ============================================================

test_that("gcomp × by × ATE recovers subgroup-specific ATE (continuous)", {
  df <- simulate_het_effect(n = 5000)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L + sex + A:L + A:sex
  )
  result <- contrast(
    fit,
    list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    type = "difference",
    ci_method = "sandwich",
    by = "sex"
  )

  ate_sex0 <- result$contrasts$estimate[result$contrasts$by == "0"]
  ate_sex1 <- result$contrasts$estimate[result$contrasts$by == "1"]
  expect_equal(ate_sex0, 3.0, tolerance = tol)
  expect_equal(ate_sex1, 4.5, tolerance = tol)
})

test_that("gcomp × by × ATT recovers subgroup-specific ATT (continuous)", {
  df <- simulate_het_effect(n = 5000)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L + sex + A:L + A:sex
  )
  result <- contrast(
    fit,
    list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    type = "difference",
    ci_method = "sandwich",
    estimand = "ATT",
    by = "sex"
  )

  att_sex0 <- result$contrasts$estimate[result$contrasts$by == "0"]
  att_sex1 <- result$contrasts$estimate[result$contrasts$by == "1"]
  expect_equal(att_sex0, 3.83, tolerance = tol)
  expect_equal(att_sex1, 5.33, tolerance = tol)

  expect_true(att_sex0 > 3.0 + 0.3)
  expect_true(att_sex1 > 4.5 + 0.3)
})

test_that("gcomp × by × ATC recovers subgroup-specific ATC (continuous)", {
  df <- simulate_het_effect(n = 5000)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L + sex + A:L + A:sex
  )
  result <- contrast(
    fit,
    list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    type = "difference",
    ci_method = "sandwich",
    estimand = "ATC",
    by = "sex"
  )

  atc_sex0 <- result$contrasts$estimate[result$contrasts$by == "0"]
  atc_sex1 <- result$contrasts$estimate[result$contrasts$by == "1"]
  expect_equal(atc_sex0, 2.17, tolerance = tol)
  expect_equal(atc_sex1, 3.67, tolerance = tol)

  expect_true(atc_sex0 < 3.0 - 0.3)
  expect_true(atc_sex1 < 4.5 - 0.3)
})

test_that("gcomp × by × ATT matches manual subset (continuous)", {
  df <- simulate_het_effect(n = 5000)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L + sex + A:L + A:sex
  )

  result_by <- contrast(
    fit,
    list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich",
    estimand = "ATT",
    by = "sex"
  )

  result_manual_0 <- contrast(
    fit,
    list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich",
    subset = quote(sex == 0 & A == 1)
  )
  result_manual_1 <- contrast(
    fit,
    list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich",
    subset = quote(sex == 1 & A == 1)
  )

  by_est_0 <- result_by$contrasts$estimate[result_by$contrasts$by == "0"]
  by_est_1 <- result_by$contrasts$estimate[result_by$contrasts$by == "1"]
  manual_0 <- result_manual_0$contrasts$estimate[1]
  manual_1 <- result_manual_1$contrasts$estimate[1]

  expect_equal(by_est_0, manual_0, tolerance = 1e-10)
  expect_equal(by_est_1, manual_1, tolerance = 1e-10)
})


# ============================================================
# G-COMP × BY × BINARY OUTCOME
# ============================================================

test_that("gcomp × by × ATE recovers subgroup-specific RD (binary)", {
  df <- simulate_het_binary(n = 5000)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L + sex + A:L + A:sex,
    family = "binomial"
  )
  result <- contrast(
    fit,
    list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    type = "difference",
    ci_method = "sandwich",
    by = "sex"
  )

  rd_sex0 <- result$contrasts$estimate[result$contrasts$by == "0"]
  rd_sex1 <- result$contrasts$estimate[result$contrasts$by == "1"]
  expect_equal(rd_sex0, 0.287, tolerance = 0.1)
  expect_equal(rd_sex1, 0.364, tolerance = 0.1)
})

test_that("gcomp × by × ATT recovers subgroup-specific ATT RD (binary)", {
  df <- simulate_het_binary(n = 5000)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L + sex + A:L + A:sex,
    family = "binomial"
  )
  result <- contrast(
    fit,
    list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    type = "difference",
    ci_method = "sandwich",
    estimand = "ATT",
    by = "sex"
  )

  att_sex0 <- result$contrasts$estimate[result$contrasts$by == "0"]
  att_sex1 <- result$contrasts$estimate[result$contrasts$by == "1"]
  expect_equal(att_sex0, 0.346, tolerance = 0.1)
  expect_equal(att_sex1, 0.415, tolerance = 0.1)

  result_ate <- contrast(
    fit,
    list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich",
    by = "sex"
  )
  ate_sex0 <- result_ate$contrasts$estimate[result_ate$contrasts$by == "0"]
  expect_true(att_sex0 > ate_sex0)
})

test_that("gcomp × by × ATC recovers subgroup-specific ATC RD (binary)", {
  df <- simulate_het_binary(n = 5000)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L + sex + A:L + A:sex,
    family = "binomial"
  )
  result <- contrast(
    fit,
    list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    type = "difference",
    ci_method = "sandwich",
    estimand = "ATC",
    by = "sex"
  )

  atc_sex0 <- result$contrasts$estimate[result$contrasts$by == "0"]
  atc_sex1 <- result$contrasts$estimate[result$contrasts$by == "1"]
  expect_equal(atc_sex0, 0.228, tolerance = 0.1)
  expect_equal(atc_sex1, 0.312, tolerance = 0.1)
})


# ============================================================
# G-COMP × BY × RATIO AND OR CONTRASTS (BINARY OUTCOME)
# ============================================================

test_that("gcomp × by × ATT × ratio contrast (binary)", {
  df <- simulate_het_binary(n = 5000)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L + sex + A:L + A:sex,
    family = "binomial"
  )
  result <- contrast(
    fit,
    list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    type = "ratio",
    ci_method = "sandwich",
    estimand = "ATT",
    by = "sex"
  )

  expect_true(all(result$contrasts$estimate > 1))
  expect_true(all(result$contrasts$se > 0))
  expect_equal(nrow(result$contrasts), 2L)
  expect_equal(sort(unique(result$contrasts$by)), c("0", "1"))
})

test_that("gcomp × by × ATT × OR contrast (binary)", {
  df <- simulate_het_binary(n = 5000)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L + sex + A:L + A:sex,
    family = "binomial"
  )
  result <- contrast(
    fit,
    list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    type = "or",
    ci_method = "sandwich",
    estimand = "ATT",
    by = "sex"
  )

  expect_true(all(result$contrasts$estimate > 1))
  expect_true(all(result$contrasts$se > 0))
})


# ============================================================
# IPW × BY × ESTIMAND
# ============================================================

test_that("ipw × by × ATT: n_by reflects treated subset", {
  df <- simulate_het_effect(n = 3000)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L + sex,
    method = "ipw",
    estimand = "ATT"
  )
  result <- contrast(
    fit,
    list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich",
    by = "sex"
  )

  n_treated_sex0 <- sum(df$A == 1 & df$sex == 0)
  n_treated_sex1 <- sum(df$A == 1 & df$sex == 1)

  n_by_0 <- unique(result$contrasts$n_by[result$contrasts$by == "0"])
  n_by_1 <- unique(result$contrasts$n_by[result$contrasts$by == "1"])

  expect_equal(n_by_0, n_treated_sex0)
  expect_equal(n_by_1, n_treated_sex1)
  expect_true(n_by_0 + n_by_1 < nrow(df))
})

test_that("ipw × by × ATE: n_by sums to total", {
  df <- simulate_het_effect(n = 3000)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L + sex,
    method = "ipw"
  )
  result <- contrast(
    fit,
    list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich",
    by = "sex"
  )

  expect_equal(
    sum(unique(result$contrasts$n_by)),
    nrow(df)
  )
})

test_that("ipw × by × ATC: n_by reflects control subset", {
  df <- simulate_het_effect(n = 3000)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L + sex,
    method = "ipw",
    estimand = "ATC"
  )
  result <- contrast(
    fit,
    list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich",
    by = "sex"
  )

  n_control_sex0 <- sum(df$A == 0 & df$sex == 0)
  n_control_sex1 <- sum(df$A == 0 & df$sex == 1)

  n_by_0 <- unique(result$contrasts$n_by[result$contrasts$by == "0"])
  n_by_1 <- unique(result$contrasts$n_by[result$contrasts$by == "1"])

  expect_equal(n_by_0, n_control_sex0)
  expect_equal(n_by_1, n_control_sex1)
})


# ============================================================
# MATCHING × BY × ESTIMAND
# ============================================================

test_that("matching × by × ATT: n_by reflects treated subset", {
  df <- simulate_het_effect(n = 3000)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L + sex,
    method = "matching",
    estimand = "ATT"
  )
  result <- contrast(
    fit,
    list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich",
    by = "sex"
  )

  n_treated_sex0 <- sum(df$A == 1 & df$sex == 0)
  n_treated_sex1 <- sum(df$A == 1 & df$sex == 1)

  n_by_0 <- unique(result$contrasts$n_by[result$contrasts$by == "0"])
  n_by_1 <- unique(result$contrasts$n_by[result$contrasts$by == "1"])

  expect_equal(n_by_0, n_treated_sex0)
  expect_equal(n_by_1, n_treated_sex1)
})


# ============================================================
# BY WITH FACTOR VARIABLE
# ============================================================

test_that("gcomp × by with factor variable works correctly", {
  df <- simulate_het_effect(n = 5000)
  df$sex <- factor(df$sex, levels = 0:1, labels = c("Male", "Female"))

  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L + sex + A:L + A:sex
  )
  result <- contrast(
    fit,
    list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich",
    estimand = "ATT",
    by = "sex"
  )

  expect_equal(sort(unique(result$contrasts$by)), c("Female", "Male"))
  expect_equal(nrow(result$contrasts), 2L)

  att_male <- result$contrasts$estimate[result$contrasts$by == "Male"]
  att_female <- result$contrasts$estimate[result$contrasts$by == "Female"]
  expect_equal(att_male, 3.83, tolerance = tol)
  expect_equal(att_female, 5.33, tolerance = tol)
})


# ============================================================
# BY + ESTIMAND + SUBSET (combined)
# ============================================================

test_that("gcomp × by + subset composes correctly", {
  df <- simulate_het_effect(n = 5000)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L + sex + A:L + A:sex
  )

  result_by <- contrast(
    fit,
    list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich",
    subset = quote(L > 0),
    by = "sex"
  )

  result_manual <- contrast(
    fit,
    list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich",
    subset = quote(L > 0 & sex == 0)
  )

  by_est_0 <- result_by$contrasts$estimate[result_by$contrasts$by == "0"]
  expect_equal(by_est_0, result_manual$contrasts$estimate[1], tolerance = 1e-10)
})


# ============================================================
# BY + VCOV STRUCTURE
# ============================================================

test_that("vcov returns per-stratum list with by + estimand", {
  df <- simulate_het_effect(n = 3000)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L + sex + A:L + A:sex
  )
  result <- contrast(
    fit,
    list(a1 = static(1), a0 = static(0)),
    ci_method = "sandwich",
    estimand = "ATT",
    by = "sex"
  )

  v <- vcov(result)
  expect_true(is.list(v))
  expect_equal(length(v), 2L)
  expect_true(all(vapply(v, is.matrix, logical(1))))
})


# ============================================================
# BY + BOOTSTRAP
# ============================================================

test_that("gcomp × by × ATT × bootstrap works", {
  df <- simulate_het_effect(n = 2000)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L + sex + A:L + A:sex
  )
  result <- contrast(
    fit,
    list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "bootstrap",
    n_boot = 100,
    estimand = "ATT",
    by = "sex"
  )

  att_sex0 <- result$contrasts$estimate[result$contrasts$by == "0"]
  att_sex1 <- result$contrasts$estimate[result$contrasts$by == "1"]
  expect_equal(att_sex0, 3.83, tolerance = 1.0)
  expect_equal(att_sex1, 5.33, tolerance = 1.0)
  expect_true(all(result$contrasts$se > 0))
})


# ============================================================
# LONGITUDINAL (ICE) × BY
# ============================================================

test_that("ICE gcomp × by recovers subgroup effects", {
  set.seed(99)
  n <- 500
  sex <- rbinom(n, 1, 0.5)
  L0 <- rnorm(n)
  A0 <- rbinom(n, 1, plogis(0.3 * L0))
  L1 <- L0 + 0.5 * A0 + rnorm(n)
  A1 <- rbinom(n, 1, plogis(0.3 * L1))
  Y <- 2 +
    1.5 * A0 +
    1.5 * A1 +
    0.5 * L0 +
    0.5 * L1 +
    0.8 * sex * (A0 + A1) +
    rnorm(n)

  t0 <- data.frame(
    id = seq_len(n),
    time = 0L,
    A = A0,
    L = L0,
    sex = sex,
    Y = NA_real_
  )
  t1 <- data.frame(
    id = seq_len(n),
    time = 1L,
    A = A1,
    L = L1,
    sex = sex,
    Y = Y
  )
  long <- rbind(t0, t1)

  fit <- causat(
    long,
    outcome = "Y",
    treatment = "A",
    confounders = ~sex,
    confounders_tv = ~L,
    id = "id",
    time = "time"
  )
  result <- contrast(
    fit,
    interventions = list(
      always = static(1),
      never = static(0)
    ),
    reference = "never",
    ci_method = "sandwich",
    by = "sex"
  )

  expect_equal(nrow(result$contrasts), 2L)
  expect_equal(sort(unique(result$contrasts$by)), c("0", "1"))
  expect_true(all(result$contrasts$se > 0))

  eff_sex0 <- result$contrasts$estimate[result$contrasts$by == "0"]
  eff_sex1 <- result$contrasts$estimate[result$contrasts$by == "1"]
  expect_true(eff_sex1 > eff_sex0)
})


# ============================================================
# IPW × BY × BINARY OUTCOME
# ============================================================

test_that("ipw × by × ATT × binary outcome", {
  df <- simulate_het_binary(n = 3000)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L + sex,
    method = "ipw",
    family = "binomial",
    estimand = "ATT"
  )
  result <- contrast(
    fit,
    list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    type = "difference",
    ci_method = "sandwich",
    by = "sex"
  )

  expect_equal(nrow(result$contrasts), 2L)
  expect_true(all(result$contrasts$se > 0))
  expect_true(all(result$contrasts$estimate > 0))
})


# ============================================================
# MATCHING × BY × BINARY OUTCOME
# ============================================================

test_that("matching × by × ATT × binary outcome", {
  df <- simulate_het_binary(n = 3000)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L + sex,
    method = "matching",
    family = "binomial",
    estimand = "ATT"
  )
  result <- contrast(
    fit,
    list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    type = "difference",
    ci_method = "sandwich",
    by = "sex"
  )

  expect_equal(nrow(result$contrasts), 2L)
  expect_true(all(result$contrasts$se > 0))
  expect_true(all(result$contrasts$estimate > 0))
})


# ============================================================
# CONTINUOUS TREATMENT × BY
# ============================================================

test_that("gcomp × by × continuous treatment × shift", {
  set.seed(42)
  n <- 3000
  L <- rnorm(n)
  sex <- rbinom(n, 1, 0.5)
  A <- 1 + 0.5 * L + rnorm(n)
  Y <- 1 + (2 + 0.5 * sex) * A + L + rnorm(n)
  df <- data.frame(Y = Y, A = A, L = L, sex = sex)

  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L + sex + A:sex
  )
  result <- contrast(
    fit,
    interventions = list(shifted = shift(-1), observed = NULL),
    reference = "observed",
    ci_method = "sandwich",
    by = "sex"
  )

  eff_sex0 <- result$contrasts$estimate[result$contrasts$by == "0"]
  eff_sex1 <- result$contrasts$estimate[result$contrasts$by == "1"]
  expect_equal(eff_sex0, -2.0, tolerance = 0.5)
  expect_equal(eff_sex1, -2.5, tolerance = 0.5)
})
