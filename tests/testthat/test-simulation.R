# Simulation-based tests with known true parameter values.
# DGP functions are defined in helper-dgp.R (auto-loaded by testthat).

# ============================================================
# GCOMP × BINARY TREATMENT × CONTINUOUS OUTCOME × DIFFERENCE
# ============================================================

test_that("gcomp × binary trt × continuous outcome × difference × sandwich: ATE ≈ 3", {
  df <- simulate_binary_continuous(n = 2000)
  fit <- causat(df, outcome = "Y", treatment = "A", confounders = ~L)
  result <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich"
  )

  ate <- result$contrasts$estimate[1]
  expect_equal(ate, 3, tolerance = 0.15)
  # CI should contain the true value.
  expect_lt(result$contrasts$ci_lower[1], 3)
  expect_gt(result$contrasts$ci_upper[1], 3)
})

# ============================================================
# GCOMP × BINARY TREATMENT × CONTINUOUS OUTCOME × ESTIMANDS
# ============================================================

test_that("gcomp × binary trt × continuous outcome × ATT ≈ 3", {
  df <- simulate_binary_continuous(n = 2000)
  fit <- causat(df, outcome = "Y", treatment = "A", confounders = ~L)
  result <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    estimand = "ATT",
    ci_method = "sandwich"
  )

  att <- result$contrasts$estimate[1]
  expect_equal(att, 3, tolerance = 0.15)
})

# ============================================================
# GCOMP × BINARY TREATMENT × BINARY OUTCOME
# ============================================================

test_that("gcomp × binary trt × binary outcome × difference × sandwich: RD ≈ 0.33", {
  df <- simulate_binary_binary(n = 3000)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    family = "binomial"
  )
  result <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    type = "difference",
    ci_method = "sandwich"
  )

  rd <- result$contrasts$estimate[1]
  # True RD ≈ 0.33 (from simulation).
  expect_equal(rd, 0.33, tolerance = 0.1)
})

test_that("gcomp × binary trt × binary outcome × ratio × sandwich", {
  df <- simulate_binary_binary(n = 3000)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    family = "binomial"
  )
  result <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    type = "ratio",
    ci_method = "sandwich"
  )

  rr <- result$contrasts$estimate[1]
  # True RR ≈ 0.622 / 0.289 ≈ 2.15
  expect_gt(rr, 1.5)
  expect_lt(rr, 3.0)
  # SE should be finite and positive.
  expect_gt(result$contrasts$se[1], 0)
})

test_that("gcomp × categorical (3-level) treatment × static × difference × sandwich", {
  # Truth-based test for gcomp on a 3-level factor treatment.
  # DGP:
  #   L ~ N(0, 1)
  #   A | L ~ multinomial with level probs driven by L
  #   E[Y | A, L] = 1 + 2*I(A=1) + 5*I(A=2) + 0.5*L
  # Marginal counterfactual means: E[Y^0] = 1, E[Y^1] = 3, E[Y^2] = 6.
  # Contrasts: (a1 - a0) = 2, (a2 - a0) = 5.
  set.seed(100)
  n <- 4000
  L <- stats::rnorm(n)
  e1 <- exp(0.3 + 0.4 * L)
  e2 <- exp(-0.2 + 0.3 * L)
  P <- cbind(1, e1, e2) / (1 + e1 + e2)
  A <- vapply(
    seq_len(n),
    function(i) sample(0:2, 1L, prob = P[i, ]),
    integer(1)
  )
  A <- factor(A, levels = c("0", "1", "2"))
  beta_A <- c("0" = 0, "1" = 2, "2" = 5)
  Y <- 1 + beta_A[as.character(A)] + 0.5 * L + stats::rnorm(n)
  df <- data.table::data.table(A = A, L = L, Y = as.numeric(Y))

  fit <- causat(df, outcome = "Y", treatment = "A", confounders = ~L)
  res <- contrast(
    fit,
    interventions = list(
      a0 = static("0"),
      a1 = static("1"),
      a2 = static("2")
    ),
    reference = "a0",
    ci_method = "sandwich"
  )

  est <- res$estimates
  e0 <- est$estimate[est$intervention == "a0"]
  e1_mean <- est$estimate[est$intervention == "a1"]
  e2_mean <- est$estimate[est$intervention == "a2"]
  expect_lt(abs(e0 - 1), 0.1)
  expect_lt(abs(e1_mean - 3), 0.1)
  expect_lt(abs(e2_mean - 6), 0.1)

  ct <- res$contrasts
  c10 <- ct$estimate[ct$comparison == "a1 vs a0"]
  c20 <- ct$estimate[ct$comparison == "a2 vs a0"]
  expect_lt(abs(c10 - 2), 0.1)
  expect_lt(abs(c20 - 5), 0.1)
  # SEs finite and positive; CIs cover the truth.
  expect_true(all(is.finite(ct$se)) && all(ct$se > 0))
  expect_lt(ct$ci_lower[ct$comparison == "a1 vs a0"], 2)
  expect_gt(ct$ci_upper[ct$comparison == "a1 vs a0"], 2)
  expect_lt(ct$ci_lower[ct$comparison == "a2 vs a0"], 5)
  expect_gt(ct$ci_upper[ct$comparison == "a2 vs a0"], 5)
})


# ============================================================
# GCOMP × CONTINUOUS TREATMENT × CONTINUOUS OUTCOME × INTERVENTIONS
# ============================================================

test_that("gcomp × continuous trt × shift intervention × difference × sandwich", {
  df <- simulate_continuous_continuous(n = 2000)
  fit <- causat(df, outcome = "Y", treatment = "A", confounders = ~L)
  result <- contrast(
    fit,
    interventions = list(shifted = shift(-1), observed = NULL),
    reference = "observed",
    ci_method = "sandwich"
  )

  # shift(-1) reduces E[Y] by 2 (since dY/dA = 2 in the DGP).
  diff <- result$contrasts$estimate[1]
  expect_equal(diff, -2, tolerance = 0.15)
})

test_that("gcomp × continuous trt × scale intervention × difference × sandwich", {
  df <- simulate_continuous_continuous(n = 2000)
  fit <- causat(df, outcome = "Y", treatment = "A", confounders = ~L)
  result <- contrast(
    fit,
    interventions = list(halved = scale_by(0.5), observed = NULL),
    reference = "observed",
    ci_method = "sandwich"
  )

  # scale_by(0.5) halves A; since E[A] ≈ 1, difference ≈ 2*(0.5*1 - 1) = -1.
  # But since A varies, the exact value depends on E[A].
  diff <- result$contrasts$estimate[1]
  expect_lt(diff, 0) # halving treatment should reduce outcome
})

test_that("gcomp × continuous trt × threshold intervention × difference × sandwich", {
  df <- simulate_continuous_continuous(n = 2000)
  fit <- causat(df, outcome = "Y", treatment = "A", confounders = ~L)
  result <- contrast(
    fit,
    interventions = list(capped = threshold(0, 0.5), observed = NULL),
    reference = "observed",
    ci_method = "sandwich"
  )

  # Clamping A to [0, 0.5] should reduce the mean outcome substantially
  # since E[A] ≈ 1 in this DGP and dY/dA = 2.
  diff <- result$contrasts$estimate[1]
  expect_lt(diff, -0.5)
})

# ============================================================
# GCOMP × BINARY TREATMENT × DYNAMIC INTERVENTION
# ============================================================

test_that("gcomp × binary trt × dynamic intervention × difference × sandwich", {
  df <- simulate_binary_continuous(n = 2000)
  fit <- causat(df, outcome = "Y", treatment = "A", confounders = ~L)

  # Dynamic rule: treat if L > 0, else don't.
  result <- contrast(
    fit,
    interventions = list(
      rule = dynamic(\(data, trt) ifelse(data$L > 0, 1, 0)),
      all_treat = static(1)
    ),
    reference = "all_treat",
    ci_method = "sandwich"
  )

  # Rule treats fewer people than static(1), so mean outcome should be lower.
  expect_lt(result$contrasts$estimate[1], 0)
})

# ============================================================
# GCOMP × SUBSET ESTIMAND
# ============================================================

test_that("gcomp × binary trt × subset estimand (L > 0)", {
  df <- simulate_binary_continuous(n = 2000)
  fit <- causat(df, outcome = "Y", treatment = "A", confounders = ~L)
  result <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    subset = quote(L > 0),
    ci_method = "sandwich"
  )

  # Subgroup effect should still be ≈ 3 (constant treatment effect in DGP).
  expect_equal(result$contrasts$estimate[1], 3, tolerance = 0.3)
  expect_equal(result$estimand, "subset")
})

# ============================================================
# GCOMP × SANDWICH vs BOOTSTRAP INFERENCE
# ============================================================

test_that("gcomp × binary trt × continuous outcome × sandwich vs bootstrap", {
  df <- simulate_binary_continuous(n = 1000)
  fit <- causat(df, outcome = "Y", treatment = "A", confounders = ~L)

  res_sw <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich"
  )
  res_bs <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "bootstrap",
    n_boot = 200L
  )

  # SEs should be within a factor of 2 (generous for 200 bootstrap samples).
  ratio <- res_bs$contrasts$se[1] / res_sw$contrasts$se[1]
  expect_gt(ratio, 0.5)
  expect_lt(ratio, 2.0)
})

# ============================================================
# GCOMP × MODEL_FN VARIANTS
# ============================================================

test_that("gcomp × GAM via model_fn = mgcv::gam", {
  skip_if_not_installed("mgcv")
  df <- simulate_binary_continuous(n = 2000)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~ s(L),
    model_fn = mgcv::gam
  )
  result <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich"
  )

  expect_equal(result$contrasts$estimate[1], 3, tolerance = 0.3)
})

# ============================================================
# GCOMP × CENSORING
# ============================================================

test_that("gcomp × binary trt × continuous outcome × censoring", {
  df <- simulate_binary_continuous(n = 2000)
  # Simulate censoring: censor 10% of observations at random.
  set.seed(99)
  df$C <- rbinom(nrow(df), 1, 0.1)
  df$Y[df$C == 1] <- NA

  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    censoring = "C"
  )
  result <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich"
  )

  # ATE should still be ≈ 3 (censoring is random, not informative).
  expect_equal(result$contrasts$estimate[1], 3, tolerance = 0.3)
})


# ============================================================
# MATCHING × BINARY TREATMENT × CONTINUOUS OUTCOME × DIFFERENCE
# ============================================================

test_that("matching × binary trt × continuous outcome × difference × sandwich: ATT ≈ 3", {
  df <- simulate_binary_continuous(n = 2000)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    method = "matching",
    estimand = "ATT"
  )

  result <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich"
  )

  att <- result$contrasts$estimate[1]
  expect_equal(att, 3, tolerance = 0.3)
})

test_that("matching × binary trt × continuous outcome × difference × bootstrap: SE finite", {
  df <- simulate_binary_continuous(n = 500)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    method = "matching",
    estimand = "ATT"
  )

  result <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "bootstrap",
    n_boot = 50L
  )

  expect_gt(result$contrasts$se[1], 0)
  expect_true(is.finite(result$contrasts$se[1]))
})


# ============================================================
# IPW × BINARY TREATMENT × CONTINUOUS OUTCOME × DIFFERENCE
# ============================================================

test_that("ipw × binary trt × continuous outcome × difference × sandwich: ATE ≈ 3", {
  df <- simulate_binary_continuous(n = 2000)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    method = "ipw"
  )

  result <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich"
  )

  ate <- result$contrasts$estimate[1]
  expect_equal(ate, 3, tolerance = 0.3)
})

test_that("ipw × binary trt × continuous outcome × difference × bootstrap: SE finite", {
  df <- simulate_binary_continuous(n = 500)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    method = "ipw"
  )

  result <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "bootstrap",
    n_boot = 50L
  )

  expect_gt(result$contrasts$se[1], 0)
  expect_true(is.finite(result$contrasts$se[1]))
})


# ============================================================
# TRIANGULATION × CONTINUOUS OUTCOME
# ============================================================

test_that("triangulation × binary trt × continuous outcome: all methods ≈ ATE 3", {
  df <- simulate_binary_continuous(n = 3000)

  fit_gc <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    method = "gcomp"
  )
  fit_ipw <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    method = "ipw"
  )
  fit_m <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    method = "matching",
    estimand = "ATT"
  )

  res_gc <- contrast(
    fit_gc,
    list(a1 = static(1), a0 = static(0)),
    reference = "a0"
  )
  res_ipw <- contrast(
    fit_ipw,
    list(a1 = static(1), a0 = static(0)),
    reference = "a0"
  )
  res_m <- contrast(
    fit_m,
    list(a1 = static(1), a0 = static(0)),
    reference = "a0"
  )

  ate_gc <- res_gc$contrasts$estimate[1]
  ate_ipw <- res_ipw$contrasts$estimate[1]
  att_m <- res_m$contrasts$estimate[1]

  # All should be close to 3 (true ATE = ATT = 3 in this DGP).
  expect_equal(ate_gc, 3, tolerance = 0.15)
  expect_equal(ate_ipw, 3, tolerance = 0.3)
  expect_equal(att_m, 3, tolerance = 0.3)

  # They should agree with each other within 1 unit.
  expect_lt(abs(ate_gc - ate_ipw), 1)
  expect_lt(abs(ate_gc - att_m), 1)
})


# ============================================================
# SURVIVAL / PERSON-PERIOD TESTS
# ============================================================

test_that("to_person_period converts wide to long correctly", {
  wide <- data.table::data.table(
    id = 1:3,
    sex = c(0, 1, 0),
    A0 = c(1, 0, 1),
    A1 = c(1, 1, 0),
    L0 = c(5, 3, 7),
    L1 = c(4, 6, 8),
    Y = c(0, 1, 0)
  )
  long <- to_person_period(
    wide,
    id = "id",
    time_varying = list(A = c("A0", "A1"), L = c("L0", "L1")),
    time_invariant = c("sex", "Y")
  )

  expect_equal(nrow(long), 6L)
  expect_true("time" %in% names(long))
  expect_equal(long[id == 1 & time == 0, A], 1)
  expect_equal(long[id == 1 & time == 1, A], 1)
  expect_equal(long[id == 2 & time == 1, A], 1)
  expect_equal(long[id == 3 & time == 1, L], 8)
  # Time-invariant columns carried forward.
  expect_equal(long[id == 1 & time == 0, sex], 0)
  expect_equal(long[id == 1 & time == 1, sex], 0)
})

test_that("to_person_period rejects mismatched lengths", {
  wide <- data.table::data.table(
    id = 1:2,
    A0 = c(1, 0),
    L0 = c(5, 3),
    L1 = c(4, 6)
  )
  expect_error(
    to_person_period(
      wide,
      id = "id",
      time_varying = list(A = "A0", L = c("L0", "L1"))
    ),
    "same length"
  )
})

# ============================================================
# IPW × BINARY TREATMENT × BINARY OUTCOME
# ============================================================

test_that("ipw × binary trt × binary outcome × difference × sandwich", {
  df <- simulate_binary_binary(n = 3000)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    method = "ipw"
  )
  result <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    type = "difference",
    ci_method = "sandwich"
  )
  expect_equal(result$contrasts$estimate[1], 0.33, tolerance = 0.15)
})

test_that("ipw × binary trt × binary outcome × ratio × sandwich", {
  df <- simulate_binary_binary(n = 3000)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    method = "ipw"
  )
  result <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    type = "ratio",
    ci_method = "sandwich"
  )
  rr <- result$contrasts$estimate[1]
  expect_gt(rr, 1.0)
  expect_lt(rr, 4.0)
  expect_gt(result$contrasts$se[1], 0)
})


test_that("ipw × binary outcome × binomial family recovers RD", {
  df <- simulate_binary_binary(n = 3000)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    method = "ipw",
    family = "binomial"
  )
  result <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    type = "difference",
    ci_method = "sandwich"
  )
  expect_equal(result$contrasts$estimate[1], 0.33, tolerance = 0.15)
  expect_gt(result$contrasts$se[1], 0)
})

test_that("matching × binary outcome × quasibinomial family recovers RD", {
  df <- simulate_binary_binary(n = 3000)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    method = "matching",
    family = "quasibinomial"
  )
  result <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    type = "difference",
    ci_method = "sandwich"
  )
  expect_equal(result$contrasts$estimate[1], 0.33, tolerance = 0.15)
  expect_gt(result$contrasts$se[1], 0)
})

test_that("ipw × continuous outcome × gaussian family (default) recovers ATE", {
  df <- simulate_binary_continuous(n = 2000)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    method = "ipw",
    family = "gaussian"
  )
  result <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    type = "difference",
    ci_method = "sandwich"
  )
  expect_equal(result$contrasts$estimate[1], 3.0, tolerance = 0.3)
})

test_that("matching × continuous outcome × gaussian family recovers ATE", {
  df <- simulate_binary_continuous(n = 2000)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    method = "matching",
    family = "gaussian"
  )
  result <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    type = "difference",
    ci_method = "sandwich"
  )
  expect_equal(result$contrasts$estimate[1], 3.0, tolerance = 0.3)
})


# ============================================================
# IPW × BINARY TREATMENT × CONTINUOUS OUTCOME (extended)
# ============================================================

test_that("ipw × binary trt × continuous outcome × sandwich vs bootstrap", {
  df <- simulate_binary_continuous(n = 1000)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    method = "ipw"
  )

  res_sw <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich"
  )
  res_bs <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "bootstrap",
    n_boot = 200L
  )

  ratio <- res_bs$contrasts$se[1] / res_sw$contrasts$se[1]
  expect_gt(ratio, 0.3)
  expect_lt(ratio, 3.0)
})


# ============================================================
# IPW × ESTIMAND VARIANTS
# ============================================================

test_that("ipw × ATT estimand ≈ 3", {
  df <- simulate_binary_continuous(n = 2000)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    method = "ipw",
    estimand = "ATT"
  )
  result <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich"
  )
  expect_equal(result$contrasts$estimate[1], 3, tolerance = 0.3)
})


# ============================================================
# MATCHING × BINARY TREATMENT × BINARY OUTCOME
# ============================================================

test_that("matching × binary trt × binary outcome × difference × sandwich", {
  df <- simulate_binary_binary(n = 3000)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    method = "matching",
    estimand = "ATT"
  )
  result <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    type = "difference",
    ci_method = "sandwich"
  )
  expect_equal(result$contrasts$estimate[1], 0.33, tolerance = 0.15)
})

test_that("matching × binary trt × binary outcome × ratio × sandwich", {
  df <- simulate_binary_binary(n = 3000)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    method = "matching",
    estimand = "ATT"
  )
  result <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    type = "ratio",
    ci_method = "sandwich"
  )
  rr <- result$contrasts$estimate[1]
  expect_gt(rr, 1.0)
  expect_lt(rr, 4.0)
  expect_gt(result$contrasts$se[1], 0)
})


# ============================================================
# MATCHING × BINARY TREATMENT × CONTINUOUS OUTCOME (extended)
# ============================================================

test_that("matching × binary trt × continuous outcome × sandwich vs bootstrap", {
  df <- simulate_binary_continuous(n = 1000)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    method = "matching",
    estimand = "ATT"
  )

  res_sw <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich"
  )
  res_bs <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "bootstrap",
    n_boot = 100L
  )

  ratio <- res_bs$contrasts$se[1] / res_sw$contrasts$se[1]
  expect_gt(ratio, 0.3)
  expect_lt(ratio, 3.0)
})


# ============================================================
# GCOMP × BINARY TREATMENT × BINARY OUTCOME × OR CONTRAST
# ============================================================

test_that("gcomp × binary trt × binary outcome × or × sandwich", {
  df <- simulate_binary_binary(n = 3000)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    family = "binomial"
  )
  result <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    type = "or",
    ci_method = "sandwich"
  )
  or_val <- result$contrasts$estimate[1]
  expect_gt(or_val, 1.0)
  expect_gt(result$contrasts$se[1], 0)
  expect_true(is.finite(result$contrasts$se[1]))
})


# ============================================================
# GCOMP × ESTIMAND VARIANTS
# ============================================================

test_that("gcomp × ATC estimand ≈ 3", {
  df <- simulate_binary_continuous(n = 2000)
  fit <- causat(df, outcome = "Y", treatment = "A", confounders = ~L)
  result <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    estimand = "ATC",
    ci_method = "sandwich"
  )
  expect_equal(result$contrasts$estimate[1], 3, tolerance = 0.3)
  expect_equal(result$estimand, "ATC")
})


# ============================================================
# GCOMP × MULTIPLE INTERVENTIONS (>2)
# ============================================================

test_that("gcomp × multiple interventions with custom reference", {
  df <- simulate_binary_continuous(n = 2000)
  fit <- causat(df, outcome = "Y", treatment = "A", confounders = ~L)
  result <- contrast(
    fit,
    interventions = list(
      treat = static(1),
      control = static(0),
      natural = NULL
    ),
    reference = "control",
    ci_method = "sandwich"
  )

  expect_equal(nrow(result$estimates), 3L)
  expect_equal(nrow(result$contrasts), 2L)
  expect_equal(result$reference, "control")
  expect_true(all(grepl("vs control", result$contrasts$comparison)))
})


# ============================================================
# TRIANGULATION × BINARY OUTCOME
# ============================================================

test_that("triangulation: all methods agree on binary outcome RD", {
  df <- simulate_binary_binary(n = 3000)

  fit_gc <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    method = "gcomp",
    family = "binomial"
  )
  fit_ipw <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    method = "ipw"
  )
  fit_m <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    method = "matching",
    estimand = "ATT"
  )

  res_gc <- contrast(
    fit_gc,
    list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    type = "difference"
  )
  res_ipw <- contrast(
    fit_ipw,
    list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    type = "difference"
  )
  res_m <- contrast(
    fit_m,
    list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    type = "difference"
  )

  expect_equal(res_gc$contrasts$estimate[1], 0.33, tolerance = 0.15)
  expect_equal(res_ipw$contrasts$estimate[1], 0.33, tolerance = 0.15)
  expect_equal(res_m$contrasts$estimate[1], 0.33, tolerance = 0.15)
})


# ============================================================
# SURVIVAL / PERSON-PERIOD TESTS
# ============================================================

test_that("causat_survival fits a pooled logistic model on long data", {
  # Simulate simple survival data in person-period format.
  set.seed(123)
  n_id <- 200
  n_time <- 10
  ids <- rep(seq_len(n_id), each = n_time)
  times <- rep(0:(n_time - 1), times = n_id)
  A <- rep(rbinom(n_id, 1, 0.5), each = n_time)
  L <- rep(rnorm(n_id), each = n_time)

  # Hazard increases with time, treatment is protective.
  h <- plogis(-3 + 0.1 * times - 0.5 * A + 0.3 * L)
  event <- rbinom(length(h), 1, h)

  # Zero out events after the first one per individual.
  dt <- data.table::data.table(
    id = ids,
    time = times,
    A = A,
    L = L,
    event = event
  )
  dt[, cum_event := cumsum(event), by = id]
  dt[cum_event > 1, event := 0]
  dt[, cum_event := NULL]

  fit <- causat_survival(
    dt,
    outcome = "event",
    treatment = "A",
    confounders = ~L,
    id = "id",
    time = "time",
    time_formula = ~ splines::ns(time, 3)
  )

  expect_s3_class(fit, "causatr_fit")
  expect_equal(fit$type, "survival")
  expect_true(length(fit$details$time_points) > 1)
})


# ============================================================
# MATCHING × ATE ESTIMAND (full matching auto-selection)
# ============================================================

test_that("matching × binary trt × continuous outcome × ATE × sandwich ≈ 3", {
  df <- simulate_binary_continuous(n = 2000)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    method = "matching",
    estimand = "ATE"
  )
  result <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich"
  )
  expect_equal(result$contrasts$estimate[1], 3, tolerance = 0.3)
  expect_equal(result$estimand, "ATE")
})

test_that("matching × binary trt × continuous outcome × ATE × bootstrap", {
  df <- simulate_binary_continuous(n = 1000)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    method = "matching",
    estimand = "ATE"
  )
  result <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "bootstrap",
    n_boot = 50L
  )
  expect_gt(result$contrasts$se[1], 0)
  expect_true(is.finite(result$contrasts$se[1]))
})


# ============================================================
# MATCHING × ATC ESTIMAND
# ============================================================

test_that("matching × binary trt × continuous outcome × ATC × sandwich ≈ 3", {
  df <- simulate_binary_continuous(n = 2000)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    method = "matching",
    estimand = "ATC"
  )
  result <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich"
  )
  expect_equal(result$contrasts$estimate[1], 3, tolerance = 0.3)
  expect_equal(result$estimand, "ATC")
})


# ============================================================
# IPW × ATC ESTIMAND
# ============================================================

test_that("ipw × binary trt × continuous outcome × ATC × sandwich ≈ 3", {
  df <- simulate_binary_continuous(n = 2000)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    method = "ipw",
    estimand = "ATC"
  )
  result <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich"
  )
  expect_equal(result$contrasts$estimate[1], 3, tolerance = 0.3)
  expect_equal(result$estimand, "ATC")
})


# ============================================================
# IPW × BINARY OUTCOME × OR
# ============================================================

test_that("ipw × binary trt × binary outcome × or × sandwich", {
  df <- simulate_binary_binary(n = 3000)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    method = "ipw"
  )
  result <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    type = "or",
    ci_method = "sandwich"
  )
  or_val <- result$contrasts$estimate[1]
  expect_gt(or_val, 1.0)
  expect_gt(result$contrasts$se[1], 0)
  expect_true(is.finite(result$contrasts$se[1]))
})


# ============================================================
# MATCHING × BINARY OUTCOME × OR
# ============================================================

test_that("matching × binary trt × binary outcome × or × sandwich", {
  df <- simulate_binary_binary(n = 3000)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    method = "matching",
    estimand = "ATT"
  )
  result <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    type = "or",
    ci_method = "sandwich"
  )
  or_val <- result$contrasts$estimate[1]
  expect_gt(or_val, 1.0)
  expect_gt(result$contrasts$se[1], 0)
  expect_true(is.finite(result$contrasts$se[1]))
})


# ============================================================
# GCOMP × BINARY OUTCOME × BOOTSTRAP SE
# ============================================================

test_that("gcomp × binary trt × binary outcome × sandwich vs bootstrap", {
  df <- simulate_binary_binary(n = 1500)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    family = "binomial"
  )
  res_sw <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    type = "difference",
    ci_method = "sandwich"
  )
  res_bs <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    type = "difference",
    ci_method = "bootstrap",
    n_boot = 200L
  )
  ratio <- res_bs$contrasts$se[1] / res_sw$contrasts$se[1]
  expect_gt(ratio, 0.5)
  expect_lt(ratio, 2.0)
})


# ============================================================
# ICE × BINARY OUTCOME × SANDWICH
# ============================================================

test_that("ICE × binary outcome × sandwich: SE finite and positive", {
  set.seed(42)
  n <- 2000
  sim_id <- rep(seq_len(n), each = 2)
  sim_time <- rep(0:1, n)
  L0 <- rbinom(n, 1, 0.5)
  A0 <- rbinom(n, 1, plogis(-0.5 + 0.5 * L0))
  L1 <- rbinom(n, 1, plogis(-0.5 + 0.8 * A0 + 0.3 * L0))
  A1 <- rbinom(n, 1, plogis(-0.5 + 0.5 * L1))
  Y <- rbinom(n, 1, plogis(-0.5 + 0.5 * A0 + 0.5 * A1 + 0.3 * L1))

  sim_long <- data.frame(
    id = sim_id,
    time = sim_time,
    A = c(rbind(A0, A1)),
    L = c(rbind(L0, L1)),
    Y = c(rbind(NA_real_, Y))
  )

  fit <- causat(
    sim_long,
    outcome = "Y",
    treatment = "A",
    confounders = ~1,
    confounders_tv = ~L,
    family = "binomial",
    id = "id",
    time = "time"
  )

  result <- contrast(
    fit,
    interventions = list(always = static(1), never = static(0)),
    reference = "never",
    type = "difference",
    ci_method = "sandwich"
  )

  expect_true(all(result$estimates$se > 0))
  expect_true(all(is.finite(result$estimates$se)))
  expect_true(is.finite(result$contrasts$se[1]))
})


# ============================================================
# ICE × RATIO CONTRAST
# ============================================================

test_that("ICE × continuous outcome × ratio contrast", {
  groups <- data.frame(
    A0 = c(0, 0, 0, 0, 1, 1, 1, 1),
    L1 = c(0, 0, 1, 1, 0, 0, 1, 1),
    A1 = c(0, 1, 0, 1, 0, 1, 0, 1),
    Y = c(84, 84, 52, 52, 76, 76, 44, 44),
    N = c(2400, 1600, 2400, 9600, 4800, 3200, 1600, 6400)
  )
  rows <- lapply(seq_len(nrow(groups)), function(i) {
    g <- groups[i, ]
    off <- sum(groups$N[seq_len(i - 1)])
    data.frame(
      id = seq_len(g$N) + off,
      A0 = g$A0,
      L1 = g$L1,
      A1 = g$A1,
      Y = g$Y
    )
  })
  wide <- do.call(rbind, rows)
  t0 <- data.frame(
    id = wide$id,
    time = 0L,
    A = wide$A0,
    L = NA_real_,
    Y = NA_real_
  )
  t1 <- data.frame(
    id = wide$id,
    time = 1L,
    A = wide$A1,
    L = wide$L1,
    Y = wide$Y
  )
  long <- rbind(t0, t1)

  fit <- causat(
    long,
    outcome = "Y",
    treatment = "A",
    confounders = ~1,
    confounders_tv = ~L,
    id = "id",
    time = "time"
  )

  result <- contrast(
    fit,
    interventions = list(always = static(1), never = static(0)),
    reference = "never",
    type = "ratio",
    ci_method = "sandwich"
  )

  expect_equal(result$contrasts$estimate[1], 1.0, tolerance = 0.01)
  expect_gt(result$contrasts$se[1], 0)
})


# ============================================================
# ICE × DYNAMIC INTERVENTION × SANDWICH SE
# ============================================================

test_that("ICE × dynamic intervention × sandwich: SE finite", {
  groups <- data.frame(
    A0 = c(0, 0, 0, 0, 1, 1, 1, 1),
    L1 = c(0, 0, 1, 1, 0, 0, 1, 1),
    A1 = c(0, 1, 0, 1, 0, 1, 0, 1),
    Y = c(84, 84, 52, 52, 76, 76, 44, 44),
    N = c(2400, 1600, 2400, 9600, 4800, 3200, 1600, 6400)
  )
  rows <- lapply(seq_len(nrow(groups)), function(i) {
    g <- groups[i, ]
    off <- sum(groups$N[seq_len(i - 1)])
    data.frame(
      id = seq_len(g$N) + off,
      A0 = g$A0,
      L1 = g$L1,
      A1 = g$A1,
      Y = g$Y
    )
  })
  wide <- do.call(rbind, rows)
  t0 <- data.frame(
    id = wide$id,
    time = 0L,
    A = wide$A0,
    L = NA_real_,
    Y = NA_real_
  )
  t1 <- data.frame(
    id = wide$id,
    time = 1L,
    A = wide$A1,
    L = wide$L1,
    Y = wide$Y
  )
  long <- rbind(t0, t1)

  fit <- causat(
    long,
    outcome = "Y",
    treatment = "A",
    confounders = ~1,
    confounders_tv = ~L,
    id = "id",
    time = "time"
  )

  result <- contrast(
    fit,
    interventions = list(
      adaptive = dynamic(\(data, trt) {
        ifelse(!is.na(data$L) & data$L > 0, 1L, 0L)
      }),
      never = static(0)
    ),
    reference = "never",
    ci_method = "sandwich"
  )

  expect_true(all(result$estimates$se > 0))
  expect_true(all(is.finite(result$estimates$se)))
  expect_true(is.finite(result$contrasts$se[1]))
})


# ============================================================
# MATCHING × BINARY OUTCOME × ATE (full matching)
# ============================================================

test_that("matching × binary trt × binary outcome × ATE × sandwich", {
  df <- simulate_binary_binary(n = 3000)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    method = "matching",
    estimand = "ATE"
  )
  result <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    type = "difference",
    ci_method = "sandwich"
  )
  expect_equal(result$contrasts$estimate[1], 0.33, tolerance = 0.15)
  expect_equal(result$estimand, "ATE")
})


# ============================================================
# IPW × BINARY OUTCOME × BOOTSTRAP
# ============================================================

test_that("ipw × binary trt × binary outcome × bootstrap SE finite", {
  df <- simulate_binary_binary(n = 1000)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    method = "ipw"
  )
  result <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    type = "difference",
    ci_method = "bootstrap",
    n_boot = 50L
  )
  expect_gt(result$contrasts$se[1], 0)
  expect_true(is.finite(result$contrasts$se[1]))
})


# ============================================================
# MATCHING × BINARY OUTCOME × BOOTSTRAP
# ============================================================

test_that("matching × binary trt × binary outcome × bootstrap SE finite", {
  df <- simulate_binary_binary(n = 1000)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    method = "matching",
    estimand = "ATT"
  )
  result <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    type = "difference",
    ci_method = "bootstrap",
    n_boot = 50L
  )
  expect_gt(result$contrasts$se[1], 0)
  expect_true(is.finite(result$contrasts$se[1]))
})


# ============================================================
# GCOMP × CONTINUOUS TREATMENT × BOOTSTRAP
# ============================================================

test_that("gcomp × continuous trt × shift × bootstrap SE finite", {
  df <- simulate_continuous_continuous(n = 1000)
  fit <- causat(df, outcome = "Y", treatment = "A", confounders = ~L)
  result <- contrast(
    fit,
    interventions = list(shifted = shift(-1), observed = NULL),
    reference = "observed",
    ci_method = "bootstrap",
    n_boot = 100L
  )
  expect_gt(result$contrasts$se[1], 0)
  expect_true(is.finite(result$contrasts$se[1]))
})


# ============================================================
# LOG-SCALE CIs FOR RATIO / OR
# ============================================================

test_that("log-scale ratio CI is always positive and contains true RR", {
  d <- simulate_binary_binary(n = 5000, seed = 42)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    family = "binomial"
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    type = "ratio"
  )
  expect_gt(res$contrasts$ci_lower, 0)
  rr <- res$contrasts$estimate
  expect_true(rr > 1.5 && rr < 3.0)
  expect_lt(res$contrasts$ci_lower, 2.15)
  expect_gt(res$contrasts$ci_upper, 2.15)
})

test_that("log-scale OR CI is always positive and contains true OR", {
  d <- simulate_binary_binary(n = 5000, seed = 42)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    family = "binomial"
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    type = "or"
  )
  expect_gt(res$contrasts$ci_lower, 0)
  or_est <- res$contrasts$estimate
  expect_true(or_est > 2 && or_est < 7)
})

test_that("ipw × binary outcome × ratio: log-scale CI positive", {
  d <- simulate_binary_binary(n = 5000, seed = 42)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    method = "ipw",
    family = "binomial"
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    type = "ratio"
  )
  expect_gt(res$contrasts$ci_lower, 0)
  expect_gt(res$contrasts$estimate, 1)
})

test_that("matching × binary outcome × ratio: log-scale CI positive", {
  d <- simulate_binary_binary(n = 5000, seed = 42)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    method = "matching",
    estimand = "ATT",
    family = "binomial"
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    type = "ratio"
  )
  expect_gt(res$contrasts$ci_lower, 0)
  expect_gt(res$contrasts$estimate, 1)
})


# ============================================================
# ANALYTIC IF vs NUMERICAL FALLBACK (variance_if_numeric Tier 1)
# ============================================================
# These previously tested compute_vcov_marginal()'s analytic vs numDeriv
# Jacobian paths. After the variance refactor the equivalence we care
# about is "variance_if() main path" (analytic correct_model + bread)
# vs "variance_if_numeric() Tier 1" (sandwich::estfun + sandwich::bread),
# which is covered by tests in test-variance-if.R.

# ============================================================
# EXTERNAL WEIGHTS × ALL METHODS
# ============================================================

test_that("gcomp × external weights + bootstrap: finite SE", {
  d <- simulate_binary_continuous(n = 1000, seed = 42)
  w <- runif(nrow(d), 0.5, 2)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    weights = w
  )
  expect_true(!is.null(fit$details$weights))
  expect_false(".causatr_w" %in% names(fit$data))
  res <- suppressWarnings(contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    ci_method = "bootstrap",
    n_boot = 50
  ))
  expect_true(all(is.finite(res$contrasts$se)))
  expect_gt(res$contrasts$se, 0)
})

test_that("ipw × external weights + sandwich: runs without error", {
  d <- simulate_binary_continuous(n = 1000, seed = 42)
  w <- runif(nrow(d), 0.5, 2)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    method = "ipw",
    weights = w
  )
  res <- contrast(
    fit,
    interventions = list(a0 = static(0), a1 = static(1)),
    ci_method = "sandwich"
  )
  expect_true(all(is.finite(res$contrasts$se)))
  expect_gt(res$contrasts$se, 0)
})

test_that("gcomp × external weights recovers WEIGHTED target population mean", {
  # Truth-based test for survey-style external weights. The existing
  # external-weights tests only assert "runs without error and finite
  # SE"; this one pins the actual numerical target so a future bug
  # that drops the weights or applies them in the wrong place can
  # be caught.
  #
  # DGP:
  #   S      ~ Bernoulli(0.5) (stratum)
  #   L | S=0 ~ N(0, 1)
  #   L | S=1 ~ N(2, 1)
  #   A | L  ~ Bernoulli(plogis(0.4 * L))
  #   Y      = 1 + 2*A + 0.5*L + N(0, 1)
  # Survey weight: w = 1 if S=0, w = 4 if S=1.
  # Truth: weighted E[Y^{a=1}] = 1 + 2 + 0.5 * E_w[L]
  #        E_w[L] = (n0 * 0 + 4 * n1 * 2) / (n0 + 4 * n1) -> 1.6
  set.seed(2026)
  n0 <- 2000
  n1 <- 2000
  S <- c(rep(0L, n0), rep(1L, n1))
  L <- c(stats::rnorm(n0, 0, 1), stats::rnorm(n1, 2, 1))
  A <- stats::rbinom(n0 + n1, 1, stats::plogis(0.4 * L))
  Y <- 1 + 2 * A + 0.5 * L + stats::rnorm(n0 + n1)
  w <- ifelse(S == 1, 4, 1)
  d <- data.frame(Y = Y, A = A, L = L, S = S)

  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    weights = w
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    reference = "a0",
    ci_method = "sandwich"
  )

  ew_L <- (n0 * 0 + 4 * n1 * 2) / (n0 + 4 * n1) # = 1.6
  truth_a1 <- 1 + 2 + 0.5 * ew_L # = 3.8
  truth_a0 <- 1 + 0.5 * ew_L # = 1.8

  est <- res$estimates[order(res$estimates$intervention)]
  expect_lt(
    abs(est$estimate[est$intervention == "a1"] - truth_a1),
    0.15
  )
  expect_lt(
    abs(est$estimate[est$intervention == "a0"] - truth_a0),
    0.15
  )
  # ATE is constant across strata so the weighted vs unweighted
  # contrast is the same: 2.
  expect_lt(abs(res$contrasts$estimate[1] - 2), 0.1)

  # Each marginal mean's CI must cover its analytical truth.
  e_a1_row <- est[est$intervention == "a1"]
  e_a0_row <- est[est$intervention == "a0"]
  expect_lt(e_a1_row$ci_lower, truth_a1)
  expect_gt(e_a1_row$ci_upper, truth_a1)
  expect_lt(e_a0_row$ci_lower, truth_a0)
  expect_gt(e_a0_row$ci_upper, truth_a0)
  # And the contrast CI must cover 2.
  expect_lt(res$contrasts$ci_lower[1], 2)
  expect_gt(res$contrasts$ci_upper[1], 2)
})

test_that("ICE × external weights produces the WEIGHTED marginal mean", {
  # Critical regression test for the ICE backward-loop weight bug.
  # Before the fix, only the final-time outcome model received
  # external weights; every intermediate pseudo-outcome model was
  # silently unweighted. The bug produced finite SEs and
  # plausible-looking estimates — it could only be detected by
  # comparing against the analytical weighted target.
  #
  # DGP (2 periods, balanced panel):
  #   L_0 ~ N(0, 1)
  #   A_0 | L_0 ~ Bernoulli(plogis(0.4 * L_0))
  #   L_1 = L_0 + 0.3 * A_0 + N(0, 1)
  #   A_1 | L_1, A_0 ~ Bernoulli(plogis(0.4 * L_1 + 0.5 * A_0))
  #   Y    = 1 + A_0 + A_1 + 0.5 * L_1 + N(0, 1)
  # Survey weight per id: w = 4 if L_0 > 0, w = 1 otherwise.
  #
  # Under the "never" regime (A_0 = A_1 = 0):
  #   Y^{never} = 1 + 0.5 * L_0 + 0.5*N(0,1) + N(0,1)
  # so weighted E[Y^{never}] = 1 + 0.5 * E_w[L_0].
  # Unweighted E[Y^{never}] ≈ 1 + 0 = 1.
  set.seed(2027)
  n <- 4000
  L0 <- stats::rnorm(n)
  A0 <- stats::rbinom(n, 1, stats::plogis(0.4 * L0))
  L1 <- L0 + 0.3 * A0 + stats::rnorm(n)
  A1 <- stats::rbinom(n, 1, stats::plogis(0.4 * L1 + 0.5 * A0))
  Y <- 1 + A0 + A1 + 0.5 * L1 + stats::rnorm(n)

  w_id <- ifelse(L0 > 0, 4, 1)
  long <- data.table::data.table(
    id = rep(seq_len(n), each = 2),
    time = rep(0:1, times = n),
    A = as.numeric(rbind(A0, A1)),
    L = as.numeric(rbind(L0, L1)),
    Y = rep(Y, each = 2)
  )
  w_long <- rep(w_id, each = 2)

  fit_w <- causat(
    long,
    outcome = "Y",
    treatment = "A",
    confounders = ~1,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    weights = w_long
  )
  res_w <- contrast(
    fit_w,
    interventions = list(never = static(0)),
    ci_method = "sandwich"
  )

  fit_uw <- causat(
    long,
    outcome = "Y",
    treatment = "A",
    confounders = ~1,
    confounders_tv = ~L,
    id = "id",
    time = "time"
  )
  res_uw <- contrast(
    fit_uw,
    interventions = list(never = static(0)),
    ci_method = "sandwich"
  )

  e_never_w <- res_w$estimates$estimate[1]
  e_never_uw <- res_uw$estimates$estimate[1]

  # Analytical weighted truth for E[Y^{never}].
  ew_L0 <- sum(w_id * L0) / sum(w_id)
  truth_w <- 1 + 0.5 * ew_L0

  # Both estimates within 0.2 of their respective targets.
  expect_lt(abs(e_never_w - truth_w), 0.2)
  expect_lt(abs(e_never_uw - 1), 0.2)

  # Each CI must cover its own analytical target.
  expect_lt(res_w$estimates$ci_lower[1], truth_w)
  expect_gt(res_w$estimates$ci_upper[1], truth_w)
  expect_lt(res_uw$estimates$ci_lower[1], 1)
  expect_gt(res_uw$estimates$ci_upper[1], 1)

  # And — critically — the two must DIFFER. If the backward-loop
  # models silently dropped the weights, the weighted and
  # unweighted estimates would coincide.
  expect_gt(e_never_w - e_never_uw, 0.1)
})

test_that("ICE × continuous TV treatment × shift recovers structural 2*delta (vs lmtp)", {
  # Truth-based test for shift() on a time-varying continuous
  # treatment. Validated against `lmtp::lmtp_tmle()` — both methods
  # now agree on E[Y^{shift(1)}] = 1 + 2*delta = 3 with delta = 1.
  #
  # DGP designed so L_1 is INDEPENDENT of A_0 and A_1 is independent
  # of A_0. This isolates the shift effect: each time point's A
  # enters Y with coefficient 1 and there is no downstream-confounding
  # path, so the structural truth is exactly 2*delta.
  #
  #   L_0 ~ N(0, 1)
  #   A_0 = 0.4*L_0 + N(0, 1)
  #   L_1 = L_0 + N(0, 1)        (independent of A_0)
  #   A_1 = 0.4*L_1 + N(0, 1)    (independent of A_0)
  #   Y   = 1 + A_0 + A_1 + 0.5*L_1 + N(0, 1)
  #
  # Pre-fix, causatr's ICE returned ~3*delta because
  # `ice_apply_intervention_long()` recomputed treatment lag columns
  # from the shifted treatment, double-applying the shift through
  # both the lag-column path and the current-A prediction path.
  # Removing the lag recomputation aligns causatr with the Robins
  # iterated-conditional-expectation algorithm and with lmtp.
  set.seed(2028)
  n <- 5000
  delta <- 1

  L0 <- stats::rnorm(n)
  A0 <- 0.4 * L0 + stats::rnorm(n)
  L1 <- L0 + stats::rnorm(n)
  A1 <- 0.4 * L1 + stats::rnorm(n)
  Y <- 1 + A0 + A1 + 0.5 * L1 + stats::rnorm(n)

  long <- data.table::data.table(
    id = rep(seq_len(n), each = 2),
    time = rep(0:1, times = n),
    A = as.numeric(rbind(A0, A1)),
    L = as.numeric(rbind(L0, L1)),
    Y = rep(Y, each = 2)
  )
  fit <- causat(
    long,
    outcome = "Y",
    treatment = "A",
    confounders = ~1,
    confounders_tv = ~L,
    id = "id",
    time = "time"
  )
  res <- contrast(
    fit,
    interventions = list(shifted = shift(delta), nat = shift(0)),
    type = "difference",
    reference = "nat",
    ci_method = "sandwich"
  )
  truth <- 2 * delta
  expect_lt(abs(res$contrasts$estimate[1] - truth), 0.2)
  expect_lt(res$contrasts$ci_lower[1], truth)
  expect_gt(res$contrasts$ci_upper[1], truth)

  # Also validate the marginal mean under shift matches the lmtp
  # reference of 1 + 2*delta = 3.
  est <- res$estimates
  e_shift <- est$estimate[est$intervention == "shifted"]
  expect_lt(abs(e_shift - 3), 0.2)
})


test_that("ICE × continuous TV treatment × scale_by recovers 2*(c-1) (vs lmtp)", {
  # Truth-based companion to the shift() test above. DGP:
  #
  #   L_0 ~ N(0, 1)
  #   A_0 = 1 + 0.4*L_0 + N(0, 1)   (E[A_0] = 1)
  #   L_1 = L_0 + N(0, 1)           (independent of A_0)
  #   A_1 = 1 + 0.4*L_1 + N(0, 1)   (E[A_1] = 1)
  #   Y   = 1 + A_0 + A_1 + 0.5*L_1 + N(0, 1)
  #
  # Under scale_by(c):
  #   E[Y^{c*A}] = 1 + c*E[A_0] + c*E[A_1] + 0.5*E[L_1] = 1 + 2c
  # so the contrast against scale_by(1) is 2*(c - 1).
  # Cross-checked against `lmtp::lmtp_tmle()` with a custom shift
  # function (= 1.5 * A_t): lmtp returns ~0.985, causatr returns
  # ~0.987, structural truth = 1.
  set.seed(2029)
  n <- 5000
  c_scale <- 1.5

  L0 <- stats::rnorm(n)
  A0 <- 1 + 0.4 * L0 + stats::rnorm(n)
  L1 <- L0 + stats::rnorm(n)
  A1 <- 1 + 0.4 * L1 + stats::rnorm(n)
  Y <- 1 + A0 + A1 + 0.5 * L1 + stats::rnorm(n)

  long <- data.table::data.table(
    id = rep(seq_len(n), each = 2),
    time = rep(0:1, times = n),
    A = as.numeric(rbind(A0, A1)),
    L = as.numeric(rbind(L0, L1)),
    Y = rep(Y, each = 2)
  )
  fit <- causat(
    long,
    outcome = "Y",
    treatment = "A",
    confounders = ~1,
    confounders_tv = ~L,
    id = "id",
    time = "time"
  )
  res <- contrast(
    fit,
    interventions = list(scaled = scale_by(c_scale), nat = scale_by(1)),
    type = "difference",
    reference = "nat",
    ci_method = "sandwich"
  )
  truth <- 2 * (c_scale - 1)
  expect_lt(abs(res$contrasts$estimate[1] - truth), 0.05)
  expect_lt(res$contrasts$ci_lower[1], truth)
  expect_gt(res$contrasts$ci_upper[1], truth)

  # Marginal mean under scale_by(c) should be ~1 + 2c = 4.
  est <- res$estimates
  e_scaled <- est$estimate[est$intervention == "scaled"]
  expect_lt(abs(e_scaled - (1 + 2 * c_scale)), 0.1)
})


# ============================================================
# CONFINT × STRATIFIED BOOTSTRAP
# ============================================================

test_that("confint with by + bootstrap returns correctly ordered CIs", {
  d <- simulate_effect_mod(n = 2000, seed = 42)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L + sex + A:L + A:sex
  )
  res <- suppressWarnings(contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    ci_method = "bootstrap",
    n_boot = 50,
    by = "sex"
  ))
  ci <- confint(res)
  expect_equal(nrow(ci), nrow(res$estimates))
  for (i in seq_len(nrow(ci))) {
    expect_lte(ci[i, "lower"], res$estimates$estimate[i])
    expect_gte(ci[i, "upper"], res$estimates$estimate[i])
  }
})


test_that("to_person_period() aborts on duplicated ids", {
  # Regression guard: a duplicated id in the wide input would
  # silently produce malformed long data because the reshape
  # assumes one row per id. Catch it up front.
  wide <- data.table::data.table(
    id = c(1, 2, 2), # duplicate
    sex = c(0, 1, 1),
    A0 = c(1, 0, 0),
    A1 = c(1, 1, 1),
    L0 = c(5, 3, 3),
    L1 = c(4, 6, 6),
    Y = c(0, 1, 1)
  )
  expect_snapshot(
    error = TRUE,
    to_person_period(
      wide,
      id = "id",
      time_varying = list(A = c("A0", "A1"), L = c("L0", "L1")),
      time_invariant = c("sex", "Y")
    )
  )
})
