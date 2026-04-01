# Simulation-based tests with known true parameter values.
#
# DGP (data-generating process):
#   L ~ N(0, 1)
#   A | L ~ Bernoulli(expit(0.5 * L))       [binary treatment]
#   Y | A, L ~ N(2 + 3*A + 1.5*L, sd = 1)   [continuous outcome]
#
# True ATE  = E[Y(1)] - E[Y(0)] = 3
# True ATT  = 3 (same because treatment effect is constant in this DGP)
# True E[Y(1)] = 2 + 3 + 1.5*E[L] = 5 (since E[L] = 0)
# True E[Y(0)] = 2 + 1.5*E[L] = 2

simulate_binary_continuous <- function(n = 2000, seed = 42) {
  set.seed(seed)
  L <- rnorm(n)
  ps <- plogis(0.5 * L)
  A <- rbinom(n, 1, ps)
  Y <- 2 + 3 * A + 1.5 * L + rnorm(n)
  data.frame(Y = Y, A = A, L = L)
}

# DGP for binary outcome:
#   L ~ N(0, 1)
#   A | L ~ Bernoulli(expit(0.5 * L))
#   Y | A, L ~ Bernoulli(expit(-1 + 1.5*A + 0.8*L))
#
# True risk under A=1: E[expit(-1 + 1.5 + 0.8*L)] ≈ 0.622 (via simulation)
# True risk under A=0: E[expit(-1 + 0.8*L)]        ≈ 0.289
# True RD ≈ 0.333

simulate_binary_binary <- function(n = 2000, seed = 42) {
  set.seed(seed)
  L <- rnorm(n)
  ps <- plogis(0.5 * L)
  A <- rbinom(n, 1, ps)
  Y <- rbinom(n, 1, plogis(-1 + 1.5 * A + 0.8 * L))
  data.frame(Y = Y, A = A, L = L)
}

# DGP for continuous treatment:
#   L ~ N(0, 1)
#   A | L ~ N(1 + 0.5*L, sd = 1)     [continuous treatment]
#   Y | A, L ~ N(1 + 2*A + L, sd = 1)
#
# True E[Y(a)] = 1 + 2*a + E[L] = 1 + 2*a
# shift(-1): E[Y(A-1)] vs E[Y(A)] → difference = -2

simulate_continuous_continuous <- function(n = 2000, seed = 42) {
  set.seed(seed)
  L <- rnorm(n)
  A <- 1 + 0.5 * L + rnorm(n)
  Y <- 1 + 2 * A + L + rnorm(n)
  data.frame(Y = Y, A = A, L = L)
}


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
  expect_equal(ate, 3, tolerance = 0.3)
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
  expect_equal(att, 3, tolerance = 0.3)
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
  expect_equal(diff, -2, tolerance = 0.3)
})

test_that("gcomp × continuous trt × scale intervention × difference × sandwich", {
  df <- simulate_continuous_continuous(n = 2000)
  fit <- causat(df, outcome = "Y", treatment = "A", confounders = ~L)
  result <- contrast(
    fit,
    interventions = list(halved = scale(0.5), observed = NULL),
    reference = "observed",
    ci_method = "sandwich"
  )

  # scale(0.5) halves A; since E[A] ≈ 1, difference ≈ 2*(0.5*1 - 1) = -1.
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
  expect_equal(result$contrasts$estimate[1], 3, tolerance = 0.5)
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

  expect_equal(result$contrasts$estimate[1], 3, tolerance = 0.5)
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
  expect_equal(result$contrasts$estimate[1], 3, tolerance = 0.5)
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
  expect_equal(att, 3, tolerance = 0.5)
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
  expect_equal(ate, 3, tolerance = 0.5)
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
  expect_equal(ate_gc, 3, tolerance = 0.5)
  expect_equal(ate_ipw, 3, tolerance = 0.5)
  expect_equal(att_m, 3, tolerance = 0.5)

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
  expect_equal(result$contrasts$estimate[1], 3, tolerance = 0.5)
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
  expect_equal(result$contrasts$estimate[1], 3, tolerance = 0.5)
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
# S3 METHODS: coef, confint, print, summary
# ============================================================

test_that("coef.causatr_result returns named numeric vector", {
  df <- simulate_binary_continuous(n = 1000)
  fit <- causat(df, outcome = "Y", treatment = "A", confounders = ~L)
  result <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich"
  )
  co <- coef(result)
  expect_named(co, c("a1", "a0"))
  expect_true(all(is.finite(co)))
})

test_that("confint.causatr_result returns matrix with lower/upper", {
  df <- simulate_binary_continuous(n = 1000)
  fit <- causat(df, outcome = "Y", treatment = "A", confounders = ~L)
  result <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich"
  )
  ci <- confint(result)
  expect_equal(nrow(ci), 2L)
  expect_equal(colnames(ci), c("lower", "upper"))
  expect_true(all(ci[, "lower"] < ci[, "upper"]))
})

test_that("print.causatr_fit outputs method and type", {
  df <- simulate_binary_continuous(n = 100)
  fit <- causat(df, outcome = "Y", treatment = "A", confounders = ~L)
  expect_output(print(fit), "causatr_fit")
  expect_output(print(fit), "gcomp")
})

test_that("print.causatr_result outputs contrast info", {
  df <- simulate_binary_continuous(n = 500)
  fit <- causat(df, outcome = "Y", treatment = "A", confounders = ~L)
  result <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0"
  )
  expect_output(print(result), "causatr_result")
  expect_output(print(result), "sandwich")
})

test_that("summary.causatr_fit outputs model details", {
  df <- simulate_binary_continuous(n = 100)
  fit <- causat(df, outcome = "Y", treatment = "A", confounders = ~L)
  expect_output(summary(fit), "causatr fit")
  expect_output(summary(fit), "Confounders")
})

test_that("summary.causatr_result outputs result details", {
  df <- simulate_binary_continuous(n = 500)
  fit <- causat(df, outcome = "Y", treatment = "A", confounders = ~L)
  result <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0"
  )
  expect_output(summary(result), "causatr result")
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
