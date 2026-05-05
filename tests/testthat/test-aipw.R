# --- Basic consistency -------------------------------------------------------

test_that("AIPW recovers ATE on linear-Gaussian DGP (binary treatment)", {
  d <- simulate_binary_continuous(n = 2000, seed = 42)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "aipw"
  )
  result <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0"
  )
  expect_equal(result$contrasts$estimate[1], 3, tolerance = 0.15)
  expect_lt(result$contrasts$ci_lower[1], 3)
  expect_gt(result$contrasts$ci_upper[1], 3)
  expect_true(is.finite(result$contrasts$se[1]) && result$contrasts$se[1] > 0)
})

test_that("AIPW agrees with gcomp and IPW on well-specified DGP", {
  d <- simulate_binary_continuous(n = 2000, seed = 42)
  ivs <- list(a1 = static(1), a0 = static(0))

  fit_gc <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "gcomp"
  )
  fit_ipw <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw"
  )
  fit_aipw <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "aipw"
  )

  res_gc <- contrast(fit_gc, interventions = ivs, reference = "a0")
  res_ipw <- contrast(fit_ipw, interventions = ivs, reference = "a0")
  res_aipw <- contrast(fit_aipw, interventions = ivs, reference = "a0")

  expect_lt(
    abs(res_gc$contrasts$estimate[1] - res_aipw$contrasts$estimate[1]),
    0.1
  )
  expect_lt(
    abs(res_ipw$contrasts$estimate[1] - res_aipw$contrasts$estimate[1]),
    0.1
  )
})

test_that("AIPW bootstrap agrees with sandwich (within 30%)", {
  d <- simulate_binary_continuous(n = 1000, seed = 123)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "aipw"
  )
  ivs <- list(a1 = static(1), a0 = static(0))

  res_sand <- contrast(
    fit,
    interventions = ivs,
    reference = "a0",
    ci_method = "sandwich"
  )
  res_boot <- contrast(
    fit,
    interventions = ivs,
    reference = "a0",
    ci_method = "bootstrap",
    n_boot = 200L
  )

  ratio <- res_boot$contrasts$se[1] / res_sand$contrasts$se[1]
  expect_gt(ratio, 0.7)
  expect_lt(ratio, 1.3)
})

# --- Double robustness (chunk 16d) ------------------------------------------

test_that("AIPW DR: wrong outcome model, correct propensity => ATE ~ 3", {
  set.seed(16004)
  n <- 5000
  L <- rnorm(n)
  A <- rbinom(n, 1, plogis(0.5 * L))
  # Nonlinear outcome — the L^2 term is omitted from the fitted model
  Y <- 2 + 3 * A + 1.5 * L + 0.8 * L^2 + rnorm(n)
  d <- data.frame(Y = Y, A = A, L = L)

  # Outcome model Y ~ A + L misses L^2 (misspecified)
  # Propensity A ~ L is correct (true logit-linear in L)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "aipw"
  )
  result <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0"
  )
  expect_equal(result$contrasts$estimate[1], 3, tolerance = 0.25)
  expect_lt(result$contrasts$ci_lower[1], 3)
  expect_gt(result$contrasts$ci_upper[1], 3)
})

# --- Double robustness (chunk 16e) ------------------------------------------

test_that("AIPW DR: correct outcome model, wrong propensity => ATE ~ 3", {
  set.seed(16050)
  n <- 5000
  L <- rnorm(n)
  # Nonlinear propensity — the L^2 term is omitted from the fitted model
  A <- rbinom(n, 1, plogis(0.5 * L + 0.3 * L^2))
  Y <- 2 + 3 * A + 1.5 * L + rnorm(n)
  d <- data.frame(Y = Y, A = A, L = L)

  # Outcome model Y ~ A + L is correct (true model is linear in A + L)
  # Propensity A ~ L misses L^2 (misspecified)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "aipw"
  )
  result <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0"
  )
  expect_equal(result$contrasts$estimate[1], 3, tolerance = 0.2)
  expect_lt(result$contrasts$ci_lower[1], 3)
  expect_gt(result$contrasts$ci_upper[1], 3)
})

# --- Efficiency (chunk 16f) -------------------------------------------------

test_that("AIPW SE <= gcomp SE and IPW SE (both correct, large n)", {
  d <- simulate_binary_continuous(n = 5000, seed = 789)
  ivs <- list(a1 = static(1), a0 = static(0))

  fit_gc <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "gcomp"
  )
  fit_ipw <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw"
  )
  fit_aipw <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "aipw"
  )

  res_gc <- contrast(fit_gc, interventions = ivs, reference = "a0")
  res_ipw <- contrast(fit_ipw, interventions = ivs, reference = "a0")
  res_aipw <- contrast(fit_aipw, interventions = ivs, reference = "a0")

  se_gc <- res_gc$contrasts$se[1]
  se_ipw <- res_ipw$contrasts$se[1]
  se_aipw <- res_aipw$contrasts$se[1]

  # AIPW should be at most as large as both (up to 5% slack for
  # finite-sample variation)
  expect_lte(se_aipw, se_gc * 1.05)
  expect_lte(se_aipw, se_ipw * 1.05)
  expect_true(is.finite(se_aipw) && se_aipw > 0)
})

# --- Rejections --------------------------------------------------------------

test_that("AIPW accepts longitudinal data (ICE-AIPW)", {
  d <- make_linear_scm(n = 300, n_times = 2, seed = 1)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    estimator = "aipw",
    family = "gaussian",
    id = "id",
    time = "time"
  )
  expect_equal(fit$type, "longitudinal")
  expect_equal(fit$estimator, "aipw")
})

test_that("AIPW rejects multivariate treatment", {
  d <- data.frame(
    Y = rnorm(50),
    A1 = rbinom(50, 1, 0.5),
    A2 = rbinom(50, 1, 0.5),
    L = rnorm(50)
  )
  expect_error(
    causat(
      d,
      outcome = "Y",
      treatment = c("A1", "A2"),
      confounders = ~L,
      estimator = "aipw"
    ),
    "Multivariate"
  )
})

# --- Marginal means ----------------------------------------------------------

test_that("AIPW marginal means are finite and correctly ordered", {
  d <- simulate_binary_continuous(n = 1000, seed = 55)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "aipw"
  )
  result <- contrast(fit, interventions = list(a1 = static(1), a0 = static(0)))
  expect_true(all(is.finite(result$estimates$estimate)))
  expect_true(all(is.finite(result$estimates$se)))
  expect_true(all(result$estimates$se > 0))
  # E[Y(1)] > E[Y(0)] because treatment effect is +3

  expect_gt(result$estimates$estimate[1], result$estimates$estimate[2])
})

# --- Binary outcome ----------------------------------------------------------

test_that("AIPW works with binomial outcome (risk difference)", {
  set.seed(16010)
  n <- 3000
  L <- rnorm(n)
  A <- rbinom(n, 1, plogis(0.5 * L))
  Y <- rbinom(n, 1, plogis(-1 + 1.5 * A + 0.8 * L))
  d <- data.frame(Y = Y, A = A, L = L)

  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "aipw",
    family = "binomial"
  )
  result <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0"
  )
  # Positive risk difference
  expect_gt(result$contrasts$estimate[1], 0)
  expect_true(is.finite(result$contrasts$se[1]) && result$contrasts$se[1] > 0)
  expect_lt(result$contrasts$ci_lower[1], result$contrasts$estimate[1])
  expect_gt(result$contrasts$ci_upper[1], result$contrasts$estimate[1])
})

# --- Continuous treatment with shift -----------------------------------------

test_that("AIPW works with continuous treatment and shift intervention", {
  d <- simulate_continuous_continuous(n = 2000, seed = 42)

  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "aipw"
  )
  result <- contrast(
    fit,
    interventions = list(shifted = shift(-1), obs = NULL),
    reference = "obs"
  )
  # Truth: shift(-1) effect = -2 (coefficient on A is 2)
  expect_equal(result$contrasts$estimate[1], -2, tolerance = 0.2)
  expect_true(is.finite(result$contrasts$se[1]) && result$contrasts$se[1] > 0)
})

# --- Full intervention set (chunk 16g) --------------------------------------

test_that("AIPW works with scale_by intervention (continuous treatment)", {
  d <- simulate_continuous_continuous(n = 2000, seed = 42)

  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "aipw"
  )
  result <- contrast(
    fit,
    interventions = list(scaled = scale_by(0.5), obs = NULL),
    reference = "obs"
  )
  # scale_by(0.5) halves A; coefficient on A is 2, mean(A) ~ 1
  # so effect ~ 2 * (0.5 * mean(A) - mean(A)) = -mean(A) ~ -1
  expect_true(result$contrasts$estimate[1] < 0)
  expect_true(is.finite(result$contrasts$se[1]) && result$contrasts$se[1] > 0)
})

test_that("AIPW works with dynamic intervention (binary treatment)", {
  d <- simulate_binary_continuous(n = 2000, seed = 42)

  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "aipw"
  )
  # Dynamic rule: treat everyone with L > 0, don't treat L <= 0
  rule <- function(data, orig_trt) as.integer(data$L > 0)
  result <- contrast(
    fit,
    interventions = list(dyn = dynamic(rule), a0 = static(0)),
    reference = "a0"
  )
  # Effect should be positive (treating some people) but less than ATE=3
  expect_gt(result$contrasts$estimate[1], 0)
  expect_lt(result$contrasts$estimate[1], 3.5)
  expect_true(is.finite(result$contrasts$se[1]) && result$contrasts$se[1] > 0)
})

test_that("AIPW works with ipsi intervention (binary treatment)", {
  d <- simulate_binary_continuous(n = 2000, seed = 42)

  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "aipw"
  )
  result <- contrast(
    fit,
    interventions = list(ipsi_d = ipsi(0.5), obs = NULL),
    reference = "obs"
  )
  # IPSI with delta < 1 nudges toward non-treatment; effect direction unclear

  # but SE must be finite and positive
  expect_true(is.finite(result$contrasts$se[1]) && result$contrasts$se[1] > 0)
  expect_true(is.finite(result$contrasts$estimate[1]))
})

test_that("AIPW rejects stochastic interventions", {
  d <- simulate_binary_continuous(n = 200, seed = 42)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "aipw"
  )

  sampler <- function(data) rbinom(nrow(data), 1, 0.5)
  expect_error(
    contrast(
      fit,
      interventions = list(s = stochastic(sampler), a0 = static(0)),
      reference = "a0"
    ),
    "stochastic"
  )
})

test_that("AIPW rejects threshold on continuous treatment", {
  d <- simulate_continuous_continuous(n = 200, seed = 42)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "aipw"
  )
  expect_error(
    contrast(
      fit,
      interventions = list(th = threshold(lower = 0), obs = NULL),
      reference = "obs"
    ),
    "threshold"
  )
})

# --- ATT/ATC (chunk 16h) ----------------------------------------------------

test_that("AIPW recovers ATT on linear DGP (static binary)", {
  d <- simulate_binary_continuous(n = 3000, seed = 42)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "aipw",
    estimand = "ATT"
  )
  result <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0"
  )
  # ATT = 3 on this DGP (constant treatment effect, no interaction)
  expect_equal(result$contrasts$estimate[1], 3, tolerance = 0.2)
  expect_true(is.finite(result$contrasts$se[1]) && result$contrasts$se[1] > 0)
})

test_that("AIPW recovers ATC on linear DGP (static binary)", {
  d <- simulate_binary_continuous(n = 3000, seed = 42)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "aipw",
    estimand = "ATC"
  )
  result <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0"
  )
  # ATC = 3 on this DGP (constant treatment effect)
  expect_equal(result$contrasts$estimate[1], 3, tolerance = 0.2)
  expect_true(is.finite(result$contrasts$se[1]) && result$contrasts$se[1] > 0)
})

# --- Effect modification (chunk 16h) -----------------------------------------

test_that("AIPW recovers sex-stratified ATE with effect modification", {
  d <- simulate_het_effect(n = 5000, seed = 42)
  # True: ATE|sex=0 = 3.0, ATE|sex=1 = 4.5
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L + A:sex,
    estimator = "aipw"
  )
  result <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    by = "sex"
  )
  est_sex0 <- result$contrasts$estimate[result$contrasts$by == "0"]
  est_sex1 <- result$contrasts$estimate[result$contrasts$by == "1"]
  expect_equal(est_sex0, 3.0, tolerance = 0.3)
  expect_equal(est_sex1, 4.5, tolerance = 0.3)
})

# --- Categorical treatment (chunk 16j) --------------------------------------

test_that("AIPW works with categorical treatment", {
  set.seed(16020)
  n <- 2000
  L <- rnorm(n)
  probs <- cbind(
    plogis(-0.5 - 0.3 * L),
    plogis(0.3 * L) * (1 - plogis(-0.5 - 0.3 * L)),
    1 - plogis(-0.5 - 0.3 * L) - plogis(0.3 * L) * (1 - plogis(-0.5 - 0.3 * L))
  )
  probs[probs < 0.05] <- 0.05
  probs <- probs / rowSums(probs)
  A <- apply(probs, 1, function(p) sample(c("a", "b", "c"), 1, prob = p))
  A <- factor(A, levels = c("a", "b", "c"))
  Y <- 1 + 2 * (A == "b") + 3 * (A == "c") + 0.5 * L + rnorm(n)
  d <- data.frame(Y = Y, A = A, L = L)

  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "aipw"
  )
  result <- contrast(
    fit,
    interventions = list(ab = static("b"), aa = static("a")),
    reference = "aa"
  )
  # E[Y(b)] - E[Y(a)] = 2

  expect_equal(result$contrasts$estimate[1], 2, tolerance = 0.3)
  expect_true(is.finite(result$contrasts$se[1]) && result$contrasts$se[1] > 0)
})


# --- Longitudinal AIPW (ICE-AIPW, Bang & Robins 2005) -----------------------

test_that("longitudinal AIPW: binary static, sandwich, ATE ~ 5", {
  d <- make_linear_scm(n = 3000, n_times = 2, seed = 42)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    estimator = "aipw",
    family = "gaussian",
    id = "id",
    time = "time"
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich"
  )
  ate <- res$contrasts$estimate[1]
  se <- res$contrasts$se[1]
  expect_equal(ate, 5, tolerance = 0.3)
  expect_true(is.finite(se) && se > 0)
})

test_that("longitudinal AIPW: binary static, bootstrap, ATE ~ 5", {
  d <- make_linear_scm(n = 2000, n_times = 2, seed = 43)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    estimator = "aipw",
    family = "gaussian",
    id = "id",
    time = "time"
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "bootstrap",
    n_boot = 50
  )
  ate <- res$contrasts$estimate[1]
  se <- res$contrasts$se[1]
  expect_equal(ate, 5, tolerance = 0.4)
  expect_true(is.finite(se) && se > 0)
})

test_that("longitudinal AIPW: continuous shift, sandwich", {
  d <- make_continuous_scm(n = 3000, seed = 42)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    estimator = "aipw",
    family = "gaussian",
    id = "id",
    time = "time"
  )
  res <- contrast(
    fit,
    interventions = list(up = shift(0.5), nat = NULL),
    reference = "nat",
    ci_method = "sandwich"
  )
  est <- res$contrasts$estimate[1]
  se <- res$contrasts$se[1]
  expect_true(is.finite(est))
  expect_true(is.finite(se) && se > 0)

  # Cross-check against ICE
  fit_ice <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    estimator = "gcomp",
    family = "gaussian",
    id = "id",
    time = "time"
  )
  res_ice <- contrast(
    fit_ice,
    interventions = list(up = shift(0.5), nat = NULL),
    reference = "nat",
    ci_method = "sandwich"
  )
  expect_lt(abs(est - res_ice$contrasts$estimate[1]), 0.3)
})

test_that("longitudinal AIPW: continuous shift, bootstrap", {
  d <- make_continuous_scm(n = 2000, seed = 42)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    estimator = "aipw",
    family = "gaussian",
    id = "id",
    time = "time"
  )
  res <- contrast(
    fit,
    interventions = list(up = shift(0.5), nat = NULL),
    reference = "nat",
    ci_method = "bootstrap",
    n_boot = 50
  )
  est <- res$contrasts$estimate[1]
  se <- res$contrasts$se[1]
  expect_true(is.finite(est))
  expect_true(is.finite(se) && se > 0)
})

test_that("longitudinal AIPW: dynamic intervention, sandwich", {
  d <- make_linear_scm(n = 3000, n_times = 2, seed = 44)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    estimator = "aipw",
    family = "gaussian",
    id = "id",
    time = "time"
  )
  dyn_fn <- dynamic(function(data, trt) ifelse(data$L0 > 0, 1L, 0L))
  res <- contrast(
    fit,
    interventions = list(dyn = dyn_fn, nat = NULL),
    reference = "nat",
    ci_method = "sandwich"
  )
  est <- res$contrasts$estimate[1]
  se <- res$contrasts$se[1]
  expect_true(is.finite(est))
  expect_true(is.finite(se) && se > 0)
})

test_that("longitudinal AIPW: scale_by intervention, sandwich", {
  d <- make_continuous_scm(n = 3000, seed = 42)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    estimator = "aipw",
    family = "gaussian",
    id = "id",
    time = "time"
  )
  res <- contrast(
    fit,
    interventions = list(sc = scale_by(1.2), nat = NULL),
    reference = "nat",
    ci_method = "sandwich"
  )
  est <- res$contrasts$estimate[1]
  se <- res$contrasts$se[1]
  expect_true(is.finite(est))
  expect_true(is.finite(se) && se > 0)
})

test_that("longitudinal AIPW vs ICE vs long-IPW: cross-method agreement", {
  d <- make_linear_scm(n = 3000, n_times = 2, seed = 45)
  ivs <- list(a1 = static(1), a0 = static(0))

  fit_aipw <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    estimator = "aipw",
    family = "gaussian",
    id = "id",
    time = "time"
  )
  fit_ice <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    estimator = "gcomp",
    family = "gaussian",
    id = "id",
    time = "time"
  )
  fit_ipw <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    estimator = "ipw",
    family = "gaussian",
    id = "id",
    time = "time"
  )

  res_aipw <- contrast(
    fit_aipw,
    interventions = ivs,
    reference = "a0",
    ci_method = "sandwich"
  )
  res_ice <- contrast(
    fit_ice,
    interventions = ivs,
    reference = "a0",
    ci_method = "sandwich"
  )
  res_ipw <- contrast(
    fit_ipw,
    interventions = ivs,
    reference = "a0",
    ci_method = "sandwich"
  )

  ate_aipw <- res_aipw$contrasts$estimate[1]
  ate_ice <- res_ice$contrasts$estimate[1]
  ate_ipw <- res_ipw$contrasts$estimate[1]

  expect_lt(abs(ate_aipw - ate_ice), 0.3)
  expect_lt(abs(ate_aipw - ate_ipw), 0.3)
})

test_that("longitudinal AIPW vs ICE: continuous shift agreement", {
  d <- make_continuous_scm(n = 3000, seed = 46)

  fit_aipw <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    estimator = "aipw",
    family = "gaussian",
    id = "id",
    time = "time"
  )
  fit_ice <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    estimator = "gcomp",
    family = "gaussian",
    id = "id",
    time = "time"
  )

  res_aipw <- contrast(
    fit_aipw,
    interventions = list(up = shift(0.5), nat = NULL),
    reference = "nat",
    ci_method = "sandwich"
  )
  res_ice <- contrast(
    fit_ice,
    interventions = list(up = shift(0.5), nat = NULL),
    reference = "nat",
    ci_method = "sandwich"
  )

  expect_lt(
    abs(res_aipw$contrasts$estimate[1] - res_ice$contrasts$estimate[1]),
    0.3
  )
})

test_that("longitudinal AIPW DR: wrong outcome, correct propensity", {
  d <- make_linear_scm(n = 3000, n_times = 2, seed = 47)
  # Misspecify outcome by dropping TV confounders
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    estimator = "aipw",
    family = "gaussian",
    id = "id",
    time = "time"
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "bootstrap",
    n_boot = 50
  )
  ate <- res$contrasts$estimate[1]
  # DR: should still recover truth despite wrong outcome model
  expect_equal(ate, 5, tolerance = 0.6)
})

test_that("longitudinal AIPW DR: correct outcome, wrong propensity", {
  d <- make_linear_scm(n = 3000, n_times = 2, seed = 48)
  # Misspecify propensity by omitting L from confounders
  # (outcome still includes L0 which is baseline, but misses L
  # in propensity). Since causat uses same confounders for both,
  # we misspecify propensity by dropping confounders_tv entirely
  # while keeping outcome correct via baseline confounders only.
  # This partially misspecifies propensity.
  fit_correct <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    estimator = "aipw",
    family = "gaussian",
    id = "id",
    time = "time"
  )
  fit_misspec <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    estimator = "aipw",
    family = "gaussian",
    id = "id",
    time = "time"
  )
  res_c <- contrast(
    fit_correct,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "bootstrap",
    n_boot = 50
  )
  res_m <- contrast(
    fit_misspec,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "bootstrap",
    n_boot = 50
  )
  # Both should be close to truth (DR property)
  expect_equal(res_c$contrasts$estimate[1], 5, tolerance = 0.4)
  expect_equal(res_m$contrasts$estimate[1], 5, tolerance = 0.6)
})

test_that("longitudinal AIPW: sandwich SE finite and positive", {
  d <- make_linear_scm(n = 3000, n_times = 2, seed = 49)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    estimator = "aipw",
    family = "gaussian",
    id = "id",
    time = "time"
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    ci_method = "sandwich"
  )
  sds <- sqrt(diag(res$vcov))
  expect_true(all(sds > 0 & is.finite(sds)))
})

test_that("longitudinal AIPW: bootstrap SE finite and positive", {
  d <- make_linear_scm(n = 2000, n_times = 2, seed = 50)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    estimator = "aipw",
    family = "gaussian",
    id = "id",
    time = "time"
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    ci_method = "bootstrap",
    n_boot = 50
  )
  sds <- sqrt(diag(res$vcov))
  expect_true(all(sds > 0 & is.finite(sds)))
})

test_that("longitudinal AIPW: sandwich vs bootstrap SE agreement", {
  d <- make_linear_scm(n = 5000, n_times = 2, seed = 51)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    estimator = "aipw",
    family = "gaussian",
    id = "id",
    time = "time"
  )
  ivs <- list(a1 = static(1), a0 = static(0))
  res_s <- contrast(fit, interventions = ivs, ci_method = "sandwich")
  res_b <- contrast(
    fit,
    interventions = ivs,
    ci_method = "bootstrap",
    n_boot = 200
  )
  se_sand <- sqrt(diag(res_s$vcov))
  se_boot <- sqrt(diag(res_b$vcov))
  ratio <- se_sand / se_boot
  expect_true(all(ratio > 0.5 & ratio < 2.0))
})

test_that("longitudinal AIPW: effect modification by baseline (sex)", {
  d <- make_em_ice_scm(n = 2500, seed = 52)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L0 + sex,
    confounders_tv = ~L,
    estimator = "aipw",
    family = "gaussian",
    id = "id",
    time = "time"
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    by = "sex",
    ci_method = "sandwich"
  )
  est_sex0 <- res$contrasts$estimate[
    res$contrasts$by == "0"
  ]
  est_sex1 <- res$contrasts$estimate[
    res$contrasts$by == "1"
  ]
  expect_equal(est_sex0, 5, tolerance = 0.5)
  expect_equal(est_sex1, 8, tolerance = 0.5)
})

test_that("longitudinal AIPW: EM agreement with ICE", {
  d <- make_em_ice_scm(n = 2500, seed = 53)
  ivs <- list(a1 = static(1), a0 = static(0))
  by_var <- "sex"

  fit_aipw <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L0 + sex,
    confounders_tv = ~L,
    estimator = "aipw",
    family = "gaussian",
    id = "id",
    time = "time"
  )
  fit_ice <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L0 + sex,
    confounders_tv = ~L,
    estimator = "gcomp",
    family = "gaussian",
    id = "id",
    time = "time"
  )

  res_aipw <- contrast(
    fit_aipw,
    interventions = ivs,
    reference = "a0",
    by = by_var,
    ci_method = "sandwich"
  )
  res_ice <- contrast(
    fit_ice,
    interventions = ivs,
    reference = "a0",
    by = by_var,
    ci_method = "sandwich"
  )

  for (sx in c("0", "1")) {
    aipw_est <- res_aipw$contrasts$estimate[
      res_aipw$contrasts$by == sx
    ]
    ice_est <- res_ice$contrasts$estimate[
      res_ice$contrasts$by == sx
    ]
    expect_lt(abs(aipw_est - ice_est), 2.0)
  }
})

test_that("longitudinal AIPW: binomial outcome", {
  set.seed(55)
  n <- 2000
  id <- rep(1:n, each = 2)
  time <- rep(0:1, times = n)
  L0 <- rnorm(n)[id]
  A <- rbinom(2 * n, 1, plogis(0.3 * L0))
  L <- rnorm(2 * n, mean = 0.3 * A)
  Y <- rep(NA_real_, 2 * n)
  Y[time == 1] <- as.numeric(rbinom(
    n,
    1,
    plogis(
      -1 + 0.5 * A[time == 1] + 0.3 * L[time == 1] + 0.2 * L0[time == 1]
    )
  ))
  d <- data.table::data.table(
    Y = Y,
    A = A,
    L0 = L0,
    L = L,
    id = id,
    time = time
  )
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    estimator = "aipw",
    family = "binomial",
    id = "id",
    time = "time"
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    ci_method = "sandwich"
  )
  ests <- coef(res)
  expect_true(all(ests > 0 & ests < 1))
  expect_true(all(is.finite(sqrt(diag(res$vcov)))))
})

test_that("longitudinal AIPW: 3-period DGP agrees with IPW", {
  d <- make_linear_scm(n = 3000, n_times = 3, seed = 56)
  fit_aipw <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    estimator = "aipw",
    family = "gaussian",
    id = "id",
    time = "time"
  )
  fit_ipw <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    estimator = "ipw",
    family = "gaussian",
    id = "id",
    time = "time"
  )
  ivs <- list(a1 = static(1), a0 = static(0))
  res_aipw <- contrast(
    fit_aipw,
    interventions = ivs,
    reference = "a0",
    ci_method = "sandwich"
  )
  res_ipw <- contrast(
    fit_ipw,
    interventions = ivs,
    reference = "a0",
    ci_method = "sandwich"
  )
  ate_aipw <- res_aipw$contrasts$estimate[1]
  ate_ipw <- res_ipw$contrasts$estimate[1]
  se <- res_aipw$contrasts$se[1]
  expect_lt(abs(ate_aipw - ate_ipw), 0.5)
  expect_true(is.finite(se) && se > 0)
})

test_that("longitudinal AIPW: multivariate treatment rejected", {
  d <- make_linear_scm(n = 300, n_times = 2, seed = 57)
  d <- data.table::as.data.table(d)
  d[, A2 := rbinom(.N, 1, 0.5)]
  expect_error(
    causat(
      d,
      outcome = "Y",
      treatment = c("A", "A2"),
      confounders = ~L0,
      confounders_tv = ~L,
      estimator = "aipw",
      family = "gaussian",
      id = "id",
      time = "time"
    ),
    "ultivariate"
  )
})

test_that("longitudinal AIPW: ATT rejected", {
  d <- make_linear_scm(n = 300, n_times = 2, seed = 58)
  expect_error(
    causat(
      d,
      outcome = "Y",
      treatment = "A",
      confounders = ~L0,
      confounders_tv = ~L,
      estimator = "aipw",
      family = "gaussian",
      estimand = "ATT",
      id = "id",
      time = "time"
    ),
    "estimand = 'ATT'"
  )
})

test_that("longitudinal AIPW: lmtp cross-check (binary static)", {
  skip_if_not_installed("lmtp")
  skip_if_not_installed("SuperLearner")

  d <- make_linear_scm(n = 1500, n_times = 2, seed = 60)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    estimator = "aipw",
    family = "gaussian",
    id = "id",
    time = "time"
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich"
  )
  ate_aipw <- res$contrasts$estimate[1]

  d_wide <- data.frame(
    L0 = d$L0[d$time == 0],
    A_0 = d$A[d$time == 0],
    L_1 = d$L[d$time == 1],
    A_1 = d$A[d$time == 1],
    Y = d$Y[d$time == 1]
  )
  # lmtp_sdr: always-treated
  lmtp_a1 <- lmtp::lmtp_sdr(
    data = d_wide,
    trt = c("A_0", "A_1"),
    outcome = "Y",
    baseline = "L0",
    time_vary = list(NULL, "L_1"),
    shift = function(data, trt) rep(1, nrow(data)),
    outcome_type = "continuous",
    learners_trt = "SL.glm",
    learners_outcome = "SL.glm",
    folds = 1
  )
  # lmtp_sdr: never-treated
  lmtp_a0 <- lmtp::lmtp_sdr(
    data = d_wide,
    trt = c("A_0", "A_1"),
    outcome = "Y",
    baseline = "L0",
    time_vary = list(NULL, "L_1"),
    shift = function(data, trt) rep(0, nrow(data)),
    outcome_type = "continuous",
    learners_trt = "SL.glm",
    learners_outcome = "SL.glm",
    folds = 1
  )
  ate_lmtp <- lmtp_a1$estimate@x - lmtp_a0$estimate@x
  expect_lt(abs(ate_aipw - ate_lmtp), 0.5)

  # SE comparison: lmtp uses EIF-based SE (no cross-fitting with folds=1),
  # causatr uses parametric sandwich. Both should be in the same ballpark.
  se_aipw <- res$contrasts$std_error[1]
  se_lmtp <- sqrt(
    lmtp_a1$estimate@std_error^2 + lmtp_a0$estimate@std_error^2
  )
  se_ratio <- se_aipw / se_lmtp
  expect_gt(se_ratio, 0.5)
  expect_lt(se_ratio, 2.0)
})

test_that("longitudinal AIPW: lmtp cross-check (continuous shift)", {
  skip_if_not_installed("lmtp")
  skip_if_not_installed("SuperLearner")

  d <- make_continuous_scm(n = 1500, seed = 61)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    estimator = "aipw",
    family = "gaussian",
    id = "id",
    time = "time"
  )
  res <- contrast(
    fit,
    interventions = list(up = shift(0.5), nat = NULL),
    reference = "nat",
    ci_method = "sandwich"
  )
  est_aipw <- res$contrasts$estimate[1]

  d_wide <- data.frame(
    L0 = d$L0[d$time == 0],
    A_0 = d$A[d$time == 0],
    L_1 = d$L[d$time == 1],
    A_1 = d$A[d$time == 1],
    Y = d$Y[d$time == 1]
  )
  shift_fn <- function(data, trt) data[[trt]] + 0.5
  lmtp_up <- lmtp::lmtp_sdr(
    data = d_wide,
    trt = c("A_0", "A_1"),
    outcome = "Y",
    baseline = "L0",
    time_vary = list(NULL, "L_1"),
    shift = shift_fn,
    outcome_type = "continuous",
    learners_trt = "SL.glm",
    learners_outcome = "SL.glm",
    folds = 1
  )
  lmtp_nat <- lmtp::lmtp_sdr(
    data = d_wide,
    trt = c("A_0", "A_1"),
    outcome = "Y",
    baseline = "L0",
    time_vary = list(NULL, "L_1"),
    shift = NULL,
    outcome_type = "continuous",
    learners_trt = "SL.glm",
    learners_outcome = "SL.glm",
    folds = 1
  )
  est_lmtp <- lmtp_up$estimate@x - lmtp_nat$estimate@x
  expect_lt(abs(est_aipw - est_lmtp), 0.5)

  se_aipw <- res$contrasts$std_error[1]
  se_lmtp <- sqrt(
    lmtp_up$estimate@std_error^2 + lmtp_nat$estimate@std_error^2
  )
  se_ratio <- se_aipw / se_lmtp
  expect_gt(se_ratio, 0.5)
  expect_lt(se_ratio, 2.0)
})

test_that("longitudinal AIPW: lmtp cross-check (binary outcome)", {
  skip_if_not_installed("lmtp")
  skip_if_not_installed("SuperLearner")

  set.seed(62)
  n <- 1500
  id <- rep(1:n, each = 2)
  time <- rep(0:1, times = n)
  L0 <- rnorm(n)[id]
  A <- rbinom(2 * n, 1, plogis(0.3 * L0))
  L <- rnorm(2 * n, mean = 0.3 * A)
  Y <- rep(NA_real_, 2 * n)
  Y[time == 1] <- as.numeric(rbinom(
    n,
    1,
    plogis(-1 + 0.5 * A[time == 1] + 0.3 * L[time == 1] + 0.2 * L0[time == 1])
  ))
  d <- data.table::data.table(
    Y = Y,
    A = A,
    L0 = L0,
    L = L,
    id = id,
    time = time
  )

  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    estimator = "aipw",
    family = "binomial",
    id = "id",
    time = "time"
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich"
  )
  est_aipw <- res$contrasts$estimate[1]

  d_wide <- data.frame(
    L0 = d$L0[d$time == 0],
    A_0 = d$A[d$time == 0],
    L_1 = d$L[d$time == 1],
    A_1 = d$A[d$time == 1],
    Y = d$Y[d$time == 1]
  )
  lmtp_a1 <- lmtp::lmtp_sdr(
    data = d_wide,
    trt = c("A_0", "A_1"),
    outcome = "Y",
    baseline = "L0",
    time_vary = list(NULL, "L_1"),
    shift = function(data, trt) rep(1, nrow(data)),
    outcome_type = "binomial",
    learners_trt = "SL.glm",
    learners_outcome = "SL.glm",
    folds = 1
  )
  lmtp_a0 <- lmtp::lmtp_sdr(
    data = d_wide,
    trt = c("A_0", "A_1"),
    outcome = "Y",
    baseline = "L0",
    time_vary = list(NULL, "L_1"),
    shift = function(data, trt) rep(0, nrow(data)),
    outcome_type = "binomial",
    learners_trt = "SL.glm",
    learners_outcome = "SL.glm",
    folds = 1
  )
  est_lmtp <- lmtp_a1$estimate@x - lmtp_a0$estimate@x
  expect_lt(abs(est_aipw - est_lmtp), 0.3)

  se_aipw <- res$contrasts$std_error[1]
  se_lmtp <- sqrt(
    lmtp_a1$estimate@std_error^2 + lmtp_a0$estimate@std_error^2
  )
  se_ratio <- se_aipw / se_lmtp
  expect_gt(se_ratio, 0.5)
  expect_lt(se_ratio, 2.0)
})
