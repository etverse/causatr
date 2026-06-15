# --- Basic consistency -------------------------------------------------------

test_that("AIPW recovers ATE on linear-Gaussian DGP (binary treatment)", {
  d <- simulate_binary_continuous(n = 2000, seed = 42)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "aipw",
    propensity_model_fn = stats::glm
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
    estimator = "aipw",
    propensity_model_fn = stats::glm
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
  skip_if_fast()
  d <- simulate_binary_continuous(n = 1000, seed = 123)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "aipw",
    propensity_model_fn = stats::glm
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
    n_boot = 100L
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
    estimator = "aipw",
    propensity_model_fn = stats::glm
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
    estimator = "aipw",
    propensity_model_fn = stats::glm
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
  d <- simulate_binary_continuous(n = 3000, seed = 789)
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
    estimator = "aipw",
    propensity_model_fn = stats::glm
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
    propensity_model_fn = stats::glm,
    family = "gaussian",
    id = "id",
    time = "time"
  )
  expect_equal(fit$type, "longitudinal")
  expect_equal(fit$estimator, "aipw")
})

test_that("AIPW rejects stabilize for univariate treatment", {
  d <- simulate_binary_continuous(n = 200, seed = 80)
  expect_error(
    causat(
      d,
      outcome = "Y",
      treatment = "A",
      confounders = ~L,
      estimator = "aipw",
      propensity_model_fn = stats::glm,
      stabilize = "marginal"
    ),
    class = "causatr_stabilize_univariate"
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
    estimator = "aipw",
    propensity_model_fn = stats::glm
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
    propensity_model_fn = stats::glm,
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
    estimator = "aipw",
    propensity_model_fn = stats::glm
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
    estimator = "aipw",
    propensity_model_fn = stats::glm
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
    estimator = "aipw",
    propensity_model_fn = stats::glm
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
    estimator = "aipw",
    propensity_model_fn = stats::glm
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
    estimator = "aipw",
    propensity_model_fn = stats::glm
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
    estimator = "aipw",
    propensity_model_fn = stats::glm
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
    propensity_model_fn = stats::glm,
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
    propensity_model_fn = stats::glm,
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
  d <- simulate_het_effect(n = 3000, seed = 42)
  # True: ATE|sex=0 = 3.0, ATE|sex=1 = 4.5
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L + A:sex,
    estimator = "aipw",
    propensity_model_fn = stats::glm
  )
  result <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    by = "sex"
  )
  est_sex0 <- result$contrasts$estimate[result$contrasts$by == "0"]
  est_sex1 <- result$contrasts$estimate[result$contrasts$by == "1"]
  expect_equal(est_sex0, 3.0, tolerance = 0.05)
  expect_equal(est_sex1, 4.5, tolerance = 0.05)
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
    estimator = "aipw",
    propensity_model_fn = stats::glm
  )
  result <- contrast(
    fit,
    interventions = list(ab = static("b"), aa = static("a")),
    reference = "aa"
  )
  # E[Y(b)] - E[Y(a)] = 2

  expect_equal(result$contrasts$estimate[1], 2, tolerance = 0.05)
  expect_true(is.finite(result$contrasts$se[1]) && result$contrasts$se[1] > 0)
})
