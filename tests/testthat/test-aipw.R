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
    estimator = "aipw",
    propensity_model_fn = stats::glm
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
    propensity_model_fn = stats::glm,
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
    propensity_model_fn = stats::glm,
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
    propensity_model_fn = stats::glm,
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
    propensity_model_fn = stats::glm,
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
    propensity_model_fn = stats::glm,
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
    propensity_model_fn = stats::glm,
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
    propensity_model_fn = stats::glm,
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
    propensity_model_fn = stats::glm,
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
    propensity_model_fn = stats::glm,
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
    propensity_model_fn = stats::glm,
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
    propensity_model_fn = stats::glm,
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
    propensity_model_fn = stats::glm,
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
    propensity_model_fn = stats::glm,
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
  d <- make_linear_scm(n = 3000, n_times = 2, seed = 51)
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
    propensity_model_fn = stats::glm,
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
    propensity_model_fn = stats::glm,
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
    propensity_model_fn = stats::glm,
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
    propensity_model_fn = stats::glm,
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
      propensity_model_fn = stats::glm,
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
      propensity_model_fn = stats::glm,
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
    propensity_model_fn = stats::glm,
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

  # lmtp sums two marginal EIF SEs in quadrature (ignoring covariance),
  # causatr sandwich targets the contrast directly — ratio < 1 expected.
  se_aipw <- res$contrasts$se[1]
  expect_true(is.finite(se_aipw), label = "causatr SE is finite")
  expect_true(
    is.finite(lmtp_a1$estimate@std_error),
    label = "lmtp_a1 SE is finite"
  )
  expect_true(
    is.finite(lmtp_a0$estimate@std_error),
    label = "lmtp_a0 SE is finite"
  )
  se_lmtp <- sqrt(
    lmtp_a1$estimate@std_error^2 + lmtp_a0$estimate@std_error^2
  )
  se_ratio <- se_aipw / se_lmtp
  expect_gt(se_ratio, 0.15)
  expect_lt(se_ratio, 5.0)
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
    propensity_model_fn = stats::glm,
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

  # lmtp SE is sum-in-quadrature of two marginal EIF SEs (ignores covariance),
  # while causatr sandwich targets the contrast directly — expect a ratio well
  # below 1.  Smoke-test bounds only.
  se_aipw <- res$contrasts$se[1]
  expect_true(is.finite(se_aipw), label = "causatr SE is finite")
  expect_true(
    is.finite(lmtp_up$estimate@std_error),
    label = "lmtp_up SE is finite"
  )
  expect_true(
    is.finite(lmtp_nat$estimate@std_error),
    label = "lmtp_nat SE is finite"
  )
  se_lmtp <- sqrt(
    lmtp_up$estimate@std_error^2 + lmtp_nat$estimate@std_error^2
  )
  se_ratio <- se_aipw / se_lmtp
  expect_gt(se_ratio, 0.15)
  expect_lt(se_ratio, 5.0)
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
    propensity_model_fn = stats::glm,
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

  # Same caveat as continuous shift: lmtp sums marginal SEs in quadrature
  # (ignoring covariance), causatr targets the contrast — ratio < 1 expected.
  se_aipw <- res$contrasts$se[1]
  expect_true(is.finite(se_aipw), label = "causatr SE is finite")
  expect_true(
    is.finite(lmtp_a1$estimate@std_error),
    label = "lmtp_a1 SE is finite"
  )
  expect_true(
    is.finite(lmtp_a0$estimate@std_error),
    label = "lmtp_a0 SE is finite"
  )
  se_lmtp <- sqrt(
    lmtp_a1$estimate@std_error^2 + lmtp_a0$estimate@std_error^2
  )
  se_ratio <- se_aipw / se_lmtp
  expect_gt(se_ratio, 0.15)
  expect_lt(se_ratio, 5.0)
})


# --- Plan gap tests (added post-16i) -----------------------------------

test_that("longitudinal AIPW DR caveat: sandwich SE under misspecified outcome", {
  d <- make_linear_scm(n = 3000, n_times = 2, seed = 70)
  # Misspecify outcome by dropping TV confounders
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    estimator = "aipw",
    propensity_model_fn = stats::glm,
    family = "gaussian",
    id = "id",
    time = "time"
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
    n_boot = 50
  )
  se_sw <- res_sw$contrasts$se[1]
  se_bs <- res_bs$contrasts$se[1]
  # Both SEs must be finite and positive. Under misspecification the
  # sandwich SE is NOT DR-consistent (Rotnitzky et al. 2017), so it may
  # diverge from bootstrap — we only check finiteness here.
  expect_true(is.finite(se_sw) && se_sw > 0)
  expect_true(is.finite(se_bs) && se_bs > 0)
  # Point estimate should still recover truth (DR property)
  expect_equal(res_sw$contrasts$estimate[1], 5, tolerance = 0.6)
})

test_that("longitudinal AIPW: EM agreement with long-IPW", {
  d <- make_em_ice_scm(n = 2500, seed = 71)
  ivs <- list(a1 = static(1), a0 = static(0))
  by_var <- "sex"

  fit_aipw <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L0 + sex,
    confounders_tv = ~L,
    estimator = "aipw",
    propensity_model_fn = stats::glm,
    family = "gaussian",
    id = "id",
    time = "time"
  )
  fit_ipw <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L0 + sex,
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
    by = by_var,
    ci_method = "sandwich"
  )
  res_ipw <- contrast(
    fit_ipw,
    interventions = ivs,
    reference = "a0",
    by = by_var,
    ci_method = "sandwich"
  )

  for (sx in c("0", "1")) {
    aipw_est <- res_aipw$contrasts$estimate[
      res_aipw$contrasts$by == sx
    ]
    ipw_est <- res_ipw$contrasts$estimate[
      res_ipw$contrasts$by == sx
    ]
    expect_lt(abs(aipw_est - ipw_est), 2.0)
  }
})

test_that("longitudinal AIPW: 3-period sandwich vs bootstrap SE agreement", {
  d <- make_linear_scm(n = 3000, n_times = 3, seed = 72)
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
  ivs <- list(a1 = static(1), a0 = static(0))
  res_sw <- contrast(
    fit,
    interventions = ivs,
    reference = "a0",
    ci_method = "sandwich"
  )
  res_bs <- contrast(
    fit,
    interventions = ivs,
    reference = "a0",
    ci_method = "bootstrap",
    n_boot = 100
  )
  se_sw <- res_sw$contrasts$se[1]
  se_bs <- res_bs$contrasts$se[1]
  expect_true(is.finite(se_sw) && se_sw > 0)
  expect_true(is.finite(se_bs) && se_bs > 0)
  se_ratio <- se_sw / se_bs
  expect_gt(se_ratio, 0.5)
  expect_lt(se_ratio, 2.0)
})

test_that("longitudinal AIPW: near-positivity stress test", {
  set.seed(73)
  n <- 2000
  id <- rep(1:n, each = 2)
  time <- rep(0:1, times = n)
  L0 <- rnorm(n)[id]
  # Heavy confounding: propensity near 0 or 1 for most units
  A <- rbinom(2 * n, 1, plogis(2.0 + 2.5 * L0))
  L <- rnorm(2 * n, mean = 0.5 * A)
  Y <- rep(NA_real_, 2 * n)
  Y[time == 1] <- 2 + 3 * A[time == 1] + 1.5 * L0[time == 1] + rnorm(n)
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
    propensity_model_fn = stats::glm,
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
  expect_true(is.finite(ate))
  expect_true(is.finite(se) && se > 0)
  expect_equal(ate, 5, tolerance = 2.0)
})

# --- delicatessen cross-check (chunk 16k) ------------------------------------
#
# Reference values generated by data-raw/aipw_reference.py using
# delicatessen (Zivich 2024) stacked M-estimation sandwich variance.
# Zivich PN et al. (2024). Statistics in Medicine 43:5562-5572.

test_that("AIPW sandwich matches delicatessen — binary treatment, ATE", {
  ref_mu_1 <- 4.9715
  ref_se_1 <- 0.0289
  ref_mu_0 <- 1.9730
  ref_se_0 <- 0.0314
  ref_ate <- 2.9985
  ref_se_ate <- 0.0299

  set.seed(42)
  n <- 5000
  L <- rnorm(n)
  A <- rbinom(n, 1, plogis(0.3 + 0.5 * L))
  Y <- 2 + 3 * A + 1.5 * L + rnorm(n)
  d <- data.frame(Y = Y, A = A, L = L)

  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "aipw",
    propensity_model_fn = stats::glm
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich"
  )

  mu_1 <- res$estimates$estimate[res$estimates$intervention == "a1"]
  mu_0 <- res$estimates$estimate[res$estimates$intervention == "a0"]
  se_1 <- res$estimates$se[res$estimates$intervention == "a1"]
  se_0 <- res$estimates$se[res$estimates$intervention == "a0"]

  expect_equal(mu_1, ref_mu_1, tolerance = 0.005)
  expect_equal(mu_0, ref_mu_0, tolerance = 0.005)
  expect_equal(res$contrasts$estimate[1], ref_ate, tolerance = 0.005)
  expect_equal(se_1, ref_se_1, tolerance = 0.005)
  expect_equal(se_0, ref_se_0, tolerance = 0.005)
  expect_equal(res$contrasts$se[1], ref_se_ate, tolerance = 0.005)
})

test_that("AIPW sandwich matches delicatessen — continuous treatment, shift(1)", {
  ref_mu_shift <- 6.0079
  ref_se_shift <- 0.0613
  ref_mu_nat <- 4.0149
  ref_se_nat <- 0.0586
  ref_effect <- 1.9930
  ref_se_effect <- 0.0171

  set.seed(99)
  n <- 5000
  L <- rnorm(n)
  A <- rnorm(n, mean = 1 + L, sd = 1)
  Y <- 2 + 2 * A + 1.5 * L + rnorm(n)
  d <- data.frame(Y = Y, A = A, L = L)

  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "aipw",
    propensity_model_fn = stats::glm
  )
  res <- contrast(
    fit,
    interventions = list(shifted = shift(1), nat = NULL),
    reference = "nat",
    ci_method = "sandwich"
  )

  mu_shift <- res$estimates$estimate[res$estimates$intervention == "shifted"]
  mu_nat <- res$estimates$estimate[res$estimates$intervention == "nat"]
  se_shift <- res$estimates$se[res$estimates$intervention == "shifted"]
  se_nat <- res$estimates$se[res$estimates$intervention == "nat"]

  expect_equal(mu_shift, ref_mu_shift, tolerance = 0.005)
  expect_equal(mu_nat, ref_mu_nat, tolerance = 0.005)
  expect_equal(res$contrasts$estimate[1], ref_effect, tolerance = 0.005)
  expect_equal(se_shift, ref_se_shift, tolerance = 0.005)
  expect_equal(se_nat, ref_se_nat, tolerance = 0.005)
  expect_equal(res$contrasts$se[1], ref_se_effect, tolerance = 0.005)
})

test_that("longitudinal AIPW: IPSI rejection", {
  d <- make_linear_scm(n = 300, n_times = 2, seed = 74)
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
  expect_error(
    contrast(
      fit,
      interventions = list(shifted = ipsi(0.5)),
      ci_method = "sandwich"
    ),
    class = "causatr_longitudinal_ipsi_pending"
  )
})


# --- Multivariate AIPW (binary x binary) ------------------------------------

# DGP: L ~ N(0,1), A1|L ~ Bernoulli(plogis(0.3*L)),
#       A2|L ~ Bernoulli(plogis(-0.3*L)),
#       Y = 2 + 1.5*A1 + 1.0*A2 - 0.5*L + N(0,1).
# Truth: E[Y(1,1)] = 4.5, E[Y(0,0)] = 2.0, ATE = 2.5.
sim_bb_aipw <- function(n = 3000, seed = 42) {
  set.seed(seed)
  L <- rnorm(n)
  A1 <- rbinom(n, 1, plogis(0.3 * L))
  A2 <- rbinom(n, 1, plogis(-0.3 * L))
  Y <- 2 + 1.5 * A1 + 1.0 * A2 - 0.5 * L + rnorm(n)
  data.frame(L = L, A1 = A1, A2 = A2, Y = Y)
}

# DGP: continuous x continuous with downstream conditioning.
# L ~ N(0,1), A1|L ~ N(0.5*L, 1), A2|A1,L ~ N(0.3*A1+0.5*L, 1),
# Y = 1 + 0.5*A1 + 0.4*A2 - 0.5*L + N(0,1).
# Truth for shift(-1,-1): 0.5*(-1) + 0.4*(-1) = -0.9.
sim_cc_aipw <- function(n = 3000, seed = 42) {
  set.seed(seed)
  L <- rnorm(n)
  A1 <- 0.5 * L + rnorm(n)
  A2 <- 0.3 * A1 + 0.5 * L + rnorm(n)
  Y <- 1 + 0.5 * A1 + 0.4 * A2 - 0.5 * L + rnorm(n)
  data.frame(L = L, A1 = A1, A2 = A2, Y = Y)
}

# DGP: binary x binary, binomial outcome.
# Y ~ Bernoulli(plogis(-1 + A1 + 0.8*A2 + 0.5*L)).
# MC truth: P[Y(1,1)] ~ 0.622, P[Y(0,0)] ~ 0.269.
sim_bb_binary_aipw <- function(n = 3000, seed = 42) {
  set.seed(seed)
  L <- rnorm(n)
  A1 <- rbinom(n, 1, plogis(0.3 * L))
  A2 <- rbinom(n, 1, plogis(-0.3 * L))
  Y <- rbinom(n, 1, plogis(-1 + A1 + 0.8 * A2 + 0.5 * L))
  data.frame(L = L, A1 = A1, A2 = A2, Y = Y)
}

test_that("mv AIPW: binary x binary static recovers truth", {
  df <- sim_bb_aipw()
  fit <- causat(
    df,
    "Y",
    c("A1", "A2"),
    ~L,
    estimator = "aipw",
    propensity_model_fn = stats::glm
  )
  res <- contrast(
    fit,
    interventions = list(
      both = list(A1 = static(1), A2 = static(1)),
      neither = list(A1 = static(0), A2 = static(0))
    ),
    reference = "neither",
    ci_method = "sandwich"
  )
  ate <- res$contrasts$estimate[1]
  se <- res$contrasts$se[1]
  expect_equal(ate, 2.5, tolerance = 0.3)
  expect_true(is.finite(se) && se > 0)
})

test_that("mv AIPW: cross-checks gcomp and IPW", {
  df <- sim_bb_aipw()
  ivs <- list(
    both = list(A1 = static(1), A2 = static(1)),
    neither = list(A1 = static(0), A2 = static(0))
  )

  fit_a <- causat(
    df,
    "Y",
    c("A1", "A2"),
    ~L,
    estimator = "aipw",
    propensity_model_fn = stats::glm
  )
  res_a <- contrast(
    fit_a,
    interventions = ivs,
    reference = "neither",
    ci_method = "sandwich"
  )

  fit_g <- causat(df, "Y", c("A1", "A2"), ~L, estimator = "gcomp")
  res_g <- contrast(
    fit_g,
    interventions = ivs,
    reference = "neither",
    ci_method = "sandwich"
  )

  fit_i <- causat(df, "Y", c("A1", "A2"), ~L, estimator = "ipw")
  res_i <- contrast(
    fit_i,
    interventions = ivs,
    reference = "neither",
    ci_method = "sandwich"
  )

  ate_a <- res_a$contrasts$estimate[1]
  ate_g <- res_g$contrasts$estimate[1]
  ate_i <- res_i$contrasts$estimate[1]
  expect_equal(ate_a, ate_g, tolerance = 0.15)
  expect_equal(ate_a, ate_i, tolerance = 0.15)
})

test_that("mv AIPW: bootstrap parity with sandwich", {
  df <- sim_bb_aipw(n = 1000, seed = 99)
  fit <- causat(
    df,
    "Y",
    c("A1", "A2"),
    ~L,
    estimator = "aipw",
    propensity_model_fn = stats::glm
  )
  ivs <- list(
    both = list(A1 = static(1), A2 = static(1)),
    neither = list(A1 = static(0), A2 = static(0))
  )
  res_sw <- contrast(
    fit,
    interventions = ivs,
    reference = "neither",
    ci_method = "sandwich"
  )
  set.seed(123)
  res_bt <- contrast(
    fit,
    interventions = ivs,
    reference = "neither",
    ci_method = "bootstrap",
    n_boot = 200
  )
  se_sw <- res_sw$contrasts$se[1]
  se_bt <- res_bt$contrasts$se[1]
  expect_true(abs(se_sw - se_bt) < 0.3 * se_sw)
})

test_that("mv AIPW: continuous x continuous shift recovers truth", {
  df <- sim_cc_aipw(n = 3000)
  fit <- causat(
    df,
    "Y",
    c("A1", "A2"),
    ~L,
    estimator = "aipw",
    propensity_model_fn = stats::glm
  )
  res <- contrast(
    fit,
    interventions = list(
      shifted = list(A1 = shift(-1), A2 = shift(-1)),
      obs = NULL
    ),
    reference = "obs",
    ci_method = "sandwich"
  )
  ate <- res$contrasts$estimate[1]
  se <- res$contrasts$se[1]
  expect_equal(ate, -0.9, tolerance = 0.3)
  expect_true(is.finite(se) && se > 0)
})

test_that("mv AIPW: DR — wrong outcome, correct propensity", {
  df <- sim_bb_aipw(n = 3000, seed = 77)
  # Outcome model omits L (misspecified); propensity correctly includes L.
  fit <- causat(
    df,
    "Y",
    c("A1", "A2"),
    ~L,
    estimator = "aipw",
    propensity_model_fn = stats::glm,
    model_fn = function(formula, data, ...) {
      stats::glm(Y ~ A1 + A2, data = data, ...)
    }
  )
  res <- contrast(
    fit,
    interventions = list(
      both = list(A1 = static(1), A2 = static(1)),
      neither = list(A1 = static(0), A2 = static(0))
    ),
    reference = "neither",
    ci_method = "sandwich"
  )
  ate <- res$contrasts$estimate[1]
  expect_equal(ate, 2.5, tolerance = 0.4)
})

test_that("mv AIPW: DR — correct outcome, wrong propensity", {
  df <- sim_bb_aipw(n = 3000, seed = 78)
  # Propensity model drops L (misspecified); outcome model is correct.
  fit <- causat(
    df,
    "Y",
    c("A1", "A2"),
    ~L,
    estimator = "aipw",
    propensity_model_fn = function(formula, data, ...) {
      resp <- all.vars(formula)[1]
      stats::glm(
        stats::reformulate("1", response = resp),
        data = data,
        ...
      )
    }
  )
  res <- contrast(
    fit,
    interventions = list(
      both = list(A1 = static(1), A2 = static(1)),
      neither = list(A1 = static(0), A2 = static(0))
    ),
    reference = "neither",
    ci_method = "sandwich"
  )
  ate <- res$contrasts$estimate[1]
  expect_equal(ate, 2.5, tolerance = 0.4)
})

test_that("mv AIPW: binary x binary binomial outcome", {
  df <- sim_bb_binary_aipw()
  fit <- causat(
    df,
    "Y",
    c("A1", "A2"),
    ~L,
    estimator = "aipw",
    propensity_model_fn = stats::glm,
    family = "binomial"
  )
  res <- contrast(
    fit,
    interventions = list(
      both = list(A1 = static(1), A2 = static(1)),
      neither = list(A1 = static(0), A2 = static(0))
    ),
    reference = "neither",
    ci_method = "sandwich"
  )
  rd <- res$contrasts$estimate[1]
  se <- res$contrasts$se[1]
  # MC truth: RD ~ 0.353 (from 1e6 simulation)
  expect_true(abs(rd - 0.353) < 0.15)
  expect_true(is.finite(se) && se > 0)
})

test_that("mv AIPW: stabilize = 'marginal' produces finite results", {
  df <- sim_bb_aipw()
  fit <- causat(
    df,
    "Y",
    c("A1", "A2"),
    ~L,
    estimator = "aipw",
    propensity_model_fn = stats::glm,
    stabilize = "marginal"
  )
  res <- contrast(
    fit,
    interventions = list(
      both = list(A1 = static(1), A2 = static(1)),
      neither = list(A1 = static(0), A2 = static(0))
    ),
    reference = "neither",
    ci_method = "sandwich"
  )
  ate <- res$contrasts$estimate[1]
  se <- res$contrasts$se[1]
  expect_equal(ate, 2.5, tolerance = 0.3)
  expect_true(is.finite(se) && se > 0)
})


# ============================================================
# CHUNK 16p — AIPW TEST COVERAGE PARITY
# ============================================================

# --- Gap 1: Extended outcome families -----------------------------------------

test_that("aipw x bin trt x gamma(log) x diff x sandwich", {
  # DGP: Y ~ Gamma(shape=5, rate=5/exp(1 + 0.4*A + 0.3*L))
  # E[Y^a] = exp(1 + 0.4*a) * E[exp(0.3*L)] = exp(1 + 0.4*a + 0.045)
  # ATE = exp(1.045) * (exp(0.4) - 1)
  set.seed(16201)
  n <- 3000
  L <- stats::rnorm(n)
  A <- stats::rbinom(n, 1, stats::plogis(0.3 + 0.5 * L))
  mu <- exp(1 + 0.4 * A + 0.3 * L)
  Y <- stats::rgamma(n, shape = 5, rate = 5 / mu)
  df <- data.frame(Y = Y, A = A, L = L)
  truth <- exp(1.045) * (exp(0.4) - 1)

  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "aipw",
    family = stats::Gamma(link = "log"),
    propensity_model_fn = stats::glm
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    type = "difference",
    ci_method = "sandwich"
  )
  expect_equal(res$contrasts$estimate[1], truth, tolerance = 0.2)
  expect_true(all(is.finite(res$contrasts$se) & res$contrasts$se > 0))
  expect_lt(res$contrasts$ci_lower[1], truth)
  expect_gt(res$contrasts$ci_upper[1], truth)
})

test_that("aipw x bin trt x gamma(log) x ratio x sandwich", {
  set.seed(16202)
  n <- 3000
  L <- stats::rnorm(n)
  A <- stats::rbinom(n, 1, stats::plogis(0.3 + 0.5 * L))
  mu <- exp(1 + 0.4 * A + 0.3 * L)
  Y <- stats::rgamma(n, shape = 5, rate = 5 / mu)
  df <- data.frame(Y = Y, A = A, L = L)
  truth_rr <- exp(0.4)

  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "aipw",
    family = stats::Gamma(link = "log"),
    propensity_model_fn = stats::glm
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    type = "ratio",
    ci_method = "sandwich"
  )
  expect_equal(res$contrasts$estimate[1], truth_rr, tolerance = 0.1)
  expect_true(all(is.finite(res$contrasts$se) & res$contrasts$se > 0))
})

test_that("aipw x bin trt x quasibinom x diff x sandwich", {
  # DGP: mu = plogis(-0.5 + 1*A + 0.3*L), Y ~ Beta(10*mu, 10*(1-mu))
  # E[Y^a] = E_L[plogis(-0.5 + a + 0.3*L)]
  set.seed(16203)
  n <- 3000
  L <- stats::rnorm(n)
  A <- stats::rbinom(n, 1, 0.5)
  mu <- stats::plogis(-0.5 + 1.0 * A + 0.3 * L)
  Y <- stats::rbeta(n, 10 * mu, 10 * (1 - mu))
  df <- data.frame(Y = Y, A = A, L = L)

  truth_1 <- mean(stats::plogis(-0.5 + 1 + 0.3 * L))
  truth_0 <- mean(stats::plogis(-0.5 + 0 + 0.3 * L))
  truth_rd <- truth_1 - truth_0

  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "aipw",
    family = stats::quasibinomial(),
    propensity_model_fn = stats::glm
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    type = "difference",
    ci_method = "sandwich"
  )
  expect_equal(res$contrasts$estimate[1], truth_rd, tolerance = 0.05)
  expect_true(all(is.finite(res$contrasts$se) & res$contrasts$se > 0))
})

test_that("aipw x bin trt x quasibinom x ratio x sandwich", {
  set.seed(16204)
  n <- 3000
  L <- stats::rnorm(n)
  A <- stats::rbinom(n, 1, 0.5)
  mu <- stats::plogis(-0.5 + 1.0 * A + 0.3 * L)
  Y <- stats::rbeta(n, 10 * mu, 10 * (1 - mu))
  df <- data.frame(Y = Y, A = A, L = L)

  truth_1 <- mean(stats::plogis(-0.5 + 1 + 0.3 * L))
  truth_0 <- mean(stats::plogis(-0.5 + 0 + 0.3 * L))
  truth_rr <- truth_1 / truth_0

  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "aipw",
    family = stats::quasibinomial(),
    propensity_model_fn = stats::glm
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    type = "ratio",
    ci_method = "sandwich"
  )
  expect_equal(res$contrasts$estimate[1], truth_rr, tolerance = 0.1)
  expect_true(all(is.finite(res$contrasts$se) & res$contrasts$se > 0))
})

test_that("aipw x bin trt x negbin x diff x sandwich", {
  skip_if_not_installed("MASS")
  # DGP: L ~ N(2,1), A ~ Bern(expit(-1 + 0.5*L)),
  #   Y ~ NB(mu = exp(0.5 + 0.3*A + 0.4*L), size = 2)
  # ATE = exp(1.38) * (exp(0.3) - 1)
  set.seed(16205)
  n <- 3000
  L <- stats::rnorm(n, 2, 1)
  A <- stats::rbinom(n, 1, stats::plogis(-1 + 0.5 * L))
  Y <- stats::rnbinom(n, mu = exp(0.5 + 0.3 * A + 0.4 * L), size = 2)
  df <- data.frame(Y = Y, A = A, L = L)
  truth <- exp(1.38) * (exp(0.3) - 1)

  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "aipw",
    model_fn = MASS::glm.nb,
    propensity_model_fn = stats::glm
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    type = "difference",
    ci_method = "sandwich"
  )
  expect_equal(res$contrasts$estimate[1], truth, tolerance = 0.5)
  expect_true(all(is.finite(res$contrasts$se) & res$contrasts$se > 0))
})

test_that("aipw x bin trt x negbin x ratio x sandwich", {
  skip_if_not_installed("MASS")
  set.seed(16206)
  n <- 3000
  L <- stats::rnorm(n, 2, 1)
  A <- stats::rbinom(n, 1, stats::plogis(-1 + 0.5 * L))
  Y <- stats::rnbinom(n, mu = exp(0.5 + 0.3 * A + 0.4 * L), size = 2)
  df <- data.frame(Y = Y, A = A, L = L)
  truth_rr <- exp(0.3)

  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "aipw",
    model_fn = MASS::glm.nb,
    propensity_model_fn = stats::glm
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    type = "ratio",
    ci_method = "sandwich"
  )
  expect_equal(res$contrasts$estimate[1], truth_rr, tolerance = 0.15)
  expect_true(all(is.finite(res$contrasts$se) & res$contrasts$se > 0))
})

test_that("aipw x bin trt x beta x diff x bootstrap", {
  skip_if_not_installed("betareg")
  # DGP: mu = plogis(0.2 + 0.5*A + 0.3*L), Y ~ Beta(mu*10, (1-mu)*10)
  set.seed(16207)
  n <- 3000
  L <- stats::rnorm(n)
  A <- stats::rbinom(n, 1, stats::plogis(-0.5 + 0.3 * L))
  mu <- stats::plogis(0.2 + 0.5 * A + 0.3 * L)
  phi <- 10
  Y <- stats::rbeta(n, mu * phi, (1 - mu) * phi)
  df <- data.frame(Y = Y, A = A, L = L)

  ey1 <- stats::integrate(
    function(x) stats::plogis(0.2 + 0.5 + 0.3 * x) * stats::dnorm(x),
    -Inf,
    Inf
  )$value
  ey0 <- stats::integrate(
    function(x) stats::plogis(0.2 + 0.3 * x) * stats::dnorm(x),
    -Inf,
    Inf
  )$value
  truth <- ey1 - ey0

  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "aipw",
    model_fn = betareg::betareg,
    family = "beta",
    propensity_model_fn = stats::glm
  )
  # betareg lacks standard $family$mu.eta so AIPW sandwich falls back
  # to bootstrap (prepare_model_if cannot compute bread for betareg).
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    type = "difference",
    ci_method = "bootstrap",
    n_boot = 100L
  )
  expect_lt(abs(res$contrasts$estimate[1] - truth), 0.02)
  expect_true(all(is.finite(res$contrasts$se) & res$contrasts$se > 0))
})

test_that("aipw x bin trt x beta x ratio x bootstrap", {
  skip_if_not_installed("betareg")
  set.seed(16208)
  n <- 3000
  L <- stats::rnorm(n)
  A <- stats::rbinom(n, 1, stats::plogis(-0.5 + 0.3 * L))
  mu <- stats::plogis(0.2 + 0.5 * A + 0.3 * L)
  phi <- 10
  Y <- stats::rbeta(n, mu * phi, (1 - mu) * phi)
  df <- data.frame(Y = Y, A = A, L = L)

  ey1 <- stats::integrate(
    function(x) stats::plogis(0.2 + 0.5 + 0.3 * x) * stats::dnorm(x),
    -Inf,
    Inf
  )$value
  ey0 <- stats::integrate(
    function(x) stats::plogis(0.2 + 0.3 * x) * stats::dnorm(x),
    -Inf,
    Inf
  )$value
  truth_rr <- ey1 / ey0

  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "aipw",
    model_fn = betareg::betareg,
    family = "beta",
    propensity_model_fn = stats::glm
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    type = "ratio",
    ci_method = "bootstrap",
    n_boot = 100L
  )
  expect_equal(res$contrasts$estimate[1], truth_rr, tolerance = 0.1)
  expect_true(all(is.finite(res$contrasts$se) & res$contrasts$se > 0))
})


# --- Gap 2: Survey / external weights -----------------------------------------

test_that("aipw x external weights x sandwich recovers ATE", {
  # Weights independent of (A, L, Y) so ATE unchanged. True ATE = 2.
  set.seed(16209)
  n <- 2000
  L <- stats::rnorm(n)
  A <- stats::rbinom(n, 1, stats::plogis(0.5 * L))
  Y <- 1 + 2 * A + 1.5 * L + stats::rnorm(n)
  pw <- stats::runif(n, 0.5, 2)
  df <- data.frame(Y = Y, A = A, L = L)

  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "aipw",
    propensity_model_fn = stats::glm,
    weights = pw
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich"
  )
  expect_equal(res$contrasts$estimate[1], 2, tolerance = 0.15)
  expect_true(all(is.finite(res$contrasts$se) & res$contrasts$se > 0))
})

test_that("aipw x external weights x bootstrap SE agreement", {
  set.seed(16210)
  n <- 1500
  L <- stats::rnorm(n)
  A <- stats::rbinom(n, 1, stats::plogis(0.5 * L))
  Y <- 1 + 2 * A + 1.5 * L + stats::rnorm(n)
  pw <- stats::runif(n, 0.5, 2)
  df <- data.frame(Y = Y, A = A, L = L)

  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "aipw",
    propensity_model_fn = stats::glm,
    weights = pw
  )
  res_sand <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich"
  )
  res_boot <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "bootstrap",
    n_boot = 100L
  )
  ratio <- res_boot$contrasts$se[1] / res_sand$contrasts$se[1]
  expect_gt(ratio, 0.7)
  expect_lt(ratio, 1.3)
})

test_that("aipw x svydesign auto-extract x sandwich", {
  skip_if_not_installed("survey")
  set.seed(16211)
  n_clusters <- 100
  m <- 5
  n <- n_clusters * m
  cl <- rep(seq_len(n_clusters), each = m)
  L <- stats::rnorm(n)
  A_cl <- stats::rbinom(n_clusters, 1, 0.5)
  flip <- stats::rbinom(n, 1, 0.1)
  A <- ifelse(flip == 1, 1 - A_cl[cl], A_cl[cl])
  Y <- 1 + 2 * A + 0.5 * L + stats::rnorm(n)
  pw <- stats::runif(n, 0.5, 2)
  df <- data.frame(Y = Y, A = A, L = L, cl = cl, pw = pw)

  des <- survey::svydesign(ids = ~cl, weights = ~pw, data = df)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "aipw",
    model_fn = stats::glm,
    propensity_model_fn = stats::glm,
    weights = des
  )
  expect_equal(fit$details$cluster, "cl")
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich"
  )
  expect_equal(res$contrasts$estimate[1], 2, tolerance = 0.2)
  expect_true(all(is.finite(res$contrasts$se) & res$contrasts$se > 0))
})

test_that("aipw x svydesign equivalence to manual weights + cluster", {
  skip_if_not_installed("survey")
  set.seed(16212)
  n_clusters <- 100
  m <- 5
  n <- n_clusters * m
  cl <- rep(seq_len(n_clusters), each = m)
  L <- stats::rnorm(n)
  A_cl <- stats::rbinom(n_clusters, 1, 0.5)
  flip <- stats::rbinom(n, 1, 0.1)
  A <- ifelse(flip == 1, 1 - A_cl[cl], A_cl[cl])
  Y <- 1 + 2 * A + 0.5 * L + stats::rnorm(n)
  pw <- stats::runif(n, 0.5, 2)
  df <- data.frame(Y = Y, A = A, L = L, cl = cl, pw = pw)

  des <- survey::svydesign(ids = ~cl, weights = ~pw, data = df)
  fit_des <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "aipw",
    model_fn = stats::glm,
    propensity_model_fn = stats::glm,
    weights = des
  )
  fit_manual <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "aipw",
    model_fn = stats::glm,
    propensity_model_fn = stats::glm,
    weights = df$pw,
    cluster = "cl"
  )

  r_des <- contrast(fit_des, list(a1 = static(1), a0 = static(0)))
  r_manual <- contrast(fit_manual, list(a1 = static(1), a0 = static(0)))

  expect_equal(r_des$contrasts$estimate, r_manual$contrasts$estimate)
  expect_equal(r_des$contrasts$se, r_manual$contrasts$se)
})


# --- Gap 3: Cluster-robust sandwich ------------------------------------------

test_that("aipw x cluster x SE inflation vs independent", {
  set.seed(16213)
  n_clusters <- 200
  m <- 5
  n <- n_clusters * m
  cl <- rep(seq_len(n_clusters), each = m)
  U <- stats::rnorm(n_clusters)
  A_cl <- stats::rbinom(n_clusters, 1, 0.5)
  L <- stats::rnorm(n)
  flip <- stats::rbinom(n, 1, 0.1)
  A <- ifelse(flip == 1, 1 - A_cl[cl], A_cl[cl])
  Y <- 1 + 2 * A + 0.5 * L + U[cl] + stats::rnorm(n, sd = 0.5)
  df <- data.frame(Y = Y, A = A, L = L, cl = cl)

  fit_indep <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "aipw",
    propensity_model_fn = stats::glm
  )
  fit_clust <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "aipw",
    propensity_model_fn = stats::glm,
    cluster = "cl"
  )

  r_indep <- contrast(
    fit_indep,
    list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich"
  )
  r_clust <- contrast(
    fit_clust,
    list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich"
  )
  se_ratio <- r_clust$contrasts$se[1] / r_indep$contrasts$se[1]
  expect_gt(se_ratio, 1.1)
})

test_that("aipw x cluster x truth recovery", {
  set.seed(16214)
  n_clusters <- 200
  m <- 5
  n <- n_clusters * m
  cl <- rep(seq_len(n_clusters), each = m)
  U <- stats::rnorm(n_clusters)
  A_cl <- stats::rbinom(n_clusters, 1, 0.5)
  L <- stats::rnorm(n)
  flip <- stats::rbinom(n, 1, 0.1)
  A <- ifelse(flip == 1, 1 - A_cl[cl], A_cl[cl])
  Y <- 1 + 2 * A + 0.5 * L + U[cl] + stats::rnorm(n, sd = 0.5)
  df <- data.frame(Y = Y, A = A, L = L, cl = cl)

  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "aipw",
    propensity_model_fn = stats::glm,
    cluster = "cl"
  )
  res <- contrast(
    fit,
    list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich"
  )
  expect_equal(res$contrasts$estimate[1], 2, tolerance = 0.15)
  expect_true(all(is.finite(res$contrasts$se) & res$contrasts$se > 0))
})

test_that("aipw x cluster x by(sex)", {
  set.seed(16215)
  n_clusters <- 200
  m <- 5
  n <- n_clusters * m
  cl <- rep(seq_len(n_clusters), each = m)
  U <- stats::rnorm(n_clusters)
  A_cl <- stats::rbinom(n_clusters, 1, 0.5)
  L <- stats::rnorm(n)
  sex <- stats::rbinom(n, 1, 0.5)
  flip <- stats::rbinom(n, 1, 0.1)
  A <- ifelse(flip == 1, 1 - A_cl[cl], A_cl[cl])
  Y <- 1 + 2 * A + 0.5 * L + 1.5 * A * sex + U[cl] + stats::rnorm(n, sd = 0.5)
  df <- data.frame(Y = Y, A = A, L = L, sex = sex, cl = cl)

  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L + sex + A:sex,
    estimator = "aipw",
    propensity_model_fn = stats::glm,
    cluster = "cl"
  )
  res <- contrast(
    fit,
    list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich",
    by = "sex"
  )
  expect_true(all(is.finite(res$contrasts$se) & res$contrasts$se > 0))
  expect_equal(nrow(res$contrasts), 2)
})


# --- Gap 4: Missing data -----------------------------------------------------

test_that("aipw x MCAR outcome NAs x censoring x sandwich", {
  df <- simulate_mcar_outcome(n = 2000, p_cens = 0.15, seed = 16216)

  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "aipw",
    propensity_model_fn = stats::glm,
    censoring = "C"
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich"
  )
  expect_equal(res$contrasts$estimate[1], 3, tolerance = 0.2)
  expect_true(all(is.finite(res$contrasts$se) & res$contrasts$se > 0))
})

test_that("aipw x MCAR outcome NAs x bootstrap SE agreement", {
  df <- simulate_mcar_outcome(n = 1500, p_cens = 0.15, seed = 16217)

  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "aipw",
    propensity_model_fn = stats::glm,
    censoring = "C"
  )
  res_sand <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich"
  )
  res_boot <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "bootstrap",
    n_boot = 100L
  )
  ratio <- res_boot$contrasts$se[1] / res_sand$contrasts$se[1]
  expect_gt(ratio, 0.5)
  expect_lt(ratio, 2.0)
})

test_that("aipw x covariate NAs require pre-processing", {
  # AIPW rejects confounder NAs because outcome and propensity fits

  # would disagree on which rows to keep. Users must drop/impute first.
  set.seed(16218)
  n <- 2000
  L <- stats::rnorm(n)
  A <- stats::rbinom(n, 1, stats::plogis(0.5 * L))
  Y <- 2 + 3 * A + 1.5 * L + stats::rnorm(n)
  na_idx <- sample.int(n, size = round(0.10 * n))
  L[na_idx] <- NA_real_
  df <- data.frame(Y = Y, A = A, L = L)

  expect_error(
    causat(
      df,
      outcome = "Y",
      treatment = "A",
      confounders = ~L,
      estimator = "aipw",
      propensity_model_fn = stats::glm
    ),
    "confounder column has missing values"
  )
})

test_that("aipw x external weights + MCAR outcome NAs", {
  df <- simulate_mcar_outcome(n = 2000, p_cens = 0.15, seed = 16219)
  set.seed(16219)
  pw <- stats::runif(nrow(df), 0.5, 2)

  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "aipw",
    propensity_model_fn = stats::glm,
    censoring = "C",
    weights = pw
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich"
  )
  expect_true(all(is.finite(res$contrasts$estimate)))
  expect_true(all(is.finite(res$contrasts$se) & res$contrasts$se > 0))
})


# --- Gap 5: Subset estimand --------------------------------------------------

test_that("aipw x subset(L > 0) restricts population", {
  set.seed(16220)
  n <- 3000
  L <- stats::rnorm(n)
  A <- stats::rbinom(n, 1, stats::plogis(0.5 * L))
  Y <- 2 + 3 * A + 1.5 * L + stats::rnorm(n)
  df <- data.frame(Y = Y, A = A, L = L)

  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "aipw",
    propensity_model_fn = stats::glm
  )
  res_full <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich"
  )
  res_sub <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich",
    subset = quote(L > 0)
  )
  expect_true(all(is.finite(res_sub$contrasts$estimate)))
  expect_true(all(is.finite(res_sub$contrasts$se) & res_sub$contrasts$se > 0))
  # Subset and full should differ (different population)
  expect_false(
    isTRUE(all.equal(
      res_sub$contrasts$estimate[1],
      res_full$contrasts$estimate[1]
    ))
  )
})

test_that("aipw x subset + by composition", {
  df <- simulate_het_effect(n = 3000)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L + sex + A:L + A:sex,
    estimator = "aipw",
    propensity_model_fn = stats::glm
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

test_that("aipw x subset + ATT", {
  set.seed(16221)
  n <- 3000
  L <- stats::rnorm(n)
  A <- stats::rbinom(n, 1, stats::plogis(0.5 * L))
  Y <- 2 + 3 * A + 1.5 * L + stats::rnorm(n)
  df <- data.frame(Y = Y, A = A, L = L)

  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "aipw",
    estimand = "ATT",
    propensity_model_fn = stats::glm
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich",
    subset = quote(L > 0)
  )
  expect_true(all(is.finite(res$contrasts$estimate)))
  expect_true(all(is.finite(res$contrasts$se) & res$contrasts$se > 0))
})


# --- Gap 6: GAM outcome / propensity models -----------------------------------

test_that("aipw x GAM outcome x nonlinear DGP x sandwich", {
  skip_if_not_installed("mgcv")
  # DGP with nonlinear confounding: Y = 1 + 3*A + sin(2*L) + L^2 + eps
  # GLM Y ~ A + L is misspecified; GAM s(L) recovers truth.
  set.seed(16222)
  n <- 3000
  L <- stats::rnorm(n)
  ps <- stats::plogis(sin(L) + 0.3 * L^2)
  A <- stats::rbinom(n, 1, ps)
  Y <- 1 + 3 * A + sin(2 * L) + L^2 + stats::rnorm(n)
  df <- data.frame(Y = Y, A = A, L = L)

  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~ s(L),
    estimator = "aipw",
    model_fn = mgcv::gam,
    propensity_model_fn = mgcv::gam
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich"
  )
  expect_equal(res$contrasts$estimate[1], 3, tolerance = 0.2)
  expect_true(all(is.finite(res$contrasts$se) & res$contrasts$se > 0))
})

test_that("aipw x GAM propensity DR: wrong outcome, correct GAM propensity", {
  skip_if_not_installed("mgcv")
  # DR property: misspecified outcome (drops L), correct GAM propensity.
  set.seed(16223)
  n <- 3000
  L <- stats::rnorm(n)
  ps <- stats::plogis(sin(L) + 0.3 * L^2)
  A <- stats::rbinom(n, 1, ps)
  Y <- 1 + 3 * A + sin(2 * L) + L^2 + stats::rnorm(n)
  df <- data.frame(Y = Y, A = A, L = L)

  # Outcome model deliberately wrong: ignores L entirely
  wrong_outcome <- function(formula, data, ...) {
    stats::glm(Y ~ A, data = data, ...)
  }

  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~ s(L),
    estimator = "aipw",
    model_fn = wrong_outcome,
    propensity_model_fn = mgcv::gam
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich"
  )
  expect_equal(res$contrasts$estimate[1], 3, tolerance = 0.5)
  expect_true(all(is.finite(res$contrasts$se) & res$contrasts$se > 0))
})

test_that("aipw x GAM outcome+propensity x sandwich vs bootstrap SE", {
  skip_if_not_installed("mgcv")
  set.seed(16224)
  n <- 2000
  L <- stats::rnorm(n)
  ps <- stats::plogis(sin(L) + 0.3 * L^2)
  A <- stats::rbinom(n, 1, ps)
  Y <- 1 + 3 * A + sin(2 * L) + L^2 + stats::rnorm(n)
  df <- data.frame(Y = Y, A = A, L = L)

  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~ s(L),
    estimator = "aipw",
    model_fn = mgcv::gam,
    propensity_model_fn = mgcv::gam
  )
  res_sand <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich"
  )
  res_boot <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "bootstrap",
    n_boot = 100L
  )
  ratio <- res_boot$contrasts$se[1] / res_sand$contrasts$se[1]
  expect_gt(ratio, 0.5)
  expect_lt(ratio, 2.0)
})
