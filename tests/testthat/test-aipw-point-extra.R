# AIPW point estimator: multivariate treatment, non-Gaussian families,
# external/survey/cluster weights, missing data, subset, GAM.
# Split from test-aipw.R for test-file-level parallelism.

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
  expect_equal(ate, 2.5, tolerance = 0.05)
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
  skip_if_fast()
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
  expect_equal(ate, -0.9, tolerance = 0.05)
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
  expect_equal(ate, 2.5, tolerance = 0.05)
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
  expect_equal(ate, 2.5, tolerance = 0.05)
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
  expect_equal(ate, 2.5, tolerance = 0.05)
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
  # NegBin has high variance; n=3000 is ~13% relative error, converges at n=30K
  expect_equal(res$contrasts$estimate[1], truth, tolerance = 0.15)
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
  skip_if_fast()
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
  skip_if_fast()
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
  skip_if_fast()
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
  skip_if_fast()
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
  expect_equal(res$contrasts$estimate[1], 3, tolerance = 0.05)
  expect_true(all(is.finite(res$contrasts$se) & res$contrasts$se > 0))
})

test_that("aipw x GAM outcome+propensity x sandwich vs bootstrap SE", {
  skip_if_fast()
  skip_if_not_installed("mgcv")
  set.seed(16224)
  n <- 2000
  L <- stats::rnorm(n)
  ps <- stats::plogis(sin(L) + 0.3 * L^2)
  A <- stats::rbinom(n, 1, ps)
  Y <- 1 + 3 * A + sin(2 * L) + L^2 + stats::rnorm(n)
  df <- data.frame(Y = Y, A = A, L = L)

  # mgcv's GAM optimiser occasionally backs off a step on a bootstrap resample
  # ("Fitting terminated with step failure"); that is mgcv-internal numerical
  # chatter, not a causatr problem, so muffle just that message (see
  # helper-external-warnings.R) around the GAM fit and the refitting bootstrap.
  fit <- muffle_gam_step_failure(causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~ s(L),
    estimator = "aipw",
    model_fn = mgcv::gam,
    propensity_model_fn = mgcv::gam
  ))
  res_sand <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich"
  )
  res_boot <- muffle_gam_step_failure(contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "bootstrap",
    n_boot = 100L
  ))
  ratio <- res_boot$contrasts$se[1] / res_sand$contrasts$se[1]
  expect_gt(ratio, 0.5)
  expect_lt(ratio, 2.0)
})
