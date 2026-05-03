# -- Chunk 1: Constructor + rejection paths ------------------------------------

test_that("stochastic() creates a causatr_intervention", {
  sampler <- \(data, trt) rbinom(nrow(data), 1, 0.5)
  iv <- stochastic(sampler, n_mc = 50L)
  expect_s3_class(iv, "causatr_intervention")
  expect_equal(iv$type, "stochastic")
  expect_true(is.function(iv$sampler))
  expect_equal(iv$n_mc, 50L)
})

test_that("stochastic() rejects non-function sampler", {
  expect_error(stochastic(42), "function")
  expect_error(stochastic("foo"), "function")
})

test_that("stochastic() rejects bad n_mc", {
  s <- \(d, t) t
  expect_error(stochastic(s, n_mc = 0L), "positive integer")
  expect_error(stochastic(s, n_mc = -1L), "positive integer")
  expect_error(stochastic(s, n_mc = "a"), "positive integer")
  expect_error(stochastic(s, n_mc = NA), "positive integer")
})

test_that("stochastic() warns when n_mc < 10", {
  s <- \(d, t) t
  expect_warning(stochastic(s, n_mc = 5L), "very low")
  expect_warning(stochastic(s, n_mc = 1L), "very low")
  expect_no_warning(stochastic(s, n_mc = 10L))
})

test_that("stochastic() coerces n_mc to integer", {
  s <- \(d, t) t
  iv <- stochastic(s, n_mc = 100)
  expect_identical(iv$n_mc, 100L)
})

test_that("print.causatr_intervention displays stochastic correctly", {
  iv <- stochastic(\(d, t) t, n_mc = 15L)
  expect_output(print(iv), "stochastic")
})

test_that("apply_single_intervention works for stochastic (numeric)", {
  dt <- data.table::data.table(A = c(0, 1, 0, 1), L = c(1, 2, 3, 4))
  sampler <- \(data, trt) rep(1, nrow(data))
  iv <- stochastic(sampler, n_mc = 10L)
  apply_single_intervention(dt, "A", iv)
  expect_equal(dt$A, rep(1, 4))
})

test_that("apply_single_intervention works for stochastic (factor)", {
  dt <- data.table::data.table(
    A = factor(c("a", "b", "a"), levels = c("a", "b", "c")),
    L = 1:3
  )
  sampler <- \(data, trt) c("b", "c", "a")
  iv <- stochastic(sampler, n_mc = 10L)
  apply_single_intervention(dt, "A", iv)
  expect_equal(as.character(dt$A), c("b", "c", "a"))
  expect_equal(levels(dt$A), c("a", "b", "c"))
})

test_that("apply_single_intervention rejects stochastic type mismatch", {
  dt <- data.table::data.table(A = c(0, 1, 0), L = 1:3)
  bad_sampler <- \(data, trt) c("a", "b", "c")
  iv <- stochastic(bad_sampler, n_mc = 10L)
  expect_error(
    apply_single_intervention(dt, "A", iv),
    "stochastic.*character.*numeric"
  )
})

test_that("apply_single_intervention rejects stochastic wrong length", {
  dt <- data.table::data.table(A = c(0, 1, 0), L = 1:3)
  bad_sampler <- \(data, trt) c(1, 2)
  iv <- stochastic(bad_sampler, n_mc = 10L)
  expect_error(
    apply_single_intervention(dt, "A", iv),
    "length 3.*got 2"
  )
})

test_that("stochastic() rejected under IPW", {
  skip_on_cran()
  set.seed(42)
  n <- 200
  L <- rnorm(n)
  A <- rbinom(n, 1, plogis(L))
  Y <- 2 + A + L + rnorm(n)
  df <- data.frame(Y = Y, A = A, L = L)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw"
  )
  sampler <- \(data, trt) rbinom(nrow(data), 1, 0.5)
  expect_error(
    contrast(fit, interventions = list(g = stochastic(sampler))),
    "stochastic.*gcomp"
  )
})

test_that("stochastic() rejected under matching", {
  skip_on_cran()
  set.seed(42)
  n <- 200
  L <- rnorm(n)
  A <- rbinom(n, 1, plogis(L))
  Y <- 2 + A + L + rnorm(n)
  df <- data.frame(Y = Y, A = A, L = L)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "matching"
  )
  sampler <- \(data, trt) rbinom(nrow(data), 1, 0.5)
  expect_error(
    contrast(fit, interventions = list(g = stochastic(sampler))),
    "Non-static.*stochastic"
  )
})

# -- Chunk 2: Point gcomp MC integration ------------------------------------

test_that("stochastic gcomp: binary treatment, gaussian outcome", {
  skip_on_cran()
  dgp <- simulate_stochastic_binary_gaussian(n = 1000, seed = 42)
  fit <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "gcomp"
  )
  set.seed(123)
  res <- contrast(
    fit,
    interventions = list(g = stochastic(dgp$sampler, n_mc = 20L)),
    type = "difference",
    ci_method = "sandwich"
  )
  mu_g <- unname(coef(res)["g"])
  expect_equal(mu_g, dgp$truth, tolerance = 0.35)
  ci <- confint(res)
  expect_true(ci["g", 1] <= dgp$truth && dgp$truth <= ci["g", 2])
})

test_that("stochastic gcomp: binary/gaussian, bootstrap CI covers truth", {
  skip_on_cran()
  dgp <- simulate_stochastic_binary_gaussian(n = 200, seed = 99)
  fit <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "gcomp"
  )
  set.seed(456)
  res <- contrast(
    fit,
    interventions = list(g = stochastic(dgp$sampler, n_mc = 5L)),
    type = "difference",
    ci_method = "bootstrap",
    n_boot = 20
  )
  ci <- confint(res)
  expect_true(ci["g", 1] <= dgp$truth && dgp$truth <= ci["g", 2])
})

test_that("stochastic gcomp: continuous treatment, gaussian outcome", {
  skip_on_cran()
  dgp <- simulate_stochastic_continuous_gaussian(n = 1000, seed = 42)
  fit <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "gcomp"
  )
  set.seed(123)
  res <- contrast(
    fit,
    interventions = list(g = stochastic(dgp$sampler, n_mc = 20L)),
    type = "difference",
    ci_method = "sandwich"
  )
  mu_g <- unname(coef(res)["g"])
  expect_equal(mu_g, dgp$truth, tolerance = 0.35)
  ci <- confint(res)
  expect_true(ci["g", 1] <= dgp$truth && dgp$truth <= ci["g", 2])
})

test_that("stochastic gcomp: binary treatment, binomial outcome", {
  skip_on_cran()
  dgp <- simulate_stochastic_binary_binomial(n = 1000, seed = 42)
  fit <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "gcomp",
    family = binomial()
  )
  set.seed(123)
  res <- contrast(
    fit,
    interventions = list(g = stochastic(dgp$sampler, n_mc = 20L)),
    type = "difference",
    ci_method = "sandwich"
  )
  mu_g <- unname(coef(res)["g"])
  expect_equal(mu_g, dgp$truth, tolerance = 0.1)
  ci <- confint(res)
  expect_true(ci["g", 1] <= dgp$truth && dgp$truth <= ci["g", 2])
})

test_that("stochastic gcomp: categorical treatment, gaussian outcome", {
  skip_on_cran()
  dgp <- simulate_stochastic_categorical_gaussian(n = 1000, seed = 42)
  fit <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "gcomp"
  )
  set.seed(123)
  res <- contrast(
    fit,
    interventions = list(g = stochastic(dgp$sampler, n_mc = 20L)),
    type = "difference",
    ci_method = "sandwich"
  )
  mu_g <- unname(coef(res)["g"])
  expect_equal(mu_g, dgp$truth, tolerance = 0.35)
  ci <- confint(res)
  expect_true(ci["g", 1] <= dgp$truth && dgp$truth <= ci["g", 2])
})

test_that("stochastic gcomp: multivariate treatment, gaussian outcome", {
  skip_on_cran()
  dgp <- simulate_stochastic_multivariate_gaussian(n = 1000, seed = 42)
  fit <- causat(
    dgp$data,
    outcome = "Y",
    treatment = c("A1", "A2"),
    confounders = ~L,
    estimator = "gcomp"
  )
  set.seed(123)
  res <- contrast(
    fit,
    interventions = list(
      g = list(
        A1 = stochastic(dgp$sampler_a1, n_mc = 20L),
        A2 = stochastic(dgp$sampler_a2, n_mc = 20L)
      )
    ),
    type = "difference",
    ci_method = "sandwich"
  )
  mu_g <- unname(coef(res)["g"])
  expect_equal(mu_g, dgp$truth, tolerance = 0.35)
})

test_that("stochastic gcomp: binomial outcome, ratio contrast", {
  skip_on_cran()
  dgp <- simulate_stochastic_binary_binomial(n = 1000, seed = 42)
  fit <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "gcomp",
    family = binomial()
  )
  # Compare stochastic vs static(0) on the ratio scale
  set.seed(123)
  res <- contrast(
    fit,
    interventions = list(
      g = stochastic(dgp$sampler, n_mc = 20L),
      a0 = static(0)
    ),
    type = "ratio",
    ci_method = "sandwich"
  )
  # ratio = mu_g / mu_a0 — just check it's positive and finite
  est <- coef(res)
  expect_true(all(is.finite(est)))
  expect_true(est["g"] > 0 && est["a0"] > 0)
})

test_that("stochastic gcomp: binomial outcome, OR contrast", {
  skip_on_cran()
  dgp <- simulate_stochastic_binary_binomial(n = 1000, seed = 42)
  fit <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "gcomp",
    family = binomial()
  )
  set.seed(123)
  res <- contrast(
    fit,
    interventions = list(
      g = stochastic(dgp$sampler, n_mc = 20L),
      a0 = static(0)
    ),
    type = "or",
    ci_method = "sandwich"
  )
  est <- coef(res)
  expect_true(all(is.finite(est)))
})

test_that("stochastic gcomp: ATT estimand", {
  skip_on_cran()
  dgp <- simulate_stochastic_binary_gaussian(n = 1000, seed = 42)
  fit <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "gcomp"
  )
  set.seed(123)
  res <- contrast(
    fit,
    interventions = list(g = stochastic(dgp$sampler, n_mc = 20L)),
    type = "difference",
    ci_method = "sandwich",
    estimand = "ATT"
  )
  mu_g <- unname(coef(res)["g"])
  expect_true(is.finite(mu_g))
  ci <- confint(res)
  expect_true(ci["g", 1] < ci["g", 2])
})

test_that("stochastic gcomp: ATC estimand", {
  skip_on_cran()
  dgp <- simulate_stochastic_binary_gaussian(n = 1000, seed = 42)
  fit <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "gcomp"
  )
  set.seed(123)
  res <- contrast(
    fit,
    interventions = list(g = stochastic(dgp$sampler, n_mc = 20L)),
    type = "difference",
    ci_method = "sandwich",
    estimand = "ATC"
  )
  mu_g <- unname(coef(res)["g"])
  expect_true(is.finite(mu_g))
  ci <- confint(res)
  expect_true(ci["g", 1] < ci["g", 2])
})

test_that("stochastic gcomp: by-stratified estimand", {
  skip_on_cran()
  dgp <- simulate_stochastic_binary_gaussian(n = 1000, seed = 42)
  dgp$data$sex <- rbinom(nrow(dgp$data), 1, 0.5)
  fit <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L + sex,
    estimator = "gcomp"
  )
  set.seed(123)
  res <- contrast(
    fit,
    interventions = list(g = stochastic(dgp$sampler, n_mc = 15L)),
    type = "difference",
    ci_method = "sandwich",
    by = "sex"
  )
  expect_s3_class(res, "causatr_result")
  expect_equal(nrow(res$estimates), 2)
})

test_that("stochastic gcomp: subset estimand", {
  skip_on_cran()
  dgp <- simulate_stochastic_binary_gaussian(n = 1000, seed = 42)
  fit <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "gcomp"
  )
  set.seed(123)
  res <- contrast(
    fit,
    interventions = list(g = stochastic(dgp$sampler, n_mc = 15L)),
    type = "difference",
    ci_method = "sandwich",
    subset = quote(L > 0)
  )
  mu_g <- unname(coef(res)["g"])
  expect_true(is.finite(mu_g))
})

test_that("stochastic gcomp: GAM model", {
  skip_on_cran()
  skip_if_not_installed("mgcv")
  dgp <- simulate_stochastic_binary_gaussian(n = 1000, seed = 42)
  fit <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "gcomp",
    model_fn = mgcv::gam
  )
  set.seed(123)
  res <- contrast(
    fit,
    interventions = list(g = stochastic(dgp$sampler, n_mc = 15L)),
    type = "difference",
    ci_method = "sandwich"
  )
  mu_g <- unname(coef(res)["g"])
  expect_equal(mu_g, dgp$truth, tolerance = 0.3)
})

test_that("stochastic gcomp: Poisson outcome, ratio contrast", {
  skip_on_cran()
  set.seed(42)
  n <- 5000
  L <- rnorm(n)
  A <- rbinom(n, 1, plogis(0.5 * L))
  Y <- rpois(n, exp(0.5 + 0.3 * A + 0.2 * L))
  df <- data.frame(Y = Y, A = A, L = L)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "gcomp",
    family = poisson()
  )
  sampler <- \(data, trt) rbinom(nrow(data), 1, 0.5)
  set.seed(123)
  res <- contrast(
    fit,
    interventions = list(g = stochastic(sampler, n_mc = 15L)),
    type = "ratio",
    ci_method = "sandwich"
  )
  mu_g <- unname(coef(res)["g"])
  expect_true(is.finite(mu_g) && mu_g > 0)
})

test_that("stochastic gcomp: n_mc = 1 works (degenerate single draw)", {
  skip_on_cran()
  dgp <- simulate_stochastic_binary_gaussian(n = 1000, seed = 42)
  fit <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "gcomp"
  )
  set.seed(99)
  expect_warning(
    iv <- stochastic(dgp$sampler, n_mc = 1L),
    "very low"
  )
  res <- contrast(
    fit,
    interventions = list(g = iv),
    type = "difference",
    ci_method = "sandwich"
  )
  mu_g <- unname(coef(res)["g"])
  expect_true(is.finite(mu_g))
})

test_that("stochastic gcomp: reproducible with set.seed()", {
  skip_on_cran()
  dgp <- simulate_stochastic_binary_gaussian(n = 1000, seed = 42)
  fit <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "gcomp"
  )
  set.seed(777)
  res1 <- contrast(
    fit,
    interventions = list(g = stochastic(dgp$sampler, n_mc = 20L)),
    type = "difference",
    ci_method = "sandwich"
  )
  set.seed(777)
  res2 <- contrast(
    fit,
    interventions = list(g = stochastic(dgp$sampler, n_mc = 20L)),
    type = "difference",
    ci_method = "sandwich"
  )
  expect_equal(coef(res1), coef(res2))
})

test_that("stochastic gcomp: continuous shift vs deterministic shift oracle", {
  skip_on_cran()
  dgp <- simulate_stochastic_continuous_gaussian(n = 500, seed = 42)
  fit <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "gcomp"
  )
  # Stochastic: A* = A + N(0.5, 0.25) — mean shift is 0.5
  set.seed(123)
  res_stoch <- contrast(
    fit,
    interventions = list(g = stochastic(dgp$sampler, n_mc = 15L)),
    type = "difference",
    ci_method = "sandwich"
  )
  # Deterministic: A* = A + 0.5 — same mean shift
  res_det <- contrast(
    fit,
    interventions = list(g = shift(0.5)),
    type = "difference",
    ci_method = "sandwich"
  )
  # For linear outcome model Y = b0 + b1*A + b2*L + e,
  # E[Y^g] depends only on E[A*] = E[A] + 0.5, so both should agree.
  expect_equal(
    unname(coef(res_stoch)["g"]),
    unname(coef(res_det)["g"]),
    tolerance = 0.3
  )
})

test_that("stochastic gcomp: mixed interventions (stochastic + static)", {
  skip_on_cran()
  dgp <- simulate_stochastic_binary_gaussian(n = 1000, seed = 42)
  fit <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "gcomp"
  )
  set.seed(123)
  res <- contrast(
    fit,
    interventions = list(
      g = stochastic(dgp$sampler, n_mc = 20L),
      a1 = static(1),
      a0 = static(0)
    ),
    type = "difference",
    ci_method = "sandwich"
  )
  mu <- coef(res)
  # mu_g should be between static(0) and static(1)
  expect_true(mu["a0"] < mu["g"] && mu["g"] < mu["a1"])
  # Static means should match known DGP truths
  expect_equal(unname(mu["a1"]), 5, tolerance = 0.35)
  expect_equal(unname(mu["a0"]), 2, tolerance = 0.35)
})

# -- Chunk 4: ICE (longitudinal) MC integration ------------------------------

test_that("stochastic ICE: binary treatment, gaussian, 2 periods", {
  skip_on_cran()
  dgp <- simulate_stochastic_ice_binary_gaussian(n = 500, seed = 42)
  fit <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    estimator = "gcomp"
  )
  set.seed(123)
  res <- contrast(
    fit,
    interventions = list(g = stochastic(dgp$sampler, n_mc = 15L)),
    type = "difference",
    ci_method = "sandwich"
  )
  mu_g <- unname(coef(res)["g"])
  expect_equal(mu_g, dgp$truth, tolerance = 0.3)
  ci <- confint(res)
  expect_true(ci["g", 1] <= dgp$truth && dgp$truth <= ci["g", 2])
})

test_that("stochastic ICE: binary/gaussian, bootstrap CI covers truth", {
  skip_on_cran()
  dgp <- simulate_stochastic_ice_binary_gaussian(n = 200, seed = 99)
  fit <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    estimator = "gcomp"
  )
  set.seed(456)
  res <- contrast(
    fit,
    interventions = list(g = stochastic(dgp$sampler, n_mc = 5L)),
    type = "difference",
    ci_method = "bootstrap",
    n_boot = 20
  )
  ci <- confint(res)
  expect_true(ci["g", 1] <= dgp$truth && dgp$truth <= ci["g", 2])
})

test_that("stochastic ICE: continuous treatment, gaussian, 2 periods", {
  skip_on_cran()
  dgp <- simulate_stochastic_ice_continuous_gaussian(n = 500, seed = 42)
  fit <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    estimator = "gcomp"
  )
  set.seed(123)
  res <- contrast(
    fit,
    interventions = list(g = stochastic(dgp$sampler, n_mc = 15L)),
    type = "difference",
    ci_method = "sandwich"
  )
  mu_g <- unname(coef(res)["g"])
  expect_equal(mu_g, dgp$truth, tolerance = 0.3)
})

test_that("stochastic ICE: reproducible with set.seed()", {
  skip_on_cran()
  dgp <- simulate_stochastic_ice_binary_gaussian(n = 500, seed = 42)
  fit <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    estimator = "gcomp"
  )
  set.seed(777)
  res1 <- contrast(
    fit,
    interventions = list(g = stochastic(dgp$sampler, n_mc = 15L)),
    type = "difference",
    ci_method = "sandwich"
  )
  set.seed(777)
  res2 <- contrast(
    fit,
    interventions = list(g = stochastic(dgp$sampler, n_mc = 15L)),
    type = "difference",
    ci_method = "sandwich"
  )
  expect_equal(coef(res1), coef(res2))
})

test_that("stochastic ICE: mixed interventions (stochastic vs static)", {
  skip_on_cran()
  dgp <- simulate_stochastic_ice_binary_gaussian(n = 500, seed = 42)
  fit <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    estimator = "gcomp"
  )
  set.seed(123)
  res <- contrast(
    fit,
    interventions = list(
      g = stochastic(dgp$sampler, n_mc = 15L),
      always = static(1),
      never = static(0)
    ),
    type = "difference",
    ci_method = "sandwich"
  )
  mu <- coef(res)
  expect_true(mu["never"] < mu["g"] && mu["g"] < mu["always"])
})

# -- Chunk 3: Sandwich vs bootstrap agreement ---------------------------------

test_that("stochastic gcomp: sandwich and bootstrap SEs agree (point)", {
  skip_on_cran()
  dgp <- simulate_stochastic_binary_gaussian(n = 2000, seed = 42)
  fit <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "gcomp"
  )
  set.seed(303)
  res_sw <- contrast(
    fit,
    interventions = list(
      g = stochastic(dgp$sampler, n_mc = 50L),
      a1 = static(1)
    ),
    type = "difference",
    ci_method = "sandwich"
  )
  set.seed(303)
  res_bs <- contrast(
    fit,
    interventions = list(
      g = stochastic(dgp$sampler, n_mc = 50L),
      a1 = static(1)
    ),
    type = "difference",
    ci_method = "bootstrap",
    n_boot = 200
  )
  est_sw <- res_sw$estimates
  est_bs <- res_bs$estimates
  # Point estimates should be identical (same seed, same n_mc)
  expect_equal(est_sw$estimate, est_bs$estimate, tolerance = 1e-10)
  # SEs should agree within a factor of 2 (bootstrap has MC noise)
  ratio_g <- est_sw$se[est_sw$intervention == "g"] /
    est_bs$se[est_bs$intervention == "g"]
  expect_gt(ratio_g, 0.5)
  expect_lt(ratio_g, 2.0)
  ratio_a1 <- est_sw$se[est_sw$intervention == "a1"] /
    est_bs$se[est_bs$intervention == "a1"]
  expect_gt(ratio_a1, 0.5)
  expect_lt(ratio_a1, 2.0)
  # CIs should both cover truth
  ci_sw <- confint(res_sw)
  ci_bs <- confint(res_bs)
  expect_true(ci_sw["g", 1] <= dgp$truth && dgp$truth <= ci_sw["g", 2])
  expect_true(ci_bs["g", 1] <= dgp$truth && dgp$truth <= ci_bs["g", 2])
})

test_that("stochastic gcomp: sandwich and bootstrap SEs agree (binomial)", {
  skip_on_cran()
  dgp <- simulate_stochastic_binary_binomial(n = 3000, seed = 42)
  fit <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "gcomp",
    family = binomial()
  )
  set.seed(404)
  res_sw <- contrast(
    fit,
    interventions = list(g = stochastic(dgp$sampler, n_mc = 50L)),
    type = "difference",
    ci_method = "sandwich"
  )
  set.seed(404)
  res_bs <- contrast(
    fit,
    interventions = list(g = stochastic(dgp$sampler, n_mc = 50L)),
    type = "difference",
    ci_method = "bootstrap",
    n_boot = 200
  )
  ratio <- res_sw$estimates$se / res_bs$estimates$se
  expect_gt(ratio, 0.5)
  expect_lt(ratio, 2.0)
})

test_that("stochastic gcomp: sandwich and bootstrap SEs agree (continuous)", {
  skip_on_cran()
  dgp <- simulate_stochastic_continuous_gaussian(n = 2000, seed = 42)
  fit <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "gcomp"
  )
  set.seed(606)
  res_sw <- contrast(
    fit,
    interventions = list(g = stochastic(dgp$sampler, n_mc = 50L)),
    type = "difference",
    ci_method = "sandwich"
  )
  set.seed(606)
  res_bs <- contrast(
    fit,
    interventions = list(g = stochastic(dgp$sampler, n_mc = 50L)),
    type = "difference",
    ci_method = "bootstrap",
    n_boot = 200
  )
  expect_equal(
    res_sw$estimates$estimate,
    res_bs$estimates$estimate,
    tolerance = 1e-10
  )
  ratio <- res_sw$estimates$se / res_bs$estimates$se
  expect_gt(ratio, 0.5)
  expect_lt(ratio, 2.0)
})

test_that("stochastic gcomp: sandwich and bootstrap SEs agree (categorical)", {
  skip_on_cran()
  dgp <- simulate_stochastic_categorical_gaussian(
    n = 2000, seed = 42
  )
  fit <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "gcomp"
  )
  set.seed(707)
  res_sw <- contrast(
    fit,
    interventions = list(g = stochastic(dgp$sampler, n_mc = 50L)),
    type = "difference",
    ci_method = "sandwich"
  )
  set.seed(707)
  res_bs <- contrast(
    fit,
    interventions = list(g = stochastic(dgp$sampler, n_mc = 50L)),
    type = "difference",
    ci_method = "bootstrap",
    n_boot = 200
  )
  expect_equal(
    res_sw$estimates$estimate,
    res_bs$estimates$estimate,
    tolerance = 1e-10
  )
  ratio <- res_sw$estimates$se / res_bs$estimates$se
  expect_gt(ratio, 0.5)
  expect_lt(ratio, 2.0)
})

test_that("stochastic gcomp: sandwich and bootstrap SEs agree (multivariate)", {
  skip_on_cran()
  dgp <- simulate_stochastic_multivariate_gaussian(
    n = 2000, seed = 42
  )
  fit <- causat(
    dgp$data,
    outcome = "Y",
    treatment = c("A1", "A2"),
    confounders = ~L,
    estimator = "gcomp"
  )
  set.seed(808)
  res_sw <- contrast(
    fit,
    interventions = list(
      g = list(
        A1 = stochastic(dgp$sampler_a1, n_mc = 50L),
        A2 = stochastic(dgp$sampler_a2, n_mc = 50L)
      )
    ),
    type = "difference",
    ci_method = "sandwich"
  )
  set.seed(808)
  res_bs <- contrast(
    fit,
    interventions = list(
      g = list(
        A1 = stochastic(dgp$sampler_a1, n_mc = 50L),
        A2 = stochastic(dgp$sampler_a2, n_mc = 50L)
      )
    ),
    type = "difference",
    ci_method = "bootstrap",
    n_boot = 200
  )
  expect_equal(
    res_sw$estimates$estimate,
    res_bs$estimates$estimate,
    tolerance = 1e-10
  )
  ratio <- res_sw$estimates$se / res_bs$estimates$se
  expect_gt(ratio, 0.5)
  expect_lt(ratio, 2.0)
})

# -- Chunk 5: Sandwich vs bootstrap agreement (ICE) ---------------------------

test_that("stochastic ICE: sandwich and bootstrap SEs agree", {
  skip_on_cran()
  dgp <- simulate_stochastic_ice_binary_gaussian(n = 1000, seed = 42)
  fit <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    estimator = "gcomp"
  )
  set.seed(505)
  res_sw <- contrast(
    fit,
    interventions = list(g = stochastic(dgp$sampler, n_mc = 30L)),
    type = "difference",
    ci_method = "sandwich"
  )
  set.seed(505)
  res_bs <- contrast(
    fit,
    interventions = list(g = stochastic(dgp$sampler, n_mc = 30L)),
    type = "difference",
    ci_method = "bootstrap",
    n_boot = 100
  )
  est_sw <- res_sw$estimates
  est_bs <- res_bs$estimates
  # Point estimates identical (same seed)
  expect_equal(est_sw$estimate, est_bs$estimate, tolerance = 1e-10)
  # SEs within a factor of 2
  ratio <- est_sw$se / est_bs$se
  expect_gt(ratio, 0.5)
  expect_lt(ratio, 2.0)
  # Both CIs cover truth
  ci_sw <- confint(res_sw)
  ci_bs <- confint(res_bs)
  expect_true(ci_sw["g", 1] <= dgp$truth && dgp$truth <= ci_sw["g", 2])
  expect_true(ci_bs["g", 1] <= dgp$truth && dgp$truth <= ci_bs["g", 2])
})

# -- Cross-package: lmtp oracle -----------------------------------------------

test_that("stochastic gcomp: agrees with lmtp_sdr (point)", {
  skip_on_cran()
  skip_if_not_installed("lmtp")
  skip_if_not_installed("SuperLearner")

  set.seed(42)
  n <- 2000
  L <- rnorm(n)
  A <- rbinom(n, 1, plogis(0.5 * L))
  Y <- 2 + 3 * A + 1.5 * L + rnorm(n)
  df <- data.frame(Y = Y, A = A, L = L)

  sampler <- function(data, trt) {
    rbinom(nrow(data), 1, plogis(0.5 + 0.3 * data$L))
  }

  # causatr
  fit <- causat(df, outcome = "Y", treatment = "A",
                confounders = ~L, estimator = "gcomp")
  set.seed(123)
  res <- contrast(
    fit,
    interventions = list(g = stochastic(sampler, n_mc = 200L)),
    type = "difference",
    ci_method = "sandwich"
  )
  est_causatr <- res$estimates$estimate

  # lmtp: the shift function for lmtp draws from the same stochastic policy
  lmtp_shift <- function(data, trt) {
    rbinom(nrow(data), 1, plogis(0.5 + 0.3 * data[["L"]]))
  }
  lmtp_res <- tryCatch(
    suppressWarnings(suppressMessages(lmtp::lmtp_sdr(
      data = df,
      trt = "A",
      outcome = "Y",
      baseline = "L",
      shift = lmtp_shift,
      outcome_type = "continuous",
      learners_trt = "SL.glm",
      learners_outcome = "SL.glm",
      folds = 1
    ))),
    error = function(e) NULL
  )
  skip_if(is.null(lmtp_res), "lmtp::lmtp_sdr() unavailable")

  est_lmtp <- tryCatch(lmtp_res$estimate@x, error = function(e) lmtp_res$theta)
  expect_lt(abs(est_causatr - est_lmtp), 0.3)
})

test_that("stochastic ICE: agrees with lmtp_sdr (longitudinal)", {
  skip_on_cran()
  skip_if_not_installed("lmtp")
  skip_if_not_installed("SuperLearner")

  set.seed(42)
  n <- 2000
  L0 <- rnorm(n)
  A0 <- rbinom(n, 1, plogis(0.5 * L0))
  L1 <- 0.5 * A0 + 0.5 * L0 + rnorm(n, 0, 0.5)
  A1 <- rbinom(n, 1, plogis(0.3 * L1))
  Y <- 1 + A0 + A1 + 0.5 * L0 + 0.5 * L1 + rnorm(n)

  d_long <- rbind(
    data.frame(id = seq_len(n), time = 0L, A = A0, L = NA_real_,
               L0 = L0, Y = NA_real_),
    data.frame(id = seq_len(n), time = 1L, A = A1, L = L1,
               L0 = L0, Y = Y)
  )

  sampler <- function(data, trt) {
    if ("L" %in% names(data) && !all(is.na(data$L))) {
      cov_col <- data$L
      cov_col[is.na(cov_col)] <- data$L0[is.na(cov_col)]
    } else {
      cov_col <- data$L0
    }
    rbinom(nrow(data), 1, plogis(0.2 + 0.3 * cov_col))
  }

  # causatr
  fit <- causat(
    d_long,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    estimator = "gcomp"
  )
  set.seed(123)
  res <- contrast(
    fit,
    interventions = list(g = stochastic(sampler, n_mc = 200L)),
    type = "difference",
    ci_method = "sandwich"
  )
  est_causatr <- res$estimates$estimate

  # lmtp wants wide format
  d_wide <- data.frame(
    L0 = L0,
    A_0 = A0,
    L_1 = L1,
    A_1 = A1,
    Y = Y
  )

  lmtp_shift <- function(data, trt) {
    if (trt == "A_0") {
      cov <- data[["L0"]]
    } else {
      cov <- data[["L_1"]]
    }
    rbinom(nrow(data), 1, plogis(0.2 + 0.3 * cov))
  }

  lmtp_res <- tryCatch(
    suppressWarnings(suppressMessages(lmtp::lmtp_sdr(
      data = d_wide,
      trt = c("A_0", "A_1"),
      outcome = "Y",
      baseline = "L0",
      time_vary = list(NULL, "L_1"),
      shift = lmtp_shift,
      outcome_type = "continuous",
      learners_trt = "SL.glm",
      learners_outcome = "SL.glm",
      folds = 1
    ))),
    error = function(e) NULL
  )
  skip_if(is.null(lmtp_res), "lmtp::lmtp_sdr() unavailable")

  est_lmtp <- tryCatch(lmtp_res$estimate@x, error = function(e) lmtp_res$theta)
  expect_lt(abs(est_causatr - est_lmtp), 0.5)
})
