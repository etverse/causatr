# Simulation-based tests for multivariate treatments with known truth.
#
# DGP (point treatment):
#   L ~ N(0, 1)
#   A1 | L ~ Bernoulli(expit(0.3 * L))
#   A2 | L ~ Bernoulli(expit(-0.3 * L))
#   Y | A1, A2, L ~ N(2 + 1.5*A1 + 1.0*A2 - 0.5*L, sd = 1)
#
# True potential outcome means (E[L] = 0 for large n):
#   E[Y(a1,a2)] = 2 + 1.5*a1 + 1.0*a2
#   E[Y(1,1)] = 4.5
#   E[Y(1,0)] = 3.5
#   E[Y(0,1)] = 3.0
#   E[Y(0,0)] = 2.0
#   ATE(both vs neither) = 2.5
#   ATE(A1 only) = 1.5
#   ATE(A2 only) = 1.0

simulate_mv_point <- function(n = 3000, seed = 42) {
  set.seed(seed)
  L <- rnorm(n)
  A1 <- rbinom(n, 1, plogis(0.3 * L))
  A2 <- rbinom(n, 1, plogis(-0.3 * L))
  Y <- 2 + 1.5 * A1 + 1.0 * A2 - 0.5 * L + rnorm(n)
  data.frame(L = L, A1 = A1, A2 = A2, Y = Y)
}

# DGP (point treatment, continuous A2):
#   L ~ N(0, 1)
#   A1 ~ Bernoulli(0.5)
#   A2 ~ N(5 + 0.5*L, sd = 1)
#   Y = 1 + 2*A1 + 0.3*A2 - 0.5*L + N(0,1)
#
# True E[Y(a1, a2)] = 1 + 2*a1 + 0.3*a2
# shift A2 by -2: E[Y(1, A2-2)] - E[Y(0, A2-2)] = 2 (A1 effect unchanged)
# shift A2 by -2 vs natural: E[Y(a1, A2-2)] - E[Y(a1, A2)] = 0.3*(-2) = -0.6

simulate_mv_mixed <- function(n = 3000, seed = 42) {
  set.seed(seed)
  L <- rnorm(n)
  A1 <- rbinom(n, 1, 0.5)
  A2 <- 5 + 0.5 * L + rnorm(n)
  Y <- 1 + 2 * A1 + 0.3 * A2 - 0.5 * L + rnorm(n)
  data.frame(L = L, A1 = A1, A2 = A2, Y = Y)
}

# DGP (point treatment, binary outcome):
#   L ~ N(0, 1)
#   A1 ~ Bernoulli(expit(0.3*L))
#   A2 ~ Bernoulli(expit(-0.3*L))
#   Y ~ Bernoulli(expit(-1 + A1 + 0.8*A2 + 0.5*L))
#
# True risks (computed via Monte Carlo with n = 1e6):
#   P[Y(1,1)=1] ≈ 0.622
#   P[Y(0,0)=1] ≈ 0.269
#   True RD(both vs neither) ≈ 0.353
#   True RR(both vs neither) ≈ 2.31

simulate_mv_binary <- function(n = 3000, seed = 42) {
  set.seed(seed)
  L <- rnorm(n)
  A1 <- rbinom(n, 1, plogis(0.3 * L))
  A2 <- rbinom(n, 1, plogis(-0.3 * L))
  Y <- rbinom(n, 1, plogis(-1 + A1 + 0.8 * A2 + 0.5 * L))
  data.frame(L = L, A1 = A1, A2 = A2, Y = Y)
}

# DGP (longitudinal multivariate):
#   L0 ~ N(0, 1)
#   A1_0, A2_0 ~ Bernoulli(0.5)
#   L1 ~ N(0.5*A1_0, 1)
#   A1_1 ~ Bernoulli(expit(0.3*L1))
#   A2_1 ~ Bernoulli(expit(0.3*L1))
#   Y = 50 + 2*A1_0 + 1*A2_0 + 2*A1_1 + 1*A2_1 - L1 + N(0, 3)
#
# True E[Y(1,1)] (both always = 1):
#   = 50 + 2 + 1 + 2 + 1 - E[L1|A1_0=1] = 56 - 0.5 = 55.5
# True E[Y(0,0)] (both always = 0):
#   = 50 + 0 + 0 + 0 + 0 - E[L1|A1_0=0] = 50 - 0 = 50
# True ATE ≈ 5.5

simulate_mv_long <- function(n = 2000, seed = 42) {
  set.seed(seed)
  L0 <- rnorm(n)
  A1_0 <- rbinom(n, 1, 0.5)
  A2_0 <- rbinom(n, 1, 0.5)
  L1 <- rnorm(n, mean = 0.5 * A1_0)
  A1_1 <- rbinom(n, 1, plogis(0.3 * L1))
  A2_1 <- rbinom(n, 1, plogis(0.3 * L1))
  Y <- 50 + 2 * A1_0 + 1 * A2_0 + 2 * A1_1 + 1 * A2_1 - L1 + rnorm(n, 0, 3)

  sim_id <- rep(seq_len(n), each = 2)
  sim_time <- rep(0:1, n)

  data.frame(
    id = sim_id,
    time = sim_time,
    A1 = c(rbind(A1_0, A1_1)),
    A2 = c(rbind(A2_0, A2_1)),
    L = c(rbind(L0, L1)),
    Y = c(rbind(NA_real_, Y))
  )
}

tol <- 0.3


test_that("mv point: static, all four combinations recover truth", {
  df <- simulate_mv_point()
  fit <- causat(df, "Y", c("A1", "A2"), ~L)

  result <- contrast(
    fit,
    interventions = list(
      both = list(A1 = static(1), A2 = static(1)),
      a1_only = list(A1 = static(1), A2 = static(0)),
      a2_only = list(A1 = static(0), A2 = static(1)),
      neither = list(A1 = static(0), A2 = static(0))
    ),
    reference = "neither"
  )

  ey <- setNames(result$estimates$estimate, result$estimates$intervention)
  expect_true(abs(ey["both"] - 4.5) < tol)
  expect_true(abs(ey["a1_only"] - 3.5) < tol)
  expect_true(abs(ey["a2_only"] - 3.0) < tol)
  expect_true(abs(ey["neither"] - 2.0) < tol)

  ate_both <- result$contrasts$estimate[
    result$contrasts$comparison == "both vs neither"
  ]
  ate_a1 <- result$contrasts$estimate[
    result$contrasts$comparison == "a1_only vs neither"
  ]
  ate_a2 <- result$contrasts$estimate[
    result$contrasts$comparison == "a2_only vs neither"
  ]
  expect_true(abs(ate_both - 2.5) < tol)
  expect_true(abs(ate_a1 - 1.5) < tol)
  expect_true(abs(ate_a2 - 1.0) < tol)
})

test_that("mv point: A1:A2 treatment-by-treatment interaction recovers truth", {
  # DGP with a multiplicative A1*A2 interaction on top of simulate_mv_point:
  #   Y = 2 + 1.5*A1 + 1.0*A2 + 0.8*A1*A2 - 0.5*L + N(0,1)
  # Truth:
  #   E[Y(0,0)] = 2.0, E[Y(1,0)] = 3.5, E[Y(0,1)] = 3.0, E[Y(1,1)] = 5.3
  #   ATE(both vs neither)   = 3.3
  #   ATE(a1_only vs neither)= 1.5
  #   ATE(a2_only vs neither)= 1.0
  #
  # Unlike `A:modifier`, an interaction between the two treatment columns
  # is legal in IPW/matching terms (no baseline moderator) AND in point
  # gcomp's outcome model — the whole A1:A2 surface is modeled. This test
  # pins that behavior so a future refactor of the multivariate path
  # cannot silently collapse the interaction.
  set.seed(42)
  n <- 3000
  L <- rnorm(n)
  A1 <- rbinom(n, 1, plogis(0.3 * L))
  A2 <- rbinom(n, 1, plogis(-0.3 * L))
  Y <- 2 + 1.5 * A1 + 1.0 * A2 + 0.8 * A1 * A2 - 0.5 * L + rnorm(n)
  df <- data.frame(L = L, A1 = A1, A2 = A2, Y = Y)

  fit <- causat(df, "Y", c("A1", "A2"), ~ L + A1:A2)

  result <- contrast(
    fit,
    interventions = list(
      both = list(A1 = static(1), A2 = static(1)),
      a1_only = list(A1 = static(1), A2 = static(0)),
      a2_only = list(A1 = static(0), A2 = static(1)),
      neither = list(A1 = static(0), A2 = static(0))
    ),
    reference = "neither"
  )

  ey <- setNames(result$estimates$estimate, result$estimates$intervention)
  expect_true(abs(ey["both"] - 5.3) < tol)
  expect_true(abs(ey["a1_only"] - 3.5) < tol)
  expect_true(abs(ey["a2_only"] - 3.0) < tol)
  expect_true(abs(ey["neither"] - 2.0) < tol)

  ates <- setNames(
    result$contrasts$estimate,
    result$contrasts$comparison
  )
  expect_true(abs(ates["both vs neither"] - 3.3) < tol)
  expect_true(abs(ates["a1_only vs neither"] - 1.5) < tol)
  expect_true(abs(ates["a2_only vs neither"] - 1.0) < tol)
})

test_that("mv point: shift intervention on continuous A2", {
  df <- simulate_mv_mixed()
  fit <- causat(df, "Y", c("A1", "A2"), ~L)

  result <- contrast(
    fit,
    interventions = list(
      treat_shift = list(A1 = static(1), A2 = shift(-2)),
      control_shift = list(A1 = static(0), A2 = shift(-2))
    ),
    reference = "control_shift"
  )

  ate <- result$contrasts$estimate[1]
  expect_true(abs(ate - 2.0) < tol)
})

test_that("mv point: scale intervention on continuous A2", {
  df <- simulate_mv_mixed()
  fit <- causat(df, "Y", c("A1", "A2"), ~L)

  result <- contrast(
    fit,
    interventions = list(
      full = list(A1 = static(1), A2 = scale_by(1)),
      scaled = list(A1 = static(1), A2 = scale_by(0.5))
    ),
    reference = "full"
  )

  expect_s3_class(result, "causatr_result")
  diff <- result$contrasts$estimate[1]
  expect_true(diff < 0)
})

test_that("mv point: threshold intervention on continuous A2", {
  df <- simulate_mv_mixed()
  fit <- causat(df, "Y", c("A1", "A2"), ~L)

  result <- contrast(
    fit,
    interventions = list(
      uncapped = list(A1 = static(1), A2 = threshold(-Inf, Inf)),
      capped = list(A1 = static(1), A2 = threshold(-Inf, 4))
    ),
    reference = "uncapped"
  )

  expect_s3_class(result, "causatr_result")
  diff <- result$contrasts$estimate[1]
  expect_true(diff < 0)
})

test_that("mv point: dynamic intervention", {
  df <- simulate_mv_point()
  fit <- causat(df, "Y", c("A1", "A2"), ~L)

  result <- contrast(
    fit,
    interventions = list(
      adaptive = list(
        A1 = dynamic(\(data, trt) ifelse(data$L > 0, 1, 0)),
        A2 = static(1)
      ),
      always = list(A1 = static(1), A2 = static(1))
    )
  )

  expect_s3_class(result, "causatr_result")
  ey_adaptive <- result$estimates$estimate[
    result$estimates$intervention == "adaptive"
  ]
  ey_always <- result$estimates$estimate[
    result$estimates$intervention == "always"
  ]
  expect_true(ey_adaptive < ey_always)
})

test_that("mv point: binary outcome, risk difference", {
  df <- simulate_mv_binary(n = 5000)
  fit <- causat(df, "Y", c("A1", "A2"), ~L, family = "binomial")

  result <- contrast(
    fit,
    interventions = list(
      both = list(A1 = static(1), A2 = static(1)),
      neither = list(A1 = static(0), A2 = static(0))
    ),
    reference = "neither",
    type = "difference"
  )

  rd <- result$contrasts$estimate[1]
  expect_true(abs(rd - 0.353) < 0.15)
})

test_that("mv point: binary outcome, risk ratio", {
  df <- simulate_mv_binary(n = 5000)
  fit <- causat(df, "Y", c("A1", "A2"), ~L, family = "binomial")

  result <- contrast(
    fit,
    interventions = list(
      both = list(A1 = static(1), A2 = static(1)),
      neither = list(A1 = static(0), A2 = static(0))
    ),
    reference = "neither",
    type = "ratio"
  )

  rr <- result$contrasts$estimate[1]
  expect_true(rr > 1.5 && rr < 3.5)
})

test_that("mv point: binary outcome, odds ratio", {
  df <- simulate_mv_binary(n = 5000)
  fit <- causat(df, "Y", c("A1", "A2"), ~L, family = "binomial")

  result <- contrast(
    fit,
    interventions = list(
      both = list(A1 = static(1), A2 = static(1)),
      neither = list(A1 = static(0), A2 = static(0))
    ),
    reference = "neither",
    type = "or"
  )

  or_est <- result$contrasts$estimate[1]
  expect_true(or_est > 2.0)
})

test_that("mv point: subset estimand", {
  df <- simulate_mv_point()
  df$sex <- rep(c(0, 1), length.out = nrow(df))
  fit <- causat(df, "Y", c("A1", "A2"), ~ L + sex)

  result <- contrast(
    fit,
    interventions = list(
      both = list(A1 = static(1), A2 = static(1)),
      neither = list(A1 = static(0), A2 = static(0))
    ),
    subset = quote(sex == 1),
    reference = "neither"
  )

  expect_equal(result$estimand, "subset")
  ate <- result$contrasts$estimate[1]
  expect_true(abs(ate - 2.5) < tol)
})

test_that("mv point: ATT/ATC rejected", {
  df <- simulate_mv_point()
  expect_snapshot(
    error = TRUE,
    causat(df, "Y", c("A1", "A2"), ~L, estimand = "ATT")
  )
})

test_that("mv point: matching rejected", {
  df <- simulate_mv_point()
  expect_snapshot(
    error = TRUE,
    causat(df, "Y", c("A1", "A2"), ~L, estimator = "matching")
  )
})

test_that("mv point: sandwich SE is positive and finite", {
  df <- simulate_mv_point()
  fit <- causat(df, "Y", c("A1", "A2"), ~L)

  result <- contrast(
    fit,
    interventions = list(
      both = list(A1 = static(1), A2 = static(1)),
      neither = list(A1 = static(0), A2 = static(0))
    ),
    ci_method = "sandwich"
  )

  expect_true(all(result$estimates$se > 0))
  expect_true(all(is.finite(result$estimates$se)))
  expect_true(all(result$contrasts$se > 0))
  expect_true(all(is.finite(result$contrasts$se)))
})

test_that("mv point: bootstrap SE agrees with sandwich", {
  df <- simulate_mv_point(n = 500)
  fit <- causat(df, "Y", c("A1", "A2"), ~L)

  ivs <- list(
    both = list(A1 = static(1), A2 = static(1)),
    neither = list(A1 = static(0), A2 = static(0))
  )

  res_sw <- contrast(fit, ivs, ci_method = "sandwich")
  res_bs <- contrast(fit, ivs, ci_method = "bootstrap", n_boot = 50L)

  ratio <- res_bs$contrasts$se[1] / res_sw$contrasts$se[1]
  expect_true(ratio > 0.3 && ratio < 3.0)
})

test_that("mv point: three+ interventions with reference", {
  df <- simulate_mv_point()
  fit <- causat(df, "Y", c("A1", "A2"), ~L)

  result <- contrast(
    fit,
    interventions = list(
      both = list(A1 = static(1), A2 = static(1)),
      a1_only = list(A1 = static(1), A2 = static(0)),
      a2_only = list(A1 = static(0), A2 = static(1)),
      neither = list(A1 = static(0), A2 = static(0))
    ),
    reference = "neither"
  )

  expect_equal(nrow(result$contrasts), 3L)
  expect_true(all(grepl("vs neither", result$contrasts$comparison)))
})

test_that("mv point: by argument produces subgroup estimates", {
  df <- simulate_mv_point()
  df$sex <- rep(c(0, 1), length.out = nrow(df))
  fit <- causat(df, "Y", c("A1", "A2"), ~ L + sex)

  result <- contrast(
    fit,
    interventions = list(
      both = list(A1 = static(1), A2 = static(1)),
      neither = list(A1 = static(0), A2 = static(0))
    ),
    by = "sex"
  )

  expect_true("by" %in% names(result$estimates))
  expect_equal(nrow(result$estimates), 4L)
  expect_equal(nrow(result$contrasts), 2L)
})


test_that("mv longitudinal ICE: static recovers truth", {
  long <- simulate_mv_long()
  fit <- causat(
    long,
    "Y",
    c("A1", "A2"),
    confounders = ~1,
    confounders_tv = ~L,
    id = "id",
    time = "time"
  )

  expect_equal(fit$type, "longitudinal")

  result <- contrast(
    fit,
    interventions = list(
      both = list(A1 = static(1), A2 = static(1)),
      neither = list(A1 = static(0), A2 = static(0))
    ),
    reference = "neither",
    ci_method = "sandwich"
  )

  ey_both <- result$estimates$estimate[result$estimates$intervention == "both"]
  ey_neither <- result$estimates$estimate[
    result$estimates$intervention == "neither"
  ]
  ate <- result$contrasts$estimate[1]

  expect_true(abs(ey_both - 55.5) < 2.0)
  expect_true(abs(ey_neither - 50.0) < 2.0)
  expect_true(abs(ate - 5.5) < 2.0)
})

test_that("mv longitudinal ICE: dynamic intervention", {
  long <- simulate_mv_long()
  fit <- causat(
    long,
    "Y",
    c("A1", "A2"),
    confounders = ~1,
    confounders_tv = ~L,
    id = "id",
    time = "time"
  )

  result <- contrast(
    fit,
    interventions = list(
      adaptive = list(
        A1 = dynamic(\(data, trt) ifelse(!is.na(data$L) & data$L > 0, 1L, 0L)),
        A2 = static(1)
      ),
      both = list(A1 = static(1), A2 = static(1)),
      neither = list(A1 = static(0), A2 = static(0))
    ),
    reference = "neither",
    ci_method = "sandwich"
  )

  expect_equal(nrow(result$estimates), 3L)
  expect_equal(nrow(result$contrasts), 2L)

  ey_adaptive <- result$estimates$estimate[
    result$estimates$intervention == "adaptive"
  ]
  ey_both <- result$estimates$estimate[result$estimates$intervention == "both"]
  ey_neither <- result$estimates$estimate[
    result$estimates$intervention == "neither"
  ]
  expect_true(ey_adaptive > ey_neither)
  expect_true(ey_adaptive < ey_both)
})

test_that("mv longitudinal ICE: sandwich SE is positive", {
  long <- simulate_mv_long()
  fit <- causat(
    long,
    "Y",
    c("A1", "A2"),
    confounders = ~1,
    confounders_tv = ~L,
    id = "id",
    time = "time"
  )

  result <- contrast(
    fit,
    interventions = list(
      both = list(A1 = static(1), A2 = static(1)),
      neither = list(A1 = static(0), A2 = static(0))
    ),
    ci_method = "sandwich"
  )

  expect_true(all(result$estimates$se > 0))
  expect_true(all(is.finite(result$estimates$se)))
})

test_that("mv longitudinal ICE: bootstrap works", {
  long <- simulate_mv_long(n = 200)
  fit <- causat(
    long,
    "Y",
    c("A1", "A2"),
    confounders = ~1,
    confounders_tv = ~L,
    id = "id",
    time = "time"
  )

  result <- contrast(
    fit,
    interventions = list(
      both = list(A1 = static(1), A2 = static(1)),
      neither = list(A1 = static(0), A2 = static(0))
    ),
    ci_method = "bootstrap",
    n_boot = 10L
  )

  expect_equal(result$ci_method, "bootstrap")
  expect_true(all(result$estimates$se > 0))
})
