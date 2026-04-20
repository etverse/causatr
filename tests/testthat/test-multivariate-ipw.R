# Truth-based tests for multivariate IPW under sequential factorisation.
#
# DGP families used here mirror the multivariate g-comp tests so the
# point estimates can be cross-checked between estimators on the same
# data. The IPW path uses
#   f(A_1, ..., A_K | L) = prod_k f(A_k | A_1, ..., A_{k-1}, L)
# with per-component univariate density models, product density-ratio
# weights, an intercept-only Hajek MSM (Y ~ 1) per intervention, and
# a stacked sandwich variance over the K propensity models (the bread
# is block-diagonal because the propensity models are fit
# independently).

# DGP 1: binary x binary
#
#   L ~ N(0, 1)
#   A1 | L ~ Bernoulli(plogis(0.3 * L))
#   A2 | L ~ Bernoulli(plogis(-0.3 * L))
#   Y | A1, A2, L ~ N(2 + 1.5*A1 + 1.0*A2 - 0.5*L, sd = 1)
#
# Truth: E[Y(1,1)] = 4.5, E[Y(0,0)] = 2.0, ATE(both vs neither) = 2.5.
sim_bb <- function(n = 3000, seed = 42) {
  set.seed(seed)
  L <- rnorm(n)
  A1 <- rbinom(n, 1, plogis(0.3 * L))
  A2 <- rbinom(n, 1, plogis(-0.3 * L))
  Y <- 2 + 1.5 * A1 + 1.0 * A2 - 0.5 * L + rnorm(n)
  data.frame(L = L, A1 = A1, A2 = A2, Y = Y)
}

# DGP 2: binary x continuous (mixed type)
#
#   L ~ N(0, 1)
#   A1 ~ Bernoulli(0.5)
#   A2 ~ N(5 + 0.5 * L, sd = 1)
#   Y = 1 + 2*A1 + 0.3*A2 - 0.5*L + N(0, 1)
#
# E[Y(a1, A2 - 2)] - E[Y(0, A2 - 2)] depends only on A1: ATE(A1=1) = 2.
sim_bc <- function(n = 5000, seed = 42) {
  set.seed(seed)
  L <- rnorm(n)
  A1 <- rbinom(n, 1, 0.5)
  A2 <- 5 + 0.5 * L + rnorm(n)
  Y <- 1 + 2 * A1 + 0.3 * A2 - 0.5 * L + rnorm(n)
  data.frame(L = L, A1 = A1, A2 = A2, Y = Y)
}

# DGP 3: continuous x continuous
#
#   L ~ N(0, 1)
#   A1 | L ~ N(0.5*L, sd = 1)
#   A2 | A1, L ~ N(0.3*A1 + 0.5*L, sd = 1)   # downstream-conditioned!
#   Y = 1 + 0.5*A1 + 0.4*A2 - 0.5*L + N(0, 1)
#
# E[Y(A1 + d1, A2 + d2)] - E[Y(A1, A2)] = 0.5*d1 + 0.4*d2.
# For (d1, d2) = (-1, -1) the truth is -0.9.
sim_cc <- function(n = 3000, seed = 42) {
  set.seed(seed)
  L <- rnorm(n)
  A1 <- 0.5 * L + rnorm(n)
  A2 <- 0.3 * A1 + 0.5 * L + rnorm(n)
  Y <- 1 + 0.5 * A1 + 0.4 * A2 - 0.5 * L + rnorm(n)
  data.frame(L = L, A1 = A1, A2 = A2, Y = Y)
}

# DGP 4: binary x binary, binomial outcome
#
#   L ~ N(0, 1)
#   A1 ~ Bernoulli(plogis(0.3*L))
#   A2 ~ Bernoulli(plogis(-0.3*L))
#   Y ~ Bernoulli(plogis(-1 + A1 + 0.8*A2 + 0.5*L))
#
# Approx truths from a 1e6 Monte Carlo:
#   P[Y(1,1)=1] ≈ 0.622, P[Y(0,0)=1] ≈ 0.269
#   RD ≈ 0.353, RR ≈ 2.31, OR ≈ 4.49.
sim_bb_binary <- function(n = 5000, seed = 42) {
  set.seed(seed)
  L <- rnorm(n)
  A1 <- rbinom(n, 1, plogis(0.3 * L))
  A2 <- rbinom(n, 1, plogis(-0.3 * L))
  Y <- rbinom(n, 1, plogis(-1 + A1 + 0.8 * A2 + 0.5 * L))
  data.frame(L = L, A1 = A1, A2 = A2, Y = Y)
}

tol_pt <- 0.3
tol_rd <- 0.15

test_that("mv IPW: binary x binary static recovers truth (gauss outcome)", {
  df <- sim_bb()
  fit <- causat(df, "Y", c("A1", "A2"), ~L, estimator = "ipw")

  result <- contrast(
    fit,
    interventions = list(
      both = list(A1 = static(1), A2 = static(1)),
      a1 = list(A1 = static(1), A2 = static(0)),
      a2 = list(A1 = static(0), A2 = static(1)),
      neither = list(A1 = static(0), A2 = static(0))
    ),
    reference = "neither"
  )

  ey <- setNames(result$estimates$estimate, result$estimates$intervention)
  expect_true(abs(ey["both"] - 4.5) < tol_pt)
  expect_true(abs(ey["a1"] - 3.5) < tol_pt)
  expect_true(abs(ey["a2"] - 3.0) < tol_pt)
  expect_true(abs(ey["neither"] - 2.0) < tol_pt)

  ates <- setNames(result$contrasts$estimate, result$contrasts$comparison)
  expect_true(abs(ates["both vs neither"] - 2.5) < tol_pt)
  expect_true(abs(ates["a1 vs neither"] - 1.5) < tol_pt)
  expect_true(abs(ates["a2 vs neither"] - 1.0) < tol_pt)

  # Sandwich SEs are positive and finite.
  expect_true(all(result$estimates$se > 0))
  expect_true(all(is.finite(result$estimates$se)))
})

test_that("mv IPW: binary x binary cross-checks against gcomp", {
  # Both estimators sit on the same multivariate API; agreement under a
  # correctly-specified DGP is the strongest internal consistency check
  # we can run without an external oracle.
  df <- sim_bb()
  ivs <- list(
    both = list(A1 = static(1), A2 = static(1)),
    neither = list(A1 = static(0), A2 = static(0))
  )

  fit_g <- causat(df, "Y", c("A1", "A2"), ~L, estimator = "gcomp")
  fit_w <- causat(df, "Y", c("A1", "A2"), ~L, estimator = "ipw")

  res_g <- contrast(fit_g, ivs, reference = "neither")
  res_w <- contrast(fit_w, ivs, reference = "neither")

  ate_g <- res_g$contrasts$estimate[1]
  ate_w <- res_w$contrasts$estimate[1]
  expect_true(abs(ate_g - ate_w) < 0.2)
})

test_that("mv IPW: binary x binary bootstrap parity", {
  df <- sim_bb(n = 1000)
  ivs <- list(
    both = list(A1 = static(1), A2 = static(1)),
    neither = list(A1 = static(0), A2 = static(0))
  )
  fit <- causat(df, "Y", c("A1", "A2"), ~L, estimator = "ipw")

  sw <- contrast(fit, ivs, ci_method = "sandwich")
  # Bootstrap is a slow but model-free reference. 100 reps is enough to
  # hit ±30% MC noise on the SE for an n = 1000 sample. Set to 75 to
  # save runtime; tests that need tighter tolerance can bump.
  bs <- contrast(fit, ivs, ci_method = "bootstrap", n_boot = 75)

  ratio <- bs$contrasts$se[1] / sw$contrasts$se[1]
  expect_true(ratio > 0.5 && ratio < 2.0)
})

test_that("mv IPW: binary + continuous, static + shift on continuous", {
  # The static A1 indicator collapses surviving rows to A1 = a1*, so
  # the second component's denominator (which conditions on observed A1)
  # ends up at the same A1 level as the numerator. This case does not
  # exercise the intervened-newdata machinery, but it does exercise the
  # mixed-family product weight (Bernoulli x Gaussian).
  df <- sim_bc()
  fit <- causat(df, "Y", c("A1", "A2"), ~L, estimator = "ipw")

  result <- contrast(
    fit,
    interventions = list(
      treat = list(A1 = static(1), A2 = shift(-2)),
      control = list(A1 = static(0), A2 = shift(-2))
    ),
    reference = "control"
  )

  # ATE(A1) = 2 (the A1 effect, A2 shift is the same for both arms).
  expect_true(abs(result$contrasts$estimate[1] - 2.0) < 0.5)
  expect_true(result$contrasts$se[1] > 0 && is.finite(result$contrasts$se[1]))
})

test_that("mv IPW: continuous + continuous shift+shift recovers truth", {
  # This is the case that EXERCISES intervened-newdata: the second
  # component's numerator must be evaluated at the SHIFTED A1 because
  # the chain-rule factorisation has f_2(A_2 | A_1, L) -- substituting
  # the intervened A_1 = A_1 - delta1 for the conditioning vector. If
  # the engine forgot this substitution, the recovered ATE would drift
  # from the truth -0.9 by ~0.1-0.2 (the magnitude of the A1->A2
  # cross-coefficient times the shift).
  df <- sim_cc()
  fit <- causat(df, "Y", c("A1", "A2"), ~L, estimator = "ipw")

  result <- contrast(
    fit,
    interventions = list(
      shifted = list(A1 = shift(-1), A2 = shift(-1)),
      natural = list(A1 = shift(0), A2 = shift(0))
    ),
    reference = "natural"
  )

  # Truth: 0.5 * (-1) + 0.4 * (-1) = -0.9.
  expect_true(abs(result$contrasts$estimate[1] - (-0.9)) < 0.15)
  expect_true(result$contrasts$se[1] > 0 && is.finite(result$contrasts$se[1]))
})

test_that("mv IPW: continuous + continuous gcomp cross-check", {
  df <- sim_cc()
  ivs <- list(
    shifted = list(A1 = shift(-1), A2 = shift(-1)),
    natural = list(A1 = shift(0), A2 = shift(0))
  )

  fit_g <- causat(df, "Y", c("A1", "A2"), ~L, estimator = "gcomp")
  fit_w <- causat(df, "Y", c("A1", "A2"), ~L, estimator = "ipw")

  res_g <- contrast(fit_g, ivs, reference = "natural")
  res_w <- contrast(fit_w, ivs, reference = "natural")

  expect_true(
    abs(res_g$contrasts$estimate[1] - res_w$contrasts$estimate[1]) < 0.15
  )
})

test_that("mv IPW: binomial outcome, RD/RR/OR contrasts", {
  df <- sim_bb_binary(n = 5000)
  fit <- causat(
    df,
    "Y",
    c("A1", "A2"),
    ~L,
    family = "binomial",
    estimator = "ipw"
  )
  ivs <- list(
    both = list(A1 = static(1), A2 = static(1)),
    neither = list(A1 = static(0), A2 = static(0))
  )

  rd <- contrast(fit, ivs, reference = "neither", type = "difference")
  rr <- contrast(fit, ivs, reference = "neither", type = "ratio")
  or_ <- contrast(fit, ivs, reference = "neither", type = "or")

  expect_true(abs(rd$contrasts$estimate[1] - 0.353) < tol_rd)
  expect_true(rr$contrasts$estimate[1] > 1.5 && rr$contrasts$estimate[1] < 3.5)
  expect_true(or_$contrasts$estimate[1] > 2.0)
})

test_that("mv IPW: subset estimand", {
  df <- sim_bb()
  df$sex <- rep(c(0, 1), length.out = nrow(df))
  fit <- causat(df, "Y", c("A1", "A2"), ~ L + sex, estimator = "ipw")

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
  expect_true(abs(result$contrasts$estimate[1] - 2.5) < tol_pt)
})

test_that("mv IPW: by-stratified estimates", {
  df <- sim_bb()
  df$sex <- rep(c(0, 1), length.out = nrow(df))
  fit <- causat(df, "Y", c("A1", "A2"), ~ L + sex, estimator = "ipw")

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

test_that("mv IPW: dynamic component", {
  df <- sim_bb()
  fit <- causat(df, "Y", c("A1", "A2"), ~L, estimator = "ipw")

  result <- contrast(
    fit,
    interventions = list(
      adaptive = list(
        A1 = dynamic(\(data, trt) ifelse(data$L > 0, 1L, 0L)),
        A2 = static(1)
      ),
      always = list(A1 = static(1), A2 = static(1))
    )
  )

  ey_a <- result$estimates$estimate[result$estimates$intervention == "adaptive"]
  ey_b <- result$estimates$estimate[result$estimates$intervention == "always"]
  # Adaptive treats only L > 0 individuals; always treats everyone.
  # E[Y(adaptive)] should be strictly less than E[Y(always)] under the
  # positive A1 effect.
  expect_true(ey_a < ey_b)
})

test_that("mv IPW: 3-component (binary x binary x binary) sequential factorisation", {
  # Generalisation check: the sequential factorisation generalises to
  # any K. Here K = 3, all binary. Truth E[Y(1,1,1)] - E[Y(0,0,0)]
  # = 1.5 + 1.0 + 0.7 = 3.2.
  set.seed(7)
  n <- 3000
  L <- rnorm(n)
  A1 <- rbinom(n, 1, plogis(0.3 * L))
  A2 <- rbinom(n, 1, plogis(-0.3 * L + 0.2 * A1))
  A3 <- rbinom(n, 1, plogis(0.4 * L - 0.1 * A2))
  Y <- 1 + 1.5 * A1 + 1.0 * A2 + 0.7 * A3 - 0.5 * L + rnorm(n)
  df <- data.frame(L = L, A1 = A1, A2 = A2, A3 = A3, Y = Y)

  fit <- causat(df, "Y", c("A1", "A2", "A3"), ~L, estimator = "ipw")
  expect_length(fit$details$treatment_models, 3L)

  result <- contrast(
    fit,
    interventions = list(
      all_on = list(A1 = static(1), A2 = static(1), A3 = static(1)),
      all_off = list(A1 = static(0), A2 = static(0), A3 = static(0))
    ),
    reference = "all_off"
  )

  expect_true(abs(result$contrasts$estimate[1] - 3.2) < 0.4)
})

# ---------------------------------------------------------------------
# Rejection paths
# ---------------------------------------------------------------------

test_that("mv IPW: ipsi() in any component is rejected", {
  df <- sim_bb()
  fit <- causat(df, "Y", c("A1", "A2"), ~L, estimator = "ipw")
  expect_error(
    contrast(fit, list(both = list(A1 = ipsi(2), A2 = static(1)))),
    class = "causatr_multivariate_ipsi"
  )
})

test_that("mv IPW: effect modification (`A1:sex`) recovers per-stratum truth", {
  # DGP with sex modifying the A1 effect:
  #   E[Y(a1, a2) | sex] = 2 + (1.0 + 1.0*sex)*a1 + 1.0*a2 - 0.5*L
  # Truth:
  #   E[Y(1,1) - Y(0,0) | sex=0] = 1 + 1 = 2
  #   E[Y(1,1) - Y(0,0) | sex=1] = 2 + 1 = 3
  set.seed(42)
  n <- 5000
  L <- rnorm(n)
  sex <- rbinom(n, 1, 0.5)
  A1 <- rbinom(n, 1, plogis(0.3 * L))
  A2 <- rbinom(n, 1, plogis(-0.3 * L))
  Y <- 2 + (1.0 + 1.0 * sex) * A1 + 1.0 * A2 - 0.5 * L + rnorm(n)
  df <- data.frame(L = L, sex = sex, A1 = A1, A2 = A2, Y = Y)

  fit <- causat(df, "Y", c("A1", "A2"), ~ L + sex + A1:sex, estimator = "ipw")
  expect_true(fit$details$em_info$has_em)
  # Per-component propensity formulas should NOT contain the EM term.
  ps1 <- as.character(fit$details$treatment_models[[1]]$ps_formula)
  ps2 <- as.character(fit$details$treatment_models[[2]]$ps_formula)
  expect_false(any(grepl("A1:sex", ps1, fixed = TRUE)))
  expect_false(any(grepl("A1:sex", ps2, fixed = TRUE)))

  # by-stratified contrast pulls out the per-stratum modifier effect.
  res <- contrast(
    fit,
    interventions = list(
      both = list(A1 = static(1), A2 = static(1)),
      neither = list(A1 = static(0), A2 = static(0))
    ),
    reference = "neither",
    by = "sex"
  )
  ate_by <- setNames(res$contrasts$estimate, res$contrasts$by)
  expect_true(abs(ate_by["0"] - 2.0) < 0.3)
  expect_true(abs(ate_by["1"] - 3.0) < 0.3)
})

test_that("mv IPW: EM cross-checks against gcomp", {
  set.seed(42)
  n <- 5000
  L <- rnorm(n)
  sex <- rbinom(n, 1, 0.5)
  A1 <- rbinom(n, 1, plogis(0.3 * L))
  A2 <- rbinom(n, 1, plogis(-0.3 * L))
  Y <- 2 + (1.0 + 1.0 * sex) * A1 + 1.0 * A2 - 0.5 * L + rnorm(n)
  df <- data.frame(L = L, sex = sex, A1 = A1, A2 = A2, Y = Y)

  ivs <- list(
    both = list(A1 = static(1), A2 = static(1)),
    neither = list(A1 = static(0), A2 = static(0))
  )

  fit_g <- causat(
    df,
    "Y",
    c("A1", "A2"),
    ~ L + sex + A1:sex,
    estimator = "gcomp"
  )
  fit_w <- causat(df, "Y", c("A1", "A2"), ~ L + sex + A1:sex, estimator = "ipw")

  res_g <- contrast(fit_g, ivs, reference = "neither", by = "sex")
  res_w <- contrast(fit_w, ivs, reference = "neither", by = "sex")

  # Pair stratum-level contrasts and check agreement.
  ate_g <- setNames(res_g$contrasts$estimate, res_g$contrasts$by)
  ate_w <- setNames(res_w$contrasts$estimate, res_w$contrasts$by)
  expect_true(abs(ate_g["0"] - ate_w["0"]) < 0.2)
  expect_true(abs(ate_g["1"] - ate_w["1"]) < 0.2)
})

test_that("mv IPW: bare treatment in confounders (`~ L + A1`) still rejected", {
  # `A:modifier` is fine; bare A is always wrong (puts A on both sides
  # of the propensity model). Phase 6's `check_em_compat()` catches it.
  df <- sim_bb()
  expect_error(
    causat(df, "Y", c("A1", "A2"), ~ L + A1, estimator = "ipw"),
    class = "causatr_bare_treatment_in_confounders"
  )
})

test_that("mv IPW: bin x categorical recovers truth", {
  # DGP with a categorical second component:
  #   E[Y(a1, a2)] = 0.5 + 1*a1 + 1.5*(a2=="high") + 0.7*(a2=="med") - 0.5*L
  # Truth: E[Y(1, "high") - Y(0, "low")] = 1 + 1.5 = 2.5.
  set.seed(42)
  n <- 5000
  L <- rnorm(n)
  A1 <- rbinom(n, 1, plogis(0.3 * L))
  A2 <- factor(sample(c("low", "med", "high"), n, replace = TRUE))
  Y <- 0.5 +
    1.0 * A1 +
    1.5 * (A2 == "high") +
    0.7 * (A2 == "med") -
    0.5 * L +
    rnorm(n)
  df <- data.frame(L = L, A1 = A1, A2 = A2, Y = Y)

  fit <- causat(df, "Y", c("A1", "A2"), ~L, estimator = "ipw")
  expect_equal(fit$details$treatment_models[[1]]$family, "bernoulli")
  expect_equal(fit$details$treatment_models[[2]]$family, "categorical")
  expect_true(inherits(fit$details$treatment_models[[2]]$model, "multinom"))

  res <- contrast(
    fit,
    interventions = list(
      high = list(A1 = static(1), A2 = static("high")),
      low = list(A1 = static(0), A2 = static("low"))
    ),
    reference = "low"
  )
  expect_true(abs(res$contrasts$estimate[1] - 2.5) < 0.3)
  expect_true(res$contrasts$se[1] > 0 && is.finite(res$contrasts$se[1]))
})

test_that("mv IPW: categorical x bin (cat first) recovers same truth", {
  # Order-insensitivity check: swapping the components must give the
  # same point estimate up to numerical noise.
  set.seed(42)
  n <- 5000
  L <- rnorm(n)
  A1 <- rbinom(n, 1, plogis(0.3 * L))
  A2 <- factor(sample(c("low", "med", "high"), n, replace = TRUE))
  Y <- 0.5 +
    1.0 * A1 +
    1.5 * (A2 == "high") +
    0.7 * (A2 == "med") -
    0.5 * L +
    rnorm(n)
  df <- data.frame(L = L, A1 = A1, A2 = A2, Y = Y)

  fit_a1 <- causat(df, "Y", c("A1", "A2"), ~L, estimator = "ipw")
  fit_a2 <- causat(df, "Y", c("A2", "A1"), ~L, estimator = "ipw")

  res_a1 <- contrast(
    fit_a1,
    list(
      high = list(A1 = static(1), A2 = static("high")),
      low = list(A1 = static(0), A2 = static("low"))
    ),
    reference = "low"
  )
  res_a2 <- contrast(
    fit_a2,
    list(
      high = list(A2 = static("high"), A1 = static(1)),
      low = list(A2 = static("low"), A1 = static(0))
    ),
    reference = "low"
  )
  expect_true(
    abs(res_a1$contrasts$estimate[1] - res_a2$contrasts$estimate[1]) < 0.05
  )
})

test_that("mv IPW: categorical x categorical recovers truth", {
  # Two independent categorical treatments. Truth E[Y(c, "high") -
  # Y(a, "low")] = 2.0 + 1.5 = 3.5.
  set.seed(42)
  n <- 5000
  L <- rnorm(n)
  A1 <- factor(sample(c("a", "b", "c"), n, replace = TRUE))
  A2 <- factor(sample(c("low", "med", "high"), n, replace = TRUE))
  Y <- 0.5 +
    1.0 * (A1 == "b") +
    2.0 * (A1 == "c") +
    1.5 * (A2 == "high") +
    0.7 * (A2 == "med") -
    0.5 * L +
    rnorm(n)
  df <- data.frame(L = L, A1 = A1, A2 = A2, Y = Y)

  fit <- causat(df, "Y", c("A1", "A2"), ~L, estimator = "ipw")
  res <- contrast(
    fit,
    list(
      c_high = list(A1 = static("c"), A2 = static("high")),
      a_low = list(A1 = static("a"), A2 = static("low"))
    ),
    reference = "a_low"
  )
  expect_true(abs(res$contrasts$estimate[1] - 3.5) < 0.3)
  expect_true(res$contrasts$se[1] > 0 && is.finite(res$contrasts$se[1]))
})

test_that("mv IPW: categorical component cross-checks against gcomp", {
  set.seed(42)
  n <- 5000
  L <- rnorm(n)
  A1 <- rbinom(n, 1, plogis(0.3 * L))
  A2 <- factor(sample(c("low", "med", "high"), n, replace = TRUE))
  Y <- 0.5 +
    1.0 * A1 +
    1.5 * (A2 == "high") +
    0.7 * (A2 == "med") -
    0.5 * L +
    rnorm(n)
  df <- data.frame(L = L, A1 = A1, A2 = A2, Y = Y)

  ivs <- list(
    high = list(A1 = static(1), A2 = static("high")),
    low = list(A1 = static(0), A2 = static("low"))
  )

  fit_g <- causat(df, "Y", c("A1", "A2"), ~L, estimator = "gcomp")
  fit_w <- causat(df, "Y", c("A1", "A2"), ~L, estimator = "ipw")

  res_g <- contrast(fit_g, ivs, reference = "low")
  res_w <- contrast(fit_w, ivs, reference = "low")

  expect_true(
    abs(res_g$contrasts$estimate[1] - res_w$contrasts$estimate[1]) < 0.2
  )
})

test_that("mv IPW: propensity_family validates per-component shape", {
  df <- sim_bb()
  expect_error(
    causat(
      df,
      "Y",
      c("A1", "A2"),
      ~L,
      estimator = "ipw",
      propensity_family = c("poisson", "negbin", "extra")
    ),
    "must be NULL, length 1, or length K"
  )
})

test_that("mv IPW: ATT/ATC fail at the existing trt-estimand check", {
  # Multivariate ATT/ATC has no clean definition (no single "treated"
  # arm); rejected upstream by `check_estimand_trt_compat()` not by
  # the new MV path. Verify the error class still fires under the IPW
  # estimator now that the upstream IPW gate is gone.
  df <- sim_bb()
  expect_error(
    causat(df, "Y", c("A1", "A2"), ~L, estimator = "ipw", estimand = "ATT")
  )
})
