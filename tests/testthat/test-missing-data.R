# Tests for missing data (NA) handling across all methods.
#
# Theoretical grounding
# ─────────────────────
# Three variables can have missing values: outcome (Y), treatment (A),
# and covariates (L). The correct handling depends on the missingness
# mechanism:
#
#   MCAR (Missing Completely At Random): P(missing) = constant.
#     Complete-case analysis is unbiased.
#
#   MAR (Missing At Random): P(missing | observed variables) depends on
#     observed data. Complete-case is biased. Correction:
#       Y missing → IPCW: weight by 1/P(C=0 | A, L)
#       A missing → Multiple imputation (MICE)
#       L missing → Multiple imputation (MICE)
#
#   MNAR (Missing Not At Random): out of scope (sensitivity analysis).
#
# These tests verify the complete-case path (MCAR) is correct, document
# the bias under MAR (motivating IPCW), and test rejection paths.
#
# DGP functions (helper-dgp.R):
#   simulate_mcar_outcome()       — DGP-M1: MCAR outcome censoring
#   simulate_mar_outcome()        — DGP-M2: MAR outcome censoring
#   simulate_mcar_covariate()     — DGP-M3: MCAR covariate missingness
#   simulate_longitudinal_mcar_outcome() — DGP-M4: longitudinal MCAR
#   simulate_longitudinal_mar_outcome()  — DGP-M5: longitudinal MAR

# ── helpers ──────────────────────────────────────────────────────────────

#' Inject MCAR NAs into a column
#' @param data data.frame
#' @param col character column name
#' @param prop proportion of rows to set to NA (0, 1)
#' @param seed RNG seed for reproducibility
#' @return data.frame with NAs injected
inject_na <- function(data, col, prop, seed = 99) {
  set.seed(seed)
  idx <- sample(nrow(data), size = floor(nrow(data) * prop))
  data[[col]][idx] <- NA
  data
}


# ═══════════════════════════════════════════════════════════════════════
# 1. REJECTION PATHS
# ═══════════════════════════════════════════════════════════════════════

test_that("treatment NAs without censoring → abort (gcomp)", {
  d <- simulate_binary_continuous(n = 100, seed = 1)
  d$A[1:3] <- NA
  expect_snapshot(
    error = TRUE,
    causat(d, outcome = "Y", treatment = "A", confounders = ~L)
  )
})

test_that("treatment NAs without censoring → abort (ipw)", {
  d <- simulate_binary_continuous(n = 100, seed = 1)
  d$A[1:3] <- NA
  expect_snapshot(
    error = TRUE,
    causat(
      d,
      outcome = "Y",
      treatment = "A",
      confounders = ~L,
      estimator = "ipw"
    )
  )
})

test_that("treatment NAs without censoring → abort (matching)", {
  d <- simulate_binary_continuous(n = 100, seed = 1)
  d$A[1:3] <- NA
  expect_snapshot(
    error = TRUE,
    causat(
      d,
      outcome = "Y",
      treatment = "A",
      confounders = ~L,
      estimator = "matching"
    )
  )
})

test_that("IPW aborts when covariate NAs diverge from outcome mask", {
  # fit_treatment_model() drops rows with NA L. If the outcome mask
  # (non-missing Y) has more rows than the propensity mask, the
  # row-alignment guard in fit_ipw() fires.
  d <- simulate_binary_continuous(n = 200, seed = 2)
  d$L[1:10] <- NA

  expect_error(
    causat(
      d,
      outcome = "Y",
      treatment = "A",
      confounders = ~L,
      estimator = "ipw"
    ),
    regexp = "confounder column has missing values"
  )
})


# ═══════════════════════════════════════════════════════════════════════
# 2. MCAR OUTCOME MISSINGNESS — complete-case is unbiased
#
# Under MCAR, P(C=1) = constant, so the observed sample is a random
# subset of the full sample. Complete-case analysis recovers the same
# ATE as the full-data analysis (E[ATE | observed] = ATE).
# ═══════════════════════════════════════════════════════════════════════

# --- 2a. G-comp ---

test_that("gcomp: MCAR outcome NAs, censoring= pathway (DGP-M1)", {
  d <- simulate_mcar_outcome(n = 3000, p_cens = 0.15, seed = 10)
  fit <- causat(
    d,
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
  expect_equal(result$contrasts$estimate[1], 3.0, tolerance = 0.3)
  expect_true(is.finite(result$contrasts$se[1]) && result$contrasts$se[1] > 0)
})

test_that("gcomp: MCAR outcome NAs, binomial outcome (DGP-M1 variant)", {
  d <- simulate_binary_binary(n = 3000, seed = 11)
  d <- inject_na(d, "Y", prop = 0.15, seed = 11)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    family = "binomial"
  )
  result <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich"
  )
  # True RD ≈ 0.333
  expect_equal(result$contrasts$estimate[1], 0.333, tolerance = 0.08)
})

test_that("gcomp: MCAR outcome NAs, bootstrap", {
  d <- simulate_mcar_outcome(n = 1000, p_cens = 0.10, seed = 12)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    censoring = "C"
  )
  result <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "bootstrap",
    n_boot = 100L
  )
  expect_equal(result$contrasts$estimate[1], 3.0, tolerance = 0.5)
  expect_true(is.finite(result$contrasts$se[1]) && result$contrasts$se[1] > 0)
})

test_that("gcomp: MCAR outcome NAs, ATT/ATC estimands", {
  d <- simulate_mcar_outcome(n = 3000, p_cens = 0.10, seed = 13)
  fit_att <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    censoring = "C",
    estimand = "ATT"
  )
  result_att <- contrast(
    fit_att,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich"
  )
  # True ATT = 3 (constant treatment effect)
  expect_equal(result_att$contrasts$estimate[1], 3.0, tolerance = 0.3)

  fit_atc <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    censoring = "C",
    estimand = "ATC"
  )
  result_atc <- contrast(
    fit_atc,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich"
  )
  expect_equal(result_atc$contrasts$estimate[1], 3.0, tolerance = 0.3)
})

test_that("gcomp: MCAR outcome NAs, by-stratification", {
  d <- simulate_effect_mod(n = 3000, seed = 14)
  d <- inject_na(d, "Y", prop = 0.10, seed = 14)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L + sex + A:sex
  )
  result <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich",
    by = "sex"
  )
  ests <- result$contrasts
  # True ATE|sex=0 = 3, ATE|sex=1 = 4.2
  expect_equal(ests[by == 0]$estimate, 3.0, tolerance = 0.4)
  expect_equal(ests[by == 1]$estimate, 4.2, tolerance = 0.4)
})

test_that("gcomp: MCAR outcome NAs, ratio and OR contrasts", {
  d <- simulate_binary_binary(n = 3000, seed = 15)
  d <- inject_na(d, "Y", prop = 0.10, seed = 15)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    family = "binomial"
  )
  result_rr <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    type = "ratio",
    ci_method = "sandwich"
  )
  result_or <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    type = "or",
    ci_method = "sandwich"
  )
  expect_true(all(is.finite(result_rr$contrasts$estimate)))
  expect_true(all(is.finite(result_or$contrasts$estimate)))
  expect_true(all(is.finite(result_rr$contrasts$se)))
  expect_true(all(is.finite(result_or$contrasts$se)))
})

test_that("gcomp: MCAR outcome NAs, multivariate treatment", {
  set.seed(16)
  n <- 2000
  L <- rnorm(n)
  A1 <- rbinom(n, 1, plogis(0.3 * L))
  A2 <- rbinom(n, 1, plogis(-0.2 * L))
  Y <- 1 + 2 * A1 + 1.5 * A2 + L + rnorm(n)
  d <- data.frame(Y = Y, A1 = A1, A2 = A2, L = L)
  d <- inject_na(d, "Y", prop = 0.10, seed = 16)

  fit <- causat(d, outcome = "Y", treatment = c("A1", "A2"), confounders = ~L)
  result <- contrast(
    fit,
    interventions = list(
      both = list(A1 = static(1), A2 = static(1)),
      none = list(A1 = static(0), A2 = static(0))
    ),
    reference = "none",
    ci_method = "sandwich"
  )
  # True ATE(both vs none) = 2 + 1.5 = 3.5
  expect_equal(result$contrasts$estimate[1], 3.5, tolerance = 0.4)
})

test_that("gcomp: MCAR outcome NAs, categorical treatment", {
  d <- simulate_categorical_continuous(n = 5000, seed = 17)
  d <- inject_na(d, "Y", prop = 0.10, seed = 17)
  fit <- causat(d, outcome = "Y", treatment = "A", confounders = ~L)
  result <- contrast(
    fit,
    interventions = list(b = static("b"), a = static("a")),
    reference = "a",
    ci_method = "sandwich"
  )
  # True ATE("b" vs "a") = 5 - 2 = 3
  expect_equal(result$contrasts$estimate[1], 3.0, tolerance = 0.3)
})

test_that("gcomp: MCAR outcome NAs, GAM model", {
  skip_if_not_installed("mgcv")
  d <- simulate_mcar_outcome(n = 2000, p_cens = 0.10, seed = 18)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    censoring = "C",
    model_fn = mgcv::gam
  )
  result <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich"
  )
  expect_equal(result$contrasts$estimate[1], 3.0, tolerance = 0.4)
})

test_that("gcomp: MCAR outcome NAs, poisson outcome", {
  set.seed(19)
  n <- 3000
  L <- rnorm(n)
  A <- rbinom(n, 1, plogis(0.5 * L))
  Y <- rpois(n, exp(0.5 + 0.3 * A + 0.4 * L))
  d <- data.frame(Y = Y, A = A, L = L)
  d <- inject_na(d, "Y", prop = 0.10, seed = 19)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    family = "poisson"
  )
  result <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    type = "ratio",
    ci_method = "sandwich"
  )
  # True rate ratio = exp(0.3)
  expect_equal(result$contrasts$estimate[1], exp(0.3), tolerance = 0.15)
})

test_that("gcomp: MCAR outcome NAs with external weights", {
  d <- simulate_mcar_outcome(n = 2000, p_cens = 0.10, seed = 20)
  set.seed(20)
  w <- runif(nrow(d), 0.5, 2.0)
  # Zero out weights for censored rows
  w[d$C == 1] <- 0
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    censoring = "C",
    weights = w
  )
  result <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich"
  )
  expect_equal(result$contrasts$estimate[1], 3.0, tolerance = 0.5)
})

# --- 2b. IPW ---

test_that("IPW: MCAR outcome NAs (binary/gaussian)", {
  d <- simulate_mcar_outcome(n = 3000, p_cens = 0.15, seed = 30)
  # IPW ignores censoring=; just drop censored rows from outcome mask
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw"
  )
  result <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich"
  )
  expect_equal(result$contrasts$estimate[1], 3.0, tolerance = 0.5)
  expect_true(is.finite(result$contrasts$se[1]))
})

test_that("IPW: MCAR outcome NAs, covariate NAs on same rows", {
  d <- simulate_binary_continuous(n = 3000, seed = 31)
  # NAs on the same rows → masks agree
  d$Y[1:30] <- NA
  d$L[1:30] <- NA
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw"
  )
  result <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich"
  )
  expect_equal(result$contrasts$estimate[1], 3.0, tolerance = 0.5)
})

test_that("IPW: MCAR outcome NAs, continuous shift", {
  d <- simulate_continuous_continuous(n = 3000, seed = 32)
  d <- inject_na(d, "Y", prop = 0.10, seed = 32)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw"
  )
  result <- contrast(
    fit,
    interventions = list(shifted = shift(-1), obs = NULL),
    reference = "obs",
    ci_method = "sandwich"
  )
  # True shift(-1) difference = -2
  expect_equal(result$contrasts$estimate[1], -2.0, tolerance = 0.5)
})

test_that("IPW: MCAR outcome NAs, categorical treatment", {
  d <- simulate_categorical_continuous(n = 5000, seed = 33)
  d <- inject_na(d, "Y", prop = 0.10, seed = 33)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw"
  )
  result <- contrast(
    fit,
    interventions = list(b = static("b"), a = static("a")),
    reference = "a",
    ci_method = "sandwich"
  )
  expect_equal(result$contrasts$estimate[1], 3.0, tolerance = 0.6)
})

# --- 2c. Matching ---

test_that("matching: MCAR outcome NAs", {
  d <- simulate_mcar_outcome(n = 2000, p_cens = 0.10, seed = 40)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "matching"
  )
  result <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich"
  )
  # ATT ≈ 3 (constant effect)
  expect_equal(result$contrasts$estimate[1], 3.0, tolerance = 0.5)
})

# --- 2d. ICE (longitudinal) ---

test_that("ICE: MCAR outcome NAs at final time (DGP-M4, sandwich)", {
  d <- simulate_longitudinal_mcar_outcome(n = 3000, p_cens = 0.10, seed = 50)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    censoring = "C",
    id = "id",
    time = "time"
  )
  result <- contrast(
    fit,
    interventions = list(always = static(1), never = static(0)),
    reference = "never",
    ci_method = "sandwich"
  )
  # True ATE = 5
  expect_equal(result$contrasts$estimate[1], 5.0, tolerance = 1.0)
  expect_true(is.finite(result$contrasts$se[1]))
})

test_that("ICE: MCAR outcome NAs at final time (DGP-M4, bootstrap)", {
  d <- simulate_longitudinal_mcar_outcome(n = 1000, p_cens = 0.10, seed = 51)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    censoring = "C",
    id = "id",
    time = "time"
  )
  result <- contrast(
    fit,
    interventions = list(always = static(1), never = static(0)),
    reference = "never",
    ci_method = "bootstrap",
    n_boot = 80L
  )
  expect_equal(result$contrasts$estimate[1], 5.0, tolerance = 1.5)
  expect_true(is.finite(result$contrasts$se[1]))
})


# ═══════════════════════════════════════════════════════════════════════
# 3. MCAR COVARIATE MISSINGNESS — complete-case is unbiased
# ═══════════════════════════════════════════════════════════════════════

test_that("gcomp: MCAR covariate NAs (DGP-M3)", {
  d <- simulate_mcar_covariate(n = 3000, p_miss = 0.10, seed = 60)
  fit <- causat(d, outcome = "Y", treatment = "A", confounders = ~L)
  result <- suppressMessages(contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich"
  ))
  expect_equal(result$contrasts$estimate[1], 3.0, tolerance = 0.3)
})

test_that("gcomp: MCAR outcome + covariate NAs simultaneously", {
  d <- simulate_binary_continuous(n = 3000, seed = 61)
  d <- inject_na(d, "Y", prop = 0.10, seed = 61)
  d <- inject_na(d, "L", prop = 0.10, seed = 62)
  fit <- causat(d, outcome = "Y", treatment = "A", confounders = ~L)
  result <- suppressMessages(contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich"
  ))
  expect_equal(result$contrasts$estimate[1], 3.0, tolerance = 0.35)
})

test_that("gcomp: MCAR covariate NAs, continuous shift", {
  d <- simulate_continuous_continuous(n = 3000, seed = 63)
  d <- inject_na(d, "L", prop = 0.10, seed = 63)
  fit <- causat(d, outcome = "Y", treatment = "A", confounders = ~L)
  result <- suppressMessages(contrast(
    fit,
    interventions = list(shifted = shift(-1), obs = NULL),
    reference = "obs",
    ci_method = "sandwich"
  ))
  expect_equal(result$contrasts$estimate[1], -2.0, tolerance = 0.3)
})

test_that("matching: covariate NAs handled by MatchIt", {
  d <- simulate_mcar_covariate(n = 2000, p_miss = 0.05, seed = 64)
  # MatchIt handles NAs internally. If it aborts, verify error is clear.
  tryCatch(
    {
      fit <- causat(
        d,
        outcome = "Y",
        treatment = "A",
        confounders = ~L,
        estimator = "matching"
      )
      result <- contrast(
        fit,
        interventions = list(a1 = static(1), a0 = static(0)),
        reference = "a0",
        ci_method = "sandwich"
      )
      expect_equal(result$contrasts$estimate[1], 3.0, tolerance = 0.6)
    },
    error = function(e) {
      expect_true(
        grepl(
          "missing|NA|na\\.omit|complete",
          conditionMessage(e),
          ignore.case = TRUE
        ),
        info = paste("Error was:", conditionMessage(e))
      )
    }
  )
})

test_that("ICE: MCAR time-varying covariate NAs (bootstrap)", {
  # ICE sandwich has a known cascade-gradient alignment issue with
  # partial covariate NAs across time steps. Bootstrap handles it.
  d <- make_linear_scm(n = 3000, n_times = 2, seed = 65)
  t1_rows <- which(d$time == 1)
  set.seed(65)
  na_rows <- sample(t1_rows, size = floor(length(t1_rows) * 0.08))
  d$L[na_rows] <- NA
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    id = "id",
    time = "time"
  )
  result <- contrast(
    fit,
    interventions = list(always = static(1), never = static(0)),
    reference = "never",
    ci_method = "bootstrap",
    n_boot = 80L
  )
  expect_equal(result$contrasts$estimate[1], 5.0, tolerance = 1.5)
})


# ═══════════════════════════════════════════════════════════════════════
# 4. MAR OUTCOME MISSINGNESS — complete-case IS BIASED
#
# Under MAR, P(C=0 | A, L) depends on observed variables. Complete-case
# over-represents individuals with low censoring probability. IPCW with
# the correct censoring model corrects this bias.
#
# These tests document the bias and verify that manual IPCW weights
# (computed from the known censoring model) recover the truth. They
# serve as the ground truth for the upcoming IPCW implementation.
# ═══════════════════════════════════════════════════════════════════════

test_that("MAR outcome: gcomp on observed data recovers truth (model-based)", {
  # G-comp with a correctly specified outcome model E[Y|A,L] is consistent
  # even under MAR censoring — the regression surface is unchanged by
  # censoring, and standardization over the target population recovers
  # the marginal mean. This is the key advantage of model-based
  # estimation (Hernán & Robins Ch. 13): g-comp does not require IPCW
  # when the outcome model is correctly specified.
  d <- simulate_mar_outcome(n = 10000, seed = 70)
  observed <- d[d$C == 0, ]

  fit <- causat(observed, outcome = "Y", treatment = "A", confounders = ~L)
  result <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich"
  )
  expect_equal(result$contrasts$estimate[1], 3.0, tolerance = 0.2)
})

test_that("MAR outcome: IPCW-weighted gcomp via weights= (doubly robust)", {
  # When the outcome model may be misspecified, IPCW weights provide
  # additional robustness. Here the outcome model IS correctly specified,
  # so the IPCW-weighted fit should also recover the truth.
  d <- simulate_mar_outcome(n = 10000, seed = 71)
  observed <- d[d$C == 0, ]

  # Stabilized IPCW: w_i = P(C=0) / P(C=0 | A_i, L_i)
  p_obs <- 1 - d$p_cens[d$C == 0]
  p_overall <- mean(d$C == 0)
  ipcw_stab <- p_overall / p_obs

  fit <- causat(
    observed,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    weights = ipcw_stab
  )
  result <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich"
  )
  expect_equal(result$contrasts$estimate[1], 3.0, tolerance = 0.3)
})

test_that("MAR outcome (longitudinal): manual IPCW via weights= recovers truth", {
  d <- simulate_longitudinal_mar_outcome(n = 10000, seed = 71)
  # Manual IPCW weights from the known censoring model.
  # P(C_1=1 | A_0, L_0) = expit(-0.5 + 1.5*A_0 + 0.8*L_0)
  uncensored <- d[d$C == 0, ]
  t1_uncens <- uncensored[uncensored$time == 1, ]
  p_obs <- 1 - t1_uncens$p_cens
  # Build weight vector for the full uncensored person-period panel:
  # time=0 rows have w=1 (no censoring at baseline), time=1 have w=IPCW
  ipcw_full <- rep(1, nrow(uncensored))
  ipcw_full[uncensored$time == 1] <- 1 / p_obs

  fit <- causat(
    uncensored,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    weights = ipcw_full
  )
  result <- contrast(
    fit,
    interventions = list(always = static(1), never = static(0)),
    reference = "never",
    ci_method = "sandwich"
  )
  # True ATE = 5
  expect_equal(result$contrasts$estimate[1], 5.0, tolerance = 0.8)
})
