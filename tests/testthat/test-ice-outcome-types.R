# Tests for non-gaussian/non-binomial outcome families with longitudinal ICE.
# Phase 13 added Poisson / Gamma / NB / betareg for point gcomp and IPW;
# this file covers their longitudinal (ICE backward-recursion) paths.
#
# Source change: `family_pseudo` switch in `fit_ice()` (`R/ice.R`) now maps
# "poisson" -> quasipoisson so that intermediate pseudo-outcome steps do not
# trigger "non-integer counts" warnings. Gamma, NB, and betareg are
# unaffected by the switch (they pass through or ignore the family arg).
#
# Oracle strategy:
#   Poisson  — analytical truth (log-normal moment formula) + Python M-estimation
#              sandwich + lmtp_sdr count
#   Gamma    — lmtp_sdr continuous
#   NB       — lmtp_sdr count (same Poisson DGP, MASS::glm.nb model_fn)
#   betareg  — smoke only (no standard oracle)

# ---------------------------------------------------------------------------
# Poisson
# ---------------------------------------------------------------------------

test_that("longitudinal ICE: Poisson — point estimates match analytical truth", {
  # Truth: E[Y^{A=1}] = exp(2.2) ~ 9.025, E[Y^{A=0}] = exp(0.6) ~ 1.822
  # (Derivation in helper-dgp.R DGP-ICE-POIS comment block)
  dat <- simulate_longitudinal_poisson(n = 3000L, seed = 1L)
  fit <- causat(
    dat,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    id = "id",
    time = "t",
    estimator = "gcomp",
    family = "poisson"
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1L), a0 = static(0L)),
    reference = "a0",
    type = "difference",
    ci_method = "sandwich"
  )
  ate <- res$contrasts$estimate[1]
  est <- res$estimates

  expect_equal(
    ate,
    exp(2.2) - exp(0.6),
    tolerance = 0.35,
    label = "ICE Poisson ATE vs analytical truth"
  )
  expect_equal(
    est$estimate[est$intervention == "a1"],
    exp(2.2),
    tolerance = 0.25,
    label = "ICE Poisson E[Y^{A=1}] vs truth"
  )
  expect_equal(
    est$estimate[est$intervention == "a0"],
    exp(0.6),
    tolerance = 0.15,
    label = "ICE Poisson E[Y^{A=0}] vs truth"
  )
  expect_true(is.finite(res$contrasts$ci_lower[1]))
  expect_true(is.finite(res$contrasts$ci_upper[1]))
  # 95% CI must cover the true ATE
  expect_true(
    res$contrasts$ci_lower[1] < exp(2.2) - exp(0.6) &&
      exp(2.2) - exp(0.6) < res$contrasts$ci_upper[1],
    label = "ICE Poisson sandwich CI covers true ATE"
  )
})

test_that("longitudinal ICE: Poisson — no non-integer counts warnings at pseudo-steps", {
  # The `family_pseudo` switch must map poisson -> quasipoisson for pseudo-steps.
  # Any "non-integer counts" warning here means the switch is broken.
  dat <- simulate_longitudinal_poisson(n = 500L, seed = 2L)
  fit <- causat(
    dat,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    id = "id",
    time = "t",
    estimator = "gcomp",
    family = "poisson"
  )
  saw_integer_warning <- FALSE
  withCallingHandlers(
    contrast(
      fit,
      interventions = list(a1 = static(1L)),
      type = "difference",
      ci_method = "sandwich"
    ),
    warning = function(w) {
      if (grepl("non-integer", conditionMessage(w), fixed = TRUE)) {
        saw_integer_warning <<- TRUE
      }
      invokeRestart("muffleWarning")
    }
  )
  expect_false(
    saw_integer_warning,
    label = "No 'non-integer counts' from Poisson ICE pseudo-steps"
  )
})

test_that("longitudinal ICE: Poisson — Python M-estimation cross-check", {
  fix_path <- testthat::test_path(
    "fixtures",
    "python",
    "ice_poisson_tau2_results.csv"
  )
  skip_if(!file.exists(fix_path), "Python ICE Poisson fixture absent")

  py <- read.csv(fix_path)
  dat <- simulate_longitudinal_poisson(n = 300L, seed = 2025L)
  fit <- causat(
    dat,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    id = "id",
    time = "t",
    estimator = "gcomp",
    family = "poisson"
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1L), a0 = static(0L)),
    reference = "a0",
    type = "difference",
    ci_method = "sandwich"
  )
  r_a1 <- res$estimates[res$estimates$intervention == "a1", ]
  r_a0 <- res$estimates[res$estimates$intervention == "a0", ]
  py_a1 <- py[py$arm == "a1", ]
  py_a0 <- py[py$arm == "a0", ]

  # Point estimates: same data, small deviation from different column handling
  expect_equal(
    r_a1$estimate,
    py_a1$estimate,
    tolerance = 1e-2,
    label = "ICE Poisson a1 estimate vs Python M-estimation"
  )
  expect_equal(
    r_a0$estimate,
    py_a0$estimate,
    tolerance = 1e-2,
    label = "ICE Poisson a0 estimate vs Python M-estimation"
  )
  # SEs: same plug-in sandwich should agree within 1%
  expect_equal(
    r_a1$se,
    py_a1$se,
    tolerance = 1e-2,
    label = "ICE Poisson a1 SE vs Python M-estimation"
  )
  expect_equal(
    r_a0$se,
    py_a0$se,
    tolerance = 1e-2,
    label = "ICE Poisson a0 SE vs Python M-estimation"
  )
})

test_that("longitudinal ICE: Poisson — lmtp cross-check (static binary)", {
  skip_if_fast()
  skip_if_not_installed("lmtp")

  dat <- simulate_longitudinal_poisson(n = 2000L, seed = 3L)
  fit <- causat(
    dat,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    id = "id",
    time = "t",
    estimator = "gcomp",
    family = "poisson"
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1L), a0 = static(0L)),
    reference = "a0",
    type = "difference",
    ci_method = "sandwich"
  )
  ate_causatr <- res$contrasts$estimate[1]

  # Reshape to wide for lmtp
  d_wide <- reshape(
    as.data.frame(dat),
    idvar = "id",
    timevar = "t",
    direction = "wide",
    v.names = c("A", "L", "Y"),
    sep = "_"
  )
  d_wide <- d_wide[, c("id", "L0", "A_1", "A_2", "L_1", "L_2", "Y_2")]

  theta_of <- function(r) tryCatch(r$estimate@x, error = function(e) r$theta)
  run_lmtp <- function(shift_fn) {
    tryCatch(
      suppressWarnings(suppressMessages(lmtp::lmtp_sdr(
        data = d_wide,
        trt = c("A_1", "A_2"),
        outcome = "Y_2",
        baseline = "L0",
        time_vary = list(c("L_1"), c("L_2")),
        shift = shift_fn,
        outcome_type = "count",
        learners_trt = "SL.glm",
        learners_outcome = "SL.glm",
        folds = 1
      ))),
      error = function(e) NULL
    )
  }
  r_always <- run_lmtp(function(data, trt) rep(1, nrow(data)))
  r_never <- run_lmtp(function(data, trt) rep(0, nrow(data)))
  skip_if(is.null(r_always) || is.null(r_never), "lmtp::lmtp_sdr() failed")

  ate_lmtp <- theta_of(r_always) - theta_of(r_never)

  expect_equal(
    ate_causatr,
    exp(2.2) - exp(0.6),
    tolerance = 0.4,
    label = "causatr ICE Poisson vs analytical truth"
  )
  expect_equal(
    ate_lmtp,
    exp(2.2) - exp(0.6),
    tolerance = 0.5,
    label = "lmtp_sdr Poisson vs analytical truth"
  )
  expect_equal(
    ate_causatr,
    ate_lmtp,
    tolerance = 0.5,
    label = "causatr ICE Poisson vs lmtp_sdr"
  )
})

# ---------------------------------------------------------------------------
# Gamma
# ---------------------------------------------------------------------------

test_that("longitudinal ICE: Gamma — lmtp cross-check (continuous shift)", {
  skip_if_fast()
  skip_if_not_installed("lmtp")

  dat <- simulate_longitudinal_gamma(n = 2000L, seed = 4L)
  fit <- causat(
    dat,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    id = "id",
    time = "t",
    estimator = "gcomp",
    family = Gamma(link = "log")
  )
  res <- contrast(
    fit,
    interventions = list(shifted = shift(0.5), natural = shift(0)),
    reference = "natural",
    type = "difference",
    ci_method = "sandwich"
  )
  ate_causatr <- res$contrasts$estimate[1]

  d_wide <- reshape(
    as.data.frame(dat),
    idvar = "id",
    timevar = "t",
    direction = "wide",
    v.names = c("A", "L", "Y"),
    sep = "_"
  )
  d_wide <- d_wide[, c("id", "L0", "A_1", "A_2", "L_1", "L_2", "Y_2")]

  theta_of <- function(r) tryCatch(r$estimate@x, error = function(e) r$theta)
  run_lmtp <- function(shift_fn) {
    tryCatch(
      suppressWarnings(suppressMessages(lmtp::lmtp_sdr(
        data = d_wide,
        trt = c("A_1", "A_2"),
        outcome = "Y_2",
        baseline = "L0",
        time_vary = list(c("L_1"), c("L_2")),
        shift = shift_fn,
        outcome_type = "continuous",
        learners_trt = "SL.glm",
        learners_outcome = "SL.glm",
        folds = 1
      ))),
      error = function(e) NULL
    )
  }
  r_shifted <- run_lmtp(function(data, trt) data[[trt]] + 0.5)
  r_natural <- run_lmtp(function(data, trt) data[[trt]])
  skip_if(is.null(r_shifted) || is.null(r_natural), "lmtp::lmtp_sdr() failed")

  ate_lmtp <- theta_of(r_shifted) - theta_of(r_natural)

  expect_true(is.finite(ate_causatr), label = "causatr ICE Gamma ATE is finite")
  expect_true(is.finite(ate_lmtp), label = "lmtp_sdr Gamma ATE is finite")
  # Both should be positive (shift increases Y on log-link Gamma)
  expect_true(ate_causatr > 0, label = "causatr ICE Gamma ATE is positive")
  expect_equal(
    ate_causatr,
    ate_lmtp,
    tolerance = 0.4,
    label = "causatr ICE Gamma vs lmtp_sdr"
  )
})

test_that("longitudinal ICE: Gamma — sandwich CI is finite and SE is positive", {
  dat <- simulate_longitudinal_gamma(n = 1000L, seed = 5L)
  fit <- causat(
    dat,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    id = "id",
    time = "t",
    estimator = "gcomp",
    family = Gamma(link = "log")
  )
  res <- contrast(
    fit,
    interventions = list(shifted = shift(0.5), natural = shift(0)),
    reference = "natural",
    type = "difference",
    ci_method = "sandwich"
  )
  expect_true(all(is.finite(res$estimates$estimate)))
  expect_true(all(is.finite(res$estimates$se)))
  expect_true(all(res$estimates$se > 0))
  expect_true(all(is.finite(res$contrasts$ci_lower)))
  expect_true(all(is.finite(res$contrasts$ci_upper)))
})

# ---------------------------------------------------------------------------
# Negative binomial (MASS::glm.nb)
# ---------------------------------------------------------------------------

test_that("longitudinal ICE: MASS::glm.nb — lmtp cross-check (static binary)", {
  skip_if_fast()
  skip_if_not_installed("lmtp")
  skip_if_not_installed("MASS")

  # NB is an over-dispersed Poisson; same marginal mean DGP; same lmtp oracle.
  dat <- simulate_longitudinal_poisson(n = 2000L, seed = 12L)
  fit <- causat(
    dat,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    id = "id",
    time = "t",
    estimator = "gcomp",
    model_fn = MASS::glm.nb
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1L), a0 = static(0L)),
    reference = "a0",
    type = "difference",
    ci_method = "sandwich"
  )
  ate_causatr <- res$contrasts$estimate[1]

  d_wide <- reshape(
    as.data.frame(dat),
    idvar = "id",
    timevar = "t",
    direction = "wide",
    v.names = c("A", "L", "Y"),
    sep = "_"
  )
  d_wide <- d_wide[, c("id", "L0", "A_1", "A_2", "L_1", "L_2", "Y_2")]

  theta_of <- function(r) tryCatch(r$estimate@x, error = function(e) r$theta)
  run_lmtp <- function(shift_fn) {
    tryCatch(
      suppressWarnings(suppressMessages(lmtp::lmtp_sdr(
        data = d_wide,
        trt = c("A_1", "A_2"),
        outcome = "Y_2",
        baseline = "L0",
        time_vary = list(c("L_1"), c("L_2")),
        shift = shift_fn,
        outcome_type = "count",
        learners_trt = "SL.glm",
        learners_outcome = "SL.glm",
        folds = 1
      ))),
      error = function(e) NULL
    )
  }
  r_always <- run_lmtp(function(data, trt) rep(1, nrow(data)))
  r_never <- run_lmtp(function(data, trt) rep(0, nrow(data)))
  skip_if(is.null(r_always) || is.null(r_never), "lmtp::lmtp_sdr() failed")

  ate_lmtp <- theta_of(r_always) - theta_of(r_never)

  expect_equal(
    ate_causatr,
    exp(2.2) - exp(0.6),
    tolerance = 0.4,
    label = "causatr ICE NB vs analytical truth"
  )
  expect_equal(
    ate_causatr,
    ate_lmtp,
    tolerance = 0.5,
    label = "causatr ICE NB vs lmtp_sdr"
  )
})

# ---------------------------------------------------------------------------
# betareg
# ---------------------------------------------------------------------------

test_that("longitudinal ICE: betareg — smoke test (point + bootstrap)", {
  skip_if_fast()
  skip_if_not_installed("betareg")

  dat <- simulate_longitudinal_betareg(n = 500L, seed = 7L)
  # betareg::betareg uses optim directly and does not handle rank-deficient
  # design matrices (unlike stats::glm which drops redundant columns). In this
  # DGP, L at t=1 equals L0, so including both L0 and L as confounders creates
  # perfect collinearity at the pseudo-step. Use only the TV confounder.
  # betareg's mu/precision coefficient list is non-conformable with the analytic
  # ICE score matrix; the supported variance path is bootstrap (matching the
  # point-gcomp Tier 1 numeric fallback for betareg).
  fit <- causat(
    dat,
    outcome = "Y",
    treatment = "A",
    confounders_tv = ~L,
    id = "id",
    time = "t",
    estimator = "gcomp",
    model_fn = betareg::betareg
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1L), a0 = static(0L)),
    reference = "a0",
    type = "difference",
    ci_method = "bootstrap",
    n_boot = 50L
  )
  ate <- res$contrasts$estimate[1]
  se <- res$estimates$se

  expect_true(is.finite(ate), label = "betareg ICE ATE is finite")
  expect_true(ate > 0, label = "betareg ICE ATE is positive (A=1 > A=0)")
  expect_true(all(is.finite(se)), label = "betareg ICE SEs are finite")
  expect_true(all(se > 0), label = "betareg ICE SEs are positive")
  expect_true(
    is.finite(res$contrasts$ci_lower[1]) &&
      is.finite(res$contrasts$ci_upper[1]),
    label = "betareg ICE CI bounds are finite"
  )
})
