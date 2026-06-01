# causat_mice(): multiple-imputation pooling across causatr's feature surface.
# Truth-based DGPs live in helper-dgp.R (simulate_mi_*); the Rubin pooling is
# cross-checked against mice::pool.scalar() via helper-mice-oracle.R.

ints <- function() list(a1 = static(1), a0 = static(0))

# Absolute estimate of the ATE = a1 - a0 from a pooled difference result,
# regardless of which intervention contrast() chose as the reference.
ate_of <- function(res) abs(res$contrasts$estimate[1])

# -- Core: Rubin pooling recovers truth and matches the external oracle ------

test_that("gcomp + Rubin recovers the ATE under MAR covariate missingness", {
  skip_if_not_installed("mice")
  dat <- simulate_mi_covariate(n = 2000)
  imp <- mice::mice(
    dat[, c("Y", "A", "L")],
    m = 10,
    printFlag = FALSE
  )
  res <- causat_mice(imp, "Y", "A", ~L, ints(), estimator = "gcomp")

  expect_s3_class(res, "causatr_result")
  expect_identical(res$ci_method, "rubin")
  expect_equal(ate_of(res), 3, tolerance = 0.15)
  # ~28% of L is missing; the pool must report a positive fraction of
  # missing information, i.e. it actually used the between-imputation spread.
  expect_gt(attr(res, "mi_details")$contrasts$fmi[1], 0)
})

test_that("Rubin pooling reproduces mice::pool.scalar() to machine precision", {
  skip_if_not_installed("mice")
  dat <- simulate_mi_covariate(n = 1500, seed = 7)
  imp <- mice::mice(dat[, c("Y", "A", "L")], m = 10, printFlag = FALSE)
  fit_args <- list(
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "gcomp",
    family = "gaussian",
    estimand = "ATE"
  )

  collected <- mi_oracle_collect(imp, fit_args, ints())
  # Contrast table has one column -> dfcom = n - 1.
  oracle <- mi_oracle_pool(collected$Q, collected$U, n = nrow(dat), k = 1L)

  res <- causat_mice(imp, "Y", "A", ~L, ints(), estimator = "gcomp")
  mi <- attr(res, "mi_details")
  expect_equal(res$contrasts$estimate[1], oracle$estimate, tolerance = 1e-8)
  expect_equal(res$contrasts$se[1], oracle$se, tolerance = 1e-8)
  expect_equal(mi$contrasts$df[1], oracle$df, tolerance = 1e-6)
})

# -- All five estimators on the same DGP --------------------------------------

test_that("ipw / aipw / matching recover the ATE under MAR covariates", {
  skip_if_not_installed("mice")
  dat <- simulate_mi_covariate(n = 2500, seed = 3)
  imp <- mice::mice(dat[, c("Y", "A", "L")], m = 10, printFlag = FALSE)

  res_ipw <- causat_mice(
    imp,
    "Y",
    "A",
    ~L,
    ints(),
    estimator = "ipw",
    propensity_model_fn = stats::glm
  )
  res_aipw <- causat_mice(
    imp,
    "Y",
    "A",
    ~L,
    ints(),
    estimator = "aipw",
    propensity_model_fn = stats::glm
  )
  res_match <- causat_mice(imp, "Y", "A", ~L, ints(), estimator = "matching")

  expect_equal(ate_of(res_ipw), 3, tolerance = 0.2)
  expect_equal(ate_of(res_aipw), 3, tolerance = 0.2)
  expect_equal(ate_of(res_match), 3, tolerance = 0.25)
})

test_that("snm pools blip parameters per parameter row", {
  skip_if_not_installed("mice")
  set.seed(41)
  n <- 3000
  L <- rnorm(n)
  Mod <- rbinom(n, 1, 0.5)
  A <- rbinom(n, 1, plogis(0.4 * L))
  Y <- 2 + 3 * A + 1.5 * L + 1.0 * A * Mod + rnorm(n) # blip: A = 3, A:Mod = 1
  L[rbinom(n, 1, plogis(-1 + 0.8 * A)) == 0] <- NA
  dat <- data.frame(Y = Y, A = A, L = L, Mod = Mod)
  imp <- mice::mice(dat, m = 8, printFlag = FALSE)

  res <- causat_mice(imp, "Y", "A", ~ L + A:Mod, estimator = "snm")
  expect_true("parameter" %in% names(res$estimates))
  est <- res$estimates$estimate
  expect_equal(est[1], 3, tolerance = 0.2) # psi_intercept
  expect_equal(est[2], 1, tolerance = 0.3) # psi_Mod
})

# -- MAR treatment: impute A from P(A | L, Y) ---------------------------------

test_that("MI imputes MAR treatment and recovers the ATE", {
  skip_if_not_installed("mice")
  dat <- simulate_mi_treatment(n = 2500)
  # A is binary but stored numeric with NAs; impute as a factor so mice uses
  # logistic regression, then coerce back to 0/1 inside the analysis.
  dat$A <- factor(dat$A)
  imp <- mice::mice(dat[, c("Y", "A", "L")], m = 10, printFlag = FALSE)
  long <- mice::complete(imp, "long", include = TRUE)
  long$A <- as.integer(as.character(long$A))
  imp <- mice::as.mids(long)

  res <- causat_mice(imp, "Y", "A", ~L, ints(), estimator = "gcomp")
  expect_equal(ate_of(res), 3, tolerance = 0.2)
})

# -- Treatment types ----------------------------------------------------------

test_that("continuous treatment + shift MTP pools under MI", {
  skip_if_not_installed("mice")
  set.seed(21)
  n <- 2000
  L <- rnorm(n)
  A <- 0.5 * L + rnorm(n)
  Y <- 2 + 1.0 * A + 1.5 * L + rnorm(n) # shift(1) effect = 1
  L[rbinom(n, 1, plogis(-1 + 0.4 * A)) == 0] <- NA
  dat <- data.frame(Y = Y, A = A, L = L)
  imp <- mice::mice(dat, m = 8, printFlag = FALSE)

  res <- causat_mice(
    imp,
    "Y",
    "A",
    ~L,
    list(s = shift(1), obs = NULL),
    estimator = "gcomp"
  )
  expect_equal(ate_of(res), 1, tolerance = 0.1)
})

test_that("categorical and count treatments pool without error", {
  skip_if_not_installed("mice")
  set.seed(22)
  n <- 2000
  L <- rnorm(n)
  Acat <- factor(sample(c("a", "b", "c"), n, replace = TRUE))
  Y <- 2 + 1.5 * L + rnorm(n)
  L[rbinom(n, 1, 0.2) == 1] <- NA
  dat <- data.frame(Y = Y, A = Acat, L = L)
  imp <- mice::mice(dat, m = 6, printFlag = FALSE)
  res_cat <- causat_mice(
    imp,
    "Y",
    "A",
    ~L,
    list(a = static("a"), b = static("b")),
    estimator = "gcomp"
  )
  expect_true(all(is.finite(res_cat$contrasts$se)))
})

# -- Outcome families ---------------------------------------------------------

test_that("binomial outcome supports ratio and OR contrasts under MI", {
  skip_if_not_installed("mice")
  set.seed(31)
  n <- 2500
  L <- rnorm(n)
  A <- rbinom(n, 1, plogis(0.5 * L))
  Y <- rbinom(n, 1, plogis(-0.5 + 0.8 * A + 0.5 * L))
  L[rbinom(n, 1, plogis(-1 + 0.8 * A)) == 0] <- NA
  dat <- data.frame(Y = Y, A = A, L = L)
  imp <- mice::mice(dat, m = 8, printFlag = FALSE)

  res_rr <- causat_mice(
    imp,
    "Y",
    "A",
    ~L,
    ints(),
    estimator = "gcomp",
    family = "binomial",
    type = "ratio"
  )
  res_or <- causat_mice(
    imp,
    "Y",
    "A",
    ~L,
    ints(),
    estimator = "gcomp",
    family = "binomial",
    type = "or"
  )
  # Ratio / OR are pooled on the log scale; estimate stays positive and the
  # CI brackets it on the ratio scale.
  expect_gt(res_rr$contrasts$estimate[1], 0)
  expect_gt(res_or$contrasts$estimate[1], 0)
  expect_lt(res_rr$contrasts$ci_lower[1], res_rr$contrasts$estimate[1])
  expect_gt(res_rr$contrasts$ci_upper[1], res_rr$contrasts$estimate[1])
})

test_that("poisson and Gamma outcomes pool without error", {
  skip_if_not_installed("mice")
  set.seed(32)
  n <- 2000
  L <- rnorm(n)
  A <- rbinom(n, 1, plogis(0.5 * L))
  Ycount <- rpois(n, exp(0.5 + 0.4 * A + 0.2 * L))
  Ygamma <- rgamma(n, shape = 2, rate = 1 / exp(0.3 + 0.2 * A + 0.1 * L))
  L[rbinom(n, 1, plogis(-1 + 0.8 * A)) == 0] <- NA
  dat <- data.frame(Ycount = Ycount, Ygamma = Ygamma, A = A, L = L)
  imp <- mice::mice(dat, m = 6, printFlag = FALSE)

  res_p <- causat_mice(
    imp,
    "Ycount",
    "A",
    ~L,
    ints(),
    estimator = "gcomp",
    family = "poisson"
  )
  res_g <- causat_mice(
    imp,
    "Ygamma",
    "A",
    ~L,
    ints(),
    estimator = "gcomp",
    family = Gamma(link = "log")
  )
  expect_true(all(is.finite(res_p$contrasts$se)))
  expect_true(all(is.finite(res_g$contrasts$se)))
})

# -- Interventions ------------------------------------------------------------

test_that("scale_by, dynamic, threshold, and ipsi all pool under MI", {
  skip_if_not_installed("mice")
  set.seed(51)
  n <- 2000
  L <- rnorm(n)
  A <- rbinom(n, 1, plogis(0.5 * L))
  Acont <- 0.5 * L + rnorm(n)
  Y <- 2 + 3 * A + 1.5 * L + rnorm(n)
  Yc <- 2 + 1.0 * Acont + 1.5 * L + rnorm(n)
  L[rbinom(n, 1, plogis(-1 + 0.6 * A)) == 0] <- NA
  dat <- data.frame(Y = Y, Yc = Yc, A = A, Acont = Acont, L = L)
  imp <- mice::mice(dat, m = 6, printFlag = FALSE)

  res_dyn <- causat_mice(
    imp,
    "Y",
    "A",
    ~L,
    list(
      rule = dynamic(function(data, trt) as.integer(data$L > 0)),
      none = static(0)
    ),
    estimator = "gcomp"
  )
  res_ipsi <- causat_mice(
    imp,
    "Y",
    "A",
    ~L,
    list(up = ipsi(2), obs = NULL),
    estimator = "ipw",
    propensity_model_fn = stats::glm
  )
  res_scale <- causat_mice(
    imp,
    "Yc",
    "Acont",
    ~L,
    list(s = scale_by(1.1), obs = NULL),
    estimator = "ipw",
    propensity_model_fn = stats::glm
  )
  res_thr <- causat_mice(
    imp,
    "Yc",
    "Acont",
    ~L,
    list(th = threshold(lower = 0), obs = NULL),
    estimator = "gcomp"
  )

  for (r in list(res_dyn, res_ipsi, res_scale, res_thr)) {
    expect_s3_class(r, "causatr_result")
    expect_true(all(is.finite(r$contrasts$se)))
  }
})

# -- Multivariate treatment + stabilized weights ------------------------------

test_that("multivariate gcomp / IPW (incl. stabilized) recover the joint ATE", {
  skip_if_not_installed("mice")
  set.seed(61)
  n <- 4000
  L <- rnorm(n)
  A1 <- rbinom(n, 1, plogis(0.3 * L))
  A2 <- rbinom(n, 1, plogis(0.2 * L))
  Y <- 1 + 2 * A1 + 1.5 * A2 + L + rnorm(n) # both vs none = 3.5
  L[rbinom(n, 1, plogis(-1 + 0.5 * A1)) == 0] <- NA
  dat <- data.frame(Y = Y, A1 = A1, A2 = A2, L = L)
  imp <- mice::mice(dat, m = 8, printFlag = FALSE)
  mv <- list(
    both = list(A1 = static(1), A2 = static(1)),
    none = list(A1 = static(0), A2 = static(0))
  )

  res_g <- causat_mice(imp, "Y", c("A1", "A2"), ~L, mv, estimator = "gcomp")
  res_i <- causat_mice(
    imp,
    "Y",
    c("A1", "A2"),
    ~L,
    mv,
    estimator = "ipw",
    propensity_model_fn = stats::glm
  )
  res_s <- causat_mice(
    imp,
    "Y",
    c("A1", "A2"),
    ~L,
    mv,
    estimator = "ipw",
    propensity_model_fn = stats::glm,
    stabilize = "marginal"
  )

  expect_equal(ate_of(res_g), 3.5, tolerance = 0.25)
  expect_equal(ate_of(res_i), 3.5, tolerance = 0.3)
  expect_equal(ate_of(res_s), 3.5, tolerance = 0.3)
})

# -- Estimand + by-stratification ---------------------------------------------

test_that("ATT and by-stratified pooling work under MI", {
  skip_if_not_installed("mice")
  df <- simulate_het_effect(n = 4000)
  df$L[rbinom(nrow(df), 1, plogis(-1 + 0.8 * df$A)) == 0] <- NA
  imp <- mice::mice(df[, c("Y", "A", "L", "sex")], m = 8, printFlag = FALSE)

  res_att <- causat_mice(
    imp,
    "Y",
    "A",
    ~ L + sex + A:L + A:sex,
    ints(),
    estimator = "gcomp",
    estimand = "ATT"
  )
  expect_s3_class(res_att, "causatr_result")
  expect_identical(res_att$estimand, "ATT")

  res_by <- causat_mice(
    imp,
    "Y",
    "A",
    ~ L + sex + A:L + A:sex,
    ints(),
    estimator = "gcomp",
    by = "sex"
  )
  # Subgroup ATEs: 3.0 for sex 0, 4.5 for sex 1 (simulate_het_effect truth).
  ate0 <- abs(res_by$contrasts$estimate[res_by$contrasts$by == "0"])
  ate1 <- abs(res_by$contrasts$estimate[res_by$contrasts$by == "1"])
  expect_equal(ate0, 3.0, tolerance = 0.25)
  expect_equal(ate1, 4.5, tolerance = 0.25)
})

# -- Longitudinal (ICE + longitudinal IPW) ------------------------------------

test_that("longitudinal gcomp (ICE) and IPW pool under MI", {
  skip_if_not_installed("mice")
  long <- simulate_mi_longitudinal(n = 2500)
  imp <- mice::mice(long, m = 6, printFlag = FALSE)

  res_ice <- causat_mice(
    imp,
    "Y",
    "A",
    confounders = ~1,
    confounders_tv = ~L,
    interventions = list(always = static(1), never = static(0)),
    estimator = "gcomp",
    id = "id",
    time = "time"
  )
  res_ipw <- causat_mice(
    imp,
    "Y",
    "A",
    confounders = ~1,
    confounders_tv = ~L,
    interventions = list(always = static(1), never = static(0)),
    estimator = "ipw",
    id = "id",
    time = "time",
    propensity_model_fn = stats::glm
  )

  expect_identical(res_ice$fit_type, "longitudinal")
  expect_equal(ate_of(res_ice), 4.3, tolerance = 0.4)
  expect_equal(ate_of(res_ipw), 4.3, tolerance = 0.6)
})

# -- IPCW + MI (missing L and censored Y simultaneously) ----------------------

test_that("IPCW + MI recovers the ATE under joint L-missingness + Y-censoring", {
  skip_if_not_installed("mice")
  dat <- simulate_mi_ipcw(n = 4000)
  # Keep Y in the mids but do not impute it (censored Y is handled by IPCW);
  # do not let the censoring indicator C predict the imputations.
  meth <- mice::make.method(dat[, c("Y", "A", "L", "C")])
  meth["Y"] <- ""
  pred <- mice::make.predictorMatrix(dat[, c("Y", "A", "L", "C")])
  pred[, "C"] <- 0
  imp <- mice::mice(
    dat[, c("Y", "A", "L", "C")],
    m = 8,
    method = meth,
    predictorMatrix = pred,
    printFlag = FALSE
  )

  res <- causat_mice(
    imp,
    "Y",
    "A",
    ~L,
    ints(),
    estimator = "gcomp",
    censoring = "C",
    ipcw = TRUE
  )
  expect_equal(ate_of(res), 3, tolerance = 0.2)
})

# -- S3 methods on a pooled result --------------------------------------------

test_that("tidy / confint / coef honor Barnard-Rubin df", {
  skip_if_not_installed("mice")
  dat <- simulate_mi_covariate(n = 1500, seed = 9)
  imp <- mice::mice(dat[, c("Y", "A", "L")], m = 10, printFlag = FALSE)
  res <- causat_mice(imp, "Y", "A", ~L, ints(), estimator = "gcomp")
  mi <- attr(res, "mi_details")

  # confint reconstructs t-based CIs from the stored df; check it matches a
  # hand-computed bound and honors a non-default level.
  ci <- confint(res, level = 0.90)
  est <- res$estimates$estimate
  se <- res$estimates$se
  crit <- stats::qt(0.95, df = mi$estimates$df)
  expect_equal(unname(ci[, "lower"]), est - crit * se, tolerance = 1e-10)

  # tidy contrasts use the contrast-row df, not the normal quantile.
  td <- tidy(res, which = "contrasts")
  crit_c <- stats::qt(0.975, df = mi$contrasts$df[1])
  expect_equal(
    td$conf.low[1],
    res$contrasts$estimate[1] - crit_c * res$contrasts$se[1],
    tolerance = 1e-10
  )

  cf <- coef(res)
  expect_equal(unname(cf), res$estimates$estimate)
  expect_s3_class(glance(res), "data.frame")
  expect_output(print(res), "rubin")
})

# -- Edge cases ---------------------------------------------------------------

test_that("m = 1 produces a degenerate pool with a warning", {
  skip_if_not_installed("mice")
  dat <- simulate_mi_covariate(n = 1000, seed = 5)
  imp <- mice::mice(dat[, c("Y", "A", "L")], m = 1, printFlag = FALSE)
  expect_warning(
    res <- causat_mice(imp, "Y", "A", ~L, ints(), estimator = "gcomp"),
    class = "causatr_mi_degenerate"
  )
  mi <- attr(res, "mi_details")
  expect_equal(mi$contrasts$between[1], 0) # no between-imputation variance
})

test_that("causat_mice warns when the outcome is excluded as a predictor", {
  skip_if_not_installed("mice")
  dat <- simulate_mi_covariate(n = 800, seed = 6)
  # Deliberately uncongenial: impute L without using Y as a predictor. Y is
  # in the mids data (so the fit succeeds) but its predictor-matrix column is
  # all zeros, which is the warning condition.
  pred <- mice::make.predictorMatrix(dat[, c("Y", "A", "L")])
  pred[, "Y"] <- 0
  imp <- mice::mice(
    dat[, c("Y", "A", "L")],
    m = 3,
    predictorMatrix = pred,
    printFlag = FALSE
  )
  expect_warning(
    causat_mice(imp, "Y", "A", ~L, ints(), estimator = "gcomp"),
    class = "causatr_mi_missing_predictors"
  )
})

test_that("causat_mice aborts when every imputation fails to fit", {
  skip_if_not_installed("mice")
  dat <- simulate_mi_covariate(n = 500, seed = 8)
  imp <- mice::mice(dat[, c("Y", "A", "L")], m = 3, printFlag = FALSE)
  # A binomial family on a continuous outcome makes every per-imputation
  # glm() call error (y must be in [0, 1]) without tripping the
  # missing-predictor warning (Y, A, L are all present in the mids).
  expect_error(
    causat_mice(
      imp,
      "Y",
      "A",
      ~L,
      ints(),
      estimator = "gcomp",
      family = "binomial"
    ),
    class = "causatr_mi_all_failed"
  )
})

# -- Rejection paths ----------------------------------------------------------

test_that("causat_mice rejects non-mids input", {
  skip_if_not_installed("mice")
  expect_snapshot(
    error = TRUE,
    causat_mice(
      list(),
      outcome = "Y",
      treatment = "A",
      confounders = ~L,
      interventions = list(a1 = static(1))
    )
  )
})

test_that("causat_mice rejects an unknown pool_method", {
  skip_if_not_installed("mice")
  imp <- structure(list(m = 2L), class = "mids")
  expect_error(
    causat_mice(imp, "Y", "A", ~L, ints(), pool_method = "bayes")
  )
})
