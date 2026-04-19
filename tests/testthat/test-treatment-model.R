# Unit tests for R/treatment_model.R (Phase 4 foundation layer).
#
# The treatment density model + density evaluator underpin the
# self-contained IPW engine. These tests pin down:
#   - family classification (binary vs continuous vs categorical)
#   - the fit pipeline end-to-end against `stats::glm` internals
#   - Bernoulli / Gaussian density evaluation against hand formulas
#   - the "categorical not yet" abort
#   - sigma extraction across GLM / LM fit shapes
#
# The categorical branch is explicitly deferred in this first Phase 4
# chunk, and the test locks in the abort message so a future PR that
# lands the multinomial branch has to remove this test alongside the
# abort.

test_that("detect_treatment_family() classifies each input type", {
  expect_identical(detect_treatment_family(c(0L, 1L, 0L, 1L)), "bernoulli")
  expect_identical(detect_treatment_family(c(0, 1, 0, 1)), "bernoulli")
  expect_identical(detect_treatment_family(c(TRUE, FALSE, TRUE)), "bernoulli")
  expect_identical(
    detect_treatment_family(rnorm(20)),
    "gaussian"
  )
  expect_identical(
    detect_treatment_family(factor(c("a", "b", "c", "a"))),
    "categorical"
  )
  expect_identical(
    detect_treatment_family(c("a", "b", "c", "a")),
    "categorical"
  )
})

test_that("detect_treatment_family() rejects unsupported types", {
  expect_error(
    detect_treatment_family(as.Date("2020-01-01") + 0:3),
    "is not supported"
  )
})

test_that("fit_treatment_model() on binary treatment returns a Bernoulli fit", {
  d <- simulate_binary_continuous(n = 500, seed = 1)
  dt <- data.table::as.data.table(d)

  tm <- fit_treatment_model(
    dt,
    treatment = "A",
    confounders = ~L,
    model_fn = stats::glm
  )

  expect_s3_class(tm, "causatr_treatment_model")
  expect_identical(tm$family, "bernoulli")
  expect_identical(tm$treatment, "A")
  expect_null(tm$sigma)
  expect_equal(sum(tm$fit_rows), nrow(d))
  # alpha_hat length matches propensity design matrix width
  expect_equal(length(tm$alpha_hat), ncol(tm$X_prop))
  # X_prop should be the fitted model's design matrix verbatim
  ref_fit <- stats::glm(A ~ L, data = d, family = stats::binomial())
  expect_equal(
    unname(tm$alpha_hat),
    unname(stats::coef(ref_fit)),
    tolerance = 1e-10
  )
})

test_that("fit_treatment_model() on continuous treatment recovers sigma", {
  d <- simulate_continuous_continuous(n = 500, seed = 2)
  dt <- data.table::as.data.table(d)

  tm <- fit_treatment_model(dt, treatment = "A", confounders = ~L)

  expect_identical(tm$family, "gaussian")
  expect_true(is.numeric(tm$sigma) && tm$sigma > 0)

  # Compare sigma and alpha against a plain glm() fit on the same data.
  # Note: `summary(glm)$sigma` is NULL for a Gaussian GLM — the
  # residual-variance estimator lives on `$dispersion`. Take sqrt to
  # get the residual SD.
  ref_fit <- stats::glm(A ~ L, data = d, family = stats::gaussian())
  expect_equal(
    tm$sigma,
    sqrt(summary(ref_fit)$dispersion),
    tolerance = 1e-10
  )
  expect_equal(
    unname(tm$alpha_hat),
    unname(stats::coef(ref_fit)),
    tolerance = 1e-10
  )
})

test_that("fit_treatment_model() on categorical treatment returns a multinomial fit", {
  d <- simulate_categorical_continuous(n = 500, seed = 10)
  dt <- data.table::as.data.table(d)

  tm <- fit_treatment_model(
    dt,
    treatment = "A",
    confounders = ~L,
    model_fn = nnet::multinom
  )

  expect_s3_class(tm, "causatr_treatment_model")
  expect_identical(tm$family, "categorical")
  expect_identical(tm$treatment, "A")
  expect_null(tm$sigma)
  expect_identical(tm$levels, c("a", "b", "c"))
  expect_equal(sum(tm$fit_rows), nrow(d))
  # alpha_hat is flattened (K-1) * p = 2 * 2 = 4 for 3-level, 2-col design
  expect_equal(length(tm$alpha_hat), 4L)
})

test_that("fit_treatment_model() on 2-level factor returns categorical", {
  set.seed(11)
  dt <- data.table::data.table(
    A = factor(sample(c("ctrl", "trt"), 200, replace = TRUE)),
    L = rnorm(200)
  )

  tm <- fit_treatment_model(
    dt,
    treatment = "A",
    confounders = ~L,
    model_fn = nnet::multinom
  )

  expect_identical(tm$family, "categorical")
  expect_identical(tm$levels, c("ctrl", "trt"))
  # 2-level factor: (K-1) * p = 1 * 2 = 2 params
  expect_equal(length(tm$alpha_hat), 2L)
})

test_that("fit_treatment_model() drops rows with NA treatment / NA confounders", {
  d <- simulate_binary_continuous(n = 200, seed = 3)
  d$A[1:5] <- NA
  d$L[10:12] <- NA
  dt <- data.table::as.data.table(d)

  tm <- fit_treatment_model(dt, treatment = "A", confounders = ~L)
  # 5 NA in A + 3 NA in L (non-overlapping) = 8 dropped rows
  expect_equal(sum(tm$fit_rows), nrow(d) - 8L)
})

test_that("evaluate_density() binary returns Bernoulli pmf per row", {
  d <- simulate_binary_continuous(n = 300, seed = 4)
  dt <- data.table::as.data.table(d)
  tm <- fit_treatment_model(dt, treatment = "A", confounders = ~L)

  fit_data <- dt[tm$fit_rows]
  a_obs <- fit_data$A
  p <- stats::predict(tm$model, newdata = fit_data, type = "response")

  # Evaluate at the observed treatment: ifelse(a == 1, p, 1 - p)
  f_obs <- evaluate_density(tm, a_obs, fit_data)
  expect_equal(f_obs, ifelse(a_obs == 1, p, 1 - p), tolerance = 1e-12)

  # Evaluate at a constant 1: f_int is p for every row
  f_int_1 <- evaluate_density(tm, rep(1, nrow(fit_data)), fit_data)
  expect_equal(f_int_1, unname(p), tolerance = 1e-12)

  # Evaluate at a constant 0: f_int is 1 - p
  f_int_0 <- evaluate_density(tm, rep(0, nrow(fit_data)), fit_data)
  expect_equal(f_int_0, unname(1 - p), tolerance = 1e-12)
})

test_that("evaluate_density() continuous returns Gaussian pdf per row", {
  d <- simulate_continuous_continuous(n = 300, seed = 5)
  dt <- data.table::as.data.table(d)
  tm <- fit_treatment_model(dt, treatment = "A", confounders = ~L)

  fit_data <- dt[tm$fit_rows]
  a_obs <- fit_data$A
  mu <- stats::predict(tm$model, newdata = fit_data, type = "response")

  f_obs <- evaluate_density(tm, a_obs, fit_data)
  expect_equal(
    f_obs,
    stats::dnorm(a_obs, mean = mu, sd = tm$sigma),
    tolerance = 1e-12
  )

  # Density at a shifted treatment matches dnorm at the shifted point
  f_shift <- evaluate_density(tm, a_obs - 0.5, fit_data)
  expect_equal(
    f_shift,
    stats::dnorm(a_obs - 0.5, mean = mu, sd = tm$sigma),
    tolerance = 1e-12
  )
})

test_that("evaluate_density() categorical returns multinomial pmf per row", {
  d <- simulate_categorical_continuous(n = 500, seed = 12)
  dt <- data.table::as.data.table(d)
  tm <- fit_treatment_model(
    dt,
    treatment = "A",
    confounders = ~L,
    model_fn = nnet::multinom
  )

  fit_data <- dt[tm$fit_rows]
  a_obs <- fit_data$A

  # Density at observed treatment should match the predict(type="probs")
  # column for each row's observed level.
  probs_raw <- predict(tm$model, newdata = fit_data, type = "probs")
  a_char <- as.character(a_obs)
  expected <- vapply(
    seq_len(nrow(fit_data)),
    function(i) probs_raw[i, a_char[i]],
    numeric(1)
  )

  f_obs <- evaluate_density(tm, a_obs, fit_data)
  expect_equal(f_obs, expected, tolerance = 1e-12)

  # Density at a fixed level "b" for all rows should be the "b" column
  f_b <- evaluate_density(tm, rep("b", nrow(fit_data)), fit_data)
  expect_equal(f_b, unname(probs_raw[, "b"]), tolerance = 1e-12)

  # All probabilities sum to ~1 per row
  f_a <- evaluate_density(tm, rep("a", nrow(fit_data)), fit_data)
  f_c <- evaluate_density(tm, rep("c", nrow(fit_data)), fit_data)
  expect_equal(f_a + f_b + f_c, rep(1, nrow(fit_data)), tolerance = 1e-12)
})

test_that("evaluate_density() rejects length mismatch", {
  d <- simulate_binary_continuous(n = 100, seed = 6)
  dt <- data.table::as.data.table(d)
  tm <- fit_treatment_model(dt, treatment = "A", confounders = ~L)

  fit_data <- dt[tm$fit_rows]
  expect_error(
    evaluate_density(tm, rep(1, 3), fit_data),
    "length"
  )
})

test_that("extract_sigma() recovers sigma from a Gaussian GLM via dispersion", {
  set.seed(7)
  df <- data.frame(A = rnorm(100), L = rnorm(100))
  glm_fit <- stats::glm(A ~ L, data = df, family = stats::gaussian())
  # `summary(glm)$sigma` is NULL on a Gaussian GLM (common gotcha);
  # `summary(glm)$dispersion` is the estimator of sigma^2.
  expect_equal(
    extract_sigma(glm_fit),
    sqrt(summary(glm_fit)$dispersion),
    tolerance = 1e-12
  )
})

test_that("extract_sigma() recovers sigma from an lm fit via $sigma", {
  set.seed(8)
  df <- data.frame(A = rnorm(100), L = rnorm(100))
  lm_fit <- stats::lm(A ~ L, data = df)
  expect_equal(extract_sigma(lm_fit), summary(lm_fit)$sigma, tolerance = 1e-12)
})

# print.causatr_treatment_model() display test ---------------------------------

# The print method is purely cosmetic but is the only public window into
# the treatment-model object. The test exercises each conditional branch
# (sigma for Gaussian, theta for NB, levels for categorical) so a
# regression that drops one of those slots is caught.
# ---- extract_sigma() residual fallbacks ----------------------------------

test_that("fit_treatment_model() aborts when negbin model is missing theta", {
  set.seed(106)
  d <- data.table::as.data.table(data.frame(
    A = rnbinom(200, mu = 3, size = 1),
    L = rnorm(200)
  ))
  # Custom propensity_model_fn that returns a glm-like object without
  # the `$theta` slot. The branch fires at line ~222 of treatment_model.R.
  bad_nb_fn <- function(formula, data, ...) {
    fit <- stats::glm(formula, data = data, family = stats::poisson())
    fit$theta <- NULL # explicitly strip theta
    fit
  }
  expect_error(
    fit_treatment_model(
      d,
      treatment = "A",
      confounders = ~L,
      model_fn = bad_nb_fn,
      propensity_family = "negbin"
    ),
    "Could not extract dispersion parameter"
  )
})

test_that("fit_treatment_model() aborts on a coefficient/design-matrix dimension mismatch", {
  # Force a mismatch by stubbing `model_fn` to return a fit whose
  # `coef()` length doesn't match the design matrix's column count
  # (the canonical symptom of an aliased / dropped column).
  set.seed(107)
  d <- data.table::as.data.table(simulate_binary_continuous(
    n = 100,
    seed = 107
  ))
  bad_fn <- function(formula, data, family, ...) {
    fit <- stats::glm(formula, data = data, family = family)
    # Splice in an extra fake coefficient so length(coef) > ncol(model.matrix).
    fit$coefficients <- c(fit$coefficients, fake = 1.23)
    fit
  }
  expect_error(
    fit_treatment_model(
      d,
      treatment = "A",
      confounders = ~L,
      model_fn = bad_fn
    ),
    "this usually means a column was aliased"
  )
})

# ---- fit_count_density() weighted path -----------------------------------

test_that("fit_count_density() forwards `weights` to the Poisson fitter", {
  # The Poisson branch (line ~492) propagates the resampled weight
  # vector via `base_args$weights <- weights_fit`. NB has the same
  # plumbing at line ~505 but `MASS::glm.nb` is brittle to non-uniform
  # weights and no production code path inside causatr exercises it
  # today, so we test only the Poisson branch.
  set.seed(108)
  n <- 200
  d <- data.table::as.data.table(data.frame(
    A = rpois(n, lambda = 3),
    L = rnorm(n)
  ))
  w <- runif(n, 0.5, 1.5)
  tm_pois <- fit_treatment_model(
    d,
    treatment = "A",
    confounders = ~L,
    propensity_family = "poisson",
    weights = w
  )
  expect_s3_class(tm_pois, "causatr_treatment_model")
  expect_identical(tm_pois$family, "poisson")
})

test_that("extract_sigma() reads `$sig2` from a mgcv::gam fit", {
  skip_if_not_installed("mgcv")
  set.seed(101)
  df <- data.frame(A = rnorm(200), L = rnorm(200))
  gam_fit <- mgcv::gam(A ~ s(L), data = df)
  # mgcv stores the residual variance as `$sig2`. summary.gam() returns
  # an object with neither `$sigma` nor `$dispersion`, so this branch
  # is the only one that fires for GAM fits.
  expect_equal(extract_sigma(gam_fit), sqrt(gam_fit$sig2), tolerance = 1e-12)
})

test_that("extract_sigma() falls back to response residuals when summary lacks sigma/dispersion", {
  # The fallback branch fires when `summary()` returns no `$sigma`
  # and no `$dispersion`, and `model$sig2` is NULL: it computes
  # `sqrt(sum(r^2) / (n - p))` directly. We exercise it with a custom
  # S3 class whose summary/residuals/coef methods are defined inside
  # this test_that block, so registration is local to the call frame
  # via `withr::defer()` cleanup.
  set.seed(102)
  r <- rnorm(50)
  cls <- "causatr_test_fakefit"
  registerS3method(
    "summary",
    cls,
    function(object, ...) list(),
    envir = .GlobalEnv
  )
  registerS3method(
    "residuals",
    cls,
    function(object, type = "response", ...) r,
    envir = .GlobalEnv
  )
  registerS3method(
    "coef",
    cls,
    function(object, ...) c(intercept = 0, slope = 0),
    envir = .GlobalEnv
  )
  withr::defer({
    s3 <- get(".__S3MethodsTable__.", envir = .GlobalEnv)
    for (nm in c("summary.", "residuals.", "coef.")) {
      key <- paste0(nm, cls)
      if (exists(key, envir = s3, inherits = FALSE)) {
        rm(list = key, envir = s3)
      }
    }
  })

  fake_fit <- structure(list(), class = cls)
  expect_equal(
    extract_sigma(fake_fit),
    sqrt(sum(r^2) / (length(r) - 2L)),
    tolerance = 1e-12
  )
})

test_that("extract_sigma() aborts when neither summary nor residuals work", {
  # All sigma extraction paths fail: `summary()` errors, and the
  # `residuals()` fallback also errors. The function aborts with a
  # clear message so a downstream `dnorm()` does not silently get
  # `sigma = NA`.
  cls <- "causatr_test_brokenfit"
  registerS3method(
    "summary",
    cls,
    function(object, ...) stop("no summary"),
    envir = .GlobalEnv
  )
  registerS3method(
    "residuals",
    cls,
    function(object, ...) stop("no residuals"),
    envir = .GlobalEnv
  )
  withr::defer({
    s3 <- get(".__S3MethodsTable__.", envir = .GlobalEnv)
    for (nm in c("summary.", "residuals.")) {
      key <- paste0(nm, cls)
      if (exists(key, envir = s3, inherits = FALSE)) {
        rm(list = key, envir = s3)
      }
    }
  })

  fake_fit <- structure(list(), class = cls)
  expect_error(extract_sigma(fake_fit), "Could not extract residual SD")
})

# ---- evaluate_density() defensive guards ---------------------------------

test_that("evaluate_density() rejects a non-treatment_model `treatment_model`", {
  expect_error(
    evaluate_density(
      list(family = "bernoulli"),
      treatment_values = c(0, 1),
      newdata = data.frame(L = c(0, 1))
    ),
    "causatr_treatment_model"
  )
})

test_that("evaluate_density() aborts on an unknown family tag", {
  d <- data.table::as.data.table(simulate_binary_continuous(n = 50, seed = 103))
  tm <- fit_treatment_model(d, treatment = "A", confounders = ~L)
  evil_tm <- tm
  evil_tm$family <- "frobnitz"
  expect_error(
    evaluate_density(
      evil_tm,
      treatment_values = d$A,
      newdata = d
    ),
    "Unknown treatment family"
  )
})

# ---- evaluate_categorical_density() K=2 + unknown-level paths ------------

test_that("evaluate_categorical_density() handles K=2 (vector predict return)", {
  # nnet::multinom returns a vector (not a matrix) for K = 2 -- the
  # function builds the n x 2 probability matrix from `cbind(1 - p, p)`.
  # This regression test pins that branch.
  set.seed(104)
  d <- data.frame(
    A = factor(sample(c("a", "b"), 200, replace = TRUE)),
    L = rnorm(200)
  )
  m <- nnet::multinom(A ~ L, data = d, trace = FALSE)
  out <- evaluate_categorical_density(
    m,
    trt_levels = c("a", "b"),
    treatment_values = d$A,
    newdata = d
  )
  expect_length(out, nrow(d))
  expect_true(all(out > 0 & out < 1))
})

test_that("evaluate_categorical_density() aborts on a treatment value not in model levels", {
  set.seed(105)
  d <- data.frame(
    A = factor(sample(c("a", "b", "c"), 200, replace = TRUE)),
    L = rnorm(200)
  )
  m <- nnet::multinom(A ~ L, data = d, trace = FALSE)
  expect_error(
    evaluate_categorical_density(
      m,
      trt_levels = c("a", "b", "c"),
      treatment_values = c("a", "z", "b"), # 'z' is not a known level
      newdata = d[1:3, ]
    ),
    "not found in model levels"
  )
})

# ---- print.causatr_treatment_model display test --------------------------

test_that("print.causatr_treatment_model prints expected fields", {
  # The print method is internal (`@noRd`, not registered in NAMESPACE),
  # so generic S3 dispatch via `print()` from outside the package would
  # hit `print.default`. We call the method by its full name to exercise
  # the body the user actually maintains.
  d_bin <- data.table::as.data.table(
    simulate_binary_continuous(n = 200, seed = 11)
  )
  tm_bin <- fit_treatment_model(
    d_bin,
    treatment = "A",
    confounders = ~L,
    model_fn = stats::glm
  )
  expect_s3_class(tm_bin, "causatr_treatment_model")
  out_bin <- capture.output(print.causatr_treatment_model(tm_bin))
  expect_match(out_bin[1], "causatr_treatment_model: bernoulli")
  expect_true(any(grepl("treatment:", out_bin)))
  expect_true(any(grepl("p_alpha:", out_bin)))

  d_cont <- data.table::as.data.table(
    simulate_continuous_continuous(n = 200, seed = 12)
  )
  tm_cont <- fit_treatment_model(
    d_cont,
    treatment = "A",
    confounders = ~L,
    model_fn = stats::glm
  )
  out_cont <- capture.output(print.causatr_treatment_model(tm_cont))
  expect_true(any(grepl("sigma:", out_cont)))

  set.seed(13)
  d_cat <- data.table::as.data.table(data.frame(
    A = factor(sample(letters[1:3], 200, replace = TRUE)),
    L = rnorm(200)
  ))
  tm_cat <- fit_treatment_model(
    d_cat,
    treatment = "A",
    confounders = ~L,
    model_fn = nnet::multinom
  )
  out_cat <- capture.output(print.causatr_treatment_model(tm_cat))
  expect_true(any(grepl("levels:", out_cat)))
})
