# Tests for the censoring model primitive (Phase 14a).
#
# Validates:
#   - fit_censoring_model() coefficient recovery on DGP-M2
#   - compute_ipcw_weights() stabilization properties
#   - fit_censoring_models_longitudinal() on DGP-M5
#   - check_ipcw_inputs() validation
#   - make_ipcw_weight_fn() closure correctness
#   - refit_censoring_weights() bootstrap replay


# ── check_ipcw_inputs() ────────────────────────────────────────────

test_that("ipcw=TRUE without censoring aborts", {
  expect_error(
    check_ipcw_inputs(
      ipcw = TRUE,
      censoring = NULL,
      call = rlang::current_env()
    ),
    class = "causatr_ipcw_no_censoring"
  )
})

test_that("ipcw=FALSE passes silently", {
  expect_silent(
    check_ipcw_inputs(
      ipcw = FALSE,
      censoring = NULL
    )
  )
})

test_that("non-binary censoring column aborts", {
  expect_error(
    check_ipcw_inputs(
      ipcw = TRUE,
      censoring = "C",
      censoring_col = c(0, 1, 2, 0)
    ),
    class = "causatr_ipcw_non_binary"
  )
})

test_that("high censoring rate warns", {
  cens_col <- c(rep(1L, 90), rep(0L, 10))
  expect_warning(
    check_ipcw_inputs(
      ipcw = TRUE,
      censoring = "C",
      censoring_col = cens_col
    ),
    class = "causatr_ipcw_high_censoring"
  )
})

test_that("valid binary censoring passes", {
  cens_col <- c(rep(0L, 70), rep(1L, 30))
  expect_silent(
    check_ipcw_inputs(
      ipcw = TRUE,
      censoring = "C",
      censoring_col = cens_col
    )
  )
})


# ── fit_censoring_model() ──────────────────────────────────────────

test_that("censoring model recovers DGP-M2 coefficients", {
  d <- simulate_mar_outcome(n = 5000, seed = 142)
  dt <- data.table::as.data.table(d)

  cm <- fit_censoring_model(
    data = dt,
    censoring = "C",
    treatment = "A",
    confounders = ~L,
    model_fn = stats::glm
  )

  expect_s3_class(cm, "causatr_censoring_model")

  # DGP-M2 censoring model: C ~ expit(-0.5 + 1.5*A + 1.0*L)
  # So the uncensoring model (1-C) ~ expit(0.5 - 1.5*A - 1.0*L)
  coefs <- cm$alpha_hat
  expect_equal(coefs[["(Intercept)"]], 0.5, tolerance = 0.2)
  expect_equal(coefs[["A"]], -1.5, tolerance = 0.3)
  expect_equal(coefs[["L"]], -1.0, tolerance = 0.2)
})

test_that("censoring model fit_rows excludes NA censoring", {
  dt <- data.table::data.table(
    Y = c(1, 2, 3, 4, 5),
    A = c(1, 0, 1, 0, 1),
    L = c(0.5, -0.3, 0.1, 0.8, -0.2),
    C = c(0L, 1L, 0L, NA_integer_, 0L)
  )
  cm <- fit_censoring_model(
    data = dt,
    censoring = "C",
    treatment = "A",
    confounders = ~L
  )
  expect_equal(sum(cm$fit_rows), 4L)
  expect_false(cm$fit_rows[4])
})

test_that("censoring model with all uncensored works", {
  dt <- data.table::data.table(
    Y = rnorm(100),
    A = rbinom(100, 1, 0.5),
    L = rnorm(100),
    C = rep(0L, 100)
  )
  cm <- fit_censoring_model(
    data = dt,
    censoring = "C",
    treatment = "A",
    confounders = ~L
  )
  # All P(uncensored) should be 1 (or close)
  expect_true(all(cm$p_uncensored > 0.99))
  expect_equal(cm$p_marginal, 1)
})


# ── compute_ipcw_weights() ─────────────────────────────────────────

test_that("stabilized IPCW weights have mean ~1 among uncensored", {
  d <- simulate_mar_outcome(n = 5000, seed = 143)
  dt <- data.table::as.data.table(d)

  cm <- fit_censoring_model(
    data = dt,
    censoring = "C",
    treatment = "A",
    confounders = ~L
  )
  w <- compute_ipcw_weights(
    cm,
    n_total = nrow(dt),
    censoring_col = dt$C,
    stabilize = TRUE
  )

  # Weights should have length == nrow(data)
  expect_length(w, nrow(dt))

  # Censored rows should have weight 0
  expect_true(all(w[dt$C == 1L] == 0))

  # Uncensored rows should have positive weights
  expect_true(all(w[dt$C == 0L] > 0))

  # Mean of stabilized weights among uncensored should be ~1
  w_uncens <- w[dt$C == 0L]
  expect_equal(mean(w_uncens), 1.0, tolerance = 0.05)
})

test_that("unstabilized weights are >= 1 for uncensored", {
  d <- simulate_mar_outcome(n = 2000, seed = 144)
  dt <- data.table::as.data.table(d)

  cm <- fit_censoring_model(
    data = dt,
    censoring = "C",
    treatment = "A",
    confounders = ~L
  )
  w <- compute_ipcw_weights(
    cm,
    n_total = nrow(dt),
    censoring_col = dt$C,
    stabilize = FALSE
  )

  w_uncens <- w[dt$C == 0L]
  expect_true(all(w_uncens >= 1 - 1e-10))
})

test_that("IPCW weights with no censoring are all 1", {
  dt <- data.table::data.table(
    Y = rnorm(100),
    A = rbinom(100, 1, 0.5),
    L = rnorm(100),
    C = rep(0L, 100)
  )
  cm <- fit_censoring_model(
    data = dt,
    censoring = "C",
    treatment = "A",
    confounders = ~L
  )
  w <- compute_ipcw_weights(
    cm,
    n_total = 100,
    censoring_col = dt$C,
    stabilize = TRUE
  )
  expect_equal(w, rep(1, 100), tolerance = 1e-6)
})


# ── make_ipcw_weight_fn() ─────────────────────────────────────────

test_that("weight closure reproduces fitted weights at alpha_hat", {
  d <- simulate_mar_outcome(n = 2000, seed = 145)
  dt <- data.table::as.data.table(d)

  cm <- fit_censoring_model(
    data = dt,
    censoring = "C",
    treatment = "A",
    confounders = ~L
  )

  wfn <- make_ipcw_weight_fn(
    cm,
    n_total = nrow(dt),
    censoring_col = dt$C,
    stabilize = TRUE
  )

  w_from_fn <- wfn(cm$alpha_hat)
  w_direct <- compute_ipcw_weights(
    cm,
    n_total = nrow(dt),
    censoring_col = dt$C,
    stabilize = TRUE
  )

  expect_equal(w_from_fn, w_direct, tolerance = 1e-10)
})

test_that("weight closure responds to parameter perturbation", {
  d <- simulate_mar_outcome(n = 1000, seed = 146)
  dt <- data.table::as.data.table(d)

  cm <- fit_censoring_model(
    data = dt,
    censoring = "C",
    treatment = "A",
    confounders = ~L
  )

  wfn <- make_ipcw_weight_fn(
    cm,
    n_total = nrow(dt),
    censoring_col = dt$C,
    stabilize = TRUE
  )

  w0 <- wfn(cm$alpha_hat)
  # Perturb intercept
  gamma_perturbed <- cm$alpha_hat
  gamma_perturbed[1] <- gamma_perturbed[1] + 0.5
  w1 <- wfn(gamma_perturbed)

  # Weights should change
  expect_false(isTRUE(all.equal(w0, w1)))
  # But censored rows should still be 0
  expect_true(all(w1[dt$C == 1L] == 0))
})


# ── fit_censoring_models_longitudinal() ────────────────────────────

test_that("longitudinal censoring models fit on DGP-M5", {
  d <- simulate_longitudinal_mar_outcome(n = 3000, seed = 147)
  dt <- data.table::as.data.table(d)

  # prepare_data() adds lag columns; must include censoring
  # so the C column is preserved
  dt <- prepare_data(
    dt,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    history = 1L,
    censoring = "C",
    cluster = NULL
  )

  time_points <- sort(unique(dt$time))

  result <- fit_censoring_models_longitudinal(
    data = dt,
    censoring = "C",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    model_fn = stats::glm,
    id = "id",
    time = "time",
    time_points = time_points,
    history = 1L
  )

  expect_type(result, "list")
  expect_named(
    result,
    c("models", "cumulative_weights", "per_period_weights")
  )

  # Cumulative weights should be length nrow(data)
  expect_length(result$cumulative_weights, nrow(dt))

  # At time 0, no censoring -> all weights should be 1
  t0_rows <- dt$time == 0
  expect_equal(
    result$per_period_weights[["0"]],
    rep(1, sum(t0_rows))
  )

  # At time 1, there IS censoring -> model should be non-NULL
  expect_s3_class(
    result$models[["1"]],
    "causatr_censoring_model"
  )
})

test_that("longitudinal IPCW cumulative weights are zero for censored", {
  d <- simulate_longitudinal_mar_outcome(n = 2000, seed = 148)
  dt <- data.table::as.data.table(d)

  dt <- prepare_data(
    dt,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    history = 1L,
    censoring = "C",
    cluster = NULL
  )

  time_points <- sort(unique(dt$time))

  result <- fit_censoring_models_longitudinal(
    data = dt,
    censoring = "C",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    model_fn = stats::glm,
    id = "id",
    time = "time",
    time_points = time_points,
    history = 1L
  )

  # Censored rows (C=1) should have cumulative weight = 0
  censored_rows <- dt$C == 1L
  expect_true(all(result$cumulative_weights[censored_rows] == 0))

  # Uncensored rows should have positive weights
  uncensored_rows <- dt$C == 0L
  expect_true(
    all(result$cumulative_weights[uncensored_rows] > 0)
  )
})


# ── print method ───────────────────────────────────────────────────

test_that("print.causatr_censoring_model produces output", {
  d <- simulate_mar_outcome(n = 500, seed = 149)
  dt <- data.table::as.data.table(d)

  cm <- fit_censoring_model(
    data = dt,
    censoring = "C",
    treatment = "A",
    confounders = ~L
  )

  expect_output(print(cm), "causatr_censoring_model")
  expect_output(print(cm), "Formula")
  expect_output(print(cm), "Marginal P")
})
