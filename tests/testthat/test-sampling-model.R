test_that("check_transport_inputs: rejects matching estimator", {
  d <- data.frame(L = 1:10, S = rep(0:1, 5))
  expect_snapshot(
    check_transport_inputs("S", d$S, "target", "matching"),
    error = TRUE
  )
})

test_that("check_transport_inputs: rejects NA in target column", {
  d <- data.frame(L = 1:10, S = c(rep(0:1, 4), NA, NA))
  expect_snapshot(
    check_transport_inputs("S", d$S, "target", "gcomp"),
    error = TRUE
  )
})

test_that("check_transport_inputs: rejects non-binary target column", {
  d <- data.frame(L = 1:10, S = c(0, 1, 2, 0, 1, 2, 0, 1, 2, 0))
  expect_snapshot(
    check_transport_inputs("S", d$S, "target", "gcomp"),
    error = TRUE
  )
})

test_that("check_transport_inputs: rejects degenerate target (all S=1)", {
  d <- data.frame(L = 1:10, S = rep(1L, 10))
  expect_snapshot(
    check_transport_inputs("S", d$S, "target", "gcomp"),
    error = TRUE
  )
})

test_that("check_transport_inputs: warns on extreme selection (>95%)", {
  s_col <- c(rep(1L, 97), rep(0L, 3))
  expect_warning(
    check_transport_inputs("S", s_col, "target", "gcomp"),
    class = "causatr_transport_extreme_selection"
  )
})

test_that("check_transport_inputs: passes silently on valid binary S", {
  s_col <- rep(0:1, 50)
  expect_invisible(
    check_transport_inputs("S", s_col, "target", "gcomp")
  )
})

test_that("fit_sampling_model: returns causatr_sampling_model", {
  d <- simulate_transport(n = 200)
  d <- data.table::as.data.table(d)
  result <- fit_sampling_model(d, "S", ~L)
  expect_s3_class(result, "causatr_sampling_model")
})

test_that("fit_sampling_model: recovers known sampling coefficients", {
  # DGP-T1: P(S=1|L) = expit(-0.5 + 1.0 * L)
  # With n=5000, coefficients should be close to (-0.5, 1.0)
  d <- simulate_transport(n = 5000, seed = 1)
  d <- data.table::as.data.table(d)
  m <- fit_sampling_model(d, "S", ~L)
  expect_equal(unname(m$gamma_hat[1]), -0.5, tolerance = 0.15)
  expect_equal(unname(m$gamma_hat[2]), 1.0, tolerance = 0.15)
})

test_that("fit_sampling_model: formula has no treatment term", {
  d <- simulate_transport(n = 200)
  d <- data.table::as.data.table(d)
  m <- fit_sampling_model(d, "S", ~L)
  rhs_vars <- all.vars(m$sampling_formula)
  expect_false("A" %in% rhs_vars)
  expect_true("L" %in% rhs_vars)
})

test_that("fit_sampling_model: p_study in (0,1), p_marginal matches mean(S)", {
  d <- simulate_transport(n = 500)
  d <- data.table::as.data.table(d)
  m <- fit_sampling_model(d, "S", ~L)
  expect_true(all(m$p_study > 0 & m$p_study < 1))
  expect_equal(m$p_marginal, mean(d$S), tolerance = 1e-10)
})

test_that("fit_sampling_model: X_fit dimensions match n_fit x n_coef", {
  d <- simulate_transport(n = 300)
  d <- data.table::as.data.table(d)
  m <- fit_sampling_model(d, "S", ~L)
  n_fit <- sum(m$fit_rows)
  n_coef <- length(m$gamma_hat)
  expect_equal(nrow(m$X_fit), n_fit)
  expect_equal(ncol(m$X_fit), n_coef)
})

test_that("fit_sampling_model: fit_rows excludes only NA confounders", {
  d <- simulate_transport(n = 300)
  d$L[c(1, 5, 10)] <- NA
  d <- data.table::as.data.table(d)
  m <- fit_sampling_model(d, "S", ~L)
  expect_equal(sum(!m$fit_rows), 3L)
})

test_that("fit_sampling_model: aborts on NA in target column", {
  d <- simulate_transport(n = 100)
  d$S[1] <- NA_integer_
  d <- data.table::as.data.table(d)
  expect_error(
    fit_sampling_model(d, "S", ~L),
    class = "causatr_transport_na_target"
  )
})

test_that("print.causatr_sampling_model: produces output", {
  d <- simulate_transport(n = 200)
  d <- data.table::as.data.table(d)
  m <- fit_sampling_model(d, "S", ~L)
  expect_output(print(m), "causatr_sampling_model")
})

test_that("causat: target = 'S' stores sampling model on fit", {
  d <- simulate_transport(n = 500)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    target = "S"
  )
  expect_s3_class(fit$details$sampling_model, "causatr_sampling_model")
  expect_true(isTRUE(fit$details$transport))
  expect_equal(fit$target, "S")
  expect_equal(fit$target_subset, "target")
})

test_that("causat: target = NULL leaves fit unchanged", {
  d <- simulate_transport(n = 500)
  d_study <- d[d$S == 1, ]
  fit <- causat(
    d_study,
    outcome = "Y",
    treatment = "A",
    confounders = ~L
  )
  expect_null(fit$target)
  expect_null(fit$details$sampling_model)
  expect_null(fit$details$transport)
})

test_that("causat: target_subset = 'all' is stored correctly", {
  d <- simulate_transport(n = 500)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    target = "S",
    target_subset = "all"
  )
  expect_equal(fit$target_subset, "all")
})

test_that("causat: target column survives in fit$data", {
  d <- simulate_transport(n = 500)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    target = "S"
  )
  expect_true("S" %in% names(fit$data))
})

test_that("causat: matching + target aborts", {
  d <- simulate_transport(n = 500)
  d_study <- d[d$S == 1, ]
  d_study$S <- sample(0:1, nrow(d_study), replace = TRUE)
  expect_snapshot(
    causat(
      d_study,
      outcome = "Y",
      treatment = "A",
      confounders = ~L,
      estimator = "matching",
      target = "S",
      method = "nearest"
    ),
    error = TRUE
  )
})
