# Unit tests for `replay_fit()` — the central `c(base_args, dots)` +
# do.call helper introduced by Alternative E in the 2026-04-15
# third-round critical-review dots audit. This file pins the
# edge-case behavior so the four refit sites (gcomp, IPW, matching,
# ICE) can trust the helper without restating the invariants.

test_that("replay_fit: base_args win over duplicate dots keys", {
  f <- function(x, method = "default") method
  r <- replay_fit(
    f,
    base_args = list(x = 1, method = "base"),
    dots = list(method = "dots")
  )
  expect_equal(r, "base")
})

test_that("replay_fit: positional (unnamed) dots are dropped", {
  f <- function(x, method = "default") method
  r <- replay_fit(
    f,
    base_args = list(x = 1),
    dots = list("positional_stray", method = "kept")
  )
  expect_equal(r, "kept")
})

test_that("replay_fit: reserved keys are stripped from dots", {
  f <- function(x, method = "default") method
  r <- replay_fit(
    f,
    base_args = list(x = 1),
    dots = list(method = "blocked"),
    reserved = "method"
  )
  expect_equal(r, "default")
})

test_that("replay_fit: NULL dots is a no-op", {
  f <- function(x, method = "default") method
  r <- replay_fit(
    f,
    base_args = list(x = 1),
    dots = NULL
  )
  expect_equal(r, "default")
})

test_that("replay_fit: non-conflicting dots pass through", {
  f <- function(x, method = "default", power = 1) c(method, power)
  r <- replay_fit(
    f,
    base_args = list(x = 1),
    dots = list(method = "user", power = 2)
  )
  expect_equal(r, c("user", "2"))
})

test_that("replay_fit: empty dots list is identical to NULL", {
  f <- function(x, method = "default") method
  expect_identical(
    replay_fit(f, list(x = 1), list()),
    replay_fit(f, list(x = 1), NULL)
  )
})

test_that("replay_fit: mix of reserved + base_args strip behaves additively", {
  f <- function(x, method = "default", link = "identity") {
    paste(method, link, sep = "/")
  }
  r <- replay_fit(
    f,
    base_args = list(x = 1, method = "base"),
    dots = list(method = "dots", link = "blocked"),
    reserved = "link"
  )
  expect_equal(r, "base/identity")
})

test_that("replay_fit: unnamed entry alongside named works (name kept)", {
  # A stray positional element should not contaminate the named keys.
  f <- function(x, method = "default") method
  r <- replay_fit(
    f,
    base_args = list(x = 1),
    dots = list("stray1", method = "kept", "stray2")
  )
  expect_equal(r, "kept")
})

test_that("gcomp bootstrap replays user `...` via replay_fit (end to end)", {
  df <- simulate_binary_continuous(n = 300, seed = 11L)
  df$Y <- pmax(0, round(df$Y + 5))
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "gcomp",
    family = stats::quasipoisson()
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "ratio",
    ci_method = "bootstrap",
    n_boot = 20L
  )
  expect_s3_class(res, "causatr_result")
  expect_true(all(res$estimates$estimate > 0))
})

test_that("IPW bootstrap replays method=cbps via replay_fit (end to end)", {
  df <- simulate_binary_continuous(n = 400, seed = 12L)
  fit <- suppressWarnings(
    causat(
      df,
      outcome = "Y",
      treatment = "A",
      confounders = ~L,
      estimator = "ipw",
      method = "cbps"
    )
  )
  expect_equal(fit$details$dots$method, "cbps")
  res <- suppressWarnings(
    contrast(
      fit,
      interventions = list(a1 = static(1), a0 = static(0)),
      type = "difference",
      ci_method = "bootstrap",
      n_boot = 10L
    )
  )
  expect_s3_class(res, "causatr_result")
})

test_that("replay_fit: R3 — base R catches unknown dots at fit time", {
  # Confirms the dots audit R3 finding: stats::glm's glm.control
  # errors on unused arguments, so bogus dots are not silent.
  df <- simulate_binary_continuous(n = 100, seed = 13L)
  expect_error(
    causat(
      df,
      outcome = "Y",
      treatment = "A",
      confounders = ~L,
      estimator = "gcomp",
      bogus_never_used = 42
    ),
    regexp = "unused argument"
  )
})

test_that("replay_fit: R5 — old fit objects without $dots still bootstrap", {
  # Backward compatibility: a fit restored from a saved RDS that
  # predates the `$dots` field must not abort at bootstrap time.
  df <- simulate_binary_continuous(n = 200, seed = 14L)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "gcomp"
  )
  fit$details$dots <- NULL
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "bootstrap",
    n_boot = 10L
  )
  expect_s3_class(res, "causatr_result")
})
