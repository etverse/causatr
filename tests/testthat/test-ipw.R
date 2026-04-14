test_that("causat(method = 'ipw') fits on simple data", {
  df <- data.frame(
    Y = c(1, 2, 3, 4, 5, 6, 7, 8),
    A = c(0, 0, 0, 0, 1, 1, 1, 1),
    L = c(1, 2, 3, 4, 1, 2, 3, 4)
  )
  fit <- suppressWarnings(causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    method = "ipw"
  ))
  expect_s3_class(fit, "causatr_fit")
  expect_equal(fit$method, "ipw")
})

test_that("IPW with external weights stores in details", {
  d <- simulate_binary_continuous(n = 500, seed = 42)
  w <- runif(nrow(d), 0.5, 2)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    method = "ipw",
    weights = w
  )
  expect_true(!is.null(fit$details$weights))
  expect_equal(length(fit$details$weights), nrow(d))
  expect_false(".causatr_w" %in% names(fit$data))
})

test_that("IPW bootstrap with external weights gives finite SE", {
  d <- simulate_binary_continuous(n = 500, seed = 42)
  w <- runif(nrow(d), 0.5, 2)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    method = "ipw",
    weights = w
  )
  res <- suppressWarnings(contrast(
    fit,
    interventions = list(a0 = static(0), a1 = static(1)),
    ci_method = "bootstrap",
    n_boot = 50
  ))
  expect_true(all(is.finite(res$contrasts$se)))
  expect_true(res$contrasts$se > 0)
})

test_that("IPW continuous treatment via WeightIt fits", {
  d <- simulate_continuous_continuous(n = 2000, seed = 42)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    method = "ipw"
  )
  expect_s3_class(fit, "causatr_fit")
  expect_equal(fit$method, "ipw")
})

test_that("causat(method = 'ipw') rejects longitudinal data", {
  df <- data.frame(
    Y = c(0, 1, 0, 1),
    A = c(0, 1, 0, 1),
    L = c(1, 1, 2, 2),
    id = c(1, 1, 2, 2),
    time = c(0, 1, 0, 1)
  )
  expect_error(
    causat(
      df,
      outcome = "Y",
      treatment = "A",
      confounders = ~L,
      method = "ipw",
      id = "id",
      time = "time"
    ),
    "longitudinal"
  )
})


test_that("IPW × non-Mparts WeightIt method × sandwich returns finite SE at contrast time", {
  # Coverage gap from FEATURE_COVERAGE_MATRIX: non-Mparts WeightIt
  # methods (gbm, super, bart, optweight, energy, npcbps) drop the
  # propensity-uncertainty correction from the sandwich SE. The
  # fit-time warning was already tested elsewhere; this asserts the
  # contrast-time behavior — sandwich still returns a finite
  # (understated) SE, and bootstrap is the recommended alternative.
  #
  # We avoid depending on the extra packages those methods need
  # (gbm, dbarts, osqp, CBPS) by stubbing `WeightIt::weightit` to
  # strip the Mparts attribute from a plain glm weight object. This
  # exercises the same code path that runs when a user picks an
  # actual non-Mparts method on their own machine.
  skip_if_not_installed("WeightIt")
  set.seed(314)
  n <- 600
  L <- stats::rnorm(n)
  A <- stats::rbinom(n, 1, stats::plogis(0.4 * L))
  Y <- 1 + 0.5 * A + 0.3 * L + stats::rnorm(n)
  df <- data.frame(Y = Y, A = A, L = L)

  real_w <- WeightIt::weightit(
    A ~ L,
    data = df,
    method = "glm",
    estimand = "ATE"
  )
  attr(real_w, "Mparts") <- NULL
  real_w$method <- "gbm"

  # Stub `WeightIt::weightit` for this test only. `local_mocked_bindings`
  # does not dispatch on namespace-qualified calls, so use
  # `assignInNamespace` + on.exit to restore the original binding.
  ns_orig <- WeightIt::weightit
  on.exit(
    assignInNamespace("weightit", ns_orig, ns = "WeightIt"),
    add = TRUE
  )
  assignInNamespace(
    "weightit",
    function(...) real_w,
    ns = "WeightIt"
  )

  expect_warning(
    fit <- causat(
      df,
      outcome = "Y",
      treatment = "A",
      confounders = ~L,
      method = "ipw"
    ),
    "does not implement the M-estimation correction"
  )

  res_sand <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    ci_method = "sandwich"
  )
  expect_true(all(is.finite(res_sand$estimates$se)))
  expect_true(all(res_sand$estimates$se > 0))
  expect_true(all(is.finite(res_sand$contrasts$se)))
  expect_true(all(res_sand$contrasts$se > 0))

  res_boot <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    ci_method = "bootstrap",
    n_boot = 200L
  )
  expect_true(all(is.finite(res_boot$contrasts$se)))
  expect_true(all(res_boot$contrasts$se > 0))
})


test_that("IPW rejects shift intervention until Phase 4", {
  # Phase-3 IPW only handles static interventions (binary/categorical).
  # Shift / MTP / IPSI / dynamic require self-contained IPW (Phase 4)
  # because WeightIt does not return density-ratio weights for
  # non-static interventions. Check the error path is not silent.
  data("nhefs", package = "causatr")
  fit <- causat(
    nhefs,
    outcome = "wt82_71",
    treatment = "qsmk",
    confounders = ~ sex + age + wt71,
    method = "ipw",
    censoring = "censored"
  )
  expect_snapshot(
    error = TRUE,
    contrast(
      fit,
      interventions = list(s = shift(1), c = static(0)),
      ci_method = "sandwich"
    )
  )
})
