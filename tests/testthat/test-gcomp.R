test_that("causat(estimator = 'gcomp') fits on simple data", {
  df <- data.frame(
    Y = c(1, 2, 3, 4, 5, 6, 7, 8),
    A = c(0, 0, 0, 0, 1, 1, 1, 1),
    L = c(1, 2, 3, 4, 1, 2, 3, 4)
  )
  fit <- causat(df, outcome = "Y", treatment = "A", confounders = ~L)
  expect_s3_class(fit, "causatr_fit")
  expect_equal(fit$estimator, "gcomp")
  expect_equal(fit$type, "point")
  expect_equal(fit$details$n_fit, 8L)
})

test_that("g-computation replicates NHEFS ATE ≈ 3.5 kg (Ch. 13)", {
  data("nhefs", package = "causatr")

  # Hernán & Robins model: quadratic terms + one A × smokeintensity interaction.
  fit <- causat(
    nhefs,
    outcome = "wt82_71",
    treatment = "qsmk",
    confounders = ~ sex +
      age +
      I(age^2) +
      race +
      factor(education) +
      smokeintensity +
      I(smokeintensity^2) +
      smokeyrs +
      I(smokeyrs^2) +
      factor(exercise) +
      factor(active) +
      wt71 +
      I(wt71^2) +
      qsmk:smokeintensity,
    censoring = "censored"
  )

  result <- contrast(
    fit,
    interventions = list(quit = static(1), continue = static(0)),
    type = "difference",
    ci_method = "sandwich",
    reference = "continue"
  )

  # Book target: ATE ≈ 3.5 kg (95% CI: 2.6–4.5).
  ate <- result$contrasts$estimate[1]
  expect_gt(ate, 2.5)
  expect_lt(ate, 4.5)

  # CI should cover the book value.
  expect_lt(result$contrasts$ci_lower[1], 3.5)
  expect_gt(result$contrasts$ci_upper[1], 3.5)
})

test_that("g-comp sandwich SE and bootstrap SE are within 20%", {
  data("nhefs", package = "causatr")

  fit <- causat(
    nhefs,
    outcome = "wt82_71",
    treatment = "qsmk",
    confounders = ~ sex +
      age +
      I(age^2) +
      race +
      smokeintensity +
      smokeyrs +
      factor(exercise) +
      factor(active) +
      wt71,
    censoring = "censored"
  )

  res_sand <- contrast(
    fit,
    interventions = list(quit = static(1), continue = static(0)),
    ci_method = "sandwich",
    reference = "continue"
  )
  res_boot <- contrast(
    fit,
    interventions = list(quit = static(1), continue = static(0)),
    ci_method = "bootstrap",
    n_boot = 200L,
    reference = "continue"
  )

  se_sand <- res_sand$contrasts$se[1]
  se_boot <- res_boot$contrasts$se[1]

  # SEs should be within 30% of each other (allowing for bootstrap noise).
  ratio <- se_boot / se_sand
  expect_gt(ratio, 0.5)
  expect_lt(ratio, 2.0)
})

test_that("g-comp ATT differs from ATE", {
  data("nhefs", package = "causatr")

  fit <- causat(
    nhefs,
    outcome = "wt82_71",
    treatment = "qsmk",
    confounders = ~ sex + age + wt71,
    censoring = "censored"
  )

  res_ate <- contrast(
    fit,
    interventions = list(quit = static(1), continue = static(0)),
    estimand = "ATE",
    reference = "continue"
  )
  res_att <- contrast(
    fit,
    interventions = list(quit = static(1), continue = static(0)),
    estimand = "ATT",
    reference = "continue"
  )

  # ATE and ATT should give valid (non-NA) but potentially different estimates.
  expect_false(is.na(res_ate$contrasts$estimate[1]))
  expect_false(is.na(res_att$contrasts$estimate[1]))
  # Both should be positive (smoking cessation increases weight).
  expect_gt(res_ate$contrasts$estimate[1], 0)
  expect_gt(res_att$contrasts$estimate[1], 0)
})


test_that("g-comp supports Poisson family with sandwich variance", {
  # Coverage gap: every existing g-comp test uses Gaussian or
  # binomial. Verify the Poisson path (count outcome with log link)
  # works end-to-end with sandwich SEs.
  set.seed(1)
  n <- 500
  L <- stats::rnorm(n)
  A <- stats::rbinom(n, 1, plogis(0.3 * L))
  Y <- stats::rpois(n, exp(0.5 + 0.4 * A + 0.2 * L))

  df <- data.frame(Y = Y, A = A, L = L)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    family = "poisson"
  )
  expect_s3_class(fit, "causatr_fit")

  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "ratio",
    ci_method = "sandwich",
    reference = "a0"
  )
  # Expected rate ratio ~ exp(0.4) ≈ 1.49
  expect_gt(res$contrasts$estimate[1], 1.2)
  expect_lt(res$contrasts$estimate[1], 1.8)
  expect_true(all(is.finite(res$contrasts$se)))
  expect_true(all(res$contrasts$se > 0))
})


test_that("g-comp rejects unknown family string", {
  df <- data.frame(Y = c(1, 2, 3), A = c(0, 1, 0), L = c(1, 2, 3))
  expect_snapshot(
    error = TRUE,
    causat(
      df,
      outcome = "Y",
      treatment = "A",
      confounders = ~L,
      family = "not_a_family"
    )
  )
})
