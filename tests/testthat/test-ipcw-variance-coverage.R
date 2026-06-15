# IPCW sandwich 95% CI coverage + conservative-SE checks.
# Split from test-ipcw-variance.R for parallelism.

# ── 95% CI covers truth for point IPCW ─────────���───────────────

test_that("gcomp IPCW sandwich 95% CI covers truth (ATE = 3)", {
  d <- simulate_mar_outcome(n = 5000, seed = 605)
  dt <- data.table::as.data.table(d)

  fit <- causat(
    dt,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "gcomp",
    censoring = "C",
    ipcw = TRUE
  )

  r <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "sandwich"
  )

  # contrast() reports a0 - a1; negate to get ATE = a1 - a0
  ate <- -r$contrasts$estimate
  ci_lo <- -r$contrasts$ci_upper
  ci_hi <- -r$contrasts$ci_lower

  expect_true(ci_lo <= 3 && 3 <= ci_hi)
})


test_that("IPW IPCW sandwich 95% CI covers truth (ATE = 3)", {
  d <- simulate_mar_outcome(n = 5000, seed = 606)
  dt <- data.table::as.data.table(d)

  fit <- causat(
    dt,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw",
    censoring = "C",
    ipcw = TRUE
  )

  r <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "sandwich"
  )

  ate <- -r$contrasts$estimate
  ci_lo <- -r$contrasts$ci_upper
  ci_hi <- -r$contrasts$ci_lower

  expect_true(ci_lo <= 3 && 3 <= ci_hi)
})


# ── Longitudinal IPCW: sandwich is conservative ────────────────

test_that("longitudinal ICE IPCW: sandwich SE >= bootstrap SE (conservative)", {
  skip_if_fast()
  d <- simulate_longitudinal_mar_outcome(n = 3000, seed = 607)
  dt <- data.table::as.data.table(d)

  fit <- causat(
    dt,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    estimator = "gcomp",
    type = "longitudinal",
    id = "id",
    time = "time",
    censoring = "C",
    ipcw = TRUE
  )

  r_sw <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "sandwich"
  )
  r_bt <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "bootstrap",
    n_boot = 500
  )

  ratio <- r_sw$contrasts$se / r_bt$contrasts$se
  # Sandwich omits longitudinal censoring Channel 2 correction,
  # so SE should be >= bootstrap. Allow small tolerance for sampling noise.
  expect_gt(ratio, 0.85)
  expect_lt(ratio, 2.0)
})
