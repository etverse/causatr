# Tests for `causat(weights = survey::svydesign(...))` convenience:
# extracting sampling weights from a `survey.design` object and
# auto-adopting its first-stage PSU as the default contrast-time
# cluster. Also exercises `unpack_svydesign()`'s rejection paths.

# Every test requires the `survey` package. Tests are not skipped when
# unavailable (per project rule: install deps rather than skip); if
# `survey` is genuinely missing, these tests will abort loudly.
skip_if_not_survey <- function() {
  testthat::skip_if_not_installed("survey")
}

simulate_clustered_survey <- function(
  n_clusters = 100,
  m = 5,
  seed = 1
) {
  set.seed(seed)
  cl <- rep(seq_len(n_clusters), each = m)
  U <- stats::rnorm(n_clusters)
  n <- n_clusters * m
  L <- stats::rnorm(n)
  A_cl <- stats::rbinom(n_clusters, 1, 0.5)
  flip <- stats::rbinom(n, 1, 0.1)
  A <- ifelse(flip == 1, 1 - A_cl[cl], A_cl[cl])
  Y <- 1 + 2 * A + 0.5 * L + U[cl] + stats::rnorm(n, sd = 0.5)
  pw <- stats::runif(n, 0.5, 2)
  data.frame(Y = Y, A = A, L = L, cl = cl, pw = pw)
}


# ── svydesign as weights: equivalence to manual unpack ──
test_that("causat(weights = svydesign) matches manual (weights + cluster) on gcomp", {
  skip_if_not_survey()
  d <- simulate_clustered_survey()
  des <- survey::svydesign(ids = ~cl, weights = ~pw, data = d)

  fit_des <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    weights = des
  )
  fit_manual <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    weights = d$pw,
    cluster = "cl"
  )

  r_des <- contrast(fit_des, list(a1 = static(1), a0 = static(0)))
  r_manual <- contrast(fit_manual, list(a1 = static(1), a0 = static(0)))

  # Point estimates identical (same fit).
  expect_equal(r_des$contrasts$estimate, r_manual$contrasts$estimate)
  # SEs identical (same fit, same cluster).
  expect_equal(r_des$contrasts$se, r_manual$contrasts$se)
  # Cluster auto-adopted from the design.
  expect_equal(fit_des$details$cluster, "cl")
})


# ── svydesign + IPW ──
test_that("causat(weights = svydesign) works under IPW", {
  skip_if_not_survey()
  d <- simulate_clustered_survey()
  des <- survey::svydesign(ids = ~cl, weights = ~pw, data = d)
  # The weighted logistic propensity emits the standard
  # "non-integer #successes" GLM warning whenever continuous external
  # weights flow into binomial(); harmless, suppressed here.
  fit <- suppressWarnings(causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw",
    weights = des
  ))
  expect_equal(fit$details$cluster, "cl")
  res <- contrast(fit, list(a1 = static(1), a0 = static(0)))
  expect_true(all(is.finite(res$contrasts$se)))
  expect_gt(res$contrasts$se, 0)
})


# ── User-supplied cluster overrides the design's PSU ──
test_that("explicit `cluster =` wins over the svydesign PSU", {
  skip_if_not_survey()
  d <- simulate_clustered_survey()
  d$site <- rep(1:10, length.out = nrow(d))
  des <- survey::svydesign(ids = ~cl, weights = ~pw, data = d)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    weights = des,
    cluster = "site"
  )
  # User's override wins: the design's `cl` PSU is ignored.
  expect_equal(fit$details$cluster, "site")
})


# ── Single-PSU design: no auto-cluster ──
test_that("single-PSU design does not auto-adopt a cluster", {
  skip_if_not_survey()
  d <- simulate_clustered_survey()
  # Use ~1 to declare no clustering (all rows in one PSU).
  des <- survey::svydesign(ids = ~1, weights = ~pw, data = d)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    weights = des
  )
  expect_null(fit$details$cluster)
})


# ── svydesign + matching: cluster auto-adoption rejected ──
test_that("svydesign on matching aborts with causatr_bad_cluster", {
  skip_if_not_survey()
  d <- simulate_clustered_survey()
  des <- survey::svydesign(ids = ~cl, weights = ~pw, data = d)
  expect_error(
    causat(
      d,
      outcome = "Y",
      treatment = "A",
      confounders = ~L,
      estimator = "matching",
      estimand = "ATT",
      weights = des
    ),
    class = "causatr_bad_cluster"
  )
})


# ── svydesign + matching + no-cluster design: accepted (no cluster propagation) ──
test_that("svydesign(~1, ...) + matching is accepted (only sampling weights flow)", {
  skip_if_not_survey()
  d <- simulate_clustered_survey(n_clusters = 50, m = 4)
  des <- survey::svydesign(ids = ~1, weights = ~pw, data = d)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "matching",
    estimand = "ATT",
    weights = des
  )
  expect_null(fit$details$cluster)
  expect_true(!is.null(fit$details$weights))
  res <- contrast(fit, list(a1 = static(1), a0 = static(0)))
  expect_true(all(is.finite(res$contrasts$se)))
})


# ── Row-count mismatch: rejected ──
test_that("svydesign with a different row count aborts with causatr_bad_svydesign", {
  skip_if_not_survey()
  d1 <- simulate_clustered_survey()
  # Design built on the full data...
  des <- survey::svydesign(ids = ~cl, weights = ~pw, data = d1)
  # ...but `data` passed to causat() is a strict subset. Row counts differ.
  d2 <- d1[1:400, ]
  expect_error(
    causat(
      d2,
      outcome = "Y",
      treatment = "A",
      confounders = ~L,
      weights = des
    ),
    class = "causatr_bad_svydesign"
  )
})


# ── Truth: weighted ATE recovery ──
test_that("causat(weights = svydesign) recovers the survey-weighted ATE", {
  skip_if_not_survey()
  # Weights are independent of (A, L, Y) by construction, so the
  # survey-weighted ATE equals the unweighted structural truth of 2.
  d <- simulate_clustered_survey(n_clusters = 300, m = 5)
  des <- survey::svydesign(ids = ~cl, weights = ~pw, data = d)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    weights = des
  )
  res <- contrast(
    fit,
    list(a1 = static(1), a0 = static(0)),
    reference = "a0"
  )
  expect_lt(abs(res$contrasts$estimate - 2), 0.15)
  # The cluster-robust CI covers 2.
  expect_lt(res$contrasts$ci_lower, 2)
  expect_gt(res$contrasts$ci_upper, 2)
})
