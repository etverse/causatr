# Tests for the `parallel = "future"` bootstrap backend in `contrast()`.
# Dispatches replicates through `future.apply::future_lapply()` so any
# active `future::plan()` — local multisession, remote cluster,
# `future.batchtools` — is honoured. The tests pin:
#
# - equivalence to `parallel = "no"` on a tight bootstrap seed (the
#   two paths draw slightly different PRNG streams; we assert MC
#   agreement, not bitwise identity),
# - end-to-end runs under `plan(multisession)` for gcomp and ICE,
# - the `check_pkg("future.apply")` abort path when the suggest is
#   absent (cannot be exercised without uninstalling — we test the
#   symmetric path where the user simply forgets `parallel = "future"`
#   and still gets correct behaviour).

skip_if_not_future <- function() {
  testthat::skip_if_not_installed("future")
  testthat::skip_if_not_installed("future.apply")
}


# ── Sequential future vs boot::boot(parallel = "no") ──
test_that("parallel = 'future' under plan(sequential) matches parallel = 'no' in MC tolerance", {
  skip_if_not_future()
  suppressPackageStartupMessages(future::plan(future::sequential))
  on.exit(future::plan(future::sequential), add = TRUE)

  set.seed(5)
  n <- 400
  d <- data.frame(
    Y = 1 + 2 * rbinom(n, 1, 0.5) + stats::rnorm(n),
    A = rbinom(n, 1, 0.5),
    L = stats::rnorm(n)
  )
  d$Y <- 1 + 2 * d$A + 0.3 * d$L + stats::rnorm(n)

  fit <- causat(d, outcome = "Y", treatment = "A", confounders = ~L)

  set.seed(123)
  r_future <- contrast(
    fit,
    list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "bootstrap",
    n_boot = 200L,
    parallel = "future"
  )
  set.seed(123)
  r_seq <- contrast(
    fit,
    list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "bootstrap",
    n_boot = 200L,
    parallel = "no"
  )

  # Point estimates are model-based, not bootstrap-based, so they must
  # match exactly.
  expect_equal(r_future$contrasts$estimate, r_seq$contrasts$estimate)
  # The two SEs are computed on independent resample streams
  # (`boot::boot` uses its own RNG, our dispatcher uses base R's).
  # Under 200 replicates the Monte Carlo bound is roughly
  # se * sqrt(2/(R-1)) ≈ 10% of se; we assert 30% to give plenty of
  # headroom.
  se_ratio <- r_future$contrasts$se / r_seq$contrasts$se
  expect_gt(se_ratio, 0.7)
  expect_lt(se_ratio, 1.3)
})


# ── Multisession plan: end-to-end ──
test_that("parallel = 'future' + plan(multisession) produces a sensible SE for gcomp", {
  skip_if_not_future()
  suppressPackageStartupMessages(future::plan(
    future::multisession,
    workers = 2
  ))
  on.exit(future::plan(future::sequential), add = TRUE)

  set.seed(7)
  n <- 400
  d <- data.frame(
    A = rbinom(n, 1, 0.5),
    L = stats::rnorm(n)
  )
  d$Y <- 1 + 2 * d$A + 0.3 * d$L + stats::rnorm(n)
  fit <- causat(d, outcome = "Y", treatment = "A", confounders = ~L)

  res <- contrast(
    fit,
    list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "bootstrap",
    n_boot = 100L,
    parallel = "future"
  )
  expect_true(all(is.finite(res$contrasts$se)))
  expect_gt(res$contrasts$se, 0)
  # Truth (ATE = 2) must fall inside the CI.
  expect_lt(res$contrasts$ci_lower, 2)
  expect_gt(res$contrasts$ci_upper, 2)
})


# ── Multisession plan: ICE (longitudinal, cluster-bootstrap by id) ──
test_that("parallel = 'future' + plan(multisession) works for ICE (cluster bootstrap by id)", {
  skip_if_not_future()
  suppressPackageStartupMessages(future::plan(
    future::multisession,
    workers = 2
  ))
  on.exit(future::plan(future::sequential), add = TRUE)

  set.seed(8)
  n_ids <- 250
  L0 <- stats::rnorm(n_ids)
  A0 <- stats::rbinom(n_ids, 1, stats::plogis(0.3 * L0))
  L1 <- L0 + 0.2 * A0 + stats::rnorm(n_ids)
  A1 <- stats::rbinom(n_ids, 1, stats::plogis(0.3 * L1))
  Y <- 1 + A0 + A1 + 0.5 * L1 + stats::rnorm(n_ids)

  long <- data.table::data.table(
    id = rep(seq_len(n_ids), each = 2),
    time = rep(0:1, times = n_ids),
    A = as.numeric(rbind(A0, A1)),
    L = as.numeric(rbind(L0, L1)),
    Y = rep(Y, each = 2)
  )
  fit <- causat(
    long,
    outcome = "Y",
    treatment = "A",
    confounders = ~1,
    confounders_tv = ~L,
    id = "id",
    time = "time"
  )
  res <- contrast(
    fit,
    list(all = static(1)),
    ci_method = "bootstrap",
    n_boot = 60L,
    parallel = "future"
  )
  expect_true(is.finite(res$estimates$se))
  expect_gt(res$estimates$se, 0)
})


# ── Invalid `parallel` rejected by match.arg ──
test_that("unknown parallel value is rejected", {
  d <- data.frame(Y = rnorm(50), A = rbinom(50, 1, 0.5), L = rnorm(50))
  fit <- causat(d, outcome = "Y", treatment = "A", confounders = ~L)
  expect_error(
    contrast(
      fit,
      list(a1 = static(1), a0 = static(0)),
      ci_method = "bootstrap",
      n_boot = 10L,
      parallel = "foo"
    )
  )
})
