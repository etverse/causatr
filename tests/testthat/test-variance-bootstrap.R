# Unit tests for R/variance_bootstrap.R helpers (process_boot_results
# and refit_model). The end-to-end bootstrap path is already exercised
# by every test that uses `ci_method = "bootstrap"`. These tests target
# the per-replicate failure-handling and dispatcher branches that the
# happy path doesn't reach.

# ---- process_boot_results() failure-handling branches --------------------

# Helper: build a minimal `boot::boot`-shaped object with a hand-rolled
# `t` matrix so we can drive the failure-counting paths without
# actually running a bootstrap.
fake_boot_result <- function(t_mat) {
  structure(
    list(t = t_mat),
    class = "boot"
  )
}

test_that("process_boot_results() warns and returns NA vcov when < 2 replicates succeed", {
  # All-NA `t` matrix simulates the case where every replicate failed
  # (e.g. factor-level mismatch in every resample).
  t_mat <- matrix(NA_real_, nrow = 100, ncol = 2)
  br <- fake_boot_result(t_mat)
  expect_warning(
    out <- process_boot_results(br, c("a1", "a0"), n_boot = 100),
    "fewer than 2 non-NA"
  )
  expect_true(all(is.na(out$vcov)))
  expect_equal(dim(out$vcov), c(2L, 2L))
  expect_equal(out$boot_info$n_ok, 0L)
  expect_equal(out$boot_info$n_fail, 100L)
  expect_equal(out$boot_info$n_fail_by_int, c(a1 = 100L, a0 = 100L))
})

test_that("process_boot_results() emits the gentle warning when failure rate <= 20%", {
  # 10 of 100 (= 10%) replicates failed -- under the 20% threshold so
  # the message is the "discarded" variant rather than "may be unreliable".
  set.seed(41)
  t_mat <- matrix(rnorm(200), nrow = 100, ncol = 2)
  t_mat[1:10, ] <- NA_real_
  br <- fake_boot_result(t_mat)
  expect_warning(
    out <- process_boot_results(br, c("a1", "a0"), n_boot = 100),
    "discarded"
  )
  expect_equal(out$boot_info$n_ok, 90L)
  expect_equal(out$boot_info$n_fail, 10L)
  # Variance is real: built from the 90 successful replicates.
  expect_true(all(is.finite(out$vcov)))
})

test_that("process_boot_results() emits the strong warning when failure rate > 20%", {
  # 30 of 100 (= 30%) replicates failed -- triggers the
  # "may be unreliable" wording that signals model instability.
  set.seed(42)
  t_mat <- matrix(rnorm(200), nrow = 100, ncol = 2)
  t_mat[1:30, ] <- NA_real_
  br <- fake_boot_result(t_mat)
  expect_warning(
    out <- process_boot_results(br, c("a1", "a0"), n_boot = 100),
    "may be unreliable"
  )
  expect_equal(out$boot_info$n_ok, 70L)
})

test_that("process_boot_results() succeeds silently when no replicates fail", {
  set.seed(43)
  t_mat <- matrix(rnorm(200), nrow = 100, ncol = 2)
  br <- fake_boot_result(t_mat)
  # No NAs -> no warning. expect_no_warning was added in testthat 3.2.
  expect_no_warning(
    out <- process_boot_results(br, c("a1", "a0"), n_boot = 100)
  )
  expect_equal(out$boot_info$n_ok, 100L)
  expect_equal(out$boot_info$n_fail, 0L)
})

# ---- refit_model() dispatcher --------------------------------------------

test_that("refit_model() aborts on an unknown estimator", {
  # The dispatcher routes to refit_gcomp / refit_ipw / refit_matching
  # for the supported estimators and aborts on anything else. The abort
  # is a defensive guard for a future estimator that ships without a
  # bootstrap refit hook.
  fake_fit <- list(estimator = "frobnitz")
  expect_error(
    refit_model(fake_fit, d_b = data.frame()),
    "Bootstrap is not supported"
  )
})
