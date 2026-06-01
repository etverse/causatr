# Boot MI (von Hippel 2020) variance engine. The decomposition is unit-tested
# exactly via pool_table_boot(); the full causat_mice(pool_method="boot_mi")
# path is checked for truth recovery, structure, reproducibility, and the
# uncongeniality correction relative to Rubin's rules.

test_that("pool_table_boot() implements V = between - within/M exactly", {
  # Build a B x M x 1 array with known per-cell values so the random-effects
  # decomposition can be checked against a hand computation.
  B <- 5
  M <- 2
  cells <- matrix(
    c(
      1.0,
      1.2,
      0.8,
      1.0,
      1.5,
      1.3,
      0.9,
      1.1,
      1.2,
      1.0
    ),
    nrow = B,
    byrow = TRUE
  )
  arr <- array(cells, dim = c(B, M, 1))
  pool <- pool_table_boot(arr, scale = "linear", conf_level = 0.95)

  boot_means <- rowMeans(cells)
  b_comp <- stats::var(boot_means)
  w_comp <- mean(apply(cells, 1, stats::var))
  v_hat <- max(b_comp - w_comp / M, .Machine$double.eps)

  expect_equal(pool$estimate[1], mean(boot_means), tolerance = 1e-12)
  expect_equal(pool$se[1], sqrt(v_hat), tolerance = 1e-12)
  expect_equal(pool$df[1], B - 1L)
  expect_equal(pool$diagnostics$between[1], b_comp, tolerance = 1e-12)
  expect_equal(pool$diagnostics$within[1], w_comp, tolerance = 1e-12)
})

test_that("pool_table_boot() pools ratios on the log scale", {
  B <- 4
  M <- 2
  cells <- matrix(
    c(2, 2.2, 1.8, 2.0, 2.4, 2.1, 1.9, 2.3),
    nrow = B,
    byrow = TRUE
  )
  arr <- array(cells, dim = c(B, M, 1))
  pool <- pool_table_boot(arr, scale = "log", conf_level = 0.95)
  expect_equal(
    pool$estimate[1],
    exp(mean(rowMeans(log(cells)))),
    tolerance = 1e-12
  )
  expect_gt(pool$ci_lower[1], 0) # ratio CI stays on the positive scale
})

test_that("boot_mi recovers the ATE and reports a finite SE", {
  skip_if_not_installed("mice")
  dat <- simulate_mi_covariate(n = 1500, seed = 12)
  imp <- mice::mice(dat[, c("Y", "A", "L")], m = 5, printFlag = FALSE)
  res <- causat_mice(
    imp,
    "Y",
    "A",
    ~L,
    list(a1 = static(1), a0 = static(0)),
    estimator = "gcomp",
    pool_method = "boot_mi",
    B = 40,
    M = 2,
    seed = 1
  )

  expect_identical(res$ci_method, "boot_mi")
  expect_equal(abs(res$contrasts$estimate[1]), 3, tolerance = 0.2)
  expect_true(is.finite(res$contrasts$se[1]))
  mi <- attr(res, "mi_details")
  expect_identical(mi$pool_method, "boot_mi")
  expect_true(mi$B <= 40L)
})

test_that("boot_mi rejects B < 2 and M < 2", {
  skip_if_not_installed("mice")
  dat <- simulate_mi_covariate(n = 400, seed = 13)
  imp <- mice::mice(dat[, c("Y", "A", "L")], m = 3, printFlag = FALSE)
  expect_error(
    causat_mice(
      imp,
      "Y",
      "A",
      ~L,
      list(a1 = static(1), a0 = static(0)),
      estimator = "gcomp",
      pool_method = "boot_mi",
      B = 1,
      M = 2
    ),
    class = "causatr_mi_bad_boot"
  )
  expect_error(
    causat_mice(
      imp,
      "Y",
      "A",
      ~L,
      list(a1 = static(1), a0 = static(0)),
      estimator = "gcomp",
      pool_method = "boot_mi",
      B = 20,
      M = 1
    ),
    class = "causatr_mi_bad_boot"
  )
})

test_that("boot_mi is reproducible across parallel backends with a fixed seed", {
  skip_if_not_installed("mice")
  skip_if_not_installed("future.apply")
  dat <- simulate_mi_covariate(n = 800, seed = 14)
  imp <- mice::mice(dat[, c("Y", "A", "L")], m = 3, printFlag = FALSE)
  args <- list(
    imp,
    "Y",
    "A",
    ~L,
    list(a1 = static(1), a0 = static(0)),
    estimator = "gcomp",
    pool_method = "boot_mi",
    B = 12,
    M = 2,
    seed = 99
  )
  r_seq <- do.call(causat_mice, args)
  r_fut <- do.call(causat_mice, c(args, list(parallel = "future")))
  # future_lapply uses parallel-safe RNG streams; the point estimate (a mean
  # over resamples) is stable even though the exact stream differs, so check
  # the estimate to a loose tolerance rather than bit-for-bit.
  expect_equal(
    r_seq$contrasts$estimate[1],
    r_fut$contrasts$estimate[1],
    tolerance = 0.3
  )
})

test_that("Boot MI corrects Rubin's conservatism under uncongeniality", {
  skip_on_cran()
  skip_if_not_installed("mice")
  # DGP-MI5 has an A*L interaction; imputing L without Y is uncongenial.
  # Rubin's between-imputation variance then overestimates the imputation
  # uncertainty, so its pooled SE is larger than Boot MI's. Average over a
  # few datasets to damp Monte Carlo noise in the comparison.
  rubin_se <- numeric(4)
  boot_se <- numeric(4)
  for (s in 1:4) {
    dat <- simulate_mi_uncongenial(n = 1200, seed = 200 + s)
    pred <- mice::make.predictorMatrix(dat[, c("Y", "A", "L")])
    pred[, "Y"] <- 0 # exclude Y -> uncongenial
    imp <- mice::mice(
      dat[, c("Y", "A", "L")],
      m = 10,
      predictorMatrix = pred,
      printFlag = FALSE
    )
    ints5 <- list(a1 = static(1), a0 = static(0))
    # Excluding Y from the predictor matrix is exactly the uncongeniality
    # condition, so the missing-predictor warning is expected here. Assert it
    # fired and keep the returned result (expect_warning returns the
    # condition, not the value, so capture the warning explicitly instead).
    expect_warn_value <- function(expr) {
      seen <- FALSE
      val <- withCallingHandlers(
        expr,
        causatr_mi_missing_predictors = function(w) {
          seen <<- TRUE
          invokeRestart("muffleWarning")
        }
      )
      expect_true(seen)
      val
    }
    rr <- expect_warn_value(
      causat_mice(imp, "Y", "A", ~L, ints5, estimator = "gcomp")
    )
    bb <- expect_warn_value(
      causat_mice(
        imp,
        "Y",
        "A",
        ~L,
        ints5,
        estimator = "gcomp",
        pool_method = "boot_mi",
        B = 60,
        M = 2,
        seed = s
      )
    )
    rubin_se[s] <- rr$contrasts$se[1]
    boot_se[s] <- bb$contrasts$se[1]
    # Both must still recover the marginal truth of 3.
    expect_equal(abs(rr$contrasts$estimate[1]), 3, tolerance = 0.3)
    expect_equal(abs(bb$contrasts$estimate[1]), 3, tolerance = 0.3)
  }
  expect_gt(mean(rubin_se), mean(boot_se))
})
