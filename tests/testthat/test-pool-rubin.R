# Rubin's-rules engine: unit tests against the mice::pool.scalar() oracle and
# the closed-form Barnard-Rubin degrees of freedom.

test_that("rubin_scalar() reproduces mice::pool.scalar() exactly", {
  skip_if_not_installed("mice")
  set.seed(1)
  m <- 12
  Q <- rnorm(m, mean = 3, sd = 0.4)
  U <- runif(m, 0.05, 0.15)
  n <- 500
  k <- 1

  ours <- rubin_scalar(Q, U, m = m, dfcom = n - k)
  oracle <- mice::pool.scalar(Q, U, n = n, k = k)

  expect_equal(ours$qbar, oracle$qbar, tolerance = 1e-12)
  expect_equal(ours$ubar, oracle$ubar, tolerance = 1e-12)
  expect_equal(ours$b, oracle$b, tolerance = 1e-12)
  expect_equal(ours$t, oracle$t, tolerance = 1e-12)
  expect_equal(ours$df, oracle$df, tolerance = 1e-8)
  expect_equal(ours$fmi, oracle$fmi, tolerance = 1e-8)
})

test_that("barnard_rubin_df() matches mice:::barnard.rubin()", {
  skip_if_not_installed("mice")
  # Exercise a grid of (m, between, total, dfcom) including the infinite
  # complete-data df (large-sample Rubin 1987) limit.
  grid <- expand.grid(
    m = c(2, 5, 20),
    b = c(0.01, 0.1),
    t = c(0.2, 0.5),
    dfcom = c(Inf, 50, 1000)
  )
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    expect_equal(
      barnard_rubin_df(g$m, g$b, g$t, g$dfcom),
      mice:::barnard.rubin(g$m, g$b, g$t, g$dfcom),
      tolerance = 1e-10,
      info = paste(g, collapse = ",")
    )
  }
})

test_that("rubin_scalar() total variance follows T = Ubar + (1 + 1/m) B", {
  set.seed(2)
  m <- 8
  Q <- rnorm(m, 1, 0.3)
  U <- rep(0.1, m)
  rs <- rubin_scalar(Q, U, m = m, dfcom = Inf)
  expect_equal(rs$ubar, 0.1)
  expect_equal(rs$b, stats::var(Q))
  expect_equal(rs$t, 0.1 + (1 + 1 / m) * stats::var(Q))
})

test_that("rubin_scalar() handles the m = 1 degenerate case", {
  rs <- rubin_scalar(Q = 2.5, U = 0.09, m = 1L, dfcom = 100)
  expect_equal(rs$qbar, 2.5)
  expect_equal(rs$b, 0)
  expect_equal(rs$t, 0.09) # total collapses to within
  expect_equal(rs$df, 100)
  expect_true(is.na(rs$fmi))
})

test_that("build_pooled_table() emits canonical columns and a by column", {
  pool <- list(
    estimate = c(1, 2),
    se = c(0.1, 0.2),
    ci_lower = c(0.8, 1.6),
    ci_upper = c(1.2, 2.4)
  )
  dt <- build_pooled_table(
    labels = c("a1", "a0"),
    label_col = "intervention",
    by = NULL,
    pool = pool
  )
  expect_s3_class(dt, "data.table")
  expect_named(dt, c("intervention", "estimate", "se", "ci_lower", "ci_upper"))

  dt_by <- build_pooled_table(
    labels = c("a1", "a0"),
    label_col = "intervention",
    by = c("F", "M"),
    pool = pool
  )
  expect_true("by" %in% names(dt_by))
  expect_equal(names(dt_by)[1:2], c("intervention", "by"))
})
