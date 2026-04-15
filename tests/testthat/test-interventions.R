test_that("static() creates a causatr_intervention", {
  iv <- static(1)
  expect_s3_class(iv, "causatr_intervention")
  expect_equal(iv$type, "static")
  expect_equal(iv$value, 1)
})

test_that("shift() creates a causatr_intervention", {
  iv <- shift(-10)
  expect_s3_class(iv, "causatr_intervention")
  expect_equal(iv$type, "shift")
  expect_equal(iv$delta, -10)
})

test_that("scale_by() creates a causatr_intervention", {
  iv <- scale_by(0.5)
  expect_s3_class(iv, "causatr_intervention")
  expect_equal(iv$type, "scale")
  expect_equal(iv$factor, 0.5)
})

test_that("threshold() creates a causatr_intervention", {
  iv <- threshold(0, 20)
  expect_s3_class(iv, "causatr_intervention")
  expect_equal(iv$type, "threshold")
  expect_equal(iv$lower, 0)
  expect_equal(iv$upper, 20)
})

test_that("dynamic() creates a causatr_intervention", {
  rule <- \(data, trt) ifelse(data$x > 0, 1, 0)
  iv <- dynamic(rule)
  expect_s3_class(iv, "causatr_intervention")
  expect_equal(iv$type, "dynamic")
  expect_true(is.function(iv$rule))
})

test_that("dynamic() rejects non-functions", {
  expect_snapshot(error = TRUE, dynamic("not a function"))
})

test_that("ipsi() creates a causatr_intervention", {
  iv <- ipsi(2)
  expect_s3_class(iv, "causatr_intervention")
  expect_equal(iv$type, "ipsi")
  expect_equal(iv$delta, 2)
})

test_that("ipsi() rejects non-positive delta", {
  expect_snapshot(error = TRUE, ipsi(0))
  expect_snapshot(error = TRUE, ipsi(-1))
})

test_that("dynamic() refuses non-numeric treatment columns", {
  # 2026-04-15 fourth-round critical review Issue #10: previously
  # `apply_single_intervention()` only validated the return type of
  # the dynamic rule, not the treatment column's own type. A rule
  # returning numeric silently coerced a character / factor column,
  # which then blew up downstream in `predict()` with a cryptic
  # factor-level mismatch (or, worse, silently predicted on a numeric
  # dummy the outcome model never saw). Repro: character treatment
  # "low"/"high", rule returning `rep(1, n)`. Pre-fix behavior: column
  # overwritten to `c(1,1,1,1)`, class numeric. Fix: abort at apply
  # time with a clean pointer to `static()`.
  d <- data.table::data.table(A = c("low", "high", "low", "high"), Y = 1:4)
  iv <- dynamic(function(data, trt) rep(1, nrow(data)))
  expect_error(
    apply_intervention(d, "A", iv),
    "numeric treatment column"
  )

  # Factor treatment — same rejection path.
  d2 <- data.table::data.table(
    A = factor(c("low", "high", "low", "high")),
    Y = 1:4
  )
  expect_error(
    apply_intervention(d2, "A", iv),
    "numeric treatment column"
  )

  # Numeric treatment still works (regression guard).
  d3 <- data.table::data.table(A = c(0.1, 0.5, 0.9, 0.3), Y = 1:4)
  iv_num <- dynamic(function(data, trt) as.numeric(trt > 0.4))
  res <- apply_intervention(d3, "A", iv_num)
  expect_equal(res$A, c(0, 1, 1, 0))
})

test_that("ipsi() rejects NA / NaN delta with a clean message", {
  # 2026-04-15 fourth-round critical review Issue #9: `is_scalar_double`
  # returns TRUE for NA_real_, so the following `delta <= 0` evaluated
  # to NA and users saw the cryptic "missing value where TRUE/FALSE
  # needed" abort. The explicit `is.na(delta)` guard mirrors the
  # defensive checks already present in shift/scale_by/threshold.
  expect_error(ipsi(NA_real_), "non-NA positive")
  expect_error(ipsi(NaN),     "non-NA positive")
})

test_that("static() rejects invalid inputs", {
  expect_snapshot(error = TRUE, static(NA))
  expect_snapshot(error = TRUE, static(c(1, 2)))
})

test_that("static() accepts character values", {
  iv <- static("A")
  expect_s3_class(iv, "causatr_intervention")
  expect_equal(iv$value, "A")
})

test_that("static() with character value works through full pipeline", {
  set.seed(99)
  n <- 500
  L <- rnorm(n)
  A <- sample(
    c("low", "med", "high"),
    n,
    replace = TRUE,
    prob = c(0.3, 0.4, 0.3)
  )
  Y <- 1 + 2 * (A == "med") + 4 * (A == "high") + 0.5 * L + rnorm(n)
  d <- data.frame(Y = Y, A = A, L = L)

  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "gcomp"
  )
  res <- contrast(
    fit,
    interventions = list(low = static("low"), high = static("high")),
    type = "difference"
  )
  expect_true(res$contrasts$estimate > 2)
  expect_true(res$contrasts$se > 0)
})

test_that("shift() rejects invalid inputs", {
  expect_snapshot(error = TRUE, shift(NA))
  expect_snapshot(error = TRUE, shift("a"))
})

test_that("scale_by() rejects invalid inputs", {
  expect_snapshot(error = TRUE, scale_by(NA))
  expect_snapshot(error = TRUE, scale_by("a"))
})

test_that("threshold() rejects invalid inputs", {
  expect_snapshot(error = TRUE, threshold(5, 3))
  expect_snapshot(error = TRUE, threshold(NA, 10))
})

test_that("print.causatr_intervention() works", {
  expect_snapshot(print(static(1)))
  expect_snapshot(print(shift(-5)))
})
