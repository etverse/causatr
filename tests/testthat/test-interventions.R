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

test_that("dynamic() enforces type compatibility with treatment column", {
  # 2026-04-15 fourth-round critical review Issue #10: previously
  # `apply_single_intervention()` only validated the return type of
  # the dynamic rule, not that it matched the treatment column's type.
  # A rule returning numeric on a character / factor column silently
  # coerced it, and downstream `predict()` either errored with a
  # cryptic factor-level mismatch or (worse) silently evaluated at a
  # numeric dummy the outcome model never saw. Fix: require
  # type-compatible returns — numeric→numeric, character→character,
  # factor accepts factor-with-matching-levels or character whose
  # values are all existing levels.

  iv_num <- dynamic(function(data, trt) rep(1, nrow(data)))

  # Numeric treatment still works.
  d_num <- data.table::data.table(A = c(0.1, 0.5, 0.9, 0.3), Y = 1:4)
  iv_rule <- dynamic(function(data, trt) as.numeric(trt > 0.4))
  expect_equal(apply_intervention(d_num, "A", iv_rule)$A, c(0, 1, 1, 0))

  # Character → numeric: rejected.
  d_chr <- data.table::data.table(
    A = c("low", "high", "low", "high"),
    Y = 1:4
  )
  expect_error(
    apply_intervention(d_chr, "A", iv_num),
    "character"
  )

  # Character → character: accepted (the real use case the earlier,
  # too-restrictive fix was blocking).
  iv_chr_rule <- dynamic(function(data, trt) {
    ifelse(seq_len(nrow(data)) %% 2L == 0L, "high", "low")
  })
  res_chr <- apply_intervention(d_chr, "A", iv_chr_rule)
  expect_equal(res_chr$A, c("low", "high", "low", "high"))
  expect_type(res_chr$A, "character")

  # Factor → numeric: rejected.
  d_fct <- data.table::data.table(
    A = factor(c("low", "high", "low", "high"), levels = c("low", "high")),
    Y = 1:4
  )
  expect_error(
    apply_intervention(d_fct, "A", iv_num),
    "factor"
  )

  # Factor → character with existing levels: accepted, coerced back to
  # factor with the original levels preserved.
  iv_fct_rule <- dynamic(function(data, trt) rep("high", nrow(data)))
  res_fct <- apply_intervention(d_fct, "A", iv_fct_rule)
  expect_s3_class(res_fct$A, "factor")
  expect_equal(levels(res_fct$A), c("low", "high"))
  expect_equal(as.character(res_fct$A), rep("high", 4))

  # Factor → character with an unknown level: rejected (rather than
  # silently becoming NA on coercion).
  iv_unknown <- dynamic(function(data, trt) rep("medium", nrow(data)))
  expect_error(
    apply_intervention(d_fct, "A", iv_unknown),
    "not present as levels"
  )
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
