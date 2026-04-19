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
  expect_error(ipsi(NaN), "non-NA positive")
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

# ---- apply_intervention() multivariate dispatch -----------------------------
#
# `apply_intervention()` has a defensive guard for multivariate
# treatments: every name in the sub-intervention list must match a
# column in `treatment`. The guard catches user typos like
# `list(X = ..., A2 = ...)` for `treatment = c("A1", "A2")` -- without
# it, data.table's `[, (X) := ...]` would silently create a spurious
# column instead of touching A1, biasing the contrast. The single-
# treatment fast path (`length(treatment) == 1L`) and the natural-
# course shortcut (`is.null(iv)`) are also exercised here so the
# branches are pinned.

test_that("apply_intervention() returns data unchanged when iv is NULL", {
  d <- data.table::data.table(A = c(0, 1, 0), L = c(1, 2, 3))
  out <- apply_intervention(d, "A", NULL)
  expect_equal(out$A, c(0, 1, 0))
  # Result must be a copy: mutating it must not touch the caller's data.
  out$A <- c(9, 9, 9)
  expect_equal(d$A, c(0, 1, 0))
})

test_that("apply_intervention() applies a static intervention to a single treatment", {
  d <- data.table::data.table(A = c(0, 1, 0, 1), L = c(1, 2, 3, 4))
  out <- apply_intervention(d, "A", static(1))
  expect_equal(out$A, c(1, 1, 1, 1))
  expect_equal(d$A, c(0, 1, 0, 1)) # original untouched
})

test_that("apply_intervention() multivariate: applies named per-treatment interventions", {
  d <- data.table::data.table(A1 = c(0, 1), A2 = c(2, 3), L = c(1, 2))
  out <- apply_intervention(
    d,
    c("A1", "A2"),
    list(A1 = static(1), A2 = shift(0.5))
  )
  expect_equal(out$A1, c(1, 1))
  expect_equal(out$A2, c(2.5, 3.5))
})

test_that("apply_intervention() multivariate aborts on missing intervention name", {
  d <- data.table::data.table(A1 = c(0, 1), A2 = c(2, 3))
  expect_error(
    apply_intervention(
      d,
      c("A1", "A2"),
      list(A1 = static(1)) # A2 missing
    ),
    "Missing"
  )
})

test_that("apply_intervention() multivariate aborts on extra (unexpected) intervention name", {
  d <- data.table::data.table(A1 = c(0, 1), A2 = c(2, 3))
  expect_error(
    apply_intervention(
      d,
      c("A1", "A2"),
      list(A1 = static(1), A2 = static(0), A3 = static(0))
    ),
    "Unexpected"
  )
})

# ---- apply_single_intervention() defensive branches ----------------------
#
# This is the data.table-mutating companion to
# `apply_intervention_to_values()`: same family of dispatch logic but
# operates on a `data.table` column in place. The dynamic-rule type
# checks here are stricter than the value-only counterpart because
# they protect downstream `predict()` from receiving a wrong-type
# treatment column. All abort branches need direct exercising.

test_that("apply_single_intervention() static / shift / scale / threshold update column in place", {
  d <- data.table::data.table(A = c(1, 2, 3, 4))
  apply_single_intervention(d, "A", static(99))
  expect_equal(d$A, c(99, 99, 99, 99))

  d <- data.table::data.table(A = c(1, 2, 3, 4))
  apply_single_intervention(d, "A", shift(0.5))
  expect_equal(d$A, c(1.5, 2.5, 3.5, 4.5))

  d <- data.table::data.table(A = c(1, 2, 3, 4))
  apply_single_intervention(d, "A", scale_by(2))
  expect_equal(d$A, c(2, 4, 6, 8))

  d <- data.table::data.table(A = c(-2, 0, 2, 5))
  apply_single_intervention(d, "A", threshold(0, 3))
  expect_equal(d$A, c(0, 0, 2, 3))
})

test_that("apply_single_intervention() dynamic on numeric A: rejects non-numeric return", {
  d <- data.table::data.table(A = c(1.0, 2.0, 3.0))
  bad_rule <- dynamic(function(d, a) c("a", "b", "c"))
  expect_error(
    apply_single_intervention(d, "A", bad_rule),
    "rule returned character but treatment column"
  )
})

test_that("apply_single_intervention() dynamic on factor A: char-with-known-level coerces back", {
  d <- data.table::data.table(
    A = factor(c("a", "b", "a"), levels = c("a", "b"))
  )
  good_rule <- dynamic(function(d, a) c("b", "a", "b"))
  apply_single_intervention(d, "A", good_rule)
  expect_s3_class(d$A, "factor")
  expect_equal(as.character(d$A), c("b", "a", "b"))
  expect_equal(levels(d$A), c("a", "b"))
})

test_that("apply_single_intervention() dynamic on factor A: char-with-unknown-level aborts", {
  d <- data.table::data.table(
    A = factor(c("a", "b", "a"), levels = c("a", "b"))
  )
  bad_rule <- dynamic(function(d, a) c("a", "z", "a"))
  expect_error(
    apply_single_intervention(d, "A", bad_rule),
    "value\\(s\\) not present as levels"
  )
})

test_that("apply_single_intervention() dynamic on factor A: factor-with-mismatched-levels aborts", {
  d <- data.table::data.table(
    A = factor(c("a", "b", "a"), levels = c("a", "b"))
  )
  wrong_factor <- factor(c("a", "b", "a"), levels = c("a", "b", "c"))
  bad_rule <- dynamic(function(d, a) wrong_factor)
  expect_error(
    apply_single_intervention(d, "A", bad_rule),
    "factor with levels that do not match"
  )
})

test_that("apply_single_intervention() dynamic on factor A: numeric return aborts", {
  d <- data.table::data.table(
    A = factor(c("a", "b", "a"), levels = c("a", "b"))
  )
  bad_rule <- dynamic(function(d, a) c(1.5, 2.5, 1.5))
  expect_error(
    apply_single_intervention(d, "A", bad_rule),
    "double but treatment column.*is a factor"
  )
})

test_that("apply_single_intervention() dynamic on character A: non-character return aborts", {
  d <- data.table::data.table(A = c("a", "b", "a"))
  bad_rule <- dynamic(function(d, a) c(1, 2, 1))
  expect_error(
    apply_single_intervention(d, "A", bad_rule),
    "double but treatment column.*is character"
  )
})

test_that("apply_single_intervention() dynamic on unsupported column type aborts", {
  # Logical treatment columns aren't a supported shape; the dynamic
  # branch falls through to the type-not-supported abort.
  d <- data.table::data.table(A = c(TRUE, FALSE, TRUE))
  rule <- dynamic(function(d, a) c(TRUE, TRUE, TRUE))
  expect_error(
    apply_single_intervention(d, "A", rule),
    "does not support treatment columns of type"
  )
})

test_that("apply_single_intervention() dynamic length-mismatch aborts", {
  d <- data.table::data.table(A = c(1.0, 2.0, 3.0))
  bad_rule <- dynamic(function(d, a) c(1, 2)) # length 2, should be 3
  expect_error(
    apply_single_intervention(d, "A", bad_rule),
    "must return a vector of length 3"
  )
})

test_that("apply_single_intervention() rejects ipsi with a clear message", {
  # IPSI is IPW-only; under gcomp / matching, apply_single_intervention()
  # is the gate that aborts. The error names the supported estimator
  # so the user immediately knows what to switch to.
  d <- data.table::data.table(A = c(0, 1, 0))
  expect_error(
    apply_single_intervention(d, "A", ipsi(2)),
    "only supported under .estimator = 'ipw'"
  )
})

test_that("apply_single_intervention() rejects an unknown intervention type", {
  d <- data.table::data.table(A = c(0, 1, 0))
  evil_iv <- structure(
    list(type = "frobnitz"),
    class = c("causatr_intervention", "list")
  )
  expect_error(
    apply_single_intervention(d, "A", evil_iv),
    "Unknown intervention type"
  )
})

test_that("apply_intervention() multivariate aborts on a typo (renamed sub-intervention)", {
  # Regression: pre-guard, this would silently create a spurious `X`
  # column via data.table's `[, (X) := ...]` and never touch A1.
  d <- data.table::data.table(A1 = c(0, 1), A2 = c(2, 3))
  expect_error(
    apply_intervention(
      d,
      c("A1", "A2"),
      list(X = static(1), A2 = shift(-10))
    ),
    "Missing.*A1.*Unexpected.*X"
  )
  # Original data must still be intact -- the guard fires before any
  # mutation reaches the data.table copy.
  expect_equal(d$A1, c(0, 1))
})
