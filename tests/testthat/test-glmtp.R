# Tests for natural-history modified treatment policies (G-LMTPs), the
# augmented-data sequential regression of Diaz, Williams, Morzywolek & Rudolph
# (2026, arXiv:2605.24167). The estimator is genuinely irreducible for policies
# that depend on the natural-value history of treatment: the standard ICE
# recursion would condition on the OBSERVED lag rather than the counterfactual
# natural value, returning a quietly-wrong number. The truth oracle is therefore
# forward Monte-Carlo simulation of the natural-history regime (lmtp's one-shot
# shift computes the standard LMTP, the wrong estimand here), backed by exact
# limiting-case equivalences against the existing engine.

# ---------------------------------------------------------------------------
# Chunk 1 -- augmented-frame infrastructure (R/glmtp_augment.R)
# ---------------------------------------------------------------------------

test_that("glmtp_support returns the sorted discrete support and drops censored rows", {
  dt <- data.table::data.table(
    A = c(0, 1, 1, 0, 2, 1),
    C = c(0, 0, 0, 0, 1, 0)
  )
  # Binary support, no censoring column.
  expect_identical(glmtp_support(dt[A %in% c(0, 1)], "A"), c(0, 1))
  # Ordinal support {0,1,2}; the censored row (C == 1, A == 2) is excluded,
  # so 2 drops out of the support entirely.
  expect_identical(glmtp_support(dt, "A", censoring = "C"), c(0, 1))
  # Without the censoring filter the full ordinal support is returned.
  expect_identical(glmtp_support(dt, "A"), c(0, 1, 2))
})

test_that("glmtp_support rejects continuous, factor, and multivariate treatment", {
  cont <- data.table::data.table(A = c(0.1, 1.7, 2.3))
  expect_error(glmtp_support(cont, "A"), class = "causatr_glmtp_continuous_trt")

  fac <- data.table::data.table(A = factor(c("lo", "hi", "lo")))
  expect_error(glmtp_support(fac, "A"), class = "causatr_glmtp_continuous_trt")

  mv <- data.table::data.table(A1 = c(0, 1), A2 = c(1, 0))
  expect_error(
    glmtp_support(mv, c("A1", "A2")),
    class = "causatr_glmtp_continuous_trt"
  )
})

test_that("glmtp_enumerate_labels builds the product support with the empty base", {
  support <- c(0, 1)
  # t = 0 -> a single empty label (the base of the recursion for q_1).
  l0 <- glmtp_enumerate_labels(support, 0L)
  expect_length(l0, 1L)
  expect_length(l0[[1L]], 0L)

  # t = 1 -> |A| labels.
  l1 <- glmtp_enumerate_labels(support, 1L)
  expect_length(l1, 2L)
  expect_setequal(vapply(l1, identity, numeric(1)), support)

  # t = 2 -> |A|^2 labels, all unique, first coordinate varies slowest.
  l2 <- glmtp_enumerate_labels(support, 2L)
  expect_length(l2, 4L)
  keys <- vapply(l2, glmtp_label_key, character(1))
  expect_length(unique(keys), 4L)

  # Ordinal support: cardinality is |A|^t.
  expect_length(glmtp_enumerate_labels(c(0, 1, 2), 3L), 27L)
})

test_that("glmtp_label_key is collision-free and maps the empty sequence to ''", {
  expect_identical(glmtp_label_key(numeric(0)), "")
  expect_identical(glmtp_label_key(c(0, 1, 0)), "0|1|0")
  # Distinct sequences that would collide under naive concatenation get
  # distinct keys (1,11) vs (11,1).
  expect_false(
    identical(glmtp_label_key(c(1, 11)), glmtp_label_key(c(11, 1)))
  )
})

test_that("glmtp_check_tractable caps the worst-step blow-up", {
  # Binary, 5 periods: worst step enumerates 2^4 = 16 labels -- within budget.
  expect_identical(glmtp_check_tractable(c(0, 1), 5L, budget = 1024L), 16)
  # A single period has no history to enumerate -> 1 label.
  expect_identical(glmtp_check_tractable(c(0, 1), 1L), 1)
  # Binary, 12 periods: 2^11 = 2048 > 1024 budget -> abort.
  expect_error(
    glmtp_check_tractable(c(0, 1), 12L, budget = 1024L),
    class = "causatr_glmtp_too_many"
  )
  # Raising the budget admits the same problem.
  expect_identical(
    glmtp_check_tractable(c(0, 1), 12L, budget = 4096L),
    2048
  )
})
