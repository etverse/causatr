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
  A <- sample(c("low", "med", "high"), n,
    replace = TRUE, prob = c(0.3, 0.4, 0.3)
  )
  Y <- 1 + 2 * (A == "med") + 4 * (A == "high") + 0.5 * L + rnorm(n)
  d <- data.frame(Y = Y, A = A, L = L)

  fit <- causat(d, outcome = "Y", treatment = "A", confounders = ~L,
    method = "gcomp"
  )
  res <- contrast(fit,
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
