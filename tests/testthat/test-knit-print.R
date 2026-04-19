# knit_print.causatr_result + .onLoad registration tests.
#
# `knit_print.causatr_result()` is a knitr S3 method registered at
# package load time via `.onLoad()` in R/zzz.R. Both functions ship at
# 0% coverage because they only fire when the package is rendered via
# knitr/Quarto. These tests exercise both paths directly so a regression
# in the table builder or the conditional registration is caught.

test_that(".onLoad registers knit_print.causatr_result with knitr", {
  skip_if_not_installed("knitr")
  # `.onLoad()` runs once when the namespace is first loaded; on a fresh
  # devtools::load_all() that already happened. We assert the end state
  # (knitr's S3 method registry knows about `knit_print.causatr_result`)
  # AND directly invoke .onLoad() so the function body itself is
  # covered by covr -- otherwise it sits at 0% because covr instruments
  # the namespace after .onLoad has already run.
  causatr:::.onLoad("dummy_libname", "causatr")

  fn <- utils::getS3method(
    "knit_print",
    "causatr_result",
    optional = TRUE,
    envir = asNamespace("knitr")
  )
  expect_true(is.function(fn))
  expect_identical(
    body(fn),
    body(getNamespace("causatr")$knit_print.causatr_result)
  )
})

test_that(".onLoad is a silent no-op when knitr is unavailable", {
  # The `requireNamespace("knitr", ...)` guard inside .onLoad() returns
  # FALSE if knitr is not installed; the function returns invisibly
  # without registering anything. We mock requireNamespace to FALSE
  # for the duration of this call to exercise the guard branch on
  # systems that DO have knitr installed.
  testthat::local_mocked_bindings(
    requireNamespace = function(package, ...) {
      if (identical(package, "knitr")) FALSE else TRUE
    },
    .package = "base"
  )
  expect_invisible(causatr:::.onLoad("dummy_libname", "causatr"))
})

test_that("knit_print.causatr_result emits header + two tinytable HTML tables", {
  skip_if_not_installed("knitr")
  skip_if_not_installed("tinytable")

  d <- simulate_binary_continuous(n = 200, seed = 21)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "gcomp"
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "sandwich"
  )
  expect_s3_class(res, "causatr_result")

  out <- knit_print.causatr_result(res)
  expect_s3_class(out, "knit_asis")
  txt <- as.character(out)
  # Metadata header surfaces the four core fields.
  expect_match(txt, "Estimator")
  expect_match(txt, "Estimand")
  expect_match(txt, "Contrast")
  expect_match(txt, "CI method")
  expect_match(txt, "N:")
  # Tinytable emits caption blocks for the two tables.
  expect_match(txt, "Intervention means")
  expect_match(txt, "Contrasts")
})

test_that("knit_print.causatr_result falls back when tinytable is missing", {
  skip_if_not_installed("knitr")

  d <- simulate_binary_continuous(n = 100, seed = 22)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "gcomp"
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference"
  )

  # Force the `requireNamespace("tinytable", ...)` guard inside
  # knit_print.causatr_result() to fail by stubbing requireNamespace
  # in the function's enclosing environment. The fallback path returns
  # `knitr::normal_print(x)`, which dispatches to print.causatr_result.
  testthat::local_mocked_bindings(
    requireNamespace = function(package, ...) {
      if (identical(package, "tinytable")) FALSE else TRUE
    },
    .package = "base"
  )
  out <- capture.output(knit_print.causatr_result(res))
  # `print.causatr_result()` prints a banner with "Estimator:" followed
  # by the estimator's display label ("G-computation" for gcomp).
  expect_true(any(grepl("Estimator", out)))
  expect_true(any(grepl("G-computation", out)))
})
