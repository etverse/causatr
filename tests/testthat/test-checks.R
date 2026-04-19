# Unit tests for the small input-validation helpers in R/checks.R and
# R/utils.R. Each one is mostly a one-line guard that's exercised
# transitively by the public API on the happy path; the abort branches
# are easy to miss in coverage but matter because they are the user's
# first signal that an input is wrong.

# ---- check_string() -----------------------------------------------------

test_that("check_string() accepts a single character string silently", {
  expect_invisible(check_string("hello"))
  expect_null(check_string("x"))
})

test_that("check_string() aborts on non-string inputs", {
  expect_error(check_string(1L), "must be a single character string")
  expect_error(check_string(c("a", "b")), "must be a single character string")
  expect_error(check_string(NA_character_), "must be a single character string")
  expect_error(check_string(NULL), "must be a single character string")
})

test_that("check_string() reports the user-supplied arg name", {
  bad <- 1L
  expect_error(check_string(bad), "`bad`")
})

# ---- check_formula() ----------------------------------------------------

test_that("check_formula() accepts one- and two-sided formulas silently", {
  expect_invisible(check_formula(~ L1 + L2))
  expect_invisible(check_formula(Y ~ A + L))
})

test_that("check_formula() aborts on non-formula inputs", {
  expect_error(check_formula("Y ~ A"), "must be a formula")
  expect_error(check_formula(list()), "must be a formula")
  expect_error(check_formula(NULL), "must be a formula")
})

test_that("check_formula() reports the user-supplied arg name", {
  bad_arg <- "not a formula"
  expect_error(check_formula(bad_arg), "`bad_arg`")
})

# ---- check_col_exists() -------------------------------------------------

test_that("check_col_exists() accepts an existing column silently", {
  d <- data.frame(A = 1:3, L = 1:3)
  expect_invisible(check_col_exists(d, "A"))
})

test_that("check_col_exists() aborts on a missing column", {
  d <- data.frame(A = 1:3)
  expect_error(check_col_exists(d, "L"), "Column `L`.*not found")
})

# ---- check_pkg() --------------------------------------------------------

test_that("check_pkg() accepts an installed package silently", {
  # `stats` is part of base R -- always installed.
  expect_invisible(check_pkg("stats"))
})

test_that("check_pkg() aborts on a missing package with installation hint", {
  # Use a deliberately nonsense package name that no one will ever
  # publish to CRAN. The abort message must include the install
  # instructions so the user can act on it without further docs.
  expect_error(
    check_pkg("ZZZNotARealPackageZZZ"),
    "ZZZNotARealPackageZZZ.*install\\.packages"
  )
})

# ---- check_intervention_list() ------------------------------------------
#
# Top-level structure validator for the `interventions` argument.
# Three valid shapes per element: NULL (natural course), bare
# causatr_intervention, named list of causatr_intervention (multivariate).
# All other shapes abort with a class-specific message.

test_that("check_intervention_list() accepts a well-formed list", {
  expect_invisible(check_intervention_list(list(
    a1 = static(1),
    a0 = static(0),
    natural = NULL
  )))
})

test_that("check_intervention_list() accepts a multivariate sub-list", {
  expect_invisible(check_intervention_list(list(
    both1 = list(A1 = static(1), A2 = shift(0.5)),
    natural = NULL
  )))
})

test_that("check_intervention_list() aborts on non-list / empty inputs", {
  expect_error(check_intervention_list(NULL), "must be a named list")
  expect_error(check_intervention_list(list()), "must be a named list")
  expect_error(check_intervention_list("static"), "must be a named list")
})

test_that("check_intervention_list() aborts on missing names", {
  expect_error(
    check_intervention_list(list(static(1), static(0))),
    "must be named"
  )
  expect_error(
    check_intervention_list(list(a1 = static(1), static(0))),
    "must be named"
  )
})

test_that("check_intervention_list() aborts on duplicate names", {
  expect_error(
    check_intervention_list(list(a = static(1), a = static(0))),
    "duplicated name"
  )
})

test_that("check_intervention_list() aborts on a bare scalar that isn't an intervention", {
  expect_error(
    check_intervention_list(list(bad = "static(1)")),
    "must be a `causatr_intervention`"
  )
  expect_error(
    check_intervention_list(list(bad = 1)),
    "must be a `causatr_intervention`"
  )
})

test_that("check_intervention_list() aborts on an unnamed multivariate sub-list", {
  expect_error(
    check_intervention_list(list(
      both = list(static(1), shift(0.5)) # missing names inside
    )),
    "all elements are named"
  )
})

test_that("check_intervention_list() aborts on a multivariate sub-list element that isn't an intervention", {
  expect_error(
    check_intervention_list(list(
      both = list(A1 = static(1), A2 = "not an intervention")
    )),
    "must be a `causatr_intervention`"
  )
})

# ---- check_causat_inputs() ---------------------------------------------
#
# The top-level validator for `causat()`. Most argument-shape checks
# delegate to small helpers (check_string / check_formula / etc.,
# tested above), so the tests here target the validator's own
# additional branches: treatment-arg shape, outcome == treatment guard,
# and the `history` integer/Inf accumulator.

test_that("check_causat_inputs() rejects non-character or empty `treatment`", {
  d <- data.frame(Y = 1, A = 1, L = 1)
  expect_error(
    check_causat_inputs(
      data = d,
      outcome = "Y",
      treatment = 1L, # not character
      confounders = ~L,
      confounders_tv = NULL,
      id = NULL,
      time = NULL,
      history = NULL,
      estimator = "gcomp",
      estimand = "ATE"
    ),
    "must be a character string or character vector"
  )
  expect_error(
    check_causat_inputs(
      data = d,
      outcome = "Y",
      treatment = character(0),
      confounders = ~L,
      confounders_tv = NULL,
      id = NULL,
      time = NULL,
      history = NULL,
      estimator = "gcomp",
      estimand = "ATE"
    ),
    "must be a character string or character vector"
  )
})

test_that("check_causat_inputs() rejects outcome == treatment", {
  d <- data.frame(Y = 1, A = 1)
  expect_error(
    check_causat_inputs(
      data = d,
      outcome = "Y",
      treatment = "Y", # typo: outcome reused as treatment
      confounders = ~A,
      confounders_tv = NULL,
      id = NULL,
      time = NULL,
      history = NULL,
      estimator = "gcomp",
      estimand = "ATE"
    ),
    "must be different columns"
  )
})

test_that("check_causat_inputs() rejects history that isn't a positive integer or Inf", {
  d <- data.frame(Y = 1, A = 1, L = 1, id = 1, t = 0)
  # First branch: not scalar double / integer / Inf at all -- pass
  # a character to trip the type check.
  expect_error(
    check_causat_inputs(
      data = d,
      outcome = "Y",
      treatment = "A",
      confounders = ~L,
      confounders_tv = ~L,
      id = "id",
      time = "t",
      history = "k", # bad type
      estimator = "gcomp",
      estimand = "ATE"
    ),
    "must be a positive integer or `Inf`"
  )
})
