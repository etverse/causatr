# Unit tests for R/prepare_data.R helpers. The main prepare_data()
# function is exercised transitively by every causat() call; these
# tests target the warn_confounder_variation() guard, which warns when
# a user mis-classifies a confounder (baseline vs time-varying) for an
# ICE fit. A regression here would silently mis-specify the model
# without surfacing a diagnostic.

# Helper: build a 3-period long-format panel with one truly time-varying
# variable (Lt) and one truly baseline variable (sex).
mk_long_panel <- function(n_id = 30, seed = 51) {
  set.seed(seed)
  ids <- seq_len(n_id)
  long <- data.table::CJ(id = ids, t = 0:2)
  # sex is constant within id
  sex <- rbinom(n_id, 1, 0.5)
  long[, sex := sex[id]]
  # Lt varies within id
  long[, Lt := rnorm(.N)]
  # A baseline-but-actually-time-varying column for the warning test
  long[, Lwrong := rnorm(.N)]
  # A time-varying-but-actually-baseline column for the inverse warning
  baseline_const <- rnorm(n_id)
  long[, Lconst := baseline_const[id]]
  long
}

test_that("warn_confounder_variation() warns when a time-varying confounder is constant within id", {
  d <- mk_long_panel(n_id = 30, seed = 51)
  # `Lconst` is declared time-varying but is constant within every
  # individual -- the warning suggests moving it to baseline.
  # `confounders` only references genuinely-baseline `sex` so the
  # baseline-varies branch is silent.
  expect_warning(
    warn_confounder_variation(
      d,
      confounders = ~sex,
      confounders_tv = ~Lconst,
      id = "id"
    ),
    "does not vary within any individual"
  )
})

test_that("warn_confounder_variation() warns when a baseline confounder varies within id", {
  d <- mk_long_panel(n_id = 30, seed = 52)
  # `Lwrong` is declared baseline but varies within individuals --
  # the warning suggests moving it to confounders_tv.
  # `confounders_tv` only references genuinely-time-varying `Lt` so
  # the tv-not-varying branch is silent.
  expect_warning(
    warn_confounder_variation(
      d,
      confounders = ~ sex + Lwrong,
      confounders_tv = ~Lt,
      id = "id"
    ),
    "varies within some individuals"
  )
})

test_that("warn_confounder_variation() is silent when classifications are correct", {
  d <- mk_long_panel(n_id = 30, seed = 53)
  expect_no_warning(
    warn_confounder_variation(
      d,
      confounders = ~sex,
      confounders_tv = ~Lt,
      id = "id"
    )
  )
})

test_that("warn_confounder_variation() ignores the treatment column in baseline_vars", {
  # Regression: `~ L + A:L` would naively flag `A` as a baseline
  # confounder via `all.vars()`. The function subtracts the treatment
  # vector from baseline_vars before the within-id variation check,
  # so a varying treatment column does not trip the baseline warning.
  d <- mk_long_panel(n_id = 30, seed = 54)
  d[, A := rbinom(.N, 1, 0.5)]
  expect_no_warning(
    warn_confounder_variation(
      d,
      confounders = ~ sex + A:sex,
      confounders_tv = ~Lt,
      id = "id",
      treatment = "A"
    )
  )
})

test_that("warn_confounder_variation() works with NULL confounders_tv", {
  d <- mk_long_panel(n_id = 30, seed = 55)
  expect_no_warning(
    warn_confounder_variation(
      d,
      confounders = ~sex,
      confounders_tv = NULL,
      id = "id"
    )
  )
})
