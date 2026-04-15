# Unit tests for `check_estimand_intervention_compat()`
#
# This gate protects users from silently getting the wrong estimand
# when they combine ATT / ATC with a non-static intervention on the
# IPW or matching engine. See `PHASE_4_INTERVENTIONS_SELF_IPW.md §8`
# for the rationale.
#
# The tests exercise the helper directly rather than going through
# `contrast()` because the Phase-3 `check_interventions_compat()` gate
# still rejects any non-static intervention under IPW / matching
# upstream of this check — so until the Phase 4 `fit_ipw()` rewrite
# lands, the new gate is only reachable via the internal helper. Once
# `fit_ipw()` is rewritten, the existing integration tests (e.g.
# `test-contrast.R`) will pick up the same gate through the full
# pipeline.

# ---- g-comp: no-op regardless of (estimand, intervention) ----------

test_that("gcomp skips the check for every (estimand, intervention)", {
  # Each sub-case would fail under IPW / matching but must pass for
  # gcomp — the outcome model is estimand-agnostic, so reweighting
  # over a different target population works for any intervention.
  for (est in c("ATE", "ATT", "ATC")) {
    expect_silent(
      check_estimand_intervention_compat(
        est,
        list(a = static(1)),
        "gcomp"
      )
    )
    expect_silent(
      check_estimand_intervention_compat(
        est,
        list(sh = shift(-0.5)),
        "gcomp"
      )
    )
    expect_silent(
      check_estimand_intervention_compat(
        est,
        list(ips = ipsi(2), sh = shift(1)),
        "gcomp"
      )
    )
  }
})

# ---- ATE + IPW / matching: no-op for every intervention ------------

test_that("ATE allows every intervention under IPW and matching", {
  # ATE is well-defined for every intervention the IPW engine
  # supports. The new gate is only concerned with ATT / ATC.
  ivs <- list(
    a1 = static(1),
    a0 = static(0),
    sh = shift(-0.5),
    sc = scale_by(0.75),
    ip = ipsi(1.8),
    dy = dynamic(function(d, a) as.numeric(d$L > 0)),
    nc = NULL
  )
  expect_silent(
    check_estimand_intervention_compat("ATE", ivs, "ipw")
  )
  expect_silent(
    check_estimand_intervention_compat("ATE", ivs, "matching")
  )
})

# ---- ATT + non-static: the core rejection case ---------------------

test_that("ATT rejects shift under IPW with `causatr_bad_estimand_intervention`", {
  # This is the primary enforcement path. The error class is how
  # the Phase 4 doc pins this test, and downstream code that wants to
  # react programmatically (e.g. a future "auto-promote ATT to ATE"
  # pathway) can trap on it.
  expect_error(
    check_estimand_intervention_compat(
      "ATT",
      list(sh = shift(-0.5)),
      "ipw"
    ),
    class = "causatr_bad_estimand_intervention"
  )
})

test_that("ATT rejects ipsi / scale_by / dynamic under IPW", {
  expect_error(
    check_estimand_intervention_compat(
      "ATT",
      list(ip = ipsi(2)),
      "ipw"
    ),
    class = "causatr_bad_estimand_intervention"
  )
  expect_error(
    check_estimand_intervention_compat(
      "ATT",
      list(sc = scale_by(0.5)),
      "ipw"
    ),
    class = "causatr_bad_estimand_intervention"
  )
  expect_error(
    check_estimand_intervention_compat(
      "ATT",
      list(dy = dynamic(function(d, a) as.numeric(d$L > 0))),
      "ipw"
    ),
    class = "causatr_bad_estimand_intervention"
  )
})

test_that("ATC rejects non-static under matching the same way", {
  # The gate is symmetric across ATT / ATC and across the two engines
  # that bake the estimand into their weights / matched sets.
  expect_error(
    check_estimand_intervention_compat(
      "ATC",
      list(sh = shift(1)),
      "matching"
    ),
    class = "causatr_bad_estimand_intervention"
  )
})

# ---- Mixed lists: only non-static entries must be named in the message

test_that("the error message lists the offending intervention names", {
  # A common user pattern: request a bundle of contrasts where some
  # are static-vs-static and one is accidentally a shift. The gate
  # should name the specific offender(s) so the user knows what to
  # drop or switch, not just that "something is wrong".
  err <- tryCatch(
    check_estimand_intervention_compat(
      "ATT",
      list(
        a1 = static(1),
        a0 = static(0),
        drop5 = shift(-5),
        half = scale_by(0.5)
      ),
      "ipw"
    ),
    error = function(e) e
  )
  expect_s3_class(err, "causatr_bad_estimand_intervention")
  msg <- conditionMessage(err)
  expect_match(msg, "drop5", fixed = TRUE)
  expect_match(msg, "half", fixed = TRUE)
  # Static entries should NOT appear in the "offending" list.
  expect_false(grepl("'a1'", msg, fixed = TRUE))
  expect_false(grepl("'a0'", msg, fixed = TRUE))
})

# ---- Natural course entries are always allowed --------------------

test_that("NULL (natural course) entries are allowed under ATT / ATC", {
  # The observed marginal mean among the treated / controls is a
  # well-defined quantity for any estimator, so a `NULL` entry in the
  # interventions list should never trip the gate.
  expect_silent(
    check_estimand_intervention_compat(
      "ATT",
      list(obs = NULL),
      "ipw"
    )
  )
  expect_silent(
    check_estimand_intervention_compat(
      "ATT",
      list(a1 = static(1), obs = NULL),
      "matching"
    )
  )
})

# ---- Multivariate interventions ------------------------------------

test_that("multivariate interventions are flagged sub-entry by sub-entry", {
  # Multivariate treatment: `iv` is a plain named list of
  # `causatr_intervention` objects, one per treatment column. Any
  # non-static sub-entry flags the whole regime.
  mv_bad <- list(A1 = static(1), A2 = shift(-1))
  expect_error(
    check_estimand_intervention_compat(
      "ATT",
      list(joint = mv_bad),
      "ipw"
    ),
    class = "causatr_bad_estimand_intervention"
  )

  mv_good <- list(A1 = static(1), A2 = static(0))
  expect_silent(
    check_estimand_intervention_compat(
      "ATT",
      list(joint = mv_good),
      "ipw"
    )
  )
})

# ---- Error message carries both hint lines ------------------------

test_that("error message points at both `estimand = 'ATE'` and `gcomp`", {
  # The hint lines are load-bearing — users who hit this error
  # need to know they have two ways out (switch estimand, or switch
  # estimator). A snapshot would be more robust but testthat skips
  # snapshots on CRAN by default, so we pin the two key strings
  # directly here.
  err <- tryCatch(
    check_estimand_intervention_compat(
      "ATT",
      list(drop5 = shift(-5)),
      "ipw"
    ),
    error = function(e) e
  )
  msg <- conditionMessage(err)
  expect_match(msg, "estimand = 'ATE'", fixed = TRUE)
  expect_match(msg, "estimator = 'gcomp'", fixed = TRUE)
  expect_match(msg, "drop5", fixed = TRUE)
})
