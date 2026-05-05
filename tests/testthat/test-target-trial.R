test_that("target_trial() returns correct S3 class", {
  tt <- target_trial(
    eligibility = "Adults aged 25-74",
    treatment_strategy = "Quit smoking vs. continue"
  )
  expect_s3_class(tt, "causatr_target_trial")
})

test_that("target_trial() stores only non-NULL fields", {
  tt <- target_trial(eligibility = "All adults", outcome = "Weight change")
  expect_named(tt, c("eligibility", "outcome"))
  expect_equal(tt$eligibility, "All adults")
  expect_equal(tt$outcome, "Weight change")
})

test_that("target_trial() with all fields populates correctly", {
  tt <- target_trial(
    eligibility = "A",
    treatment_strategy = "B",
    assignment = "C",
    follow_up = "D",
    outcome = "E",
    causal_contrast = "F",
    model = "G",
    censoring = "H"
  )
  expect_length(tt, 8)
})

test_that("target_trial() rejects non-string arguments", {
  expect_snapshot(target_trial(eligibility = 42), error = TRUE)
  expect_snapshot(target_trial(outcome = c("a", "b")), error = TRUE)
})

test_that("print.causatr_target_trial() snapshot", {
  tt <- target_trial(
    eligibility = "Adults aged 25-74 who smoke >= 5 cig/day",
    treatment_strategy = "Quit smoking vs. continue smoking",
    follow_up = "Baseline (1971) to end of follow-up (1982)",
    outcome = "Weight change in kg at end of follow-up",
    causal_contrast = "ATE: E[Y(quit)] - E[Y(continue)]"
  )
  expect_snapshot(print(tt))
})

test_that("print omits NULL fields", {
  tt <- target_trial(eligibility = "Everyone")
  out <- capture.output(print(tt))
  expect_true(any(grepl("Eligibility", out)))
  expect_false(any(grepl("Outcome", out)))
  expect_false(any(grepl("Censoring", out)))
})

test_that("target_trial() with no arguments returns empty object", {
  tt <- target_trial()
  expect_s3_class(tt, "causatr_target_trial")
  expect_length(tt, 0)
})
