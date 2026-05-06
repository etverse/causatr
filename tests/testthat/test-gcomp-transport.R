test_that("gcomp transport: fit_rows restricted to study rows (S=1)", {
  d <- simulate_transport(n = 500, seed = 1)
  fit <- causat(
    data = d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "gcomp",
    target = "S"
  )
  # Fitting rows must be S=1 and non-missing Y
  expected <- which(d$S == 1L & !is.na(d$Y))
  actual <- which(fit$details$fit_rows)
  expect_equal(sort(actual), sort(expected))
  expect_equal(fit$details$n_fit, length(expected))
  expect_equal(fit$details$target, "S")
})

test_that("gcomp transport: target_subset stored on fit", {
  d <- simulate_transport(n = 200, seed = 2)
  fit_tgt <- causat(
    d,
    "Y",
    "A",
    ~L,
    estimator = "gcomp",
    target = "S",
    target_subset = "target"
  )
  fit_all <- causat(
    d,
    "Y",
    "A",
    ~L,
    estimator = "gcomp",
    target = "S",
    target_subset = "all"
  )
  expect_equal(fit_tgt$target_subset, "target")
  expect_equal(fit_all$target_subset, "all")
})

test_that("gcomp transport: check_transport_inputs rejects ATT/ATC", {
  d <- simulate_transport(n = 200, seed = 3)
  expect_snapshot(
    causat(
      d,
      "Y",
      "A",
      ~L,
      estimator = "gcomp",
      target = "S",
      estimand = "ATT"
    ),
    error = TRUE
  )
  expect_snapshot(
    causat(
      d,
      "Y",
      "A",
      ~L,
      estimator = "gcomp",
      target = "S",
      estimand = "ATC"
    ),
    error = TRUE
  )
})

test_that("gcomp transport (transportability): recovers target ATE ≈ 3", {
  # DGP-T1: Y = 2 + 3A + 1.5L + A*L + eps
  # Target ATE = 3 + E[L | S=0]. Because P(S=1|L) = expit(-0.5 + L),
  # S=0 rows over-represent L<0, so E[L|S=0] < 0, making target ATE < 3.
  # But the MARGINAL ATE (over L ~ N(0,1)) = 3 + E[L] = 3.
  # Transportability here targets S=0 pop whose ATE != 3 unless we compute
  # the truth from the data.
  set.seed(42)
  d <- simulate_transport(n = 6000, seed = 42)

  fit <- causat(
    data = d,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L + A:L,
    estimator = "gcomp",
    target = "S",
    target_subset = "target"
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "sandwich"
  )

  # Truth: 3 + E[L | S=0] from large-n simulation
  truth_ate <- 3 + mean(d$L[d$S == 0])
  est <- coef(res)["a1"] - coef(res)["a0"]
  expect_lt(abs(est - truth_ate), 0.15)
})

test_that("gcomp transport (generalizability): recovers marginal ATE ≈ 3", {
  d <- simulate_transport(n = 6000, seed = 99)

  fit <- causat(
    data = d,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L + A:L,
    estimator = "gcomp",
    target = "S",
    target_subset = "all"
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "sandwich"
  )

  # Target = all rows; marginal ATE = 3 + E[L] ≈ 3 (since E[L] = 0 for N(0,1))
  truth_ate <- 3 + mean(d$L)
  est <- coef(res)["a1"] - coef(res)["a0"]
  expect_lt(abs(est - truth_ate), 0.15)
})

test_that("gcomp transport corrects study bias vs. naive study-only estimate", {
  d <- simulate_transport(n = 6000, seed = 7)

  # Transport: averages over target (S=0) rows
  fit_transport <- causat(
    d,
    "Y",
    "A",
    ~ L + A:L,
    estimator = "gcomp",
    target = "S",
    target_subset = "target"
  )
  res_transport <- contrast(
    fit_transport,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "sandwich"
  )

  # Naive study-only: fit and average over S=1 rows only
  d_study <- d[d$S == 1, ]
  fit_naive <- causat(d_study, "Y", "A", ~ L + A:L, estimator = "gcomp")
  res_naive <- contrast(
    fit_naive,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "sandwich"
  )

  truth_target <- 3 + mean(d$L[d$S == 0])
  truth_study <- 3 + mean(d$L[d$S == 1])

  est_transport <- coef(res_transport)["a1"] - coef(res_transport)["a0"]
  est_naive <- coef(res_naive)["a1"] - coef(res_naive)["a0"]

  # Transport estimate closer to target truth than naive
  expect_lt(abs(est_transport - truth_target), abs(est_naive - truth_target))
  # Naive is biased toward study truth (E[L|S=1] > 0)
  expect_lt(abs(est_naive - truth_study), 0.15)
})

test_that("gcomp transport: sandwich SE is plausible (ratio to bootstrap in (0.5, 2))", {
  d <- simulate_transport(n = 2000, seed = 11)
  fit <- causat(
    d,
    "Y",
    "A",
    ~ L + A:L,
    estimator = "gcomp",
    target = "S",
    target_subset = "target"
  )
  ivs <- list(a1 = static(1), a0 = static(0))
  res_sw <- contrast(
    fit,
    interventions = ivs,
    type = "difference",
    ci_method = "sandwich"
  )
  res_bt <- contrast(
    fit,
    interventions = ivs,
    type = "difference",
    ci_method = "bootstrap",
    n_boot = 200
  )

  v_sw <- vcov(res_sw)
  se_sw <- sqrt(v_sw["a1", "a1"] + v_sw["a0", "a0"] - 2 * v_sw["a1", "a0"])
  v_bt <- vcov(res_bt)
  se_bt <- sqrt(v_bt["a1", "a1"] + v_bt["a0", "a0"] - 2 * v_bt["a1", "a0"])
  ratio <- se_sw / se_bt
  expect_gt(ratio, 0.5)
  expect_lt(ratio, 2.0)
})

test_that("gcomp transport: target rows with NA outcome/treatment are handled", {
  # S=0 rows naturally have NA Y and NA A in simulate_transport()
  # The estimator should work without errors
  d <- simulate_transport(n = 500, seed = 5)
  expect_true(all(is.na(d$Y[d$S == 0])))
  expect_true(all(is.na(d$A[d$S == 0])))

  fit <- causat(d, "Y", "A", ~L, estimator = "gcomp", target = "S")
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "sandwich"
  )
  expect_s3_class(res, "causatr_result")
  expect_false(anyNA(coef(res)))
})

test_that("gcomp without transport is unaffected (target = NULL default)", {
  d <- simulate_transport(n = 500, seed = 6)
  d_study <- d[d$S == 1, ]

  # No target arg -- should fit and average over all study rows normally
  fit <- causat(d_study, "Y", "A", ~L, estimator = "gcomp")
  expect_null(fit$target)
  expect_null(fit$details$target)

  res <- contrast(
    fit,
    list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "sandwich"
  )
  expect_s3_class(res, "causatr_result")
})
