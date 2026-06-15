test_that("IPW transport: fit_rows restricted to study rows (S=1)", {
  d <- simulate_transport(n = 500, seed = 1)
  fit <- causat(
    data = d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw",
    target = "S"
  )
  expected <- which(d$S == 1L & !is.na(d$Y))
  actual <- which(fit$details$fit_rows)
  expect_equal(sort(actual), sort(expected))
  expect_equal(fit$details$n_fit, length(expected))
  expect_equal(fit$target, "S")
})

test_that("IPW transport: target_subset stored on fit", {
  d <- simulate_transport(n = 200, seed = 2)
  fit_tgt <- causat(
    d,
    "Y",
    "A",
    ~L,
    estimator = "ipw",
    target = "S",
    target_subset = "target"
  )
  fit_all <- causat(
    d,
    "Y",
    "A",
    ~L,
    estimator = "ipw",
    target = "S",
    target_subset = "all"
  )
  expect_equal(fit_tgt$target_subset, "target")
  expect_equal(fit_all$target_subset, "all")
})

test_that("IPW transport (transportability): recovers target ATE", {
  # DGP: Y = 2 + 3A + 1.5L + A*L + eps
  # Target ATE = 3 + E[L | S=0], computed from the drawn sample.
  # IPW is noisier than g-comp (SE ~ 0.07 at n=20000), so we use a
  # larger sample and tolerance of ~2 SE.
  d <- simulate_transport(n = 20000, seed = 42)

  fit <- causat(
    data = d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw",
    target = "S",
    target_subset = "target"
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "sandwich"
  )

  truth_ate <- 3 + mean(d$L[d$S == 0])
  est <- coef(res)["a1"] - coef(res)["a0"]
  expect_lt(abs(est - truth_ate), 0.15)
})

test_that("IPW transport (generalizability): recovers marginal ATE", {
  d <- simulate_transport(n = 20000, seed = 99)

  fit <- causat(
    data = d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw",
    target = "S",
    target_subset = "all"
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "sandwich"
  )

  # Target = all rows; marginal ATE = 3 + E[L] ~ 3
  truth_ate <- 3 + mean(d$L)
  est <- coef(res)["a1"] - coef(res)["a0"]
  expect_lt(abs(est - truth_ate), 0.15)
})

test_that("IPW transport corrects study bias vs. naive study-only estimate", {
  d <- simulate_transport(n = 20000, seed = 7)

  fit_transport <- causat(
    d,
    "Y",
    "A",
    ~L,
    estimator = "ipw",
    target = "S",
    target_subset = "target"
  )
  res_transport <- contrast(
    fit_transport,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "sandwich"
  )

  d_study <- d[d$S == 1, ]
  fit_naive <- causat(d_study, "Y", "A", ~L, estimator = "ipw")
  res_naive <- contrast(
    fit_naive,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "sandwich"
  )

  truth_target <- 3 + mean(d$L[d$S == 0])
  est_transport <- coef(res_transport)["a1"] - coef(res_transport)["a0"]
  est_naive <- coef(res_naive)["a1"] - coef(res_naive)["a0"]

  expect_lt(abs(est_transport - truth_target), abs(est_naive - truth_target))
})

test_that("IPW transport: sandwich SE plausible (ratio to bootstrap)", {
  skip_if_fast()
  d <- simulate_transport(n = 3000, seed = 11)
  fit <- causat(
    d,
    "Y",
    "A",
    ~L,
    estimator = "ipw",
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
  expect_gt(ratio, 0.6)
  expect_lt(ratio, 1.8)
})

test_that("IPW transport: target rows with NA outcome/treatment handled", {
  d <- simulate_transport(n = 500, seed = 5)
  expect_true(all(is.na(d$Y[d$S == 0])))
  expect_true(all(is.na(d$A[d$S == 0])))

  fit <- causat(d, "Y", "A", ~L, estimator = "ipw", target = "S")
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "sandwich"
  )
  expect_s3_class(res, "causatr_result")
  expect_false(anyNA(coef(res)))
})

test_that("IPW transport cross-check: matches TransportHealth::transportIP()", {
  skip_if_fast()
  # TransportHealth::transportIP() fits the same models: P(A=1|L) logistic for
  # confounding, P(S=1|L) logistic for sampling, then a weighted MSM Y ~ A on
  # study rows. For binary A the MSM coefficient equals the Hajek ATE, matching
  # causatr's per-intervention intercept-only MSM approach.
  # Install: pak::pkg_install("CoreClinicalSciences/TransportHealth")
  skip_if_not_installed("TransportHealth")
  d <- simulate_transport(n = 10000, seed = 7)
  d_th <- data.frame(d)
  d_th$A <- factor(d_th$A)

  invisible(capture.output(
    res_ip <- TransportHealth::transportIP(
      msmFormula = Y ~ A,
      propensityScoreModel = A ~ L,
      participationModel = S ~ L,
      treatment = "A",
      participation = "S",
      response = "Y",
      data = d_th,
      transport = TRUE,
      bootstrapNum = 0
    )
  ))
  th_ate <- unname(stats::coef(res_ip$msm)["A1"])

  fit <- causat(
    d,
    "Y",
    "A",
    ~L,
    estimator = "ipw",
    target = "S",
    target_subset = "target"
  )
  res <- contrast(
    fit,
    list(a1 = static(1), a0 = static(0)),
    type = "difference",
    ci_method = "sandwich"
  )
  c_ate <- unname(coef(res)["a1"] - coef(res)["a0"])

  expect_equal(c_ate, th_ate, tolerance = 1e-4)
})

test_that("IPW without transport is unaffected (target = NULL default)", {
  d <- simulate_transport(n = 500, seed = 6)
  d_study <- d[d$S == 1, ]

  fit <- causat(d_study, "Y", "A", ~L, estimator = "ipw")
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
