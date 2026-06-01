# Triangulation: longitudinal SNM vs IPW-MSM vs ICE on shared DGPs
#
# Under correct specification and no time-varying EM, all three
# approaches — SNM blip sum, longitudinal IPW Hájek MSM, and ICE
# g-comp — must agree on the total causal effect. Disagreement is
# a bug, not a method difference (Phase 18 design doc invariant).
#
# Continuous-treatment shift is NOT included: the IPW sequential MTP
# estimand accounts for the A_0 -> L_1 -> A_1 cascade, making it
# structurally larger than the SNM per-stage blip sum. They estimate
# different quantities by design.

# --- Binary treatment: SNM vs IPW vs ICE (always=1 vs never=0) ----------

test_that("longitudinal SNM blip sum matches IPW ATE (binary, no EM)", {
  dgp <- simulate_snm_longitudinal(n = 5000, seed = 101)

  fit_snm <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders = ~1,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    type = "longitudinal",
    family = "gaussian",
    estimator = "snm",
    propensity_model_fn = stats::glm
  )
  res_snm <- contrast(fit_snm, ci_method = "sandwich")

  psi_snm <- res_snm$estimates$estimate
  names(psi_snm) <- res_snm$estimates$parameter
  snm_total <- sum(psi_snm)

  fit_ipw <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders = ~1,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    estimator = "ipw",
    propensity_model_fn = stats::glm
  )
  res_ipw <- contrast(
    fit_ipw,
    interventions = list(always = static(1), never = static(0)),
    type = "difference",
    reference = "never",
    ci_method = "sandwich"
  )
  ipw_ate <- res_ipw$contrasts$estimate

  # Both should recover the analytical truth: 3.15 + 3 = 6.15
  expect_equal(snm_total, 6.15, tolerance = 0.15)
  expect_equal(ipw_ate, 6.15, tolerance = 0.15)
  expect_equal(snm_total, ipw_ate, tolerance = 0.2)
})


test_that("longitudinal SNM blip sum matches ICE ATE (binary, no EM)", {
  dgp <- simulate_snm_longitudinal(n = 5000, seed = 102)

  fit_snm <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders = ~1,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    type = "longitudinal",
    family = "gaussian",
    estimator = "snm",
    propensity_model_fn = stats::glm
  )
  res_snm <- contrast(fit_snm, ci_method = "sandwich")
  snm_total <- sum(res_snm$estimates$estimate)

  fit_ice <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders = ~1,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    estimator = "gcomp",
    model_fn = stats::glm
  )
  res_ice <- contrast(
    fit_ice,
    interventions = list(always = static(1), never = static(0)),
    type = "difference",
    reference = "never",
    ci_method = "sandwich"
  )
  ice_ate <- res_ice$contrasts$estimate

  expect_equal(snm_total, 6.15, tolerance = 0.15)
  expect_equal(ice_ate, 6.15, tolerance = 0.15)
  expect_equal(snm_total, ice_ate, tolerance = 0.2)
})


test_that("3-way triangulation: SNM vs IPW vs ICE all agree (binary)", {
  dgp <- simulate_snm_longitudinal(n = 8000, seed = 103)

  fit_snm <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders = ~1,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    type = "longitudinal",
    family = "gaussian",
    estimator = "snm",
    propensity_model_fn = stats::glm
  )
  res_snm <- contrast(fit_snm, ci_method = "sandwich")
  snm_total <- sum(res_snm$estimates$estimate)

  fit_ipw <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders = ~1,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    estimator = "ipw",
    propensity_model_fn = stats::glm
  )
  res_ipw <- contrast(
    fit_ipw,
    interventions = list(always = static(1), never = static(0)),
    type = "difference",
    reference = "never",
    ci_method = "sandwich"
  )
  ipw_ate <- res_ipw$contrasts$estimate

  fit_ice <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders = ~1,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    estimator = "gcomp",
    model_fn = stats::glm
  )
  res_ice <- contrast(
    fit_ice,
    interventions = list(always = static(1), never = static(0)),
    type = "difference",
    reference = "never",
    ci_method = "sandwich"
  )
  ice_ate <- res_ice$contrasts$estimate

  expect_equal(snm_total, 6.15, tolerance = 0.1)
  expect_equal(ipw_ate, 6.15, tolerance = 0.1)
  expect_equal(ice_ate, 6.15, tolerance = 0.1)

  pairwise_diffs <- c(
    abs(snm_total - ipw_ate),
    abs(snm_total - ice_ate),
    abs(ipw_ate - ice_ate)
  )
  expect_true(all(pairwise_diffs < 0.15))
})


# --- Negative control: IPW-MSM biased with TV modifier, SNM is not ------
#
# On the TV-EM DGP, IPW-MSM conditioning on post-treatment M_1 is
# biased (collider bias; Robins 2000). SNM recovers the correct blip.

test_that("SNM correct but IPW biased with time-varying EM (negative control)", {
  dgp <- simulate_snm_longitudinal_tv_em(n = 8000, seed = 301)

  # SNM uses confounders_tv for the treatment model (no baseline needed)
  fit_snm <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders_outcome = ~ A:M,
    confounders_tv = ~ L + M,
    id = "id",
    time = "time",
    type = "longitudinal",
    family = "gaussian",
    estimator = "snm",
    propensity_model_fn = stats::glm,
    treatment_free = ~L
  )
  res_snm <- contrast(fit_snm, ci_method = "sandwich")
  psi_snm <- res_snm$estimates$estimate
  names(psi_snm) <- res_snm$estimates$parameter

  expect_equal(
    psi_snm[["stage0_psi_intercept"]],
    dgp$truth_psi[["stage0_psi_intercept"]],
    tolerance = 0.15
  )
  expect_equal(
    psi_snm[["stage1_psi_intercept"]],
    dgp$truth_psi[["stage1_psi_intercept"]],
    tolerance = 0.15
  )
  expect_equal(
    psi_snm[["stage1_psi_M"]],
    dgp$truth_psi[["stage1_psi_M"]],
    tolerance = 0.2
  )

  # IPW-MSM with A:M in confounders (conditioning on post-treatment M_1)
  # is biased — collider bias from M_1 being a descendant of A_0. `M` enters
  # `confounders` (baseline) via the interaction yet varies within individuals,
  # so causatr correctly warns; assert it rather than let it leak.
  expect_warning(
    fit_ipw <- causat(
      dgp$data,
      outcome = "Y",
      treatment = "A",
      confounders = ~ A:M,
      confounders_tv = ~ L + M,
      id = "id",
      time = "time",
      estimator = "ipw",
      propensity_model_fn = stats::glm
    ),
    class = "causatr_baseline_confounder_varying"
  )
  res_ipw <- contrast(
    fit_ipw,
    interventions = list(always = static(1), never = static(0)),
    type = "difference",
    reference = "never",
    ci_method = "sandwich"
  )

  # SNM averaged blip: stage0 + stage1 (with modifier averaged)
  snm_stage0 <- psi_snm[["stage0_psi_intercept"]]
  if ("stage0_psi_M" %in% names(psi_snm)) {
    snm_stage0 <- snm_stage0 +
      psi_snm[["stage0_psi_M"]] *
        mean(
          dgp$data[dgp$data$time == 0, ]$M
        )
  }
  snm_stage1 <- psi_snm[["stage1_psi_intercept"]] +
    psi_snm[["stage1_psi_M"]] * mean(dgp$data[dgp$data$time == 1, ]$M)

  snm_total <- snm_stage0 + snm_stage1
  ipw_ate <- res_ipw$contrasts$estimate

  # IPW estimate should differ from SNM total — structural collider bias
  expect_true(
    abs(ipw_ate - snm_total) > 0.1,
    label = paste0(
      "IPW (",
      round(ipw_ate, 3),
      ") should be biased away from SNM total (",
      round(snm_total, 3),
      ")"
    )
  )
})


# --- delicatessen cross-language reference ------------------------------------
# Reference values from data-raw/snm_longitudinal_reference.py using
# delicatessen stacked M-estimation (Zivich 2022) on the SAME shared
# fixture data (data-raw/snm_longitudinal_fixture.csv, n=5000, seed=101).
# Since both R and Python operate on the exact same dataset, tight
# tolerances (1e-4 for estimates, 1e-3 for SEs) are appropriate.

test_that("longitudinal SNM matches delicatessen (binary, no EM)", {
  ref <- read.csv(test_path("fixtures", "snm_longitudinal_delicatessen.csv"))

  # Load the shared fixture data (generated by R, read by Python).
  # Both R and Python operate on the exact same dataset.
  # Copy lives in fixtures/ for R CMD check; canonical source is data-raw/.
  wide <- read.csv(test_path(
    "fixtures",
    "snm_longitudinal_fixture.csv"
  ))
  d <- data.table::rbindlist(lapply(seq_len(nrow(wide)), function(i) {
    data.table::data.table(
      id = c(wide$id[i], wide$id[i]),
      time = c(0L, 1L),
      Y = c(NA_real_, wide$Y[i]),
      A = c(wide$A0[i], wide$A1[i]),
      L = c(wide$L0[i], wide$L1[i])
    )
  }))

  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~1,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    type = "longitudinal",
    family = "gaussian",
    estimator = "snm",
    propensity_model_fn = stats::glm
  )
  res <- contrast(fit, ci_method = "sandwich")
  psi <- res$estimates$estimate
  names(psi) <- res$estimates$parameter
  se <- res$estimates$se
  names(se) <- res$estimates$parameter

  ref_s0 <- ref[ref$parameter == "stage0_psi_intercept", ]
  ref_s1 <- ref[ref$parameter == "stage1_psi_intercept", ]

  expect_equal(psi[["stage0_psi_intercept"]], ref_s0$estimate, tolerance = 1e-4)
  expect_equal(psi[["stage1_psi_intercept"]], ref_s1$estimate, tolerance = 1e-4)
  # SE tolerance slightly wider (0.3%) because R's cluster-robust sandwich
  # and delicatessen's stacked sandwich use different bread/meat aggregation
  expect_equal(se[["stage0_psi_intercept"]], ref_s0$se, tolerance = 0.003)
  expect_equal(se[["stage1_psi_intercept"]], ref_s1$se, tolerance = 0.003)
})


test_that("longitudinal SNM + TF matches delicatessen (binary, no EM)", {
  ref <- read.csv(test_path(
    "fixtures",
    "snm_longitudinal_tf_delicatessen.csv"
  ))

  wide <- read.csv(test_path(
    "fixtures",
    "snm_longitudinal_fixture.csv"
  ))
  d <- data.table::rbindlist(lapply(seq_len(nrow(wide)), function(i) {
    data.table::data.table(
      id = c(wide$id[i], wide$id[i]),
      time = c(0L, 1L),
      Y = c(NA_real_, wide$Y[i]),
      A = c(wide$A0[i], wide$A1[i]),
      L = c(wide$L0[i], wide$L1[i])
    )
  }))

  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~1,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    type = "longitudinal",
    family = "gaussian",
    estimator = "snm",
    propensity_model_fn = stats::glm,
    treatment_free = ~L
  )
  res <- contrast(fit, ci_method = "sandwich")
  psi <- res$estimates$estimate
  names(psi) <- res$estimates$parameter
  se <- res$estimates$se
  names(se) <- res$estimates$parameter

  ref_s0 <- ref[ref$parameter == "stage0_psi_intercept", ]
  ref_s1 <- ref[ref$parameter == "stage1_psi_intercept", ]

  expect_equal(psi[["stage0_psi_intercept"]], ref_s0$estimate, tolerance = 1e-4)
  expect_equal(psi[["stage1_psi_intercept"]], ref_s1$estimate, tolerance = 1e-4)
  expect_equal(se[["stage0_psi_intercept"]], ref_s0$se, tolerance = 0.003)
  expect_equal(se[["stage1_psi_intercept"]], ref_s1$se, tolerance = 0.003)
})


test_that("longitudinal SNM + TV-EM + TF matches delicatessen", {
  ref <- read.csv(test_path(
    "fixtures",
    "snm_longitudinal_tv_em_delicatessen.csv"
  ))

  wide <- read.csv(test_path(
    "fixtures",
    "snm_longitudinal_tv_em_fixture.csv"
  ))
  d <- data.table::rbindlist(lapply(seq_len(nrow(wide)), function(i) {
    data.table::data.table(
      id = c(wide$id[i], wide$id[i]),
      time = c(0L, 1L),
      Y = c(NA_real_, wide$Y[i]),
      A = c(wide$A0[i], wide$A1[i]),
      L = c(wide$L0[i], wide$L1[i]),
      M = c(wide$M0[i], wide$M1[i])
    )
  }))

  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders_outcome = ~ A:M,
    confounders_tv = ~ L + M,
    id = "id",
    time = "time",
    type = "longitudinal",
    family = "gaussian",
    estimator = "snm",
    propensity_model_fn = stats::glm,
    treatment_free = ~L
  )
  res <- contrast(fit, ci_method = "sandwich")
  psi <- res$estimates$estimate
  names(psi) <- res$estimates$parameter
  se <- res$estimates$se
  names(se) <- res$estimates$parameter

  ref_s1_int <- ref[ref$parameter == "stage1_psi_intercept", ]
  ref_s1_M <- ref[ref$parameter == "stage1_psi_M", ]
  ref_s0_int <- ref[ref$parameter == "stage0_psi_intercept", ]
  ref_s0_M <- ref[ref$parameter == "stage0_psi_M", ]

  expect_equal(
    psi[["stage1_psi_intercept"]],
    ref_s1_int$estimate,
    tolerance = 1e-4
  )
  expect_equal(
    psi[["stage1_psi_M"]],
    ref_s1_M$estimate,
    tolerance = 1e-4
  )
  expect_equal(
    psi[["stage0_psi_intercept"]],
    ref_s0_int$estimate,
    tolerance = 1e-4
  )
  expect_equal(
    psi[["stage0_psi_M"]],
    ref_s0_M$estimate,
    tolerance = 1e-4
  )

  # SE tolerance wider (1%) for TV-EM + TF because the backward-transform
  # chain amplifies the bread/meat aggregation difference between R's
  # cluster-robust sandwich and delicatessen's stacked sandwich
  expect_equal(
    se[["stage1_psi_intercept"]],
    ref_s1_int$se,
    tolerance = 0.01
  )
  expect_equal(
    se[["stage1_psi_M"]],
    ref_s1_M$se,
    tolerance = 0.01
  )
  expect_equal(
    se[["stage0_psi_intercept"]],
    ref_s0_int$se,
    tolerance = 0.01
  )
  expect_equal(
    se[["stage0_psi_M"]],
    ref_s0_M$se,
    tolerance = 0.01
  )
})


# --- history = 0 tests -------------------------------------------------------
# history = 0 means no lag columns in the per-period treatment model.

test_that("longitudinal IPW works with history = 0 (no lags)", {
  dgp <- simulate_snm_longitudinal(n = 5000, seed = 201)

  fit <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders = ~1,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    estimator = "ipw",
    propensity_model_fn = stats::glm,
    history = 0
  )
  res <- contrast(
    fit,
    interventions = list(always = static(1), never = static(0)),
    type = "difference",
    reference = "never"
  )

  # IPW should still recover a reasonable ATE even without lags.
  # The model is slightly misspecified (A1 depends on A0 in truth),
  # but the bias is small.
  expect_equal(res$contrasts$estimate, 6.15, tolerance = 0.05)
})


test_that("longitudinal ICE (gcomp) works with history = 0 (no lags)", {
  dgp <- simulate_snm_longitudinal(n = 5000, seed = 202)

  fit <- causat(
    dgp$data,
    outcome = "Y",
    treatment = "A",
    confounders = ~1,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    estimator = "gcomp",
    history = 0
  )
  res <- contrast(
    fit,
    interventions = list(always = static(1), never = static(0)),
    type = "difference",
    reference = "never"
  )

  # ICE with history=0: outcome model at t=1 doesn't condition on lagged
  # treatment A0, so the indirect A0 -> L1 -> Y path is partially
  # captured but the direct A0 -> Y path at t=1 is not conditioned on.
  # This produces a valid but different estimate than the full-history
  # model. The key assertion: it runs without error and produces finite
  # results.
  expect_true(is.finite(res$contrasts$estimate))
  expect_true(res$contrasts$estimate > 0)
})
