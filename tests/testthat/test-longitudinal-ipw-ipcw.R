# Longitudinal IPW under IPCW (missing Y): the id-level analytic sandwich now
# propagates the per-period censoring model's estimation uncertainty through the
# MSM (the gamma -> beta cross-term, `compute_ipw_ipcw_correction_longitudinal()`).
# Unlike the doubly-robust AIPW case -- where double-robust orthogonality makes
# the censoring cross-term ~0.03% -- this term is large and load-bearing for IPW
# (treating the estimated IPCW weights as fixed over-states the treated arm SE).
#
# The headline oracle is the delicatessen M-estimation stack carrying the
# censoring gamma block (fixtures/python/longitudinal_ipw_ipcw_delicatessen.py),
# which the chunk-1 (AIPW) work reserved for this chunk. The companion
# "known-weights" sandwich in the same fixture (gamma held fixed) pins the
# efficiency gain the cross-term recovers (Robins-Rotnitzky-Zhao 1994).
#
# The estimator itself was also fixed here: built-in IPCW no longer weights the
# treatment-density models by the censoring weight (the standard IPW+IPCW
# construction estimates the treatment and censoring models separately and
# multiplies their weights -- Hernan & Robins 2020, Ch. 12.6 & 17). The
# point-estimate and propensity anchors below pin that.

ipcw_oracle <- function() {
  utils::read.csv(test_path(
    "fixtures",
    "python",
    "longitudinal_ipw_ipcw_delicatessen_results.csv"
  ))
}

ipcw_fit <- function(n = 2000, seed = 2025) {
  d <- data.table::as.data.table(
    simulate_longitudinal_mar_outcome(n = n, seed = seed)
  )
  causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    estimator = "ipw",
    type = "longitudinal",
    id = "id",
    time = "time",
    censoring = "C",
    ipcw = TRUE
  )
}

# --- headline: tight delicatessen oracle match ------------------------------

# The data is `simulate_longitudinal_mar_outcome(n = 2000, seed = 2025)`, the
# same call the fixture generator reshapes to wide for the Python oracle, so the
# two see byte-identical data. Point estimates are the Hajek IPW marginal means
# (match to ~1e-6); the analytic sandwich solves the identical stacked
# M-estimation system delicatessen does, so the SEs agree to ~1e-4.
test_that("longitudinal IPW + IPCW sandwich matches delicatessen", {
  skip_if(
    !file.exists(test_path(
      "fixtures",
      "python",
      "longitudinal_ipw_ipcw_delicatessen_results.csv"
    ))
  )
  r <- ipcw_oracle()
  fit <- ipcw_fit()
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    type = "difference",
    ci_method = "sandwich"
  )
  mu_1 <- res$estimates$estimate[res$estimates$intervention == "a1"]
  mu_0 <- res$estimates$estimate[res$estimates$intervention == "a0"]
  se_1 <- res$estimates$se[res$estimates$intervention == "a1"]
  se_0 <- res$estimates$se[res$estimates$intervention == "a0"]

  expect_equal(mu_1, r$mu_1, tolerance = 1e-6)
  expect_equal(mu_0, r$mu_0, tolerance = 1e-6)
  expect_equal(res$contrasts$estimate[1], r$ate, tolerance = 1e-6)
  expect_equal(se_1, r$se_1, tolerance = 1e-4)
  expect_equal(se_0, r$se_0, tolerance = 1e-4)
  expect_equal(res$contrasts$se[1], r$se_ate, tolerance = 1e-4)
})

# --- the censoring cross-term is load-bearing (two-sided vs the oracle) ------

# Disabling the gamma block in the variance engine (treating the IPCW weights as
# KNOWN) reproduces delicatessen's "known-weights" sandwich. The full sandwich
# is strictly smaller -- the gamma -> beta correction recovers the censoring
# estimation efficiency gain. Both the full and the known-weights SEs are pinned
# to the oracle to ~1e-4, so the gap is an exact magnitude, not a one-sided
# "it changed" check. For IPW the treated-arm gap is ~5% (vs ~0.03% for AIPW).
test_that("longitudinal IPW + IPCW censoring cross-term is load-bearing", {
  skip_if(
    !file.exists(test_path(
      "fixtures",
      "python",
      "longitudinal_ipw_ipcw_delicatessen_results.csv"
    ))
  )
  r <- ipcw_oracle()
  ints <- list(a1 = static(1), a0 = static(0))
  fit <- ipcw_fit()

  res_full <- contrast(
    fit,
    interventions = ints,
    reference = "a0",
    type = "difference",
    ci_method = "sandwich"
  )

  # Known-weights sandwich: drop the gamma correction by flagging the fit
  # non-IPCW for the variance engine only. The point path reads
  # `details$weights` (still the folded IPCW weights), so the marginal means and
  # MSM are identical -- only the censoring correction is skipped.
  fit_known <- fit
  fit_known$details$ipcw <- FALSE
  res_known <- contrast(
    fit_known,
    interventions = ints,
    reference = "a0",
    type = "difference",
    ci_method = "sandwich"
  )

  se1_full <- res_full$estimates$se[res_full$estimates$intervention == "a1"]
  se0_full <- res_full$estimates$se[res_full$estimates$intervention == "a0"]
  se1_known <- res_known$estimates$se[res_known$estimates$intervention == "a1"]
  se0_known <- res_known$estimates$se[res_known$estimates$intervention == "a0"]

  # Point estimates are identical with/without the variance-only flag.
  expect_equal(res_full$estimates$estimate, res_known$estimates$estimate)

  # Full SE matches the oracle's full sandwich; known-weights SE matches the
  # oracle's gamma-fixed sandwich -- both tight.
  expect_equal(se1_known, r$se_1_known, tolerance = 1e-4)
  expect_equal(se0_known, r$se_0_known, tolerance = 1e-4)
  expect_equal(
    res_known$contrasts$se[1],
    r$se_ate_known,
    tolerance = 1e-4
  )

  # Load-bearing: the gamma correction strictly reduces the SE, materially so on
  # the treated arm (where censoring depends on treatment).
  expect_lt(se1_full, se1_known)
  expect_equal(se1_full / se1_known, r$se_1 / r$se_1_known, tolerance = 1e-3)
  expect_lt(r$se_1 / r$se_1_known, 0.96) # >4% efficiency gain on the treated arm
})

# --- estimator-fix anchor: propensity is NOT IPCW-weighted ------------------

# Built-in IPCW must reweight the OUTCOME side (MSM) only. The per-period
# treatment-density models are ordinary regressions on all observed rows -- not
# IPCW-weighted, and censored rows are not dropped. The fitted prior weights are
# therefore all 1 (no external weights in this DGP), and the period-1 propensity
# equals an unweighted GLM on every period-1 row. Before the fix the prior
# weights were the per-row IPCW weights (0 on censored rows, up to ~10).
test_that("longitudinal IPW + IPCW leaves the propensity unweighted", {
  fit <- ipcw_fit()
  tp <- fit$details$time_points
  n_periods <- length(tp)

  for (k in seq_len(n_periods)) {
    tm <- fit$details$treatment_models_by_time[[k]]
    pm <- tm$model
    # No IPCW reweighting: every fitted prior weight is exactly 1.
    expect_equal(unname(pm$prior.weights), rep(1, length(pm$prior.weights)))
    # All period-k rows enter the fit (censored rows are not dropped).
    expect_equal(length(pm$y), sum(fit$data[[fit$time]] == tp[k]))
  }

  # Period-1 propensity equals an unweighted GLM on all period-1 rows.
  data1 <- fit$details$fit_data_by_time[[2]]
  ps_form <- remove_response(fit$details$per_period_formula[[2]])
  ref <- stats::glm(
    stats::reformulate(attr(stats::terms(ps_form), "term.labels"), "A"),
    family = stats::binomial(),
    data = data1
  )
  expect_equal(
    unname(coef(fit$details$treatment_models_by_time[[2]]$model)),
    unname(coef(ref)),
    tolerance = 1e-8
  )
})

# --- non-IPCW longitudinal IPW is unchanged by the propensity-weight split ---

# The propensity-weight separation only diverges from `weights` under IPCW; a
# plain longitudinal IPW fit (external weights, no censoring) is byte-identical.
test_that("longitudinal IPW without IPCW is unaffected", {
  d <- data.table::as.data.table(make_linear_scm(
    n = 1500,
    n_times = 2,
    seed = 9
  ))
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    estimator = "ipw",
    type = "longitudinal",
    id = "id",
    time = "time"
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    type = "difference",
    ci_method = "sandwich"
  )
  # True ATE (always vs never, 2 periods) = 5.
  expect_equal(res$contrasts$estimate[1], 5, tolerance = 0.25)
  expect_true(all(is.finite(res$estimates$se)))
})

# --- Monte-Carlo calibration: sandwich SE tracks the empirical sampling SD ---

# The variance-targeted oracle: under repeated sampling the per-arm IPCW sandwich
# SE must track the empirical SD of the point estimate. This detects a missing or
# mis-scaled censoring cross-term independently of the delicatessen fixture.
test_that("longitudinal IPW + IPCW sandwich SE tracks the empirical SD", {
  skip_if_fast()
  ints <- list(a1 = static(1), a0 = static(0))
  R <- 400L
  est <- matrix(NA_real_, R, 2L)
  se <- matrix(NA_real_, R, 2L)
  for (rr in seq_len(R)) {
    fit <- tryCatch(
      ipcw_fit(n = 3000, seed = 50000 + rr),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      next
    }
    sw <- tryCatch(
      contrast(
        fit,
        interventions = ints,
        reference = "a0",
        type = "difference",
        ci_method = "sandwich"
      ),
      error = function(e) NULL
    )
    if (is.null(sw)) {
      next
    }
    est[rr, ] <- sw$estimates$estimate
    se[rr, ] <- sw$estimates$se
  }
  ok <- stats::complete.cases(est)
  ratio <- colMeans(se[ok, ]) / apply(est[ok, ], 2L, stats::sd)
  expect_equal(mean(ratio), 1, tolerance = 0.05)
  expect_true(all(ratio > 0.9 & ratio < 1.1))
})

# --- bootstrap parity --------------------------------------------------------

test_that("longitudinal IPW + IPCW sandwich SE agrees with the bootstrap", {
  skip_if_fast()
  fit <- ipcw_fit(n = 4000, seed = 71)
  ints <- list(a1 = static(1), a0 = static(0))
  rs <- contrast(
    fit,
    interventions = ints,
    reference = "a0",
    type = "difference",
    ci_method = "sandwich"
  )
  rb <- contrast(
    fit,
    interventions = ints,
    reference = "a0",
    type = "difference",
    ci_method = "bootstrap",
    n_boot = 600
  )
  expect_equal(rs$estimates$estimate, rb$estimates$estimate, tolerance = 1e-8)
  expect_equal(rs$contrasts$se / rb$contrasts$se, 1, tolerance = 0.1)
})
