simulate_binary_continuous <- function(n = 2000, seed = 42) {
  set.seed(seed)
  L <- rnorm(n)
  ps <- plogis(0.5 * L)
  A <- rbinom(n, 1, ps)
  Y <- 2 + 3 * A + 1.5 * L + rnorm(n)
  data.frame(Y = Y, A = A, L = L)
}

# ============================================================
# BASIC STRUCTURE
# ============================================================

test_that("diagnose() rejects non-causatr_fit input", {
  expect_error(diagnose("not_a_fit"), "causatr_fit")
})


test_that("diagnose() rejects longitudinal fits with a clear error", {
  # Coverage gap: running diagnose() on a longitudinal fit used to
  # produce a degenerate single-row positivity table (logistic on all
  # time points pooled) and then crash inside cobalt's print method.
  # Reject explicitly until per-period diagnostics land in a future
  # phase. Snapshot locks the message so a future refactor can't
  # silently downgrade the error or accidentally re-enable the broken
  # code path.
  set.seed(1)
  n <- 60
  L0 <- stats::rnorm(n)
  A0 <- stats::rbinom(n, 1, plogis(0.4 * L0))
  L1 <- L0 + 0.3 * A0 + stats::rnorm(n)
  A1 <- stats::rbinom(n, 1, plogis(0.4 * L1 + 0.5 * A0))
  Y <- 1 + A0 + A1 + 0.5 * L1 + stats::rnorm(n)
  long <- data.table::data.table(
    id = rep(seq_len(n), each = 2),
    time = rep(0:1, times = n),
    A = as.numeric(rbind(A0, A1)),
    L = as.numeric(rbind(L0, L1)),
    Y = rep(Y, each = 2)
  )
  fit <- causat(
    long,
    outcome = "Y",
    treatment = "A",
    confounders = ~1,
    confounders_tv = ~L,
    id = "id",
    time = "time"
  )
  expect_snapshot(error = TRUE, diagnose(fit))
})

test_that("diagnose() returns causatr_diag for gcomp", {
  df <- simulate_binary_continuous(n = 500)
  fit <- causat(df, outcome = "Y", treatment = "A", confounders = ~L)
  diag <- diagnose(fit)

  expect_s3_class(diag, "causatr_diag")
  expect_equal(diag$estimator, "gcomp")
  expect_false(is.null(diag$positivity))
  expect_false(is.null(diag$balance))
  expect_null(diag$weights)
  expect_null(diag$match_quality)
})

test_that("diagnose() returns causatr_diag for IPW", {
  df <- simulate_binary_continuous(n = 500)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw"
  )
  diag <- diagnose(fit)

  expect_s3_class(diag, "causatr_diag")
  expect_equal(diag$estimator, "ipw")
  expect_false(is.null(diag$positivity))
  expect_false(is.null(diag$balance))
  expect_false(is.null(diag$weights))
  expect_null(diag$match_quality)
})

test_that("diagnose() returns causatr_diag for matching", {
  df <- simulate_binary_continuous(n = 500)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "matching",
    estimand = "ATT"
  )
  diag <- diagnose(fit)

  expect_s3_class(diag, "causatr_diag")
  expect_equal(diag$estimator, "matching")
  expect_false(is.null(diag$positivity))
  expect_false(is.null(diag$balance))
  expect_null(diag$weights)
  expect_false(is.null(diag$match_quality))
})

# ============================================================
# POSITIVITY
# ============================================================

test_that("positivity table has expected rows for binary treatment", {
  df <- simulate_binary_continuous(n = 500)
  fit <- causat(df, outcome = "Y", treatment = "A", confounders = ~L)
  diag <- diagnose(fit)

  pos <- diag$positivity
  expect_s3_class(pos, "data.table")
  expect_true("statistic" %in% names(pos))
  expect_true("value" %in% names(pos))
  expect_true("min" %in% pos$statistic)
  expect_true("max" %in% pos$statistic)
  expect_true("n_violations" %in% pos$statistic)
})

test_that("positivity uses propensity scores from weightit for IPW", {
  df <- simulate_binary_continuous(n = 500)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw"
  )
  diag <- diagnose(fit)

  pos <- diag$positivity
  ps_min <- pos$value[pos$statistic == "min"]
  ps_max <- pos$value[pos$statistic == "max"]
  expect_gt(ps_min, 0)
  expect_lt(ps_max, 1)
})

# ============================================================
# BALANCE (cobalt)
# ============================================================

test_that("balance uses cobalt::bal.tab for IPW when available", {
  skip_if_not_installed("cobalt")
  df <- simulate_binary_continuous(n = 500)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw"
  )
  diag <- diagnose(fit)

  expect_s3_class(diag$balance, "bal.tab")
})

test_that("balance uses cobalt::bal.tab for matching when available", {
  skip_if_not_installed("cobalt")
  df <- simulate_binary_continuous(n = 500)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "matching",
    estimand = "ATT"
  )
  diag <- diagnose(fit)

  expect_s3_class(diag$balance, "bal.tab")
})

test_that("balance uses cobalt::bal.tab for gcomp (unadjusted)", {
  skip_if_not_installed("cobalt")
  df <- simulate_binary_continuous(n = 500)
  fit <- causat(df, outcome = "Y", treatment = "A", confounders = ~L)
  diag <- diagnose(fit)

  expect_s3_class(diag$balance, "bal.tab")
})

# ============================================================
# WEIGHT DISTRIBUTION (IPW)
# ============================================================

test_that("weight summary has treated/control/overall for IPW", {
  df <- simulate_binary_continuous(n = 500)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw"
  )
  diag <- diagnose(fit)

  w <- diag$weights
  expect_s3_class(w, "data.table")
  expect_equal(nrow(w), 3L)
  expect_true(all(c("treated", "control", "overall") %in% w$group))
  expect_true("ess" %in% names(w))
  expect_true(all(w$ess > 0))
})

# ============================================================
# MATCH QUALITY (matching)
# ============================================================

test_that("match_quality has expected fields for matching", {
  df <- simulate_binary_continuous(n = 500)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "matching",
    estimand = "ATT"
  )
  diag <- diagnose(fit)

  mq <- diag$match_quality
  expect_s3_class(mq, "data.table")
  expect_true("n_total" %in% mq$statistic)
  expect_true("n_matched" %in% mq$statistic)
  expect_true("n_discarded" %in% mq$statistic)
  expect_true("pct_retained" %in% mq$statistic)
})

# ============================================================
# PRINT AND SUMMARY
# ============================================================

test_that("print.causatr_diag outputs method info", {
  df <- simulate_binary_continuous(n = 500)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw"
  )
  diag <- diagnose(fit)

  expect_output(print(diag), "causatr_diag")
  expect_output(print(diag), "ipw")
})

test_that("summary.causatr_diag works", {
  df <- simulate_binary_continuous(n = 500)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw"
  )
  diag <- diagnose(fit)

  expect_output(summary(diag), "causatr_diag")
})

# ============================================================
# PLOT (cobalt love.plot)
# ============================================================

test_that("plot.causatr_diag produces a love plot for IPW", {
  skip_if_not_installed("cobalt")
  skip_if_not_installed("ggplot2")
  df <- simulate_binary_continuous(n = 500)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw"
  )
  diag <- diagnose(fit)

  p <- plot(diag)
  expect_true(inherits(p, "gg") || inherits(p, "gtable"))
})

test_that("plot.causatr_diag produces a love plot for matching", {
  skip_if_not_installed("cobalt")
  skip_if_not_installed("ggplot2")
  df <- simulate_binary_continuous(n = 500)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "matching",
    estimand = "ATT"
  )
  diag <- diagnose(fit)

  p <- plot(diag)
  expect_true(inherits(p, "gg") || inherits(p, "gtable"))
})

test_that("plot.causatr_diag errors for gcomp (no love plot available)", {
  df <- simulate_binary_continuous(n = 500)
  fit <- causat(df, outcome = "Y", treatment = "A", confounders = ~L)
  diag <- diagnose(fit)

  expect_error(plot(diag), "IPW or matching")
})

# ============================================================
# CUSTOM PARAMETERS
# ============================================================

test_that("diagnose() respects custom thresholds", {
  skip_if_not_installed("cobalt")
  df <- simulate_binary_continuous(n = 500)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw"
  )
  diag <- diagnose(fit, thresholds = c(m = 0.05, v = 1.5))

  expect_s3_class(diag$balance, "bal.tab")
})

test_that("diagnose() respects custom ps_bounds", {
  df <- simulate_binary_continuous(n = 500)
  fit <- causat(df, outcome = "Y", treatment = "A", confounders = ~L)

  diag_wide <- diagnose(fit, ps_bounds = c(0.01, 0.99))
  diag_narrow <- diagnose(fit, ps_bounds = c(0.1, 0.9))

  n_viol_wide <- diag_wide$positivity$value[
    diag_wide$positivity$statistic == "n_violations"
  ]
  n_viol_narrow <- diag_narrow$positivity$value[
    diag_narrow$positivity$statistic == "n_violations"
  ]
  expect_lte(n_viol_wide, n_viol_narrow)
})


# Fifth-round critical review, R1: diagnose() must exclude censored rows
# when computing positivity, balance, and simple balance for gcomp fits
# with a censoring column. Before the fix, diagnose() used all non-NA-
# outcome rows, ignoring censoring. The main gcomp pipeline excludes
# censored rows via get_fit_rows(data, outcome, censoring), so
# diagnostics computed on a different sample were misleading.
# Repro script: /tmp/causatr_repro_r1_diagnose_censoring.R
test_that("diagnose() excludes censored rows for gcomp with censoring", {
  set.seed(42)
  n <- 200
  d <- data.table::data.table(
    Y = rnorm(n),
    A = rbinom(n, 1, 0.5),
    L = rnorm(n),
    C = c(rep(0L, 150), rep(1L, 50))
  )

  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "gcomp",
    censoring = "C"
  )

  # The main pipeline should have fit on 150 rows (excluding C == 1)
  expect_equal(sum(fit$details$fit_rows), 150L)

  # Positivity diagnostics should also use 150 rows
  pos <- compute_positivity(fit, ps_bounds = c(0.025, 0.975))
  # The PS model is fit on the same 150 rows, so the number of
  # observations in the positivity table equals the fit-row count.
  ps_formula <- build_ps_formula(fit$confounders, fit$treatment)
  fit_rows <- get_fit_rows(fit$data, fit$outcome, fit$censoring)
  ps_model <- stats::glm(
    ps_formula,
    data = fit$data[fit_rows],
    family = stats::binomial()
  )
  expect_equal(length(stats::fitted(ps_model)), 150L)

  # Balance should also use 150 rows
  bal <- compute_balance_simple(fit)
  expect_false(is.null(bal))
  # The balance table should be based on uncensored rows only.
  # We check that the SMD is computed on the right subset by
  # verifying directly:
  d_uncensored <- d[d$C == 0]
  smd_manual <- (mean(d_uncensored$L[d_uncensored$A == 1]) -
    mean(d_uncensored$L[d_uncensored$A == 0])) /
    sqrt(
      (var(d_uncensored$L[d_uncensored$A == 1]) +
        var(d_uncensored$L[d_uncensored$A == 0])) /
        2
    )
  expect_equal(bal$smd[bal$variable == "L"], smd_manual, tolerance = 1e-10)
})

# compute_positivity() short-circuits ------------------------------------
#
# The PS-quantile table is only meaningful for a single binary
# treatment. Multivariate, non-{0,1}, non-bernoulli IPW, and matching
# fits with a missing distance vector all need to return NULL silently
# so `diagnose()` can fall back to its other diagnostics.

test_that("compute_positivity() returns NULL for multivariate treatment", {
  fake_fit <- list(
    treatment = c("A1", "A2"),
    estimator = "gcomp",
    data = data.table::data.table(),
    confounders = ~L
  )
  expect_null(compute_positivity(fake_fit, ps_bounds = c(0.05, 0.95)))
})

test_that("compute_positivity() returns NULL for non-{0,1} treatment values", {
  d <- data.table::data.table(A = c(2, 5, 7), L = rnorm(3), Y = rnorm(3))
  fake_fit <- list(
    treatment = "A",
    estimator = "gcomp",
    data = d,
    confounders = ~L,
    outcome = "Y"
  )
  expect_null(compute_positivity(fake_fit, ps_bounds = c(0.05, 0.95)))
})

test_that("compute_positivity() returns NULL for IPW with non-bernoulli treatment_model", {
  # Build a fake fit shaped like an IPW result whose treatment_model
  # has family != "bernoulli". The early NULL return at line ~171
  # protects diagnose() from trying to threshold a continuous PS.
  d <- data.table::data.table(A = c(0, 1, 0, 1), L = rnorm(4), Y = rnorm(4))
  fake_tm <- list(family = "gaussian")
  fake_fit <- list(
    treatment = "A",
    estimator = "ipw",
    data = d,
    confounders = ~L,
    outcome = "Y",
    details = list(treatment_model = fake_tm)
  )
  expect_null(compute_positivity(fake_fit, ps_bounds = c(0.05, 0.95)))
})

# ============================================================
# CHUNK 11A: PER-INTERVENTION DISPATCH
# ============================================================

# The chunk 11a rewrite of diagnose() introduced a per-intervention
# panel architecture. Default `diagnose(fit)` produces a single panel
# keyed `obs` using the observed-treatment Horvitz-Thompson view (the
# pre-chunk-11a default behaviour, preserved here as a regression
# anchor). User-supplied `interventions =` produces one panel per key,
# each with a per-intervention density-ratio weight summary.

test_that("diagnose() default returns a single `obs` panel matching the legacy view", {
  df <- simulate_binary_continuous(n = 500)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw"
  )
  diag <- diagnose(fit)

  expect_s3_class(diag, "causatr_diag")
  expect_named(diag$per_intervention, "obs")
  expect_equal(diag$interventions, "obs")
  # Backward-compat top-level slots must point to the obs panel.
  expect_identical(diag$positivity, diag$per_intervention$obs$positivity)
  expect_identical(diag$balance, diag$per_intervention$obs$balance)
  expect_identical(diag$weights, diag$per_intervention$obs$weights)

  # The default obs panel for binary IPW reports the observed-arm
  # HT weights (treated -> 1/p, control -> 1/(1-p)). Reconstruct
  # those manually and check the weight summary matches.
  ps <- as.numeric(stats::predict(
    fit$details$propensity_model,
    newdata = fit$data[fit$details$fit_rows],
    type = "response"
  ))
  a <- fit$data[fit$details$fit_rows][[fit$treatment]]
  w_obs <- ifelse(a == 1, 1 / ps, 1 / (1 - ps))
  ess_treated_manual <- sum(w_obs[a == 1])^2 / sum(w_obs[a == 1]^2)
  panel_w <- diag$per_intervention$obs$weights
  expect_equal(
    panel_w$ess[panel_w$group == "treated"],
    ess_treated_manual,
    tolerance = 1e-10
  )
})

test_that("diagnose() with interventions = list(a1, a0) produces two panels", {
  df <- simulate_binary_continuous(n = 500)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw"
  )
  diag <- diagnose(
    fit,
    interventions = list(a1 = static(1), a0 = static(0))
  )

  expect_named(diag$per_intervention, c("a1", "a0"))
  expect_equal(diag$interventions, c("a1", "a0"))

  # Each panel carries the (intervention-independent) shared
  # positivity / balance plus a per-intervention weight table.
  expect_false(is.null(diag$per_intervention$a1$positivity))
  expect_false(is.null(diag$per_intervention$a0$positivity))
  expect_false(is.null(diag$per_intervention$a1$weights))
  expect_false(is.null(diag$per_intervention$a0$weights))

  # Manual reconstruction: the a1 panel's treated-arm weights are
  # `1{A=1}/p` (=> mean(w[treated]) = mean(1/p[treated])); the a0
  # panel's control-arm weights are `1{A=0}/(1-p)`.
  ps <- as.numeric(stats::predict(
    fit$details$propensity_model,
    newdata = fit$data[fit$details$fit_rows],
    type = "response"
  ))
  a <- fit$data[fit$details$fit_rows][[fit$treatment]]
  w_a1 <- as.numeric(a == 1) / ps
  w_a0 <- as.numeric(a == 0) / (1 - ps)

  panel_a1 <- diag$per_intervention$a1$weights
  panel_a0 <- diag$per_intervention$a0$weights

  ess <- function(w) {
    s2 <- sum(w^2)
    if (s2 == 0) 0 else sum(w)^2 / s2
  }
  expect_equal(
    panel_a1$ess[panel_a1$group == "treated"],
    ess(w_a1[a == 1]),
    tolerance = 1e-10
  )
  expect_equal(
    panel_a0$ess[panel_a0$group == "control"],
    ess(w_a0[a == 0]),
    tolerance = 1e-10
  )
  # The a1 panel's control arm has all-zero weights (HT indicator
  # zeros the off-arm), so ESS = 0 by the divide-by-zero guard.
  expect_equal(panel_a1$ess[panel_a1$group == "control"], 0)
  expect_equal(panel_a0$ess[panel_a0$group == "treated"], 0)
})

test_that("diagnose() with `interventions = list(natural = NULL)` reports unit weights", {
  # Distinct from the auto-injected `obs` panel (which uses the
  # observed-arm HT view): a literal NULL intervention requests
  # `compute_density_ratio_weights(tm, data, NULL)` which returns
  # the all-1 natural-course weight vector. Useful as a sanity
  # check that the sentinel branch is wired correctly.
  df <- simulate_binary_continuous(n = 500)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw"
  )
  diag <- diagnose(fit, interventions = list(natural = NULL))

  panel <- diag$per_intervention$natural$weights
  expect_equal(panel$mean[panel$group == "overall"], 1, tolerance = 1e-10)
  expect_equal(panel$min[panel$group == "overall"], 1, tolerance = 1e-10)
  expect_equal(panel$max[panel$group == "overall"], 1, tolerance = 1e-10)
})

test_that("diagnose() rejects un-named `interventions =` lists", {
  df <- simulate_binary_continuous(n = 500)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw"
  )
  expect_error(
    diagnose(fit, interventions = list(static(1), static(0))),
    "must be named"
  )
})

test_that("print.causatr_diag handles multi-panel output", {
  df <- simulate_binary_continuous(n = 500)
  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw"
  )
  diag <- diagnose(
    fit,
    interventions = list(a1 = static(1), a0 = static(0))
  )
  out <- capture.output(print(diag))
  expect_true(any(grepl("Intervention: a1", out)))
  expect_true(any(grepl("Intervention: a0", out)))
})

# ============================================================
# CHUNK 11A: PENDING REJECTION GATES
# ============================================================

# Each unsupported combination must abort with its `causatr_diag_*_pending`
# class so future chunks can lift the rejection by removing the gate.

test_that("diagnose() aborts with causatr_diag_em_pending when by= is supplied", {
  df <- simulate_binary_continuous(n = 200)
  fit <- causat(df, outcome = "Y", treatment = "A", confounders = ~L)
  expect_error(
    diagnose(fit, by = "L"),
    class = "causatr_diag_em_pending"
  )
})

test_that("diagnose() aborts with causatr_diag_estimand_pending for IPW under ATT", {
  set.seed(7)
  n <- 600
  L <- rnorm(n)
  A <- rbinom(n, 1, plogis(0.5 * L))
  Y <- 2 + 3 * A + 1.5 * L + rnorm(n)
  d <- data.table::data.table(Y = Y, A = A, L = L)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw",
    estimand = "ATT"
  )
  expect_error(
    diagnose(fit),
    class = "causatr_diag_estimand_pending"
  )
})

test_that("diagnose() aborts with causatr_diag_estimand_pending for IPW under ATC", {
  set.seed(7)
  n <- 600
  L <- rnorm(n)
  A <- rbinom(n, 1, plogis(0.5 * L))
  Y <- 2 + 3 * A + 1.5 * L + rnorm(n)
  d <- data.table::data.table(Y = Y, A = A, L = L)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw",
    estimand = "ATC"
  )
  expect_error(
    diagnose(fit),
    class = "causatr_diag_estimand_pending"
  )
})

test_that("diagnose() aborts with causatr_diag_continuous_pending for IPW + continuous", {
  set.seed(11)
  n <- 400
  L <- rnorm(n)
  A <- 0.5 * L + rnorm(n)
  Y <- 1 + 2 * A + L + rnorm(n)
  d <- data.table::data.table(Y = Y, A = A, L = L)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw"
  )
  expect_error(
    diagnose(fit),
    class = "causatr_diag_continuous_pending"
  )
})

test_that("diagnose() aborts with causatr_diag_categorical_pending for IPW + categorical", {
  skip_if_not_installed("nnet")
  set.seed(13)
  n <- 600
  L <- rnorm(n)
  pr <- exp(cbind(0, 0.4 * L, -0.3 * L))
  pr <- pr / rowSums(pr)
  A_idx <- apply(pr, 1, function(p) sample(seq_along(p), 1, prob = p))
  A <- factor(c("low", "med", "hi")[A_idx])
  Y <- 1 + 2 * (A == "med") + 3 * (A == "hi") + L + rnorm(n)
  d <- data.table::data.table(Y = Y, A = A, L = L)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw"
  )
  expect_error(
    diagnose(fit),
    class = "causatr_diag_categorical_pending"
  )
})

test_that("diagnose() aborts with causatr_diag_count_pending for IPW + count", {
  skip_if_not_installed("MASS")
  set.seed(17)
  n <- 600
  L <- rnorm(n)
  A <- rpois(n, lambda = exp(0.3 + 0.2 * L))
  Y <- 1 + 2 * A + L + rnorm(n)
  d <- data.table::data.table(Y = Y, A = A, L = L)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw",
    propensity_family = "poisson"
  )
  expect_error(
    diagnose(fit),
    class = "causatr_diag_count_pending"
  )
})

test_that("diagnose() aborts with causatr_diag_multivariate_pending for multivariate IPW", {
  set.seed(19)
  n <- 500
  L <- rnorm(n)
  A1 <- rbinom(n, 1, plogis(0.4 * L))
  A2 <- rbinom(n, 1, plogis(0.3 * L + 0.5 * A1))
  Y <- 1 + A1 + A2 + L + rnorm(n)
  d <- data.table::data.table(Y = Y, A1 = A1, A2 = A2, L = L)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = c("A1", "A2"),
    confounders = ~L,
    estimator = "ipw"
  )
  expect_error(
    diagnose(fit),
    class = "causatr_diag_multivariate_pending"
  )
})
