# Longitudinal AIPW under IPCW (missing Y): the stacked-EE sandwich carries a
# per-period censoring (gamma) block so the numerical bread captures the
# censoring-model estimation cross-terms. These tests pin the new behaviour
# against independent ground truth (Monte-Carlo empirical SD, bootstrap), check
# the estimating equation is faithful with the gamma block, and document that
# the AIPW censoring cross-term is near-zero by double-robust orthogonality
# (so the SE barely moves, but the block is wired and correct).
#
# A delicatessen cross-language oracle is reserved for the longitudinal IPW
# IPCW sandwich, where the censoring cross-term is large (~15% of the SE) and
# load-bearing; for AIPW the cross-term is negligible (orthogonality) and the
# MC-empirical-SD + faithfulness + bootstrap checks below are the tight oracle.

# --- shared per-period censoring EE block -----------------------------------

test_that("make_ipcw_weight_fn_longitudinal reproduces weights + vanishing scores", {
  d <- simulate_longitudinal_mar_outcome(n = 600, seed = 31)
  dt <- data.table::as.data.table(d)
  fit <- causat(
    dt,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    estimator = "aipw",
    type = "longitudinal",
    id = "id",
    time = "time",
    censoring = "C",
    ipcw = TRUE
  )

  cb <- make_ipcw_weight_fn_longitudinal(fit)
  # The full-length weight closure at gamma_hat must reproduce the fitted
  # cumulative IPCW weights exactly (per-period stabilized weight per row).
  expect_equal(
    cb$weight_full_fn(cb$gamma_hat),
    fit$details$ipcw_weights,
    tolerance = 1e-10
  )
  # Each per-period censoring logistic score must vanish at gamma_hat (the MLE),
  # which is the faithfulness the stacked sandwich relies on.
  for (b in cb$blocks) {
    s <- glm_score_matrix(b$model, cb$gamma_hat[b$gamma_cols])
    expect_lt(max(abs(colSums(s))), 1e-6)
  }
})

# --- stacked EE faithfulness + gamma block wiring ---------------------------

# Rebuild the stacked estimating equation for one arm and confirm (a) the
# summed score vanishes at the fitted theta with the gamma block present, and
# (b) the bread's (beta, gamma) cross-block is non-zero -- i.e. the censoring
# cross-term is genuinely captured, not silently dropped.
test_that("AIPW IPCW stacked EE is faithful and the gamma block is wired", {
  d <- simulate_longitudinal_mar_outcome(n = 1200, seed = 7)
  dt <- data.table::as.data.table(d)
  fit <- causat(
    dt,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    estimator = "aipw",
    type = "longitudinal",
    id = "id",
    time = "time",
    censoring = "C",
    ipcw = TRUE
  )

  ir <- ice_aipw_iterate(fit, static(1))
  tp <- fit$details$time_points
  rows_first <- fit$data[[fit$time]] == tp[1]
  all_ids <- as.character(fit$data[rows_first][[fit$id]])
  n <- length(all_ids)
  id_to_idx <- stats::setNames(seq_len(n), all_ids)
  target <- rep(TRUE, n)
  w_t <- rep(1, n)
  built <- build_aipw_long_psi(fit, ir, all_ids, n, id_to_idx, target, w_t, n)

  expect_gt(built$total_gamma, 0L) # IPCW => a gamma block exists
  psi_hat <- built$psi_fn(built$theta_hat)
  col_scale <- sqrt(colSums(psi_hat^2))
  expect_lt(max(abs(colSums(psi_hat)) / pmax(col_scale, 1e-8)), 1e-4)

  B <- -numDeriv::jacobian(
    function(th) colSums(built$psi_fn(th)),
    built$theta_hat
  ) /
    n
  ga <- built$total_alpha + built$total_beta + seq_len(built$total_gamma)
  beta_rows <- built$total_alpha + seq_len(built$total_beta)
  # The outcome score depends on gamma through the IPCW weights => non-zero.
  expect_gt(max(abs(B[beta_rows, ga])), 1e-4)
})

# Without IPCW the gamma block is absent (total_gamma == 0), so the stacked EE
# is byte-identical to the pre-existing (validated) non-IPCW assembly.
test_that("AIPW non-IPCW stacked EE carries no gamma block", {
  d <- make_linear_scm(n = 500, n_times = 2, seed = 4)
  dt <- data.table::as.data.table(d)
  fit <- causat(
    dt,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    estimator = "aipw",
    type = "longitudinal",
    id = "id",
    time = "time"
  )
  ir <- ice_aipw_iterate(fit, static(1))
  tp <- fit$details$time_points
  rf <- fit$data[[fit$time]] == tp[1]
  ids <- as.character(fit$data[rf][[fit$id]])
  n <- length(ids)
  built <- build_aipw_long_psi(
    fit,
    ir,
    ids,
    n,
    stats::setNames(seq_len(n), ids),
    rep(TRUE, n),
    rep(1, n),
    n
  )
  expect_equal(built$total_gamma, 0L)
})

# --- double-robust orthogonality of the censoring cross-term ----------------

# The AIPW (ICE-AIPW) estimator is first-order insensitive to the censoring
# nuisance, so estimating gamma vs treating it as known changes the
# marginal-mean SE only negligibly -- yet the cross-term is non-zero and must be
# present for the SE to be exactly correct. Two-sided: the with-gamma SE equals
# the known-gamma SE to within 0.5%, confirming orthogonality while the wiring
# test above confirms the block is active.
test_that("AIPW IPCW censoring cross-term is orthogonal (tiny but present)", {
  d <- simulate_longitudinal_mar_outcome(n = 1500, seed = 7)
  dt <- data.table::as.data.table(d)
  fit <- causat(
    dt,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    estimator = "aipw",
    type = "longitudinal",
    id = "id",
    time = "time",
    censoring = "C",
    ipcw = TRUE
  )
  ir <- ice_aipw_iterate(fit, static(1))
  tp <- fit$details$time_points
  rf <- fit$data[[fit$time]] == tp[1]
  ids <- as.character(fit$data[rf][[fit$id]])
  n <- length(ids)
  built <- build_aipw_long_psi(
    fit,
    ir,
    ids,
    n,
    stats::setNames(seq_len(n), ids),
    rep(TRUE, n),
    rep(1, n),
    n
  )
  psi_hat <- built$psi_fn(built$theta_hat)
  B <- -numDeriv::jacobian(
    function(th) colSums(built$psi_fn(th)),
    built$theta_hat
  ) /
    n
  mu <- built$mu_index
  se_of <- function(Bm, psi, mu_i) {
    em <- numeric(ncol(Bm))
    em[mu_i] <- 1
    as.numeric(sqrt(sum((psi %*% solve(t(Bm), em))^2))) / nrow(psi)
  }
  se_full <- se_of(B, psi_hat, mu)
  ga <- built$total_alpha + built$total_beta + seq_len(built$total_gamma)
  keep <- setdiff(seq_along(built$theta_hat), ga)
  se_known <- se_of(B[keep, keep], psi_hat[, keep], match(mu, keep))
  expect_equal(se_full, se_known, tolerance = 0.005)
})

# --- Monte-Carlo calibration (the test that detects a missing cross-term) ----

# Under repeated sampling, the per-arm IPCW marginal-mean sandwich SE must track
# the empirical sampling SD. This is the variance-targeted oracle: it isolates
# the SE without conflating it with finite-sample point bias.
test_that("AIPW IPCW sandwich SE tracks the empirical sampling SD", {
  skip_if_fast()
  ints <- list(a1 = static(1), a0 = static(0))
  R <- 300L
  est <- matrix(NA_real_, R, 2L)
  se <- matrix(NA_real_, R, 2L)
  for (r in seq_len(R)) {
    d <- simulate_longitudinal_mar_outcome(n = 2000, seed = 20000 + r)
    dt <- data.table::as.data.table(d)
    fit <- tryCatch(
      causat(
        dt,
        outcome = "Y",
        treatment = "A",
        confounders = ~L0,
        confounders_tv = ~L,
        estimator = "aipw",
        type = "longitudinal",
        id = "id",
        time = "time",
        censoring = "C",
        ipcw = TRUE
      ),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      next
    }
    sw <- tryCatch(
      contrast(
        fit,
        interventions = ints,
        type = "difference",
        ci_method = "sandwich"
      ),
      error = function(e) NULL
    )
    if (is.null(sw)) {
      next
    }
    est[r, ] <- sw$estimates$estimate
    se[r, ] <- sw$estimates$se
  }
  ok <- stats::complete.cases(est)
  ratio <- colMeans(se[ok, ]) / apply(est[ok, ], 2L, stats::sd)
  expect_equal(mean(ratio), 1, tolerance = 0.05)
  expect_true(all(ratio > 0.9 & ratio < 1.1))
})

test_that("AIPW IPCW sandwich SE agrees with the bootstrap", {
  skip_if_fast()
  d <- simulate_longitudinal_mar_outcome(n = 3000, seed = 41)
  dt <- data.table::as.data.table(d)
  fit <- causat(
    dt,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    estimator = "aipw",
    type = "longitudinal",
    id = "id",
    time = "time",
    censoring = "C",
    ipcw = TRUE
  )
  ints <- list(a1 = static(1), a0 = static(0))
  rs <- contrast(
    fit,
    interventions = ints,
    type = "difference",
    ci_method = "sandwich"
  )
  rb <- contrast(
    fit,
    interventions = ints,
    type = "difference",
    ci_method = "bootstrap",
    n_boot = 400
  )
  expect_equal(rs$estimates$estimate, rb$estimates$estimate, tolerance = 1e-8)
  expect_equal(rs$contrasts$se / rb$contrasts$se, 1, tolerance = 0.12)
})
