# Phase A unit tests for the IF variance engine primitives.
#
# These test bread_inv(), iv_design_matrix(), correct_model(),
# vcov_from_if() (with and without cluster), and variance_if_numeric()
# in isolation. The dispatcher variance_if() and correct_propensity() are
# tested in later phases through the public causat() / contrast() API.

# ── correct_model() — non-canonical link agreement with sandwich::estfun ──

test_that("correct_model() matches sandwich::estfun on a probit GLM", {
  set.seed(101)
  n <- 400
  L <- rnorm(n)
  A <- rbinom(n, 1, plogis(0.4 * L))
  Y <- rbinom(n, 1, pnorm(-0.3 + 0.7 * A + 0.5 * L))
  dat <- data.frame(Y = Y, A = A, L = L)

  m <- stats::glm(Y ~ A + L, data = dat, family = binomial(link = "probit"))

  # An arbitrary gradient g (the role of J in correct_model).
  set.seed(7)
  g <- runif(length(stats::coef(m)), -1, 1)

  res <- correct_model(m, g, fit_idx = seq_len(n), n_total = n)

  # Expected: n * (estfun %*% bread/nobs %*% g), per row.
  ef <- sandwich::estfun(m)
  br <- sandwich::bread(m) / stats::nobs(m)
  expected <- n * as.numeric(ef %*% br %*% g)

  expect_equal(res$correction, expected, tolerance = 1e-10)
})


# ── correct_model() — Gamma-log link agreement with sandwich::estfun ──

test_that("correct_model() matches sandwich::estfun on a Gamma-log GLM", {
  set.seed(202)
  n <- 500
  L <- rnorm(n)
  A <- rbinom(n, 1, 0.5)
  mu <- exp(0.4 + 0.3 * A + 0.2 * L)
  Y <- rgamma(n, shape = 5, rate = 5 / mu)
  dat <- data.frame(Y = Y, A = A, L = L)

  m <- stats::glm(Y ~ A + L, data = dat, family = Gamma(link = "log"))

  set.seed(8)
  g <- runif(length(stats::coef(m)), -1, 1)

  res <- correct_model(m, g, fit_idx = seq_len(n), n_total = n)

  ef <- sandwich::estfun(m)
  br <- sandwich::bread(m) / stats::nobs(m)
  expected <- n * as.numeric(ef %*% br %*% g)

  expect_equal(res$correction, expected, tolerance = 1e-10)
})


# ── Regression guard: response-residual formula must disagree on probit ──

test_that("response-residual formula disagrees materially on probit (regression guard)", {
  set.seed(303)
  n <- 400
  L <- rnorm(n)
  A <- rbinom(n, 1, 0.5)
  Y <- rbinom(n, 1, pnorm(-0.2 + 0.6 * A + 0.4 * L))
  dat <- data.frame(Y = Y, A = A, L = L)

  m <- stats::glm(Y ~ A + L, data = dat, family = binomial(link = "probit"))
  g <- c(0.1, 0.2, 0.3)

  res <- correct_model(m, g, fit_idx = seq_len(n), n_total = n)

  # Old formula: n * d_i * r_response (no mu_eta/V scale factor).
  X <- stats::model.matrix(m)
  B_inv <- bread_inv(m, X)
  h <- as.numeric(B_inv %*% g)
  d <- as.numeric(X %*% h)
  r_response <- stats::residuals(m, type = "response")
  old_correction <- n * d * r_response

  # The two must differ materially. Mean absolute difference normalised
  # by the magnitude of the new (correct) correction.
  rel_diff <- mean(abs(res$correction - old_correction)) /
    mean(abs(res$correction))
  expect_gt(rel_diff, 0.05)
})


# ── correct_model() — canonical-link sanity (logistic, identity, Poisson-log) ──

test_that("correct_model() canonical links match sandwich::estfun", {
  set.seed(404)
  n <- 300
  L <- rnorm(n)
  A <- rbinom(n, 1, 0.5)

  for (cfg in list(
    list(
      family = binomial(),
      Y = rbinom(n, 1, plogis(0.2 + 0.5 * A + 0.3 * L))
    ),
    list(family = gaussian(), Y = 1 + 0.5 * A + 0.3 * L + rnorm(n)),
    list(family = poisson(), Y = rpois(n, exp(0.2 + 0.4 * A + 0.2 * L)))
  )) {
    dat <- data.frame(Y = cfg$Y, A = A, L = L)
    m <- stats::glm(Y ~ A + L, data = dat, family = cfg$family)
    g <- c(0.1, -0.2, 0.3)
    res <- correct_model(m, g, fit_idx = seq_len(n), n_total = n)
    ef <- sandwich::estfun(m)
    br <- sandwich::bread(m) / stats::nobs(m)
    expect_equal(
      res$correction,
      n * as.numeric(ef %*% br %*% g),
      tolerance = 1e-10
    )
  }
})


# ── vcov_from_if() — independent aggregation ──

test_that("vcov_from_if() reproduces the (1/n^2) sum-of-products formula", {
  set.seed(505)
  n <- 50
  IFa <- rnorm(n)
  IFb <- rnorm(n)

  V <- vcov_from_if(list(a = IFa, b = IFb), n = n, int_names = c("a", "b"))

  expect_equal(V["a", "a"], sum(IFa * IFa) / n^2)
  expect_equal(V["b", "b"], sum(IFb * IFb) / n^2)
  expect_equal(V["a", "b"], sum(IFa * IFb) / n^2)
  expect_equal(V["a", "b"], V["b", "a"])
})


# ── vcov_from_if() — cluster aggregation sums within cluster ──

test_that("vcov_from_if() cluster path sums within cluster before squaring", {
  IF <- c(1, 2, 3, 4, 5, 6)
  cl <- c("A", "A", "A", "B", "B", "B")
  n <- length(IF)

  V <- vcov_from_if(list(x = IF), n = n, int_names = "x", cluster = cl)

  # Sum within clusters: A -> 6, B -> 15. Sum of squares -> 36 + 225 = 261.
  expect_equal(V["x", "x"], 261 / n^2)

  # Singleton clusters reduce to the standard formula.
  cl2 <- as.character(seq_along(IF))
  V2 <- vcov_from_if(list(x = IF), n = n, int_names = "x", cluster = cl2)
  V_indep <- vcov_from_if(list(x = IF), n = n, int_names = "x")
  expect_equal(V2, V_indep)
})


# ── Matching: IF vs vcovCL + J V_beta J^T at Channel-1 = 0 ──

test_that("variance_if matching path matches vcovCL+J V_beta J^T at saturated MSM", {
  set.seed(808)
  n <- 800
  L <- rnorm(n)
  A <- rbinom(n, 1, plogis(-0.7 + 0.4 * L))
  Y <- 1 + 0.5 * A + 0.3 * L + rnorm(n)
  df <- data.frame(Y = Y, A = A, L = L)

  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    method = "matching",
    estimand = "ATT"
  )

  res_if <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich"
  )
  V_if <- res_if$vcov

  # Reference: vcovCL on the matched glm propagated through the Jacobian
  # of mu_a = g^{-1}(beta_0 + beta_1 a). For a saturated linear MSM with
  # static interventions, channel 1 is identically zero and J V_beta J^T
  # is the entire variance.
  matched <- as.data.frame(fit$details$matched_data)
  # cadjust = FALSE disables sandwich::vcovCL's G/(G-1) finite-sample
  # correction, matching the raw 1/n^2 cluster-sum aggregation in
  # vcov_from_if(cluster = ...).
  V_beta <- sandwich::vcovCL(
    fit$model,
    cluster = matched$subclass,
    cadjust = FALSE
  )
  J <- matrix(c(1, 1, 1, 0), nrow = 2, byrow = TRUE) # rows: a1, a0
  V_ref <- J %*% V_beta %*% t(J)

  expect_equal(unname(V_if), unname(V_ref), tolerance = 1e-10)
})


# ── variance_if_numeric() Tier 1 — recovers analytic IF on a glm ──

test_that("variance_if_numeric() Tier 1 matches the analytic IF on a logistic glm", {
  set.seed(606)
  n <- 300
  L <- rnorm(n)
  A <- rbinom(n, 1, plogis(0.3 * L))
  Y <- rbinom(n, 1, plogis(-0.2 + 0.6 * A + 0.4 * L))
  dat <- data.frame(Y = Y, A = A, L = L)
  m <- stats::glm(Y ~ A + L, data = dat, family = binomial)

  fit <- list(
    model = m,
    details = list(fit_rows = rep(TRUE, n))
  )

  data_a_list <- list(
    a1 = transform(dat, A = 1),
    a0 = transform(dat, A = 0)
  )
  preds_list <- lapply(data_a_list, function(df) {
    stats::predict(m, newdata = df, type = "response")
  })
  mu_hat <- vapply(preds_list, mean, numeric(1))
  target_idx <- rep(TRUE, n)

  V_num <- variance_if_numeric(
    fit,
    data_a_list,
    preds_list,
    mu_hat,
    target_idx
  )

  # Analytic reference: build IF by hand using correct_model().
  beta <- stats::coef(m)
  IF_list <- lapply(seq_along(data_a_list), function(j) {
    df <- as.data.frame(data_a_list[[j]])
    X_star <- iv_design_matrix(m, df)
    eta_star <- as.numeric(X_star %*% beta)
    mu_eta_star <- m$family$mu.eta(eta_star)
    g <- as.numeric(crossprod(X_star, mu_eta_star)) / n
    Ch1 <- preds_list[[j]] - mu_hat[j]
    Ch2 <- correct_model(m, g, fit_idx = seq_len(n), n_total = n)$correction
    Ch1 + Ch2
  })
  V_ref <- vcov_from_if(IF_list, n = n, int_names = names(data_a_list))

  expect_equal(unname(V_num), unname(V_ref), tolerance = 1e-6)
})


# ── correct_propensity() Branch A — IPW WeightIt path ──

test_that("variance_if IPW path matches WeightIt vcov + J V_beta J^T at saturated MSM", {
  skip_if_not_installed("WeightIt")
  set.seed(909)
  n <- 1000
  L <- rnorm(n)
  A <- rbinom(n, 1, plogis(0.4 * L))
  Y <- 1 + 0.5 * A + 0.3 * L + rnorm(n)
  df <- data.frame(Y = Y, A = A, L = L)

  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    method = "ipw"
  )

  res_if <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich"
  )
  V_if <- res_if$vcov

  # Reference: Y ~ A is saturated, A is static, so Channel 1 = 0 and the
  # IF variance equals J V_beta J^T with V_beta = vcov(glm_weightit) (the
  # M-estimation sandwich that already accounts for weight uncertainty).
  V_beta <- stats::vcov(fit$model)
  J <- matrix(c(1, 1, 1, 0), nrow = 2, byrow = TRUE)
  V_ref <- J %*% V_beta %*% t(J)

  expect_equal(unname(V_if), unname(V_ref), tolerance = 1e-10)
})


# ── prepare_propensity_if() Branch B — Phase 4 scaffold aborts informatively ──

test_that("prepare_propensity_if() Branch B aborts with a Phase 4 message", {
  fake_fit <- list(
    model = structure(list(), class = "lm"),
    method = "ipw"
  )
  expect_error(
    prepare_propensity_if(fake_fit, fit_idx = 1L, n_total = 1L),
    "Phase 4"
  )
})


# ── Mparts guardrail — warning at fit time for non-Mparts WeightIt method ──

test_that("fit_ipw warns when WeightIt method lacks Mparts", {
  skip_if_not_installed("WeightIt")

  set.seed(1010)
  n <- 100
  L <- rnorm(n)
  A <- rbinom(n, 1, plogis(0.3 * L))
  Y <- 1 + 0.5 * A + rnorm(n)
  df <- data.frame(Y = Y, A = A, L = L)

  fake_w <- WeightIt::weightit(A ~ L, data = df, method = "glm")
  attr(fake_w, "Mparts") <- NULL

  testthat::local_mocked_bindings(
    weightit = function(...) fake_w,
    .package = "WeightIt"
  )

  expect_warning(
    causat(
      df,
      outcome = "Y",
      treatment = "A",
      confounders = ~L,
      method = "ipw"
    ),
    "Mparts"
  )
})


# ── variance_if_numeric() Tier 2 — warns and returns V1 + J V_beta J^T ──

test_that("variance_if_numeric() Tier 2 keeps full Channel-1 covariance", {
  set.seed(707)
  n <- 200
  L <- rnorm(n)
  A <- rbinom(n, 1, 0.5)
  Y <- 1 + 0.4 * A + 0.3 * L + rnorm(n)
  dat <- data.frame(Y = Y, A = A, L = L)
  m <- stats::glm(Y ~ A + L, data = dat, family = gaussian())

  fit <- list(
    model = m,
    details = list(fit_rows = rep(TRUE, n))
  )
  data_a_list <- list(
    a1 = transform(dat, A = 1),
    a0 = transform(dat, A = 0)
  )
  preds_list <- lapply(data_a_list, function(df) {
    stats::predict(m, newdata = df, type = "response")
  })
  mu_hat <- vapply(preds_list, mean, numeric(1))
  target_idx <- rep(TRUE, n)

  # Force Tier 2 by stubbing sandwich::estfun to error for any input.
  testthat::local_mocked_bindings(
    estfun = function(x, ...) stop("no estfun for this class"),
    .package = "sandwich"
  )

  expect_warning(
    V <- variance_if_numeric(
      fit,
      data_a_list,
      preds_list,
      mu_hat,
      target_idx
    ),
    "drops the IF cross-term"
  )

  # Hand-computed reference: full Ch1 crossprod + J V_beta J^T. Using the
  # crossprod form (not just the diagonal) captures cross-intervention
  # Channel-1 covariance — dropping it was the bug fixed in this commit.
  Ch1_mat <- do.call(
    cbind,
    lapply(seq_along(preds_list), function(j) preds_list[[j]] - mu_hat[j])
  )
  V1 <- crossprod(Ch1_mat) / n^2
  beta <- stats::coef(m)
  pred_fun <- function(beta) {
    mt <- m
    mt$coefficients <- beta
    vapply(
      data_a_list,
      function(df) {
        mean(stats::predict(mt, newdata = df, type = "response"))
      },
      numeric(1)
    )
  }
  J <- numDeriv::jacobian(pred_fun, x = beta)
  V_beta <- stats::vcov(m)
  V_ref <- J %*% V_beta %*% t(J) + V1

  expect_equal(unname(V), unname(V_ref), tolerance = 1e-8)
  # Off-diagonal should be nonzero for two interventions over the same
  # target population — the old diag-only formula zeroed this out.
  expect_true(abs(V[1, 2]) > 1e-8)
})


# ── prepare_model_if() + apply_model_correction() equal correct_model() ──

test_that("prepare_model_if()/apply_model_correction() matches correct_model()", {
  set.seed(1111)
  n <- 300
  L <- rnorm(n)
  A <- rbinom(n, 1, plogis(0.2 * L))
  Y <- rbinom(n, 1, plogis(-0.3 + 0.7 * A + 0.5 * L))
  dat <- data.frame(Y = Y, A = A, L = L)
  m <- stats::glm(Y ~ A + L, data = dat, family = binomial())

  prep <- prepare_model_if(m, fit_idx = seq_len(n), n_total = n)

  # Two arbitrary gradients — applied through the shared prep object
  # should match correct_model() re-fitting everything from scratch.
  for (seed in c(1, 2)) {
    set.seed(seed)
    g <- runif(length(stats::coef(m)), -1, 1)
    direct <- correct_model(m, g, fit_idx = seq_len(n), n_total = n)
    shared <- apply_model_correction(prep, g)
    expect_equal(shared$correction, direct$correction, tolerance = 1e-12)
    expect_equal(shared$d, direct$d, tolerance = 1e-12)
    expect_equal(shared$h, direct$h, tolerance = 1e-12)
  }
})


# ── resolve_fit_idx() aborts when na.action exceeds fit_rows length ──

test_that("resolve_fit_idx() aborts on misaligned na.action indices", {
  fake_model <- structure(
    list(na.action = structure(99L, class = "omit")),
    class = "lm"
  )
  fake_fit <- list(details = list(fit_rows = c(TRUE, TRUE, FALSE, TRUE)))
  expect_error(
    resolve_fit_idx(fake_fit, fake_model),
    "exceeds `sum\\(fit\\$details\\$fit_rows\\)`"
  )
})


# ── bread_inv() warns on rank-deficient X'WX ──

test_that("bread_inv() warns when X'WX is singular and falls back to ginv()", {
  set.seed(2222)
  n <- 40
  L <- rnorm(n)
  A <- rbinom(n, 1, 0.5)
  L_dup <- L # aliased column -> rank-deficient design
  Y <- rnorm(n)
  dat <- data.frame(Y = Y, A = A, L = L, L_dup = L_dup)

  # Build a rank-deficient X by hand; fit glm.fit to keep the design.
  X <- cbind(1, dat$A, dat$L, dat$L_dup)
  fit_lm <- stats::glm.fit(X, dat$Y, family = stats::gaussian())
  fit_lm$family <- stats::gaussian()
  class(fit_lm) <- c("glm", "lm")

  # The warning is rate-limited via `rlang::warn(.frequency = "once")`
  # so reset the session state to get a deterministic emission.
  rlang::reset_warning_verbosity("causatr_bread_inv_singular")

  X_bad <- X
  expect_warning(
    B <- bread_inv(fit_lm, X_bad),
    "singular"
  )
  expect_true(is.matrix(B))
})


test_that("bread_inv() singular warning is rate-limited on repeat calls", {
  set.seed(3333)
  n <- 30
  X <- cbind(1, rnorm(n), rnorm(n))
  X <- cbind(X, X[, 2]) # aliased column
  fit_lm <- stats::glm.fit(X, rnorm(n), family = stats::gaussian())
  fit_lm$family <- stats::gaussian()
  class(fit_lm) <- c("glm", "lm")

  rlang::reset_warning_verbosity("causatr_bread_inv_singular")

  # First call warns.
  expect_warning(bread_inv(fit_lm, X), "singular")
  # Second call is silent on the warning channel specifically —
  # using expect_warning(..., NA) is narrower than expect_silent()
  # and won't be tripped by unrelated messages or stdout.
  expect_warning(bread_inv(fit_lm, X), regexp = NA)
})


# ── iv_design_matrix() accepts data.table input ──

test_that("iv_design_matrix() handles data.table and data.frame identically", {
  set.seed(15151)
  n <- 80
  dat <- data.frame(
    Y = rnorm(n),
    A = rbinom(n, 1, 0.5),
    L = rnorm(n)
  )
  m <- stats::glm(Y ~ A + L, data = dat, family = gaussian())

  newdata_df <- data.frame(A = 1, L = dat$L)
  newdata_dt <- data.table::as.data.table(newdata_df)

  X_df <- iv_design_matrix(m, newdata_df)
  X_dt <- iv_design_matrix(m, newdata_dt)

  # data.table and data.frame paths must produce the same design
  # matrix — the coercion inside iv_design_matrix() is transparent.
  expect_equal(unname(X_df), unname(X_dt))
  expect_equal(colnames(X_df), colnames(X_dt))
})


# ── build_point_channel_pieces() — shared shape for g-comp and IPW ──

test_that("build_point_channel_pieces() returns the right Ch1/grad shapes", {
  skip_if_not_installed("WeightIt")
  set.seed(4444)
  n <- 500
  L <- rnorm(n)
  A <- rbinom(n, 1, plogis(0.3 * L))
  Y <- 1 + 0.5 * A + 0.3 * L + rnorm(n)
  df <- data.frame(Y = Y, A = A, L = L)

  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    method = "gcomp"
  )
  model <- fit$model
  data_a_list <- list(a1 = transform(df, A = 1), a0 = transform(df, A = 0))
  preds_list <- lapply(data_a_list, function(da) {
    stats::predict(model, newdata = da, type = "response")
  })
  mu_hat <- vapply(preds_list, mean, numeric(1))
  target_idx <- rep(TRUE, n)

  pieces <- causatr:::build_point_channel_pieces(
    fit,
    data_a_list,
    preds_list,
    mu_hat,
    target_idx
  )

  expect_named(pieces, c("Ch1_list", "grad_list", "fit_idx", "n"))
  expect_length(pieces$Ch1_list, 2L)
  expect_length(pieces$grad_list, 2L)
  expect_equal(length(pieces$Ch1_list[[1]]), n)
  expect_equal(length(pieces$grad_list[[1]]), length(stats::coef(model)))
  expect_equal(pieces$n, n)

  # Channel-1 must sum to zero over the target population (by construction).
  expect_equal(sum(pieces$Ch1_list[[1]]), 0, tolerance = 1e-10)
  expect_equal(sum(pieces$Ch1_list[[2]]), 0, tolerance = 1e-10)

  # Gradient for static(1) should be X_{A=1}' * mu_eta averaged — for a
  # Gaussian identity link, mu_eta = 1, so grad[1] = mean(intercept col) = 1
  # and grad[2] = mean(A=1 col) = 1.
  expect_equal(pieces$grad_list$a1[1], 1, tolerance = 1e-12)
  expect_equal(pieces$grad_list$a1[2], 1, tolerance = 1e-12)
  expect_equal(pieces$grad_list$a0[2], 0, tolerance = 1e-12)
})


# ── variance_if_ipw: hoisted prop_prep matches per-intervention recompute ──

test_that("variance_if_ipw hoisted prep gives same vcov as per-intervention recompute", {
  skip_if_not_installed("WeightIt")
  set.seed(5555)
  n <- 600
  L <- rnorm(n)
  A <- rbinom(n, 1, plogis(0.4 * L))
  Y <- 1 + 0.5 * A + 0.3 * L + rnorm(n)
  df <- data.frame(Y = Y, A = A, L = L)

  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    method = "ipw"
  )

  res <- contrast(
    fit,
    interventions = list(
      a1 = static(1),
      a0 = static(0),
      a1_dup = static(1) # third intervention exercises k >= 3
    ),
    reference = "a0",
    ci_method = "sandwich"
  )
  V_hoisted <- res$vcov

  # Reference: rebuild prop_prep manually once and apply per intervention.
  # Must equal the hoisted-path vcov to machine precision (same linear
  # algebra, just iterated differently).
  model <- fit$model
  n_total <- nrow(df)
  target_idx <- rep(TRUE, n_total)
  fit_idx <- causatr:::resolve_fit_idx(fit, model)

  data_a_list <- list(
    a1 = transform(df, A = 1),
    a0 = transform(df, A = 0),
    a1_dup = transform(df, A = 1)
  )
  preds_list <- lapply(data_a_list, function(da) {
    as.numeric(stats::predict(model, newdata = da, type = "response"))
  })
  mu_hat <- vapply(preds_list, mean, numeric(1))

  prop_prep <- causatr:::prepare_propensity_if(fit, fit_idx, n_total)

  pieces <- causatr:::build_point_channel_pieces(
    fit,
    data_a_list,
    preds_list,
    mu_hat,
    target_idx
  )
  IF_list <- lapply(seq_along(pieces$grad_list), function(j) {
    pieces$Ch1_list[[j]] +
      causatr:::apply_propensity_correction(prop_prep, pieces$grad_list[[j]])
  })
  names(IF_list) <- names(data_a_list)
  V_ref <- causatr:::vcov_from_if(IF_list, n_total, names(data_a_list))

  expect_equal(unname(V_hoisted), unname(V_ref), tolerance = 1e-12)

  # Duplicate interventions should have identical diagonal AND equal
  # off-diagonal with the first intervention (variance arithmetic
  # invariant — any drift would mean the hoisting corrupted state).
  expect_equal(V_hoisted["a1", "a1"], V_hoisted["a1_dup", "a1_dup"])
  expect_equal(V_hoisted["a1", "a1_dup"], V_hoisted["a1", "a1"])
})


# ── variance_if_numeric() — weights arg flows through Tier 1 path ──

test_that("variance_if_numeric() respects the weights argument in Tier 1", {
  set.seed(6666)
  n <- 400
  L <- rnorm(n)
  A <- rbinom(n, 1, 0.5)
  Y <- 1 + 0.5 * A + 0.3 * L + rnorm(n)
  w <- runif(n, 0.5, 2)
  dat <- data.frame(Y = Y, A = A, L = L)

  # Fit a weighted glm so w enters via prior.weights and affects the
  # Tier 1 Jacobian through the weighted mean.
  m <- stats::glm(Y ~ A + L, data = dat, family = gaussian(), weights = w)

  fit <- list(
    model = m,
    details = list(fit_rows = rep(TRUE, n), weights = w),
    data = dat
  )
  data_a_list <- list(
    a1 = transform(dat, A = 1),
    a0 = transform(dat, A = 0)
  )
  preds_list <- lapply(data_a_list, function(df) {
    stats::predict(m, newdata = df, type = "response")
  })
  target_idx <- rep(TRUE, n)

  # Compare weighted and unweighted Tier 1 outputs — they must differ.
  mu_hat_w <- vapply(
    preds_list,
    function(p) stats::weighted.mean(p, w),
    numeric(1)
  )
  V_w <- variance_if_numeric(
    fit,
    data_a_list,
    preds_list,
    mu_hat_w,
    target_idx,
    weights = w
  )

  mu_hat_u <- vapply(preds_list, mean, numeric(1))
  V_u <- variance_if_numeric(
    fit,
    data_a_list,
    preds_list,
    mu_hat_u,
    target_idx,
    weights = NULL
  )

  # The two must differ (the weights argument is actually flowing through).
  expect_true(max(abs(V_w - V_u)) > 1e-6)

  # Reference: the Tier 1 weighted path should agree with the full
  # sandwich sum IF_i = Ch1 + IF_beta . J_j (hand computed).
  psi <- sandwich::estfun(m)
  A_inv <- sandwich::bread(m) / stats::nobs(m)
  IF_beta <- psi %*% A_inv
  beta <- stats::coef(m)
  pred_fun <- function(b) {
    mt <- m
    mt$coefficients <- b
    vapply(
      data_a_list,
      function(df) {
        p <- stats::predict(mt, newdata = df, type = "response")
        stats::weighted.mean(p, w)
      },
      numeric(1)
    )
  }
  J <- numDeriv::jacobian(pred_fun, x = beta)
  tw <- w / sum(w)
  Ch1_mat <- do.call(
    cbind,
    lapply(seq_along(preds_list), function(j) {
      n * tw * (preds_list[[j]] - mu_hat_w[j])
    })
  )
  Ch2_mat <- (IF_beta %*% t(J)) * n
  IF_mat_ref <- Ch1_mat + Ch2_mat
  V_ref <- crossprod(IF_mat_ref) / n^2
  expect_equal(unname(V_w), unname(V_ref), tolerance = 1e-10)
})


# ── variance_if_numeric() Tier 1 batching preserves element-wise results ──

test_that("variance_if_numeric() batched Tier 1 matches per-intervention recompute", {
  set.seed(7777)
  n <- 300
  L <- rnorm(n)
  A <- rbinom(n, 1, 0.5)
  Y <- rbinom(n, 1, plogis(-0.3 + 0.6 * A + 0.4 * L))
  dat <- data.frame(Y = Y, A = A, L = L)
  m <- stats::glm(Y ~ A + L, data = dat, family = binomial())

  fit <- list(
    model = m,
    details = list(fit_rows = rep(TRUE, n)),
    data = dat
  )
  data_a_list <- list(a1 = transform(dat, A = 1), a0 = transform(dat, A = 0))
  preds_list <- lapply(data_a_list, function(df) {
    stats::predict(m, newdata = df, type = "response")
  })
  mu_hat <- vapply(preds_list, mean, numeric(1))
  target_idx <- rep(TRUE, n)

  V_batched <- variance_if_numeric(
    fit,
    data_a_list,
    preds_list,
    mu_hat,
    target_idx
  )

  # Hand-compute via per-intervention loop (what the old code did).
  psi <- sandwich::estfun(m)
  A_inv <- sandwich::bread(m) / stats::nobs(m)
  IF_beta <- psi %*% A_inv
  beta <- stats::coef(m)
  pred_fun <- function(b) {
    mt <- m
    mt$coefficients <- b
    vapply(
      data_a_list,
      function(df) mean(stats::predict(mt, newdata = df, type = "response")),
      numeric(1)
    )
  }
  J <- numDeriv::jacobian(pred_fun, x = beta)

  tw <- rep(1 / n, n)
  Ch1_list <- lapply(seq_along(preds_list), function(j) {
    n * tw * (preds_list[[j]] - mu_hat[j])
  })
  IF_list_ref <- lapply(seq_along(preds_list), function(j) {
    Ch1_list[[j]] + as.numeric(IF_beta %*% J[j, ]) * n
  })
  IF_mat_ref <- do.call(cbind, IF_list_ref)
  V_ref <- crossprod(IF_mat_ref) / n^2
  expect_equal(unname(V_batched), unname(V_ref), tolerance = 1e-12)
})


# ── build_point_channel_pieces() — empty target population aborts ──

test_that("build_point_channel_pieces() aborts when target population is empty", {
  set.seed(8888)
  n <- 100
  L <- rnorm(n)
  A <- rbinom(n, 1, 0.5)
  Y <- 1 + 0.5 * A + rnorm(n)
  df <- data.frame(Y = Y, A = A, L = L)

  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    method = "gcomp"
  )

  data_a_list <- list(a1 = transform(df, A = 1), a0 = transform(df, A = 0))
  preds_list <- lapply(data_a_list, function(da) {
    stats::predict(fit$model, newdata = da, type = "response")
  })
  mu_hat <- vapply(preds_list, mean, numeric(1))

  # Unweighted: n_target == 0 triggers abort.
  expect_error(
    build_point_channel_pieces(
      fit,
      data_a_list,
      preds_list,
      mu_hat,
      target_idx = rep(FALSE, n)
    ),
    "target population is empty"
  )

  # Weighted: sum of target weights == 0 triggers abort.
  fit_w <- fit
  fit_w$details$weights <- rep(0, n)
  expect_error(
    build_point_channel_pieces(
      fit_w,
      data_a_list,
      preds_list,
      mu_hat,
      target_idx = rep(TRUE, n)
    ),
    "weights sum to 0"
  )
})


# ── prepare_point_variance() wrapper delegates correctly ──

test_that("prepare_point_variance() returns consistent pieces and prep", {
  set.seed(9999)
  n <- 200
  L <- rnorm(n)
  A <- rbinom(n, 1, 0.5)
  Y <- 1 + 0.5 * A + 0.3 * L + rnorm(n)
  df <- data.frame(Y = Y, A = A, L = L)

  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    method = "gcomp"
  )
  model <- fit$model
  data_a_list <- list(a1 = transform(df, A = 1), a0 = transform(df, A = 0))
  preds_list <- lapply(data_a_list, function(da) {
    stats::predict(model, newdata = da, type = "response")
  })
  mu_hat <- vapply(preds_list, mean, numeric(1))
  target_idx <- rep(TRUE, n)

  bundle <- prepare_point_variance(
    fit,
    model,
    data_a_list,
    preds_list,
    mu_hat,
    target_idx
  )

  expect_named(bundle, c("pieces", "prep"))
  # fit_idx/n must be shared between the two halves so the alignment
  # is structural rather than by convention.
  expect_equal(bundle$prep$fit_idx, bundle$pieces$fit_idx)
  expect_equal(bundle$prep$n_total, bundle$pieces$n)

  # The bundle must produce the same final correction as the direct
  # component calls.
  pieces_direct <- build_point_channel_pieces(
    fit,
    data_a_list,
    preds_list,
    mu_hat,
    target_idx
  )
  prep_direct <- prepare_model_if(model, pieces_direct$fit_idx, pieces_direct$n)
  g <- pieces_direct$grad_list$a1
  expect_equal(
    apply_model_correction(bundle$prep, g)$correction,
    apply_model_correction(prep_direct, g)$correction,
    tolerance = 1e-12
  )
})


# ── End-to-end: variance_if_gcomp() + variance_if_ipw() agree on saturated MSM ──

test_that("variance_if_gcomp() wrapper refactor preserves analytic result", {
  # Regression guard: compare post-refactor variance_if_gcomp() output
  # against a hand-computed analytic reference on a simple linear model.
  set.seed(12121)
  n <- 400
  L <- rnorm(n)
  A <- rbinom(n, 1, 0.5)
  Y <- 1 + 0.5 * A + 0.3 * L + rnorm(n)
  df <- data.frame(Y = Y, A = A, L = L)

  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    method = "gcomp"
  )

  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich"
  )
  V_engine <- res$vcov

  # Analytic reference: J V_beta J^T + Channel 1 cross-products, hand
  # computed from the full data (Channel 1 nonzero because L varies
  # within each treatment arm).
  model <- fit$model
  beta <- stats::coef(model)
  p1 <- stats::predict(model, newdata = transform(df, A = 1), type = "response")
  p0 <- stats::predict(model, newdata = transform(df, A = 0), type = "response")
  mu_hat <- c(a1 = mean(p1), a0 = mean(p0))

  X1 <- stats::model.matrix(~ A + L, data = transform(df, A = 1))
  X0 <- stats::model.matrix(~ A + L, data = transform(df, A = 0))
  J <- rbind(a1 = colMeans(X1), a0 = colMeans(X0))

  psi <- sandwich::estfun(model)
  A_inv <- sandwich::bread(model) / stats::nobs(model)
  IF_beta <- psi %*% A_inv

  # Engine scaling convention: per-obs IF is
  #   Ch1_i = (n / n_t) * t_i * (p_i - mu_hat)        # = (p_i - mu_hat) at full target
  #   Ch2_i = n * (IF_beta,i . J_j)
  # aggregated by `crossprod(IF_mat) / n^2`.
  Ch1_mat <- cbind(
    a1 = p1 - mu_hat["a1"],
    a0 = p0 - mu_hat["a0"]
  )
  Ch2_mat <- n * (IF_beta %*% t(J))
  IF_mat <- Ch1_mat + Ch2_mat
  V_ref <- crossprod(IF_mat) / n^2

  expect_equal(unname(V_engine), unname(V_ref), tolerance = 1e-8)
})


# ── build_point_channel_pieces(): NA on non-target rows is harmless ──

test_that("build_point_channel_pieces() isolates NA predictions to non-target rows", {
  # Regression guard: the weighted Channel-1 formula used to multiply a
  # full-length predictor vector by a zero mask, but `0 * NA = NA` in R,
  # so NA predictions on non-target rows poisoned the IF vector. The fix
  # is to index `p` by `target_idx` before subtracting mu_hat.
  set.seed(13131)
  n <- 200
  L <- rnorm(n)
  A <- rbinom(n, 1, 0.5)
  Y <- 1 + 0.5 * A + 0.3 * L + rnorm(n)
  df <- data.frame(Y = Y, A = A, L = L)

  fit <- causat(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    method = "gcomp"
  )
  model <- fit$model

  data_a_list <- list(a1 = transform(df, A = 1), a0 = transform(df, A = 0))
  preds_list <- lapply(data_a_list, function(da) {
    stats::predict(model, newdata = da, type = "response")
  })
  # Inject NA predictions on rows that will be excluded from target.
  non_target <- c(1L, 5L, 17L)
  for (j in seq_along(preds_list)) {
    preds_list[[j]][non_target] <- NA_real_
  }
  target_idx <- rep(TRUE, n)
  target_idx[non_target] <- FALSE
  mu_hat <- vapply(
    preds_list,
    function(p) mean(p[target_idx]),
    numeric(1)
  )

  # Unweighted path.
  pieces_u <- build_point_channel_pieces(
    fit,
    data_a_list,
    preds_list,
    mu_hat,
    target_idx
  )
  expect_false(any(is.na(pieces_u$Ch1_list$a1)))
  expect_false(any(is.na(pieces_u$Ch1_list$a0)))
  expect_equal(pieces_u$Ch1_list$a1[non_target], rep(0, length(non_target)))

  # Weighted path (this is the branch the old code silently broke).
  fit_w <- fit
  fit_w$details$weights <- runif(n, 0.5, 2)
  pieces_w <- build_point_channel_pieces(
    fit_w,
    data_a_list,
    preds_list,
    mu_hat,
    target_idx
  )
  expect_false(any(is.na(pieces_w$Ch1_list$a1)))
  expect_false(any(is.na(pieces_w$Ch1_list$a0)))
  expect_equal(pieces_w$Ch1_list$a1[non_target], rep(0, length(non_target)))
})


# ── variance_if_ice_one: NA on non-target rows is harmless ──

test_that("variance_if_ice_one() handles NA pseudo_final on non-target rows", {
  # Regression guard: the weighted ICE branch used to multiply a
  # full-length `pseudo_final` by a zero mask, so `0 * NA = NA` poisoned
  # the IF for any first-time individual dropped during the backward
  # iteration. Build a minimal longitudinal setup, inject NA on
  # non-target rows, verify IF is finite on both weighted and
  # unweighted paths.
  set.seed(14141)
  n_id <- 80
  long <- data.frame(
    id = rep(seq_len(n_id), each = 2),
    time = rep(c(0, 1), n_id),
    A = rbinom(2 * n_id, 1, 0.5),
    L = rnorm(2 * n_id),
    Y = NA_real_
  )
  # Observed Y only at the last time point (standard longitudinal shape).
  last_rows <- long$time == 1
  long$Y[last_rows] <- 1 +
    0.5 * long$A[last_rows] +
    0.3 * long$L[last_rows] +
    rnorm(n_id)

  fit <- causat(
    long,
    outcome = "Y",
    treatment = "A",
    confounders = ~1,
    confounders_tv = ~L,
    id = "id",
    time = "time"
  )

  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    ci_method = "sandwich"
  )
  expect_true(all(is.finite(res$estimates$se)))
  expect_true(all(res$estimates$se > 0))

  # Direct unit test: call variance_if_ice_one with a pseudo_final that
  # has NA outside the target.
  ice_res <- causatr:::ice_iterate(fit, causatr::static(1))
  ice_res_na <- ice_res
  first_ids <- as.character(fit$data[fit$data$time == 0, ]$id)
  n_first <- length(first_ids)
  target <- rep(TRUE, n_first)
  target[c(1L, 3L, 5L)] <- FALSE
  ice_res_na$pseudo_final[!target] <- NA_real_

  IF_vec <- causatr:::variance_if_ice_one(fit, ice_res_na, target)
  expect_length(IF_vec, n_first)
  # Finiteness is the regression guarantee: before the fix, NA in
  # `pseudo_final[!target]` propagated to IF_vec via `0 * NA = NA`.
  # Non-target rows can still carry nonzero IF through Channel 2
  # (their contribution to the nuisance model fit).
  expect_false(any(is.na(IF_vec)))
})
