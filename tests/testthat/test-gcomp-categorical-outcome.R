# Multinomial-outcome g-computation (Phase 23a-1).
#
# The estimand is the K-vector P(Y = k | do(A = a)) per intervention. Tests
# validate every step of the pipeline (fit -> per-class marginal means ->
# per-class contrasts -> bootstrap variance -> result schema -> S3 layer)
# against two oracles: the large-n softmax g-computation truth and
# `marginaleffects::avg_predictions()` run on causatr's own fitted model
# (an exact point oracle). Complex study designs (ATT, by-strata, weights,
# subset, continuous/categorical treatment, K = 4) are exercised explicitly.

# --- Point estimates: truth + marginaleffects parity -----------------------

test_that("multinomial point estimates match softmax truth and marginaleffects", {
  skip_if_not_installed("marginaleffects")
  d <- sim_multinom_binary(n = 8000, seed = 11)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "gcomp",
    model_fn = nnet::multinom,
    trace = FALSE
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "bootstrap",
    n_boot = 40
  )

  # Schema: class labels recorded, one vcov block per class, K rows per int.
  expect_equal(res$class_labels, mo_labels3)
  expect_setequal(names(res$vcov), mo_labels3)
  expect_equal(nrow(res$estimates), 2L * 3L)
  expect_true(all(c("intervention", "class") %in% names(res$estimates)))

  # Probabilities sum to 1 within each intervention.
  s <- tapply(res$estimates$estimate, res$estimates$intervention, sum)
  expect_equal(as.vector(s), c(1, 1), tolerance = 1e-8)

  mu1 <- res$estimates[intervention == "a1"]
  mu0 <- res$estimates[intervention == "a0"]
  est1 <- stats::setNames(mu1$estimate, mu1$class)[mo_labels3]
  est0 <- stats::setNames(mu0$estimate, mu0$class)[mo_labels3]

  # Exact parity with marginaleffects on causatr's own multinom fit.
  me1 <- me_multinom_means(fit$model, fit$data, "A", 1)
  me0 <- me_multinom_means(fit$model, fit$data, "A", 0)
  expect_equal(est1, me1[mo_labels3], tolerance = 1e-6)
  expect_equal(est0, me0[mo_labels3], tolerance = 1e-6)

  # Large-n softmax g-computation truth.
  t1 <- mo_gcomp_truth(
    mo_gen_al_binary,
    mo_eta3,
    mo_labels3,
    function(A, L) rep(1, length(L))
  )
  t0 <- mo_gcomp_truth(
    mo_gen_al_binary,
    mo_eta3,
    mo_labels3,
    function(A, L) rep(0, length(L))
  )
  expect_equal(est1, t1[mo_labels3], tolerance = 0.03)
  expect_equal(est0, t0[mo_labels3], tolerance = 0.03)
})

test_that("per-class difference contrasts match marginaleffects avg_comparisons", {
  skip_if_not_installed("marginaleffects")
  d <- sim_multinom_binary(n = 8000, seed = 11)
  fit <- causat(
    d,
    "Y",
    "A",
    confounders = ~L,
    estimator = "gcomp",
    model_fn = nnet::multinom,
    trace = FALSE
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    type = "difference",
    ci_method = "bootstrap",
    n_boot = 40
  )

  # marginaleffects RD for A: 0 -> 1, per class.
  me <- marginaleffects::avg_comparisons(
    fit$model,
    variables = list(A = c(0, 1)),
    type = "probs"
  )
  me_rd <- stats::setNames(me$estimate, as.character(me$group))[mo_labels3]
  con <- stats::setNames(res$contrasts$estimate, res$contrasts$class)[
    mo_labels3
  ]
  # causatr contrast is "a1 vs a0" = mu1 - mu0, matching A: 0 -> 1.
  expect_equal(con, me_rd, tolerance = 1e-6)
})

# --- Contrast scales: ratio and OR -----------------------------------------

test_that("per-class ratio and OR contrasts are exact functions of the means", {
  d <- sim_multinom_binary(n = 8000, seed = 11)
  fit <- causat(
    d,
    "Y",
    "A",
    confounders = ~L,
    estimator = "gcomp",
    model_fn = nnet::multinom,
    trace = FALSE
  )
  ivs <- list(a1 = static(1), a0 = static(0))

  res_d <- contrast(
    fit,
    ivs,
    reference = "a0",
    type = "difference",
    ci_method = "bootstrap",
    n_boot = 40
  )
  est <- res_d$estimates
  mu1 <- stats::setNames(
    est[intervention == "a1"]$estimate,
    est[intervention == "a1"]$class
  )[mo_labels3]
  mu0 <- stats::setNames(
    est[intervention == "a0"]$estimate,
    est[intervention == "a0"]$class
  )[mo_labels3]

  res_r <- contrast(
    fit,
    ivs,
    reference = "a0",
    type = "ratio",
    ci_method = "bootstrap",
    n_boot = 40
  )
  rr <- stats::setNames(res_r$contrasts$estimate, res_r$contrasts$class)[
    mo_labels3
  ]
  expect_equal(rr, mu1 / mu0, tolerance = 1e-8)

  res_or <- contrast(
    fit,
    ivs,
    reference = "a0",
    type = "or",
    ci_method = "bootstrap",
    n_boot = 40
  )
  or <- stats::setNames(res_or$contrasts$estimate, res_or$contrasts$class)[
    mo_labels3
  ]
  expect_equal(
    or,
    (mu1 / (1 - mu1)) / (mu0 / (1 - mu0)),
    tolerance = 1e-8
  )
})

# --- Treatment-type composition --------------------------------------------

test_that("multinomial composes with a continuous treatment + shift MTP", {
  skip_if_not_installed("marginaleffects")
  # A shift MTP extrapolates the fitted multinom to A + delta, so finite-sample
  # error against the softmax truth is larger than for a static intervention
  # (verified to converge: max gap 0.018 at n = 8000 -> 0.001 at n = 150000).
  # n = 20000 keeps the truth check comfortably inside tolerance; the exact
  # marginaleffects parity below is the primary point oracle.
  d <- sim_multinom_continuous(n = 20000, seed = 12)
  fit <- causat(
    d,
    "Y",
    "A",
    confounders = ~L,
    estimator = "gcomp",
    model_fn = nnet::multinom,
    trace = FALSE
  )
  res <- contrast(
    fit,
    interventions = list(up = shift(1), base = shift(0)),
    reference = "base",
    ci_method = "bootstrap",
    n_boot = 40
  )

  est_up <- res$estimates[intervention == "up"]
  up <- stats::setNames(est_up$estimate, est_up$class)[mo_labels3]

  # Oracle: marginaleffects on the shifted newdata (A := A + 1).
  nd <- as.data.frame(fit$data)
  nd$A <- nd$A + 1
  me <- marginaleffects::avg_predictions(
    fit$model,
    newdata = nd,
    type = "probs"
  )
  me_up <- stats::setNames(me$estimate, as.character(me$group))[mo_labels3]
  expect_equal(up, me_up, tolerance = 1e-6)

  # Large-n truth marginalises over the joint (A, L) under the shift.
  t_up <- mo_gcomp_truth(
    mo_gen_al_continuous,
    mo_eta3,
    mo_labels3,
    function(A, L) A + 1
  )
  expect_equal(up, t_up[mo_labels3], tolerance = 0.03)
})

test_that("multinomial composes with a categorical (3-level) treatment", {
  skip_if_not_installed("marginaleffects")
  d <- sim_multinom_cat_trt(n = 7000, seed = 13)
  fit <- causat(
    d,
    "Y",
    "A",
    confounders = ~L,
    estimator = "gcomp",
    model_fn = nnet::multinom,
    trace = FALSE
  )
  res <- contrast(
    fit,
    interventions = list(hi = static("hi"), lo = static("lo")),
    reference = "lo",
    ci_method = "bootstrap",
    n_boot = 40
  )
  est_hi <- res$estimates[intervention == "hi"]
  hi <- stats::setNames(est_hi$estimate, est_hi$class)[mo_labels3]
  me_hi <- me_multinom_means(fit$model, fit$data, "A", "hi")
  expect_equal(hi, me_hi[mo_labels3], tolerance = 1e-6)
  # Probabilities still sum to 1.
  s <- tapply(res$estimates$estimate, res$estimates$intervention, sum)
  expect_equal(as.vector(s), c(1, 1), tolerance = 1e-8)
})

test_that("the result schema generalises to K = 4 outcome classes", {
  skip_if_not_installed("marginaleffects")
  d <- sim_multinom_binary(
    n = 7000,
    seed = 14,
    eta_fn = mo_eta4,
    labels = mo_labels4
  )
  fit <- causat(
    d,
    "Y",
    "A",
    confounders = ~L,
    estimator = "gcomp",
    model_fn = nnet::multinom,
    trace = FALSE
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "bootstrap",
    n_boot = 40
  )
  expect_equal(res$class_labels, mo_labels4)
  expect_equal(nrow(res$estimates), 2L * 4L)
  expect_equal(length(res$vcov), 4L)
  s <- tapply(res$estimates$estimate, res$estimates$intervention, sum)
  expect_equal(as.vector(s), c(1, 1), tolerance = 1e-8)

  mu1 <- res$estimates[intervention == "a1"]
  est1 <- stats::setNames(mu1$estimate, mu1$class)[mo_labels4]
  me1 <- me_multinom_means(fit$model, fit$data, "A", 1)
  expect_equal(est1, me1[mo_labels4], tolerance = 1e-6)
})

# --- Complex study designs --------------------------------------------------

test_that("ATT standardises a multinomial outcome over the treated", {
  skip_if_not_installed("marginaleffects")
  d <- sim_multinom_binary(n = 9000, seed = 11)
  fit <- causat(
    d,
    "Y",
    "A",
    confounders = ~L,
    estimator = "gcomp",
    model_fn = nnet::multinom,
    trace = FALSE,
    estimand = "ATT"
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "bootstrap",
    n_boot = 40
  )
  expect_equal(res$estimand, "ATT")
  # Oracle: marginaleffects standardised over the treated rows only.
  treated <- fit$data[fit$data$A == 1, ]
  me1 <- me_multinom_means(fit$model, treated, "A", 1)
  mu1 <- res$estimates[intervention == "a1"]
  est1 <- stats::setNames(mu1$estimate, mu1$class)[mo_labels3]
  expect_equal(est1, me1[mo_labels3], tolerance = 1e-6)
})

test_that("by-stratified multinomial yields per-stratum per-class tables", {
  skip_if_not_installed("marginaleffects")
  d <- sim_multinom_binary(n = 9000, seed = 11)
  d$G <- factor(ifelse(d$L > 0, "hi", "lo"))
  fit <- causat(
    d,
    "Y",
    "A",
    confounders = ~ L + G,
    estimator = "gcomp",
    model_fn = nnet::multinom,
    trace = FALSE
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    by = "G",
    ci_method = "bootstrap",
    n_boot = 40
  )
  expect_true(all(c("by", "class") %in% names(res$estimates)))
  expect_equal(nrow(res$estimates), 2L * 2L * 3L)
  expect_equal(res$class_labels, mo_labels3)

  for (g in c("hi", "lo")) {
    sub <- fit$data[fit$data$G == g, ]
    me1 <- me_multinom_means(fit$model, sub, "A", 1)
    est_g <- res$estimates[by == g & intervention == "a1"]
    got <- stats::setNames(est_g$estimate, est_g$class)[mo_labels3]
    expect_equal(got, me1[mo_labels3], tolerance = 1e-6)
  }
})

test_that("by-stratified multinomial bootstrap result prints and confints", {
  # The by-stratified boot_info is a list-of-lists, and its `n_requested`
  # carries the user's `n_boot = 40` (a double); the print collapse must
  # coerce it before an integer-typed vapply. confint() over the same result
  # must emit one CI row per (by, intervention, class).
  d <- sim_multinom_binary(n = 6000, seed = 11)
  d$G <- factor(ifelse(d$L > 0, "hi", "lo"))
  fit <- causat(
    d,
    "Y",
    "A",
    confounders = ~ L + G,
    estimator = "gcomp",
    model_fn = nnet::multinom,
    trace = FALSE
  )
  res <- contrast(
    fit,
    list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    by = "G",
    ci_method = "bootstrap",
    n_boot = 40
  )
  expect_no_error(capture.output(print(res)))
  ci <- confint(res)
  expect_equal(nrow(ci), 2L * 2L * 3L)
  expect_true(all(ci[, "lower"] <= ci[, "upper"]))
})

test_that("confint handles a degenerate by-stratum for a multinomial result", {
  # A by-stratum with fewer than two successful bootstrap replicates takes the
  # NA fallback, which must emit the per-stratum row count (interventions x
  # classes), not the intervention count, or the CI rownames assignment aborts.
  d <- sim_multinom_binary(n = 6000, seed = 11)
  d$G <- factor(ifelse(d$L > 0, "hi", "lo"))
  fit <- causat(
    d,
    "Y",
    "A",
    confounders = ~ L + G,
    estimator = "gcomp",
    model_fn = nnet::multinom,
    trace = FALSE
  )
  res <- contrast(
    fit,
    list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    by = "G",
    ci_method = "bootstrap",
    n_boot = 40
  )
  # Force one stratum degenerate (a single retained replicate).
  res$boot_t[["hi"]] <- res$boot_t[["hi"]][1, , drop = FALSE]
  ci <- confint(res)
  expect_equal(nrow(ci), nrow(res$estimates))
  expect_equal(
    rownames(ci),
    paste(res$estimates$intervention, res$estimates$class, sep = ":")
  )
  expect_true(all(is.na(ci[res$estimates$by == "hi", ])))
})

test_that("external weights give a weighted multinomial g-formula", {
  skip_if_not_installed("marginaleffects")
  set.seed(7)
  d <- sim_multinom_binary(n = 8000, seed = 11)
  w <- stats::runif(nrow(d), 0.5, 2)
  fit <- causat(
    d,
    "Y",
    "A",
    confounders = ~L,
    estimator = "gcomp",
    model_fn = nnet::multinom,
    trace = FALSE,
    weights = w
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "bootstrap",
    n_boot = 40
  )
  nd <- as.data.frame(fit$data)
  nd$A <- 1
  me <- marginaleffects::avg_predictions(
    fit$model,
    newdata = nd,
    type = "probs",
    wts = w
  )
  me1 <- stats::setNames(me$estimate, as.character(me$group))[mo_labels3]
  mu1 <- res$estimates[intervention == "a1"]
  est1 <- stats::setNames(mu1$estimate, mu1$class)[mo_labels3]
  expect_equal(est1, me1, tolerance = 1e-6)
})

test_that("IPCW de-biases a MAR-censored multinomial outcome", {
  d <- sim_multinom_binary(n = 10000, seed = 11)
  # Cens = 1 marks a censored (missing) outcome; missingness depends on L
  # (MAR), so a complete-case fit would be biased and IPCW must repair it.
  d$Cens <- stats::rbinom(nrow(d), 1L, stats::plogis(-0.5 - 0.8 * d$L))
  d$Yobs <- d$Y
  d$Yobs[d$Cens == 1L] <- NA
  fit <- causat(
    d,
    "Yobs",
    "A",
    confounders = ~L,
    estimator = "gcomp",
    model_fn = nnet::multinom,
    trace = FALSE,
    censoring = "Cens",
    ipcw = TRUE,
    confounders_censoring = ~L
  )
  res <- contrast(
    fit,
    list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "bootstrap",
    n_boot = 40
  )
  # IPCW targets the full-population g-formula, so the estimate recovers the
  # softmax truth despite the MAR censoring.
  t1 <- mo_gcomp_truth(
    mo_gen_al_binary,
    mo_eta3,
    mo_labels3,
    function(A, L) rep(1, length(L))
  )
  mu1 <- res$estimates[intervention == "a1"]
  est1 <- stats::setNames(mu1$estimate, mu1$class)[mo_labels3]
  expect_equal(est1, t1[mo_labels3], tolerance = 0.03)
  s <- tapply(res$estimates$estimate, res$estimates$intervention, sum)
  expect_equal(as.vector(s), c(1, 1), tolerance = 1e-8)
})

test_that("subset estimand restricts the multinomial standardisation set", {
  skip_if_not_installed("marginaleffects")
  d <- sim_multinom_binary(n = 9000, seed = 11)
  fit <- causat(
    d,
    "Y",
    "A",
    confounders = ~L,
    estimator = "gcomp",
    model_fn = nnet::multinom,
    trace = FALSE
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    subset = quote(L > 0),
    ci_method = "bootstrap",
    n_boot = 40
  )
  expect_equal(res$estimand, "subset")
  sub <- fit$data[fit$data$L > 0, ]
  me1 <- me_multinom_means(fit$model, sub, "A", 1)
  mu1 <- res$estimates[intervention == "a1"]
  est1 <- stats::setNames(mu1$estimate, mu1$class)[mo_labels3]
  expect_equal(est1, me1[mo_labels3], tolerance = 1e-6)
})

test_that("more than two interventions produce K rows each", {
  d <- sim_multinom_continuous(n = 7000, seed = 12)
  fit <- causat(
    d,
    "Y",
    "A",
    confounders = ~L,
    estimator = "gcomp",
    model_fn = nnet::multinom,
    trace = FALSE
  )
  res <- contrast(
    fit,
    interventions = list(
      base = shift(0),
      up = shift(1),
      down = shift(-1)
    ),
    reference = "base",
    ci_method = "bootstrap",
    n_boot = 40
  )
  expect_equal(nrow(res$estimates), 3L * 3L)
  expect_equal(nrow(res$contrasts), 2L * 3L)
  s <- tapply(res$estimates$estimate, res$estimates$intervention, sum)
  expect_equal(as.vector(s), c(1, 1, 1), tolerance = 1e-8)
})

test_that("spline confounders fit and standardise for a multinomial outcome", {
  skip_if_not_installed("marginaleffects")
  d <- sim_multinom_binary(n = 8000, seed = 11)
  fit <- causat(
    d,
    "Y",
    "A",
    confounders = ~ splines::ns(L, 3),
    estimator = "gcomp",
    model_fn = nnet::multinom,
    trace = FALSE
  )
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "bootstrap",
    n_boot = 40
  )
  me1 <- me_multinom_means(fit$model, fit$data, "A", 1)
  mu1 <- res$estimates[intervention == "a1"]
  est1 <- stats::setNames(mu1$estimate, mu1$class)[mo_labels3]
  expect_equal(est1, me1[mo_labels3], tolerance = 1e-6)
})

# --- Bootstrap variance -----------------------------------------------------

test_that("bootstrap variance is reproducible and near the delta-method SE", {
  skip_if_fast()
  skip_if_not_installed("marginaleffects")
  d <- sim_multinom_binary(n = 6000, seed = 11)
  fit <- causat(
    d,
    "Y",
    "A",
    confounders = ~L,
    estimator = "gcomp",
    model_fn = nnet::multinom,
    trace = FALSE
  )
  ivs <- list(a1 = static(1), a0 = static(0))

  set.seed(123)
  res1 <- contrast(
    fit,
    ivs,
    reference = "a0",
    ci_method = "bootstrap",
    n_boot = 300
  )
  set.seed(123)
  res2 <- contrast(
    fit,
    ivs,
    reference = "a0",
    ci_method = "bootstrap",
    n_boot = 300
  )
  expect_equal(res1$estimates$se, res2$estimates$se)
  expect_true(all(is.finite(res1$estimates$se) & res1$estimates$se > 0))

  # Bootstrap SE within a loose ratio of the marginaleffects delta-method SE
  # (the tight equality is the job of the analytic sandwich chunk).
  me <- marginaleffects::avg_predictions(
    fit$model,
    newdata = transform(as.data.frame(fit$data), A = 1),
    type = "probs"
  )
  me_se <- stats::setNames(me$std.error, as.character(me$group))[mo_labels3]
  mu1 <- res1$estimates[intervention == "a1"]
  boot_se <- stats::setNames(mu1$se, mu1$class)[mo_labels3]
  ratio <- boot_se / me_se
  expect_true(all(ratio > 0.6 & ratio < 1.6))
})

test_that("confint() returns per-(intervention,class) percentile intervals", {
  d <- sim_multinom_binary(n = 6000, seed = 11)
  fit <- causat(
    d,
    "Y",
    "A",
    confounders = ~L,
    estimator = "gcomp",
    model_fn = nnet::multinom,
    trace = FALSE
  )
  res <- contrast(
    fit,
    list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "bootstrap",
    n_boot = 100
  )
  ci <- confint(res)
  expect_equal(nrow(ci), 2L * 3L)
  expect_equal(
    rownames(ci),
    paste(res$estimates$intervention, res$estimates$class, sep = ":")
  )
  expect_true(all(ci[, "lower"] < ci[, "upper"]))
})

# --- S3 layer ---------------------------------------------------------------

test_that("S3 methods render the multinomial result", {
  d <- sim_multinom_binary(n = 6000, seed = 11)
  fit <- causat(
    d,
    "Y",
    "A",
    confounders = ~L,
    estimator = "gcomp",
    model_fn = nnet::multinom,
    trace = FALSE
  )
  res <- contrast(
    fit,
    list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "bootstrap",
    n_boot = 60
  )

  expect_named(
    coef(res),
    paste(res$estimates$intervention, res$estimates$class, sep = ":")
  )

  td <- tidy(res, which = "all")
  expect_true("class" %in% names(td))
  expect_equal(sum(td$type == "mean"), 2L * 3L)
  expect_equal(sum(td$type == "contrast"), 3L)

  g <- glance(res)
  expect_equal(g$n_interventions, 2L)
  expect_equal(g$n_classes, 3L)

  out <- capture.output(print(res))
  expect_true(any(grepl("Class probabilities", out)))
  expect_true(any(grepl("Contrasts \\(per class\\)", out)))

  skip_if_not_installed("forrest")
  expect_no_error(plot(res, which = "means"))
  expect_no_error(plot(res, which = "contrasts"))
})

# --- Byte-identical guard for the scalar path ------------------------------

test_that("a scalar-outcome gcomp result is unchanged by the multinomial path", {
  set.seed(3)
  n <- 2000
  L <- stats::rnorm(n)
  A <- stats::rbinom(n, 1L, stats::plogis(0.3 * L))
  Y <- stats::rbinom(n, 1L, stats::plogis(-0.2 + 0.8 * A + 0.5 * L))
  d <- data.frame(Y = Y, A = A, L = L)
  fit <- causat(
    d,
    "Y",
    "A",
    confounders = ~L,
    estimator = "gcomp",
    family = "binomial"
  )
  res <- contrast(
    fit,
    list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich"
  )
  expect_null(res$class_labels)
  expect_false("class" %in% names(res$estimates))
  expect_true(is.matrix(res$vcov))
  expect_named(coef(res), c("a1", "a0"))
  expect_equal(nrow(res$estimates), 2L)
  expect_equal(nrow(res$contrasts), 1L)
})

# --- Rejection gates --------------------------------------------------------

test_that("categorical outcome rejects unsupported estimators at causat()", {
  d <- sim_multinom_binary(n = 500, seed = 11)

  expect_error(
    causat(d, "Y", "A", confounders = ~L, estimator = "ipw"),
    class = "causatr_categorical_outcome_unsupported"
  )
  expect_error(
    causat(d, "Y", "A", confounders = ~L, estimator = "aipw"),
    class = "causatr_categorical_outcome_unsupported"
  )
  expect_error(
    causat(d, "Y", "A", confounders = ~L, estimator = "matching"),
    class = "causatr_categorical_outcome_unsupported"
  )
  expect_error(
    causat(
      d,
      "Y",
      "A",
      confounders = ~L,
      estimator = "snm",
      treatment_free = ~L
    ),
    class = "causatr_snm_categorical_outcome"
  )
})

test_that("categorical outcome rejects longitudinal and transport designs", {
  set.seed(5)
  np <- 400
  dl <- data.frame(
    id = rep(seq_len(np), each = 2L),
    time = rep(1:2, np),
    A = stats::rbinom(np * 2L, 1L, 0.5),
    L = stats::rnorm(np * 2L)
  )
  ylevels <- c("none", "mild", "severe")
  yobs <- factor(
    sample(ylevels, np, replace = TRUE),
    levels = ylevels
  )
  dl$Y <- c(rbind(rep(NA_character_, np), as.character(yobs)))
  dl$Y <- factor(dl$Y, levels = ylevels)
  expect_error(
    causat(
      dl,
      "Y",
      "A",
      confounders = ~L,
      estimator = "gcomp",
      id = "id",
      time = "time"
    ),
    class = "causatr_categorical_outcome_unsupported"
  )

  d <- sim_multinom_binary(n = 600, seed = 11)
  d$S <- stats::rbinom(nrow(d), 1L, 0.5)
  expect_error(
    causat(
      d,
      "Y",
      "A",
      confounders = ~L,
      estimator = "gcomp",
      target = "S"
    ),
    class = "causatr_categorical_outcome_unsupported"
  )
})

test_that("categorical outcome rejects sandwich and stochastic at contrast()", {
  d <- sim_multinom_binary(n = 1500, seed = 11)
  fit <- causat(
    d,
    "Y",
    "A",
    confounders = ~L,
    estimator = "gcomp",
    model_fn = nnet::multinom,
    trace = FALSE
  )
  expect_error(
    contrast(
      fit,
      list(a1 = static(1), a0 = static(0)),
      ci_method = "sandwich"
    ),
    class = "causatr_categorical_outcome_sandwich"
  )
  expect_error(
    contrast(
      fit,
      list(
        s = stochastic(function(n) stats::rbinom(n, 1L, 0.5)),
        a0 = static(0)
      ),
      ci_method = "bootstrap",
      n_boot = 10
    ),
    class = "causatr_categorical_outcome_unsupported"
  )
})
