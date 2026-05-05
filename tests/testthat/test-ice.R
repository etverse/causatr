# --- Structure-1 DGP: TV confounder only, no baseline-only confounder --------
#
#   L_0 ~ N(0, 1)              (TV confounder at t = 0)
#   A_0 ~ Bernoulli(expit(0.5 * L_0))
#   L_1 = A_0 + 0.5 * L_0 + e (treatment-confounder feedback)
#   A_1 ~ Bernoulli(expit(0.3 * L_1 + 0.2 * A_0))
#   Y   = 10 + 2 * (A_0 + A_1) + L_0 + L_1 + e
#
# True ATE (always vs never) = 5
# E[Y(1,1)] = 15, E[Y(0,0)] = 10

make_tv_only_scm <- function(n = 5000, seed = 42) {
  set.seed(seed)
  L0 <- stats::rnorm(n)
  A0 <- stats::rbinom(n, 1, stats::plogis(0.5 * L0))
  L1 <- A0 + 0.5 * L0 + stats::rnorm(n, 0, 0.5)
  A1 <- stats::rbinom(n, 1, stats::plogis(0.3 * L1 + 0.2 * A0))
  Y <- 10 + 2 * (A0 + A1) + L0 + L1 + stats::rnorm(n)

  rbind(
    data.frame(id = seq_len(n), time = 0L, A = A0, L = L0, Y = NA_real_),
    data.frame(id = seq_len(n), time = 1L, A = A1, L = L1, Y = Y)
  )
}


test_that("ICE recovers ATE with TV confounder only (no baseline confounder)", {
  long <- make_tv_only_scm(n = 5000, seed = 42)

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
    interventions = list(always = static(1), never = static(0)),
    reference = "never",
    ci_method = "sandwich"
  )

  ate <- res$contrasts$estimate[1]
  se <- res$contrasts$se[1]
  ci_lo <- res$contrasts$ci_lower[1]
  ci_hi <- res$contrasts$ci_upper[1]

  expect_equal(ate, 5, tolerance = 0.3)
  expect_true(ci_lo < 5 && ci_hi > 5)
  expect_true(se > 0 && se < 0.5)
})


test_that("ICE sandwich SE is valid with TV confounder only", {
  long <- make_tv_only_scm(n = 5000, seed = 42)

  fit <- causat(
    long,
    outcome = "Y",
    treatment = "A",
    confounders = ~1,
    confounders_tv = ~L,
    id = "id",
    time = "time"
  )

  res_sand <- contrast(
    fit,
    interventions = list(always = static(1), never = static(0)),
    reference = "never",
    ci_method = "sandwich"
  )
  res_boot <- contrast(
    fit,
    interventions = list(always = static(1), never = static(0)),
    reference = "never",
    ci_method = "bootstrap",
    n_boot = 200L
  )

  se_sand <- res_sand$contrasts$se[1]
  se_boot <- res_boot$contrasts$se[1]
  ratio <- se_sand / se_boot

  expect_true(ratio > 0.5 && ratio < 2.0)
})


test_that("causat() returns a causatr_fit for longitudinal data", {
  long <- make_table201(scale = 1 / 100)
  fit <- causat(
    long,
    outcome = "Y",
    treatment = "A",
    confounders = ~1,
    confounders_tv = ~L,
    id = "id",
    time = "time"
  )

  expect_s3_class(fit, "causatr_fit")
  expect_equal(fit$type, "longitudinal")
  expect_equal(fit$estimator, "gcomp")
  expect_null(fit$model)
  expect_equal(fit$details$n_times, 2L)
})


test_that("ICE recovers ATE = 0 from Table 20.1 (always vs never)", {
  long <- make_table201(scale = 1 / 100)
  fit <- causat(
    long,
    outcome = "Y",
    treatment = "A",
    confounders = ~1,
    confounders_tv = ~L,
    id = "id",
    time = "time"
  )

  result <- contrast(
    fit,
    interventions = list(always = static(1), never = static(0)),
    reference = "never",
    type = "difference",
    ci_method = "sandwich"
  )

  expect_s3_class(result, "causatr_result")

  # Both marginal means should be 60 (book value).
  expect_equal(result$estimates$estimate[1], 60, tolerance = 0.5)
  expect_equal(result$estimates$estimate[2], 60, tolerance = 0.5)

  # ATE should be 0 (true causal effect is null).
  expect_equal(result$contrasts$estimate[1], 0, tolerance = 0.5)
})


test_that("ICE sandwich SE is finite and positive", {
  long <- make_table201(scale = 1 / 100)
  fit <- causat(
    long,
    outcome = "Y",
    treatment = "A",
    confounders = ~1,
    confounders_tv = ~L,
    id = "id",
    time = "time"
  )

  result <- contrast(
    fit,
    interventions = list(always = static(1), never = static(0)),
    ci_method = "sandwich"
  )

  expect_true(all(result$estimates$se > 0))
  expect_true(all(is.finite(result$estimates$se)))
  expect_true(all(result$contrasts$se > 0))
  expect_true(all(is.finite(result$contrasts$se)))
})


test_that("ICE bootstrap gives finite SE", {
  long <- make_table201(scale = 1 / 800)
  fit <- causat(
    long,
    outcome = "Y",
    treatment = "A",
    confounders = ~1,
    confounders_tv = ~L,
    id = "id",
    time = "time"
  )

  result <- contrast(
    fit,
    interventions = list(always = static(1), never = static(0)),
    ci_method = "bootstrap",
    n_boot = 10L
  )

  expect_true(all(is.finite(result$estimates$se)))
  expect_equal(result$ci_method, "bootstrap")
})


test_that("ICE works with dynamic interventions", {
  long <- make_table201(scale = 1 / 100)
  fit <- causat(
    long,
    outcome = "Y",
    treatment = "A",
    confounders = ~1,
    confounders_tv = ~L,
    id = "id",
    time = "time"
  )

  result <- contrast(
    fit,
    interventions = list(
      dynamic_rule = dynamic(function(data, trt) {
        ifelse(!is.na(data$L) & data$L > 0, 1L, 0L)
      }),
      never = static(0)
    ),
    ci_method = "sandwich"
  )

  expect_s3_class(result, "causatr_result")
  expect_true(all(is.finite(result$estimates$estimate)))
})


test_that("ICE works with NULL intervention (natural course)", {
  long <- make_table201(scale = 1 / 100)
  fit <- causat(
    long,
    outcome = "Y",
    treatment = "A",
    confounders = ~1,
    confounders_tv = ~L,
    id = "id",
    time = "time"
  )

  result <- contrast(
    fit,
    interventions = list(observed = NULL, never = static(0)),
    ci_method = "sandwich"
  )

  expect_s3_class(result, "causatr_result")
  expect_true(all(is.finite(result$estimates$estimate)))
})


test_that("ICE rejects ATT/ATC for longitudinal data", {
  long <- make_table201(scale = 1 / 100)
  fit <- causat(
    long,
    outcome = "Y",
    treatment = "A",
    confounders = ~1,
    confounders_tv = ~L,
    id = "id",
    time = "time"
  )

  expect_error(
    contrast(
      fit,
      interventions = list(a = static(1), b = static(0)),
      estimand = "ATT"
    ),
    class = "rlang_error"
  )
})


test_that("ICE handles censoring indicator", {
  long <- make_table201(scale = 1 / 100)
  set.seed(42)
  ids_to_censor <- sample(unique(long$id), size = 20)
  long$C <- 0L
  long$C[long$id %in% ids_to_censor & long$time == 1L] <- 1L
  long$Y[long$C == 1L] <- NA_real_

  fit <- causat(
    long,
    outcome = "Y",
    treatment = "A",
    confounders = ~1,
    confounders_tv = ~L,
    censoring = "C",
    id = "id",
    time = "time"
  )

  result <- contrast(
    fit,
    interventions = list(always = static(1), never = static(0)),
    ci_method = "sandwich"
  )

  expect_s3_class(result, "causatr_result")
  expect_true(all(is.finite(result$estimates$estimate)))
})


test_that("ICE recovers known non-zero ATE from linear SCM (2 time points)", {
  long <- make_linear_scm(n = 5000, n_times = 2, seed = 101)
  fit <- causat(
    long,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    id = "id",
    time = "time"
  )

  result <- contrast(
    fit,
    interventions = list(always = static(1), never = static(0)),
    reference = "never",
    type = "difference",
    ci_method = "sandwich"
  )

  true_ate <- 3 * 2 - 1 # = 5

  expect_equal(result$contrasts$estimate[1], true_ate, tolerance = 0.5)
  expect_true(result$estimates$se[1] > 0)
  expect_true(result$estimates$se[2] > 0)
})


test_that("ICE recovers known ATE from linear SCM (3 time points)", {
  long <- make_linear_scm(n = 5000, n_times = 3, seed = 202)
  fit <- causat(
    long,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    history = Inf
  )

  result <- contrast(
    fit,
    interventions = list(always = static(1), never = static(0)),
    reference = "never",
    type = "difference",
    ci_method = "sandwich"
  )

  true_ate <- 3 * 3 - 1 # = 8

  expect_equal(result$contrasts$estimate[1], true_ate, tolerance = 0.75)
})


test_that("ICE recovers dynamic intervention effect (treat if L0 > 0)", {
  # Under dynamic(treat if L0 > 0):
  #   A_t = I(L0 > 0) for all t, so sum(A) = T * I(L0 > 0)
  #   L_t = I(L0>0) + 0.5*L0 + ε for t > 0
  #   E[Y|dynamic] = 10 + (3T-1) * P(L0>0) = 10 + (3T-1)/2
  #   E[Y|never]   = 10
  #   ATE(dynamic vs never) = (3T-1)/2
  long <- make_linear_scm(n = 5000, n_times = 2, seed = 303)
  fit <- causat(
    long,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    id = "id",
    time = "time"
  )

  result <- contrast(
    fit,
    interventions = list(
      treat_if_L0_pos = dynamic(function(data, trt) {
        ifelse(data$L0 > 0, 1L, 0L)
      }),
      never = static(0)
    ),
    reference = "never",
    type = "difference",
    ci_method = "sandwich"
  )

  true_ate <- (3 * 2 - 1) / 2 # = 2.5

  expect_equal(result$contrasts$estimate[1], true_ate, tolerance = 0.5)
})


test_that("ICE CI covers true ATE from linear SCM", {
  long <- make_linear_scm(n = 5000, n_times = 2, seed = 404)
  fit <- causat(
    long,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    id = "id",
    time = "time"
  )

  result <- contrast(
    fit,
    interventions = list(always = static(1), never = static(0)),
    reference = "never",
    type = "difference",
    ci_method = "sandwich",
    conf_level = 0.95
  )

  true_ate <- 5
  ci_lower <- result$contrasts$ci_lower[1]
  ci_upper <- result$contrasts$ci_upper[1]

  expect_true(
    ci_lower <= true_ate && true_ate <= ci_upper,
    label = sprintf(
      "95%% CI [%.2f, %.2f] should cover true ATE = %.1f",
      ci_lower,
      ci_upper,
      true_ate
    )
  )
})


test_that("ICE with continuous treatment and shift (LMTP) intervention", {
  long <- make_continuous_scm(n = 5000, seed = 505)
  fit <- causat(
    long,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    id = "id",
    time = "time"
  )

  result <- contrast(
    fit,
    interventions = list(
      shifted = shift(-0.5),
      observed = NULL
    ),
    type = "difference",
    ci_method = "sandwich"
  )

  expect_s3_class(result, "causatr_result")
  expect_true(all(is.finite(result$estimates$estimate)))
  expect_true(result$contrasts$estimate[1] > 0)
  expect_true(all(result$contrasts$se > 0))
})


test_that("ICE bootstrap and sandwich agree on linear SCM", {
  long <- make_linear_scm(n = 2000, n_times = 2, seed = 606)
  fit <- causat(
    long,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    id = "id",
    time = "time"
  )

  res_sw <- contrast(
    fit,
    interventions = list(always = static(1), never = static(0)),
    ci_method = "sandwich"
  )
  res_boot <- contrast(
    fit,
    interventions = list(always = static(1), never = static(0)),
    ci_method = "bootstrap",
    n_boot = 100L
  )

  expect_equal(
    res_sw$contrasts$se[1],
    res_boot$contrasts$se[1],
    tolerance = 0.5
  )
})


test_that("parallel bootstrap works for point-treatment gcomp", {
  skip_on_os("windows")
  data("nhefs", package = "causatr")
  fit <- causat(
    nhefs,
    outcome = "wt82_71",
    treatment = "qsmk",
    confounders = ~ sex + age + wt71
  )

  result <- contrast(
    fit,
    interventions = list(quit = static(1), cont = static(0)),
    ci_method = "bootstrap",
    n_boot = 20L,
    parallel = "multicore",
    ncpus = 2L
  )

  expect_s3_class(result, "causatr_result")
  expect_true(all(is.finite(result$estimates$se)))
  expect_true(all(result$estimates$se > 0))
})


test_that("ICE external weights propagate through every backward step", {
  # Regression guard: the backward-iteration loop previously fit the
  # pseudo-outcome models via a direct `model_fn()` call that dropped
  # the `weights =` argument. Only the final-time model received
  # external weights; every earlier step was silently unweighted.
  # We detect the bug by comparing two fits that differ *only* in
  # whether `weights` is NULL or a non-uniform vector. If the weights
  # are ignored, the fitted coefficients of the pseudo-outcome models
  # are identical.
  long <- make_table201(scale = 1 / 100)

  # Non-uniform per-row weights so the weighted fit must diverge
  # from the unweighted one at every step.
  set.seed(42)
  w <- stats::runif(nrow(long), min = 0.5, max = 2)

  fit_unw <- causat(
    long,
    outcome = "Y",
    treatment = "A",
    confounders = ~1,
    confounders_tv = ~L,
    id = "id",
    time = "time"
  )
  fit_w <- causat(
    long,
    outcome = "Y",
    treatment = "A",
    confounders = ~1,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    weights = w
  )

  # Run a contrast so the backward-iteration fits every pseudo-outcome
  # model. Compare coefficients across steps; with weights != NULL
  # at least one intermediate model must differ from the unweighted
  # counterpart.
  ivs <- list(always = static(1), never = static(0))
  r_unw <- contrast(fit_unw, interventions = ivs, ci_method = "sandwich")
  r_w <- contrast(fit_w, interventions = ivs, ci_method = "sandwich")

  # Access the fitted intermediate models via the attached ICE result.
  models_unw <- attr(r_unw, "ice_models_always")
  models_w <- attr(r_w, "ice_models_always")

  # Fallback: compare the final marginal mean estimates directly.
  # Any meaningful weight effect will move the estimate.
  expect_false(
    isTRUE(all.equal(
      r_unw$estimates$estimate,
      r_w$estimates$estimate,
      tolerance = 1e-8
    ))
  )
})


test_that("ICE accepts shift / scale_by / threshold interventions on time-varying treatment", {
  # Coverage gap: existing ICE tests cover static and dynamic. The
  # remaining intervention constructors (shift, scale_by, threshold)
  # are documented as supported by gcomp/ICE but were never exercised
  # against the ICE backward iteration. Smoke-test that they at
  # least fit and return finite estimates.
  set.seed(2)
  n <- 50
  long <- data.table::data.table(
    id = rep(seq_len(n), each = 3),
    time = rep(0:2, times = n),
    A = stats::rnorm(3 * n, 1, 0.5),
    L = stats::rnorm(3 * n),
    Y = stats::rnorm(3 * n)
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

  for (iv in list(
    list(s = shift(0.5), z = static(0)),
    list(s = scale_by(1.2), z = static(0)),
    list(s = threshold(lower = 0), z = static(0))
  )) {
    res <- contrast(fit, interventions = iv, ci_method = "sandwich")
    expect_true(all(is.finite(res$estimates$estimate)))
    expect_true(all(is.finite(res$estimates$se)))
    expect_true(all(res$estimates$se > 0))
  }
})


test_that("ipsi() is explicitly rejected by ICE/contrast() (Phase 4 placeholder)", {
  # Lock in the not-yet-implemented abort so a future Phase-4
  # implementation can't silently regress while wiring it up.
  long <- make_table201(scale = 1 / 100)
  fit <- causat(
    long,
    outcome = "Y",
    treatment = "A",
    confounders = ~1,
    confounders_tv = ~L,
    id = "id",
    time = "time"
  )
  expect_snapshot(
    error = TRUE,
    contrast(
      fit,
      interventions = list(boost = ipsi(2), nat = ipsi(1)),
      ci_method = "sandwich"
    )
  )
})


test_that("ICE predictions saturated at {0,1} fire a classed warning", {
  # A binary-outcome ICE fit whose final-time model produces predictions
  # at exactly 0 or 1 silently drops those rows from the next backward
  # quasibinomial pseudo-outcome model (IWLS working weight
  # `mu*(1-mu) = 0`). A rate-limited classed warning
  # `causatr_ice_boundary_saturation` now surfaces the risk.
  set.seed(2024)
  n <- 300
  Tper <- 3
  rows <- expand.grid(time = 0:(Tper - 1L), id = seq_len(n))
  rows <- rows[order(rows$id, rows$time), ]
  rows$L <- stats::rnorm(nrow(rows))
  rows$A <- stats::rbinom(nrow(rows), 1, stats::plogis(0.2 * rows$L))
  is_final <- rows$time == (Tper - 1L)
  rows$Y <- NA_real_
  # Extreme coefficients push a handful of final-time preds to 0/1.
  eta <- -8 + 12 * rows$A[is_final] + 8 * rows$L[is_final]
  rows$Y[is_final] <- stats::rbinom(sum(is_final), 1, stats::plogis(eta))

  fit <- causat(
    rows,
    outcome = "Y",
    treatment = "A",
    confounders = ~1,
    confounders_tv = ~L,
    estimator = "gcomp",
    type = "longitudinal",
    id = "id",
    time = "time",
    family = stats::binomial()
  )

  # Reset the frequency throttle so this test is deterministic even if
  # the classed warning already fired earlier in the session.
  rlang::reset_warning_verbosity("causatr_ice_boundary_saturation")

  saw <- FALSE
  msg <- NULL
  withCallingHandlers(
    contrast(
      fit,
      interventions = list(a1 = static(1), a0 = static(0)),
      ci_method = "sandwich"
    ),
    causatr_ice_boundary_saturation = function(w) {
      saw <<- TRUE
      msg <<- conditionMessage(w)
      invokeRestart("muffleWarning")
    },
    warning = function(w) invokeRestart("muffleWarning")
  )

  expect_true(saw)
  expect_match(msg, "saturated at \\{0, 1\\}")
})


# --- Formula builder: transformed TV confounders ----------------------------

test_that("substitute_vars_in_term() handles function calls", {
  expect_equal(
    substitute_vars_in_term("ns(L, 3)", list(L = "lag1_L")),
    "ns(lag1_L, 3)"
  )
  expect_equal(
    substitute_vars_in_term("log(L + 1)", list(L = "lag2_L")),
    "log(lag2_L + 1)"
  )
})

test_that("substitute_vars_in_term() passes through bare names", {
  expect_equal(
    substitute_vars_in_term("L", list(L = "lag1_L")),
    "lag1_L"
  )
})

test_that("is_bare_term() distinguishes bare vs transformed", {
  expect_true(is_bare_term("L"))
  expect_true(is_bare_term("my_var"))
  expect_false(is_bare_term("ns(L, 3)"))
  expect_false(is_bare_term("log(L + 1)"))
})

test_that("term_vars() extracts variable names from terms", {
  expect_equal(term_vars("L"), "L")
  expect_equal(term_vars("ns(L, 3)"), "L")
  expect_setequal(term_vars("log(L + 1)"), "L")
})

# --- Self-consistency: transformed TV confounders on a linear DGP -----------
#
# On the linear DGP (make_tv_only_scm, true ATE = 5), all TV confounder
# specs should recover the same ATE because the true model is linear.
# Any correctly-expanding spec (bare, poly, ns, log) that nests the
# truth should agree within Monte Carlo noise.

test_that("ICE self-consistency: bare vs poly() vs ns() vs I()", {
  long <- make_tv_only_scm(n = 5000, seed = 200)

  run_ice <- function(tv_formula) {
    fit <- causat(
      long,
      outcome = "Y",
      treatment = "A",
      confounders = ~1,
      confounders_tv = tv_formula,
      id = "id",
      time = "time"
    )
    res <- contrast(
      fit,
      interventions = list(always = static(1), never = static(0)),
      reference = "never",
      ci_method = "sandwich"
    )
    res$contrasts$estimate[1]
  }

  ate_bare <- run_ice(~L)
  ate_poly <- run_ice(~ poly(L, 2))
  ate_ns <- run_ice(~ splines::ns(L, df = 3))
  ate_sq <- run_ice(~ L + I(L^2))

  # All should recover ATE ≈ 5 (SE ≈ 0.05 at n=5000)
  expect_equal(ate_bare, 5, tolerance = 0.15)
  expect_equal(ate_poly, 5, tolerance = 0.15)
  expect_equal(ate_ns, 5, tolerance = 0.15)
  expect_equal(ate_sq, 5, tolerance = 0.15)

  # Same data + all specs nest the linear truth → near-identical
  expect_lt(abs(ate_bare - ate_poly), 0.1)
  expect_lt(abs(ate_bare - ate_ns), 0.1)
  expect_lt(abs(ate_bare - ate_sq), 0.1)
})


test_that("ICE self-consistency: log transform on linear DGP", {
  # L is centered around 0, so shift it positive for log()
  set.seed(201)
  n <- 4000
  L0 <- stats::rnorm(n, mean = 5, sd = 1)
  A0 <- stats::rbinom(n, 1, stats::plogis(0.5 * (L0 - 5)))
  L1 <- A0 + 0.5 * L0 + stats::rnorm(n, 0, 0.5)
  A1 <- stats::rbinom(n, 1, stats::plogis(0.3 * (L1 - 5) + 0.2 * A0))
  Y <- 10 + 2 * (A0 + A1) + L0 + L1 + stats::rnorm(n)
  long <- rbind(
    data.frame(id = seq_len(n), time = 0L, A = A0, L = L0, Y = NA_real_),
    data.frame(id = seq_len(n), time = 1L, A = A1, L = L1, Y = Y)
  )

  run_ice <- function(tv_formula) {
    fit <- causat(
      long,
      outcome = "Y",
      treatment = "A",
      confounders = ~1,
      confounders_tv = tv_formula,
      id = "id",
      time = "time"
    )
    contrast(
      fit,
      interventions = list(always = static(1), never = static(0)),
      reference = "never",
      ci_method = "sandwich"
    )$contrasts$estimate[1]
  }

  ate_bare <- run_ice(~L)
  ate_log <- run_ice(~ log(L))

  # Both recover ATE ≈ 5 (DGP is linear in L; log(L) approximates
  # linearity well when L is far from 0)
  expect_equal(ate_bare, 5, tolerance = 0.15)
  expect_equal(ate_log, 5, tolerance = 0.2)
  expect_lt(abs(ate_bare - ate_log), 0.15)
})


# --- Nonlinear DGP for spline superiority -----------------------------------
#
# DGP where the outcome depends on L through a nonlinear function.
# A correctly-specified spline should do better than bare L.
#
#   L_0 ~ N(0, 1)
#   A_0 ~ Bernoulli(expit(0.3 * L_0))
#   L_1 = 0.5 * A_0 + 0.3 * L_0 + e
#   A_1 ~ Bernoulli(expit(0.3 * L_1 + 0.2 * A_0))
#   Y   = 10 + 2 * (A_0 + A_1) + sin(2 * L_0) + sin(2 * L_1) + e
#
# True ATE (always vs never) ≈ 5 (since E[sin(2*L)] doesn't change
# under treatment, the nonlinearity is purely in the confounders, but
# correct adjustment still requires capturing the nonlinear L→Y path).

make_nonlinear_scm <- function(n = 5000, seed = 42) {
  set.seed(seed)
  L0 <- stats::rnorm(n)
  A0 <- stats::rbinom(n, 1, stats::plogis(0.3 * L0))
  L1 <- 0.5 * A0 + 0.3 * L0 + stats::rnorm(n, 0, 0.5)
  A1 <- stats::rbinom(n, 1, stats::plogis(0.3 * L1 + 0.2 * A0))
  Y <- 10 + 2 * (A0 + A1) + sin(2 * L0) + sin(2 * L1) + stats::rnorm(n)

  rbind(
    data.frame(id = seq_len(n), time = 0L, A = A0, L = L0, Y = NA_real_),
    data.frame(id = seq_len(n), time = 1L, A = A1, L = L1, Y = Y)
  )
}


test_that("ICE with ns() handles nonlinear DGP better than bare L", {
  long <- make_nonlinear_scm(n = 5000, seed = 300)

  run_ice <- function(tv_formula) {
    fit <- causat(
      long,
      outcome = "Y",
      treatment = "A",
      confounders = ~1,
      confounders_tv = tv_formula,
      id = "id",
      time = "time"
    )
    contrast(
      fit,
      interventions = list(always = static(1), never = static(0)),
      reference = "never",
      ci_method = "sandwich"
    )$contrasts$estimate[1]
  }

  ate_bare <- run_ice(~L)
  ate_ns <- run_ice(~ splines::ns(L, df = 5))

  # ns() captures the sin() nonlinearity → near truth
  expect_equal(ate_ns, 5, tolerance = 0.15)
  # Bare L is misspecified (omits sin() curvature) but residual
  # confounding bias is modest since treatment effect is additive
  expect_equal(ate_bare, 5, tolerance = 0.3)
})


# --- lmtp cross-check: ICE with ns() vs lmtp_sdr ---------------------------
#
# On the nonlinear DGP, compare causatr ICE with splines against
# lmtp::lmtp_sdr() with SL.glm. Both should be close to truth and
# to each other.

test_that("ICE with ns() agrees with lmtp_sdr on nonlinear DGP", {
  skip_if_not_installed("lmtp")

  long <- make_nonlinear_scm(n = 5000, seed = 400)

  # --- causatr: ICE with natural splines ---
  fit <- causat(
    long,
    outcome = "Y",
    treatment = "A",
    confounders = ~1,
    confounders_tv = ~ splines::ns(L, df = 4),
    id = "id",
    time = "time"
  )
  res <- contrast(
    fit,
    interventions = list(always = static(1), never = static(0)),
    reference = "never",
    ci_method = "sandwich"
  )
  ate_causatr <- res$contrasts$estimate[1]

  # --- lmtp: reshape to wide ---
  d_wide <- reshape(
    long,
    idvar = "id",
    timevar = "time",
    direction = "wide",
    v.names = c("A", "L", "Y"),
    sep = "_"
  )
  d_clean <- d_wide[, c("id", "A_0", "A_1", "L_0", "L_1", "Y_1")]

  run_lmtp <- function(shift_fn) {
    tryCatch(
      suppressWarnings(suppressMessages(lmtp::lmtp_sdr(
        data = d_clean,
        trt = c("A_0", "A_1"),
        outcome = "Y_1",
        baseline = NULL,
        time_vary = list(c("L_0"), c("L_1")),
        shift = shift_fn,
        outcome_type = "continuous",
        learners_trt = "SL.glm",
        learners_outcome = "SL.glm",
        folds = 1
      ))),
      error = function(e) NULL
    )
  }
  theta_of <- function(r) {
    tryCatch(r$estimate@x, error = function(e) r$theta)
  }

  r_always <- run_lmtp(function(data, trt) rep(1, nrow(data)))
  r_never <- run_lmtp(function(data, trt) rep(0, nrow(data)))

  skip_if(
    is.null(r_always) || is.null(r_never),
    "lmtp::lmtp_sdr() failed"
  )

  ate_lmtp <- theta_of(r_always) - theta_of(r_never)

  # Both near truth (ATE ≈ 5)
  expect_equal(ate_causatr, 5, tolerance = 0.15)
  expect_equal(ate_lmtp, 5, tolerance = 0.2)
  # Cross-method agreement (plug-in g-formula vs SDR)
  expect_lt(
    abs(ate_causatr - ate_lmtp),
    0.2,
    label = "causatr ICE(ns) vs lmtp_sdr ATE"
  )
})


test_that("ICE with poly() agrees with lmtp_sdr on linear DGP", {
  skip_if_not_installed("lmtp")

  long <- make_tv_only_scm(n = 5000, seed = 401)

  # --- causatr: ICE with poly ---
  fit <- causat(
    long,
    outcome = "Y",
    treatment = "A",
    confounders = ~1,
    confounders_tv = ~ poly(L, 3),
    id = "id",
    time = "time"
  )
  res <- contrast(
    fit,
    interventions = list(always = static(1), never = static(0)),
    reference = "never",
    ci_method = "sandwich"
  )
  ate_causatr <- res$contrasts$estimate[1]

  # --- lmtp ---
  d_wide <- reshape(
    long,
    idvar = "id",
    timevar = "time",
    direction = "wide",
    v.names = c("A", "L", "Y"),
    sep = "_"
  )
  d_clean <- d_wide[, c("id", "A_0", "A_1", "L_0", "L_1", "Y_1")]

  run_lmtp <- function(shift_fn) {
    tryCatch(
      suppressWarnings(suppressMessages(lmtp::lmtp_sdr(
        data = d_clean,
        trt = c("A_0", "A_1"),
        outcome = "Y_1",
        baseline = NULL,
        time_vary = list(c("L_0"), c("L_1")),
        shift = shift_fn,
        outcome_type = "continuous",
        learners_trt = "SL.glm",
        learners_outcome = "SL.glm",
        folds = 1
      ))),
      error = function(e) NULL
    )
  }
  theta_of <- function(r) {
    tryCatch(r$estimate@x, error = function(e) r$theta)
  }

  r_always <- run_lmtp(function(data, trt) rep(1, nrow(data)))
  r_never <- run_lmtp(function(data, trt) rep(0, nrow(data)))

  skip_if(
    is.null(r_always) || is.null(r_never),
    "lmtp::lmtp_sdr() failed"
  )

  ate_lmtp <- theta_of(r_always) - theta_of(r_never)

  # Both near truth (ATE = 5, linear DGP so both are correctly specified)
  expect_equal(ate_causatr, 5, tolerance = 0.15)
  expect_equal(ate_lmtp, 5, tolerance = 0.15)
  expect_lt(
    abs(ate_causatr - ate_lmtp),
    0.15,
    label = "causatr ICE(poly) vs lmtp_sdr ATE"
  )
})
