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
  expect_equal(fit$method, "gcomp")
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
