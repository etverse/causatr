# Stratified ICE (Phase 22a)
#
# `causat(..., stratified = "G")` fits a separate per-step outcome /
# pseudo-outcome model for each level of a baseline stratum G, then merges
# per-stratum predictions back into the backward ICE recursion. Because the
# stratification is on a baseline (time-invariant) variable and individuals
# never cross strata in the recursion, stratified ICE on the full panel is
# *exactly* equivalent to running pooled ICE separately on each stratum
# subset and concatenating the individual counterfactual pseudo-outcomes.
# That equivalence is the primary internal oracle below; truth-based DGPs
# and an lmtp cross-check confirm the causal target.

# Per-individual counterfactual pseudo-outcomes at baseline (one per id,
# ordered by first-period id). Computed directly so the point-estimate
# tests avoid the bootstrap (stratified ICE rejects the sandwich).
ice_pf <- function(fit, iv) {
  ice_iterate(fit, iv)$pseudo_final
}


# Test A — exact internal oracle ------------------------------------------
#
# Stratified ICE == pooled ICE run separately per stratum subset, to
# machine precision, for a gaussian outcome (per-stratum OLS fitted values
# equal the subset OLS fitted values exactly). The stratifying column `sex`
# is also placed in `confounders` so the per-stratum formula exercises
# `strip_stratum_terms()` (it must drop the constant `sex` term, leaving the
# same formula the subset fit builds).
test_that("stratified ICE equals per-stratum pooled ICE exactly (gaussian)", {
  d <- make_em_ice_scm(n = 2000, n_times = 2, seed = 7)
  sexb <- d$sex[d$time == 0]

  fit_s <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L0 + sex,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    stratified = "sex"
  )
  fit0 <- causat(
    d[d$sex == 0, ],
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    id = "id",
    time = "time"
  )
  fit1 <- causat(
    d[d$sex == 1, ],
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    id = "id",
    time = "time"
  )

  for (iv in list(static(1), static(0), shift(0))) {
    pf_s <- ice_pf(fit_s, iv)
    expect_equal(pf_s[sexb == 0], ice_pf(fit0, iv), tolerance = 1e-8)
    expect_equal(pf_s[sexb == 1], ice_pf(fit1, iv), tolerance = 1e-8)
  }

  # The stored model list is one model per stratum at each time point.
  ic <- ice_iterate(fit_s, static(1))
  expect_named(ic$models[[1]], c("0", "1"))
  expect_named(ic$models[[2]], c("0", "1"))
  expect_s3_class(ic$models[[2]][["0"]], "glm")
})


# Test B — truth-based recovery, binary outcome (non-collapsible) ----------
#
# make_em_ice_binom_scm has a sex-modified treatment effect on the logit
# scale. Documented Monte-Carlo truths:
#   RD | sex = 0 ~ 0.495,  RD | sex = 1 ~ 0.663,  marginal RD ~ 0.58.
# Stratified ICE is correctly specified within each stratum and recovers
# these; the naive pooled additive model (omitting the sex x A interaction)
# is misspecified and, under the non-collapsible logit link, its
# standardized marginal RD is biased away from truth.
test_that("stratified ICE recovers sex-modified RD; pooled additive is biased", {
  d <- make_em_ice_binom_scm(n = 8000, seed = 42)
  sexb <- d$sex[d$time == 0]

  fit_s <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L0 + sex,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    family = "binomial",
    stratified = "sex"
  )
  pf1 <- ice_pf(fit_s, static(1))
  pf0 <- ice_pf(fit_s, static(0))

  rd_marg <- mean(pf1) - mean(pf0)
  rd0 <- mean(pf1[sexb == 0]) - mean(pf0[sexb == 0])
  rd1 <- mean(pf1[sexb == 1]) - mean(pf0[sexb == 1])

  expect_equal(rd_marg, 0.58, tolerance = 0.03)
  expect_equal(rd0, 0.495, tolerance = 0.04)
  expect_equal(rd1, 0.663, tolerance = 0.04)

  # Naive pooled additive ICE (no per-stratum models, no sex:A term).
  fit_p <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L0 + sex,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    family = "binomial"
  )
  rd_pooled <- mean(ice_pf(fit_p, static(1))) - mean(ice_pf(fit_p, static(0)))

  # Stratified is strictly closer to the marginal truth than the
  # misspecified pooled additive model.
  expect_gt(abs(rd_pooled - 0.58), abs(rd_marg - 0.58))
})


# Test C — external oracle: lmtp per-stratum cross-check ------------------
#
# Validates the stratified-ICE per-stratum counterfactual means against
# lmtp::lmtp_sdr() fit separately on each stratum subset (the semiparametric
# analogue of what stratified ICE does parametrically). Tier-2 + lmtp guard.
test_that("stratified ICE per-stratum RD agrees with lmtp", {
  skip_if_not_installed("lmtp")
  skip_if_fast()

  d <- make_em_ice_binom_scm(n = 4000, seed = 42)
  sexb <- d$sex[d$time == 0]

  fit_s <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    family = "binomial",
    stratified = "sex"
  )
  pf1 <- ice_pf(fit_s, static(1))
  pf0 <- ice_pf(fit_s, static(0))

  d_wide <- reshape(
    d,
    idvar = "id",
    timevar = "time",
    direction = "wide",
    v.names = c("A", "L", "Y"),
    sep = "_"
  )
  d_clean <- d_wide[, c("id", "L0", "A_0", "A_1", "L_1", "Y_1")]
  sex_wide <- d_wide$sex

  run_lmtp <- function(data_sub, val) {
    suppressWarnings(suppressMessages(lmtp::lmtp_sdr(
      data = data_sub,
      trt = c("A_0", "A_1"),
      outcome = "Y_1",
      baseline = "L0",
      time_vary = list(NULL, "L_1"),
      shift = function(data, trt) rep(val, nrow(data)),
      outcome_type = "binomial",
      learners_trt = "SL.glm",
      learners_outcome = "SL.glm",
      folds = 1
    )))
  }

  for (s in c(0, 1)) {
    rd_causatr <- mean(pf1[sexb == s]) - mean(pf0[sexb == s])
    ds <- d_clean[sex_wide == s, ]
    rd_lmtp <- run_lmtp(ds, 1)$estimate@x - run_lmtp(ds, 0)$estimate@x
    expect_lt(
      abs(rd_causatr - rd_lmtp),
      0.06,
      label = paste0("causatr vs lmtp RD | sex=", s)
    )
  }
})


# Test D — bootstrap variance ---------------------------------------------
#
# Stratified ICE variance is bootstrap-only; the refit path must thread
# `stratified` through `fit_ice()` on every resample. Tier-1 asserts a
# finite, positive, well-ordered interval; Tier-2 checks coverage at n.
test_that("stratified ICE bootstrap variance is finite and well-formed", {
  d <- make_em_ice_binom_scm(n = 2500, seed = 11)
  fit_s <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    family = "binomial",
    stratified = "sex"
  )
  res <- contrast(
    fit_s,
    interventions = list(always = static(1), never = static(0)),
    type = "difference",
    reference = "never",
    ci_method = "bootstrap",
    n_boot = 50
  )
  se <- res$contrasts$se[1]
  expect_true(is.finite(se) && se > 0)
  expect_lt(res$contrasts$ci_lower[1], res$contrasts$ci_upper[1])
})

test_that("stratified ICE bootstrap CI covers the marginal truth", {
  skip_if_fast()
  d <- make_em_ice_binom_scm(n = 6000, seed = 3)
  fit_s <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    family = "binomial",
    stratified = "sex"
  )
  res <- contrast(
    fit_s,
    interventions = list(always = static(1), never = static(0)),
    type = "difference",
    reference = "never",
    ci_method = "bootstrap",
    n_boot = 300
  )
  expect_gt(res$contrasts$ci_upper[1], 0.58)
  expect_lt(res$contrasts$ci_lower[1], 0.58)
})


# Test E — composition with other ICE features ---------------------------

test_that("stratified ICE composes with a continuous-treatment shift", {
  d <- make_em_ice_cont_scm(n = 2000, seed = 9)
  fit_s <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    stratified = "sex"
  )
  pf <- ice_pf(fit_s, shift(1))
  pf0 <- ice_pf(fit_s, shift(0))
  expect_true(all(is.finite(pf)) && all(is.finite(pf0)))
  # shift(1) increases the (positive-effect) outcome on average.
  expect_gt(mean(pf), mean(pf0))
})

test_that("stratified ICE composes with a dynamic covariate rule", {
  d <- make_em_ice_scm(n = 2000, n_times = 2, seed = 5)
  fit_s <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    stratified = "sex"
  )
  rule <- dynamic(function(data, trt) as.integer(data$L > 0))
  pf <- ice_pf(fit_s, rule)
  expect_length(pf, 2000)
  expect_true(all(is.finite(pf)))
})

test_that("stratified ICE composes with external weights", {
  d <- make_em_ice_scm(n = 1500, n_times = 2, seed = 8)
  set.seed(8)
  w <- runif(nrow(d), 0.5, 2)
  fit_s <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    weights = w,
    stratified = "sex"
  )
  pf <- ice_pf(fit_s, static(1))
  expect_true(all(is.finite(pf)))
})

test_that("stratified ICE composes with a censoring row-filter and factor strata", {
  d <- make_em_ice_scm(n = 2000, n_times = 2, seed = 6)
  # Censor 8% of individuals at the final period: drop their outcome row
  # from fitting. C = 1 marks a censored person-period.
  d$C <- 0L
  set.seed(6)
  cens_ids <- sample(unique(d$id), size = 160)
  d$C[d$id %in% cens_ids & d$time == 1] <- 1L
  d$Y[d$id %in% cens_ids & d$time == 1] <- NA_real_
  # Factor-coded stratum exercises the non-integer split path.
  d$sexf <- factor(d$sex, labels = c("F", "M"))

  fit_s <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    censoring = "C",
    stratified = "sexf"
  )
  pf <- ice_pf(fit_s, static(1))
  expect_true(all(is.finite(pf)))
  ic <- ice_iterate(fit_s, static(1))
  expect_named(ic$models[[2]], c("F", "M"))
})


# Test F — rejection paths (classed errors) ------------------------------
#
# expect_error(class = ) rather than snapshots: several of these abort at
# causat() time with `call = call`, which embeds the (long, data-dependent)
# user call in the condition and makes message snapshots brittle.
test_that("stratified rejects non-ICE configurations", {
  dp <- data.frame(
    id = 1:200,
    G = rbinom(200, 1, 0.5),
    L = rnorm(200),
    A = rbinom(200, 1, 0.5),
    Y = rnorm(200)
  )
  # Point g-computation.
  expect_error(
    causat(dp, "Y", "A", ~G, estimator = "gcomp", stratified = "G"),
    class = "causatr_stratified_not_ice"
  )
  d <- make_em_ice_scm(n = 400, n_times = 2, seed = 1)
  # IPW (longitudinal, but not gcomp).
  expect_error(
    causat(
      d,
      "Y",
      "A",
      ~L0,
      confounders_tv = ~L,
      id = "id",
      time = "time",
      estimator = "ipw",
      stratified = "sex"
    ),
    class = "causatr_stratified_not_ice"
  )
})

test_that("stratified rejects a missing, time-varying, or continuous column", {
  d <- make_em_ice_scm(n = 400, n_times = 2, seed = 1)
  base_args <- list(
    data = d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    id = "id",
    time = "time"
  )
  expect_error(
    do.call(causat, c(base_args, list(stratified = "nope"))),
    class = "causatr_stratified_not_found"
  )
  # `L` is time-varying within individuals.
  expect_error(
    do.call(causat, c(base_args, list(stratified = "L"))),
    class = "causatr_stratified_not_baseline"
  )
  # `L0` is a continuous baseline covariate.
  expect_error(
    do.call(causat, c(base_args, list(stratified = "L0"))),
    class = "causatr_stratified_too_many"
  )
})

test_that("stratified ICE rejects the analytic sandwich", {
  d <- make_em_ice_scm(n = 600, n_times = 2, seed = 1)
  fit_s <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    stratified = "sex"
  )
  expect_error(
    contrast(
      fit_s,
      interventions = list(a1 = static(1), a0 = static(0)),
      ci_method = "sandwich"
    ),
    class = "causatr_stratified_ice_sandwich"
  )
})
