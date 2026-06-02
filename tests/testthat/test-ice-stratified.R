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
# ordered by first-period id). Computed directly so the point-estimate tests
# inspect the recursion output without going through a variance path.
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
    # Parametric ICE vs lmtp's semiparametric SDR: a two-sided agreement
    # bound on the per-stratum RD (|diff| < absolute tolerance set by the
    # SDR's cross-fit Monte-Carlo error).
    expect_lt(
      abs(rd_causatr - rd_lmtp),
      0.06,
      label = paste0("causatr vs lmtp RD | sex=", s)
    )
  }
})


# Test D — bootstrap variance ---------------------------------------------
#
# Stratified ICE variance has both a bootstrap and (Test G) a sandwich path;
# the bootstrap refit must thread `stratified` through `fit_ice()` on every
# resample. Both tests pin the bootstrap SE against an analytic oracle (the
# Test-G sandwich / delicatessen) with a two-sided tolerance, rather than
# checking only that the interval is finite and ordered.
test_that("stratified ICE bootstrap SE matches the analytic sandwich (binomial)", {
  skip_if_fast()
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
  ivs <- list(always = static(1), never = static(0))
  sw <- contrast(
    fit_s,
    interventions = ivs,
    type = "difference",
    reference = "never",
    ci_method = "sandwich"
  )
  bt <- contrast(
    fit_s,
    interventions = ivs,
    type = "difference",
    reference = "never",
    ci_method = "bootstrap",
    n_boot = 300
  )
  # Same point estimate exactly; SE agreement up to bootstrap MC noise.
  expect_equal(
    bt$contrasts$estimate[1],
    sw$contrasts$estimate[1],
    tolerance = 1e-8
  )
  expect_equal(bt$contrasts$se[1], sw$contrasts$se[1], tolerance = 0.15)
})

test_that("stratified ICE bootstrap SE matches delicatessen at n = 6000", {
  skip_if_fast()
  # Same DGP+seed as the binomial delicatessen fixture (Test G4): the
  # ID-cluster bootstrap SE must reproduce the M-estimation sandwich SE.
  ref_se_ate <- 0.014747
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
  expect_equal(res$contrasts$se[1], ref_se_ate, tolerance = 0.15)
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

# Test G — analytic sandwich variance --------------------------------------
#
# The per-stratum x per-time stacked-EE sandwich. Because G is baseline and
# individuals never cross strata, the system is block-diagonal across strata:
# Channel 2 (nuisance correction) runs once per stratum on disjoint rows, while
# Channel 1 (sampling) stays global. The two oracles below pin this down.

# G1 — exact reduction to the pooled sandwich at S = 1. A constant stratum
# column (one level) must give byte-identical SEs to the unstratified fit, a
# tripwire on the refactor and the stratum driver.
test_that("stratified ICE sandwich equals pooled sandwich at S = 1", {
  d <- make_em_ice_scm(n = 2000, n_times = 2, seed = 7)
  d$one <- 1L

  ivs <- list(a1 = static(1), a0 = static(0))
  fit_s <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    stratified = "one"
  )
  fit_p <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    id = "id",
    time = "time"
  )
  rs <- contrast(
    fit_s,
    interventions = ivs,
    type = "difference",
    ci_method = "sandwich"
  )
  rp <- contrast(
    fit_p,
    interventions = ivs,
    type = "difference",
    ci_method = "sandwich"
  )

  expect_equal(rs$contrasts$se[1], rp$contrasts$se[1], tolerance = 1e-8)
  expect_equal(
    rs$contrasts$estimate[1],
    rp$contrasts$estimate[1],
    tolerance = 1e-8
  )
})


# G2 — cross-package: delicatessen M-estimation sandwich (PRIMARY oracle).
# delicatessen's MEstimator computes the *same* plug-in M-estimation sandwich
# causatr does, so for a gaussian (identity-link) DGP the two coincide to
# numerical precision — not a loose semiparametric cross-check. Reference
# values from data-raw/ice_stratified_reference.py on the wide snapshot
# data-raw/ice_fixture_stratified_2t.csv, regenerated here from the SAME
# seeded DGP (make_em_ice_scm, seed = 2026) so both engines fit identical rows.
test_that("stratified ICE sandwich matches delicatessen (gaussian, 2 strata)", {
  # delicatessen ice_gcomp_2t_stratified(family="gaussian") reference
  # (deli_venv, v4.2). Both means, both per-arm SEs, the ATE and the ATE SE
  # are pinned — the entire stacked sandwich, not just one number.
  ref_mu_always <- 16.571339
  ref_se_always <- 0.048470
  ref_mu_never <- 9.958901
  ref_se_never <- 0.040256
  ref_ate <- 6.612438
  ref_se_ate <- 0.062145

  d <- make_em_ice_scm(n = 3000, n_times = 2, seed = 2026)
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
  res <- contrast(
    fit_s,
    interventions = list(always = static(1), never = static(0)),
    type = "difference",
    reference = "never",
    ci_method = "sandwich"
  )
  est <- res$estimates
  ia <- which(est$intervention == "always")
  ino <- which(est$intervention == "never")

  expect_equal(est$estimate[ia], ref_mu_always, tolerance = 1e-4)
  expect_equal(est$se[ia], ref_se_always, tolerance = 1e-4)
  expect_equal(est$estimate[ino], ref_mu_never, tolerance = 1e-4)
  expect_equal(est$se[ino], ref_se_never, tolerance = 1e-4)
  expect_equal(res$contrasts$estimate[1], ref_ate, tolerance = 1e-4)
  expect_equal(res$contrasts$se[1], ref_se_ate, tolerance = 1e-4)
})


# G3 — sandwich vs ID-cluster bootstrap parity on a genuine 2-stratum DGP.
test_that("stratified ICE sandwich agrees with ID-cluster bootstrap", {
  skip_if_fast()
  d <- make_em_ice_scm(n = 3000, n_times = 2, seed = 21)
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
  ivs <- list(always = static(1), never = static(0))
  rs <- contrast(
    fit_s,
    interventions = ivs,
    type = "difference",
    reference = "never",
    ci_method = "sandwich"
  )
  rb <- contrast(
    fit_s,
    interventions = ivs,
    type = "difference",
    reference = "never",
    ci_method = "bootstrap",
    n_boot = 200
  )
  expect_equal(rs$contrasts$se[1], rb$contrasts$se[1], tolerance = 0.2)
})


# G4 — cross-package: delicatessen sandwich on the non-collapsible BINOMIAL
# DGP (logit final / quasibinomial pseudo). Same tight M-estimation oracle as
# G2; reference from data-raw/ice_stratified_reference.py(family="binomial")
# on make_em_ice_binom_scm(n = 6000, seed = 3). The recovered marginal RD
# (~0.577) sits at the documented MC truth 0.58, and every per-arm mean / SE
# matches delicatessen, so the whole binomial stacked sandwich is pinned.
test_that("stratified ICE sandwich matches delicatessen (binomial, 2 strata)", {
  ref_mu_always <- 0.848290
  ref_se_always <- 0.007363
  ref_mu_never <- 0.271522
  ref_se_never <- 0.010392
  ref_ate <- 0.576768
  ref_se_ate <- 0.014747

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
    ci_method = "sandwich"
  )
  est <- res$estimates
  ia <- which(est$intervention == "always")
  ino <- which(est$intervention == "never")

  expect_equal(est$estimate[ia], ref_mu_always, tolerance = 1e-3)
  expect_equal(est$se[ia], ref_se_always, tolerance = 1e-3)
  expect_equal(est$estimate[ino], ref_mu_never, tolerance = 1e-3)
  expect_equal(est$se[ino], ref_se_never, tolerance = 1e-3)
  expect_equal(res$contrasts$estimate[1], ref_ate, tolerance = 1e-3)
  expect_equal(res$contrasts$se[1], ref_se_ate, tolerance = 1e-3)
})


# G5 — factor strata vs integer strata: relabelling 0/1 -> F/M cannot change
# the fit, so the per-arm sandwich SEs must be byte-identical. A tight equality
# oracle (1e-10) for the factor-coded split path -- not a finite/positive check.
test_that("stratified ICE sandwich is invariant to factor vs integer stratum coding", {
  d <- make_em_ice_scm(n = 2000, n_times = 2, seed = 6)
  d$sexf <- factor(d$sex, labels = c("F", "M"))
  ivs <- list(a1 = static(1), a0 = static(0))

  fit_i <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    stratified = "sex"
  )
  fit_f <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    stratified = "sexf"
  )
  ri <- contrast(
    fit_i,
    interventions = ivs,
    type = "difference",
    ci_method = "sandwich"
  )
  rf <- contrast(
    fit_f,
    interventions = ivs,
    type = "difference",
    ci_method = "sandwich"
  )

  expect_equal(rf$contrasts$se[1], ri$contrasts$se[1], tolerance = 1e-10)
  expect_equal(rf$estimates$se, ri$estimates$se, tolerance = 1e-10)
})


# G6 — sandwich vs ID-cluster bootstrap parity for continuous-treatment shift
# and external weights (the composition paths delicatessen does not cover here).
# A two-sided numerical comparison against the bootstrap variance, not a
# finite/positive check; tolerance set by the bootstrap's own Monte-Carlo noise
# at n_boot = 300.
test_that("stratified ICE sandwich matches bootstrap under shift and weights", {
  skip_if_fast()
  dc <- make_em_ice_cont_scm(n = 3000, seed = 9)
  fc <- causat(
    dc,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    stratified = "sex"
  )
  ivs_c <- list(s1 = shift(1), s0 = shift(0))
  rc_sw <- contrast(
    fc,
    interventions = ivs_c,
    type = "difference",
    ci_method = "sandwich"
  )
  rc_bt <- contrast(
    fc,
    interventions = ivs_c,
    type = "difference",
    ci_method = "bootstrap",
    n_boot = 300
  )
  expect_equal(rc_sw$contrasts$se[1], rc_bt$contrasts$se[1], tolerance = 0.15)

  dw <- make_em_ice_scm(n = 3000, n_times = 2, seed = 8)
  set.seed(8)
  w <- runif(nrow(dw), 0.5, 2)
  fw <- causat(
    dw,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    weights = w,
    stratified = "sex"
  )
  ivs_w <- list(a1 = static(1), a0 = static(0))
  rw_sw <- contrast(
    fw,
    interventions = ivs_w,
    type = "difference",
    ci_method = "sandwich"
  )
  rw_bt <- contrast(
    fw,
    interventions = ivs_w,
    type = "difference",
    ci_method = "bootstrap",
    n_boot = 300
  )
  expect_equal(rw_sw$contrasts$se[1], rw_bt$contrasts$se[1], tolerance = 0.15)
})


# G7 — stochastic intervention through the per-stratum sandwich driver.
#
# Critical-review round 2026-06-02: the `is_stoch` branch of
# variance_if_ice_chain() is the one path no other Test-G case exercises
# through the stratified driver. Repro /tmp/causatr_repro_stoch_parity.R
# confirmed the per-draw frame and the deterministic frame share the same
# row set (intervention application is row-preserving), so the chain's
# Monte-Carlo gradient is well-defined even with censoring (valid_target not
# all TRUE). Pinned here by sandwich-vs-ID-cluster-bootstrap parity, a
# two-sided numerical comparison.
test_that("stratified ICE sandwich matches bootstrap for a stochastic intervention", {
  skip_if_fast()
  d <- make_em_ice_scm(n = 2500, n_times = 2, seed = 11)
  # Censor 10% at the final period so target rows are absent at fit time and
  # the chain's `valid_target` / `vt_m` masks are non-trivial.
  d$C <- 0L
  set.seed(11)
  cens <- sample(unique(d$id), 250)
  d$C[d$id %in% cens & d$time == 1] <- 1L
  d$Y[d$id %in% cens & d$time == 1] <- NA_real_

  samp <- function(data, treatment) stats::rbinom(nrow(data), 1, 0.5)
  ivs <- list(s = stochastic(samp, n_mc = 30L), n = static(0))
  fit_s <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    censoring = "C",
    stratified = "sex"
  )
  sw <- contrast(
    fit_s,
    interventions = ivs,
    type = "difference",
    reference = "n",
    ci_method = "sandwich"
  )
  bt <- contrast(
    fit_s,
    interventions = ivs,
    type = "difference",
    reference = "n",
    ci_method = "bootstrap",
    n_boot = 200
  )
  expect_equal(sw$contrasts$se[1], bt$contrasts$se[1], tolerance = 0.2)
})
