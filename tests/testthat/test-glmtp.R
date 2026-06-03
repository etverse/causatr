# Tests for natural-history modified treatment policies (G-LMTPs), the
# augmented-data sequential regression of Diaz, Williams, Morzywolek & Rudolph
# (2026, arXiv:2605.24167). The estimator is genuinely irreducible for policies
# that depend on the natural-value history of treatment: the standard ICE
# recursion would condition on the OBSERVED lag rather than the counterfactual
# natural value, returning a quietly-wrong number. The truth oracle is therefore
# forward Monte-Carlo simulation of the natural-history regime (lmtp's one-shot
# shift computes the standard LMTP, the wrong estimand here), backed by exact
# limiting-case equivalences against the existing engine.

# ---------------------------------------------------------------------------
# Chunk 1 -- augmented-frame infrastructure (R/glmtp_augment.R)
# ---------------------------------------------------------------------------

test_that("glmtp_support returns the sorted discrete support and drops censored rows", {
  dt <- data.table::data.table(
    A = c(0, 1, 1, 0, 2, 1),
    C = c(0, 0, 0, 0, 1, 0)
  )
  # Binary support, no censoring column.
  expect_identical(glmtp_support(dt[A %in% c(0, 1)], "A"), c(0, 1))
  # Ordinal support {0,1,2}; the censored row (C == 1, A == 2) is excluded,
  # so 2 drops out of the support entirely.
  expect_identical(glmtp_support(dt, "A", censoring = "C"), c(0, 1))
  # Without the censoring filter the full ordinal support is returned.
  expect_identical(glmtp_support(dt, "A"), c(0, 1, 2))
})

test_that("glmtp_support rejects continuous, factor, and multivariate treatment", {
  cont <- data.table::data.table(A = c(0.1, 1.7, 2.3))
  expect_error(glmtp_support(cont, "A"), class = "causatr_glmtp_continuous_trt")

  fac <- data.table::data.table(A = factor(c("lo", "hi", "lo")))
  expect_error(glmtp_support(fac, "A"), class = "causatr_glmtp_continuous_trt")

  mv <- data.table::data.table(A1 = c(0, 1), A2 = c(1, 0))
  expect_error(
    glmtp_support(mv, c("A1", "A2")),
    class = "causatr_glmtp_continuous_trt"
  )
})

test_that("glmtp_enumerate_labels builds the product support with the empty base", {
  support <- c(0, 1)
  # t = 0 -> a single empty label (the base of the recursion for q_1).
  l0 <- glmtp_enumerate_labels(support, 0L)
  expect_length(l0, 1L)
  expect_length(l0[[1L]], 0L)

  # t = 1 -> |A| labels.
  l1 <- glmtp_enumerate_labels(support, 1L)
  expect_length(l1, 2L)
  expect_setequal(vapply(l1, identity, numeric(1)), support)

  # t = 2 -> |A|^2 labels, all unique, first coordinate varies slowest.
  l2 <- glmtp_enumerate_labels(support, 2L)
  expect_length(l2, 4L)
  keys <- vapply(l2, glmtp_label_key, character(1))
  expect_length(unique(keys), 4L)

  # Ordinal support: cardinality is |A|^t.
  expect_length(glmtp_enumerate_labels(c(0, 1, 2), 3L), 27L)
})

test_that("glmtp_label_key is collision-free and maps the empty sequence to ''", {
  expect_identical(glmtp_label_key(numeric(0)), "")
  expect_identical(glmtp_label_key(c(0, 1, 0)), "0|1|0")
  # Distinct sequences that would collide under naive concatenation get
  # distinct keys (1,11) vs (11,1).
  expect_false(
    identical(glmtp_label_key(c(1, 11)), glmtp_label_key(c(11, 1)))
  )
})

test_that("glmtp_check_tractable caps the worst-step blow-up", {
  # Binary, 5 periods: worst step enumerates 2^4 = 16 labels -- within budget.
  expect_identical(glmtp_check_tractable(c(0, 1), 5L, budget = 1024L), 16)
  # A single period has no history to enumerate -> 1 label.
  expect_identical(glmtp_check_tractable(c(0, 1), 1L), 1)
  # Binary, 12 periods: 2^11 = 2048 > 1024 budget -> abort.
  expect_error(
    glmtp_check_tractable(c(0, 1), 12L, budget = 1024L),
    class = "causatr_glmtp_too_many"
  )
  # Raising the budget admits the same problem.
  expect_identical(
    glmtp_check_tractable(c(0, 1), 12L, budget = 4096L),
    2048
  )
})

# ---------------------------------------------------------------------------
# Chunk 2/3 -- the augmented engine: constructors, truth, limiting cases
# ---------------------------------------------------------------------------

# Fit a standard longitudinal-ICE object on a glmtp DGP. The augmented engine
# reuses fit_ice metadata, so this is the same fit for both ice_iterate() and
# glmtp_iterate().
glmtp_fit <- function(data, family = "gaussian") {
  causat(
    data,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    id = "id",
    time = "t",
    estimator = "gcomp",
    family = family,
    history = Inf
  )
}

# Point estimate for the ATE target (mean baseline pseudo-outcome), bypassing
# the bootstrap so the truth checks are deterministic and fast.
glmtp_point <- function(fit, iv) {
  mean(glmtp_iterate(fit, iv)$pseudo_final)
}

test_that("grace_period / carry_forward validate args and carry the policy", {
  expect_error(grace_period(-1L), class = "causatr_bad_intervention_arg")
  expect_error(grace_period(1.5), class = "causatr_bad_intervention_arg")
  expect_error(
    grace_period(1L, budget = 0L),
    class = "causatr_bad_intervention_arg"
  )
  expect_error(
    carry_forward(seed = "later"),
    class = "causatr_bad_intervention_arg"
  )
  expect_error(
    carry_forward(seed = c(1, 2)),
    class = "causatr_bad_intervention_arg"
  )

  g <- grace_period(2L)
  expect_s3_class(g, "causatr_glmtp")
  expect_s3_class(g, "causatr_intervention")
  expect_identical(g$window, 2L)
  expect_true(is.function(g$policy))

  cf <- carry_forward()
  expect_identical(cf$seed, "baseline")
  expect_true(any_glmtp(list(a = g, b = NULL)))
  expect_false(any_glmtp(list(a = static(1), b = NULL)))
  expect_output(print(g), "grace_period")
})

# Truth oracle: forward Monte-Carlo of the natural-history regime (the paper's
# Proposition 1). The constants below were computed at n_mc = 2e6, seed = 1 via
# the helper-glmtp-dgp.R forward-simulators and are reconfirmed in a Tier-2 test
# so they cannot silently drift.
TRUTH_GAUSS_W1 <- 2.30313
TRUTH_GAUSS_NAT <- 2.98045
TRUTH_PAPER_W1 <- 0.55472

test_that("glmtp grace_period recovers the forward-MC truth (gaussian)", {
  # Linear-gaussian absorbing-treatment SCM: every augmented per-step
  # regression is exactly linear, so the additive-GLM plug-in is consistent.
  d <- glmtp_delay_data(n = 12000L, seed = 1L)$data
  fit <- glmtp_fit(d)
  est <- glmtp_point(fit, grace_period(1L))
  # Two-sided agreement with the Proposition-1 truth (gap ~1 SE at this n).
  expect_lt(abs(est - TRUTH_GAUSS_W1), 0.025)
})

test_that("the augmented engine is necessary: naive dynamic(lag1_A) is far off", {
  # The standard ICE recursion conditions on the OBSERVED lag, so a dynamic
  # rule that reads `lag1_A` targets a different (wrong) estimand. The gap to
  # truth is enormous (~1.3 on a value of 2.3) while glmtp is within ~0.02.
  d <- glmtp_delay_data(n = 12000L, seed = 1L)$data
  fit <- glmtp_fit(d)
  mu_glmtp <- glmtp_point(fit, grace_period(1L))
  naive_rule <- dynamic(function(data, trt) {
    v <- data$lag1_A
    v[is.na(v)] <- 0L
    v
  })
  mu_naive <- mean(ice_iterate(fit, naive_rule)$pseudo_final)
  expect_lt(abs(mu_glmtp - TRUTH_GAUSS_W1), 0.03)
  expect_gt(abs(mu_naive - TRUTH_GAUSS_W1), 0.5)
})

test_that("carry_forward equals the equivalent standard-ICE baseline regime", {
  # carry_forward(baseline) sets every period's treatment to the baseline
  # natural value A_1, which is a baseline-observed quantity -- so it must equal
  # a standard ICE dynamic rule that sets A_t to the baseline-A column, exactly.
  dat <- glmtp_delay_data(n = 3000L, seed = 7L)$data
  base_A <- dat$A[dat$t == 1L]
  dat$baseA <- rep(base_A, each = max(dat$t))
  fit_glmtp_obj <- glmtp_fit(dat)
  fit_dyn <- causat(
    dat,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L0 + baseA,
    confounders_tv = ~L,
    id = "id",
    time = "t",
    estimator = "gcomp",
    history = Inf
  )
  mu_cf <- glmtp_point(fit_glmtp_obj, carry_forward())
  mu_dyn <- mean(
    ice_iterate(fit_dyn, dynamic(function(d, trt) d$baseA))$pseudo_final
  )
  expect_equal(mu_cf, mu_dyn, tolerance = 1e-8)
})

test_that("grace_period(0) reduces to the natural course exactly (absorbing DGP)", {
  # With window 0 the delay policy is "absorbing initiation at the natural
  # time", which in an absorbing-treatment DGP is the observed treatment. The
  # augmented identity-vs-grace(0) means coincide to numerical precision.
  d <- glmtp_delay_data(n = 4000L, seed = 3L)$data
  fit <- glmtp_fit(d)
  mu_g0 <- glmtp_point(fit, grace_period(0L))
  mu_nat <- glmtp_point(fit, NULL)
  expect_equal(mu_g0, mu_nat, tolerance = 1e-9)
})

test_that("glmtp replicates the Diaz et al. (2026) Section-6 delay result", {
  # Paper mechanism (binary absorbing treatment, propensity logit(-1.5 + 0.3 L),
  # L_t = 0.5 L_{t-1} + eps, tau = 5) with a single end-of-study binary outcome.
  # glmtp's augmented plug-in recovers the Proposition-1 forward-MC truth for
  # the one-period-delay policy.
  d <- glmtp_paper_data(n = 8000L, seed = 1L)$data
  fit <- glmtp_fit(d, family = "binomial")
  est <- glmtp_point(fit, grace_period(1L))
  expect_lt(abs(est - TRUTH_PAPER_W1), 0.01)
})

test_that("bootstrap CI covers the forward-MC truth (gaussian)", {
  d <- glmtp_delay_data(n = 1500L, seed = 11L)$data
  fit <- glmtp_fit(d)
  res <- contrast(
    fit,
    interventions = list(delay1 = grace_period(1L)),
    ci_method = "bootstrap",
    n_boot = 200L
  )
  row <- res$estimates[res$estimates$intervention == "delay1", ]
  expect_true(is.finite(row$se) && row$se > 0)
  expect_gte(TRUTH_GAUSS_W1, row$ci_lower)
  expect_lte(TRUTH_GAUSS_W1, row$ci_upper)
})

test_that("glmtp composes with a censoring row-filter and external weights", {
  d <- glmtp_delay_data(n = 2500L, seed = 5L)$data
  # Light random censoring at the final period (row filter), plus survey weights.
  set.seed(99)
  d$cens <- 0L
  tau <- max(d$t)
  fin <- d$t == tau
  d$cens[fin] <- stats::rbinom(sum(fin), 1L, 0.1)
  w <- runif(nrow(d), 0.5, 1.5)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    id = "id",
    time = "t",
    estimator = "gcomp",
    censoring = "cens",
    weights = w,
    history = Inf
  )
  res <- contrast(
    fit,
    interventions = list(g = grace_period(1L), nat = NULL),
    ci_method = "bootstrap",
    n_boot = 40L
  )
  expect_true(all(is.finite(res$estimates$estimate)))
  expect_true(all(is.finite(res$estimates$se) & res$estimates$se > 0))
})

# Ordinal (k > 2) treatment exercises the |A|^t label machinery, which the
# binary tests above never reach. carry_forward() reduces to a baseline regime
# for ANY discrete treatment, so it must equal the equivalent standard-ICE
# dynamic rule exactly on a 3-level treatment.
# (2026-06-03 critical review Issue #1: ordinal support was verified to 1e-14
# but untested; /tmp/causatr_repro_glmtp_adversarial.R.)
test_that("carry_forward equals the baseline regime for an ordinal (3-level) treatment", {
  set.seed(20)
  n <- 2500L
  tau <- 3L
  L0 <- stats::rnorm(n)
  A <- matrix(0L, n, tau)
  L <- matrix(0, n, tau)
  Lprev <- L0
  Aprev <- integer(n)
  for (t in seq_len(tau)) {
    Lt <- if (t == 1L) {
      0.5 * L0 + stats::rnorm(n)
    } else {
      0.5 * Lprev + 0.3 * Aprev + stats::rnorm(n)
    }
    lam <- stats::plogis(-0.2 + 0.4 * Lt + 0.5 * Aprev) * 1.2
    At <- pmin(2L, pmax(0L, stats::rpois(n, lambda = lam)))
    A[, t] <- At
    L[, t] <- Lt
    Lprev <- Lt
    Aprev <- At
  }
  # Confirm the treatment is genuinely 3-level (the |A|^t labels are exercised).
  expect_setequal(sort(unique(as.vector(A))), c(0L, 1L, 2L))
  Y <- 1 + 0.5 * rowSums(A) + 0.4 * L[, tau] - 0.3 * L0 + stats::rnorm(n)
  d <- data.frame(
    id = rep(seq_len(n), each = tau),
    t = rep(seq_len(tau), n),
    L0 = rep(L0, each = tau),
    A = as.vector(t(A)),
    L = as.vector(t(L)),
    Y = NA_real_
  )
  d$Y[d$t == tau] <- Y
  d$baseA <- rep(d$A[d$t == 1L], each = tau)
  fit <- glmtp_fit(d)
  fit2 <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L0 + baseA,
    confounders_tv = ~L,
    id = "id",
    time = "t",
    estimator = "gcomp",
    history = Inf
  )
  mu_cf <- glmtp_point(fit, carry_forward())
  mu_dyn <- mean(
    ice_iterate(fit2, dynamic(function(x, trt) x$baseA))$pseudo_final
  )
  expect_equal(mu_cf, mu_dyn, tolerance = 1e-8)
})

# cap_escalation() is internal/deferred (PHASE_22): the parametric plug-in is
# consistent only under a flexible dose model causatr's ICE does not expose
# (a factor-dose model recovers the forward-MC truth to ~1e-3; the bare-numeric
# model is ~3% biased when the cap binds). The ENGINE, however, must be bug-free
# on the ordinal natural-value-dependent path. A from-scratch hand-coded
# recursion (different code structure, bare-numeric glm) is the independent
# oracle: the two agree to numerical precision, proving the recursion -- not the
# method -- is correct, before any "plug-in bias" attribution.
# (2026-06-03 critical review; /tmp/causatr_glmtp_bugcheck.R.)
test_that("the engine matches an independent hand-coded recursion for an ordinal dose-cap", {
  d <- glmtp_ordinal_cap_data(n = 6000L, seed = 1L)$data
  fit <- glmtp_fit(d)
  for (delta in c(1, 10)) {
    mu_engine <- mean(glmtp_iterate(fit, cap_escalation(delta))$pseudo_final)
    mu_hand <- glmtp_handcode_cap(fit, delta)
    expect_equal(mu_engine, mu_hand, tolerance = 1e-9)
  }
})

# Flexible-treatment ICE term (Phase 22b-5): the bare-numeric per-step model
# misspecifies the kinked capped-dose pseudo-response (~0.8-1.1% asymptotic bias
# when the cap binds), whereas `treatment_form = ~ factor(A)` lets each dose
# level carry its own coefficient and recovers the forward-MC truth to sampling
# error. The constant below is the Proposition-1 forward-MC truth at n_mc = 2e6,
# seed = 1 (reconfirmed in the Tier-2 test so it cannot silently drift).
TRUTH_CAP_W1 <- 4.0966

test_that("factor(A) treatment term recovers the cap forward-MC truth; bare-numeric is biased", {
  # n = 40000 keeps the test affordable on CI. The headline claim is the
  # seed-robust *comparison*: factor(A) is far closer to the forward-MC truth
  # than the bare-numeric plug-in. Across seeds 1-4 at this n the factor gap is
  # 0.0015-0.0101 and the bare gap is 0.028-0.044 (ratio 3.1-22x), so the bounds
  # below hold independent of seed. (factor(A) is consistent -- the gap settles
  # to ~0.001-0.0025 by n = 80000 -- but a tight absolute band at n = 40000 is
  # seed-fragile, hence the comparative assertion.)
  d <- glmtp_ordinal_cap_data(n = 40000L, seed = 1L)$data
  fit_bare <- glmtp_fit(d)
  fit_factor <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    id = "id",
    time = "t",
    estimator = "gcomp",
    history = Inf,
    treatment_form = ~ factor(A)
  )
  gap_bare <- abs(mean(glmtp_iterate(fit_bare, cap_escalation(1))$pseudo_final) - TRUTH_CAP_W1)
  gap_factor <- abs(mean(glmtp_iterate(fit_factor, cap_escalation(1))$pseudo_final) - TRUTH_CAP_W1)
  # factor(A) lands close to the truth; the bare-numeric plug-in keeps a
  # seed-stable ~0.03-0.04 cap-kink bias; factor(A) is at least ~2x closer.
  expect_lt(gap_factor, 0.015)
  expect_gt(gap_bare, 0.02)
  expect_lt(gap_factor, gap_bare / 2)
})

test_that("factor(A) composes with grace_period() end-to-end (bootstrap CI)", {
  d <- glmtp_delay_data(n = 1500L, seed = 13L)$data
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    id = "id",
    time = "t",
    estimator = "gcomp",
    history = Inf,
    treatment_form = ~ factor(A)
  )
  res <- contrast(
    fit,
    interventions = list(delay1 = grace_period(1L), nat = NULL),
    ci_method = "bootstrap",
    n_boot = 40L
  )
  expect_true(all(is.finite(res$estimates$estimate)))
  expect_true(all(is.finite(res$estimates$se) & res$estimates$se > 0))
  # Binary support: factor(A) and bare numeric span the same two-point design,
  # so the grace-period point estimate is unchanged (sanity, not a new claim).
  mu_factor <- mean(glmtp_iterate(fit, grace_period(1L))$pseudo_final)
  mu_bare <- mean(glmtp_iterate(glmtp_fit(d), grace_period(1L))$pseudo_final)
  expect_equal(mu_factor, mu_bare, tolerance = 1e-8)
})

test_that("cap_escalation with delta >= max increase reduces to the natural course", {
  # When the cap never binds the augmented recursion reproduces the observed
  # treatment exactly -- a tight check that the policy/prediction plumbing is
  # correct independent of any model-misspecification bias.
  d <- glmtp_ordinal_cap_data(n = 4000L, seed = 2L)$data
  fit <- glmtp_fit(d)
  mu_cap_big <- mean(glmtp_iterate(fit, cap_escalation(10))$pseudo_final)
  mu_nat <- mean(glmtp_iterate(fit, NULL)$pseudo_final)
  expect_equal(mu_cap_big, mu_nat, tolerance = 1e-10)
})

test_that("a per-label fit failure surfaces a classed warning, not a silent estimate", {
  # 2026-06-03 critical review Issue #3: a per-label model failure degrades the
  # affected rows to NA and drops them from later steps. Force one failure with
  # a flaky model_fn (errors on its 3rd call -- the base fit is call 1, so this
  # hits a per-label pseudo-outcome fit) and confirm the classed warning fires.
  d <- glmtp_delay_data(n = 800L, seed = 8L)$data
  call_n <- 0L
  # The base outcome model is fit first (call 1); per-label pseudo-outcome fits
  # follow. Failing call 3 forces exactly one per-label failure. Extra args
  # (incl. `weights` when present) flow through `...` as values, avoiding glm's
  # non-standard evaluation of a re-passed `weights` symbol.
  flaky <- function(formula, data, family, ...) {
    call_n <<- call_n + 1L
    if (call_n == 3L) {
      rlang::abort("forced per-label fit failure")
    }
    stats::glm(formula, data = data, family = family, ...)
  }
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    id = "id",
    time = "t",
    estimator = "gcomp",
    model_fn = flaky,
    history = Inf
  )
  expect_warning(
    glmtp_iterate(fit, grace_period(1L)),
    class = "causatr_glmtp_label_fit_failed"
  )
})

# ---------------------------------------------------------------------------
# Rejections at the contrast() boundary
# ---------------------------------------------------------------------------

test_that("glmtp is rejected outside longitudinal g-computation", {
  # Point treatment.
  dpt <- data.frame(Y = rnorm(80), A = rbinom(80, 1, 0.5), L = rnorm(80))
  fit_pt <- causat(dpt, outcome = "Y", treatment = "A", confounders = ~L)
  expect_error(
    contrast(fit_pt, list(g = grace_period(1L))),
    class = "causatr_glmtp_not_ice"
  )
  # Longitudinal IPW.
  d <- glmtp_delay_data(n = 800L, seed = 2L)$data
  fit_ipw <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    id = "id",
    time = "t",
    estimator = "ipw"
  )
  expect_error(
    contrast(fit_ipw, list(g = grace_period(1L))),
    class = "causatr_glmtp_not_ice"
  )
})

test_that("glmtp rejects mixing, sandwich, and a continuous treatment", {
  d <- glmtp_delay_data(n = 800L, seed = 4L)$data
  fit <- glmtp_fit(d)
  expect_error(
    contrast(fit, list(g = grace_period(1L), s = static(1))),
    class = "causatr_glmtp_mixed"
  )
  expect_error(
    contrast(fit, list(g = grace_period(1L)), ci_method = "sandwich"),
    class = "causatr_glmtp_sandwich"
  )
  # Continuous treatment: make A non-integer so the discreteness gate fires.
  dc <- d
  dc$A <- dc$A + stats::runif(nrow(dc), 0, 0.3)
  fit_c <- causat(
    dc,
    outcome = "Y",
    treatment = "A",
    confounders = ~L0,
    confounders_tv = ~L,
    id = "id",
    time = "t",
    estimator = "gcomp",
    history = Inf
  )
  expect_error(
    contrast(fit_c, list(g = grace_period(1L))),
    class = "causatr_glmtp_continuous_trt"
  )
})

# ---------------------------------------------------------------------------
# Tier-2: confirm the hard-coded forward-MC truth constants do not drift
# ---------------------------------------------------------------------------

test_that("hard-coded forward-MC truth constants match a fresh 2e6 simulation", {
  skip_on_cran()
  expect_equal(
    glmtp_delay_forward_truth(1L, n_mc = 2e6, seed = 1L),
    TRUTH_GAUSS_W1,
    tolerance = 1e-3
  )
  expect_equal(
    glmtp_paper_forward_truth(1L, n_mc = 2e6, seed = 1L),
    TRUTH_PAPER_W1,
    tolerance = 1e-3
  )
  expect_equal(
    glmtp_cap_forward_truth(1, n_mc = 2e6, seed = 1L),
    TRUTH_CAP_W1,
    tolerance = 1e-3
  )
})
