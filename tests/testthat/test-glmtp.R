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
})
