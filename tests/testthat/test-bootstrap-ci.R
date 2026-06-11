# User-selectable bootstrap CI flavour (`boot_ci = "percentile" / "normal"`).
#
# Both intervals come from the same stored replicate matrix `boot_t`, so the
# tests recompute each definition directly from `boot_t` (the truth) and assert
# the stored CIs / accessors match exactly. `boot::boot.ci()` is an external
# cross-check for the percentile interval (its order-statistic interpolation
# differs slightly from quantile type-7, so the tolerance is loosened at large
# R). The point estimate, SE, and vcov must be invariant to the flavour.

# Small binomial DGP: probabilities in (0, 1) make percentile and normal
# intervals genuinely differ (the normal interval is symmetric, the percentile
# interval is not).
sim_boot_ci_dgp <- function(n = 1500, seed = 42) {
  set.seed(seed)
  L <- stats::rnorm(n)
  A <- stats::rbinom(n, 1L, stats::plogis(0.3 * L))
  Y <- stats::rbinom(n, 1L, stats::plogis(-0.3 + 0.8 * A + 0.5 * L))
  data.frame(Y = Y, A = A, L = L)
}

z975 <- stats::qnorm(0.975)

test_that("percentile bootstrap CIs equal the empirical replicate quantiles", {
  d <- sim_boot_ci_dgp()
  fit <- causat(
    d,
    "Y",
    "A",
    confounders = ~L,
    estimator = "gcomp",
    family = "binomial"
  )
  ivs <- list(a1 = static(1), a0 = static(0))
  set.seed(1)
  res <- contrast(
    fit,
    ivs,
    reference = "a0",
    ci_method = "bootstrap",
    n_boot = 500,
    boot_ci = "percentile"
  )

  expect_equal(res$boot_ci, "percentile")
  bt <- res$boot_t
  est <- as.data.frame(res$estimates)
  con <- as.data.frame(res$contrasts)

  # Means: stored CI == quantiles of each intervention's replicate column.
  for (nm in c("a1", "a0")) {
    q <- stats::quantile(bt[, nm], c(0.025, 0.975), names = FALSE)
    row <- est$intervention == nm
    expect_equal(c(est$ci_lower[row], est$ci_upper[row]), q, tolerance = 1e-12)
  }
  # Difference contrast: stored CI == quantiles of the per-replicate diff.
  diffs <- bt[, "a1"] - bt[, "a0"]
  q_d <- stats::quantile(diffs, c(0.025, 0.975), names = FALSE)
  expect_equal(c(con$ci_lower, con$ci_upper), q_d, tolerance = 1e-12)
})

test_that("normal bootstrap CIs equal the Wald interval from the bootstrap SE", {
  d <- sim_boot_ci_dgp()
  fit <- causat(
    d,
    "Y",
    "A",
    confounders = ~L,
    estimator = "gcomp",
    family = "binomial"
  )
  ivs <- list(a1 = static(1), a0 = static(0))
  set.seed(1)
  rn <- contrast(
    fit,
    ivs,
    reference = "a0",
    ci_method = "bootstrap",
    n_boot = 500,
    boot_ci = "normal"
  )
  set.seed(1)
  rp <- contrast(
    fit,
    ivs,
    reference = "a0",
    ci_method = "bootstrap",
    n_boot = 500,
    boot_ci = "percentile"
  )

  est <- as.data.frame(rn$estimates)
  bt <- rn$boot_t
  # Means: estimate +/- z * sd(replicates), and equals estimate +/- z * se.
  for (nm in c("a1", "a0")) {
    row <- est$intervention == nm
    wald <- est$estimate[row] + c(-1, 1) * z975 * stats::sd(bt[, nm])
    expect_equal(
      c(est$ci_lower[row], est$ci_upper[row]),
      wald,
      tolerance = 1e-12
    )
    se_wald <- est$estimate[row] + c(-1, 1) * z975 * est$se[row]
    expect_equal(
      c(est$ci_lower[row], est$ci_upper[row]),
      se_wald,
      tolerance = 1e-12
    )
  }
  # Difference contrast normal == estimate +/- z * se (the delta interval).
  con <- as.data.frame(rn$contrasts)
  delta <- con$estimate + c(-1, 1) * z975 * con$se
  expect_equal(c(con$ci_lower, con$ci_upper), delta, tolerance = 1e-12)

  # Point estimate, SE, and vcov are invariant to the CI flavour.
  expect_equal(as.data.frame(rp$estimates)$estimate, est$estimate)
  expect_equal(as.data.frame(rp$estimates)$se, est$se)
  expect_equal(rp$vcov, rn$vcov)
})

test_that("ratio / OR percentile CIs are the replicate quantiles on each scale", {
  d <- sim_boot_ci_dgp()
  fit <- causat(
    d,
    "Y",
    "A",
    confounders = ~L,
    estimator = "gcomp",
    family = "binomial"
  )
  ivs <- list(a1 = static(1), a0 = static(0))

  for (ty in c("ratio", "or")) {
    set.seed(5)
    res <- contrast(
      fit,
      ivs,
      reference = "a0",
      type = ty,
      ci_method = "bootstrap",
      n_boot = 500,
      boot_ci = "percentile"
    )
    bt <- res$boot_t
    ra <- bt[, "a1"]
    rr <- bt[, "a0"]
    reps_c <- if (ty == "ratio") {
      ra / rr
    } else {
      (ra / (1 - ra)) / (rr / (1 - rr))
    }
    q <- stats::quantile(
      reps_c[is.finite(reps_c)],
      c(0.025, 0.975),
      names = FALSE
    )
    con <- as.data.frame(res$contrasts)
    expect_equal(c(con$ci_lower, con$ci_upper), q, tolerance = 1e-12)
    # Percentile ratio/OR bounds are strictly positive (support-respecting).
    expect_true(con$ci_lower > 0)
  }
})

test_that("percentile CIs cross-check against boot::boot.ci", {
  skip_if_not_installed("boot")
  d <- sim_boot_ci_dgp(n = 2000, seed = 3)
  fit <- causat(
    d,
    "Y",
    "A",
    confounders = ~L,
    estimator = "gcomp",
    family = "binomial"
  )
  ivs <- list(a1 = static(1), a0 = static(0))
  set.seed(2)
  res <- contrast(
    fit,
    ivs,
    reference = "a0",
    ci_method = "bootstrap",
    n_boot = 2000,
    boot_ci = "percentile"
  )
  est <- as.data.frame(res$estimates)
  bt <- res$boot_t
  for (nm in c("a1", "a0")) {
    row <- est$intervention == nm
    bo <- structure(
      list(
        t0 = est$estimate[row],
        t = matrix(bt[, nm], ncol = 1),
        R = nrow(bt),
        sim = "ordinary",
        stype = "i",
        call = quote(boot())
      ),
      class = "boot"
    )
    bci <- boot::boot.ci(bo, type = "perc", conf = 0.95)$percent[4:5]
    # boot.ci's order-statistic interpolation differs slightly from
    # quantile() type 7; at R = 2000 the gap is well under 0.005.
    expect_equal(
      c(est$ci_lower[row], est$ci_upper[row]),
      bci,
      tolerance = 0.005
    )
  }
})

test_that("confint() and tidy() honour the stored boot_ci convention", {
  d <- sim_boot_ci_dgp()
  fit <- causat(
    d,
    "Y",
    "A",
    confounders = ~L,
    estimator = "gcomp",
    family = "binomial"
  )
  ivs <- list(a1 = static(1), a0 = static(0))
  set.seed(1)
  rp <- contrast(
    fit,
    ivs,
    reference = "a0",
    ci_method = "bootstrap",
    n_boot = 500,
    boot_ci = "percentile"
  )
  ep <- as.data.frame(rp$estimates)

  # confint default == stored percentile CIs.
  ci <- confint(rp)
  expect_equal(
    unname(ci["a1", ]),
    c(
      ep$ci_lower[ep$intervention == "a1"],
      ep$ci_upper[ep$intervention == "a1"]
    ),
    tolerance = 1e-12
  )
  # confint(boot_ci = "normal") override == Wald from the bootstrap SE.
  ci_n <- confint(rp, boot_ci = "normal")
  wald_a1 <- ep$estimate[ep$intervention == "a1"] +
    c(-1, 1) * z975 * ep$se[ep$intervention == "a1"]
  expect_equal(unname(ci_n["a1", ]), wald_a1, tolerance = 1e-12)

  # tidy() reads the stored (percentile) CIs for a bootstrap result.
  td <- tidy(rp, which = "all")
  td_m <- td[td$type == "mean", ]
  expect_equal(td_m$conf.low, ep$ci_lower, tolerance = 1e-12)
  expect_equal(td_m$conf.high, ep$ci_upper, tolerance = 1e-12)
})

test_that("boot_ci defaults to percentile and is recorded on the result", {
  d <- sim_boot_ci_dgp()
  fit <- causat(
    d,
    "Y",
    "A",
    confounders = ~L,
    estimator = "gcomp",
    family = "binomial"
  )
  ivs <- list(a1 = static(1), a0 = static(0))
  set.seed(1)
  res_default <- contrast(
    fit,
    ivs,
    reference = "a0",
    ci_method = "bootstrap",
    n_boot = 300
  )
  set.seed(1)
  res_perc <- contrast(
    fit,
    ivs,
    reference = "a0",
    ci_method = "bootstrap",
    n_boot = 300,
    boot_ci = "percentile"
  )
  expect_equal(res_default$boot_ci, "percentile")
  # Default == explicit percentile.
  expect_equal(
    as.data.frame(res_default$estimates)$ci_lower,
    as.data.frame(res_perc$estimates)$ci_lower
  )
})

test_that("boot_ci composes with the IPW and AIPW bootstrap paths", {
  d <- sim_boot_ci_dgp()
  ivs <- list(a1 = static(1), a0 = static(0))
  for (est_kind in c("ipw", "aipw")) {
    fit <- causat(
      d,
      "Y",
      "A",
      confounders = ~L,
      estimator = est_kind,
      family = "binomial"
    )
    set.seed(4)
    res <- contrast(
      fit,
      ivs,
      reference = "a0",
      ci_method = "bootstrap",
      n_boot = 400,
      boot_ci = "percentile"
    )
    bt <- res$boot_t
    est <- as.data.frame(res$estimates)
    for (nm in c("a1", "a0")) {
      row <- est$intervention == nm
      q <- stats::quantile(bt[, nm], c(0.025, 0.975), names = FALSE)
      expect_equal(
        c(est$ci_lower[row], est$ci_upper[row]),
        q,
        tolerance = 1e-12
      )
    }
  }
})

test_that("SNM honours boot_ci for the blip parameters (Path A)", {
  set.seed(7)
  n <- 1500
  L <- stats::rnorm(n)
  A <- stats::rbinom(n, 1L, stats::plogis(0.2 * L))
  Y <- 0.5 + 1.0 * A + 0.8 * A * L + 0.4 * L + stats::rnorm(n)
  d <- data.frame(Y = Y, A = A, L = L)
  fit <- causat(
    d,
    "Y",
    "A",
    confounders = ~ L + A:L,
    estimator = "snm",
    propensity_model_fn = stats::glm
  )
  set.seed(2)
  rp <- contrast(
    fit,
    ci_method = "bootstrap",
    n_boot = 400,
    boot_ci = "percentile"
  )
  set.seed(2)
  rn <- contrast(fit, ci_method = "bootstrap", n_boot = 400, boot_ci = "normal")
  ep <- as.data.frame(rp$estimates)
  en <- as.data.frame(rn$estimates)

  expect_equal(rp$boot_ci, "percentile")
  # Normal blip CIs are the Wald interval; point and SE are invariant.
  wald <- en$estimate + c(-1, 1) * z975 * en$se
  expect_equal(
    c(en$ci_lower[1], en$ci_upper[1]),
    en$estimate[1] + c(-1, 1) * z975 * en$se[1],
    tolerance = 1e-12
  )
  expect_equal(ep$estimate, en$estimate)
  expect_equal(ep$se, en$se)
  # Percentile and normal differ (the percentile interval is asymmetric).
  expect_false(isTRUE(all.equal(ep$ci_lower, en$ci_lower)))
})

test_that("an invalid boot_ci value is rejected", {
  d <- sim_boot_ci_dgp(n = 400)
  fit <- causat(
    d,
    "Y",
    "A",
    confounders = ~L,
    estimator = "gcomp",
    family = "binomial"
  )
  expect_error(
    contrast(
      fit,
      list(a1 = static(1), a0 = static(0)),
      ci_method = "bootstrap",
      n_boot = 50,
      boot_ci = "bogus"
    ),
    class = "rlang_error"
  )
})
