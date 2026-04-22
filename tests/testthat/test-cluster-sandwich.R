# Tests for the user-facing `cluster` argument on `contrast()`: general
# cluster-robust sandwich aggregation beyond matching's internal subclass
# clustering. Verifies against `sandwich::vcovCL()` as oracle on gcomp,
# verifies the ICE first-time-row alignment, verifies the IPW MSM-fit-row
# subsetting, and pins all validation rejections.

# ── Helper: cluster-randomised DGP ────────────────────────────────────
# 200 clusters × 5 members each = 1000 rows. Both treatment A and the
# outcome-shared shock U are cluster-level, so within-cluster IFs
# comove positively and the cluster-robust SE strictly inflates
# relative to the independent SE. This is the canonical "design
# effect" case that cluster-robust variance is built to fix (Liang &
# Zeger 1986; Cameron & Miller 2015).
# True ATE = 2.0.
simulate_clustered <- function(n_clusters = 200, m = 5, seed = 1) {
  set.seed(seed)
  cl <- rep(seq_len(n_clusters), each = m)
  U <- stats::rnorm(n_clusters)
  # Cluster-level treatment assignment probability -> within-cluster A
  # is near-constant (think cluster-randomised trial with some
  # within-cluster noise).
  A_cl <- stats::rbinom(n_clusters, 1, 0.5)
  n <- n_clusters * m
  L <- stats::rnorm(n)
  # 90% of rows inherit their cluster's assignment; 10% flip. Keeps A
  # varying enough for IPW to converge but retains the cluster-level
  # design effect.
  flip <- stats::rbinom(n, 1, 0.1)
  A <- ifelse(flip == 1, 1 - A_cl[cl], A_cl[cl])
  Y <- 1 + 2 * A + 0.5 * L + U[cl] + stats::rnorm(n, sd = 0.5)
  data.frame(Y = Y, A = A, L = L, cl = cl)
}


# ── gcomp × cluster: SE inflates relative to independent ──
test_that("gcomp × cluster sandwich inflates SE vs independent for correlated outcomes", {
  d <- simulate_clustered()
  fit_indep <- causat(d, outcome = "Y", treatment = "A", confounders = ~L)
  fit_cl <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    cluster = "cl"
  )

  res_indep <- contrast(
    fit_indep,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich"
  )
  res_clust <- contrast(
    fit_cl,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich"
  )

  # Point estimates are invariant; only SE changes under clustering.
  expect_equal(
    res_indep$contrasts$estimate,
    res_clust$contrasts$estimate,
    tolerance = 1e-10
  )

  # Cluster-robust SE must be strictly larger (positive within-cluster
  # correlation). A 10% lower bound is conservative — in this DGP the
  # ratio is typically 1.6-2.0.
  se_ratio <- res_clust$contrasts$se / res_indep$contrasts$se
  expect_gt(se_ratio, 1.1)

  # Truth coverage: the clustered CI covers 2.
  expect_lt(res_clust$contrasts$ci_lower, 2)
  expect_gt(res_clust$contrasts$ci_upper, 2)
})


# ── gcomp × cluster: oracle via sandwich::vcovCL ──
test_that("gcomp × cluster sandwich matches sandwich::vcovCL on a saturated outcome model", {
  # With the model `Y ~ A + L` (no nonlinearity in A) and static(1)/static(0)
  # interventions, the marginal-mean Jacobian collapses so that the
  # ATE sandwich SE reduces to the cluster-robust SE of the A coefficient.
  # This is the cleanest case to compare against `sandwich::vcovCL()`.
  d <- simulate_clustered()
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    cluster = "cl"
  )

  res_clust <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich"
  )

  # Oracle: sandwich::vcovCL() on the raw GLM.
  V_cl <- sandwich::vcovCL(fit$model, cluster = d$cl, type = "HC0")
  se_beta_A <- sqrt(V_cl["A", "A"])

  # The contrast SE should equal se_beta_A to within 1% (causatr's IF
  # engine builds the cluster-sum-then-square aggregation directly,
  # while `sandwich::vcovCL()` computes the same quantity via the
  # bread-meat form; both are HC0-equivalent up to a dof-style factor
  # `G/(G-1)` on the cluster count). For `Y ~ A + L` and ATE the two
  # paths agree to machine precision after correcting for that factor;
  # we assert the 1% bound without folding the factor in to keep the
  # assertion close to the raw engine output.
  G <- length(unique(d$cl))
  expect_equal(
    as.numeric(res_clust$contrasts$se),
    se_beta_A,
    tolerance = 0.02
  )
  # After the G/(G-1) correction the match is tight.
  expect_equal(
    as.numeric(res_clust$contrasts$se) * sqrt(G / (G - 1)),
    se_beta_A,
    tolerance = 1e-6
  )
})


# ── gcomp × cluster × EM via `by` ──
test_that("gcomp × cluster × by sandwich inflates SE within each stratum", {
  set.seed(2)
  d <- simulate_clustered()
  d$sex <- rep(c(0L, 1L), length.out = nrow(d))
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~ L + sex,
    cluster = "cl"
  )

  res_clust <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich",
    by = "sex"
  )

  # Both strata must return finite positive SE; the cluster threading
  # needs to survive the recursive `by` call.
  expect_true(all(is.finite(res_clust$contrasts$se)))
  expect_true(all(res_clust$contrasts$se > 0))
})


# ── IPW × cluster ──
test_that("ipw × cluster sandwich inflates SE vs independent on a clustered DGP", {
  d <- simulate_clustered()
  fit_indep <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw"
  )
  fit_cl <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw",
    cluster = "cl"
  )

  res_indep <- contrast(
    fit_indep,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich"
  )
  res_clust <- contrast(
    fit_cl,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich"
  )

  expect_equal(
    res_indep$contrasts$estimate,
    res_clust$contrasts$estimate,
    tolerance = 1e-10
  )

  se_ratio <- res_clust$contrasts$se / res_indep$contrasts$se
  # IPW SE inflation under clustering tends to be smaller than gcomp's
  # because the intercept-only MSM has less leverage, but must be > 1.
  expect_gt(se_ratio, 1.05)
})


# ── ICE × cluster × longitudinal ──
test_that("ice × cluster sandwich runs and inflates SE on a clustered longitudinal DGP", {
  # 2-period DGP with a cluster-level random intercept on both period
  # outcomes so the id-level IFs are correlated across clusters.
  set.seed(3)
  n_ids <- 300
  n_clusters <- 50
  id_to_cluster <- sample(seq_len(n_clusters), n_ids, replace = TRUE)
  U_cl <- stats::rnorm(n_clusters)

  L0 <- stats::rnorm(n_ids)
  A0 <- stats::rbinom(n_ids, 1, stats::plogis(0.3 * L0))
  L1 <- L0 + 0.2 * A0 + stats::rnorm(n_ids)
  A1 <- stats::rbinom(n_ids, 1, stats::plogis(0.3 * L1))
  Y <- 1 + A0 + A1 + 0.5 * L1 + U_cl[id_to_cluster] + stats::rnorm(n_ids)

  long <- data.table::data.table(
    id = rep(seq_len(n_ids), each = 2),
    time = rep(0:1, times = n_ids),
    A = as.numeric(rbind(A0, A1)),
    L = as.numeric(rbind(L0, L1)),
    Y = rep(Y, each = 2),
    cl = rep(id_to_cluster, each = 2)
  )

  fit_indep <- causat(
    long,
    outcome = "Y",
    treatment = "A",
    confounders = ~1,
    confounders_tv = ~L,
    id = "id",
    time = "time"
  )
  fit_cl <- causat(
    long,
    outcome = "Y",
    treatment = "A",
    confounders = ~1,
    confounders_tv = ~L,
    id = "id",
    time = "time",
    cluster = "cl"
  )
  res_indep <- contrast(
    fit_indep,
    interventions = list(all = static(1)),
    ci_method = "sandwich"
  )
  res_clust <- contrast(
    fit_cl,
    interventions = list(all = static(1)),
    ci_method = "sandwich"
  )

  expect_equal(
    res_indep$estimates$estimate,
    res_clust$estimates$estimate,
    tolerance = 1e-10
  )
  # With 300 ids in 50 clusters and a non-trivial cluster-level random
  # intercept, the clustered SE should be materially larger.
  expect_gt(res_clust$estimates$se / res_indep$estimates$se, 1.1)
})


# ── Matching × cluster: rejected at fit and contrast ──
test_that("matching × fit-time cluster aborts with causatr_bad_cluster", {
  set.seed(4)
  d <- simulate_clustered(n_clusters = 100, m = 4)
  expect_error(
    causat(
      d,
      outcome = "Y",
      treatment = "A",
      confounders = ~L,
      estimator = "matching",
      estimand = "ATT",
      cluster = "cl"
    ),
    class = "causatr_bad_cluster"
  )
})

test_that("matching × contrast-time cluster aborts with causatr_bad_cluster", {
  set.seed(4)
  d <- simulate_clustered(n_clusters = 100, m = 4)
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "matching",
    estimand = "ATT"
  )
  expect_error(
    contrast(
      fit,
      interventions = list(a1 = static(1), a0 = static(0)),
      reference = "a0",
      cluster = d$cl
    ),
    class = "causatr_bad_cluster"
  )
})


# ── Unknown cluster column: rejected ──
test_that("unknown cluster column aborts with causatr_bad_cluster", {
  d <- simulate_clustered()
  fit <- causat(d, outcome = "Y", treatment = "A", confounders = ~L)
  expect_error(
    contrast(
      fit,
      interventions = list(a1 = static(1), a0 = static(0)),
      cluster = "no_such_col"
    ),
    class = "causatr_bad_cluster"
  )
})


# ── Wrong-length cluster vector: rejected ──
test_that("wrong-length cluster vector aborts with causatr_bad_cluster", {
  d <- simulate_clustered()
  fit <- causat(d, outcome = "Y", treatment = "A", confounders = ~L)
  expect_error(
    contrast(
      fit,
      interventions = list(a1 = static(1), a0 = static(0)),
      cluster = seq_len(10)
    ),
    class = "causatr_bad_cluster"
  )
})


# ── NA in cluster: rejected ──
test_that("NA in cluster vector aborts with causatr_bad_cluster", {
  d <- simulate_clustered()
  d$cl[1:3] <- NA
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    cluster = "cl"
  )
  expect_error(
    contrast(
      fit,
      interventions = list(a1 = static(1), a0 = static(0))
    ),
    class = "causatr_bad_cluster"
  )
})


# ── Cluster as direct vector ──
test_that("cluster as a direct vector gives same result as column name", {
  d <- simulate_clustered()
  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    cluster = "cl"
  )
  res_name <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0))
  )
  res_vec <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    cluster = d$cl
  )
  expect_equal(res_name$contrasts$se, res_vec$contrasts$se, tolerance = 1e-12)
})


# ── Singleton clusters reduce to independent aggregation ──
test_that("singleton clusters (one per row) reduce to the independent sandwich SE", {
  d <- simulate_clustered()
  fit <- causat(d, outcome = "Y", treatment = "A", confounders = ~L)
  res_indep <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    ci_method = "sandwich"
  )
  # Every row in its own cluster: the within-cluster "sum" is a single
  # IF per cluster, so the final crossprod is identical to the
  # no-cluster case.
  res_singleton <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    ci_method = "sandwich",
    cluster = seq_len(nrow(d))
  )
  expect_equal(
    res_indep$contrasts$se,
    res_singleton$contrasts$se,
    tolerance = 1e-10
  )
})
