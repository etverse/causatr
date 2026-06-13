# Truth + external-oracle helpers for multinomial-outcome g-computation.
#
# The outcome is a single K-level factor drawn from a known softmax. Because
# the data-generating softmax is known, the g-computation estimand under
# do(A = a*) -- E_L[ P(Y = k | A = a*, L) ] -- has a closed form we evaluate by
# large-n Monte Carlo over the same DGP. `marginaleffects::avg_predictions()`
# on causatr's own fitted `nnet::multinom` object gives a second, exact oracle
# for the point estimates (same fitted model, same standardisation), and
# `avg_comparisons()` covers the per-class difference contrasts.

# Softmax of an n x K matrix of class linear predictors (reference column 0).
mo_softmax <- function(eta) {
  # Subtract the row max for numerical stability before exponentiating.
  e <- exp(eta - apply(eta, 1L, max))
  e / rowSums(e)
}

# Draw a K-level factor outcome from per-row class probabilities.
mo_draw_y <- function(probs, labels) {
  idx <- apply(probs, 1L, function(p) sample.int(length(p), 1L, prob = p))
  factor(labels[idx], levels = labels)
}

# K = 3 class linear predictors with the treatment A and confounder L entering
# linearly. eta_1 = 0 fixes the reference class ("none"). Works for binary or
# continuous A.
mo_eta3 <- function(A, L) {
  cbind(
    0,
    -0.4 + 0.9 * A + 0.6 * L,
    -0.8 + 0.5 * A - 0.4 * L
  )
}
mo_labels3 <- c("none", "mild", "severe")

# K = 4 class linear predictors, for the schema-generalisation tests.
mo_eta4 <- function(A, L) {
  cbind(
    0,
    -0.3 + 0.8 * A + 0.5 * L,
    -0.6 + 0.4 * A - 0.3 * L,
    -1.1 + 1.1 * A + 0.2 * L
  )
}
mo_labels4 <- c("c0", "c1", "c2", "c3")

# Binary-treatment K = 3 DGP. A ~ Bernoulli(plogis(0.3 L)); Y ~ softmax(eta3).
sim_multinom_binary <- function(
  n = 6000,
  seed = 11,
  eta_fn = mo_eta3,
  labels = mo_labels3
) {
  set.seed(seed)
  L <- stats::rnorm(n)
  A <- stats::rbinom(n, 1L, stats::plogis(0.3 * L))
  Y <- mo_draw_y(mo_softmax(eta_fn(A, L)), labels)
  data.frame(Y = Y, A = A, L = L)
}

# Continuous-treatment K = 3 DGP. A ~ N(0.5 L, 1); Y ~ softmax(eta3).
sim_multinom_continuous <- function(
  n = 6000,
  seed = 12,
  eta_fn = mo_eta3,
  labels = mo_labels3
) {
  set.seed(seed)
  L <- stats::rnorm(n)
  A <- 0.5 * L + stats::rnorm(n)
  Y <- mo_draw_y(mo_softmax(eta_fn(A, L)), labels)
  data.frame(Y = Y, A = A, L = L)
}

# Categorical-treatment (3-level) K = 3 outcome DGP. The treatment is an
# unordered factor entering the outcome softmax through level effects.
sim_multinom_cat_trt <- function(n = 7000, seed = 13) {
  set.seed(seed)
  L <- stats::rnorm(n)
  A <- factor(
    sample(c("lo", "mid", "hi"), n, replace = TRUE),
    levels = c("lo", "mid", "hi")
  )
  a_eff <- c(lo = 0, mid = 0.8, hi = 1.6)[as.character(A)]
  eta <- cbind(
    0,
    -0.4 + a_eff + 0.6 * L,
    -0.8 + 0.5 * a_eff - 0.4 * L
  )
  Y <- mo_draw_y(mo_softmax(eta), mo_labels3)
  data.frame(Y = Y, A = A, L = L)
}

# G-computation truth E_L[ softmax(eta(intervene(A, L), L)) ] by large-n MC.
# `gen_al(n)` draws (A, L) from the DGP joint (so shift MTPs marginalise over
# the observed treatment), and `intervene(A, L)` returns the counterfactual
# treatment column. Returns a named per-class probability vector.
mo_gcomp_truth <- function(
  gen_al,
  eta_fn,
  labels,
  intervene,
  n = 2e6,
  seed = 99
) {
  set.seed(seed)
  al <- gen_al(n)
  a_star <- intervene(al$A, al$L)
  stats::setNames(colMeans(mo_softmax(eta_fn(a_star, al$L))), labels)
}

# (A, L) joint draws matching the binary / continuous DGPs above.
mo_gen_al_binary <- function(n) {
  L <- stats::rnorm(n)
  list(A = stats::rbinom(n, 1L, stats::plogis(0.3 * L)), L = L)
}
mo_gen_al_continuous <- function(n) {
  L <- stats::rnorm(n)
  list(A = 0.5 * L + stats::rnorm(n), L = L)
}

# External oracle: per-class g-computation marginal means from
# `marginaleffects::avg_predictions()` on a fitted model, standardising over
# `data` with the treatment column set to `value`. Returns a named vector
# keyed by class label. Run on causatr's own `fit$model` for an exact match.
me_multinom_means <- function(model, data, trt_col, value) {
  nd <- as.data.frame(data)
  nd[[trt_col]] <- value
  res <- marginaleffects::avg_predictions(model, newdata = nd, type = "probs")
  stats::setNames(res$estimate, as.character(res$group))
}
