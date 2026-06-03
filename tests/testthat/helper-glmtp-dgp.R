# Data-generating processes for natural-history modified treatment policies
# (G-LMTPs; Diaz, Williams, Morzywolek & Rudolph 2026, arXiv:2605.24167).
#
# The defining feature is treatment-state feedback: the *natural value* of
# treatment at period t depends on the (intervened) treatment at t-1 through an
# absorbing exposure process, so the natural history under a delay regime
# differs from the observed history. The standard ICE recursion (which
# conditions on the OBSERVED lag) therefore targets the wrong estimand; only the
# augmented-data sequential regression recovers the true counterfactual mean.
#
# Truth oracle: forward Monte-Carlo simulation of the natural-history regime
# (the paper's Proposition 1 with correct models on a large super-population).
# `lmtp`'s one-shot shift is NOT a valid oracle here -- it computes the standard
# LMTP, which differs from the natural-history target whenever the policy has
# perturbed the trajectory.

#' Delay-by-window policy on a binary absorbing treatment
#'
#' Returns A^d_t = 1 once natural initiation has occurred by period t - window,
#' given the natural-value history `s_hist = (s_1, ..., s_t)`. The indicator is
#' monotone in t, so treatment is absorbing once it switches on. Matches the
#' policy closure built by `grace_period(window)`.
glmtp_delay_policy <- function(s_hist, window) {
  t <- length(s_hist)
  k <- t - window
  if (k <= 0L) {
    return(0)
  }
  as.numeric(any(s_hist[seq_len(k)] == 1))
}

# Default structural coefficients for the linear-gaussian delay SCM. Kept in one
# place so the observed-data generator and the forward-MC truth use identical
# parameters (any drift between them would invalidate the truth check).
glmtp_delay_params <- function() {
  list(
    tau = 4L,
    rho = 0.5, # L1 = rho * L0 + eps
    aL = 0.5, # L_t = aL * L_{t-1} + gamma * A^actual_{t-1} + eps
    gamma = 0.8, # treatment -> covariate feedback (irreducibility driver)
    sdL = 1.0,
    a0 = -1.0, # natural propensity intercept (when not yet absorbed)
    a1 = 0.5, # natural propensity covariate slope
    b0 = 1.0, # outcome intercept
    bA = 0.8, # per-period treatment effect on Y
    bL = 0.4, # final-covariate effect on Y
    bL0 = -0.5, # baseline-covariate effect on Y
    sdY = 1.0
  )
}

#' Linear-gaussian delay SCM with binary absorbing treatment (observed data)
#'
#' @description
#' Generates a longitudinal person-period data set whose binary treatment is
#' absorbing (`A_t = 1` once `A_{t-1} = 1`) and whose time-varying covariate
#' feeds back from the prior treatment. Because the SCM is linear-gaussian and
#' the grace/delay policy value is a constant given the natural-history label,
#' every per-step augmented regression is exactly linear -- so the additive-GLM
#' augmented plug-in is correctly specified and consistent for the
#' forward-MC truth.
#'
#' @param n Integer. Number of individuals.
#' @param seed Integer. RNG seed.
#' @param tau Integer. Number of periods.
#'
#' @return A list with `data` (long data.frame: `id`, `t`, `L0`, `A`, `L`, `Y`
#'   observed only at the final period) and `params`.
glmtp_delay_data <- function(n = 4000L, seed = 1L, tau = NULL) {
  p <- glmtp_delay_params()
  if (!is.null(tau)) {
    p$tau <- as.integer(tau)
  }
  tau <- p$tau
  set.seed(seed)

  L0 <- stats::rnorm(n)
  A <- matrix(0L, n, tau)
  L <- matrix(0, n, tau)
  Lprev <- L0
  Aprev <- integer(n)
  for (t in seq_len(tau)) {
    eps <- stats::rnorm(n, sd = p$sdL)
    Lt <- if (t == 1L) {
      p$rho * L0 + eps
    } else {
      p$aL * Lprev + p$gamma * Aprev + eps
    }
    pt <- stats::plogis(p$a0 + p$a1 * Lt)
    # Absorbing natural treatment: stay on once initiated.
    At <- ifelse(Aprev == 1L, 1L, stats::rbinom(n, 1L, pt))
    A[, t] <- At
    L[, t] <- Lt
    Lprev <- Lt
    Aprev <- At
  }
  cumA <- rowSums(A)
  Y <- p$b0 +
    p$bA * cumA +
    p$bL * L[, tau] +
    p$bL0 * L0 +
    stats::rnorm(n, sd = p$sdY)

  long <- data.frame(
    id = rep(seq_len(n), each = tau),
    t = rep(seq_len(tau), n),
    L0 = rep(L0, each = tau),
    A = as.vector(t(A)),
    L = as.vector(t(L)),
    Y = NA_real_
  )
  long$Y[long$t == tau] <- Y
  list(data = long, params = p)
}

#' Forward Monte-Carlo truth for the delay regime (Proposition 1)
#'
#' @description
#' Simulates the natural-history regime forward under the delay-by-`window`
#' policy on a large super-population with the *correct* structural models, and
#' returns the true counterfactual mean \eqn{E[Y^d]}. At each period the natural
#' value is drawn from its conditional law given the **intervened** history
#' (so the absorbing process is anchored to A^d_{t-1}), the policy maps the
#' natural history to the intervened value, and the next covariate is generated
#' from the **actual** (intervened) treatment. This is exactly the paper's
#' Proposition-1 truth computation.
#'
#' @param window Integer delay.
#' @param n_mc Integer super-population size.
#' @param seed Integer RNG seed.
#' @param params Structural coefficients (defaults to [glmtp_delay_params()]).
#'
#' @return Numeric scalar, the true \eqn{E[Y^d]}.
glmtp_delay_forward_truth <- function(
  window,
  n_mc = 2e6,
  seed = 1L,
  params = NULL
) {
  p <- params %||% glmtp_delay_params()
  tau <- p$tau
  set.seed(seed)

  L0 <- stats::rnorm(n_mc)
  Lprev <- L0
  Ad_prev <- integer(n_mc)
  s_hist <- matrix(0L, n_mc, tau) # natural values under the regime
  Ad <- matrix(0L, n_mc, tau) # intervened values
  L_at_tau <- numeric(n_mc)
  for (t in seq_len(tau)) {
    eps <- stats::rnorm(n_mc, sd = p$sdL)
    Lt <- if (t == 1L) {
      p$rho * L0 + eps
    } else {
      p$aL * Lprev + p$gamma * Ad_prev + eps
    }
    if (t == tau) {
      L_at_tau <- Lt
    }
    # Natural value: absorbing given the INTERVENED prior treatment.
    p_nat <- stats::plogis(p$a0 + p$a1 * Lt)
    s_t <- ifelse(Ad_prev == 1L, 1L, stats::rbinom(n_mc, 1L, p_nat))
    s_hist[, t] <- s_t
    # Delay policy: A^d_t = 1[ natural initiation by t - window ].
    k <- t - window
    Ad_t <- if (k <= 0L) {
      integer(n_mc)
    } else {
      as.integer(glmtp_row_any_initiated(s_hist[, seq_len(k), drop = FALSE]))
    }
    Ad[, t] <- Ad_t
    Lprev <- Lt
    Ad_prev <- Ad_t
  }
  cumAd <- rowSums(Ad)
  mean(p$b0 + p$bA * cumAd + p$bL * L_at_tau + p$bL0 * L0)
}

# Row-wise "ever initiated by k": TRUE if any of the first k natural values in
# the row equals 1. Base-R `rowSums(mat == 1L) > 0` over the history slice.
glmtp_row_any_initiated <- function(mat) {
  rowSums(mat == 1L) > 0L
}

# ---------------------------------------------------------------------------
# Paper replication: Diaz et al. (2026) Section-6 mechanism, adapted to a
# single end-of-study outcome (causatr's ICE targets one final outcome, not the
# paper's time-varying absorbing survival process). The treatment propensity
# logit(-1.5 + 0.3 L_t) (absorbing once initiated) and covariate
# L_t = 0.5 L_{t-1} + eps are taken verbatim from the paper; the delay-by-one
# policy is the paper's intervention. The outcome is a single binary survival
# indicator at tau driven by cumulative treatment and the final covariate.
# ---------------------------------------------------------------------------

glmtp_paper_params <- function() {
  list(
    tau = 5L,
    aL = 0.5, # L_t = 0.5 L_{t-1} + eps          (paper)
    sdL = 1.0,
    a0 = -1.5, # logit(-1.5 + 0.3 L_t)            (paper propensity)
    a1 = 0.3,
    gamma = 0.0, # paper has no A -> L feedback; irreducibility is via the
    # absorbing propensity anchored to the intervened prior treatment
    b0 = -0.5, # single end-of-study binary outcome (survival-flavoured)
    bA = 0.5,
    bL = -0.4
  )
}

#' Paper Section-6 mechanism with a single end-of-study binary outcome
#'
#' @param n Integer number of individuals.
#' @param seed Integer RNG seed.
#' @return A list with `data` (long data.frame) and `params`.
glmtp_paper_data <- function(n = 5000L, seed = 1L) {
  p <- glmtp_paper_params()
  tau <- p$tau
  set.seed(seed)

  L0 <- stats::rnorm(n)
  A <- matrix(0L, n, tau)
  L <- matrix(0, n, tau)
  Lprev <- L0
  Aprev <- integer(n)
  for (t in seq_len(tau)) {
    Lt <- p$aL * Lprev + stats::rnorm(n, sd = p$sdL)
    pt <- stats::plogis(p$a0 + p$a1 * Lt)
    At <- ifelse(Aprev == 1L, 1L, stats::rbinom(n, 1L, pt))
    A[, t] <- At
    L[, t] <- Lt
    Lprev <- Lt
    Aprev <- At
  }
  cumA <- rowSums(A)
  py <- stats::plogis(p$b0 + p$bA * cumA + p$bL * L[, tau])
  Y <- stats::rbinom(n, 1L, py)
  long <- data.frame(
    id = rep(seq_len(n), each = tau),
    t = rep(seq_len(tau), n),
    L0 = rep(L0, each = tau),
    A = as.vector(t(A)),
    L = as.vector(t(L)),
    Y = NA_real_
  )
  long$Y[long$t == tau] <- Y
  list(data = long, params = p)
}

#' Forward-MC truth for the paper mechanism under delay-by-window
#'
#' @param window Integer delay.
#' @param n_mc Integer super-population size.
#' @param seed Integer RNG seed.
#' @return Numeric scalar, true \eqn{E[Y^d]} on the probability scale.
glmtp_paper_forward_truth <- function(window, n_mc = 2e6, seed = 1L) {
  p <- glmtp_paper_params()
  tau <- p$tau
  set.seed(seed)
  L0 <- stats::rnorm(n_mc)
  Lprev <- L0
  Ad_prev <- integer(n_mc)
  s_hist <- matrix(0L, n_mc, tau)
  Ad <- matrix(0L, n_mc, tau)
  L_tau <- numeric(n_mc)
  for (t in seq_len(tau)) {
    Lt <- p$aL * Lprev + stats::rnorm(n_mc, sd = p$sdL)
    if (t == tau) {
      L_tau <- Lt
    }
    p_nat <- stats::plogis(p$a0 + p$a1 * Lt)
    s_t <- ifelse(Ad_prev == 1L, 1L, stats::rbinom(n_mc, 1L, p_nat))
    s_hist[, t] <- s_t
    k <- t - window
    Ad_t <- if (k <= 0L) {
      integer(n_mc)
    } else {
      as.integer(glmtp_row_any_initiated(s_hist[, seq_len(k), drop = FALSE]))
    }
    Ad[, t] <- Ad_t
    Lprev <- Lt
    Ad_prev <- Ad_t
  }
  cumAd <- rowSums(Ad)
  mean(stats::plogis(p$b0 + p$bA * cumAd + p$bL * L_tau))
}

# ---------------------------------------------------------------------------
# Dose-escalation cap (ordered / ordinal treatment), tau = 3, support {0,1,2}.
# A dose escalates from the ACTUAL (capped) prior dose; the cap policy limits the
# natural per-period increase at delta, comparing the natural dose at t to the
# natural dose at t-1. tau = 3 makes the policy genuinely irreducible. Used only
# to validate the engine's ordinal natural-value-dependent path against an
# independent hand-coded recursion (the public cap_escalation() release is
# deferred -- see PHASE_22). The dose enters the outcome model numerically.
# ---------------------------------------------------------------------------

#' Ordered-dose escalation SCM, tau = 3, support {0,1,2} (observed data)
#'
#' @param n Integer number of individuals.
#' @param seed Integer RNG seed.
#' @return A list with `data` (long data.frame; dose column `A`) and `support`.
glmtp_ordinal_cap_data <- function(n = 6000L, seed = 1L) {
  tau <- 3L
  m <- 2L
  set.seed(seed)
  L0 <- stats::rnorm(n)
  A <- matrix(0L, n, tau)
  L <- matrix(0, n, tau)
  Lprev <- L0
  Aprev <- integer(n)
  for (t in seq_len(tau)) {
    Lt <- 0.5 * Lprev + 0.3 * Aprev + stats::rnorm(n)
    step <- stats::rbinom(n, 2L, stats::plogis(-0.2 + 0.6 * Lt))
    At <- pmin(m, Aprev + step)
    A[, t] <- At
    L[, t] <- Lt
    Lprev <- Lt
    Aprev <- At
  }
  Y <- 1 + 0.7 * rowSums(A) + 0.4 * L[, tau] - 0.5 * L0 + stats::rnorm(n)
  long <- data.frame(
    id = rep(seq_len(n), each = tau),
    t = rep(seq_len(tau), n),
    L0 = rep(L0, each = tau),
    A = as.vector(t(A)),
    L = as.vector(t(L)),
    Y = NA_real_
  )
  long$Y[long$t == tau] <- Y
  list(data = long, support = c(0, 1, 2))
}

#' Forward-MC truth for the dose-escalation cap regime (Proposition 1)
#'
#' @description
#' Simulates the natural-history regime forward under `cap_escalation(delta)` on
#' the `glmtp_ordinal_cap_data()` SCM (tau = 3, support {0,1,2}) with the
#' *correct* structural models on a large super-population, returning the true
#' counterfactual mean \eqn{E[Y^d]}. At each period the natural dose escalates
#' from the **intervened** prior dose (the absorbing/escalating process is
#' anchored to \eqn{A^d_{t-1}}), the cap policy limits the natural per-period
#' increase by comparing the natural dose at \eqn{t} to the natural dose at
#' \eqn{t-1}, and the covariate feeds back from the **actual** (intervened)
#' treatment. This is the paper's Proposition-1 truth computation for the cap
#' policy and the oracle that the flexible-treatment (`~ factor(A)`) plug-in is
#' validated against (the bare-numeric plug-in is misspecified for the kinked
#' capped dose).
#'
#' @param delta Numeric cap on the per-period natural increase.
#' @param n_mc Integer super-population size.
#' @param seed Integer RNG seed.
#'
#' @return Numeric scalar, the true \eqn{E[Y^d]}.
glmtp_cap_forward_truth <- function(delta, n_mc = 2e6, seed = 1L) {
  tau <- 3L
  m <- 2L
  set.seed(seed)
  L0 <- stats::rnorm(n_mc)
  Lprev <- L0
  Ad_prev <- integer(n_mc) # intervened dose at t-1
  s_prev <- integer(n_mc) # natural dose at t-1
  Ad <- matrix(0L, n_mc, tau) # intervened doses under the regime
  L_tau <- numeric(n_mc)
  for (t in seq_len(tau)) {
    # Covariate feeds back from the intervened prior dose (matches the
    # observed SCM, where it feeds from the actual prior dose).
    Lt <- 0.5 * Lprev + 0.3 * Ad_prev + stats::rnorm(n_mc)
    if (t == tau) {
      L_tau <- Lt
    }
    # Natural dose at t escalates from the intervened prior dose.
    step <- stats::rbinom(n_mc, 2L, stats::plogis(-0.2 + 0.6 * Lt))
    s_t <- pmin(m, Ad_prev + step)
    # Cap the increase of the natural dose at t over the natural dose at t-1.
    Ad_t <- if (t == 1L) {
      s_t
    } else {
      ifelse(s_t - s_prev > delta, s_prev + delta, s_t)
    }
    Ad[, t] <- Ad_t
    Lprev <- Lt
    Ad_prev <- Ad_t
    s_prev <- s_t
  }
  mean(1 + 0.7 * rowSums(Ad) + 0.4 * L_tau - 0.5 * L0)
}

#' Independent hand-coded augmented recursion for the dose cap (tau = 3)
#'
#' @description
#' A from-scratch re-implementation of the augmented-data sequential regression
#' for the dose-escalation cap on a 3-period, 3-level treatment, structured
#' differently from [glmtp_iterate()] (explicit nested label loops, bare-numeric
#' `stats::glm`). It is the independent oracle proving the engine has no bug on
#' the ordinal natural-value-dependent path: the two must agree to numerical
#' precision. The cap compares the natural dose at \eqn{t} (here the value
#' predicted at) with the natural dose at \eqn{t-1} (the last label entry).
#'
#' @param fit A `causatr_fit` from `glmtp_ordinal_cap_data()` (lag columns built).
#' @param delta Numeric cap on the per-period natural increase.
#' @return Numeric scalar, the plug-in estimate \eqn{\hat\theta}.
glmtp_handcode_cap <- function(fit, delta) {
  support <- c(0, 1, 2)
  cap <- function(a_now, a_prev) {
    ifelse(a_now - a_prev > delta, a_prev + delta, a_now)
  }
  d <- as.data.frame(fit$data)
  by_t <- function(tt) d[d$t == tt, , drop = FALSE]
  # Same per-step formulas the ICE builder produces for tau = 3, history = Inf.
  f_base <- stats::as.formula(
    "Y ~ A + lag1_A + lag2_A + L + lag1_L + lag2_L + L0"
  )
  f_s2 <- stats::as.formula(".py ~ A + lag1_A + L + lag1_L + L0")
  f_s1 <- stats::as.formula(".py ~ A + L + L0")

  d3 <- by_t(3)
  m3 <- stats::glm(f_base, data = d3[!is.na(d3$Y), ])
  q3 <- list()
  for (s1 in support) {
    for (s2 in support) {
      nd <- d3
      nd$A <- cap(d3$A, s2) # cap natural A3 (=obs) vs natural A2 (=s2)
      q3[[paste(s1, s2)]] <- stats::setNames(
        as.numeric(stats::predict(m3, nd)),
        d3$id
      )
    }
  }
  d2 <- by_t(2)
  m2 <- list()
  for (s1 in support) {
    for (s2 in support) {
      k <- paste(s1, s2)
      dd <- d2
      dd$.py <- q3[[k]][as.character(d2$id)]
      m2[[k]] <- stats::glm(f_s2, data = dd)
    }
  }
  d1 <- by_t(1)
  q2 <- list()
  for (s1 in support) {
    pr <- rep(NA_real_, nrow(d2))
    for (a2 in support) {
      ix <- which(d2$A == a2)
      if (length(ix) == 0L) {
        next
      }
      nd <- d2[ix, ]
      nd$A <- cap(a2, s1) # cap natural A2 (=a2) vs natural A1 (=s1)
      pr[ix] <- as.numeric(stats::predict(m2[[paste(s1, a2)]], nd))
    }
    q2[[as.character(s1)]] <- stats::setNames(pr, d2$id)
  }
  m1 <- list()
  for (s1 in support) {
    dd <- d1
    dd$.py <- q2[[as.character(s1)]][as.character(d1$id)]
    m1[[as.character(s1)]] <- stats::glm(f_s1, data = dd)
  }
  q1 <- rep(NA_real_, nrow(d1))
  for (a1 in support) {
    ix <- which(d1$A == a1)
    if (length(ix) == 0L) {
      next
    }
    # d_1 = A_1 (no prior dose to cap against at baseline).
    q1[ix] <- as.numeric(stats::predict(m1[[as.character(a1)]], d1[ix, ]))
  }
  mean(q1)
}
