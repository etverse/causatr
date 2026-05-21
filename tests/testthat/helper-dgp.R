# Shared data-generating processes for simulation tests.
#
# DGP 1: Binary treatment, continuous outcome
#   L ~ N(0, 1)
#   A | L ~ Bernoulli(expit(0.5 * L))
#   Y | A, L ~ N(2 + 3*A + 1.5*L, sd = 1)
#
# True ATE  = E[Y(1)] - E[Y(0)] = 3
# True ATT  = 3 (constant treatment effect)
# True E[Y(1)] = 5, E[Y(0)] = 2

simulate_binary_continuous <- function(n = 2000, seed = 42) {
  set.seed(seed)
  L <- rnorm(n)
  ps <- plogis(0.5 * L)
  A <- rbinom(n, 1, ps)
  Y <- 2 + 3 * A + 1.5 * L + rnorm(n)
  data.frame(Y = Y, A = A, L = L)
}

# DGP 2: Binary treatment, binary outcome
#   L ~ N(0, 1)
#   A | L ~ Bernoulli(expit(0.5 * L))
#   Y | A, L ~ Bernoulli(expit(-1 + 1.5*A + 0.8*L))
#
# True risk under A=1: ≈ 0.622
# True risk under A=0: ≈ 0.289
# True RD ≈ 0.333

simulate_binary_binary <- function(n = 2000, seed = 42) {
  set.seed(seed)
  L <- rnorm(n)
  ps <- plogis(0.5 * L)
  A <- rbinom(n, 1, ps)
  Y <- rbinom(n, 1, plogis(-1 + 1.5 * A + 0.8 * L))
  data.frame(Y = Y, A = A, L = L)
}

# DGP 3: Continuous treatment, continuous outcome
#   L ~ N(0, 1)
#   A | L ~ N(1 + 0.5*L, sd = 1)
#   Y | A, L ~ N(1 + 2*A + L, sd = 1)
#
# True E[Y(a)] = 1 + 2*a
# shift(-1): E[Y(A-1)] vs E[Y(A)] → difference = -2

simulate_continuous_continuous <- function(n = 2000, seed = 42) {
  set.seed(seed)
  L <- rnorm(n)
  A <- 1 + 0.5 * L + rnorm(n)
  Y <- 1 + 2 * A + L + rnorm(n)
  data.frame(Y = Y, A = A, L = L)
}

# DGP 4: Binary treatment with effect modification by sex
#   L ~ N(0, 1); sex ~ Bernoulli(0.5)
#   A | L ~ Bernoulli(expit(0.5*L))
#   Y = 2 + 3*A + 1.5*L + 1.2*sex*A + N(0, 1)
#
# True ATE       = 3 + 1.2*E[sex] = 3.6
# True ATE|sex=0 = 3
# True ATE|sex=1 = 4.2

simulate_effect_mod <- function(n = 3000, seed = 42) {
  set.seed(seed)
  L <- rnorm(n)
  sex <- rbinom(n, 1, 0.5)
  A <- rbinom(n, 1, plogis(0.5 * L))
  Y <- 2 + 3 * A + 1.5 * L + 1.2 * sex * A + rnorm(n)
  data.frame(Y = Y, A = A, L = L, sex = sex)
}

# DGP 5: Binary treatment, continuous outcome, HETEROGENEOUS treatment effect
#   L ~ N(0, 1); sex ~ Bernoulli(0.5)
#   A | L ~ Bernoulli(expit(L))          [strong confounding → ATT ≠ ATE]
#   Y = 2 + (3 + 2*L)*A + L + 1.5*sex*A + N(0, 1)
#
# Treatment effect τ(L, sex) = 3 + 2*L + 1.5*sex
# Correct outcome model: Y ~ A + L + sex + A:L + A:sex
#
# Monte Carlo truth (n = 10^7, seed = 1):
#   ATE|sex=0 ≈ 3.00,  ATE|sex=1 ≈ 4.50
#   ATT|sex=0 ≈ 3.83,  ATT|sex=1 ≈ 5.33
#   ATC|sex=0 ≈ 2.17,  ATC|sex=1 ≈ 3.67

simulate_het_effect <- function(n = 5000, seed = 42) {
  set.seed(seed)
  L <- rnorm(n)
  sex <- rbinom(n, 1, 0.5)
  A <- rbinom(n, 1, plogis(L))
  Y <- 2 + (3 + 2 * L) * A + L + 1.5 * sex * A + rnorm(n)
  data.frame(Y = Y, A = A, L = L, sex = sex)
}

# DGP 6: Binary treatment, BINARY outcome, heterogeneous treatment effect
#   L ~ N(0, 1); sex ~ Bernoulli(0.5)
#   A | L ~ Bernoulli(expit(L))
#   Y ~ Bernoulli(expit(-1 + (1.5 + L)*A + 0.8*L + 0.5*sex*A))
#
# Correct outcome model: Y ~ A + L + sex + A:L + A:sex (logistic)
#
# Monte Carlo truth (n = 10^7, seed = 1):
#   ATE|sex=0 ≈ 0.287,  ATE|sex=1 ≈ 0.364
#   ATT|sex=0 ≈ 0.346,  ATT|sex=1 ≈ 0.415
#   ATC|sex=0 ≈ 0.228,  ATC|sex=1 ≈ 0.312

simulate_het_binary <- function(n = 5000, seed = 42) {
  set.seed(seed)
  L <- rnorm(n)
  sex <- rbinom(n, 1, 0.5)
  A <- rbinom(n, 1, plogis(L))
  Y <- rbinom(n, 1, plogis(-1 + (1.5 + L) * A + 0.8 * L + 0.5 * sex * A))
  data.frame(Y = Y, A = A, L = L, sex = sex)
}

# DGP 7: Categorical (3-level) treatment, continuous outcome
#   L ~ N(0, 1)
#   A | L ~ Multinomial("a","b","c") via softmax(0, 0.5*L, -0.3*L)
#   Y | A, L ~ N(mu_A + 1.5*L, sd = 1)
#     mu_A = 2 if A="a", 5 if A="b", 3 if A="c"
#
# True counterfactual means (marginalise over L ~ N(0,1)):
#   E[Y("a")] = 2,  E[Y("b")] = 5,  E[Y("c")] = 3
# True ATE("b" vs "a") = 3

simulate_categorical_continuous <- function(n = 5000, seed = 42) {
  set.seed(seed)
  L <- rnorm(n)
  # Multinomial logit: log-odds vs reference "a"
  eta_b <- 0.5 * L
  eta_c <- -0.3 * L
  denom <- 1 + exp(eta_b) + exp(eta_c)
  p_a <- 1 / denom
  p_b <- exp(eta_b) / denom
  p_c <- exp(eta_c) / denom
  probs <- cbind(p_a, p_b, p_c)
  A <- factor(
    apply(probs, 1, function(p) sample(c("a", "b", "c"), 1, prob = p)),
    levels = c("a", "b", "c")
  )
  mu <- ifelse(A == "a", 2, ifelse(A == "b", 5, 3))
  Y <- mu + 1.5 * L + rnorm(n)
  data.frame(Y = Y, A = A, L = L)
}


# Longitudinal DGPs

# Table 20.1 from Hernán & Robins (2025).
# 2 time points, true causal effect = 0.
# Demonstrates treatment-confounder feedback where naive methods give biased
# estimates but ICE g-computation correctly recovers ATE = 0.
#
# By default produces the full 32,000-individual dataset. The `scale` argument
# shrinks it proportionally for faster bootstrap tests.
make_table201 <- function(scale = 1) {
  groups <- data.frame(
    A0 = c(0, 0, 0, 0, 1, 1, 1, 1),
    L1 = c(0, 0, 1, 1, 0, 0, 1, 1),
    A1 = c(0, 1, 0, 1, 0, 1, 0, 1),
    Y = c(84, 84, 52, 52, 76, 76, 44, 44),
    N = as.integer(c(2400, 1600, 2400, 9600, 4800, 3200, 1600, 6400) * scale)
  )
  rows <- lapply(seq_len(nrow(groups)), function(i) {
    g <- groups[i, ]
    off <- sum(groups$N[seq_len(i - 1)])
    data.frame(
      id = seq_len(g$N) + off,
      A0 = g$A0,
      L1 = g$L1,
      A1 = g$A1,
      Y = g$Y
    )
  })
  wide <- do.call(rbind, rows)
  t0 <- data.frame(
    id = wide$id,
    time = 0L,
    A = wide$A0,
    L = NA_real_,
    Y = NA_real_
  )
  t1 <- data.frame(
    id = wide$id,
    time = 1L,
    A = wide$A1,
    L = wide$L1,
    Y = wide$Y
  )
  rbind(t0, t1)
}

# Linear SCM with treatment-confounder feedback and known analytical ATE.
#
# DGP:
#   L0 ~ N(0, 1)                              (time-invariant baseline confounder)
#   L_0 = 0.5 * L0 + ε_L                     (TV confounder at baseline)
#   A_0 ~ Bern(expit(0.5 * L0))               (treatment at t=0)
#   For t > 0:
#     L_t = A_{t-1} + 0.5 * L0 + ε_L          (treatment-confounder feedback)
#     A_t ~ Bern(expit(0.3 * L_t + 0.2 * A_{t-1}))
#   Y = 10 + 2 * sum(A_t) + sum(L_t) + ε_Y
#
# True ATE (always vs never) = 3 * n_times - 1
# Under always: E[Y] = 10 + 2*T + (T-1) = 9 + 3T
# Under never:  E[Y] = 10
# (L_0 does not depend on treatment, so cancels in the contrast.)
make_linear_scm <- function(n = 5000, n_times = 2, seed = 42) {
  set.seed(seed)

  L0 <- stats::rnorm(n)

  A <- matrix(NA_real_, n, n_times)
  L <- matrix(NA_real_, n, n_times)

  for (t in seq_len(n_times)) {
    if (t == 1) {
      L[, t] <- 0.5 * L0 + stats::rnorm(n, 0, 0.5)
      A[, t] <- stats::rbinom(n, 1, stats::plogis(0.5 * L0))
    } else {
      L[, t] <- A[, t - 1] + 0.5 * L0 + stats::rnorm(n, 0, 0.5)
      A[, t] <- stats::rbinom(
        n,
        1,
        stats::plogis(0.3 * L[, t] + 0.2 * A[, t - 1])
      )
    }
  }

  Y <- 10 + 2 * rowSums(A) + rowSums(L) + stats::rnorm(n)

  rows <- vector("list", n_times)
  for (t in seq_len(n_times)) {
    rows[[t]] <- data.frame(
      id = seq_len(n),
      time = t - 1L,
      A = A[, t],
      L = L[, t],
      L0 = L0,
      Y = if (t == n_times) Y else NA_real_
    )
  }
  do.call(rbind, rows)
}

# Continuous-treatment version of the linear SCM.
#
# DGP:
# Continuous-treatment version of the linear SCM.
#
# DGP:
#   L0 ~ N(0, 1)                        (time-invariant baseline confounder)
#   L_0 = 0.5 * L0 + ε_L               (TV confounder at baseline)
#   A_0 = 1 + 0.5 * L0 + ε_A           (continuous treatment)
#   L_1 = A_0 + 0.5 * L0 + ε_L         (treatment-confounder feedback)
#   A_1 = 1 + 0.3 * L_1 + 0.2 * A_0 + ε_A
#   Y = 10 + 2 * (A_0 + A_1) + L_0 + L_1 + ε_Y
make_continuous_scm <- function(n = 5000, seed = 42) {
  set.seed(seed)

  L0 <- stats::rnorm(n)
  Ltv0 <- 0.5 * L0 + stats::rnorm(n, 0, 0.5)
  A0 <- 1 + 0.5 * L0 + stats::rnorm(n, 0, 0.5)
  L1 <- A0 + 0.5 * L0 + stats::rnorm(n, 0, 0.5)
  A1 <- 1 + 0.3 * L1 + 0.2 * A0 + stats::rnorm(n, 0, 0.5)
  Y <- 10 + 2 * (A0 + A1) + Ltv0 + L1 + stats::rnorm(n)

  rbind(
    data.frame(
      id = seq_len(n),
      time = 0L,
      A = A0,
      L = Ltv0,
      L0 = L0,
      Y = NA_real_
    ),
    data.frame(
      id = seq_len(n),
      time = 1L,
      A = A1,
      L = L1,
      L0 = L0,
      Y = Y
    )
  )
}

# DGP 8: Count (Poisson) treatment, continuous outcome
#   L ~ N(0, 1)
#   A | L ~ Poisson(exp(0.5 * L))
#   Y | A, L ~ N(2 + 1.5*A + L, sd = 1)
#
# True E[Y(shift(delta))] = 2 + 1.5*(E[A] + delta)
#   where E[A] = E[exp(0.5*L)] = exp(0.125)  (MGF of N(0,1) at t=0.5)
# True shift(1) vs shift(0) difference = 1.5
# True scale_by(2) vs natural course = 1.5 * E[A] = 1.5 * exp(0.125)

simulate_count_treatment <- function(n = 3000, seed = 42) {
  set.seed(seed)
  L <- rnorm(n)
  A <- rpois(n, exp(0.5 * L))
  Y <- 2 + 1.5 * A + L + rnorm(n)
  data.frame(Y = Y, A = A, L = L)
}

# Missing-data DGPs
#
# These DGPs produce data with explicit missingness mechanisms (MCAR, MAR)
# and known censoring models so the analytical truth is available under
# each mechanism.

# DGP-M1: MCAR outcome censoring.
#   Same SCM as DGP 1 (binary trt, continuous outcome, ATE = 3).
#   C ~ Bernoulli(p_cens), independent of everything.
#   Y set to NA when C = 1.
#
# Truth: ATE = 3 (complete-case is unbiased under MCAR).
simulate_mcar_outcome <- function(n = 2000, p_cens = 0.15, seed = 42) {
  set.seed(seed)
  L <- rnorm(n)
  A <- rbinom(n, 1, plogis(0.5 * L))
  Y_full <- 2 + 3 * A + 1.5 * L + rnorm(n)
  C <- rbinom(n, 1, p_cens)
  Y <- ifelse(C == 1, NA_real_, Y_full)
  data.frame(Y = Y, A = A, L = L, C = C)
}

# DGP-M2: MAR outcome censoring (informative).
#   Same SCM as DGP 1 for the structural model.
#   Censoring depends on A and L:
#     C | A, L ~ Bernoulli(expit(-2 + 0.8*A + 0.6*L))
#
# Under MAR censoring, complete-case analysis IS biased because the
# observed sample over-represents individuals with low A and low L.
# IPCW with the correct censoring model recovers the truth.
#
# Truth: ATE = 3 (the structural effect is unchanged; censoring only
# affects which rows are observed, not the SCM).
#
# Manual IPCW weights: w_i = 1 / P(C=0 | A_i, L_i)
#   = 1 / (1 - expit(-2 + 0.8*A_i + 0.6*L_i))
simulate_mar_outcome <- function(n = 5000, seed = 42) {
  set.seed(seed)
  L <- rnorm(n)
  A <- rbinom(n, 1, plogis(0.5 * L))
  Y_full <- 2 + 3 * A + 1.5 * L + rnorm(n)
  # Strong differential censoring: treated (A=1) and high-L individuals
  # are censored much more often. The intercept (-0.5) keeps overall
  # censoring at ~30-40% so the bias is detectable at n=10000.
  p_cens <- plogis(-0.5 + 1.5 * A + 1.0 * L)
  C <- rbinom(n, 1, p_cens)
  Y <- ifelse(C == 1, NA_real_, Y_full)
  data.frame(Y = Y, Y_full = Y_full, A = A, L = L, C = C, p_cens = p_cens)
}

# DGP-M2b: MAR outcome censoring with non-linear outcome model.
#   Two confounders (L1 continuous, L2 binary), interaction A*L1,
#   and a quadratic L1^2 term. Censoring depends on A, L1, L2.
#   A misspecified linear outcome model (Y ~ A + L1 + L2) will be
#   biased by ~0.25 under complete-case analysis but corrected by IPCW.
#
#   True SCM:
#     L1 ~ N(0,1), L2 ~ Bernoulli(0.4)
#     A ~ Bernoulli(expit(0.3*L1 - 0.5*L2 + 0.2*L1*L2))
#     Y = 1 + 3*A + 1.5*L1 + 2*L2 - 0.8*A*L1 + 0.5*L1^2 + eps
#     C ~ Bernoulli(expit(-0.3 + 1.2*A + 0.8*L1 - 0.5*L2))
#
#   True ATE: E[3 - 0.8*L1] = 3 (since E[L1] = 0).
#   ~50% overall censoring; ~65% among treated.
simulate_mar_outcome_complex <- function(n = 5000, seed = 42) {
  set.seed(seed)
  L1 <- rnorm(n)
  L2 <- rbinom(n, 1, 0.4)
  A <- rbinom(n, 1, plogis(0.3 * L1 - 0.5 * L2 + 0.2 * L1 * L2))
  Y_full <- 1 + 3 * A + 1.5 * L1 + 2 * L2 - 0.8 * A * L1 + 0.5 * L1^2 + rnorm(n)
  C <- rbinom(n, 1, plogis(-0.3 + 1.2 * A + 0.8 * L1 - 0.5 * L2))
  Y <- ifelse(C == 1, NA_real_, Y_full)
  data.frame(
    Y = Y,
    Y_full = Y_full,
    A = A,
    L1 = L1,
    L2 = L2,
    C = C
  )
}

# DGP-M3: MCAR covariate missingness.
#   Same SCM as DGP 1.
#   L_obs = L with some entries set to NA at random.
#   The FULL L is available for truth computation.
#
# Truth: ATE = 3.
simulate_mcar_covariate <- function(n = 2000, p_miss = 0.10, seed = 42) {
  set.seed(seed)
  L <- rnorm(n)
  A <- rbinom(n, 1, plogis(0.5 * L))
  Y <- 2 + 3 * A + 1.5 * L + rnorm(n)
  L_obs <- L
  miss_idx <- sample(n, size = floor(n * p_miss))
  L_obs[miss_idx] <- NA
  data.frame(Y = Y, A = A, L = L_obs, L_full = L)
}

# DGP-M4: Longitudinal MCAR outcome censoring.
#   Linear SCM (same as make_linear_scm), with random dropout at final time.
#   True ATE (always vs never, 2 periods) = 5.
simulate_longitudinal_mcar_outcome <- function(
  n = 3000,
  p_cens = 0.10,
  seed = 42
) {
  d <- make_linear_scm(n = n, n_times = 2, seed = seed)
  # Censor some individuals at the final time step
  set.seed(seed + 1L)
  ids_to_censor <- sample(n, size = floor(n * p_cens))
  final_mask <- d$time == 1 & d$id %in% ids_to_censor
  d$C <- 0L
  d$C[final_mask] <- 1L
  d$Y[final_mask] <- NA_real_
  d
}

# DGP-M5: Longitudinal MAR outcome censoring (informative dropout).
#   Linear SCM with dropout depending on treatment and baseline covariate.
#   C_1 | A_0, L_0 ~ Bernoulli(expit(-2 + A_0 + 0.5*L_0))
#   True ATE = 5 (structural, censoring does not change the SCM).
simulate_longitudinal_mar_outcome <- function(n = 5000, seed = 42) {
  d <- make_linear_scm(n = n, n_times = 2, seed = seed)
  # Informative censoring at t = 1 depends on A_0 and L_0.
  # Strong coefficients to produce detectable bias under complete-case.
  t0 <- d[d$time == 0, ]
  set.seed(seed + 1L)
  p_cens <- plogis(-0.5 + 1.5 * t0$A + 0.8 * t0$L0)
  C1 <- rbinom(n, 1, p_cens)
  d$C <- 0L
  d$C[d$time == 1 & d$id %in% t0$id[C1 == 1L]] <- 1L
  d$Y[d$C == 1L] <- NA_real_
  d$p_cens <- 0
  d$p_cens[d$time == 1] <- p_cens
  d
}

# DGP-EM-ICE: Longitudinal effect modification (binary trt x binary modifier).
#
# Linear SCM with treatment-confounder feedback and sex-specific treatment
# effects. The treatment effect at each period is (2 + 1.5 * sex) per unit
# of A_k, making the sex-specific contrasts analytically computable.
#
# DGP:
#   L0 ~ N(0, 1); sex ~ Bern(0.5)
#   L_0 = 0.5 * L0 + eps_L             (TV confounder at baseline)
#   A_0 ~ Bern(expit(0.5 * L0))
#   For t > 0:
#     L_t = A_{t-1} + 0.5 * L0 + eps_L  (treatment-confounder feedback)
#     A_t ~ Bern(expit(0.3 * L_t))
#   Y = 10 + (2 + 1.5*sex) * sum(A_t) + sum(L_t) + eps_Y
#
# True ATE (always vs never):
#   2 periods: ATE|sex=0 = 5,  ATE|sex=1 = 8
#   3 periods: ATE|sex=0 = 8,  ATE|sex=1 = 12.5
#
# Derivation (2-period, always vs never):
#   Under always: E[L_0]=0, E[L_1]=1,
#     E[Y|sex=s] = 10 + 2*(2+1.5s) + 0 + 1 = 15 + 3s
#   Under never:  E[L_0]=0, E[L_1]=0, E[Y|sex=s] = 10
#   ATE|sex=s = 5 + 3s => sex=0: 5, sex=1: 8
#
# Derivation (3-period, always vs never):
#   Under always: E[L_0]=0, E[L_1]=1, E[L_2]=1,
#     E[Y|sex=s] = 10+3*(2+1.5s)+0+1+1 = 18+4.5s
#   Under never:  E[L_0]=0, E[L_1]=0, E[L_2]=0, E[Y|sex=s] = 10
#   ATE|sex=s = 8 + 4.5s => sex=0: 8, sex=1: 12.5
make_em_ice_scm <- function(n = 5000, n_times = 2, seed = 42) {
  set.seed(seed)

  L0 <- stats::rnorm(n)
  sex <- stats::rbinom(n, 1, 0.5)

  A <- matrix(NA_real_, n, n_times)
  L <- matrix(NA_real_, n, n_times)

  for (t in seq_len(n_times)) {
    if (t == 1) {
      L[, t] <- 0.5 * L0 + stats::rnorm(n, 0, 0.5)
      A[, t] <- stats::rbinom(n, 1, stats::plogis(0.5 * L0))
    } else {
      L[, t] <- A[, t - 1] + 0.5 * L0 + stats::rnorm(n, 0, 0.5)
      A[, t] <- stats::rbinom(n, 1, stats::plogis(0.3 * L[, t]))
    }
  }

  # Treatment effect is (2 + 1.5*sex) per period of treatment.
  trt_effect <- (2 + 1.5 * sex) * rowSums(A)
  Y <- 10 + trt_effect + rowSums(L) + stats::rnorm(n)

  rows <- vector("list", n_times)
  for (t in seq_len(n_times)) {
    rows[[t]] <- data.frame(
      id = seq_len(n),
      time = t - 1L,
      A = A[, t],
      L = L[, t],
      L0 = L0,
      sex = sex,
      Y = if (t == n_times) Y else NA_real_
    )
  }
  do.call(rbind, rows)
}

# DGP-EM-ICE-CONT: Longitudinal EM with continuous treatment.
#
# Continuous-treatment version of the EM ICE DGP. Treatment effect at each
# period is (2 + 1.5 * sex) per unit of A_k, making the sex-specific
# shift contrasts analytically computable.
#
# DGP (2 periods):
#   L0 ~ N(0, 1); sex ~ Bern(0.5)
#   A_0 = 1 + 0.5 * L0 + eps_A
#   L_1 = A_0 + 0.5 * L0 + eps_L
#   A_1 = 1 + 0.3 * L_1 + 0.2 * A_0 + eps_A
#   Y = 10 + (2 + 1.5*sex) * (A_0 + A_1) + L_0 + L_1 + eps_Y
#
# True shift(delta) effect (g-formula, includes indirect path via L_1):
#   MC truth (n = 5*10^6): shift(1)|sex=0 ~ 6.0, shift(1)|sex=1 ~ 9.76
#   (larger than the direct coefficient effect 4/7 because L_1 = A_0 + ...)
make_em_ice_cont_scm <- function(n = 5000, seed = 42) {
  set.seed(seed)

  L0 <- stats::rnorm(n)
  sex <- stats::rbinom(n, 1, 0.5)

  Ltv0 <- 0.5 * L0 + stats::rnorm(n, 0, 0.5)
  A0 <- 1 + 0.5 * L0 + stats::rnorm(n, 0, 0.5)
  L1 <- A0 + 0.5 * L0 + stats::rnorm(n, 0, 0.5)
  A1 <- 1 + 0.3 * L1 + 0.2 * A0 + stats::rnorm(n, 0, 0.5)

  trt_effect <- (2 + 1.5 * sex) * (A0 + A1)
  Y <- 10 + trt_effect + Ltv0 + L1 + stats::rnorm(n)

  rbind(
    data.frame(
      id = seq_len(n),
      time = 0L,
      A = A0,
      L = Ltv0,
      L0 = L0,
      sex = sex,
      Y = NA_real_
    ),
    data.frame(
      id = seq_len(n),
      time = 1L,
      A = A1,
      L = L1,
      L0 = L0,
      sex = sex,
      Y = Y
    )
  )
}

# DGP-EM-ICE-BINOM: Longitudinal EM with binary outcome.
#
# Same structure as make_em_ice_scm but with a binary outcome via logistic
# link. Treatment effect is on the log-odds scale.
#
# DGP (2 periods):
#   L0 ~ N(0, 1); sex ~ Bern(0.5)
#   A_0 ~ Bern(expit(0.5 * L0))
#   L_1 = A_0 + 0.5 * L0 + eps_L
#   A_1 ~ Bern(expit(0.3 * L_1))
#   Y ~ Bern(expit(-1 + (1 + 0.8*sex)*(A_0 + A_1) + 0.5*L0 + 0.3*L_1))
#
# No closed-form truth on the probability scale due to nonlinear link,
# but MC truth (n = 10^6, seed = 1):
#   P[Y(always)|sex=0] ~ 0.775, P[Y(never)|sex=0] ~ 0.279
#   P[Y(always)|sex=1] ~ 0.942, P[Y(never)|sex=1] ~ 0.279
#   RD|sex=0 ~ 0.495, RD|sex=1 ~ 0.663
make_em_ice_binom_scm <- function(n = 5000, seed = 42) {
  set.seed(seed)

  L0 <- stats::rnorm(n)
  sex <- stats::rbinom(n, 1, 0.5)
  Ltv0 <- 0.5 * L0 + stats::rnorm(n, 0, 0.5)
  A0 <- stats::rbinom(n, 1, stats::plogis(0.5 * L0))
  L1 <- A0 + 0.5 * L0 + stats::rnorm(n, 0, 0.5)
  A1 <- stats::rbinom(n, 1, stats::plogis(0.3 * L1))

  eta <- -1 + (1 + 0.8 * sex) * (A0 + A1) + 0.5 * Ltv0 + 0.3 * L1
  Y <- stats::rbinom(n, 1, stats::plogis(eta))

  rbind(
    data.frame(
      id = seq_len(n),
      time = 0L,
      A = A0,
      L = Ltv0,
      L0 = L0,
      sex = sex,
      Y = NA_real_
    ),
    data.frame(
      id = seq_len(n),
      time = 1L,
      A = A1,
      L = L1,
      L0 = L0,
      sex = sex,
      Y = Y
    )
  )
}

# -- Stochastic intervention DGPs -------------------------------------------
#
# Each DGP returns a list with:
#   $data  — the observed data
#   $truth — the analytical E[Y^g] under the stochastic policy
#   $sampler — the sampler function for stochastic()

# Stochastic DGP 1: Binary treatment, Gaussian outcome
#   L ~ N(0, 1)
#   A | L ~ Bernoulli(expit(0.5 * L))
#   Y = 2 + 3 * A + 1.5 * L + \epsilon,  \epsilon ~ N(0, 1)
#   Policy: A* ~ Bernoulli(expit(0.5 + 0.3 * L))
#   Truth: E[Y^g] = 2 + 3 * E_L[expit(0.5 + 0.3 * L)]
#          (since E[L] = 0 under L ~ N(0,1))
simulate_stochastic_binary_gaussian <- function(n = 5000, seed = 42) {
  set.seed(seed)
  L <- rnorm(n)
  A <- rbinom(n, 1, plogis(0.5 * L))
  Y <- 2 + 3 * A + 1.5 * L + rnorm(n)
  data <- data.frame(Y = Y, A = A, L = L)

  # Analytical truth via numerical integration
  integrand <- function(l) plogis(0.5 + 0.3 * l) * dnorm(l)
  e_pg <- integrate(integrand, lower = -10, upper = 10)$value
  truth <- 2 + 3 * e_pg

  sampler <- function(data, trt) {
    rbinom(nrow(data), 1, plogis(0.5 + 0.3 * data$L))
  }

  list(data = data, truth = truth, sampler = sampler)
}

# Stochastic DGP 2: Continuous treatment, Gaussian outcome
#   L ~ N(0, 1)
#   A | L ~ N(L, 1)
#   Y = 2 + A + 1.5 * L + \epsilon,  \epsilon ~ N(0, 1)
#   Policy: A* = A + \epsilon_{shift}, \epsilon_{shift} ~ N(0.5, 0.25)
#   Truth: E[Y^g] = 2 + E[A] + 0.5 + 1.5 * E[L]
#        = 2 + E[L] + 0.5 + 1.5 * E[L] = 2 + 0 + 0.5 + 0 = 2.5
simulate_stochastic_continuous_gaussian <- function(n = 5000, seed = 42) {
  set.seed(seed)
  L <- rnorm(n)
  A <- L + rnorm(n)
  Y <- 2 + A + 1.5 * L + rnorm(n)
  data <- data.frame(Y = Y, A = A, L = L)

  truth <- 2.5

  sampler <- function(data, trt) {
    trt + rnorm(length(trt), mean = 0.5, sd = 0.5)
  }

  list(data = data, truth = truth, sampler = sampler)
}

# Stochastic DGP 3: Binary treatment, Binomial outcome
#   L ~ N(0, 1)
#   A | L ~ Bernoulli(expit(0.5 * L))
#   Y | A, L ~ Bernoulli(expit(-1 + 2 * A + L))
#   Policy: A* ~ Bernoulli(expit(0.5 + 0.3 * L))
#   Truth: E[Y^g] = E_L[ p_g(L) * expit(1 + L) + (1 - p_g(L)) * expit(-1 + L) ]
#   where p_g(L) = expit(0.5 + 0.3 * L)
simulate_stochastic_binary_binomial <- function(n = 5000, seed = 42) {
  set.seed(seed)
  L <- rnorm(n)
  A <- rbinom(n, 1, plogis(0.5 * L))
  Y <- rbinom(n, 1, plogis(-1 + 2 * A + L))
  data <- data.frame(Y = Y, A = A, L = L)

  integrand <- function(l) {
    pg <- plogis(0.5 + 0.3 * l)
    ey <- pg * plogis(-1 + 2 + l) + (1 - pg) * plogis(-1 + l)
    ey * dnorm(l)
  }
  truth <- integrate(integrand, lower = -10, upper = 10)$value

  sampler <- function(data, trt) {
    rbinom(nrow(data), 1, plogis(0.5 + 0.3 * data$L))
  }

  list(data = data, truth = truth, sampler = sampler)
}

# Stochastic DGP 4: Categorical treatment (k=3), Gaussian outcome
#   L ~ N(0, 1)
#   A | L ~ Categorical(softmax(0, 0.5*L, -0.3*L)), levels = c("a","b","c")
#   Y = 1 + 2 * I(A=="b") + 3 * I(A=="c") + L + \epsilon
#   Policy: A* ~ Uniform("a", "b", "c")
#   Truth: E[Y^g] = 1 + 2 * (1/3) + 3 * (1/3) + E[L]
#        = 1 + 2/3 + 1 + 0 = 8/3
simulate_stochastic_categorical_gaussian <- function(n = 5000, seed = 42) {
  set.seed(seed)
  L <- rnorm(n)
  logits <- cbind(0, 0.5 * L, -0.3 * L)
  probs <- exp(logits) / rowSums(exp(logits))
  A_idx <- apply(probs, 1, function(p) sample(1:3, 1, prob = p))
  A <- factor(c("a", "b", "c")[A_idx], levels = c("a", "b", "c"))
  Y <- 1 + 2 * (A == "b") + 3 * (A == "c") + L + rnorm(n)
  data <- data.frame(Y = Y, A = A, L = L)

  truth <- 1 + 2 / 3 + 3 / 3

  sampler <- function(data, trt) {
    factor(
      sample(c("a", "b", "c"), nrow(data), replace = TRUE),
      levels = levels(trt)
    )
  }

  list(data = data, truth = truth, sampler = sampler)
}

# Stochastic DGP 5: Multivariate treatment (binary + continuous), Gaussian
#   L ~ N(0, 1)
#   A1 | L ~ Bernoulli(expit(L))
#   A2 | A1, L ~ N(L + A1, 1)
#   Y = 1 + A1 + 0.5 * A2 + L + \epsilon
#   Policy: A1* ~ Bernoulli(0.7), A2* = A2 + \epsilon, \epsilon ~ N(0.5, 0.25)
#   Truth: E[Y^g] = 1 + 0.7 + 0.5 * (E[A2] + 0.5) + E[L]
#   E[A2] = E[L + A1] = E[L] + E[A1] = 0 + E_L[expit(L)]
#   E_L[expit(L)] with L ~ N(0,1) = 0.5 (symmetry of expit around 0)
#   E[A2] = 0.5, E[A2*] = 0.5 + 0.5 = 1.0
#   Truth = 1 + 0.7 + 0.5 * 1.0 + 0 = 2.2
simulate_stochastic_multivariate_gaussian <- function(n = 5000, seed = 42) {
  set.seed(seed)
  L <- rnorm(n)
  A1 <- rbinom(n, 1, plogis(L))
  A2 <- L + A1 + rnorm(n)
  Y <- 1 + A1 + 0.5 * A2 + L + rnorm(n)
  data <- data.frame(Y = Y, A1 = A1, A2 = A2, L = L)

  # Exact E_L[expit(L)] via integration
  e_expit_l <- integrate(
    function(l) plogis(l) * dnorm(l),
    lower = -10,
    upper = 10
  )$value
  e_a2 <- e_expit_l
  truth <- 1 + 0.7 + 0.5 * (e_a2 + 0.5)

  sampler_a1 <- function(data, trt) {
    rbinom(nrow(data), 1, 0.7)
  }
  sampler_a2 <- function(data, trt) {
    trt + rnorm(length(trt), mean = 0.5, sd = 0.5)
  }

  list(
    data = data,
    truth = truth,
    sampler_a1 = sampler_a1,
    sampler_a2 = sampler_a2
  )
}

# -- Longitudinal stochastic DGPs -------------------------------------------

# Longitudinal stochastic DGP 1: Binary treatment, 2-period, Gaussian outcome
#   L_0 ~ N(0, 1)
#   A_0 ~ Bernoulli(expit(0.5 * L_0))
#   L_1 = 0.5 * A_0 + 0.5 * L_0 + \epsilon_L
#   A_1 ~ Bernoulli(expit(0.3 * L_1))
#   Y = 1 + A_0 + A_1 + 0.5 * L_0 + 0.5 * L_1 + \epsilon_Y
#   Policy: A*_k ~ Bernoulli(expit(0.2 + 0.3 * L_k))
#
# Truth via large-n simulation oracle (n = 10^6) since analytical
# nested integration is complex with treatment-confounder feedback.
simulate_stochastic_ice_binary_gaussian <- function(
  n = 5000,
  seed = 42
) {
  set.seed(seed)
  L0 <- rnorm(n)
  Ltv0 <- 0.5 * L0 + rnorm(n, 0, 0.5)
  A0 <- rbinom(n, 1, plogis(0.5 * L0))
  L1 <- 0.5 * A0 + 0.5 * L0 + rnorm(n, 0, 0.5)
  A1 <- rbinom(n, 1, plogis(0.3 * L1))
  Y <- 1 + A0 + A1 + 0.5 * Ltv0 + 0.5 * L1 + rnorm(n)

  data <- rbind(
    data.frame(
      id = seq_len(n),
      time = 0L,
      A = A0,
      L = Ltv0,
      L0 = L0,
      Y = NA_real_
    ),
    data.frame(id = seq_len(n), time = 1L, A = A1, L = L1, L0 = L0, Y = Y)
  )

  # Monte Carlo truth via large-n simulation (n = 10^6)
  set.seed(1)
  n_mc_truth <- 1e6
  L0_mc <- rnorm(n_mc_truth)
  Ltv0_mc <- 0.5 * L0_mc + rnorm(n_mc_truth, 0, 0.5)
  A0_mc <- rbinom(n_mc_truth, 1, plogis(0.2 + 0.3 * Ltv0_mc))
  L1_mc <- 0.5 * A0_mc + 0.5 * L0_mc + rnorm(n_mc_truth, 0, 0.5)
  A1_mc <- rbinom(n_mc_truth, 1, plogis(0.2 + 0.3 * L1_mc))
  Y_mc <- 1 + A0_mc + A1_mc + 0.5 * Ltv0_mc + 0.5 * L1_mc
  truth <- mean(Y_mc)

  sampler <- function(data, trt) {
    rbinom(nrow(data), 1, plogis(0.2 + 0.3 * data$L))
  }

  list(data = data, truth = truth, sampler = sampler)
}

# Longitudinal stochastic DGP 2: Continuous treatment, 2-period, Gaussian
#   L_0 ~ N(0, 1)
#   A_0 ~ N(L_0, 1)
#   L_1 = 0.5 * A_0 + 0.5 * L_0 + \epsilon_L
#   A_1 ~ N(L_1, 1)
#   Y = 1 + 0.5 * A_0 + 0.5 * A_1 + 0.5 * L_0 + 0.5 * L_1 + \epsilon_Y
#   Policy: A*_k = A_k + 0.5 (deterministic shift encoded as stochastic)
#   Truth via simulation oracle.
simulate_stochastic_ice_continuous_gaussian <- function(
  n = 5000,
  seed = 42
) {
  set.seed(seed)
  L0 <- rnorm(n)
  Ltv0 <- 0.5 * L0 + rnorm(n, 0, 0.5)
  A0 <- L0 + rnorm(n)
  L1 <- 0.5 * A0 + 0.5 * L0 + rnorm(n, 0, 0.5)
  A1 <- L1 + rnorm(n)
  Y <- 1 + 0.5 * A0 + 0.5 * A1 + 0.5 * Ltv0 + 0.5 * L1 + rnorm(n)

  data <- rbind(
    data.frame(
      id = seq_len(n),
      time = 0L,
      A = A0,
      L = Ltv0,
      L0 = L0,
      Y = NA_real_
    ),
    data.frame(id = seq_len(n), time = 1L, A = A1, L = L1, L0 = L0, Y = Y)
  )

  # Additive shift of +0.5 at each time point.
  # Use simulation oracle.
  set.seed(1)
  n_mc_truth <- 1e6
  L0_mc <- rnorm(n_mc_truth)
  Ltv0_mc <- 0.5 * L0_mc + rnorm(n_mc_truth, 0, 0.5)
  A0_mc <- L0_mc + rnorm(n_mc_truth) + 0.5
  L1_mc <- 0.5 * A0_mc + 0.5 * L0_mc + rnorm(n_mc_truth, 0, 0.5)
  A1_mc <- L1_mc + rnorm(n_mc_truth) + 0.5
  Y_mc <- 1 + 0.5 * A0_mc + 0.5 * A1_mc + 0.5 * Ltv0_mc + 0.5 * L1_mc
  truth <- mean(Y_mc)

  sampler <- function(data, trt) {
    trt + 0.5
  }

  list(data = data, truth = truth, sampler = sampler)
}


# DGP-T1: Binary treatment, interaction outcome, transportability
#   L ~ N(0, 1)
#   P(S = 1 | L) = expit(-0.5 + 1.0 * L)   [study under-represents L < 0]
#   A | L, S = 1 ~ Bernoulli(expit(0.2 + 0.3 * L))
#   Y | A, L ~ N(2 + 3 * A + 1.5 * L + 1.0 * A * L, 1)
#
# Target ATE (over L ~ N(0,1)): E[Y^1 - Y^0] = 3 + 1.0 * E[L] = 3
# Study ATE (over L | S = 1):  3 + 1.0 * E[L | S = 1] != 3
#   (S = 1 oversamples high-L individuals via the selection mechanism)
#
# The A*L interaction makes the study and target ATEs diverge, which is
# necessary for truth-based tests in chunks 17b/17c.
simulate_transport <- function(n = 3000, seed = 42) {
  set.seed(seed)
  L <- rnorm(n)
  ps_sampling <- plogis(-0.5 + 1.0 * L)
  S <- rbinom(n, 1, ps_sampling)

  ps_treat <- plogis(0.2 + 0.3 * L)
  A <- ifelse(S == 1L, rbinom(n, 1, ps_treat), NA_integer_)

  Y <- ifelse(
    S == 1L,
    2 + 3 * A + 1.5 * L + 1.0 * A * L + rnorm(n),
    NA_real_
  )

  data.frame(Y = Y, A = A, L = L, S = S)
}


# DGP for multivariate IPW × transport (chunk 17j).
#
# Y = 1 + 2*A1 + 1.5*A2 + L + A1*L + eps,  eps ~ N(0,1)
# A1 | L, S=1 ~ Bernoulli(expit(0.2 + 0.3*L))
# A2 | L, S=1 ~ Bernoulli(expit(-0.1 + 0.2*L))
# P(S=1 | L) = expit(-0.5 + L)
#
# Target ATE [(A1=1,A2=1) vs (A1=0,A2=0)]:
#   Y^{1,1} - Y^{0,0} = (1+2+1.5+L+L) - (1+L) = 3.5 + L
#   => E_target[ATE] = 3.5 + E_target[L]
#
# Generalizability (target_subset="all"): E[L] = 0, so ATE = 3.5.
# Transportability (target_subset="target"): E[L|S=0] < 0, so ATE < 3.5.
# Study ATE: E[L|S=1] > 0, so naive study estimator is upward-biased.
simulate_mv_transport <- function(n = 6000, seed = 42) {
  set.seed(seed)
  L <- rnorm(n)
  ps_sampling <- plogis(-0.5 + L)
  S <- rbinom(n, 1, ps_sampling)

  ps_a1 <- plogis(0.2 + 0.3 * L)
  ps_a2 <- plogis(-0.1 + 0.2 * L)
  A1 <- ifelse(S == 1L, rbinom(n, 1, ps_a1), NA_integer_)
  A2 <- ifelse(S == 1L, rbinom(n, 1, ps_a2), NA_integer_)

  Y <- ifelse(
    S == 1L,
    1 + 2 * A1 + 1.5 * A2 + L + A1 * L + rnorm(n),
    NA_real_
  )

  data.frame(Y = Y, A1 = A1, A2 = A2, L = L, S = S)
}


# DGP for IPCW × transport (chunk 17k).
#
# Combines differential MAR censoring with sampling selection.
# Study rows (S=1) have treatment, outcome (possibly censored), and covariates.
# Target rows (S=0) have covariates only (A=NA, Y=NA).
#
#   L ~ N(0, 1)
#   P(S = 1 | L) = expit(-0.5 + L)
#   A | L, S=1 ~ Bernoulli(expit(0.2 + 0.3*L))
#   C | A, L, S=1 ~ Bernoulli(expit(-1.5 + 0.5*A + 0.3*L))   [~20% censoring]
#   Y | A, L, S=1, C=0 ~ N(2 + 3*A + 1.5*L + A*L, 1)
#
# True target ATE (transportability): 3 + E[L|S=0]
# True target ATE (generalizability): 3 + E[L] = 3
# Study ATE: 3 + E[L|S=1] (biased by sampling)
#
# The interaction A*L makes the ATE depend on the covariate distribution,
# so the naive study estimator differs from the target estimand.
# The differential censoring (C depends on A and L) biases complete-case
# analysis; IPCW corrects for this.
simulate_ipcw_transport <- function(n = 6000, seed = 42) {
  set.seed(seed)
  L <- rnorm(n)
  ps_sampling <- plogis(-0.5 + 1.0 * L)
  S <- rbinom(n, 1, ps_sampling)

  ps_treat <- plogis(0.2 + 0.3 * L)
  A <- ifelse(S == 1L, rbinom(n, 1, ps_treat), NA_integer_)

  Y_full <- ifelse(
    S == 1L,
    2 + 3 * A + 1.5 * L + 1.0 * A * L + rnorm(n),
    NA_real_
  )

  # Differential censoring: treated + high-L individuals censored more.
  # Intercept -1.5 keeps overall censoring at ~20% so study rows remain usable.
  p_cens <- ifelse(S == 1L, plogis(-1.5 + 0.5 * A + 0.3 * L), 0)
  C <- ifelse(S == 1L, rbinom(n, 1, p_cens), 0L)
  Y <- ifelse(C == 1L, NA_real_, Y_full)

  data.frame(Y = Y, Y_full = Y_full, A = A, L = L, S = S, C = C)
}


# DGP for MTP + transport tests (continuous treatment).
#
#   L ~ N(0, 1)
#   P(S = 1 | L) = expit(-0.5 + L)
#   A | L, S=1 ~ N(0.5 + 0.3L, 1)   (continuous treatment on study rows)
#   Y | A, L, S=1 ~ N(2 + 3A + 1.5L + AL, 1)
#
# Truth for shift(delta) vs natural course on the target population:
#   E[Y^{A+delta}] - E[Y^A] = 3*delta + delta*E[L|target]
#
# The shift effect depends on L through the A*L interaction term.
# Target rows (S=0) have A = NA (not observed).
simulate_mtp_transport <- function(n = 5000, seed = 42) {
  set.seed(seed)
  L <- rnorm(n)
  p_study <- plogis(-0.5 + 1.0 * L)
  S <- rbinom(n, 1, p_study)

  A <- ifelse(S == 1L, rnorm(n, 0.5 + 0.3 * L, 1), NA_real_)
  Y <- ifelse(
    S == 1L,
    2 + 3 * A + 1.5 * L + A * L + rnorm(n),
    NA_real_
  )

  data.frame(Y = Y, A = A, L = L, S = S)
}


# Kang & Schafer (2007) adversarial DGP for double-robustness testing.
#
# True covariates Z1..Z4 ~ iid N(0, 1).
# Propensity: \pi(Z) = expit(-Z1 + 0.5 Z2 - 0.25 Z3 - 0.1 Z4).
# Treatment: A ~ Bernoulli(\pi(Z)).
# Outcome (does not depend on A):
#   Y = 210 + 27.4 Z1 + 13.7 Z2 + 13.7 Z3 + 13.7 Z4 + \epsilon.
#
# The analyst's observable (misspecified) covariates:
#   X1 = exp(Z1 / 2)
#   X2 = Z2 / (1 + exp(Z1)) + 10
#   X3 = (Z1 Z3 / 25 + 0.6)^3
#   X4 = (Z2 + Z4 + 20)^2
#
# True population mean E[Y] = 210 (since E[Z_j] = 0).
# ATE = 0 (Y does not depend on A).
#
# Four scenarios:
#   S1: outcome ~ Z, PS ~ Z         => all methods recover 210.
#   S2: outcome ~ X, PS ~ Z         => g-comp biased; IPW/AIPW recover 210.
#   S3: outcome ~ Z, PS ~ X         => IPW biased; g-comp/AIPW recover 210.
#   S4: outcome ~ X, PS ~ X         => all biased (negative control).
#
# Returns a list:
#   $data — data.frame with Y, A, Z1..Z4, X1..X4.
#   $truth — 210.
simulate_kang_schafer <- function(n = 5000, seed = 42) {
  set.seed(seed)

  Z1 <- rnorm(n)
  Z2 <- rnorm(n)
  Z3 <- rnorm(n)
  Z4 <- rnorm(n)

  ps <- plogis(-Z1 + 0.5 * Z2 - 0.25 * Z3 - 0.1 * Z4)
  A <- rbinom(n, 1, ps)

  Y <- 210 + 27.4 * Z1 + 13.7 * Z2 + 13.7 * Z3 + 13.7 * Z4 + rnorm(n)

  # Misspecified covariates (nonlinear transformations of Z).
  X1 <- exp(Z1 / 2)
  X2 <- Z2 / (1 + exp(Z1)) + 10
  X3 <- (Z1 * Z3 / 25 + 0.6)^3
  X4 <- (Z2 + Z4 + 20)^2

  list(
    data = data.frame(
      Y = Y,
      A = A,
      Z1 = Z1,
      Z2 = Z2,
      Z3 = Z3,
      Z4 = Z4,
      X1 = X1,
      X2 = X2,
      X3 = X3,
      X4 = X4
    ),
    truth = 210
  )
}


# Naimi, Cole & Kennedy (2017)-inspired longitudinal DGP.
#
# Two time points, binary treatment, binary time-varying confounder,
# continuous outcome (CD4 count). Treatment-confounder feedback at t=1.
# A continuous baseline confounder (age) makes the propensity model
# well-identified at both time points.
#
# DGP:
#   age ~ N(40, 10)                              (baseline, continuous)
#   L0 ~ Bernoulli(expit(-1 + 0.02 * age))       (TV confounder at t=0)
#   A0 ~ Bernoulli(expit(-1 + L0 + 0.01 * age))  (ART at t=0)
#   L1 ~ Bernoulli(expit(-1 + 2*A0 + 3*L0))      (TV confounder at t=1)
#   A1 ~ Bernoulli(expit(-1 + 2*A0 + L1))         (ART at t=1)
#   Y = 200*A0 + 200*A1 + 500*L0 + 100*L1 + 2*age + epsilon
#
# MC truth (n = 10^6, seed = 1): computed below.
simulate_naimi_longitudinal <- function(n = 5000, seed = 42) {
  set.seed(seed)

  age <- rnorm(n, 40, 10)
  L0 <- rbinom(n, 1, plogis(-1 + 0.02 * age))
  A0 <- rbinom(n, 1, plogis(-1 + L0 + 0.01 * age))
  L1 <- rbinom(n, 1, plogis(-1 + 2 * A0 + 3 * L0))
  A1 <- rbinom(n, 1, plogis(-1 + 2 * A0 + L1))
  Y <- 200 * A0 + 200 * A1 + 500 * L0 + 100 * L1 + 2 * age + rnorm(n, 0, 50)

  data <- rbind(
    data.frame(
      id = seq_len(n),
      time = 0L,
      A = A0,
      Ltv = as.numeric(L0),
      age = age,
      Y = NA_real_
    ),
    data.frame(
      id = seq_len(n),
      time = 1L,
      A = A1,
      Ltv = as.numeric(L1),
      age = age,
      Y = Y
    )
  )

  # MC truth via large-n simulation.
  set.seed(1)
  n_mc <- 1e6
  age_mc <- rnorm(n_mc, 40, 10)
  L0_mc <- rbinom(n_mc, 1, plogis(-1 + 0.02 * age_mc))

  # Always treated: A0=1, A1=1.
  L1_always <- rbinom(n_mc, 1, plogis(-1 + 2 * 1 + 3 * L0_mc))
  Y_always <- 200 + 200 + 500 * L0_mc + 100 * L1_always + 2 * age_mc

  # Never treated: A0=0, A1=0.
  L1_never <- rbinom(n_mc, 1, plogis(-1 + 2 * 0 + 3 * L0_mc))
  Y_never <- 500 * L0_mc + 100 * L1_never + 2 * age_mc

  truth_always <- mean(Y_always)
  truth_never <- mean(Y_never)
  truth_ate <- truth_always - truth_never

  list(
    data = data,
    truth_ate = truth_ate,
    truth_always = truth_always,
    truth_never = truth_never
  )
}


# --- SNM DGPs (Phase 18) ---------------------------------------------------

#' Point-treatment SNM DGP with effect modification
#'
#' From PHASE_18_SNMS.md design doc:
#'   L ~ N(0, 1), M = I(L > 0)
#'   A | L ~ N(0.5 * L, 1)   (continuous treatment)
#'   Y = 2 + 3*A + 1.5*L + 2*A*M + eps, eps ~ N(0, 1)
#'
#' Linear blip: gamma(a, l; psi) = a * (psi_0 + psi_M * m)
#' Truth: psi_0 = 3, psi_M = 2
#'
#' @param n Sample size.
#' @param seed RNG seed.
#' @return List with `data`, `truth_psi` (named vector).
#' @noRd
simulate_snm_point <- function(n = 2000, seed = 42) {
  set.seed(seed)
  L <- stats::rnorm(n)
  M <- as.numeric(L > 0)
  A <- stats::rnorm(n, mean = 0.5 * L, sd = 1)
  eps <- stats::rnorm(n)
  Y <- 2 + 3 * A + 1.5 * L + 2 * A * M + eps

  data <- data.table::data.table(Y = Y, A = A, L = L, M = M)
  list(
    data = data,
    truth_psi = c(psi_intercept = 3, psi_M = 2)
  )
}


#' Point-treatment SNM DGP without effect modification (constant ATE)
#'
#'   L ~ N(0, 1)
#'   A | L ~ N(0.5 * L, 1)
#'   Y = 2 + 3*A + 1.5*L + eps, eps ~ N(0, 1)
#'
#' Linear blip: gamma(a; psi) = a * psi_0
#' Truth: psi_0 = 3
#'
#' @param n Sample size.
#' @param seed RNG seed.
#' @return List with `data`, `truth_psi`.
#' @noRd
simulate_snm_point_no_em <- function(n = 2000, seed = 42) {
  set.seed(seed)
  L <- stats::rnorm(n)
  A <- stats::rnorm(n, mean = 0.5 * L, sd = 1)
  eps <- stats::rnorm(n)
  Y <- 2 + 3 * A + 1.5 * L + eps

  data <- data.table::data.table(Y = Y, A = A, L = L)
  list(
    data = data,
    truth_psi = c(psi_intercept = 3)
  )
}


#' Point-treatment SNM DGP with binary treatment and effect modification
#'
#'   L ~ N(0, 1), M = I(L > 0)
#'   A | L ~ Bernoulli(expit(0.5 * L))
#'   Y = 2 + 3*A + 1.5*L + 2*A*M + eps, eps ~ N(0, 1)
#'
#' Truth: psi_0 = 3, psi_M = 2
#'
#' @param n Sample size.
#' @param seed RNG seed.
#' @return List with `data`, `truth_psi`.
#' @noRd
simulate_snm_point_binary <- function(n = 2000, seed = 42) {
  set.seed(seed)
  L <- stats::rnorm(n)
  M <- as.numeric(L > 0)
  A <- stats::rbinom(n, 1, stats::plogis(0.5 * L))
  eps <- stats::rnorm(n)
  Y <- 2 + 3 * A + 1.5 * L + 2 * A * M + eps

  data <- data.table::data.table(Y = Y, A = A, L = L, M = M)
  list(
    data = data,
    truth_psi = c(psi_intercept = 3, psi_M = 2)
  )
}


#' Point-treatment SNM DGP with time-varying (post-treatment) modifier
#'
#' The modifier M depends on treatment A, making it genuinely
#' post-treatment. SNMs identify the blip under treatment-model
#' correctness alone — M can be post-treatment. IPW-MSM with M as
#' a modifier would condition on a descendant of A (biased).
#'
#'   L ~ N(0, 1)
#'   A | L ~ N(0.5 * L, 1)
#'   M = 0.3 * A + 0.5 * L + eps_M, eps_M ~ N(0, 0.5)
#'   Y = 2 + 3*A + 1.5*L + 2*A*M + eps_Y, eps_Y ~ N(0, 1)
#'
#' Linear blip: gamma(a, m; psi) = a * (psi_0 + psi_M * m)
#' With treatment-free model ~L: psi_0 = 3, psi_M = 2.
#' Without treatment-free model, moment-condition estimates differ
#' because the blip absorbs the A -> M -> Y indirect path.
#'
#' @param n Sample size.
#' @param seed RNG seed.
#' @return List with `data`, `truth_psi_with_tf` (under TF model).
#' @noRd
simulate_snm_point_tv_modifier <- function(n = 2000, seed = 42) {
  set.seed(seed)
  L <- stats::rnorm(n)
  A <- stats::rnorm(n, mean = 0.5 * L, sd = 1)
  M <- 0.3 * A + 0.5 * L + stats::rnorm(n, sd = sqrt(0.5))
  eps <- stats::rnorm(n)
  Y <- 2 + 3 * A + 1.5 * L + 2 * A * M + eps

  data <- data.table::data.table(Y = Y, A = A, L = L, M = M)
  list(
    data = data,
    truth_psi_with_tf = c(psi_intercept = 3, psi_M = 2)
  )
}


#' Point-treatment SNM DGP with binary treatment and TV modifier
#'
#'   L ~ N(0, 1)
#'   A | L ~ Bernoulli(expit(0.5 * L))
#'   M = 0.5 * A + 0.5 * L + eps_M, eps_M ~ N(0, 0.5)
#'   Y = 2 + 3*A + 1.5*L + 2*A*M + eps_Y, eps_Y ~ N(0, 1)
#'
#' With treatment-free model ~L: psi_0 = 3, psi_M = 2.
#'
#' @param n Sample size.
#' @param seed RNG seed.
#' @return List with `data`, `truth_psi_with_tf`.
#' @noRd
simulate_snm_point_tv_modifier_binary <- function(n = 2000, seed = 42) {
  set.seed(seed)
  L <- stats::rnorm(n)
  A <- stats::rbinom(n, 1, stats::plogis(0.5 * L))
  M <- 0.5 * A + 0.5 * L + stats::rnorm(n, sd = sqrt(0.5))
  eps <- stats::rnorm(n)
  Y <- 2 + 3 * A + 1.5 * L + 2 * A * M + eps

  data <- data.table::data.table(Y = Y, A = A, L = L, M = M)
  list(
    data = data,
    truth_psi_with_tf = c(psi_intercept = 3, psi_M = 2)
  )
}


#' Longitudinal SNM DGP: 2-period, binary treatment, no EM
#'
#' Per-stage backward-sequential estimation (Robins 1994). Each
#' stage has its own blip parameter estimated via backward induction.
#'
#'   L_0 ~ N(0, 1)
#'   A_0 | L_0 ~ Bernoulli(expit(0.3 * L_0))
#'   L_1 = 0.5 * L_0 + 0.3 * A_0 + eps_L,  eps_L ~ N(0, sqrt(0.5))
#'   A_1 | L_1, A_0 ~ Bernoulli(expit(0.3 * L_1 + 0.2 * A_0))
#'   Y = 2 + 3 * A_0 + 3 * A_1 + 1.5 * L_0 + 0.5 * L_1 + eps_Y
#'
#' The stage-0 blip captures the total causal effect of \eqn{A_0} on
#' \eqn{Y} holding \eqn{A_1} fixed, which includes the mediated path
#' \eqn{A_0 \to L_1 \to Y} (coefficient \eqn{0.3 \times 0.5 = 0.15}).
#' DTRreg confirms convergence to these values at n = 500k.
#'
#' Truth (per-stage): stage0_psi_intercept = 3.15,
#'   stage1_psi_intercept = 3
#'
#' @param n Number of individuals.
#' @param seed RNG seed.
#' @return List with `data`, `truth_psi` (named vector of per-stage
#'   blip parameters).
#' @noRd
simulate_snm_longitudinal <- function(n = 2000, seed = 42) {
  set.seed(seed)

  L0 <- stats::rnorm(n)
  A0 <- stats::rbinom(n, 1, stats::plogis(0.3 * L0))

  # A_0 -> L_1 coefficient = 0.3; L_1 -> Y coefficient = 0.5
  # So the mediated A_0 -> L_1 -> Y path contributes 0.3 * 0.5 = 0.15
  L1 <- 0.5 * L0 + 0.3 * A0 + stats::rnorm(n, sd = sqrt(0.5))
  A1 <- stats::rbinom(n, 1, stats::plogis(0.3 * L1 + 0.2 * A0))

  Y <- 2 + 3 * A0 + 3 * A1 + 1.5 * L0 + 0.5 * L1 + stats::rnorm(n)

  # Interleave per-individual rows: (time 0, time 1) pairs
  data <- data.table::data.table(
    id = rep(seq_len(n), each = 2L),
    time = rep(0:1, n),
    Y = as.numeric(rbind(NA, Y)),
    A = as.numeric(rbind(A0, A1)),
    L = as.numeric(rbind(L0, L1))
  )

  # Stage-1 blip = 3 (direct A_1 -> Y)
  # Stage-0 blip = 3 + 0.3 * 0.5 = 3.15 (direct + mediated A_0 -> L_1 -> Y)
  list(
    data = data,
    truth_psi = c(
      stage0_psi_intercept = 3.15,
      stage1_psi_intercept = 3
    )
  )
}


#' Longitudinal SNM DGP: 2-period, continuous treatment, no EM
#'
#' Same structure as `simulate_snm_longitudinal()` but with Gaussian
#' treatment models. Same mediation channel applies: the stage-0 blip
#' includes the \eqn{A_0 \to L_1 \to Y} path.
#'
#'   L_0 ~ N(0, 1)
#'   A_0 | L_0 ~ N(0.5 * L_0, 1)
#'   L_1 = 0.5 * L_0 + 0.3 * A_0 + eps_L
#'   A_1 | L_1, A_0 ~ N(0.3 * L_1 + 0.2 * A_0, 1)
#'   Y = 2 + 3 * A_0 + 3 * A_1 + 1.5 * L_0 + 0.5 * L_1 + eps_Y
#'
#' Truth (per-stage): stage0_psi_intercept = 3.15,
#'   stage1_psi_intercept = 3
#'
#' @param n Number of individuals.
#' @param seed RNG seed.
#' @return List with `data`, `truth_psi`.
#' @noRd
simulate_snm_longitudinal_continuous <- function(n = 2000, seed = 42) {
  set.seed(seed)

  L0 <- stats::rnorm(n)
  A0 <- stats::rnorm(n, mean = 0.5 * L0, sd = 1)

  L1 <- 0.5 * L0 + 0.3 * A0 + stats::rnorm(n, sd = sqrt(0.5))
  A1 <- stats::rnorm(n, mean = 0.3 * L1 + 0.2 * A0, sd = 1)

  Y <- 2 + 3 * A0 + 3 * A1 + 1.5 * L0 + 0.5 * L1 + stats::rnorm(n)

  # Interleave per-individual rows: (time 0, time 1) pairs
  data <- data.table::data.table(
    id = rep(seq_len(n), each = 2L),
    time = rep(0:1, n),
    Y = as.numeric(rbind(NA, Y)),
    A = as.numeric(rbind(A0, A1)),
    L = as.numeric(rbind(L0, L1))
  )

  list(
    data = data,
    truth_psi = c(
      stage0_psi_intercept = 3.15,
      stage1_psi_intercept = 3
    )
  )
}


#' Longitudinal SNM DGP: 2-period, binary treatment, time-varying EM
#'
#' Headline Phase 18 demonstration: the stage-1 blip includes a
#' time-varying modifier \eqn{M_1 = 1\{L_1 > 0\}} that is
#' post-treatment (L_1 depends on A_0). SNMs handle this correctly;
#' IPW-MSM conditioning on M_1 introduces collider bias.
#'
#' M_0 = 1{L_0 > 0} is baseline (not affected by treatment).
#' Both stages share the same blip specification (intercept + M),
#' so the `build_blip_design_matrix()` applies uniformly.
#'
#'   L_0 ~ N(0, 1),  M_0 = 1{L_0 > 0}
#'   A_0 | L_0 ~ Bernoulli(expit(0.3 * L_0))
#'   L_1 = 0.5 * L_0 + 0.3 * A_0 + eps_L,  eps_L ~ N(0, sqrt(0.5))
#'   M_1 = 1{L_1 > 0}                       (post-treatment!)
#'   A_1 | L_1, A_0 ~ Bernoulli(expit(0.3 * L_1 + 0.2 * A_0))
#'   Y = 2 + 1*A_0 + 2*A_0*M_0 + 2*A_1 + 2*A_1*M_1
#'       + 1.5*L_0 + 0.5*L_1 + eps_Y
#'
#' Stage-0 blip: gamma_0 = a_0 * (psi_00 + psi_0M * M_0)
#'   psi_00 = 1 + 0.3*0.5 = 1.15 (direct + mediated A_0 -> L_1 -> Y)
#'   psi_0M = 2 (direct A_0*M_0 coefficient; M_0 is baseline,
#'     no mediation through it)
#'
#' Stage-1 blip: gamma_1 = a_1 * (psi_10 + psi_1M * M_1)
#'   psi_10 = 2, psi_1M = 2  (direct coefficients)
#'
#' @param n Number of individuals.
#' @param seed RNG seed.
#' @return List with `data` (person-period data.table with columns
#'   id, time, Y, A, L, M) and `truth_psi` (named numeric vector).
#' @noRd
simulate_snm_longitudinal_tv_em <- function(n = 2000, seed = 42) {
  set.seed(seed)

  L0 <- stats::rnorm(n)
  M0 <- as.numeric(L0 > 0)
  A0 <- stats::rbinom(n, 1, stats::plogis(0.3 * L0))

  L1 <- 0.5 * L0 + 0.3 * A0 + stats::rnorm(n, sd = sqrt(0.5))
  M1 <- as.numeric(L1 > 0)
  A1 <- stats::rbinom(n, 1, stats::plogis(0.3 * L1 + 0.2 * A0))

  Y <- 2 +
    1 * A0 +
    2 * A0 * M0 +
    2 * A1 +
    2 * A1 * M1 +
    1.5 * L0 +
    0.5 * L1 +
    stats::rnorm(n)

  data <- data.table::data.table(
    id = rep(seq_len(n), each = 2L),
    time = rep(0:1, n),
    Y = as.numeric(rbind(NA, Y)),
    A = as.numeric(rbind(A0, A1)),
    L = as.numeric(rbind(L0, L1)),
    M = as.numeric(rbind(M0, M1))
  )

  list(
    data = data,
    truth_psi = c(
      stage0_psi_intercept = 1.15,
      stage0_psi_M = 2,
      stage1_psi_intercept = 2,
      stage1_psi_M = 2
    )
  )
}
