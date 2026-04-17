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
#   L0 ~ N(0, 1)                              (baseline confounder)
#   A_0 ~ Bern(expit(0.5 * L0))               (treatment at t=0)
#   For t > 0:
#     L_t = A_{t-1} + 0.5 * L0 + ε_L          (treatment-confounder feedback)
#     A_t ~ Bern(expit(0.3 * L_t + 0.2 * A_{t-1}))
#   Y = 10 + 2 * sum(A_t) + L0 + sum(L_t) + ε_Y
#
# True ATE (always vs never) = 3 * n_times - 1
# Under always: E[Y] = 10 + 2*T + (T-1) = 9 + 3T
# Under never:  E[Y] = 10
make_linear_scm <- function(n = 5000, n_times = 2, seed = 42) {
  set.seed(seed)

  L0 <- stats::rnorm(n)

  A <- matrix(NA_real_, n, n_times)
  L <- matrix(NA_real_, n, n_times)

  for (t in seq_len(n_times)) {
    if (t == 1) {
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

  Y <- 10 + 2 * rowSums(A) + L0 + rowSums(L, na.rm = TRUE) + stats::rnorm(n)

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
#   L0 ~ N(0, 1)
#   A_0 = 1 + 0.5 * L0 + ε_A           (continuous treatment)
#   L_1 = A_0 + 0.5 * L0 + ε_L         (treatment-confounder feedback)
#   A_1 = 1 + 0.3 * L_1 + 0.2 * A_0 + ε_A
#   Y = 10 + 2 * (A_0 + A_1) + L0 + L_1 + ε_Y
make_continuous_scm <- function(n = 5000, seed = 42) {
  set.seed(seed)

  L0 <- stats::rnorm(n)
  A0 <- 1 + 0.5 * L0 + stats::rnorm(n, 0, 0.5)
  L1 <- A0 + 0.5 * L0 + stats::rnorm(n, 0, 0.5)
  A1 <- 1 + 0.3 * L1 + 0.2 * A0 + stats::rnorm(n, 0, 0.5)
  Y <- 10 + 2 * (A0 + A1) + L0 + L1 + stats::rnorm(n)

  rbind(
    data.frame(
      id = seq_len(n),
      time = 0L,
      A = A0,
      L = NA_real_,
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
#   A_0 ~ Bern(expit(0.5 * L0))
#   For t > 0:
#     L_t = A_{t-1} + 0.5 * L0 + eps_L  (treatment-confounder feedback)
#     A_t ~ Bern(expit(0.3 * L_t))
#   Y = 10 + (2 + 1.5*sex) * sum(A_t) + L0 + sum(L_t) + eps_Y
#
# True ATE (always vs never):
#   2 periods: ATE|sex=0 = 5,  ATE|sex=1 = 8
#   3 periods: ATE|sex=0 = 8,  ATE|sex=1 = 12.5
#
# Derivation (2-period, always vs never):
#   Under always: E[L_1] = 1, E[Y|sex=s] = 10 + 2*(2+1.5s) + 0 + 1 = 15 + 3s
#   Under never:  E[L_1] = 0, E[Y|sex=s] = 10
#   ATE|sex=s = 5 + 3s => sex=0: 5, sex=1: 8
#
# Derivation (3-period, always vs never):
#   Under always: E[L_1]=1, E[L_2]=1, E[Y|sex=s] = 10+3*(2+1.5s)+0+1+1 = 18+4.5s
#   Under never:  E[L_1]=0, E[L_2]=0, E[Y|sex=s] = 10
#   ATE|sex=s = 8 + 4.5s => sex=0: 8, sex=1: 12.5
make_em_ice_scm <- function(n = 5000, n_times = 2, seed = 42) {
  set.seed(seed)

  L0 <- stats::rnorm(n)
  sex <- stats::rbinom(n, 1, 0.5)

  A <- matrix(NA_real_, n, n_times)
  L <- matrix(NA_real_, n, n_times)

  for (t in seq_len(n_times)) {
    if (t == 1) {
      A[, t] <- stats::rbinom(n, 1, stats::plogis(0.5 * L0))
    } else {
      L[, t] <- A[, t - 1] + 0.5 * L0 + stats::rnorm(n, 0, 0.5)
      A[, t] <- stats::rbinom(n, 1, stats::plogis(0.3 * L[, t]))
    }
  }

  # Treatment effect is (2 + 1.5*sex) per period of treatment.
  trt_effect <- (2 + 1.5 * sex) * rowSums(A)
  Y <- 10 + trt_effect + L0 + rowSums(L, na.rm = TRUE) + stats::rnorm(n)

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
#   Y = 10 + (2 + 1.5*sex) * (A_0 + A_1) + L0 + L_1 + eps_Y
#
# True shift(delta) effect (g-formula, includes indirect path via L_1):
#   MC truth (n = 5*10^6): shift(1)|sex=0 ~ 6.0, shift(1)|sex=1 ~ 9.76
#   (larger than the direct coefficient effect 4/7 because L_1 = A_0 + ...)
make_em_ice_cont_scm <- function(n = 5000, seed = 42) {
  set.seed(seed)

  L0 <- stats::rnorm(n)
  sex <- stats::rbinom(n, 1, 0.5)

  A0 <- 1 + 0.5 * L0 + stats::rnorm(n, 0, 0.5)
  L1 <- A0 + 0.5 * L0 + stats::rnorm(n, 0, 0.5)
  A1 <- 1 + 0.3 * L1 + 0.2 * A0 + stats::rnorm(n, 0, 0.5)

  trt_effect <- (2 + 1.5 * sex) * (A0 + A1)
  Y <- 10 + trt_effect + L0 + L1 + stats::rnorm(n)

  rbind(
    data.frame(
      id = seq_len(n),
      time = 0L,
      A = A0,
      L = NA_real_,
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
#   P[Y(always)|sex=0] ~ 0.767, P[Y(never)|sex=0] ~ 0.287
#   P[Y(always)|sex=1] ~ 0.938, P[Y(never)|sex=1] ~ 0.287
#   RD|sex=0 ~ 0.480, RD|sex=1 ~ 0.651
make_em_ice_binom_scm <- function(n = 5000, seed = 42) {
  set.seed(seed)

  L0 <- stats::rnorm(n)
  sex <- stats::rbinom(n, 1, 0.5)
  A0 <- stats::rbinom(n, 1, stats::plogis(0.5 * L0))
  L1 <- A0 + 0.5 * L0 + stats::rnorm(n, 0, 0.5)
  A1 <- stats::rbinom(n, 1, stats::plogis(0.3 * L1))

  eta <- -1 + (1 + 0.8 * sex) * (A0 + A1) + 0.5 * L0 + 0.3 * L1
  Y <- stats::rbinom(n, 1, stats::plogis(eta))

  rbind(
    data.frame(
      id = seq_len(n),
      time = 0L,
      A = A0,
      L = NA_real_,
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
