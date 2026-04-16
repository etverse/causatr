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
