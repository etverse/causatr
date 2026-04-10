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
