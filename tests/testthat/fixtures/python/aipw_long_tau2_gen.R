# Generate the shared wide-format data for the longitudinal-AIPW
# delicatessen cross-check (tau = 2, binary treatment, continuous outcome,
# static(1) vs static(0)). Writes one balanced and one unbalanced
# (monotone-dropout) data set. Both nuisance models are correctly
# specified, so the analytic sandwich SE is consistent and must match
# delicatessen's M-estimation sandwich.
#
# Person-period encoding for causatr (confounders = ~1, confounders_tv = ~L):
#   time 0 row: L = L0,  A = A0,  Y = NA
#   time 1 row: L = L1,  A = A1,  Y = Y
# Per-period propensity:  A ~ 1 + L (+ lag1_A + lag1_L at t = 1)
# Final outcome model:    Y ~ A + L + lag1_A + lag1_L
#
# Run from the repository root:
#   Rscript tests/testthat/fixtures/python/aipw_long_tau2_gen.R

set.seed(2025)
n <- 4000

L0 <- rnorm(n)
A0 <- rbinom(n, 1, plogis(0.2 + 0.6 * L0))
L1 <- rnorm(n, mean = 0.4 + 0.5 * A0 + 0.3 * L0)
A1 <- rbinom(n, 1, plogis(0.1 + 0.5 * L1 + 0.4 * A0))
Y <- rnorm(n, mean = 2 + 1.5 * A0 + 2 * A1 + 1.0 * L1 + 0.5 * L0)

# Monotone dropout: treated-at-baseline units drop out before t = 1 with
# probability 0.35 (informative w.r.t. A0). R1 = 1 iff present at t = 1.
R1 <- rep(1L, n)
drop <- A0 == 1L & runif(n) < 0.35
R1[drop] <- 0L

wide <- data.frame(
  id = seq_len(n),
  L0 = L0,
  A0 = A0,
  L1 = L1,
  A1 = A1,
  Y = Y,
  R1 = R1
)

here <- "tests/testthat/fixtures/python"

# Balanced data set: every unit observed at both periods.
utils::write.csv(
  wide[, c("id", "L0", "A0", "L1", "A1", "Y")],
  file.path(here, "aipw_long_tau2_balanced_data.csv"),
  row.names = FALSE
)

# Unbalanced data set: drop the t = 1 covariate / treatment / outcome for
# units with R1 == 0 (NA in the wide file; Python gates them out, R drops
# the corresponding person-period rows).
wide_unb <- wide
wide_unb$L1[wide_unb$R1 == 0L] <- NA
wide_unb$A1[wide_unb$R1 == 0L] <- NA
wide_unb$Y[wide_unb$R1 == 0L] <- NA
utils::write.csv(
  wide_unb,
  file.path(here, "aipw_long_tau2_unbalanced_data.csv"),
  row.names = FALSE
)

cat("wrote balanced (n =", n, ") and unbalanced (dropout =",
  sum(wide$R1 == 0L), ") data sets\n")
