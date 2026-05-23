#!/usr/bin/env Rscript
# Generate fixture CSV for the TV-EM longitudinal SNM delicatessen cross-check.
# Replicates the simulate_snm_longitudinal_tv_em() DGP from helper-dgp.R
# with n=5000, seed=101 and writes wide-format CSV for Python consumption.

set.seed(101)
n <- 5000L

L0 <- rnorm(n)
M0 <- as.numeric(L0 > 0)
A0 <- rbinom(n, 1, plogis(0.3 * L0))

L1 <- 0.5 * L0 + 0.3 * A0 + rnorm(n, sd = sqrt(0.5))
M1 <- as.numeric(L1 > 0)
A1 <- rbinom(n, 1, plogis(0.3 * L1 + 0.2 * A0))

Y <- 2 +
  1 * A0 +
  2 * A0 * M0 +
  2 * A1 +
  2 * A1 * M1 +
  1.5 * L0 +
  0.5 * L1 +
  rnorm(n)

wide <- data.frame(
  id = seq_len(n),
  A0 = A0, L0 = L0, M0 = M0,
  A1 = A1, L1 = L1, M1 = M1,
  Y = Y
)

write.csv(wide, "data-raw/snm_longitudinal_tv_em_fixture.csv", row.names = FALSE)
write.csv(
  wide,
  "tests/testthat/fixtures/snm_longitudinal_tv_em_fixture.csv",
  row.names = FALSE
)

cat("Wrote", nrow(wide), "rows to fixture CSVs\n")
