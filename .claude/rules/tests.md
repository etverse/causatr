---
paths:
  - "tests/testthat/**"
---

# Tests

- Run a file with `Rscript -e 'devtools::load_all(quiet = TRUE); testthat::test_file("tests/testthat/test-<name>.R")'`. Slow blocks call `skip_if_fast()` (`helper-parallel.R`); `CAUSATR_TEST_TIER=fast` skips them, and the default tier is `full`.
- Shared data-generating processes are in `helper-dgp.R` and `helper-glmtp-dgp.R`; oracle wrappers are in the `helper-*-oracle.R` files; Python oracle scripts and their outputs are in `fixtures/python/`.
- In a hand-built longitudinal DGP with interleaved rows (`id = rep(1:n, each = 2)`, `time = rep(1:2, n)`), the end-of-follow-up outcome goes on the final-period rows: `Y = c(rbind(rep(NA_real_, n), Y_obs))`, matching `A = c(rbind(A_1, A_2))`. `Y = c(rep(NA_real_, n), Y_obs)` scrambles the outcome against the covariates, and a test with only loose assertions still passes on the meaningless estimate. When reviewing a longitudinal test, check the outcome construction first.
- External oracles are compared at the contrast level, not the influence-function level.
- `lmtp::lmtp_sdr()` is an oracle for multivariate longitudinal shift interventions too (`trt = list(c("A1_t", "A2_t"), ...)` with a shift over those columns). causatr's multivariate longitudinal IPW shift agrees with it to about 3%.
