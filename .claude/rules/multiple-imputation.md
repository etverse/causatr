---
paths:
  - "R/causat_mice.R"
  - "R/pool_*.R"
  - "tests/testthat/test-causat-mice.R"
  - "tests/testthat/test-pool-*.R"
  - "tests/testthat/test-mi-coverage.R"
---

# Multiple imputation

- `causat_mice()` is the analysis step for covariates or treatment missing at random. It loops `causat()` and `contrast()` over the completed datasets of a `mice` `mids` object and pools, so it works for any `causat()` configuration. It does not impute, and it never imputes the outcome: the outcome is a predictor in the imputation model, and a missing outcome is IPCW's job.
- Pooling is per row: means and differences on the natural scale, ratios and odds ratios on the log scale, SNM per blip parameter. Diagnostics per row are in the `"mi_details"` attribute, which `confint()` and `tidy()` read to build t-based intervals.
- `pool_method = "rubin"` (Barnard-Rubin degrees of freedom, as `mice::pool.scalar()`) is exact under congeniality; under uncongeniality its variance can be biased in either direction (Meng 1994; Bartlett and Hughes 2020). `pool_method = "boot_mi"` (von Hippel 2020, bootstrap then impute) reaches nominal coverage only while the point estimator stays consistent: a misspecified imputation that biases the estimate defeats both rules.
- Transport with multiple imputation has not been validated.
