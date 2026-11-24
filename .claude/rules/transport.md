---
paths:
  - "R/sampling_model.R"
  - "R/variance_if_ipw_sampling.R"
  - "R/diagnose_sampling.R"
  - "tests/testthat/test-*transport*.R"
  - "tests/testthat/test-sampling-model.R"
---

# Transport

- Transport fits a sampling model `P(S = 1 | L)` on the combined study and target data and weights study rows by the sampling odds `(1 − p̂)/p̂`. `target_subset = "target"` transports to the target population; `"all"` generalises to study plus target.
- A modified treatment policy with transport marginalises over `P(A | L, S = 1)` by Monte Carlo, because target rows have no observed treatment: exact enumeration for a binary treatment, Monte Carlo integration for a continuous one. Only the bootstrap is supported for this combination.
- Matching with transport is rejected (see the matching rule).
