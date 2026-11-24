---
paths:
  - "R/matching.R"
  - "tests/testthat/test-matching*.R"
---

# Matching

These are deliberate. Do not report them as bugs without a numerical reproducer.

- The matching MSM is `Y ~ A`, or `Y ~ A + modifier + A:modifier` with effect modification (`build_matching_msm_formula()`). A bare treatment term in the confounders (`~ L + A`) is rejected by `check_em_compat()`.
- Matching is binary-treatment only. MatchIt rejects other treatments; `fit_matching()` stops earlier with an error that points to `gcomp` or `ipw`.
- The matching variance clusters on `subclass` to aggregate within matched subclasses; the cluster adjustment is needed.
- Matching with transport (`target =`) is rejected: the matching estimand is fixed when the fit is made and cannot take sampling-odds weights.
- A saturated MSM hides some variance errors. Check against a non-saturated MSM before reporting a matching variance problem.
