---
paths:
  - "R/**"
  - "tests/testthat/**"
---

# Combination axes

A combination audit runs over these axes. The values are the ones the code accepts; whether a combination is supported is decided by the rejection paths in `R/checks.R` and the fitters, and recorded per combination in `FEATURE_COVERAGE_MATRIX.md`. Re-read the argument lists when an axis may have changed.

| Axis | Values | Where |
|---|---|---|
| estimator | `gcomp`, `ipw`, `aipw`, `matching`, `snm` | `causat(estimator =)` |
| treatment timing | point; longitudinal when `id` and `time` are given (ICE g-computation, longitudinal IPW, AIPW and SNM); transport when `target` is given | `causat()` |
| treatment type | binary, continuous, categorical, count, multivariate (a vector `treatment`) | `causat(treatment =, propensity_family =)` |
| outcome model | any `stats::glm` family, `MASS::glm.nb`, `betareg::betareg`, `mgcv::gam`, `nnet::multinom` (categorical outcome) | `causat(model_fn =, family =)` |
| intervention | `static()`, `shift()`, `scale_by()`, `threshold()`, `dynamic()`, `ipsi()`, `stochastic()`, and the natural-history policies `grace_period()`, `carry_forward()`, `cap_escalation()` | `R/interventions.R`, `R/glmtp_interventions.R` |
| estimand | `ATE`, `ATT`, `ATC`, and `by`-stratified results | `causat(estimand =)`, `contrast(by =)` |
| contrast scale | `difference`, `ratio`, `or` | `contrast(type =)` |
| variance | `sandwich` (with a numerical fallback), `bootstrap` with `boot_ci` `percentile` or `normal` | `contrast(ci_method =, boot_ci =)` |
| weights | none, external or survey `weights`, `stabilize` `none` or `marginal` | `causat()` |
| missing data | complete case, `censoring` as a row filter, IPCW for missing outcome (`ipcw = TRUE`), multiple imputation through `causat_mice()` | `causat()`, `R/causat_mice.R` |
| nuisance options | `treatment_free` (SNM), `stratified` and `treatment_form` (longitudinal g-computation) | `causat()` |
