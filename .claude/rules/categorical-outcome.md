---
paths:
  - "R/variance_if_multinom.R"
  - "R/contrast.R"
  - "R/family.R"
  - "tests/testthat/test-gcomp-categorical-outcome.R"
  - "tests/testthat/helper-multinom-oracle.R"
---

# Categorical outcomes

- A K-level factor outcome fitted with `model_fn = nnet::multinom` has a vector estimand, `P(Y = k | do(A = a))` per intervention, summing to 1 over k. It is one outcome, not several.
- Detection is after the fit (`is_multinom_outcome()`), and `compute_contrast()` branches early to `compute_contrast_multinom()` and returns before the scalar code, so a scalar result is unchanged by this path. The categorical result adds a `class` column, stores `vcov` as a list per class and records `class_labels`; both paths share `compute_pairwise_contrasts()`, and the odds ratio is one class against the rest.
- Variance is the analytic sandwich per class (`variance_if_gcomp_multinom()`, built on the stacked multinomial score of `prepare_model_if_multinom()`) for complete-case, weighted and IPCW fits, or the bootstrap (`variance_bootstrap_multinom()`). With weights, `w ≡ 1` reproduces the unweighted path exactly.
- Supported for point g-computation only: `check_categorical_outcome()` rejects SNM (`causatr_snm_categorical_outcome`) and IPW, AIPW, matching, longitudinal and transport fits (`causatr_categorical_outcome_unsupported`), and `contrast()` rejects stochastic interventions.
- Oracles: a large-n softmax truth, `marginaleffects::avg_predictions()` on the same fit (agreement near 1e-15), and the Python M-estimation stacks in `fixtures/python/multinom_gcomp_*.py`. `marginaleffects` is a point oracle only: its delta-method SE conditions on the standardisation sample and leaves out the empirical-distribution term that the influence-function sandwich carries, so the SE oracle is the M-estimation stack.
