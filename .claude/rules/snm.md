---
paths:
  - "R/snm*.R"
  - "R/variance_if_snm*.R"
  - "R/variance_bootstrap_snm.R"
  - "tests/testthat/test-snm*.R"
---

# Structural nested models

These are deliberate. Do not report them as bugs without a numerical reproducer.

- An SNM estimates blip parameters ψ, not counterfactual means, so `contrast()` with `interventions =` is rejected for SNM fits.
- `treatment_free = ~ L` estimates (β, ψ) jointly (Vansteelandt and Joffe 2014). It absorbs the L → Y variance and leaves the point estimate unchanged under a correctly specified blip.
- With `treatment_free` in the longitudinal SNM, `variance_if_snm_longitudinal_tf()` builds the joint per-stage `(βₖ, ψₖ)` system with its cross-stage derivatives and then takes the ψ margin. The residual is `Hₖ − Zₖβₖ − γₖ`; the ψ-only residual `Hₖ − γₖ` is wrong. The error cancels when `E[R·Z] ≈ 0`, so a test with a correctly specified treatment model does not catch it: test with a treatment-free formula that includes covariates absent from the propensity model.
- Bootstrap intervals for the by-stratified averaged blip are Wald intervals, because there are no per-stratum replicates; its SE is the delta method on the global vcov.
