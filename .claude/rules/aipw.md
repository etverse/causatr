---
paths:
  - "R/aipw*.R"
  - "R/variance_if_aipw*.R"
  - "tests/testthat/test-aipw*.R"
  - "tests/testthat/helper-aipw-long-oracle.R"
---

# AIPW

These are deliberate. Do not report them as bugs without a numerical reproducer.

- The longitudinal AIPW sandwich is the full stacked-estimating-equation sandwich, `V = B⁻¹ M B⁻ᵀ / n`, in `variance_if_aipw_longitudinal()` and `R/variance_if_aipw_longitudinal_ee.R`. `ψ(θ)` stacks the per-period propensity scores, the per-step outcome and pseudo-outcome scores and the marginal-mean equation, and the bread is the numerical Jacobian (`numDeriv`) of the summed score. That captures every cross-term, including the baseline pseudo-regression block `B⁻¹[μ, β₁]`, so the SE is consistent on balanced and unbalanced panels (monotone dropout, censoring row filter). A forward-cascade assembly that drops that block underestimates the SE by about 50% under dropout; do not reintroduce one, or an abort for unbalanced panels.
- `aipw_long_stacked_if()` checks that the reconstructed score vanishes (`max|colSums(ψ(θ̂))|`) and otherwise raises `causatr_longitudinal_aipw_sandwich_unfaithful` and points to the bootstrap.
- The analytic sandwich covers GLM nuisance models, `MASS::glm.nb` and `nnet::multinom` propensities. `mgcv::gam` and `betareg::betareg` nuisance models raise `causatr_longitudinal_aipw_sandwich_model` and use the bootstrap, because their bread is not the score Jacobian.
- Validation: the summed score vanishes to 1e-11; the SE matches the `delicatessen` M-estimation oracle to 1e-7 (`aipw_long_tau2_delicatessen.py`, balanced and unbalanced) and the bootstrap and delete-one-id jackknife to a ratio of about 1.0.
- Longitudinal AIPW rejects a covariate missing within an observed person-period with `causatr_longitudinal_aipw_missing_covariate`: the recursion has no complete-case fallback there, unlike pooled ICE. Missing time-varying covariates go through `causat_mice()` and missing outcomes through IPCW; do not add complete-case handling to silence the error.
- `stabilize = "marginal"` is rejected for longitudinal AIPW (`causatr_stabilize_longitudinal_aipw`): the recursion uses Hájek-normalised single-period weights. Transport (`target =`) is not supported for longitudinal AIPW.
- Longitudinal IPSI is rejected under AIPW: there is no counterfactual treatment stream for the recursion.
