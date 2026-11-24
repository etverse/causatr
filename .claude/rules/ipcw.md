---
paths:
  - "R/censoring.R"
  - "R/ipw.R"
  - "R/aipw.R"
  - "R/variance_if.R"
  - "R/variance_if_ipw*.R"
  - "R/variance_if_aipw*.R"
  - "R/variance_if_multinom.R"
  - "tests/testthat/test-*ipcw*.R"
  - "tests/testthat/test-censoring.R"
  - "tests/testthat/test-delicatessen-nhefs.R"
---

# Censoring weights

These are deliberate unless marked as a known bug. Do not report them without a numerical reproducer.

- Built-in IPCW reweights the outcome side only. `causat()` keeps the external, pre-IPCW weights for the treatment models through `propensity_weights` (`causat()` → `fit_ipw()` → `fit_longitudinal_ipw()`, and `refit_ipw()` in the bootstrap); only the folded weights (external × IPCW) reach the MSM and outcome models. This is the separate-then-multiply construction of Hernán and Robins (2020, chapters 12.6 and 17), `W = 1/P(A|H) × 1/P(C=0|H)`. The propensity-unweighted anchor in `test-longitudinal-ipw-ipcw.R` checks it (`prior.weights ≡ 1`, coefficients equal to the unweighted GLM to 1e-8).
- Longitudinal IPW: the censoring term of the sandwich is large, about 5% of the treated-arm SE on an informatively censored process, and is subtracted, like the propensity and transport corrections (Robins, Rotnitzky and Zhao 1994). `compute_ipw_ipcw_correction_longitudinal()` propagates the censoring model through the MSM (γ → β) and adds no γ → α term, because the propensity is not IPCW-weighted. It matches `delicatessen` to 1e-13 for the estimate and 1e-4 for the SE (`longitudinal_ipw_ipcw_delicatessen.py`).
- Longitudinal AIPW: the censoring block (γ) in the stacked estimating equations is required for an exact sandwich even though it moves the SE by only about 0.03%. Double robustness makes `B[μ, γ] = 0` while `B[β, γ] ≠ 0`; `test-aipw-longitudinal-ipcw.R` checks that the block is wired.
- Point g-computation with IPCW carries both the indirect path (γ → outcome β → mean) and the direct path through the IPCW-weighted average, `∂μ/∂γ = (1/Σw) Σ (dwᵢ/dγ)(pᵢ − μ)`, via `ipcw_direct_grad_setup()` and `ipcw_direct_grad()`. Point IPW and AIPW need no direct term, because their marginal mean is a model parameter.
- Known bug: point IPW and point AIPW still fit the propensity IPCW-weighted on the uncensored rows (`fit_ipw()` passes `weights[fit_rows]`). For point IPW this leaves out the γ → α cross-term and the sandwich SE is 3 to 7% off per arm against `delicatessen`; the point estimate is consistent. The fix is to fit the propensity unweighted on all rows with a stacked-estimating-equation sandwich, then tighten the NHEFS IPCW tolerance in `test-delicatessen-nhefs.R` (now 0.1). Adding a γ → α term to the current sandwich is not that fix. Point AIPW with IPCW is non-standard but orthogonal, and its SE is acceptable.
