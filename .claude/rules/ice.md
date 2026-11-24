---
paths:
  - "R/ice*.R"
  - "R/glmtp*.R"
  - "R/variance_if_ice*.R"
  - "R/variance_if_glmtp.R"
  - "R/to_person_period.R"
  - "tests/testthat/test-ice*.R"
  - "tests/testthat/test-glmtp.R"
  - "tests/testthat/helper-glmtp-dgp.R"
---

# ICE g-computation and natural-history policies

These are deliberate. Do not report them as bugs without a numerical reproducer.

## ICE

- The intervention applies to the current-time treatment only. Lag columns hold the observed `A_{k-1}, A_{k-2}, …` at every backward step; recomputing them from the intervened treatment double-counts shift, scale, threshold and dynamic interventions.
- ICE fits its sequential outcome models in `contrast()`, not in `causat()`, because they depend on the intervention.
- Under non-uniform weights, the step > 0 cascade gradient in the ICE sandwich carries the `w_{k-1,j}` factor from the previous weighted fit; dropping it halves the SE. The fit-row bread times the target-row gradient is correct by the delta method. Run sandwich against bootstrap before calling it a scope mismatch.
- `variance_if_ice_chain()` restricts each step's counterfactual design to model-complete rows (`ice_model_complete_rows()`), because `model.matrix()` drops rows with a missing term and would desynchronise the design from the id-indexed `match()` vectors. The dropped ids are exactly those outside step k's gated score (`I(C_{i,k} = 0)`, Zivich et al. 2024); do not add them back. The baseline pseudo-outcome is defined for everyone, so the point estimate and the sampling term are unchanged, and the correction is valid under MCAR only. The missing-covariate sandwich matches bootstrap and jackknife to a ratio of about 1.0 (`test-missing-data.R`).
- Before reporting an ICE variance bug, derive the influence-function decomposition on paper and compare sandwich with bootstrap numerically.
- Stratified ICE (`stratified = "G"`) fits one model per level of a baseline column at each step. Since `G` is time-invariant and people never change stratum, it equals pooled ICE run within each stratum; the estimand is unchanged. It is longitudinal g-computation only (`causatr_stratified_not_ice`).
- `treatment_form = ~ factor(A)` or `~ splines::ns(A, df)` changes only the treatment's design term in the per-step models; the intervention still writes the numeric treatment column. `NULL` reproduces the bare numeric term exactly. Longitudinal g-computation only (`causatr_treatment_form_not_ice`).

## Natural-history policies

- `grace_period()`, `carry_forward()` and `cap_escalation()` are an intervention type (class `causatr_glmtp`), not an estimator; the estimator stays `gcomp`. `contrast()` routes them to `glmtp_iterate()` instead of `ice_iterate()`. There is deliberately no `estimator = "glmtp"`.
- `glmtp_iterate()` conditions on the observed lag columns and uses the natural-history label `s̄ₜ`, not the lags, as the policy input. That decoupling is the method (Díaz, Williams, Morzywołek and Rudolph 2026, arXiv:2605.24167): the tower property integrates over the observed distribution and the intervention enters only through the current-period policy value.
- `lmtp` is not an oracle for history-dependent policies: its one-shot shift computes the standard LMTP, a different estimand. The oracle is a forward Monte Carlo of the natural-history regime (Proposition 1 of the paper; `helper-glmtp-dgp.R`).
- Both variance paths are valid: the id-cluster bootstrap and the analytic sandwich in `variance_if_glmtp()`.
- A multivariate natural-history policy is rejected by design (`causatr_glmtp_multivariate`), not missing. The paper develops the augmented recursion for a scalar discrete exposure only, and a joint support enumerates `|A⁽¹⁾ × … × A⁽ᴷ⁾|^(τ−1)` labels. Do not add it by joint coding.
- `glmtp_iterate()` reads its final `Q[[1L]]` positionally behind a `length(Q) == 1L` guard because the element's name is `""`, and `list[[""]]` returns `NULL` in R. Do not replace it with a name lookup.

## Tests

- `test-glmtp.R` checks `cap_escalation()` with `treatment_form = ~ factor(A)` at n = 40000 by comparing the bare and factor estimators (`gap_factor < 0.015`, `gap_bare > 0.02`, `gap_factor < gap_bare / 2`). The factor estimator's gap to the truth is 0.0015 to 0.010 at that size depending on the seed and settles near 0.001 to 0.0025 only from n = 80000, while the bare estimator's bias is a stable 0.03 to 0.04. An absolute 0.005 band at n = 40000 fails on some seeds; verify at n ≥ 80000 before reading a factor gap as bias.
- The `cap_escalation()` bootstrap test at n = 3000 checks plumbing: its load-bearing assertion is the 1e-8 agreement between `contrast()` and a direct `glmtp_iterate()` call. At that size a CI half-width near 0.05 cannot separate the consistent estimator from the biased one, so coverage or `cap < natural` assertions there test nothing.
