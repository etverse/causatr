# Phase 8 — Multivariate Treatment IPW

> **Status: DONE (2026-04-20)**

## Scope

Extend the self-contained IPW density-ratio engine to support multivariate (joint) treatments `treatment = c("A1", "A2")`. Point g-comp already handles multivariate treatments; this phase closes the IPW gap.

## Current state

- **G-comp (point):** `fit_gcomp()` accepts `treatment = c("A1", "A2")` — the outcome model `Y ~ A1 + A2 + L` handles it natively. `apply_intervention()` loops over each treatment variable. Sandwich and bootstrap work.
- **IPW:** `fit_ipw()` hard-aborts at `length(treatment) > 1L` (R/ipw.R) with "Multivariate treatments are not supported under `estimator = 'ipw'`."
- **Matching:** `fit_matching()` also blocks multivariate — MatchIt itself only handles binary treatments.

## Design

### Joint treatment density

The key question: how to model `f(A1, A2 | L)`.

**Option A: Product of conditionals (sequential factorisation)**
```
f(A1, A2 | L) = f(A1 | L) · f(A2 | A1, L)
```
Each factor is a standard univariate density model that `fit_treatment_model()` already handles. `A1` enters the second model as an additional covariate. This is the simplest approach and mirrors how longitudinal IPW chains per-period models.

**Option B: Multivariate normal (continuous case)**
```
f(A1, A2 | L) = MVN(μ(L), Σ(L))
```
Requires a multivariate regression model. More statistically efficient but harder to implement and less flexible (can't mix binary + continuous treatments).

**Recommendation:** Option A. It reuses the existing `fit_treatment_model()` machinery, supports mixed treatment types (e.g. one binary, one continuous), and the product factorisation is valid under the full model.

### Joint density-ratio weights

For a joint intervention `(a1*, a2*)`:
```
w_i = f(a1*_i | L_i) · f(a2*_i | a1*_i, L_i) / [f(A1_i | L_i) · f(A2_i | A1_i, L_i)]
```

Each factor ratio is computed via the existing `compute_density_ratio_weights()` on the corresponding univariate model. The product is the joint weight.

### Sandwich variance

The stacked M-estimation system includes parameters from both density models + the MSM. The variance engine extends `compute_ipw_if_self_contained_one()` to chain through multiple propensity models — each contributes a block to the bread matrix and a column block to the cross-derivative.

### Interventions

All intervention types that work on univariate IPW (`static`, `shift`, `scale_by`, `dynamic`, `ipsi`) should work per-component. The user specifies a list of interventions, one per treatment variable:
```r
contrast(fit,
  interventions = list(
    both_up = list(shift(1), shift(0.5)),
    control = list(static(0), static(0))
  ),
  reference = "control"
)
```

## Items

- [x] Remove `length(treatment) > 1L` gate in `fit_ipw()` and `check_causat_inputs()`.
- [x] Add `fit_treatment_models()` (plural) in `R/treatment_model.R` to fit the sequential factorisation `f(A_k | A_{1..k-1}, L)` per component.
- [x] Add `compute_density_ratio_weights_mv()` in `R/ipw_weights.R`. The k-th factor's denominator stays at observed conditioning; the numerator substitutes the inverse-map values of upstream components into the conditioning columns. Reuses `evaluate_density()` per component but does NOT delegate to the univariate `compute_density_ratio_weights()` (which would use one newdata for both halves and wipe out the cross-component conditioning shift).
- [x] Add `make_weight_fn_mv()` building a stacked-alpha closure across K models for the variance engine; per-component sub-closures via `mv_natural_course_closure()` / `mv_ht_closure()` / `mv_pushforward_closure()`. `force()` every captured arg to avoid the R for-loop promise gotcha.
- [x] Add `intervention_inverse_map()` helper that returns d_k^{-1}(A_k) per intervention type (the value plugged into downstream conditioning).
- [x] Add `compute_ipw_if_self_contained_mv_one()` in `R/variance_if.R`. Computes the stacked cross-derivative `[A_{β,α₁}, ..., A_{β,α_K}]` via `numDeriv::jacobian` on the product-weight closure, then sums K block-diagonal propensity corrections (one `apply_model_correction()` per propensity model).
- [x] Wire dispatch in `fit_ipw()`, `compute_ipw_contrast_point()`, `variance_if_ipw()` on `length(treatment) > 1L` (stored as `fit$details$is_multivariate`).
- [x] Bootstrap path (`refit_ipw()`) replays the same `treatment` slot — multivariate just works through the existing `replay_fit()` plumbing.
- [x] Reject IPSI / categorical / count components and effect modification under multivariate IPW with classed errors (`causatr_multivariate_*`).
- [x] Truth-based tests in `tests/testthat/test-multivariate-ipw.R`: 12 tests covering binary × binary, binary × continuous, continuous × continuous, K = 3 binary, binom outcome with diff/ratio/OR, by-stratified, subset, dynamic, gcomp cross-check, sandwich-vs-bootstrap parity. Plus 5 rejection tests.
- [x] Update `FEATURE_COVERAGE_MATRIX.md`, `CLAUDE.md`, `NEWS.md`.

## Dependencies

Phase 4 (self-contained IPW engine). Independent of Phases 6, 9, 10.

## Out of scope

- Multivariate matching (MatchIt limitation)
- Longitudinal multivariate IPW (combine with Phase 10)
- Effect modification for multivariate treatment (Phase 6 + Phase 8 interaction; deferred)
