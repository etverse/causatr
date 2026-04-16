# Phase 8 — Multivariate Treatment IPW

> **Status: PENDING (design doc)**

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

- [ ] Remove `length(treatment) > 1L` gate in `fit_ipw()`
- [ ] Extend `fit_treatment_model()` to fit sequential factorisation: `f(A1|L)`, `f(A2|A1,L)`, ...
- [ ] Extend `compute_density_ratio_weights()` for product weights
- [ ] Extend `make_weight_fn()` for product weight closure (for variance engine)
- [ ] Extend `compute_ipw_if_self_contained_one()` for multi-model propensity sandwich
- [ ] Extend `apply_intervention()` for per-component intervention lists
- [ ] Truth-based tests: 2-treatment DGP validated against multivariate g-comp
- [ ] Mixed-type test: one binary + one continuous treatment

## Dependencies

Phase 4 (self-contained IPW engine). Independent of Phases 6, 7, 9, 10.

## Out of scope

- Multivariate matching (MatchIt limitation)
- Longitudinal multivariate IPW (combine with Phase 10)
- Effect modification for multivariate treatment (Phase 6 + Phase 8 interaction; deferred)
