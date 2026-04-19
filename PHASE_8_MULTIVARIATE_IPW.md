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

### Survival composition (Phase 7 Track A under IPW)

Multivariate + point-survival composes cleanly and is an in-scope deliverable of Phase 8, not a Phase 7 deliverable. The composition pattern is fixed now:

- **Joint density is a baseline property.** For point survival, $A_1, A_2$ are measured at baseline and $L$ is baseline. The joint density $f(A_1, A_2 \mid L)$ is fit **once on the original-row data** (one row per individual), using the same sequential factorisation as in the non-survival case. Fitting on the expanded person-period data would double-count individuals and distort the likelihood.
- **Weight broadcast.** The product density-ratio weight
  $$
  w_i = \frac{f(a_1^* \mid L_i) \cdot f(a_2^* \mid a_1^*, L_i)}{f(A_{1,i} \mid L_i) \cdot f(A_{2,i} \mid A_{1,i}, L_i)}
  $$
  is computed once per individual and broadcast onto every person-period row for that individual. $w_i$ is constant in $k$ because neither the numerator nor the denominator depends on time for point treatments.
- **Hazard MSM.** The weighted pooled-logistic hazard is fit on person-period data with the broadcast weights:
  $$
  \text{logit}\, \hat{h}(k) = \gamma_0 + \gamma(k)
  $$
  (intercept + flexible time basis, no $A$ or $L$ — IPW has already absorbed covariate adjustment into the weights, following the same Hájek $Y \sim 1$ pattern used by univariate IPW). Survival curve: $\hat{S}^{(a_1^*, a_2^*)}(t) = \prod_{k \leq t} (1 - \hat{h}(k))$.
- **Variance.** The stacked EE system has **two** propensity-model blocks (for $A_1$ and $A_2 \mid A_1$) plus the hazard-MSM block. Channel 2 in `compute_ipw_if_self_contained_one()` extends to chain through both propensity models (already required for non-survival multivariate; see § "Sandwich variance"). The survival-specific addition is the cross-time delta-method aggregation from Phase 7 § "Cross-time variance aggregation": the per-individual IF is an $n \times |t\text{-grid}|$ matrix, with row $i$'s weight-derivative broadcast across every period at-risk for $i$. The IF stacking across time for a given individual is just the delta-method chain on the cumulative product.
- **Mixed-type survival DGP.** One binary + one continuous baseline treatment with a Weibull or piecewise-exponential event time is the canonical test. External oracle: none directly — `lmtp::lmtp_tmle` does not support joint multivariate baseline treatment. Cross-check instead against multivariate g-comp + pooled logistic (Track A under gcomp), which already supports multivariate treatment for non-survival.

**Longitudinal multivariate + survival is out of scope for Phase 8.** Combining Phase 8 (multivariate) with Phase 10 (longitudinal) with Phase 7 (survival) is a three-way composition documented in Phase 10 § "Survival composition"; it is explicitly deferred past Phase 8.

## Items

- [ ] Remove `length(treatment) > 1L` gate in `fit_ipw()`
- [ ] Extend `fit_treatment_model()` to fit sequential factorisation: `f(A1|L)`, `f(A2|A1,L)`, ...
- [ ] Extend `compute_density_ratio_weights()` for product weights
- [ ] Extend `make_weight_fn()` for product weight closure (for variance engine)
- [ ] Extend `compute_ipw_if_self_contained_one()` for multi-model propensity sandwich
- [ ] Extend `apply_intervention()` for per-component intervention lists
- [ ] Truth-based tests: 2-treatment DGP validated against multivariate g-comp
- [ ] Mixed-type test: one binary + one continuous treatment
- [ ] **Survival composition:** broadcast joint weight onto person-period rows; weighted pooled-logistic hazard MSM; cross-time delta-method IF aggregation on $\hat{S}^{(a_1^*, a_2^*)}(t)$; test against multivariate gcomp + pooled logistic on a Weibull/piecewise-exponential DGP. *Depends on Phase 7 Track A (chunks 7a–7d).*

## Dependencies

Phase 4 (self-contained IPW engine). Independent of Phases 6, 7, 9, 10.

## Out of scope

- Multivariate matching (MatchIt limitation)
- Longitudinal multivariate IPW (combine with Phase 10)
- Effect modification for multivariate treatment (Phase 6 + Phase 8 interaction; deferred)
