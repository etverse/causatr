# Phase 10 — Longitudinal IPW

> **Status: PENDING (design doc)**
>
> **Depends on:** Phase 4 (self-contained IPW engine), Phase 5 (ICE data structures/conventions)

## Scope

Extend the self-contained IPW engine to support time-varying treatments. This is the IPW analogue of what ICE g-comp (Phase 5) does for the outcome-model approach: sequential density-ratio weights for time-varying exposures under a marginal structural model.

## Current state

- **ICE g-comp (Phase 5):** Fully handles longitudinal data via backward-iteration outcome models. Sandwich variance via stacked EE.
- **IPW:** `fit_ipw()` hard-aborts on `type == "longitudinal"` with "Longitudinal IPW is not supported."

## Design

### Sequential treatment density models

At each time point k, fit:
```
f(A_k | A̅_{k-1}, L̅_k) via fit_treatment_model()
```

This reuses the existing `fit_treatment_model()` infrastructure on per-period subsets, with treatment history and time-varying confounders as covariates.

### Cumulative product weights

For a longitudinal intervention d = (d_0, d_1, ..., d_K):
```
W_i = ∏_k  f(d_k(A_{i,k}, L_{i,k}) | H_{i,k}) / f(A_{i,k} | H_{i,k})
```

where H_{i,k} = (A̅_{i,k-1}, L̅_{i,k}).

Each per-period ratio is computed via the existing density-ratio machinery. The cumulative product is the longitudinal weight.

### Stabilized weights

Stabilized weights use a simpler numerator model:
```
SW_i = ∏_k  f(d_k | A̅_{k-1}) / f(A_k | A̅_{k-1}, L̅_k)
```

The `numerator` formula parameter (already in `causat()`) specifies what enters the numerator model.

### Time-varying MSM

The weighted MSM is fit on the full person-period data with cumulative weights:
```
E[Y^d | V] = g^{-1}(β₀ + β₁·cum_A + β₂·V + ...)
```

where V is baseline covariates and cum_A summarises the treatment history under intervention d.

For the intercept-only Hájek estimator (matching the point IPW pattern), the per-intervention MSM is simply `Y ~ 1` on the final-period outcomes with cumulative weights.

### Sandwich variance

The stacked M-estimation system includes:
- K propensity model parameter blocks (one per time point)
- 1 MSM parameter block

The bread matrix is block-triangular (each propensity model is fit independently). The cross-derivative chains through all K models, extending `compute_ipw_if_self_contained_one()` to handle the cumulative product.

### Interventions

All longitudinal-compatible interventions from ICE should work:
- `static(a)` — set treatment to a at every time point
- `shift(delta)` — shift by delta at every time point
- `scale_by(factor)` — scale by factor at every time point
- `dynamic(rule)` — apply deterministic rule at every time point
- `threshold(lower, upper)` — **rejected** under IPW (same as point: pushforward has point masses)

### Survival composition (Phase 7 Track B under IPW)

Longitudinal IPW + time-to-event outcome is the IPW analogue of ICE-survival (Phase 7 Track B). It is an in-scope deliverable of Phase 10, not of Phase 7. The composition pattern is fixed now:

- **Per-period weights indexed by $(i, k)$.** The cumulative density-ratio weight
  $$
  W_{i,k} = \prod_{j \leq k} \frac{f(d_j(A_{i,j}, L_{i,j}) \mid H_{i,j})}{f(A_{i,j} \mid H_{i,j})}
  $$
  is already indexed by $(i, k)$ in the non-survival design. Each person-period row receives the cumulative weight up to its own time $k$. This is the standard longitudinal MSM weight (Robins 1999; Cole & Hernán 2008).
- **Weighted hazard MSM on person-period data.** Restrict to uncensored rows and fit
  $$
  \text{logit}\, \hat{h}(k) = \gamma_0 + \gamma(k), \quad \text{weights} = W_{i,k}
  $$
  (intercept + flexible time basis; no $A$ or $L$, following the Hájek pattern). Stabilized weights follow the same pattern as non-survival longitudinal IPW via a `numerator` formula.
- **Survival curve.** Standard cumulative product: $\hat{S}^d(t) = \prod_{k \leq t}(1 - \hat{h}(k))$. Because the hazard MSM has only an intercept + time basis, $\hat{h}(k)$ is constant across individuals at a given $k$, and the cumulative product is a population-level survival curve (covariate adjustment is already baked into $W_{i,k}$).
- **Variance.** The stacked EE system grows to $K$ propensity-model blocks (one per time period) + the hazard-MSM block. Channel 2's chain-rule cross-derivative through the $K$ propensity models is already required for non-survival longitudinal IPW (§ "Sandwich variance"). The survival addition is the cross-time delta-method aggregation from Phase 7 § "Cross-time variance aggregation": per-individual IF matrix of shape $n \times |t\text{-grid}|$; pointwise SEs from the diagonal; RMST SEs from quadratic forms.
- **Oracle.** `lmtp::lmtp_tmle(outcome_type = "survival", cens = <col>, trt = <list>, shift = ...)` is directly applicable and well-validated. Cross-check point estimates against it; cross-method triangulation against ICE-hazards (Track B under gcomp) on a DGP where both estimators are consistent.
- **Censoring.** Phase 10 + survival uses the `censoring =` row-filter pathway as today (restricts the hazard MSM to uncensored person-period rows). Built-in IPCW + longitudinal survival is owned by Phase 14 (survival is Phase 14's motivating use case).

**Three-way composition with Phase 8 (multivariate longitudinal survival IPW) is out of scope for Phase 10.** It requires joint per-period density models ($f(A_{1,k}, A_{2,k} \mid H_k)$) and lands after both Phase 8 and Phase 10 ship. Mentioned here for completeness only.

## Items

- [ ] Remove `type == "longitudinal"` gate in `fit_ipw()`
- [ ] `fit_longitudinal_ipw()`: fit per-period treatment density models
- [ ] `compute_longitudinal_weights()`: cumulative product density-ratio weights
- [ ] `compute_ipw_contrast_longitudinal()`: weighted MSM on final-period outcomes
- [ ] Stabilized weights via `numerator` formula
- [ ] `variance_if_ipw_longitudinal()`: stacked sandwich for K propensity models + MSM
- [ ] Bootstrap: resample individuals (all person-period rows together)
- [ ] Sequential positivity warnings (redistributed from old Phase 7)
- [ ] Truth-based tests: linear-Gaussian DGP validated against ICE g-comp and/or `lmtp::lmtp_ipw()`
- [ ] Cross-method agreement: longitudinal IPW vs ICE g-comp on same DGP
- [ ] **Survival composition:** cumulative $(i,k)$-indexed weights, weighted pooled-logistic hazard MSM on uncensored person-period rows, cumulative-product survival curve, stacked sandwich across $K$ propensity models + hazard MSM with cross-time delta; cross-check against `lmtp::lmtp_tmle(outcome_type = "survival")` and triangulate with ICE-hazards (Phase 7 Track B under gcomp). *Depends on Phase 7 chunks 7a–7b and 7e.*

## Dependencies

Phase 4 (self-contained IPW engine), Phase 5 (ICE data structures, person-period conventions).

## Out of scope

- Longitudinal multivariate IPW (combine with Phase 8 multivariate)
- Grace period / visit-process interventions (future enhancement)
- Stratified ICE option (future enhancement)
