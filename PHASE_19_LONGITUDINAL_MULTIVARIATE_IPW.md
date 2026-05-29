# Phase 19 — Longitudinal Multivariate IPW

> **Status: COMPLETE — 19-trim, 19a, 19b, 19c all shipped.**
>
> **Depends on:** Phase 8 (multivariate IPW), Phase 10 (longitudinal IPW).

## Scope

Combine Phase 8 (multivariate IPW via sequential factorisation
$\prod_k f_k(A_k \mid A_{1..k-1}, L)$) with Phase 10 (longitudinal IPW
via per-period density chain $\prod_t f_t(A_t \mid \bar A_{t-1}, \bar L_t)$)
to support **joint time-varying treatments**: `treatment = c("A1", "A2")`
combined with `id =` and `time =`. Each individual contributes one
$(A_1, A_2)$ vector at every period.

The full factorisation chains over BOTH axes (time × component):

$$
f(\bar A_1, \bar A_2 \mid \bar L) =
\prod_{t=1}^{T} \prod_{k=1}^{K}
f_{t,k}\bigl(A_{t,k} \;\Big|\; A_{t,1..k-1},\, \bar A_{t-1},\, \bar L_t\bigr).
$$

## Current state

The chunk 10a `causatr_longitudinal_multivariate_pending` rejection
in `R/longitudinal_ipw.R` hard-aborts when `length(treatment) > 1L`
under `type = "longitudinal"`. Both Phase 8 (multivariate point IPW)
and Phase 10 (univariate longitudinal IPW) ship the building blocks
this phase composes.

## Design

### Sequential factorisation across time × component

At each `(time_points[t], treatment[k])`:

1. Subset to rows where `time == time_points[t]` (one row per id).
2. Build the conditioning formula
   `A_k ~ A_{1..k-1} + lag1_A_1 + ... + lag1_A_K + L + lag1_L + baseline`.
   Treatment lags from the prior period AND prior components within
   the current period both enter the conditioning set.
3. Fit `f_{t,k}` via `fit_treatment_model()` on that subset.

This produces $T \times K$ density models stored as a 2D list
indexed by `(time, component)`. Number of independent propensity
parameters: $\sum_{t,k} p_{t,k}$.

### Cumulative product weight

$$
W_i = \prod_{t=1}^{T} \prod_{k=1}^{K}
\frac{f_{t,k}\bigl(d_k(A_{t,k}) \mid \bar A_{t,1..k-1}^{\mathrm{obs}}, \bar A_{t-1}^{\mathrm{obs}}, \bar L_t\bigr) \cdot \lvert \mathrm{Jac}\,d_k^{-1} \rvert}{f_{t,k}\bigl(A_{t,k} \mid \bar A_{t,1..k-1}^{\mathrm{obs}}, \bar A_{t-1}^{\mathrm{obs}}, \bar L_t\bigr)}.
$$

Both numerator and denominator condition on observed upstream values
(sequential MTP semantics, Díaz et al. 2023). Per-component intervention
broadcasts across periods.

### Stacked sandwich

$T \times K$ propensity blocks fit on **disjoint** time subsets but
**potentially overlapping** within-period components, so the bread is:
- Block-diagonal across time (independent fits at each period).
- Within a period, block-diagonal across components (each component
  fit independently via the chain rule).

Net: the full $\sum_{t,k} p_{t,k} \times \sum_{t,k} p_{t,k}$ propensity
bread is block-diagonal, and the variance engine sums $T \times K$
single-model `apply_model_correction()` calls. Same shape as Phase 8e
multivariate variance, just looped over time periods as well.

### Stabilization (extends Phase 10b)

Per-period AND per-component numerator models
$g_{t,k}(A_{t,k} \mid \bar A_{t-1}, V)$ drop $L$ from the conditioning
set. Default behaviour drops time-varying confounders only; custom
`numerator = ~ V` adds baseline conditioning.

### Effect modification (extends Phase 10c)

`A_k:modifier` interactions (where `A_k` is any of the joint treatment
components) strip from per-period × per-component propensity formulas
via the chunk 10c wiring. MSM expands to `Y ~ 1 + modifier` via
`build_ipw_msm_formula()` (one MSM total, not per period or per component).

## Items

### 19-trim — cross-cutting weight truncation ✅

- [x] `truncate_weights()` helper in `R/ipw_weights.R`.
- [x] `trim =` parameter on `contrast()`, threaded through all weight paths: point IPW, multivariate IPW, longitudinal IPW, point AIPW, longitudinal AIPW.
- [x] Sandwich variance: truncation threshold fixed at alpha_hat under numDeriv perturbation.
- [x] Bootstrap: each replicate recomputes its own quantile threshold.
- [x] lmtp cross-checks for multivariate and longitudinal with matching `.trim` values.
- [x] Vignette updates (ipw.qmd, longitudinal.qmd, diagnostics.qmd).

### 19a — core engine ✅

- [x] Lift the `causatr_longitudinal_multivariate_pending` rejection in `fit_longitudinal_ipw()`. Replaced with conditional rejections for stabilize+MV and EM+MV.
- [x] Refactor `fit_longitudinal_ipw()` to dispatch on `length(treatment)`: univariate path (chunk 10a) vs. multivariate path. Multivariate path nests `fit_treatment_models()` (Phase 8) inside the per-period loop so each period gets a list of K models.
- [x] Extend `compute_longitudinal_weights()` to call `compute_density_ratio_weights_mv()` per period and multiply across both axes.
- [x] Extend `make_weight_fn_longitudinal()` to build a stacked-α closure with $T \times K$ blocks (delegate the within-period part to `make_weight_fn_mv()`).
- [x] Extend `compute_ipw_if_self_contained_long_one()` to sum $T \times K$ single-model corrections (period-loop nests component-loop).
- [x] Bootstrap unchanged: id-clustered resampling already preserves both axes.
- [x] IPSI check in `compute_ipw_contrast_longitudinal()` handles MV (per-component check).
- [x] `build_longitudinal_ps_formula()` handles vector `treatment` (placeholder response for MV).

### 19b — stabilization ✅

- [x] Extend the per-period numerator-model sweep to be per-period × per-component when `stabilize = "marginal"` AND `length(treatment) > 1L`. Reuse the denominator `fit_treatment_models()` path with `confounders = remove_response(build_longitudinal_numerator_ps_formula(...))` so the marginal numerator drops time-varying L while keeping lags + V; the per-component chain rule then layers on the prior components $A_{1..k-1}$ automatically.
- [x] Stash $T \times K$ numerator models alongside $T \times K$ denominator models via the `numerator_models` / `stabilize` attributes on each period's denominator `causatr_treatment_models` list (read by the MV weight + variance closures).
- [x] Lift the `causatr_longitudinal_mv_stabilize_pending` rejection.

### 19c — effect modification

- [x] Reuse the chunk 10c per-period propensity stripping via `parse_effect_mod()`'s `confounder_terms` slot. The Phase 8b multivariate EM already handles per-component stripping; Phase 10c handles per-period. Composition is automatic if both wirings are honoured. Confirmed: the `causatr_longitudinal_mv_em_pending` rejection was the only blocker; lifting it makes the existing wiring (`em_info` → `details`, per-component `fit_treatment_models()` stripping, final-period `build_ipw_msm_formula()` expansion) compose with no MV-specific EM code.
- [x] Truth test on a 2-period × 2-component DGP with EM (`make_em_mv_long_scm()`, a binary MV analogue of `make_em_ice_scm()` with analytical static truths 9 and 15 by `sex`).

## Tests

- Truth-based: 2-period × 2-component continuous DGP. Cross-check IPW point estimate against ICE g-comp on the same DGP.
- Static + static recovers identical estimate as multivariate point IPW with `type = "point"` on the per-period collapsed data (sanity invariant).
- Stabilized + static identical to unstabilized + static (Phase 8e / 10b invariant carries over).
- EM cross-method agreement vs ICE on the 2-period × 2-component sex-stratified DGP.
- Bootstrap parity within MC tolerance.

## Out of scope

- Time-varying effect modification: deferred to Phase 18 (SNMs).
- IPSI per period (separate Phase 20).

## Downstream: longitudinal AIPW composition

When Phase 19 ships, longitudinal AIPW (Phase 16i, `fit_aipw_longitudinal()`)
will need a corresponding multivariate composition update:
`fit_aipw_longitudinal()` must dispatch to the multivariate per-period
density chain, and `ice_aipw_iterate()` must use multivariate cumulative
weights. The `causatr_longitudinal_multivariate_pending` rejection in
`R/aipw.R` gates this.

## References

Díaz I, Williams N, Hoffman KL, Schenck EJ (2023). Non-parametric
causal effects based on longitudinal modified treatment policies.
*JASA* 118:846–857.

Robins JM, Hernán MA, Brumback B (2000). Marginal structural models
and causal inference in epidemiology. *Epidemiology* 11(5):550–560.
