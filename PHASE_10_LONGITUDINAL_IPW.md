# Phase 10 — Longitudinal IPW

> **Status: chunk 10a DONE (2026-04-26); chunks 10b stabilization + 10c effect modification PENDING.**
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

## Items (chunks 10a, 10b, 10c)

### 10a — core engine (DONE 2026-04-26)

- [x] Remove `type == "longitudinal"` gate in `fit_ipw()` and `check_causat_inputs()` (matching stays rejected).
- [x] `fit_longitudinal_ipw()` in `R/longitudinal_ipw.R`: loop over `time_points` and fit one `causatr_treatment_model` per period via `fit_treatment_model()` on the period-k subset; per-period formula via `build_longitudinal_ps_formula()` mirrors `ice_build_formula()`'s lag handling.
- [x] `compute_longitudinal_weights()` in `R/ipw_weights.R`: per-id cumulative product weight via per-period reuse of `compute_density_ratio_weights()`. Sequential-MTP semantics — both numerator and denominator condition on observed upstream treatments.
- [x] `compute_ipw_contrast_longitudinal()` in `R/longitudinal_ipw.R`: per-intervention final-period weighted Hájek MSM (`Y ~ 1`).
- [x] `make_weight_fn_longitudinal()` in `R/ipw_weights.R`: stacked-α closure for the variance engine's cross-derivative.
- [x] `variance_if_ipw_longitudinal()` + `compute_ipw_if_self_contained_long_one()` in `R/variance_if.R`: K+1-block stacked sandwich with block-diagonal propensity bread (parallel to `compute_ipw_if_self_contained_mv_one()`'s shape, periods replacing components). Natural-course short-circuit when `alpha_hat_stacked = numeric(0)`.
- [x] `ipw_longitudinal_variance_bootstrap()` in `R/variance_bootstrap.R`: id-clustered bootstrap mirroring `ice_variance_bootstrap()`'s clone-and-reassign pattern. `refit_ipw()` is type-aware and threads `id` / `time` to the longitudinal fitter.
- [x] Sequential positivity warning `causatr_longitudinal_seq_positivity` (per-period max > 100; rate-limited).
- [x] Truth-based tests in `tests/testthat/test-longitudinal-ipw.R`: cross-method agreement vs ICE on continuous shift + binary static DGPs, sandwich vs bootstrap parity, lmtp::lmtp_sdr point-estimate cross-check (skipped if SuperLearner unavailable), positivity warning fires/doesn't-fire, natural course recovers observed marginal mean, all rejection paths under `_pending` classed errors.
- [x] Coverage matrix + NEWS + CLAUDE.md updates.

### 10b — stabilized weights (PENDING)

- [ ] Lift `causatr_longitudinal_stabilize_pending` and `causatr_longitudinal_numerator_pending` rejections.
- [ ] Per-period numerator models $g_k(A_k \mid \bar A_{k-1})$ (drops $L$); fit via a parallel sweep through `time_points` in `fit_longitudinal_ipw()` and stash in `attr(fit$details$treatment_models_by_time, "numerator_models")`.
- [ ] Extend `compute_longitudinal_weights()` to swap the per-period numerator under `stabilize = "marginal"`.
- [ ] Extend `make_weight_fn_longitudinal()` to precompute fixed-γ numerator vectors per period and route to a stabilized closure (γ held fixed under `numDeriv` perturbation; same nuisance-fixed convention as multivariate Phase 8e).
- [ ] `refit_ipw()` already replays `stabilize`; just confirm it threads through under `type = "longitudinal"`.
- [ ] Tests: numerator-model structure, stabilized + static recovers same point estimate as unstabilized + static, stabilized + shift truth recovery, bootstrap captures γ uncertainty.

### 10c — effect modification (PENDING)

- [ ] Lift `causatr_longitudinal_em_pending` rejection.
- [ ] Per-period propensity formula stripping via `parse_effect_mod()`'s `confounder_terms` slot in `fit_longitudinal_ipw()`.
- [ ] MSM expansion `Y ~ 1 + modifier` via the existing `build_ipw_msm_formula()`; nothing to change in the variance engine (the MSM is already built generically).
- [ ] Carry the Phase 6 baseline-modifier constraint as a doc-level note (modifier MUST be baseline; time-varying modifier under MSM violates Robins 2000).
- [ ] Truth-based test on `make_em_ice_scm()`-style DGP (sex-stratified ATE recovers analytical truth).

## Dependencies

Phase 4 (self-contained IPW engine), Phase 5 (ICE data structures, person-period conventions).

## Out of scope

- Longitudinal multivariate IPW (combine with Phase 8 multivariate)
- Grace period / visit-process interventions (future enhancement)
- Stratified ICE option (future enhancement)
