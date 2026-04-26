# Phase 10 — Longitudinal IPW

> **Status: chunks 10a + 10b + 10c DONE (2026-04-26). Phase 10 complete.**
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

### 10b — stabilized weights (DONE 2026-04-26)

- [x] Lift the chunk 10a `_stabilize_pending` / `_numerator_pending` rejections in `fit_longitudinal_ipw()`. Replaced with a `causatr_longitudinal_numerator_without_stabilize` guard for the (custom numerator + `stabilize = "none"`) combination.
- [x] Per-period numerator models $g_k(A_k \mid \bar A_{k-1}, V)$ via a second sweep through `time_points` in `fit_longitudinal_ipw()`. Stashed in `fit$details$numerator_models_by_time` (parallel to `treatment_models_by_time`). `build_longitudinal_numerator_ps_formula()` is the per-period formula builder; default behaviour drops time-varying confounders and keeps treatment lags (chain-rule validity), and `numerator = ~ V` adds the user-supplied baseline conditioning set on top of the lags.
- [x] `compute_longitudinal_weights()` accepts `numerator_models_by_time =` and routes per-period weights through new `compute_stabilized_period_weight()` helper when stabilized. Reuses the per-family branches of the unstabilized engine; only the density-evaluation model is swapped.
- [x] `make_weight_fn_longitudinal()` accepts `numerator_models_by_time =` and builds per-period stabilized closures via new `make_long_stabilized_period_closure()` helper, which delegates to the existing `mv_stabilized_closure()` for the actual fixed-γ / alpha-dependent weight (same primitive used by Phase 8e multivariate IPW).
- [x] `refit_ipw()` already replays `stabilize` and `numerator`; chunk 10a's type-aware threading covers the longitudinal case.
- [x] Tests in `test-longitudinal-ipw.R`: numerator-model structure under default and custom `numerator = ~ L0`, stabilized + static binary recovers identical estimate **and** SE as unstabilized (T-long-ipw-stab2), stabilized + shift cross-method agreement vs ICE on continuous DGP (T-long-ipw-stab3), bootstrap captures γ uncertainty (T-long-ipw-stab4).
- [x] Coverage matrix + NEWS + CLAUDE.md updates.

### 10c — effect modification (DONE 2026-04-26)

- [x] Lift the chunk 10a `causatr_longitudinal_em_pending` rejection in `fit_longitudinal_ipw()`. `check_em_compat()` upstream still rejects bare treatment in confounders (`~ L + A`) with `causatr_bare_treatment_in_confounders` — the only invalid use of treatment-touching terms.
- [x] Per-period propensity formulas strip `A:modifier` via the chunk 10a wiring (`em_info$confounder_terms` → `baseline_terms`); modifier main effects (`sex`) are retained in the propensity. No new code -- the existing chunk 10a path was already EM-aware.
- [x] MSM expansion via `build_ipw_msm_formula(outcome, em_info)` in `compute_ipw_contrast_longitudinal()`: replaces the hard-coded `Y ~ 1` with `Y ~ 1 + modifier` whenever `em_info$has_em` is `TRUE`. Variance engine unchanged -- the expanded MSM design matrix flows through `iv_design_matrix()` and `apply_model_correction()` transparently.
- [x] Doc-level baseline-modifier constraint (Robins 2000) carried in NEWS, PHASE_10 doc, and the longitudinal vignette.
- [x] Truth-based tests in `test-longitudinal-ipw.R`: T-long-ipw-em1 recovers ATE|sex=0 = 5 and ATE|sex=1 = 8 from `make_em_ice_scm()`; T-long-ipw-em2 verifies per-period propensity strips `A:sex` and the MSM expands to `Y ~ 1 + sex`; T-long-ipw-em3 cross-method agreement vs ICE EM.
- [x] Coverage matrix + NEWS + CLAUDE.md updates.

## Dependencies

Phase 4 (self-contained IPW engine), Phase 5 (ICE data structures, person-period conventions).

## Out of scope

- Longitudinal multivariate IPW (combine with Phase 8 multivariate)
- Grace period / visit-process interventions (future enhancement)
- Stratified ICE option (future enhancement)
