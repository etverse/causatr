# Phase 5 — Longitudinal ICE G-Computation

> **Status: DONE**

## Scope

ICE (iterated conditional expectation) g-computation for time-varying treatments, sandwich variance via stacked estimating equations, censoring, parallel bootstrap.

## Key algorithm: ICE backward iteration

```
For strategy Ā* = (a₀*, a₁*, ..., a*_K):

1. Fit outcome model at final time τ:
   E[Y_τ | Ā_{τ-1}, L̄_{τ-1}] among uncensored at τ

2. Predict under intervention Ā* (keep observed L̄):
   → get Ŷ*_{τ-1}

3. Fit pseudo-outcome model at τ-1:
   E[Ŷ*_{τ-1} | Ā_{τ-2}, L̄_{τ-2}] among uncensored at τ-1

4. Predict under intervention → get Ŷ*_{τ-2}

5. Repeat backward to time 0

6. μ̂_τ(Ā*) = mean(Ŷ*_0)
```

## Items

- [x] `R/ice.R` — ICE g-computation engine (backward iteration via `fit_ice()` + `ice_iterate()`)
- [x] `causat()` longitudinal path: detect `id` + `time`, dispatch to `fit_ice()`
- [x] Sandwich variance for ICE via stacked estimating equations (manual influence functions, no `geex` dependency)
- [x] Censoring handling within ICE (restrict to uncensored at each backward step)
- [x] External IPCW weights via `weights` argument for longitudinal
- [x] Bootstrap variance for ICE (resample individuals, parallel via `boot::boot`)
- [x] Dynamic interventions for longitudinal data (static, shift, scale, threshold, dynamic, NULL)
- [x] `parallel` / `ncpus` arguments for `contrast()` bootstrap
- [x] Vignette: `longitudinal.qmd` (Table 20.1 treatment-confounder feedback demo)

## Deferred to other phases

- IPW for time-varying treatments (`WeightIt::weightitMSM()`) $\to$ Phase 4
- Sequential positivity warnings $\to$ Phase 11 (diagnose rewrite)
- Stratified ICE option $\to$ Phase 15
