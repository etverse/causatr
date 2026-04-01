# Phase 5 — Longitudinal ICE G-Computation

> **Status: PENDING**
> Book chapters: 19, 20, 21
> Key paper: Zivich et al. (2024), Statistics in Medicine 43:5562–5572

## Scope

ICE (iterated conditional expectation) g-computation for time-varying treatments, sandwich variance via stacked estimating equations, censoring, time-varying IPW delegation.

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

- [ ] `R/ice.R` — ICE g-computation engine (backward iteration)
- [ ] `causat()` longitudinal path: detect `id` + `time`, call `fit_ice()`
- [ ] Sandwich variance for ICE via stacked estimating equations (`geex` package)
- [ ] Censoring handling within ICE (restrict to uncensored at each backward step)
- [ ] External IPCW weights via `weights` argument for longitudinal
- [ ] IPW for time-varying treatments — delegate to `WeightIt::weightitMSM()`
- [ ] Sequential positivity warnings in `causat()`
- [ ] Stratified ICE option (`causat(..., stratified = TRUE)`)
- [ ] Vignette: `longitudinal.Rmd` (include Table 20.1 treatment-confounder feedback demo)

## Implementation guide

> To be created from claude.ai summarizing Ch. 19–21 when this phase begins.
