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
- [x] **Fix: ICE sandwich under time-varying covariate missingness.** With MCAR
  missingness in a time-varying confounder, the sandwich aborted with a raw
  `subscript out of bounds`: the cascade gradient built each step's
  counterfactual design via `iv_design_matrix()` → `model.matrix()`, which
  silently `na.omit`-drops missing-term rows, desyncing the design from the
  id-indexed `match()` vectors in `variance_if_ice_chain()`. Fix:
  `ice_model_complete_rows()` restricts each step's prediction frame to the rows
  the step model can score — exactly the individuals present in that step's
  gated estimating equation (the `I(C_{i,k}=0)` factor of the ICE stacked score,
  Zivich et al. 2024, *Stat. Med.* 43:5562). The baseline pseudo-outcome stays
  defined for everyone (point estimate and Channel 1 unchanged), and with no
  missing data the frame is byte-identical to before. Validity is **MCAR-only**;
  under MAR use `causat_mice()` or IPCW. Shared by the stratified-ICE chain.
  Validated: sandwich matches bootstrap + delete-one-id jackknife (~1.0) on
  2-period and 3-period (intermediate-NA) panels (`test-missing-data.R`).

## Deferred to other phases

- IPW for time-varying treatments (`WeightIt::weightitMSM()`) $\to$ Phase 4
- Sequential positivity warnings $\to$ Phase 11 (diagnose rewrite)
- Stratified ICE option $\to$ Phase 15
