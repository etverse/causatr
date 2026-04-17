# Phase 2 — Point-Treatment G-Computation + Inference

> **Status: DONE**
> Book chapters: 11, 13

## Scope

G-computation (parametric g-formula) for point treatments with sandwich and bootstrap inference. Binary + continuous outcomes, binary + continuous treatments, all intervention types.

## What was implemented

1. `causat(estimator = "gcomp")` $\to$ `fit_gcomp_point()` in `R/gcomp.R`
   - Builds formula from `treatment` + `confounders` term labels
   - Fits on uncensored + non-NA-outcome rows
   - Calls `model_fn(formula, data, family, weights, ...)` — user-pluggable (default `stats::glm`)
   - Stores `fit_rows` logical vector for variance alignment

2. `contrast()` $\to$ `compute_contrast()` in `R/contrast.R`
   - Creates counterfactual datasets via `apply_intervention()`
   - Predicts under each intervention, averages over target population
   - Excludes NA predictions from standardisation average
   - Pairwise contrasts: difference, ratio (delta method gradient), odds ratio (delta method gradient)
   - `get_target_idx()`: ATE (all rows), ATT (A==1), ATC (A==0), subset (quoted expression)

3. Sandwich variance in `R/variance_sandwich.R`
   - G-comp: `sandwich::sandwich(model)` for $V_\beta$ (Huber--White HC0)
   - Propagated via `compute_vcov_marginal()`: $J V_\beta J^\top$ using `numDeriv::jacobian()`

4. Bootstrap variance in `R/variance_bootstrap.R`
   - Full-pipeline resampling via `boot::boot()`
   - `refit_gcomp()`: resample $\to$ filter censored/NA $\to$ refit model_fn

## Key decisions made during implementation

- **`model_fn` parameter** instead of hardcoded GLM/GAM detection
- **ci_method: "sandwich" and "bootstrap" only** — delta method applied internally for ratio/OR contrasts, not a separate ci_method
- **Sandwich via $J V_\beta J^\top$** (marginaleffects approach) rather than stacked estimating equations (Zivich et al.). Asymptotically equivalent; simpler to implement. Stacked EE via `geex` planned for Phase 5 (ICE).
- **NA handling**: rows with missing confounders produce NA predictions $\to$ excluded from target population automatically
- **`numDeriv` added to Imports** for Jacobian computation

## Testing targets

| Test | Expected | Status |
|---|---|---|
| NHEFS g-formula ATE | $\approx$ 3.5 kg (CI: 2.6–4.5) | Passing |
| Simulated DGP ATE = 3 | $\approx$ 3 | Passing |
| Binary outcome RD $\approx$ 0.33 | $\approx$ 0.33 | Passing |
| Risk ratio (binary outcome) | > 1 | Passing |
| Shift intervention (dY/dA = 2) | $\approx -2$ | Passing |
| Scale, threshold, dynamic | Correct direction | Passing |
| Sandwich SE $\approx$ bootstrap SE | Within factor of 2 | Passing |
| GAM via model_fn = mgcv::gam | $\approx$ 3 | Passing |
| Censoring handling | $\approx$ 3 (random censoring) | Passing |
| ATT estimand | Valid, non-NA | Passing |
| Subset estimand | $\approx$ 3 (constant effect DGP) | Passing |
