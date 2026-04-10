# Phase 4 — Extended Interventions + Self-Contained IPW Engine

> **Status: PENDING**
> Book chapters: 12 (density ratio weights), 19 (dynamic strategies), Díaz et al. 2023 (MTPs), Kennedy 2019 (IPSI)

## Scope

Self-contained IPW engine supporting dynamic/MTP/IPSI interventions, categorical treatment support, and the interventions vignette.

## Why self-contained IPW?

WeightIt estimates standard propensity-score weights `1/f(A|L)` — valid for static interventions only. For non-static interventions we need **density ratio weights**:

```
w_i = f(a_intervened_i | L_i) / f(a_observed_i | L_i)
```

This requires:
1. Fitting a treatment density model `f(A | L)` ourselves
2. Evaluating that density at both the observed and intervened treatment values
3. Computing the ratio

WeightIt doesn't provide this capability.

## Plan

### 1. Treatment density model (`R/treatment_model.R` — new file)

Fit a treatment model and expose the density evaluation:

```r
# For binary treatment: logistic regression → Bernoulli density
# For continuous treatment: GLM for mean + residual SD → normal density
# For categorical treatment: multinomial → categorical density
fit_treatment_model <- function(data, treatment, confounders, family, ...)
evaluate_density <- function(treatment_model, treatment_values, newdata)
```

### 2. Density ratio weights (`R/ipw_weights.R` — new file)

```r
compute_density_ratio_weights <- function(
  treatment_model,
  data,                # original data
  data_intervened,     # data after applying intervention
  treatment            # treatment column name
)
# Returns: w_i = f(a_intervened_i | L_i) / f(a_observed_i | L_i)
```

### 3. Self-contained IPW estimation

Update `fit_ipw()` to support a `self_contained = TRUE` mode (or detect automatically when interventions are non-static):

```r
# Static: delegate to WeightIt (current behavior)
# Non-static: use our density ratio engine
```

For the MSM with non-static interventions, the weighted estimator is:
```
Ê[Y^{d(a,l)}] = Σᵢ wᵢ Yᵢ / Σᵢ wᵢ
```
where wᵢ = f(d(Aᵢ,Lᵢ) | Lᵢ) / f(Aᵢ | Lᵢ).

### 4. Supported intervention × method combinations after Phase 4

| Intervention | G-comp | IPW (WeightIt) | IPW (self-contained) | Matching |
|---|---|---|---|---|
| `static()` | ✓ | ✓ | ✓ | ✓ |
| `shift()` | ✓ | — | ✓ | — |
| `scale_by()` | ✓ | — | ✓ | — |
| `threshold()` | ✓ | — | ✓ | — |
| `dynamic()` | ✓ | — | ✓ | — |
| `ipsi()` | — | — | ✓ | — |

### 5. Categorical treatment support

- G-comp: already works (factor treatment column in the GLM)
- IPW: `WeightIt::weightit()` handles categorical via multinomial PS
- Matching: `MatchIt::matchit()` handles categorical (multi-arm)
- Need: update `check_estimand_trt_compat()` to allow categorical + ATE

### 6. IPSI (incremental propensity score interventions)

The `ipsi(delta)` constructor already exists. Implementation:
- Multiply the odds of treatment by δ: `p_new = (δ·p) / (1 - p + δ·p)`
- Requires the propensity score → treatment density model
- Compute density ratio weights at the modified probability

### 7. Variance for self-contained IPW

- **Sandwich:** Stacked estimating equations: treatment model score + weighted outcome estimating equation. Can use `geex` or manual bread/meat.
- **Bootstrap:** Same full-pipeline resample as current.

## Items

- [ ] `R/treatment_model.R` — fit treatment density model + density evaluation
- [ ] `R/ipw_weights.R` — density ratio weight computation
- [ ] Update `R/ipw.R` — self-contained mode for non-static interventions
- [ ] Categorical treatment support in checks + across all methods
- [ ] IPSI implementation using treatment density model
- [ ] Sandwich variance for self-contained IPW (treatment model uncertainty)
- [ ] Vignette: `interventions.qmd` — shift, scale, threshold, dynamic, IPSI examples
- [ ] IPW for time-varying treatments — delegate to `WeightIt::weightitMSM()` (deferred from Phase 5)
