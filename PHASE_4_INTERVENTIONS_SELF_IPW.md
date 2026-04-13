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

The variance engine (`R/variance_if.R`) was already refactored in the Phase 2/3/4 variance pass (`VARIANCE_REFACTOR.qmd`) to plug self-contained IPW in without touching `variance_if()` itself. The dispatcher is live; only the body of `correct_propensity_self_contained()` needs to be filled in.

**Where to implement.** `R/variance_if.R`, function `correct_propensity_self_contained(fit, J, fit_idx, n_total)`. It currently aborts with a Phase 4 message. The surrounding dispatcher `correct_propensity()` already routes here when `fit$model` is not a `glm_weightit` — so the moment Phase 4's `fit_ipw()` produces a plain weighted GLM instead of calling `WeightIt::glm_weightit()`, this branch takes over automatically.

**What Phase 4 must add to `fit$details`.** The body needs three things beyond what Phase 3 stores:

1. `fit$details$propensity_model` — the treatment density model from `R/treatment_model.R` (a `glm` / multinom / etc.).
2. `fit$details$alpha_hat` — the propensity parameter vector \(\hat\alpha\).
3. `fit$details$weight_fn` — an R closure \(\alpha \mapsto w_i(\alpha)\) that reconstructs the length-`n_fit` weight vector given a candidate \(\alpha\). This closure captures the intervention, the observed treatment, and the treatment density functions so `numDeriv::jacobian()` can perturb \(\alpha\) alone.

Without these three slots, Branch B cannot be implemented.

**Implementation steps.** (Copied verbatim from `VARIANCE_REFACTOR.qmd` — the engine side is already designed around this shape.)

1. **MSM correction term.** Call
   ```r
   msm <- correct_model(fit$model, J, fit_idx, n_total)
   ```
   to get both `msm$correction` (the MSM Channel-2 contribution) *and* `msm$h = A_{ββ}^{-1} J`. The `$h` slot is why `correct_model()` was designed to return `h` alongside `$correction` in Phase A.

2. **Cross-derivative \(A_{\beta\alpha}\) via `numDeriv::jacobian`.** Build
   ```r
   psi_beta_bar <- function(alpha) {
     w      <- fit$details$weight_fn(alpha)           # length n_fit
     X      <- model.matrix(fit$model)
     fam    <- fit$model$family
     eta    <- as.numeric(X %*% beta_hat)
     mu     <- fam$linkinv(eta)
     mu_eta <- fam$mu.eta(eta)
     r      <- Y_fit - mu
     # Weighted average score at the MSM
     as.numeric(crossprod(X, w * mu_eta * r / fam$variance(mu))) / n_fit
   }
   # A_{βα} = -(1/n) sum_i d psi_β,i / d alpha^T; numDeriv gives +d(average)/d alpha,
   # so flip the sign.
   A_beta_alpha <- -numDeriv::jacobian(psi_beta_bar, x = fit$details$alpha_hat)
   ```
   The numerical jacobian is the right tool here (not an analytic formula) because `weight_fn` changes per intervention type — shift, MTP, IPSI, dynamic rule — and re-deriving the closed form for each is tedious and error-prone. `numDeriv::jacobian` handles all of them with one code path.

3. **Propensity correction via a derived gradient.** Form the derived gradient
   \[
   g^{\mathrm{prop}} = A_{\beta\alpha}^T h_{\mathrm{msm}}
   \]
   (a `q`-vector) and feed it back into the same `correct_model()` primitive, this time on the **propensity** model:
   ```r
   g_prop <- as.numeric(crossprod(A_beta_alpha, msm$h))
   prop <- correct_model(
     fit$details$propensity_model,
     g_prop,
     fit_idx_prop,
     n_total
   )
   Ch2 <- msm$correction + prop$correction
   ```
   `fit_idx_prop` is the row map for the propensity fit (usually identical to `fit_idx` if the MSM and propensity share a fit frame).

4. **Return.** `Ch2` is the full combined MSM + propensity correction per individual — the same shape returned by Branch A.

**Tests to add** (from `VARIANCE_REFACTOR.qmd`, "correct_propensity() tests"):

- **T1 — Branch A ≡ Branch B** on a static binary IPW with logistic propensity. Fit the same setup two ways: (a) via `WeightIt::weightit(method = "glm")` → Branch A; (b) via a manually fitted logistic propensity + manually computed weights + weighted `stats::glm` MSM → Branch B. Both paths must yield the same `IF_mu`, SE, and vcov to `~1e-6` absolute tolerance (the `numDeriv` Richardson accuracy floor).
- **T2 — Numerical vs analytic \(A_{\beta\alpha}\)** on the same binary static setup. Derive \(\partial\psi_{\beta,i}/\partial\alpha^T\) by hand (closed form is tractable for logistic propensity + binary static weights) and compare against `numDeriv::jacobian(psi_beta_bar, ...)` elementwise at `~1e-6`.
- **T3 — End-to-end stacked sandwich** on a tiny simulated dataset (`n = 200`, `q = p = 2`, Gaussian outcome, logistic propensity). Assemble the full \((q+p+1)\times(q+p+1)\) stacked bread/meat by hand, invert, read off the \((\mu,\mu)\) entry, and compare against `variance_if()` at `~1e-10`.
- **T4 — Non-static IPW** (shift / MTP) on NHEFS where Channel 1 ≠ 0: IF-based variance must materially exceed `J V_\beta J^T` (the delta-only shortcut) — this proves Phase 4 is actually capturing the extra uncertainty.

**Why `numDeriv::jacobian` here and not analytic.** Each intervention type (`static`, `shift`, `scale_by`, `threshold`, `dynamic`, `ipsi`) has a different `weight_fn`, so the analytic derivative would be a different ~10-line block per type, all of which need to agree with the engine's expectations. One `numDeriv::jacobian` call on the closure handles every intervention shape uniformly, at Richardson-extrapolation precision (`~1e-7`). Test T2 verifies the numerical derivative against a closed form in the one case we can fully hand-derive.

- **Bootstrap:** unchanged — `variance_bootstrap()` refits the whole pipeline on resampled data and never touches `variance_if()` / `correct_propensity()`.

## Items

- [ ] `R/treatment_model.R` — fit treatment density model + density evaluation
- [ ] `R/ipw_weights.R` — density ratio weight computation
- [ ] Update `R/ipw.R` — self-contained mode for non-static interventions. Must populate `fit$details$propensity_model`, `fit$details$alpha_hat`, and `fit$details$weight_fn` (a closure `alpha -> w_i(alpha)`) so the existing `correct_propensity()` Branch B dispatcher can take over.
- [ ] Fill in the body of `correct_propensity_self_contained()` in `R/variance_if.R` per the steps in §7 above. The dispatcher in `correct_propensity()` already routes non-`glm_weightit` fits here — no engine-level changes required.
- [ ] Categorical treatment support in checks + across all methods
- [ ] IPSI implementation using treatment density model
- [ ] Tests T1–T4 in `tests/testthat/test-variance-if.R` — Branch A ≡ Branch B, analytic vs numerical A_βα, end-to-end stacked sandwich, non-static IPW where Ch.1 ≠ 0
- [ ] Vignette: `interventions.qmd` — shift, scale, threshold, dynamic, IPSI examples
- [ ] IPW for time-varying treatments — delegate to `WeightIt::weightitMSM()` (deferred from Phase 5)
