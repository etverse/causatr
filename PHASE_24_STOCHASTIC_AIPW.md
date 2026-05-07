# Phase 24 — Stochastic Interventions under IPW and AIPW

> **Status: PENDING** (design doc)
>
> **Depends on:** Phase 12 (stochastic gcomp — done), Phase 16 (AIPW — done),
> Phase 4 (self-contained IPW — done), Phase 8 (multivariate IPW — done),
> Phase 10 (longitudinal IPW — done)
>
> **Composes with:** Phase 14 (IPCW, future), Phase 17 (transportability, future)

## Motivation

Phase 12 delivered stochastic interventions under g-computation only.
The `stochastic()` constructor stores a `sampler` — a function that
draws random treatment values — but no density function. Without a
density, the IPW and AIPW engines cannot compute the density-ratio
weight $w_i = g^*(A_i \mid L_i) \,/\, \hat g(A_i \mid L_i)$ needed
for importance weighting.

Adding an optional `density` parameter to `stochastic()` enables two
new estimators:

1. **Stochastic IPW** — importance-weighted Hájek mean under a
   user-specified stochastic policy. Originally deferred as "Phase 12b"
   in `PHASE_12_STOCHASTIC.md`.
2. **Stochastic AIPW** — doubly-robust estimation with MC-averaged
   outcome predictions and density-ratio augmentation weights. Originally
   planned as chunk 16n in `PHASE_16_AIPW.md`.

The shared `density` parameter, the overlapping weight-computation
machinery, and the testing surface justify a unified phase.

### Why the user must supply the density

The density-ratio weight $w_i = g^*(A_i \mid L_i) / \hat g(A_i \mid L_i)$
requires evaluating the intervention density $g^*$ at the observed
treatment. The propensity model $\hat g$ is already fitted, but $g^*$
is defined by the user's stochastic policy. Unlike MTPs (where the
density ratio has an analytical Jacobian from the transformation) or
IPSI (where the closed-form weight collapses directly), general
stochastic interventions have no automatic mechanism to recover the
density from the sampler alone. The user must therefore supply both:

- `sampler`: draws counterfactual treatments (for MC g-formula)
- `density`: evaluates the intervention density at observed values
  (for IPW / AIPW weights)

These are two views of the same distribution. If they are
inconsistent, the estimator is biased — this is a user responsibility
that cannot be validated automatically.

### Double robustness

Díaz & van der Laan (2012, §3) derive the augmented IPTW (A-IPTW)
estimator for stochastic interventions and prove it solves the
efficient influence curve equation, yielding:

- **Consistency** if either the outcome model or the density ratio is
  correct (double robustness).
- **Semiparametric efficiency** when both are correct.

**Caveat (Wen, Marcus & Young 2023):** double robustness can fail
when the intervention density depends on the *observed* treatment
process. For `stochastic()` where $g^*(a \mid L)$ does not depend on
observed $A$, standard DR holds. The design doc assumes this condition.

## Scope

1. `stochastic()` constructor extension with optional `density` parameter.
2. Stochastic IPW (point): density-ratio weights, all treatment types.
3. Stochastic AIPW (point): MC-averaged outcome predictions + density-ratio
   augmentation. DR tests. Efficiency tests.
4. Stochastic IPW (multivariate): per-component densities, product weight.
5. Stochastic AIPW (multivariate): combines point AIPW + multivariate weights.
6. Stochastic AIPW (longitudinal): per-period densities, cumulative weights,
   ICE-AIPW backward recursion.
7. Stochastic IPW (longitudinal): per-period density-ratio weights.
8. Combined multivariate + longitudinal stochastic.
9. Sandwich variance: MC-averaged outcome gradient + numDeriv on
   density-ratio weight.
10. Bootstrap variance.
11. Documentation and test coverage.

## Non-scope

- **TMLE / SDR with cross-fitting.** `lmtp` covers this with a
  fundamentally different estimation strategy (ML nuisances,
  cross-fitting folds, targeting step).
- **Machine-learning nuisances.** AIPW with ML nuisances requires
  cross-fitting to recover $\sqrt{n}$-consistency.
- **Matching + stochastic.** Architecturally incompatible.
- **Auto-validation of sampler/density consistency.** User responsibility.
  We validate basic properties (non-negative, finite, correct length)
  but cannot check distributional equivalence.
- **Stochastic interventions without `density` under IPW/AIPW.**
  Density-free stochastic remains gcomp-only (Phase 12).
- **Optimal treatment regime estimation (DTR).** Different problem class.

## Caveats and disadvantages

1. **User burden.** The user must supply both `sampler` and `density`
   for the same distribution. If they are inconsistent, the estimator
   is silently biased with no built-in diagnostic.
2. **Tractability.** Not all stochastic policies have tractable
   densities. Rejection-sampled policies, mixture policies with
   intractable normalisation, or policies defined by simulation do
   not have closed-form densities. In those cases, AIPW/IPW are
   unavailable and the user must use gcomp (Phase 12).
3. **Niche audience.** Stochastic interventions are already
   specialised; the density requirement narrows the user base further.
4. **DR caveat.** DR fails if $g^*$ depends on observed $A$
   (Wen et al. 2023, Theorem 1). Document this clearly in roxygen
   and vignette.

---

## Design

### 1. Constructor extension

```r
stochastic <- function(sampler, n_mc = 100L, density = NULL) {
  # ... existing sampler / n_mc validation ...
  if (!is.null(density)) {
    if (!is.function(density)) {
      rlang::abort("`density` must be a function with signature function(a, data).")
    }
  }
  new_causatr_intervention(
    "stochastic",
    list(sampler = sampler, n_mc = n_mc, density = density)
  )
}
```

**Stored fields:** `type = "stochastic"`, `sampler`, `n_mc`, `density`.

**Helper:** `has_stochastic_density(iv)` returns `TRUE` if `iv$type == "stochastic"` and `iv$density` is non-NULL.

**Density function signature:** `function(a, data)` where:
- `a`: numeric vector of treatment values (length = `nrow(data)`)
- `data`: `data.table` with all covariates available

Returns a numeric vector of the same length as `a` with the
density/PMF values $g^*(a_i \mid L_i)$. Must be non-negative and
finite.

### 2. Lifting the IPW/AIPW rejection

Modify `check_intervention_family_compat()` in `R/ipw_weights_mv.R`
(currently lines 1051–1072). Instead of unconditionally aborting on
`iv_type == "stochastic"`, check for density:

```r
if (iv_type == "stochastic") {
  if (is.null(intervention$density)) {
    rlang::abort(c(
      paste0(
        "`stochastic()` interventions without a `density` function ",
        "are only supported under `estimator = 'gcomp'`."
      ),
      i = paste0(
        "Supply `density = function(a, data) ...` in `stochastic()` ",
        "to enable IPW and AIPW estimation."
      ),
      i = paste0(
        "The density function must evaluate the intervention density/PMF ",
        "g*(a | L) at the observed treatment values."
      )
    ))
  }
  # Fall through to new stochastic density-ratio branch below
}
```

### 3. Density-ratio weight computation

New branch in `compute_density_ratio_weights()` (R/ipw_weights.R),
after the existing `scale` branch:

```r
if (iv_type == "stochastic") {
  density_fn <- intervention$density
  f_int <- density_fn(a_obs, fit_data)   # user-supplied numerator
  # Validate: non-negative, finite, correct length
  stopifnot(length(f_int) == n_fit)
  stopifnot(all(is.finite(f_int)))
  stopifnot(all(f_int >= 0))
  f_obs <- evaluate_density(treatment_model, a_obs, fit_data)
  check_density_positivity(f_obs, "stochastic density ratio")
  return(f_int / f_obs)
}
```

Corresponding closure for sandwich in `make_weight_fn()`:

```r
if (iv_type == "stochastic") {
  density_fn <- intervention$density
  f_star <- density_fn(a_obs, fit_data)  # fixed numerator
  return(function(alpha) {
    eta <- as.numeric(X_prop %*% alpha)
    f_hat <- density_at_alpha(family_tag, a_obs, eta, sigma)
    f_star / f_hat
  })
}
```

The numerator `f_star` is held fixed (user-supplied, not estimated),
so the sandwich cross-derivative $\partial w_i / \partial \alpha$
acts only through the denominator — exactly the same structure as
the existing shift/scale branches.

### 4. The stochastic AIPW functional

$$
\hat\psi_{\mathrm{AIPW}}(g^*) = \frac{1}{n}\sum_{i=1}^{n}\left[
  \underbrace{\frac{1}{M}\sum_{m=1}^{M}\hat{m}(A_{i,m}^*, L_i)}_{\text{MC-averaged gcomp}}
  + \underbrace{\frac{g^*(A_i \mid L_i)}{\hat g(A_i \mid L_i)}}_{\text{density-ratio weight}}
  \cdot \underbrace{(Y_i - \hat{m}(A_i, L_i))}_{\text{residual at observed trt}}
\right]
$$

- **First term** — MC-averaged outcome-model predictions under the
  stochastic policy. Same as Phase 12 gcomp. Uses `predict_under_intervention()`
  which already handles the MC loop.
- **Second term** — density-ratio weight evaluated at the *observed*
  treatment. No MC averaging needed — the weight is deterministic
  given the propensity model and user-supplied density.
- **Third term** — outcome residual at observed treatment. Shared
  across interventions, computed once.

Implementation in `compute_aipw_contrast_point()`: replace the
single-shot prediction path (lines 438–444 of R/aipw.R) with
`predict_under_intervention()` when the intervention is stochastic:

```r
if (is_ipsi) {
  data_a <- fit_data
  preds_g <- preds_obs
} else if (has_stochastic_component(iv)) {
  pui <- predict_under_intervention(outcome_model, fit_data, treatment, iv)
  data_a <- pui$data_a   # last MC draw (for design matrix in variance)
  preds_g <- pui$preds   # MC-averaged predictions
} else {
  data_a <- apply_intervention(fit_data, treatment, iv)
  preds_g <- stats::predict(outcome_model, newdata = data_a, type = "response")
}
```

The density-ratio weight computation proceeds unchanged — the existing
`compute_density_ratio_weights()` call now falls through to the new
stochastic branch instead of aborting.

### 5. Sandwich variance

The stacked estimating-equation system is the same as standard AIPW
(outcome model + propensity model + plug-in), but Channel 2a uses
MC-averaged gradients.

**Channel 1 (direct term):**
$$
\mathrm{Ch1}_i = \frac{1}{M}\sum_m \hat m(A_{i,m}^*, L_i) + w_i \cdot (Y_i - \hat m(A_i, L_i)) - \hat\psi
$$
The MC-averaged prediction is already computed in `preds_g`, so
Channel 1 is structurally identical to the deterministic case:
`aipw_contrib = preds_g + w_iv * resid_obs`.

**Channel 2a (outcome model correction):**

The gradient $J_\beta = \partial \hat\psi / \partial \beta$ has two
parts. For the gcomp prediction term:
$$
\frac{\partial}{\partial\beta}\left[\frac{1}{n}\sum_i \frac{1}{M}\sum_m \hat m(A_{i,m}^*, L_i)\right]
= \frac{1}{nM}\sum_i\sum_m X_{i,m}^{*\top} \mu'(\eta_{i,m}^*)
$$
This requires MC averaging the design-matrix contribution — the same
pattern already implemented for gcomp stochastic in
`R/variance_if.R:490–517`. For the residual term (evaluated at
observed treatment), no MC averaging is needed.

Implementation in `variance_if_aipw.R` Channel 2a: when
`has_stochastic_component(iv)`, loop over `n_mc` draws to compute
`grad_a` (MC-averaged). `grad_b` (residual correction at observed
treatment) is unchanged.

**Channel 2b (propensity model correction):**

The `aug_mean_alpha` closure calls `weight_fn(alpha)` internally.
The new `make_weight_fn()` stochastic branch handles the
alpha-parameterised weight correctly (fixed numerator, parameterised
denominator). `numDeriv::jacobian` on this closure gives the
cross-derivative. No structural change needed.

**MC-averaged IF linearity (Stefanski & Boos 2002):**

The MC-averaged estimator $\hat\psi = (1/M) \sum_m \hat\psi_m$ has
IF equal to the average of per-draw IFs by linearity:
$\mathrm{IF}_i = (1/M) \sum_m \mathrm{IF}_{m,i}$. The MC-averaged
gradient is the gradient of the MC-averaged prediction, so Channel 2
reduces to a single gradient computation on the averaged quantities.
No fundamental new machinery.

### 6. Bootstrap variance

No special treatment needed. Each bootstrap replicate independently:
1. Resamples individuals.
2. Refits outcome and propensity models.
3. Runs the MC g-formula (fresh draws from sampler).
4. Computes density-ratio weights (density evaluated at resampled observed A).
5. Returns the AIPW functional.

Bootstrap captures both model uncertainty and MC noise.

### 7. Multivariate stochastic

For multivariate treatment $(A_1, A_2)$, each component gets its own
density:
```r
list(
  A1 = stochastic(sampler1, n_mc = 100, density = density1),
  A2 = stochastic(sampler2, n_mc = 100, density = density2)
)
```

The product weight factorises as:
$$
w_i = \prod_{k=1}^{K} \frac{g_k^*(A_{k,i} \mid A_{<k,i}, L_i)}{\hat g_k(A_{k,i} \mid A_{<k,i}, L_i)}
$$

This fits the existing `compute_density_ratio_weights_mv()` structure,
which applies per-component weights and multiplies. The per-component
stochastic branch returns `density_k(a_k, data) / evaluate_density(tm_k, a_k, data)`.

**Note:** the user-supplied density for component $k$ may condition
on upstream treatments $A_{<k}$. The `data` argument passed to the
density function includes all columns, so the user can access
`data$A1` when evaluating the density for `A2`. This mirrors the
conditioning structure in the existing sequential factorisation.

### 8. Longitudinal stochastic

At each time period $t$, the density-ratio weight is:
$$
w_{i,t} = \frac{g_t^*(A_{i,t} \mid \bar{L}_{i,t})}{\hat g_t(A_{i,t} \mid \bar{L}_{i,t})}
$$

The cumulative weight is $W_{i,t} = \prod_{s=1}^{t} w_{i,s}$.

The ICE-AIPW backward recursion uses:
- MC-averaged predictions at each step (already implemented for ICE
  gcomp in Phase 12).
- Cumulative density-ratio weights in the augmentation term at each step.

**Per-period density functions:** the `stochastic()` intervention at each
time point carries its own density. The longitudinal engine iterates
over periods and evaluates each period's density at the observed
treatment for that period.

### 9. Treatment types

| Treatment type | Density form | Weight computation | AIPW? | IPW? |
|---|---|---|---|---|
| Binary | $g^*(1 \mid L)$ for $A=1$, $1 - g^*(1 \mid L)$ for $A=0$ | `density(a, data) / dbinom(a, 1, p_hat)` | ✓ | ✓ |
| Continuous | $g^*(a \mid L)$ (PDF) | `density(a, data) / dnorm(a, mu_hat, sigma_hat)` | ✓ | ✓ |
| Categorical | $g^*(a \mid L)$ (PMF over $K$ levels) | `density(a, data) / predict(multinom, data)[a]` | ✓ | ✓ |
| Count | $g^*(a \mid L)$ (PMF, e.g. Poisson) | `density(a, data) / dpois(a, lambda_hat)` | ✓ | ✓ |

### 10. Estimand compatibility

- **ATE:** yes.
- **ATT/ATC:** `check_estimand_intervention_compat()` already rejects
  non-static interventions for ATT/ATC. Stochastic is non-static, so
  ATT/ATC are rejected. No change needed.
- **`by`-stratified / `subset`:** yes — restricts which rows are averaged.

---

## Chunks

| Chunk | Scope | Depends on | Status |
|---|---|---|---|
| 24a | `stochastic()` constructor: add optional `density` param, validation, `has_stochastic_density()` helper, update roxygen | Phase 12 | pending |
| 24b | Stochastic IPW (point): new branch in `compute_density_ratio_weights()` + `make_weight_fn()`, lift rejection in `check_intervention_family_compat()` when density present. Binary + continuous + categorical + count. Truth-based tests. | 24a, Phase 4 | pending |
| 24c | Stochastic AIPW (point): MC-averaged preds via `predict_under_intervention()` + density-ratio weights in `compute_aipw_contrast_point()`. DR tests (wrong outcome / wrong propensity). Efficiency test (AIPW SE ≤ gcomp SE, AIPW SE ≤ IPW SE). | 24a, 24b, Phase 16 | pending |
| 24d | Stochastic IPW (multivariate): per-component densities in `compute_density_ratio_weights_mv()`. Product weight. Truth-based tests. | 24b, Phase 8 | pending |
| 24e | Stochastic AIPW (multivariate): composes 24c + 24d. DR tests on multivariate DGP. | 24c, 24d | pending |
| 24f | Stochastic AIPW (longitudinal): per-period densities in `ice_aipw_iterate()`, cumulative weights, MC-averaged predictions at each ICE step. | 24c, Phase 10 | pending |
| 24g | Stochastic IPW (longitudinal): per-period density-ratio weights in `compute_density_ratio_weights_longitudinal()`. Cumulative product. | 24b, Phase 10 | pending |
| 24h | Combined multivariate + longitudinal stochastic (IPW + AIPW). | 24d, 24f or 24g | pending |
| 24i | Sandwich variance: MC-averaged outcome gradient in `variance_if_aipw.R` Channel 2a; numDeriv on density-ratio `make_weight_fn()` in Channel 2b. Point + multivariate + longitudinal. | 24c, 24f | pending |
| 24j | Bootstrap variance: verify both nuisances refit, MC draws fresh. Point + multivariate + longitudinal. | 24c | pending |
| 24k | Documentation: update `stochastic()` roxygen, `FEATURE_COVERAGE_MATRIX.md`, `CLAUDE.md`, vignette examples. | 24a–24j | pending |

## Invariants

- When `density = NULL`, `stochastic()` MUST be rejected by IPW and
  AIPW (existing behaviour, preserved). Error message updated to
  mention the `density` parameter.
- When `density` is non-NULL, `stochastic()` MUST be accepted by IPW
  and AIPW; the density-ratio weight
  $w_i = \texttt{density}(A_i, \texttt{data}_i) / \hat g(A_i \mid L_i)$
  MUST be computed without error.
- The user-supplied density MUST be validated for non-negativity and
  finiteness. The denominator (fitted density) is guarded by
  `check_density_positivity()`.
- Under correct specification of both nuisances AND a consistent
  sampler/density pair, AIPW, IPW, and gcomp point estimates MUST
  agree up to MC noise. This is a triangulation regression invariant.
- DR property: with a correct propensity model, AIPW MUST be
  consistent even when the outcome model is wrong. With a correct
  outcome model, AIPW MUST be consistent even when the propensity
  model is wrong. (Assumes $g^*$ does not depend on observed $A$.)
- MC integration applies ONLY to the outcome-model prediction term,
  NOT to the density-ratio weight. The density ratio uses the observed
  treatment value.
- `set.seed()` before `contrast()` MUST make the MC integration
  reproducible.

## DGP for truth-based tests

### Point (binary treatment, Gaussian outcome)

Reuse the Phase 12 stochastic DGP:
```
L ~ N(0, 1)
A | L ~ Bernoulli(expit(0.3 + 0.5 * L))     (observational)
g*(A | L): A* ~ Bernoulli(expit(0.5 + 0.3 * L))  (intervention)
density: function(a, data) dbinom(a, 1, plogis(0.5 + 0.3 * data$L))
Y = 2 + 3 * A + 1.5 * L + N(0, 1)
E[Y^{g*}] = 2 + 3 * E_L[expit(0.5 + 0.3 * L)] ≈ 4.136
```

### Point (continuous treatment, Gaussian outcome)

```
L ~ N(0, 1)
A | L ~ N(L, 1)                               (observational)
g*(A | L): A* ~ N(L + 0.5, 1)                 (shifted intervention)
density: function(a, data) dnorm(a, mean = data$L + 0.5, sd = 1)
Y = 2 + A + 1.5 * L + N(0, 1)
E[Y^{g*}] = 2 + E[L + 0.5] + 1.5 * E[L] = 2.5
```

### DR test DGP

Same as Phase 16 chunks 16d/16e — linear-Gaussian with binary
treatment. Deliberately misspecify one nuisance at a time:
- Wrong outcome: `Y ~ A` (drops L). AIPW recovers truth via correct propensity.
- Wrong propensity: `A ~ 1` (drops L). AIPW recovers truth via correct outcome.

## Key files for implementation

| File | Function(s) | Change |
|---|---|---|
| `R/interventions.R` | `stochastic()` | Add `density` param, validation |
| `R/ipw_weights_mv.R` | `check_intervention_family_compat()` | Gate stochastic on density presence |
| `R/ipw_weights.R` | `compute_density_ratio_weights()`, `make_weight_fn()` | New stochastic branches |
| `R/aipw.R` | `compute_aipw_contrast_point()` | MC-averaged preds for stochastic |
| `R/variance_if_aipw.R` | `variance_if_aipw()` | MC-averaged Channel 2a gradient |
| `R/variance_bootstrap.R` | `refit_aipw()` | No change needed (MC handled by contrast) |
| `R/utils.R` | — | Add `has_stochastic_density()` helper |

## References

- Díaz I, van der Laan MJ (2012). Population intervention causal
  effects based on stochastic interventions. *Biometrics* 68:541–549.
  **(foundational: A-IPTW = AIPW for stochastic interventions)**
- Díaz I, van der Laan MJ (2013). Assessing the causal effect of
  policies: an example using stochastic interventions. *International
  Journal of Biostatistics* 9:161–174.
  **(longitudinal A-IPTW extension)**
- Haneuse S, Rotnitzky A (2013). Estimation of the effect of
  interventions that modify the received treatment. *Statistics in
  Medicine* 32:5260–5277. **(introduced MTP terminology; DR estimators)**
- Young JG, Hernán MA, Robins JM (2014). Identification, estimation
  and approximation of risk under interventions that depend on the
  natural value of treatment. *Epidemiologic Methods* 3(1):1–19.
  **(MTP ↔ stochastic equivalence)**
- Kennedy EH (2019). Nonparametric causal effects based on incremental
  propensity score interventions. *JASA* 114:645–656.
  **(IPSI; avoids positivity via propensity-score multiplication)**
- Díaz I, Williams N, Hoffman KL, Schenck EJ (2023). Nonparametric
  causal effects based on longitudinal modified treatment policies.
  *JASA* 118:846–857. **(density ratio classification trick; SDR)**
- Wen L, Marcus JL, Young JG (2023). Intervention treatment
  distributions that depend on the observed treatment process and model
  double robustness in causal survival analysis. *Statistical Methods
  in Medical Research* 32:509–523.
  **(caveat: DR fails when g* depends on observed A)**
- Shook-Sa BE, Zivich PN, Lee C, Xue K, Ross RK, Edwards JK,
  Stringer JSA, Cole SR (2025). Double robust variance estimation
  with parametric working models. *Biometrics* 81:ujaf054.
  **(sandwich SE is itself doubly robust)**
- Stefanski LA, Boos DD (2002). The calculus of M-estimation.
  *American Statistician* 56:29–38.
  **(general stacked M-estimation theory)**
- Hernán MA, Robins JM (2025). *Causal Inference: What If*. Chapman &
  Hall/CRC. Chapter 21. **(stochastic strategies textbook reference)**
