# Phase 12 — Stochastic Interventions under G-computation

> **Status: PENDING (design doc)**
>
> **Depends on:** Phases 2 (point gcomp) and 5 (ICE) — both done.
> Stochastic interventions are gcomp-only and do not touch the IPW engine.

## Motivation

causatr currently supports **deterministic** interventions (`static`,
`shift`, `scale_by`, `threshold`, `dynamic`) and one **stochastic**
intervention (`ipsi`) limited to IPW on binary treatments.

General stochastic interventions — where the counterfactual treatment is
drawn from a user-specified distribution `g(a | L)` rather than set to a
fixed value — are an important class of causal estimands studied by
Díaz & van der Laan (2012) and used throughout the `lmtp` ecosystem.
They naturally model real-world policies where compliance is imperfect,
treatment assignment is randomised, or the intervention is a distribution
shift rather than a point mass.

**G-computation is the natural home for stochastic interventions** in
causatr. The parametric g-formula computes `E[Y^g] = E_L[ E_g[ E[Y | A, L] ] ]`
by averaging model predictions over both the covariate distribution
(outer expectation) and the stochastic intervention distribution (inner
expectation). The inner expectation is evaluated via Monte Carlo
integration — sample `M` treatment draws from `g(· | L_i)` per
individual, predict under each, and average.

IPW and matching do **not** support general stochastic interventions:
the density-ratio engine requires either a deterministic map (for the
Jacobian) or a closed-form weight (IPSI only). Matching has no
mechanism for stochastic treatment assignment. These rejection paths
are already in place.

## References

- Díaz I, van der Laan MJ (2012). Population intervention causal
  effects based on stochastic interventions. *Biometrics* 68:541–549.
- Young JG, Hernán MA, Robins JM (2014). Identification, estimation
  and approximation of risk under interventions that depend on the
  natural value of treatment using observational data.
  *Epidemiologic Methods* 3(1):1–19.
- Hernán MA, Robins JM (2025). *Causal Inference: What If*. Chapman &
  Hall/CRC. Chapter 21 (stochastic strategies).

## Scope

1. New `stochastic()` intervention constructor.
2. Monte Carlo g-formula integration for **point** gcomp.
3. Monte Carlo g-formula integration for **longitudinal** ICE gcomp.
4. Sandwich variance via averaged IFs (finite-MC consistent).
5. Bootstrap variance (natural — MC sampling is inside each replicate).
6. Rejection paths for IPW, matching.
7. Truth-based simulation tests validated against `lmtp::lmtp_tmle()`.

## Non-scope

- Stochastic interventions under IPW (beyond `ipsi()`) — would require
  the user to supply the intervention density, and the variance engine
  would need a fundamentally different cross-derivative. Out of scope.
- Stochastic interventions under matching — architecturally incompatible.
- Optimal treatment regime estimation (DTR) — different problem class.
- TMLE/SDR for stochastic interventions — `lmtp` covers this.

---

## Design

### 1. `stochastic()` constructor

```r
#' Stochastic intervention (user-supplied sampling function)
#'
#' Creates a stochastic intervention where the counterfactual treatment
#' for each individual is drawn from a user-supplied distribution that
#' may depend on the individual's covariates. Each call to `sampler`
#' should return an independent draw.
#'
#' Only supported under `estimator = "gcomp"` (point and longitudinal).
#' The g-formula evaluates E[Y^g] via Monte Carlo integration: for each
#' of `n_mc` draws, the sampler assigns counterfactual treatments, the
#' outcome model predicts, and the predictions are averaged across draws.
#'
#' @param sampler A function with signature `function(data, treatment)`
#'   that returns a vector of treatment values of length `nrow(data)`.
#'   `data` is the full counterfactual data.table (same as `dynamic()`);
#'   `treatment` is the observed treatment vector. Each call must return
#'   an independent random draw from the stochastic policy.
#' @param n_mc Positive integer. Number of Monte Carlo draws for the
#'   g-formula integration. Default 100. Larger values reduce MC noise
#'   at the cost of computation time. For sandwich variance, `n_mc`
#'   should be large enough that MC noise is negligible relative to
#'   the estimation error (100–500 is typical).
#'
#' @return A `causatr_intervention` object.
#'
#' @references
#' Díaz I, van der Laan MJ (2012). Population intervention causal
#' effects based on stochastic interventions. *Biometrics* 68:541–549.
stochastic <- function(sampler, n_mc = 100L) { ... }
```

**Stored fields:** `type = "stochastic"`, `sampler`, `n_mc`.

**Convenience wrappers** (future, not Phase 10):
- `bernoulli_iv(prob_fn)` — binary treatment drawn Bernoulli with
  `P(A* = 1 | L) = prob_fn(data)`. Sugar for
  `stochastic(\(d, t) rbinom(nrow(d), 1, prob_fn(d)))`.
- `mixture(interventions, weights)` — two-point (or k-point) mixture.
  Sugar for `stochastic(\(d, t) { ... sample index, apply ... })`.

These are deferred to avoid scope creep; the raw `stochastic()` is the
primitive, and the convenience wrappers are trivial once it works.

### 2. `apply_single_intervention()` — stochastic branch

Add a `stochastic` case to the `switch()` in `apply_single_intervention()`:

```r
stochastic = {
  new_trt <- iv$sampler(data, data[[trt_col]])
  # same type-checking as `dynamic` branch
  data[, (trt_col) := new_trt]
}
```

Each call produces a **different** counterfactual dataset because the
sampler draws from a distribution. This is the key difference from
`dynamic()`, which is deterministic.

### 3. Point gcomp — MC integration in `compute_contrast()`

The current flow (deterministic):
```
data_a = apply_intervention(data, trt, iv)        # one dataset
preds  = predict(model, newdata = data_a)          # one pred vector
mu_hat = mean(preds[target])                       # one scalar
```

The stochastic flow:
```
preds_mat = matrix(0, nrow = n, ncol = n_mc)
for (m in 1:n_mc) {
  data_a_m  = apply_intervention(data, trt, iv)    # each call samples
  preds_mat[, m] = predict(model, newdata = data_a_m)
}
preds_avg = rowMeans(preds_mat)                    # per-individual MC avg
mu_hat    = mean(preds_avg[target])                # population mean
```

**Implementation strategy:** wrap the existing `lapply(interventions, ...)`
block in `compute_contrast()` with a helper that detects stochastic
interventions and runs the MC loop, returning the averaged predictions.
Deterministic interventions go through the existing single-shot path
(no performance penalty).

Helper signature:
```r
#' Predict under an intervention (deterministic or stochastic MC)
#'
#' For deterministic interventions, returns a single prediction vector.
#' For stochastic interventions, returns the row-wise MC average.
#' Also returns the per-draw prediction matrix (n × M) for sandwich
#' variance computation.
predict_under_intervention <- function(model, data, treatment, iv) { ... }
```

### 4. ICE — MC integration in `ice_iterate()`

At each backward step k, the current flow predicts under a single
deterministic intervention. For stochastic interventions:

```
At step k:
  for m in 1:n_mc:
    sample A_k from g(· | L_k) via apply_single_intervention()
    predict pseudo_outcome_k_m from model_k
  pseudo_outcome_k = rowMeans(pseudo_outcome_k_m)  # MC average
  pass pseudo_outcome_k to step k-1
```

The MC averaging happens **within** each backward step, so the
pseudo-outcome passed to the next step is already integrated over the
stochastic policy at time k. This is consistent with the ICE algorithm's
backward substitution: `E_g_k[ E[Y | A_k, H_k] ]` at each step.

### 5. Sandwich variance — averaged IFs

The estimator under a stochastic intervention `g` is:
```
μ̂_g = (1/n) Σ_i (1/M) Σ_m Ê[Y | A_{i,m}, L_i]
```
where `A_{i,m} ~ g(· | L_i)` are iid draws.

The IF of this estimator is the **MC average** of the per-draw IFs:
```
IF_i = (1/M) Σ_m IF_m(i)
```
where `IF_m(i)` is the standard gcomp IF (Channel 1 + Channel 2)
evaluated at the m-th counterfactual dataset.

**Channel 1 (direct term):**
```
Ch1_m(i) = Ê[Y | A_{i,m}, L_i] - μ̂_g
Ch1(i) = (1/M) Σ_m Ch1_m(i)
       = (1/M) Σ_m Ê[Y | A_{i,m}, L_i] - μ̂_g
       = preds_avg[i] - μ̂_g
```
This is just the averaged prediction minus the mean — identical to
what we'd compute from `preds_avg`.

**Channel 2 (model correction):**
```
Ch2_m(i) = Σ_j ∂Ê[Y | A_{j,m}, L_j]/∂β  ·  bread^{-1}  ·  score_i
Ch2(i)   = (1/M) Σ_m Ch2_m(i)
         = (Σ_j (1/M) Σ_m ∂Ê[Y | A_{j,m}, L_j]/∂β)  ·  bread^{-1}  ·  score_i
```
The MC average of the gradient is the gradient of the MC-averaged
predictions, so Channel 2 reduces to:
```
Ch2(i) = (Σ_j d_avg_j) · bread^{-1} · score_i
```
where `d_avg_j = (1/M) Σ_m ∂predict(model, A_{j,m}, L_j)/∂β`.

**Implementation:** `apply_model_correction()` already takes a gradient
vector `d_fit`. For stochastic interventions, we pass the MC-averaged
gradient instead of the single-draw gradient. The variance engine
needs no structural change — only the gradient input changes.

The **per-draw gradient matrix** `d_mat[j, p, m]` can be computed as
part of the MC loop in `predict_under_intervention()`:
```r
d_avg <- rowMeans_over_m(d_mat)  # n × p matrix, MC-averaged
```

For ICE, the same principle applies: at each backward step, the
sensitivity matrix `A_{k,k-1}` uses the MC-averaged gradient from
step k's stochastic predictions.

### 6. Bootstrap variance

No special treatment needed. Each bootstrap replicate independently:
1. Resamples individuals.
2. Refits the outcome model.
3. Runs the MC g-formula (sampling from the stochastic policy).
4. Returns the marginal mean.

The bootstrap variance naturally captures both model uncertainty and
MC noise. If `n_mc` is large, MC noise is negligible and the bootstrap
variance converges to the model variance alone.

### 7. Rejection paths

Add `"stochastic"` to the rejection set in:
- `check_interventions_compat()` for matching (already rejects all
  non-static; stochastic falls into this bucket).
- `check_intervention_family_compat()` in `R/ipw_weights.R` for IPW.
- `apply_single_intervention()` under the `ipsi` case already rejects
  non-IPW; add a parallel `stochastic` case that rejects non-gcomp.

Error message template:
```
`stochastic()` interventions are only supported under `estimator = 'gcomp'`.
ℹ Stochastic interventions require Monte Carlo integration over the
  counterfactual treatment distribution, which is only available via the
  g-formula (outcome-model prediction).
ℹ Use `causat(..., estimator = 'gcomp')`, or rewrite the intervention
  as a deterministic `dynamic()` rule if the policy is actually
  deterministic.
```

### 8. Estimand compatibility

Stochastic interventions are compatible with:
- **ATE**: yes (average over full population).
- **ATT/ATC under gcomp**: yes (average over treated/controls).
  The stochastic policy draws are independent of the estimand — we just
  change which rows we average over.
- **`subset`**: yes.
- **`by`**: yes.

No additional gating needed beyond the IPW/matching rejection.

---

## Chunks

### Chunk 1: Constructor + rejection paths

1. Add `stochastic()` constructor to `R/interventions.R`.
2. Add `stochastic` case to `apply_single_intervention()`.
3. Add rejection in IPW path (`check_intervention_family_compat()`).
4. Add tests for constructor validation and rejection paths.
5. `devtools::document()` + `devtools::test()`.

### Chunk 2: Point gcomp MC integration

1. Add `predict_under_intervention()` helper (MC loop for stochastic,
   single-shot for deterministic).
2. Modify the point gcomp block in `compute_contrast()` to use the
   helper.
3. Truth-based simulation test: binary treatment stochastic policy
   `P(A* = 1 | L) = expit(γ₀ + γ₁ L)` on a linear-Gaussian DGP
   where the analytical truth is `E[Y^g] = β₀ + β_A · E_L[p_g(L)] + β_L · E[L]`.
4. Cross-validate against `lmtp::lmtp_tmle()` on the same DGP.
5. `devtools::test()`.

### Chunk 3: Point gcomp sandwich for stochastic

1. Modify the sandwich path to compute MC-averaged gradients.
2. Pass averaged gradient to `apply_model_correction()`.
3. Truth-based test: compare sandwich SE to the Monte Carlo SE from
   repeated simulations (should agree within MC noise).
4. Cross-validate SE against bootstrap SE (should agree closely for
   large `n_mc`).
5. `devtools::test()`.

### Chunk 4: ICE MC integration

1. Modify `ice_iterate()` to handle stochastic interventions via
   per-step MC averaging.
2. Truth-based simulation test: 2-period longitudinal DGP with
   stochastic policy, validated against `lmtp::lmtp_tmle()`.
3. `devtools::test()`.

### Chunk 5: ICE sandwich for stochastic

1. Modify `variance_if_ice_one()` to use MC-averaged gradients in
   the forward sensitivity recursion.
2. Cross-validate against bootstrap SE.
3. `devtools::test()`.

### Chunk 6: Documentation + matrix update

1. Add `stochastic()` to the `@seealso` lists in all intervention
   constructors.
2. Add stochastic examples to the `contrast()` documentation.
3. Update `FEATURE_COVERAGE_MATRIX.md`.
4. Update `CLAUDE.md`.
5. `devtools::document()` + `devtools::check()`.

---

## Invariants

- `stochastic()` MUST be rejected by IPW and matching with informative
  error messages pointing to `estimator = "gcomp"`.
- With `n_mc = 1`, `stochastic()` degenerates to a single random draw
  (equivalent to `dynamic()` with a stochastic rule). This is valid
  but has high MC variance — warn if `n_mc < 10`.
- The MC integration MUST be reproducible when `set.seed()` is called
  before `contrast()`. The sampler function is called `n_mc` times
  per intervention per individual — the RNG state advances
  deterministically.
- For sandwich variance: MC-averaged IFs MUST produce the same vcov
  (up to MC noise) as the full IF computed from the averaged
  predictions. This is the "linearity of the IF average" property.
- Bootstrap replicates MUST draw fresh MC samples independently
  (not reuse the same draws across replicates).

## DGP for truth-based tests

### Point treatment (Chunk 2–3)

```
L ~ N(0, 1)
g(A | L): A* ~ Bernoulli(expit(0.5 + 0.3 * L))
Y = 2 + 3 * A + 1.5 * L + ε,  ε ~ N(0, 1)

E[Y^g] = 2 + 3 * E_L[expit(0.5 + 0.3 * L)] + 1.5 * E[L]
       = 2 + 3 * E_L[expit(0.5 + 0.3 * L)]   (since E[L] = 0)
```
The inner expectation `E_L[expit(0.5 + 0.3 * L)]` can be computed
numerically via `integrate()` against the N(0,1) density, giving an
exact analytical truth for the test.

### Longitudinal (Chunk 4–5)

```
L_0 ~ N(0, 1)
g_0: A_0 ~ Bernoulli(expit(0.2 + 0.3 * L_0))
L_1 = 0.5 * A_0 + 0.5 * L_0 + ε_L
g_1: A_1 ~ Bernoulli(expit(0.2 + 0.3 * L_1))
Y = 1 + A_0 + A_1 + 0.5 * L_1 + ε_Y

E[Y^g] computed via nested numerical integration or
validated against lmtp::lmtp_tmle() with n = 50000.
```
