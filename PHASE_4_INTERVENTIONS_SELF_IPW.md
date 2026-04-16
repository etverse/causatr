# Phase 4 — Extended Interventions + Self-Contained IPW Engine

> **Status: PENDING**
> Book chapters: 12 (density ratio weights), 19 (dynamic strategies), Díaz et al. 2023 (MTPs), Kennedy 2019 (IPSI)

## Scope

A single self-contained IPW engine that handles `static`, `shift`, `scale_by`, `dynamic`, and `ipsi` via density-ratio weights and an explicit MSM. `threshold()` is explicitly **not** supported under IPW — the pushforward of a continuous density under a boundary clamp is a mixed measure (point masses at `lo`/`hi`), so the density ratio w.r.t. the fitted `f(a|l)` is not well-defined. Users who want the clamped-treatment counterfactual should switch to `estimator = "gcomp"`, where it is handled cleanly via predict-then-average on the outcome model. WeightIt is removed from the runtime path and kept only as a test oracle for the static binary case. Adds categorical treatment support and the interventions vignette.

## Why a single engine instead of "WeightIt for static, ours for dynamic"

The original Phase 4 plan branched at fit time: WeightIt for static interventions, a new density-ratio engine for everything else. That split is a maintenance trap:

- Two propensity-score code paths to keep in sync (one M-estimation framework via `glm_weightit`, one closure-based with `numDeriv`).
- Two variance branches (Branch A = WeightIt-based with `sandwich::estfun`, Branch B = self-contained) that both have to give the same answer on static binary fits — which means we'd need T1-style equivalence tests forever, not just once during migration.
- Static binary IPW is the simplest case — a logistic propensity, density-ratio weights that collapse to `1/p` and `1/(1-p)`, an MSM that collapses to `Y ~ A`. If the unified engine can't handle it, it can't handle anything.

So Phase 4 **unifies** the engine: one `R/treatment_model.R`, one weight builder, one MSM step, one variance branch. WeightIt becomes a dev dependency — **Suggests, not Imports** — invoked only by tests as an independent reference.

## The MSM stays explicit

Even though static binary IPW could collapse to weighted means (`Σ w_i Y_i / Σ w_i` per arm), the engine keeps an explicit marginal structural model fit on the weighted pseudo-population:

```
g(E[Y^{d}]) = β₀ + β₁ s(d)      # MSM, fit by weighted GLM
```

Reasons:

1. **Non-saturated MSMs are the whole point of MTPs.** A family of shifts `shift(δ)` for δ ∈ {-10, ..., 10} only becomes a dose-response curve when you fit `Y ~ ns(δ, 3)` across the pooled pseudo-populations. Collapsing to weighted means throws that away.
2. **Effect modification (Phase 8).** The MSM is the natural slot to carry `Y ~ A + V + A:V`. Hard-coding weighted means would re-create the exact `check_confounders_no_treatment()` abort we're trying to remove in Phase 8.
3. **Architectural symmetry.** `compute_contrast()` already predicts-then-averages from a fitted model for every estimator (gcomp, IPW, matching) — keeping IPW on that pattern means one contrast machinery, not two.

For static binary the MSM degenerates to the saturated `Y ~ A` we have today, so nothing breaks for current users.

## `model_fn` covers the propensity backend

`causat()` already exposes `model_fn` (default `stats::glm`) for the outcome model in gcomp / ICE / survival. Phase 4 reuses that convention: the treatment density model is fit with the same `model_fn` contract `function(formula, data, family, weights, ...)`, so users get GAMs, splines, fractional polynomials, `MASS::glm.nb`, etc. for free without us reimplementing a propensity backend zoo.

Two API choices, to decide during implementation:

- **Option A (single fitter):** `model_fn` fits both the outcome MSM and the propensity model in the same `causat()` call. Simpler API; rules out "GLM outcome, GAM propensity."
- **Option B (split fitter):** add a sibling `propensity_model_fn` (default = `model_fn`) so the two stages can use different fitters. Slightly more surface, but matches how triangulation users actually think about it (different model classes for outcome vs treatment).

Recommendation: **Option B**, default `propensity_model_fn = model_fn`, so the simple case stays one argument and the flexible case stays one extra argument.

This deliberately drops WeightIt's eight non-GLM propensity backends (`gbm`, `super`, `bart`, `cbps`, `energy`, `optweight`, `npcbps`, `cfd`). Users who need those can still fit weights externally and pass them via `causat()`'s `weights =` slot (validated by `check_weights()`); they'll lose the propensity-uncertainty correction in the sandwich, which is exactly the `non-Mparts WeightIt method → rlang::warn()` situation we already document.

## Plan

### 1. Treatment density model — `R/treatment_model.R` (new)

Fit a treatment model and expose density evaluation under any treatment value:

```r
fit_treatment_model <- function(data, treatment, confounders,
                                family, model_fn = stats::glm, ...)
# Returns: a `causatr_treatment_model` carrying the fitted model,
# the family, the residual scale (for continuous), the original
# design info, and a closure d_fn(a, newdata) -> density value.

evaluate_density <- function(treatment_model, treatment_values, newdata)
# - binary: Bernoulli, p = predict(type = "response")
# - categorical: multinomial, p_k from predicted class probabilities
# - continuous: Normal(mean = predict, sd = residual sigma); other
#   continuous families documented as "user can swap via model_fn"
```

For continuous treatments the Gaussian assumption is the standard book treatment (Ch. 12); document it loudly and provide an escape hatch via `model_fn` (e.g. `mgcv::gam` with a non-Gaussian family) for users who need it.

### 2. Density-ratio weights — `R/ipw_weights.R` (new)

```r
compute_density_ratio_weights <- function(treatment_model,
                                          data,
                                          data_intervened,
                                          treatment)
# Returns w_i = f(a_intervened_i | L_i) / f(a_observed_i | L_i)
```

For static binary this collapses to `1/p` (treated arm) and `1/(1-p)` (control arm) — verified against `WeightIt::weightit(method = "glm")` in tests.

### 3. MSM fit on the weighted pseudo-population

`fit_ipw()` becomes a thin pipeline:

1. Fit treatment density model via `propensity_model_fn`.
2. For each intervention requested in `contrast()`, build density-ratio weights.
3. Pool weighted observations across interventions (or fit one MSM per intervention for the saturated case) and fit the MSM via `model_fn` on `(Y, design from the MSM formula, weights = w_i)`.
4. Hand the fitted MSM to `compute_contrast()`'s standard predict-then-average machinery — same code path as gcomp.

The MSM formula defaults to `Y ~ A` (saturated), with a Phase 8 hook for `Y ~ A + V + A:V` once `parse_effect_mod()` lands.

### 4. Supported intervention × method combinations after Phase 4

| Intervention | G-comp | IPW | Matching |
|---|---|---|---|
| `static()`    | ✓ | ✓ (binary / categorical) | ✓ (binary) |
| `shift()`     | ✓ | ✓ (continuous / count — integer shifts only on count) | — |
| `scale_by()`  | ✓ | ✓ (continuous / count — integer-preserving scales only on count) | — |
| `threshold()` | ✓ | ⛔ (use gcomp — pushforward has point masses) | — |
| `dynamic()`   | ✓ | ✓ (binary / categorical only) | — |
| `ipsi()`      | — | ✓ (binary) | — |

There is no longer a separate "IPW (WeightIt)" vs "IPW (self-contained)" column — there's just IPW.

The two ⛔ / "family restriction" notes reflect a genuine architectural limit of the density-ratio framework, not a scope decision:

- **`threshold()` under IPW** is rejected because the pushforward of a continuous `f(a|l)` under a boundary clamp is a mixed measure — continuous density on `(lo, hi)` plus point masses at the boundaries — so the density ratio w.r.t. the fitted Gaussian is not well-defined. `gcomp` handles it cleanly via `predict(outcome_model, newdata = clamped)` and is the correct tool for this intervention.
- **`dynamic()` under IPW** is only defined on binary / categorical treatments, where the HT indicator weight `I(A_obs = rule_i) / f(A_obs | L)` is well-posed. Deterministic rules on continuous treatments are a Dirac per individual, which has the same pushforward-degeneracy problem as `threshold()` — users should either rewrite the rule as a smooth `shift()` / `scale_by()`, or switch to `gcomp`.
- **`static(v)` under IPW** is rejected on continuous treatments (nobody is observed exactly at `v`, so the HT indicator is zero almost surely).
- **`shift()` / `scale_by()`** require a numeric treatment with a well-defined density. On continuous (Gaussian) treatments the density ratio uses the Gaussian pdf. On count treatments (Poisson / negative binomial, opt-in via `propensity_family`), the density ratio uses `dpois()` / `dnbinom()` — but only for shifts/scales that preserve the integer support (non-integer shifts on count treatments are rejected because `dpois()` at non-integer values returns 0). These interventions have no meaning on binary or categorical treatments.
- **`ipsi()`** is binary-only — Kennedy (2019)'s closed form is Bernoulli-specific.

All of these are enforced at contrast time by `check_intervention_family_compat()` in `R/ipw_weights.R`, which aborts with an actionable pointer.

### 5. Categorical treatment support

- G-comp: already works (factor treatment column in the outcome GLM).
- IPW: multinomial treatment density model via `nnet::multinom` or `VGAM::vglm` plugged in through `propensity_model_fn`; categorical density evaluation in `evaluate_density()`.
- Matching: stays binary-only (MatchIt limitation; we already error out cleanly in `fit_matching()`).
- Update `check_estimand_trt_compat()` to allow categorical + ATE.

### 5a. Count treatment support (Poisson / negative binomial)

Integer-valued treatments (dose levels, event counts, visit counts) currently fall through `detect_treatment_family()` to `"gaussian"` and get a Normal density model. This works when the integer support is dense enough that the Gaussian approximation is harmless (age in years, cigarettes/day). It breaks down on sparse count data (0–5 doses) where the Gaussian pdf evaluates between integers that carry zero true probability mass, producing badly calibrated density ratios.

**Design decisions (settled):**

1. **Explicit opt-in, not auto-detection.** A new `propensity_family` argument on `causat()` lets the user declare `propensity_family = "poisson"` or `propensity_family = "negbin"`. Auto-detecting count data from the values is fragile (age in years, income in thousands are all non-negative integers but should not get a Poisson density). The existing auto-detect for bernoulli / gaussian / categorical stays unchanged — those are unambiguous from the data. `propensity_family = NULL` (default) preserves the current behaviour.

2. **Both Poisson and negative binomial.** Poisson is a one-parameter model (`lambda = exp(X %*% alpha)`, no dispersion). Negative binomial adds a `theta` (size) parameter estimated by `MASS::glm.nb()`. Theta is treated as fixed in the variance engine — same convention as the Gaussian `sigma` — because perturbing `theta` inside `numDeriv::jacobian` would require a joint `(alpha, theta)` M-estimation setup. `MASS::glm.nb` is a recommended package (ships with R, like `nnet`).

3. **Non-integer shifts/scales are rejected.** `shift(delta)` on a count treatment is only valid when `delta` is an integer, because `dpois(a - delta, lambda)` returns 0 for non-integer arguments. Similarly `scale_by(factor)` is only valid when `a_obs * factor` is integer for every observed value. `check_intervention_family_compat()` enforces this with a clear error pointing at `estimator = "gcomp"` as the fallback for fractional shifts on count data (g-comp handles it via predict-then-average on the outcome model, no density ratio needed).

**Implementation plan:**

Touch points (5 functions modified, 1 new):

| File | Function | Change |
|---|---|---|
| `R/treatment_model.R` | `fit_treatment_model()` | Add `propensity_family` parameter; when `"poisson"`, call `fit_count_density(family = poisson())`; when `"negbin"`, call `fit_count_density(family = NULL)` with `model_fn = MASS::glm.nb` |
| `R/treatment_model.R` | `fit_count_density()` | **New.** Thin wrapper like `fit_bernoulli_density()`. For Poisson: `model_fn(formula, data, family = poisson(), weights, ...)`. For NB: `MASS::glm.nb(formula, data, weights, ...)` (no `family` arg, same pattern as `nnet::multinom`). Extracts `theta` from the fitted model for NB (`model$theta`). |
| `R/treatment_model.R` | `evaluate_density()` | Add `"poisson"` branch: `lambda = predict(model, type = "response"); dpois(treatment_values, lambda)`. Add `"negbin"` branch: `lambda = predict(model, type = "response"); dnbinom(treatment_values, mu = lambda, size = theta)`. |
| `R/ipw_weights.R` | `make_weight_fn()` | Add count closure parallel to the Gaussian one. Key difference: `lambda = exp(X_prop %*% alpha)` (log link) instead of `mu = X_prop %*% alpha` (identity link). Density: `dpois(a, lambda)` or `dnbinom(a, mu = lambda, size = theta)`. Same pushforward structure: `f(d^{-1}(A_obs) | L) / f(A_obs | L)` with `|Jac| = 1` for integer shift. |
| `R/ipw_weights.R` | `check_intervention_family_compat()` | Allow `shift()` / `scale_by()` on `"poisson"` / `"negbin"`, with an additional guard: reject non-integer delta (for shift) and reject scale factors that produce non-integer `a_obs * factor` (for scale_by). |
| `R/ipw.R` | `fit_ipw()` | Thread `propensity_family` through to `fit_treatment_model()`. When `propensity_family = "negbin"` and `propensity_model_fn` is NULL, auto-select `MASS::glm.nb` (same pattern as categorical auto-selecting `nnet::multinom`). |
| `R/causat.R` | `causat()` | Add `propensity_family = NULL` argument, plumbed through to `fit_ipw()`. |

**Variance engine:** no changes needed. The Poisson and NB GLMs are standard GLMs with `$family`, `$linear.predictors`, working weights, etc. — `prepare_model_if()` and `bread_inv()` work out of the box. (`MASS::glm.nb` returns a `negbin` object inheriting from `glm`.)

**Tests (chunk 3j):**

- Truth-based DGP: count treatment `A ~ Poisson(exp(0.5*L))`, continuous outcome `Y = 2 + 1.5*A + L + eps`. True `E[Y(shift(1))] - E[Y(shift(0))]` derivable analytically from the Poisson MTP.
- NB parity: same DGP fit with `propensity_family = "negbin"` should agree with Poisson within tolerance (NB nests Poisson).
- Rejection tests: `shift(0.5)` on count treatment aborts; `scale_by(1.5)` on integer data with odd values aborts.

**References:**

- Díaz I, van der Laan MJ (2012). Population intervention causal effects based on stochastic interventions. *Biometrics* 68:541–549. — density-ratio framework for parametric treatment models, applicable to any exponential-family density.
- Haneuse S, Rotnitzky A (2013). Estimation of the effect of interventions that modify the received treatment. *Statistics in Medicine* 32:5260–5277. — modified treatment policies with explicit density-ratio weights under parametric treatment models.
- Cameron AC, Trivedi PK (2013). *Regression Analysis of Count Data.* Cambridge University Press. — Poisson and NB regression as parametric density models for count treatments.

### 6. IPSI (incremental propensity score interventions)

The `ipsi(delta)` constructor already exists. Implementation is a special case of the unified engine:

- Fit the propensity score (binary treatment model) once.
- The "intervened density" is `p_new = (δ·p) / (1 - p + δ·p)` — i.e. the same propensity model, evaluated under shifted odds.
- Density-ratio weights then come from the same `compute_density_ratio_weights()` path.

Because IPSI never realizes a counterfactual treatment value, its only natural estimator is IPW (no g-comp column).

### 7. Variance for the unified IPW engine

The variance engine (`R/variance_if.R`) was already refactored in Phase 3 to plug self-contained IPW in without touching `variance_if()` itself. The dispatcher `prepare_propensity_if()` routes to `prepare_propensity_if_self_contained()` whenever `fit$model` is **not** a `glm_weightit`. Because Phase 4 removes the `glm_weightit` path entirely, the dispatcher's Branch A falls dormant on the runtime side and survives only as part of the test oracle (see "WeightIt as test oracle" below).

**What `fit_ipw()` must populate in `fit$details`:**

1. `fit$details$propensity_model` — the treatment density model from `R/treatment_model.R` (a `glm` / `gam` / `multinom` / etc., produced via `propensity_model_fn`).
2. `fit$details$alpha_hat` — the propensity parameter vector \(\hat\alpha\).
3. `fit$details$weight_fn` — an R closure \(\alpha \mapsto w_i(\alpha)\) that reconstructs the length-`n_fit` weight vector under a candidate \(\alpha\). The closure captures the intervention, the observed treatment, and the density functions so `numDeriv::jacobian()` can perturb \(\alpha\) alone.

Without these three slots, the self-contained branch cannot be implemented.

**Implementation steps** (the engine is already designed around this shape — see the `prepare_model_if()` / `apply_model_correction()` split in `R/variance_if.R`):

1. **MSM correction term.** Prepare once, apply per intervention:
   ```r
   msm_prep <- prepare_model_if(fit$model, fit_idx, n_total)
   msm      <- apply_model_correction(msm_prep, J)
   ```
   `msm` carries both `msm$correction` (the MSM Channel-2 contribution) *and* `msm$h = A_{ββ}^{-1} J`. The `$h` slot is why `apply_model_correction()` returns `h` alongside `$correction`.

2. **Cross-derivative \(A_{\beta\alpha}\) via `numDeriv::jacobian`.** Build
   ```r
   psi_beta_bar <- function(alpha) {
     w      <- fit$details$weight_fn(alpha)
     X      <- model.matrix(fit$model)
     fam    <- fit$model$family
     eta    <- as.numeric(X %*% beta_hat)
     mu     <- fam$linkinv(eta)
     mu_eta <- fam$mu.eta(eta)
     r      <- Y_fit - mu
     as.numeric(crossprod(X, w * mu_eta * r / fam$variance(mu))) / n_fit
   }
   A_beta_alpha <- -numDeriv::jacobian(psi_beta_bar, x = fit$details$alpha_hat)
   ```
   The numerical jacobian is the right tool here (not an analytic formula) because `weight_fn` changes per intervention type — shift, MTP, IPSI, dynamic rule — and re-deriving the closed form for each is tedious and error-prone.

3. **Propensity correction via a derived gradient.** Form
   \[
   g^{\mathrm{prop}} = A_{\beta\alpha}^T h_{\mathrm{msm}}
   \]
   and feed it back into the same primitive on the **propensity** model:
   ```r
   g_prop    <- as.numeric(crossprod(A_beta_alpha, msm$h))
   prop_prep <- prepare_model_if(fit$details$propensity_model, fit_idx_prop, n_total)
   prop      <- apply_model_correction(prop_prep, g_prop)
   Ch2       <- msm$correction + prop$correction
   ```

4. **Return.** Wrap everything in a `prop_prep` list whose `apply_propensity_correction()` returns the combined MSM + propensity correction per individual.

### 8. Estimand support under the unified engine

ATE is well-defined for every intervention the engine handles. ATT / ATC are only well-defined when "the treated" / "the controls" is an unambiguous subpopulation — i.e. static interventions on a binary treatment. For `shift` / `scale_by` / `threshold` / `dynamic` / `ipsi` the notion of "among the treated" is either undefined or exotic, the literature doesn't use it, and silently falling back to ATE weights under an ATT request would be wrong.

**Rule:**

| Intervention | ATE | ATT | ATC |
|---|---|---|---|
| `static()` on binary 0/1 treatment | ✓ | ✓ | ✓ |
| `static()` on categorical / continuous treatment | ✓ | ⛔ | ⛔ |
| `shift()` / `scale_by()` / `threshold()` / `dynamic()` | ✓ | ⛔ | ⛔ |
| `ipsi()` | ✓ | ⛔ | ⛔ |

Implementation: add `check_estimand_intervention_compat(estimand, intervention)` alongside the existing `check_estimand_trt_compat()`. The new check runs at **contrast time**, not fit time, because the intervention is not known at fit time for IPW. Phase 3 already fixes the estimand at fit time for IPW via `check_estimand_compat()` — that stays. The new check just adds: when the fit-time estimand is ATT/ATC, every intervention in the `contrast()` call must be `static` (on a binary 0/1 treatment).

Error class: `causatr_bad_estimand_intervention`. Message points users to `estimator = "gcomp"` for non-static ATT/ATC, or to switch to ATE.

### 9. Oracle testing strategy: contrast-level comparisons, no runtime branching

**Decision: Branch A is removed from R/ entirely.** The Phase 3 variance engine carried a `prepare_propensity_if()` dispatcher that routed between Branch A (`sandwich::estfun(asympt = TRUE)` on `glm_weightit`) and Branch B (self-contained, first-principles). The dispatcher was scaffolding — a way for Phase 3 to ship with Phase 4's slot names already in place. Once Phase 4 flips the runtime, keeping Branch A alive in `R/` is dead-code-waiting-to-rot:

- Two code paths to keep in sync indefinitely even though only one runs.
- Two variance branches that must give the same answer on static binary forever, maintained via a brittle cross-branch test.
- A runtime dependency on `WeightIt::glm_weightit` for a path nobody takes in production.

**Chunk 3c deletes** `prepare_propensity_if()`, `prepare_propensity_if_weightit()`, `prepare_propensity_if_self_contained()`, and `apply_propensity_correction()` from `R/variance_if.R`. The new `variance_if_ipw()` is a straight loop over interventions calling `compute_ipw_if_self_contained_one()` (committed in Chunk 3b, 6e3d42b) and aggregating via `vcov_from_if()`. One code path, no branching. The Phase 3 IPW tests either adapt to the new architecture or get deleted (see §12 below for the tests list).

**Oracle strategy: contrast-level, not IF-level.** The WeightIt and lmtp oracle tests compare **final contrast estimates and their standard errors**, not per-individual influence functions. Two reasons:

1. **Architecture independence.** The causatr self-contained engine uses HT indicator weights with intercept-only MSMs (`w_i = I(A_i = a*) / f(A_i|L_i)`, `glm(Y ~ 1, weights = w)`) for static binary, while WeightIt's `glm_weightit` uses full-IPW weights with a saturated MSM (`w_i = 1/f(A_i|L_i)`, `glm(Y ~ A, weights = w)`). These are different architectures that estimate the same functional — their per-observation IFs are not elementwise comparable, but the final point estimate and SE of `mu_1 - mu_0` agree asymptotically (and, for the same data, numerically at ~1e-6). Comparing at the contrast level lets the test be agnostic to the internal decomposition.

2. **Robustness to internal refactoring.** IF-level tests couple the test suite to implementation details of the variance engine. If we ever change the IF decomposition (say, to hoist some more computation out of the per-intervention loop for perf), an IF-level oracle breaks even though the final answer is correct. Contrast-level oracles test what the user actually consumes: `$estimates$estimate`, `$estimates$se`, `$contrasts$estimate`, `$contrasts$se`.

**WeightIt — contrast-level oracle for binary × {ATE, ATT, ATC}.** `WeightIt::glm_weightit()` supports all three estimands on binary treatment and has been field-tested for a decade. Ideal reference for the static binary case:

- **T-oracle1:** `causat(estimator = "ipw", estimand = "ATE")` + `contrast(static(1) vs static(0))` point estimate and SE ≈ `WeightIt::glm_weightit(weightit(..., estimand = "ATE"))` + manual `Y ~ A` MSM contrast, tolerance ~1e-6.
- **T-oracle2:** Same pipeline, `estimand = "ATT"`.
- **T-oracle3:** Same pipeline, `estimand = "ATC"`.
- **T-oracle4:** `estimand = "ATE"` but with `propensity_model_fn = mgcv::gam` — proves the engine handles non-GLM propensity models correctly.

Each test is wrapped in `skip_if_not_installed("WeightIt")`. All oracle code (the WeightIt reference fit + contrast extraction helper) lives in `tests/testthat/helper-ipw-weightit-oracle.R`, **not in `R/`**. WeightIt is a test-only dependency in the final architecture.

**lmtp — contrast-level oracle for non-static interventions.** `lmtp::lmtp_sdr()` with parametric learners (`learners_trt = "SL.glm"`, `learners_outcome = "SL.glm"`, `folds = 1`) is the only publicly-available package that supports the same shift / MTP / IPSI estimands as our engine. (`lmtp_ipw()` was removed in lmtp >= 1.5.0.) For a given correctly specified parametric model the SDR point estimate is consistent for the same target parameter as causatr's density-ratio IPW. It is **not** a variance oracle because its SE uses the EIF rather than the M-estimation sandwich — different objects, only asymptotically comparable.

- **T-oracle5:** `causat(shift(1))` on a linear-Gaussian DGP (`simulate_continuous_continuous`, n=3000) point estimate ≈ `lmtp::lmtp_sdr(shift = function(d, t) d[[t]] + 1, learners_trt = "SL.glm", learners_outcome = "SL.glm", folds = 1)` point estimate, tolerance ~0.15. `lmtp_ipw()` was removed in lmtp >= 1.5.0; SDR with correctly specified parametric learners is consistent for the same target parameter. The wider tolerance (vs ~1e-6 for WeightIt oracles) reflects that SDR is a different estimator (doubly robust vs pure IPW), not a calibration issue.
- **T-oracle6:** `causat(ipsi(2))` on a linear-Gaussian DGP (`simulate_binary_continuous`, n=3000) point estimate ≈ a manual IPSI weight + Hájek weighted mean comparison, tolerance ~1e-6. The manual oracle fits the same logistic propensity model and computes the Kennedy (2019) closed-form weight `w_i = (δ·A_i + (1−A_i)) / (δ·p_i + (1−p_i))` from first principles, with no causatr code involved.
- **T-oracle7 (stretch):** categorical static ATE once the categorical branch lands — compare against `lmtp_sdr` with a multinomial propensity.

All `skip_if_not_installed("lmtp")` — lmtp is also test-only.

**Variance correctness stays anchored to first-principles tests**, not to WeightIt or lmtp:

- **T-A_β_α** (landed in Chunk 3a, commit 4090e51) — hand-derived cross-derivative vs `numDeriv::jacobian` for static(1), static(0), shift, IPSI, NULL. 12 assertions at ~1e-6.
- **T-end-to-end stacked sandwich by hand** (landed in Chunk 3b, commit 6e3d42b) — hand-assembled (alpha, beta) 3×3 block bread + meat, inverted manually, compared to `sum(compute_ipw_if_self_contained_one(...)^2) / n^2` on a synthetic binary static setup. 4 assertions at 1e-6.
- **T-non-static propensity correction** (Chunk 3g) — on shift / IPSI DGPs where Channel 1 ≠ 0, the full IF-based sandwich SE materially differs from the `J V_β Jᵀ` delta-only shortcut. The propensity correction typically *reduces* the per-intervention SE (the negative cross-term from the stacked M-estimation outweighs the added propensity-uncertainty term — same mechanism as Lunceford & Davidian 2004 for static binary). For shift the per-intervention SE drops ~5-8%; for IPSI the per-intervention effect is small but the *contrast* SE drops ~90% because the shared propensity model induces large positive covariance between the IPSI and natural-course marginal means. Proves the propensity-correction channel is non-trivially contributing to the SE under non-static interventions.
- **Bootstrap parity** (Chunk 3g) — on the same DGPs, `ci_method = "bootstrap"` agrees with `ci_method = "sandwich"` within 30% Monte Carlo tolerance (n_boot=299). `variance_bootstrap()` doesn't touch the IF engine, so this is an end-to-end sanity check on the Phase 4 runtime pipeline.

### 10. Stabilization, diagnostics, pushforward sign, IPSI shortcut

Four implementation notes that would otherwise get rediscovered mid-PR:

- **Stabilization: start nonstabilized.** Weights are `f_d(A|L) / f(A|L)`, no marginal numerator. Hájek normalization (`sum(wY)/sum(w)`) already controls finite-sample variance for the saturated MSM cases. A stabilization option (`stabilize = TRUE`) is a follow-up, not Phase 4 scope.
- **Weight diagnostics (minimal shim only in Phase 4).** `diagnose()` on a Phase 4 IPW fit ships a **minimal shim** that reads `fit$details$propensity_model` (instead of `fit$weights_obj`), computes default weights via `compute_density_ratio_weights(tm, data, static(1))` for binary fits or `compute_density_ratio_weights(tm, data, NULL)` for non-binary fits, and dispatches between treatment-type-specific panels via `detect_treatment_family()`. The shim covers the common binary static ATE case so existing `causat(estimator = "ipw") |> diagnose()` calls keep working. The **full intervention-aware / treatment-type-aware / estimand-aware / longitudinal-aware rewrite**, including a new `intervention =` argument mirroring `contrast()`, is scoped in [`PHASE_9_DIAGNOSE.md`](PHASE_9_DIAGNOSE.md) and deferred to a later phase. The deferral is intentional: a proper `diagnose()` rewrite touches every estimator and every treatment type and interacts with Phase 8 effect modification — too large to fit into the Phase 4 IPW engine push.
- **Pushforward sign + Jacobian (critical).** For continuous MTPs, the weight is `f_d(A_obs | L) / f(A_obs | L)` where `f_d` is the **pushforward** of `f` under the intervention — *not* `f(d(A_obs) | L) / f(A_obs | L)`. For `shift(δ)` this means evaluating the fitted density at `A_obs − δ` (because `d⁻¹(y) = y − δ`), and for `scale_by(c)` it means evaluating at `A_obs / c` and multiplying by the Jacobian `|1/c|`. The naive "evaluate at the intervened value" formula gives `E[Y^{shift(−δ)}]` instead of `E[Y^{shift(δ)}]` — it's a sign trap that is easy to write and hard to notice without the truth-based Hájek test. The `compute_density_ratio_weights()` and `make_weight_fn()` bodies both carry a comment block pointing at this derivation.
- **IPSI closed-form shortcut.** For `ipsi(δ)` the density ratio collapses to `w_i = (δ·A_i + (1 − A_i)) / (δ·p_i + (1 − p_i))`. Use this closed form inside the unified engine instead of evaluating the density at two transformed points — it is faster and avoids numerical near-1 ratios on both sides of the ratio. The `weight_fn` closure used by the variance engine evaluates the same closed form at candidate `α` values.

### 11. Chunk plan and test status

Phase 4 lands as a sequence of focused commits rather than one big-bang rewrite, so each step can be reviewed and rolled back independently:

| Chunk | Scope | Status |
|---|---|---|
| Foundation | `R/treatment_model.R` + `R/ipw_weights.R` + `R/interventions.R` IPSI constructor cleanup + `check_estimand_intervention_compat()` + threshold/continuous-dynamic rejection | **done** (commits `2bea1ae`, `128340b`, `8d79ceb`, `28eefb6`) |
| 3a | T-A_β_α hand-derived cross-derivative test on `make_weight_fn()` for static(1), static(0), shift, IPSI, NULL | **done** (commit `4090e51`, 12 assertions at 1e-6) |
| 3b | `compute_ipw_if_self_contained_one()` helper + stacked-sandwich oracle on a hand-assembled (alpha, beta) 3×3 bread/meat system | **done** (commit `6e3d42b`, 4 assertions at 1e-6) |
| 3c.i | `fit_ipw()` rewrite on `fit_treatment_model()` + `make_weight_fn()`; `propensity_model_fn` arg on `causat()`; `compute_contrast()` IPW path rewrite via density-ratio weights + per-intervention `Y ~ 1` Hájek MSM; `variance_if_ipw()` straight loop via `compute_ipw_if_self_contained_one()`; Phase 3 IPW test suite adapted; `ipw.md` / `s3-methods.md` snapshots regenerated | **done** (commit `d9732bf`) |
| 3c.ii | Delete `prepare_propensity_if()` / `prepare_propensity_if_weightit()` / `prepare_propensity_if_self_contained()` / `apply_propensity_correction()` from `R/variance_if.R`; unify the `variance_if()` dispatcher so IPW routes through a single live `variance_if_ipw()`; route `compute_contrast()`'s IPW path through the dispatcher uniformly; sweep Branch-A narration from roxygen | **done** (commit `b23660b`, ≈213 lines deleted) |
| 3c.iii | Unify bootstrap IPW path: rewrite `refit_ipw()` to replay `fit_ipw()` on resampled data, delete `ipw_boot_replicate()` special case, route bootstrap through `compute_ipw_contrast_point()` uniformly with the sandwich path; extract `diagnose()` IPW shim into `diagnose_ipw_point()` helper for the binary static ATE case; regenerate `_snaps/diagnose.md` (no-op); sweep orphaned Phase-3-era narration | **done** (commit `9055f1f`) |
| 3d | `tests/testthat/helper-ipw-weightit-oracle.R` + contrast-level oracle tests T-oracle1..4 (binary ATE/ATT/ATC + GAM propensity, `skip_if_not_installed("WeightIt")`); DESCRIPTION move WeightIt `Imports:` → `Suggests:`; sweep stale `WeightIt` roxygen/comments from `R/`. Also surfaced and fixed a latent correctness bug: the chunk 3c.i runtime silently returned the ATE for ATT / ATC fits because `compute_density_ratio_weights()` / `make_weight_fn()` did not consume `fit$estimand`. The fix threads `estimand` through both weight builders and adds `ht_bayes_numerator(estimand, tm, fit_data, family_tag)` — the unified Bayes-rule numerator `f*(L) = f(A* \| L)` derived in the `compute_density_ratio_weights()` roxygen header. ATE/ATT/ATC now match `WeightIt::glm_weightit()` to ~1e-6 on point estimates and sandwich SEs on the same propensity model; ATT bootstrap SE tracks the sandwich SE within Monte Carlo error (verified manually before committing). | **done** |
| 3e | Categorical (multinomial) branch in `fit_treatment_model()` + `make_weight_fn()` + `evaluate_density()`; multinomial-specific variance engine (`prepare_model_if_multinom()`); `nnet::multinom` as default categorical fitter; truth-based test (static) + smoke test (dynamic) | **done** |
| 3f | lmtp contrast-level oracles T-oracle5 (shift via `lmtp::lmtp_sdr()`, `skip_if_not_installed("lmtp")`) + T-oracle6 (manual IPSI weight parity). lmtp already in `Suggests:`. Also fixed a latent bug where `apply_intervention()` aborted for `ipsi()` in `compute_ipw_contrast_point()` and `variance_if_ipw()` — IPSI does not materialize a counterfactual treatment value, so the intercept-only MSM path now skips the treatment-column modification. | **done** |
| 3g | Non-static variance regression tests — T-non-static (propensity correction materially changes SE) + bootstrap parity for shift and IPSI | **done** |
| 3h | User-facing vignette `vignettes/interventions.qmd` — intervention-type tour with estimator-by-estimator examples; also rewrote `vignettes/ipw.qmd` to remove stale WeightIt references and reflect the Phase 4 self-contained engine | **done** |
| 3i | Theory vignette `vignettes/ipw-variance-theory.qmd` — density-ratio derivation, pushforward sign + Jacobian, HT indicator form, IPSI closed form, `make_weight_fn` closure design, numerical `A_{β,α}` verification. Cross-reference from `vignettes/variance-theory.qmd` §4.2. Also updated stale `WeightIt::glm_weightit()` references in `variance-theory.qmd` §4.2 to reflect the self-contained engine. | **done** |
| 3j | Count treatment (Poisson + negative binomial) density branch via explicit `propensity_family = "poisson"` / `"negbin"` opt-in on `causat()`; `fit_count_density()` + `evaluate_density()` Poisson/NB branches; count closure in `make_weight_fn()`; non-integer shift/scale rejection in `check_intervention_family_compat()`; truth-based Poisson DGP test + NB parity + rejection snapshot tests. See §5a for full design. | **done** |

The chunk boundary is deliberately where the runtime architecture flips (3c). Chunks 3a and 3b validate the new variance machinery against first-principles references **before** any runtime code changes, so chunk 3c can focus on plumbing without also debugging the math.

**Test file status after each chunk:**

- **Foundation** and **3a/3b**: add `test-treatment-model.R`, `test-ipw-weights.R`, `test-ipw-cross-derivative.R`, `test-ipw-branch-b.R`, `test-estimand-intervention-compat.R`. All existing Phase 3 tests still pass because nothing in runtime R/ code changed at the IPW path yet.
- **3c**: all existing Phase 3 IPW tests are either rewritten for the new architecture or deleted (see §12). Snapshots regenerated. Simulation/by-estimand tolerances audited and widened where necessary.
- **3d**: new `test-ipw-weightit-oracle.R` added. WeightIt no longer on the runtime path.
- **3e**: categorical abort test in `test-treatment-model.R` replaced with positive multinomial tests; categorical abort in `test-simulation.R` replaced with truth-based static ATE test + dynamic smoke test. New DGP `simulate_categorical_continuous()` in `helper-dgp.R`.
- **3f**: new `test-ipw-lmtp-oracle.R` + `helper-ipw-lmtp-oracle.R` added. T-oracle5 uses `lmtp::lmtp_sdr()` (SDR replaced the defunct `lmtp_ipw()` in lmtp >= 1.5.0). T-oracle6 uses a manual Kennedy (2019) closed-form weight oracle. Also fixed IPSI `apply_intervention()` abort in `R/ipw.R` and `R/variance_if.R` — IPSI skips treatment-column modification since the `Y ~ 1` MSM doesn't use it.
- **3g**: `test-ipw-variance-regression.R` added. T-non-static (shift + IPSI) and bootstrap parity (shift + IPSI). 4 test_that blocks, 10 assertions.
- **3h**: new `vignettes/interventions.qmd` created; `vignettes/ipw.qmd` rewritten; `vignettes/introduction.qmd` updated. No test file changes.
- **3i**: new `vignettes/ipw-variance-theory.qmd` created; `vignettes/variance-theory.qmd` cross-references added and stale WeightIt references updated. No test file changes.
- **3j**: new `test-ipw-count.R` added. `simulate_count_treatment()` DGP in `helper-dgp.R`. Truth-based Poisson shift test, NB parity test, 6 rejection snapshot tests, unit tests for `fit_treatment_model()` and `evaluate_density()` Poisson/NB branches, `scale_by(1)` no-op positive test.

### 12. Phase 3 tests that need attention in Chunk 3c

The reconnaissance pass (between commit 6e3d42b and Chunk 3c) identified these test touchpoints:

**WeightIt-specific tests to DELETE (behavior no longer applies):**

- `test-ipw.R`: "forwards `method` to `WeightIt::weightit()`" (tests `method = "cbps"` routing through the Phase 3 WeightIt path; the whole `method = ...` mechanism is gone under Phase 4's `propensity_model_fn` approach).
- `test-ipw.R`: non-Mparts warning test ("WeightIt method 'gbm' does not implement the M-estimation correction"); warning no longer exists because WeightIt is no longer called at runtime.
- `test-critical-review-2026-04.R`: B2 test ("bootstrap refit replays user's `method = 'cbps'`"); entire test is about WeightIt cbps bootstrap replay.

**Tests to REWRITE (same coverage, new architecture):**

- `test-ipw.R`: rewrite basic fit / contrast sanity tests to use the new `fit$details$propensity_model` + `propensity_model_fn` shape. Add new tests for the `propensity_model_fn = mgcv::gam` path.
- `test-diagnose.R`: update every test that reads `fit$weights_obj$ps` / `$weights` / `$treat` / `treat.type` to use the new slots and the shimmed `diagnose()`.
- `test-s3-methods.R`: update `print.causatr_fit` / `summary.causatr_fit` IPW path tests (print output text changes because the fit no longer stores a `weightit` object).
- `test-replay-fit.R`: IPW sections; gcomp sections stay unchanged because `replay_fit()` is still used by `refit_gcomp()`.
- `test-weights-edge-cases.R`: IPW path — external weights now enter the propensity density fit via `weights = `, not WeightIt's `s.weights`.

**Tests to AUDIT (expect to pass, may need tolerance widening):**

- `test-simulation.R`, `test-by-estimand.R`, `test-variance-if.R`, `test-complex-dgp.R`: truth-based simulation tests for binary/continuous outcomes × ATE/ATT/ATC × sandwich/bootstrap. Both old and new paths estimate the same functional; point estimates should agree at ~1e-6 and SEs should agree at ~1e-4 (finite-sample drift from numDeriv central differences vs WeightIt's analytic Jacobian). Any test needing wider tolerance than ~1e-3 gets flagged and investigated.
- `test-multivariate.R`: IPW rejection tests for multivariate treatment — should still fire, no changes needed.
- `test-contrast.R`: IPW rejection tests for non-static interventions — the existing `check_interventions_compat()` gate stays in place, so these tests still fire. Update expected error messages if they reference WeightIt.

**Snapshots to REGENERATE:**

- `_snaps/ipw.md` — WeightIt error text that appears in snapshot assertions.
- `_snaps/diagnose.md` — `fit$weights_obj`-based output format.
- `_snaps/s3-methods.md` — print/summary text for IPW fits.

**Independent oracle cross-checks** (see §9):

- **T-oracle1 / T-oracle2**: `WeightIt` equivalence on static binary, wrapped in `skip_if_not_installed("WeightIt")`.
- **T-oracle3 / T-oracle4**: `lmtp::lmtp_ipw()` point-estimate equivalence on shift / IPSI, wrapped in `skip_if_not_installed("lmtp")`.

**Estimand gating:**

- `contrast()` with an ATT or ATC fit + a non-static intervention aborts with error class `causatr_bad_estimand_intervention` (snapshot).

**Intervention × family rejection:**

- `threshold()` under IPW aborts at contrast time with a pointer to `estimator = "gcomp"` (snapshot).
- `dynamic()` on continuous treatment under IPW aborts with a pointer to `gcomp` or smooth MTPs (snapshot).
- `static(v)` on continuous treatment under IPW aborts with a pointer to `shift()` / `scale_by()` (snapshot).
- `shift()` / `scale_by()` / `ipsi()` on non-matching treatment families abort with a family-specific pointer (snapshot).

## Items

**Core engine**

- [x] `R/treatment_model.R` — fit treatment density model via `propensity_model_fn`; `evaluate_density()` for any treatment value. Binary + continuous landed in the foundation chunk; categorical (multinomial) landed in chunk 3e via `nnet::multinom` + `evaluate_categorical_density()` + `prepare_model_if_multinom()` for sandwich variance.
- [x] `R/ipw_weights.R` — density-ratio weight computation + `make_weight_fn()` closure factory. Three branches: HT indicator (discrete point-mass interventions), smooth pushforward with correct sign + Jacobian (continuous MTPs), IPSI closed form. `check_intervention_family_compat()` enforces the intervention × family compatibility matrix in §4.
- [x] Rewrite `R/ipw.R` — single self-contained engine handling static + shift + scale_by + dynamic + ipsi via density-ratio weights and an explicit weighted MSM. Populates `fit$details$propensity_model`, `fit$details$treatment_model`, `fit$details$weight_fn`, `fit$details$propensity_model_fn`. WeightIt runtime call removed (chunk 3c.i, commit `d9732bf`).
- [x] Add `propensity_model_fn` argument to `causat()` (default = `model_fn`, per the **Option B** decision in §4); plumbed through to `fit_ipw()` (chunk 3c.i).
- [x] IPW Channel-2 IF lives in `compute_ipw_if_self_contained_one()` (chunk 3b, commit `6e3d42b`); the `variance_if()` dispatcher routes IPW to a single straight loop `variance_if_ipw()` and the Branch A primitives (`prepare_propensity_if*`, `apply_propensity_correction()`) are deleted (chunk 3c.ii, commit `b23660b`).
- [x] Unblocked `apply_single_intervention()` for `ipsi` (`R/interventions.R`) — the constructor `rlang::inform()` was dropped during foundation cleanup; the IPW engine handles IPSI via the closed-form weight branch in `compute_density_ratio_weights()` and never calls `apply_single_intervention()` on the treatment column itself. The remaining ipsi branch in `apply_single_intervention()` aborts with a flat error pointing at `estimator = "ipw"` for callers that try to materialize a counterfactual treatment column under gcomp / matching.

**Estimand + safety checks**

- [x] `check_estimand_intervention_compat()` — rejects ATT/ATC for non-static interventions at contrast time (error class `causatr_bad_estimand_intervention`); landed in the foundation chunk (commit `28eefb6`).

**Categorical treatment + IPSI**

- [x] Categorical treatment support across checks + IPW path via multinomial density model. Landed in chunk 3e: `fit_categorical_density()`, `evaluate_categorical_density()`, categorical HT closure in `make_weight_fn()`, `prepare_model_if_multinom()` in the variance engine, `nnet` added to `Imports:`.
- [x] Count treatment (Poisson / negative binomial) density branch via explicit `propensity_family` opt-in. **Sub-chunk 3j.** `fit_count_density()` for Poisson (`stats::poisson()`) and NB (`MASS::glm.nb`), `evaluate_density()` Poisson/NB branches, count closure in `make_weight_fn()` with log-link `lambda = exp(X %*% alpha)`, non-integer shift/scale rejection in `check_intervention_family_compat()`, `propensity_family` parameter threaded from `causat()` through `fit_ipw()` to `fit_treatment_model()` with MASS::glm.nb auto-selection for negbin.
- [x] IPSI implementation wiring through `fit_ipw()` — the closed-form weight in `R/ipw_weights.R` is consumed by the per-intervention MSM refit in `compute_contrast()`'s IPW path (chunk 3c.i).

**Dependencies**

- [x] Move `WeightIt` from `Imports:` to `Suggests:` in `DESCRIPTION`; update `R/causatr-package.R` `@importFrom` tags accordingly. (Chunk 3d — there were no `@importFrom WeightIt` tags in causatr-package.R; the stale roxygen prose in `R/causat.R`, `R/contrast.R`, `R/utils.R`, `R/checks.R`, and the `IPW (WeightIt)` print label in `R/print.R` were swept in the same chunk.)
- [x] `lmtp` already in `Suggests:` (added prior to chunk 3f).

**Diagnostics + docs**

- [x] `diagnose()` weight summaries for the new engine — minimal shim (`diagnose_ipw_point()` in `R/diagnose.R`) covering the binary static ATE case (PS distribution, treated/control/overall HT weight summary with mean / sd / min / max / ESS) shipped in chunk 3c.iii (commit `9055f1f`). Non-binary treatment families abort with `causatr_diag_unsupported_family`. Extreme-weight count + percentile truncation are part of the full Phase 9 rewrite.
- [x] Update `FEATURE_COVERAGE_MATRIX.md` — collapse the two IPW columns into one, add rows for shift / scale_by / dynamic / ipsi under IPW, add estimand rejection rows, add `threshold()` rejection row under IPW, add T-oracle3/4 rows. (Initial pass in the foundation chunk; the final pass lands alongside the `fit_ipw()` rewrite.)
- [x] Vignette: `interventions.qmd` — shift, scale, dynamic (binary), IPSI examples. `threshold()` is documented under the gcomp vignette only. Also rewrote `ipw.qmd` to remove all stale WeightIt references and reflect the Phase 4 architecture.

**Tests** (see §11 for the full list)

- [x] Unit tests for `R/treatment_model.R` and `R/ipw_weights.R` (foundation chunk — 79 assertions across both files).
- [x] T-oracle1, T-oracle2, T-oracle3 (WeightIt static binary ATE/ATT/ATC), `skip_if_not_installed("WeightIt")` (chunk 3d, `test-ipw-weightit-oracle.R`).
- [x] T-oracle4 (GAM propensity + WeightIt GLM-PS soft match), `skip_if_not_installed("WeightIt")` + `skip_if_not_installed("mgcv")` (chunk 3d).
- [x] T-oracle5 (lmtp SDR shift point-estimate parity, `skip_if_not_installed("lmtp")`), T-oracle6 (manual IPSI weight parity). Chunk 3f. `lmtp_ipw()` was removed in lmtp >= 1.5.0; T-oracle5 uses `lmtp_sdr()` with parametric learners as the point-estimate oracle.
- [x] T-A_β_α hand-derived vs numerical (chunk 3a, commit `4090e51`, `test-ipw-cross-derivative.R`, 12 assertions at 1e-6).
- [x] T-end-to-end stacked sandwich by hand (chunk 3b, commit `6e3d42b`, `test-ipw-branch-b.R`, 4 assertions at 1e-6).
- [x] T-non-static propensity correction materially changes SE. **Chunk 3g.** Shift: per-intervention SE drops ~5-8%. IPSI: contrast SE drops ~90% via off-diagonal covariance through shared propensity model. `test-ipw-variance-regression.R`.
- [x] Bootstrap parity test (sandwich vs bootstrap on the same shift / IPSI DGP). **Chunk 3g.** `ci_method = "bootstrap"` (n_boot=299) agrees with `ci_method = "sandwich"` within 30% Monte Carlo tolerance. `test-ipw-variance-regression.R`.
- [x] Estimand × intervention rejection snapshot (`test-estimand-intervention-compat.R`, foundation chunk).
- [x] T-count-poisson: truth-based Poisson count treatment DGP + integer `shift(1)` recovers ATE = 1.5. **Chunk 3j.** `test-ipw-count.R`.
- [x] T-count-negbin: NB fit on same DGP agrees with Poisson within tolerance (NB nests Poisson). **Chunk 3j.** `test-ipw-count.R`.
- [x] T-count-reject: non-integer `shift(0.5)`, non-integer-preserving `scale_by(2)`, `static()`, `dynamic()`, `threshold()`, `ipsi()` on count treatment all abort with snapshot. **Chunk 3j.** `test-ipw-count.R`.

**Deferred (explicitly not Phase 4 scope)**

- [ ] `threshold()` under IPW — rejected architecturally, not deferred. The intervention is well-defined under `gcomp`; there is no meaningful density-ratio path.
- [ ] `threshold()` on count treatments — same mixed-measure problem as continuous. Rejected.
- [ ] `dynamic()` on count treatments — deterministic rule on a discrete-but-large support. Could work via HT indicators if the support is enumerable, but the cost–benefit vs g-comp is poor. Rejected for now; revisit if demand arises.
- [ ] Stabilized weights (`stabilize = TRUE` option).
- [ ] IPW for time-varying treatments — extend the same density-ratio machinery to longitudinal IPW (pooled-over-time MSM); deferred from Phase 5.
