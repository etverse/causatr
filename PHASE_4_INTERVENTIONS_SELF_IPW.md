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
| `shift()`     | ✓ | ✓ (continuous) | — |
| `scale_by()`  | ✓ | ✓ (continuous) | — |
| `threshold()` | ✓ | ⛔ (use gcomp — pushforward has point masses) | — |
| `dynamic()`   | ✓ | ✓ (binary / categorical only) | — |
| `ipsi()`      | — | ✓ (binary) | — |

There is no longer a separate "IPW (WeightIt)" vs "IPW (self-contained)" column — there's just IPW.

The two ⛔ / "family restriction" notes reflect a genuine architectural limit of the density-ratio framework, not a scope decision:

- **`threshold()` under IPW** is rejected because the pushforward of a continuous `f(a|l)` under a boundary clamp is a mixed measure — continuous density on `(lo, hi)` plus point masses at the boundaries — so the density ratio w.r.t. the fitted Gaussian is not well-defined. `gcomp` handles it cleanly via `predict(outcome_model, newdata = clamped)` and is the correct tool for this intervention.
- **`dynamic()` under IPW** is only defined on binary / categorical treatments, where the HT indicator weight `I(A_obs = rule_i) / f(A_obs | L)` is well-posed. Deterministic rules on continuous treatments are a Dirac per individual, which has the same pushforward-degeneracy problem as `threshold()` — users should either rewrite the rule as a smooth `shift()` / `scale_by()`, or switch to `gcomp`.
- **`static(v)` under IPW** is rejected on continuous treatments (nobody is observed exactly at `v`, so the HT indicator is zero almost surely).
- **`shift()` / `scale_by()`** are continuous-only because their density-ratio formulas rely on the Gaussian pdf; they have no meaning on a binary treatment.
- **`ipsi()`** is binary-only — Kennedy (2019)'s closed form is Bernoulli-specific.

All of these are enforced at contrast time by `check_intervention_family_compat()` in `R/ipw_weights.R`, which aborts with an actionable pointer.

### 5. Categorical treatment support

- G-comp: already works (factor treatment column in the outcome GLM).
- IPW: multinomial treatment density model via `nnet::multinom` or `VGAM::vglm` plugged in through `propensity_model_fn`; categorical density evaluation in `evaluate_density()`.
- Matching: stays binary-only (MatchIt limitation; we already error out cleanly in `fit_matching()`).
- Update `check_estimand_trt_compat()` to allow categorical + ATE.

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

### 9. WeightIt and lmtp as independent test oracles (both Suggests)

After Phase 4 lands, `WeightIt` moves from `Imports:` to `Suggests:` in `DESCRIPTION`. `lmtp` is **added** to `Suggests:` in the same PR — it serves as a second independent oracle for the non-static interventions that WeightIt cannot reach. Tests that use either package are wrapped in `skip_if_not_installed("...")` so the package builds and tests cleanly without them.

**WeightIt — variance oracle for static binary.** The classical M-estimation correction for parametric IPW is what `WeightIt::glm_weightit()` + `sandwich::estfun(asympt = TRUE)` gives you. For the static binary case this is asymptotically equivalent to our unified engine's sandwich, so it's a tight numerical oracle:

- **T-oracle1:** `compute_density_ratio_weights()` on a static binary setup ≡ `WeightIt::weightit(method = "glm")$weights` to ~1e-12.
- **T-oracle2:** full-pipeline `causat(estimator = "ipw") |> contrast()` ≡ `WeightIt::glm_weightit()` + MSM by hand on the same fit, **point estimate / SE / vcov** to ~1e-6.

**lmtp — point-estimate oracle for non-static interventions.** `lmtp::lmtp_ipw()` uses density-ratio weights and accepts parametric learners (`learners_trt = "SL.glm"`, `folds = 1` for deterministic behaviour). For a given parametric propensity model it agrees with our engine on the **Hájek point estimate**. It does **not** serve as a variance oracle because its SE is built from the EIF rather than the M-estimation sandwich — different object, only asymptotically comparable on saturated binary.

- **T-oracle3:** `causat(shift(-5))` on NHEFS continuous treatment ≡ `lmtp::lmtp_ipw(shift = function(d, t) d[[t]] - 5, learners_trt = "SL.glm", folds = 1)` point estimate to ~1e-6.
- **T-oracle4:** `causat(ipsi(2))` on NHEFS binary treatment ≡ a manual IPSI weight + weighted mean comparison (no lmtp-side IPSI helper; we build the shift closure ourselves from the closed form) to ~1e-6.

The causatr variance story stays anchored to (a) T-A_β_α elementwise hand-derivation, (b) T-end-to-end stacked sandwich by hand, (c) T-non-static IF > delta-only shortcut, (d) bootstrap parity — none of which depend on lmtp.

### 10. Stabilization, diagnostics, pushforward sign, IPSI shortcut

Four implementation notes that would otherwise get rediscovered mid-PR:

- **Stabilization: start nonstabilized.** Weights are `f_d(A|L) / f(A|L)`, no marginal numerator. Hájek normalization (`sum(wY)/sum(w)`) already controls finite-sample variance for the saturated MSM cases. A stabilization option (`stabilize = TRUE`) is a follow-up, not Phase 4 scope.
- **Weight diagnostics.** `diagnose()` on an IPW fit adds weight-distribution summaries: min / max / quartiles / 99th percentile / count of weights > 99th percentile. Optional user-side truncation at a given percentile is exposed as an argument to `diagnose()`, not baked into the fit.
- **Pushforward sign + Jacobian (critical).** For continuous MTPs, the weight is `f_d(A_obs | L) / f(A_obs | L)` where `f_d` is the **pushforward** of `f` under the intervention — *not* `f(d(A_obs) | L) / f(A_obs | L)`. For `shift(δ)` this means evaluating the fitted density at `A_obs − δ` (because `d⁻¹(y) = y − δ`), and for `scale_by(c)` it means evaluating at `A_obs / c` and multiplying by the Jacobian `|1/c|`. The naive "evaluate at the intervened value" formula gives `E[Y^{shift(−δ)}]` instead of `E[Y^{shift(δ)}]` — it's a sign trap that is easy to write and hard to notice without the truth-based Hájek test. The `compute_density_ratio_weights()` and `make_weight_fn()` bodies both carry a comment block pointing at this derivation.
- **IPSI closed-form shortcut.** For `ipsi(δ)` the density ratio collapses to `w_i = (δ·A_i + (1 − A_i)) / (δ·p_i + (1 − p_i))`. Use this closed form inside the unified engine instead of evaluating the density at two transformed points — it is faster and avoids numerical near-1 ratios on both sides of the ratio. The `weight_fn` closure used by the variance engine evaluates the same closed form at candidate `α` values.

### 11. Tests to add

**Variance-engine hand-derived truth:**

- **T-A_β_α numerical vs analytic** on static binary: derive `∂ψ_β,i/∂α^T` by hand for logistic propensity + binary static weights, compare against `numDeriv::jacobian(psi_beta_bar, ...)` elementwise at ~1e-6.
- **T-end-to-end stacked sandwich** on a tiny simulated dataset (`n = 200`, `q = p = 2`, Gaussian outcome, logistic propensity): assemble the full `(q+p+1) × (q+p+1)` stacked bread/meat by hand, invert, read off the `(μ,μ)` entry, compare against `variance_if()` at ~1e-10.
- **T-non-static** (shift / MTP) on NHEFS where Channel 1 ≠ 0: IF-based variance must materially exceed `J V_β Jᵀ` (the delta-only shortcut) — proves Phase 4 is actually capturing the extra uncertainty.
- **Bootstrap parity**: `ci_method = "bootstrap"` on the same setups must agree with sandwich within Monte Carlo error — bootstrap is unchanged because `variance_bootstrap()` refits the whole pipeline and never touches `variance_if()`.

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

- [x] `R/treatment_model.R` — fit treatment density model via `propensity_model_fn`; `evaluate_density()` for any treatment value. **Binary + continuous landed in the foundation chunk**; categorical (multinomial) branch is the next sub-chunk and currently hard-aborts with `causatr_phase4_categorical_pending`.
- [x] `R/ipw_weights.R` — density-ratio weight computation + `make_weight_fn()` closure factory. Three branches: HT indicator (discrete point-mass interventions), smooth pushforward with correct sign + Jacobian (continuous MTPs), IPSI closed form. `check_intervention_family_compat()` enforces the intervention × family compatibility matrix in §4.
- [ ] Rewrite `R/ipw.R` — single self-contained engine handling static + shift + scale_by + dynamic + ipsi via density-ratio weights and an explicit weighted MSM. Must populate `fit$details$propensity_model`, `fit$details$alpha_hat`, and `fit$details$weight_fn`. Drop the WeightIt runtime call.
- [ ] Add `propensity_model_fn` argument to `causat()` (default = `model_fn`, per the **Option B** decision in §4); plumb through to `fit_ipw()`.
- [ ] Fill in the body of `prepare_propensity_if_self_contained()` in `R/variance_if.R` per §7.
- [ ] Unblock `apply_single_intervention()` for `ipsi` (`R/interventions.R`) — it currently hard-aborts — and drop the `ipsi()` "dead-end" `rlang::inform()` from the constructor.

**Estimand + safety checks**

- [ ] `check_estimand_intervention_compat()` — reject ATT/ATC for non-static interventions at contrast time (error class `causatr_bad_estimand_intervention`).

**Categorical treatment + IPSI**

- [ ] Categorical treatment support across checks + IPW path via multinomial density model. Unblocks the `categorical` branch of `fit_treatment_model()` / `evaluate_density()` / `make_weight_fn()`.
- [ ] IPSI implementation wiring through `fit_ipw()` (the closed-form weight itself already exists in `R/ipw_weights.R`).

**Dependencies**

- [ ] Move `WeightIt` from `Imports:` to `Suggests:` in `DESCRIPTION`; update `R/causatr-package.R` `@importFrom` tags accordingly.
- [ ] Add `lmtp` to `Suggests:` for the non-static point-estimate oracle tests.

**Diagnostics + docs**

- [ ] `diagnose()` weight summaries for the new engine (min / max / quartiles / 99th percentile / extreme-weight count + optional user truncation at a percentile).
- [x] Update `FEATURE_COVERAGE_MATRIX.md` — collapse the two IPW columns into one, add rows for shift / scale_by / dynamic / ipsi under IPW, add estimand rejection rows, add `threshold()` rejection row under IPW, add T-oracle3/4 rows. (Initial pass in the foundation chunk; the final pass lands alongside the `fit_ipw()` rewrite.)
- [ ] Vignette: `interventions.qmd` — shift, scale, dynamic (binary), IPSI examples. `threshold()` is documented under the gcomp vignette only.

**Tests** (see §11 for the full list)

- [x] Unit tests for `R/treatment_model.R` and `R/ipw_weights.R` (foundation chunk — 79 assertions across both files).
- [ ] T-oracle1, T-oracle2 (WeightIt, `skip_if_not_installed`)
- [ ] T-oracle3, T-oracle4 (lmtp, `skip_if_not_installed`)
- [ ] T-A_β_α hand-derived vs numerical
- [ ] T-end-to-end stacked sandwich by hand
- [ ] T-non-static IF > delta-only shortcut
- [ ] Bootstrap parity
- [ ] Estimand × intervention rejection snapshot

**Deferred (explicitly not Phase 4 scope)**

- [ ] `threshold()` under IPW — rejected architecturally, not deferred. The intervention is well-defined under `gcomp`; there is no meaningful density-ratio path.
- [ ] Stabilized weights (`stabilize = TRUE` option).
- [ ] IPW for time-varying treatments — extend the same density-ratio machinery to longitudinal IPW (pooled-over-time MSM); deferred from Phase 5.
