# Phase 23 — Categorical (Multinomial & Ordinal) Outcome Models

> **Status: IN PROGRESS** — 23a-1 shipped (point g-computation, multinomial,
> bootstrap variance, full S3 layer). Remaining chunks below are PENDING.
>
> **Depends on:** Phase 2 (point gcomp), Phase 13 (extended outcome types)

## Scope

Support a **single categorical outcome column** `Y` — one factor with K > 2
levels (e.g. `Y ∈ {none, mild, severe}`). This is **not** multivariate
outcomes (not several outcome variables): the result becomes a K-vector per
intervention only because one categorical outcome has K class probabilities
`P(Y = k | do(A = a))` that sum to 1, exactly as a binary logistic outcome
implicitly reports `P(Y = 1)` and `P(Y = 0)`. It is the outcome-side mirror of
the multinomial-*treatment* machinery causatr already has (`nnet::multinom`
propensity, `prepare_model_if_multinom()`).

- **23a — Multinomial** outcomes via `nnet::multinom`. Estimand:
  `P(Y = k | do(A = a))` per class; per-class difference / ratio / OR contrasts.
- **23b — Ordinal** outcomes via `MASS::polr` (`ordinal::clm` optional).
  Estimand: cumulative `P(Y ≥ k | do(A = a))`. `predict(polr, type = "probs")`
  returns the same n×K shape, so 23b reuses the entire 23a result schema + S3
  layer; only the detection predicate and the variance Hessian differ.

### Out of scope

- Survival / time-to-event outcomes (owned by `survatr`).
- Zero-inflated models (custom predict + variance).
- Multivariate (multiple) outcome variables.

## What is orthogonal (free) vs what needs work

A categorical outcome touches only three things: the **prediction shape**
(`predict(model, type = "probs")` is an n×K matrix, not a length-n vector), the
**result schema** (per-class rows), and the **variance + S3 layer**.

**Orthogonal — rides inside the per-class loop with no extra design:** treatment
type (binary / continuous / categorical / count / multivariate — only changes
`apply_intervention()` + the design matrix), intervention family (static / shift
/ scale_by / threshold / dynamic — only changes the counterfactual treatment
column), estimand (ATE / ATT / ATC / by / subset — only changes target-row
selection), contrast scale (difference / ratio / OR — per-class delta on the
within-class vcov block), IPCW (a fit-row filter + weight; prediction is still
n×K), and external weights / clustering.

**Needs genuine new work:** the K-valued result schema, the variance, and the
S3 layer.

## Combination classification

| Estimator | Timing | Decision | Chunk |
|---|---|---|---|
| gcomp | point | **Supported** — bootstrap (23a-1), sandwich (23a-2) | 23a-1 / 23a-2 |
| gcomp | longitudinal (ICE) | Deferred (per-class pseudo-outcome recursion) | 23a-5 |
| ipw | point | Deferred (per-class Hájek mean of `I(Y=k)`) | 23a-3 |
| ipw | longitudinal | Deferred | 23a-6 |
| aipw | point | Deferred (probs + per-class residual correction) | 23a-4 |
| aipw | longitudinal | Deferred | 23a-7 |
| gcomp / ipw / aipw | transport (`target=`) | Deferred | 23a-8 |
| point `stochastic()` | — | Deferred (MC-marginalisation over matrix preds) | 23a-4b |
| matching | point | Deferred (K parallel matched MSMs) | 23a-9 |
| **snm** | any | **Rejected by design** (`causatr_snm_categorical_outcome`) | — |
| ordinal (`polr`/`clm`) | any of the above | Deferred (reuses 23a schema) | 23b-* |

## Chunk plan (dependency-ordered)

- **23a-1 (SHIPPED)** — vector result schema + point-gcomp multinomial probs +
  **bootstrap** variance + full S3 layer + classed gates. Composes (validated)
  with binary / continuous / categorical treatment, static / shift
  interventions, K = 3 and K = 4, difference / ratio / OR, ATT, by-strata,
  external weights, subset, IPCW, ≥ 3 interventions, spline confounders.
- **23a-2** — analytic IF sandwich for point-gcomp multinomial. Reuse
  `prepare_model_if_multinom()` (`R/variance_if_core.R`) as the outcome-score
  template; add the softmax marginal-mean gradient. Oracle: the exact
  `marginaleffects::avg_predictions()` delta-method SE (already the point oracle
  in 23a-1) and a Python `delicatessen` M-estimation stack. Lift the
  `causatr_categorical_outcome_sandwich` gate.
- **23a-3** — IPW point, per-class Hájek (`I(Y=k) ~ 1` MSM per class). [23a-1]
- **23a-4** — AIPW point, per class (DR: probs + `I(Y=k) − p_k` correction).
  [23a-1, 23a-3]
- **23a-4b** — point `stochastic()` × multinomial (MC-marginalisation over the
  matrix predictions). [23a-1]
- **23a-5** — ICE / longitudinal gcomp categorical Y (per-class pseudo-outcome
  chains). [23a-1]
- **23a-6 / 23a-7** — longitudinal IPW / AIPW. [23a-3/4 + 23a-5]
- **23a-8** — transport. [23a-1]
- **23a-9** — matching (per-class matched MSM); lift the gate. [23a-1]
- **23b-1** — point-gcomp ordinal (`MASS::polr`; `ordinal::clm` as Suggests).
  Cumulative `P(Y ≥ k | do(a))`. Inherits the 23a-1 schema + S3 layer. [23a-1]
- **23b-2** — ordinal sandwich (`polr` Hessian). [23b-1, 23a-2]
- **23b-3** — ordinal under IPW / AIPW / ICE. [corresponding 23a + 23b-1]

## 23a-1 design (shipped)

**Additive, conditional generalisation.** The scalar `estimates` / `contrasts`
/ `vcov` path is untouched. A multinomial outcome is detected post-fit
(`is_multinom_outcome(model)` → `inherits(model, "multinom")`, `R/family.R`) and
`compute_contrast()` branches early to `compute_contrast_multinom()`
(`R/contrast.R`), so the scalar code object is never executed for a multinomial
fit. A scalar `causatr_result` stays byte-identical: no `class` column, `vcov` a
matrix, `class_labels = NULL`.

- **Prediction:** `predict_outcome()` (`R/contrast.R`) returns the n×K
  probability matrix for a multinom model (mirroring the K = 2 / single-row
  reshape of `evaluate_categorical_density()`), and the length-n vector
  otherwise.
- **Result schema:** `estimates` / `contrasts` gain a `class` column (K rows per
  intervention / contrast); `vcov` is a per-class named list of `k_int × k_int`
  matrices; `class_labels` records the K levels.
- **Contrasts:** the scalar delta-method block was extracted into
  `compute_pairwise_contrasts()` (`R/contrast.R`), called by both the scalar and
  the per-class paths so the math cannot diverge. OR semantics are one-vs-rest
  per class.
- **Variance:** `variance_bootstrap_multinom()` (`R/variance_bootstrap.R`)
  flattens the `k_int × K` surface in class-major order, reuses
  `process_boot_results()`, and slices per-class vcov blocks. `refit_gcomp()`
  replays `nnet::multinom` unchanged.
- **S3:** `print` / `summary` / `tidy` / `coef` / `confint` / `plot` /
  `knit_print` branch on `class` column presence / `class_labels`. `plot()`
  facets the forest plot by class (`forrest` `section = "class"`); when `by` is
  also present it facets by class with a warning (combined two-way faceting +
  nested vcov are deferred polish).
- **Gates** (`R/checks.R`, `check_categorical_outcome()`, run in `causat()`
  before `prepare_data()`):
  `causatr_snm_categorical_outcome` (SNM, by design),
  `causatr_categorical_outcome_unsupported` (matching / ipw / aipw /
  longitudinal / transport, deferred), and at `contrast()`
  `causatr_categorical_outcome_sandwich` (sandwich, until 23a-2) and
  `causatr_categorical_outcome_unsupported` (stochastic intervention, until
  23a-4b).

### Oracle / test strategy (`test-gcomp-categorical-outcome.R`)

1. **Known-softmax DGP** (`helper-multinom-oracle.R`): `Y ~ Categorical(softmax)`,
   with the g-computation truth `E_L[softmax_k(eta(a, L))]` computed by large-n
   Monte Carlo. The shift-MTP truth marginalises over the observed `(A, L)`.
2. **External oracle — `marginaleffects::avg_predictions()`** on causatr's *own*
   fitted `nnet::multinom` (exact point parity, ~1e-15), plus
   `avg_comparisons()` for the difference contrasts. This same delta-method SE
   is the analytic-SE target the deferred 23a-2 sandwich will match.
3. Per-class ratio / OR are exact functions of the means; bootstrap SE within a
   loose ratio of the delta-method SE (the tight equality is 23a-2's job).
4. Byte-identical guard for a scalar binomial result; all rejection gates.

### Rejection list (classed errors)

- **`causatr_snm_categorical_outcome`** — SNM estimates an additive
  structural-mean blip on a real-valued outcome; a K-simplex of class
  probabilities has no single additive blip. Rejected at `causat()`.
- **`causatr_categorical_outcome_unsupported`** — categorical outcome under an
  estimator / timing not yet shipped (matching / IPW / AIPW / longitudinal /
  transport, and stochastic interventions). Each later chunk lifts its slice.
- **`causatr_categorical_outcome_sandwich`** — analytic sandwich is 23a-2; use
  bootstrap for now (the betareg-ICE / G-LMTP 22b-core bootstrap-first
  precedent).
