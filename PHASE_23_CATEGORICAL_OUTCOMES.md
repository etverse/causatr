# Phase 23 — Categorical (Multinomial & Ordinal) Outcome Models

> **Status: IN PROGRESS** — 23a-1 shipped (point g-computation, multinomial,
> bootstrap variance, full S3 layer). The analytic sandwich (23a-2) is split
> into three dependency-ordered sub-chunks that together cover **every**
> variance scenario the multinomial point-gcomp surface can present:
> **23a-2a** (complete-case, SHIPPED), **23a-2b** (survey/external weights,
> SHIPPED), **23a-2c** (IPCW, PENDING). All other chunks below are PENDING.
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
| gcomp | point | **Supported** — bootstrap (23a-1); sandwich complete-case (23a-2a) + weighted (23a-2b); IPCW sandwich pending (23a-2c) | 23a-1 / 23a-2a–c |
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
- **23a-2 — analytic IF sandwich for point-gcomp multinomial.** Reuse
  `prepare_model_if_multinom()` (`R/variance_if_core.R`) as the outcome-score
  template; add the softmax marginal-mean gradient. Split into three
  dependency-ordered sub-chunks so the weighted and IPCW derivations land
  separately and nothing is dropped (see the "23a-2 design" section below for
  the full per-scenario coverage):
  - **23a-2a (complete-case)** — `variance_if_gcomp_multinom()`
    (`R/variance_if_multinom.R`): per-class softmax marginal-mean gradient +
    Channel-1 / Channel-2 assembly. Covers binary / continuous / categorical
    treatment, static / shift interventions, K = 3 and K = 4, difference /
    ratio / OR, ATE / ATT / by-strata / subset, ≥ 3 interventions, spline
    confounders. Lifts the `causatr_categorical_outcome_sandwich` gate for the
    complete-case path; keeps narrow transitional gates for weighted and IPCW
    sandwiches. Oracle: a Python `delicatessen` M-estimation stack (tight SE)
    and `marginaleffects::avg_predictions()` (tight point; its delta-method SE
    omits the Channel-1 empirical-distribution term, so it is a point oracle,
    not the SE oracle).
  - **23a-2b (survey / external weights) — SHIPPED.** Generalises the
    multinomial score / bread to prior weights `w_i`: `prepare_model_if_multinom(weights=)`
    builds the weighted information `H = sum_i w_i (diag(p_i) - p_i p_i^T) kron
    X_i X_i^T` and sets `r_score = w_i`; `variance_if_gcomp_multinom(weights=)`
    weights Channel 1 (`(n/sum_w) w_i (p_k - mu)`) and the marginal-mean gradient
    (`(1/sum_w) sum_i w_i p_k (delta - p_{l+1}) X*`). `weights = NULL` (or `w == 1`)
    reproduces the 23a-2a complete-case sandwich byte-for-byte. Lifts the weights
    slice of the sandwich gate (keyed now on `fit$details$ipcw` only). Oracle: a
    weighted `delicatessen`-style stack (`multinom_gcomp_weighted_sandwich.py`),
    matched to ~1e-7 (estimates) / ~1e-8 (SEs), plus weighted sandwich-vs-bootstrap
    parity. [23a-2a]
  - **23a-2c (IPCW)** — multinomial-aware censoring cross-term: a
    `phi_bar_cens(gamma)` closure built on the stacked weighted multinomial
    score, its Jacobian in `gamma`, projected through the stacked outcome
    `h`. Lifts the IPCW + sandwich gate. Oracle: Monte-Carlo coverage (no
    clean external M-estimation oracle stacks a multinomial outcome with an
    IPCW censoring model — documented in the test). [23a-2b]
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
   `avg_comparisons()` for the difference contrasts. `avg_predictions()` is the
   tight *point* oracle; its delta-method SE conditions on the standardisation
   sample (coefficient uncertainty only) and so omits the Channel-1
   empirical-distribution term that causatr's full IF sandwich carries. The
   tight *SE* oracle for the sandwich is therefore the full `delicatessen`
   M-estimation stack (23a-2), not the `marginaleffects` SE.
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
- **`causatr_categorical_outcome_sandwich`** — analytic sandwich is 23a-2. From
  23a-2a the **complete-case** path is supported, and from 23a-2b the
  **survey/external-weighted** path; the error remains only as a transitional
  gate for the still-unlanded IPCW sandwich slice (→ 23a-2c, keyed on
  `fit$details$ipcw`), where it points the user to `ci_method = "bootstrap"`.
  The gate is removed entirely when 23a-2c ships.

## 23a-2 design (multinomial sandwich)

The estimand surface is the `(k_int x K)` matrix of `mu_{k,a} = P(Y = k | do(A = a))`.
The sandwich produces, **per outcome class `k`**, a `k_int x k_int` vcov over the
interventions — the exact shape the bootstrap already returns (`boot_res$vcov`
is a per-class named list), so the entire contrast / S3 layer is reused
unchanged. Only the variance source swaps.

### The two reused primitives

`nnet::multinom` parametrises K-1 non-reference linear predictors
`eta_l = X beta_l` (the reference class `l = 1` has `eta_1 = 0`). Stack the
parameter as `beta = (beta_2, ..., beta_K)`, length `(K-1) * p`.

1. **Score + bread — already built.** `prepare_model_if_multinom()`
   (`R/variance_if_core.R`) returns the stacked per-observation score
   `s_i = (I(Y_i = l) - p_{l,i})_l kron X_i` and the inverse information
   `B_inv = H^{-1}` with `H = sum_i (diag(p_i^{nr}) - p_i^{nr} p_i^{nr T}) kron X_i X_i^T`
   over the K-1 non-reference classes. It was written for the multinomial
   *propensity* model but the model class is identical, so the outcome model
   reuses it verbatim (complete-case; 23a-2b adds the `w_i` weighting).
2. **Correction — already built.** `apply_model_correction(prep, g)` forms
   `correction_i = n * (X_stacked[i,] %*% B_inv %*% g)` — the Channel-2
   nuisance correction for any sensitivity gradient `g`.

### The one new derivation — the softmax marginal-mean gradient

For class `k` (one of all K classes, **including** the reference) under
intervention `a`, the marginal mean is `mu_{k,a} = (1/n_t) sum_{i in target} p_{k,i}(a)`.
Its gradient with respect to the stacked `beta` has, in the block for
non-reference class `l` (l = 1..K-1, i.e. level `class_labels[l+1]`):

    g_{k,a}^{(l)} = (1/n_t) sum_{i in target} p_{k,i}(a) * (delta_{k,l+1} - p_{l+1,i}(a)) * X_i^*

where `p_{.,i}(a)` are the softmax probabilities at the counterfactual design
`X_i^*` (already computed as `preds_list[[a]]`), `delta_{k,l+1} = 1` iff class
`k` is the non-reference class `l+1`, and `X_i^*` is the counterfactual design
row. This is the softmax Jacobian `dp_k/deta_l = p_k(delta_{kl} - p_l)` chained
through `deta_l/dbeta_l = X_i^*`. The reference class (`k = 1`) is the
`delta = 0` case (`g_{1,a}^{(l)} = -(1/n_t) sum_i p_{1,i} p_{l+1,i} X_i^*`); it
falls straight out of the same formula, so no class is special-cased.

Per class `k` and intervention `a`:

- **Channel 1** (empirical-distribution sampling, the term `marginaleffects`
  omits): `Ch1_{k,a,i} = (n / n_t) (p_{k,i}(a) - mu_{k,a})` on target rows,
  0 elsewhere — the per-class analogue of `build_point_channel_pieces()`.
- **Channel 2** (nuisance): `apply_model_correction(prep, g_{k,a})$correction`.
- `IF_{k,a,i} = Ch1_{k,a,i} + Channel2_{k,a,i}`; then
  `vcov_from_if(IF_list_over_a, n, int_names)` gives the per-class
  `k_int x k_int` block.

`variance_if_gcomp_multinom()` (new, `R/variance_if_multinom.R`) loops this over
classes and returns the per-class vcov list. `compute_contrast_multinom()`
selects this list when `ci_method == "sandwich"` (and sets `boot_t = NULL`,
`use_perc = FALSE` so the Wald/delta bounds are used); the per-class contrast
math stays in the shared `compute_pairwise_contrasts()`.

### Full scenario coverage (so nothing is forgotten)

Every cell below is a multinomial point-gcomp **sandwich** scenario. "Rides"
means it is carried by an existing orthogonal mechanism and needs only a test,
not new variance code.

**23a-2a — complete-case (no weights, no IPCW).**

| Scenario | Mechanism | Test / oracle |
|---|---|---|
| binary trt, static, ATE, difference | core gradient | delicatessen SE + marginaleffects point |
| ratio / OR contrast | `compute_pairwise_contrasts()` delta on per-class vcov | delicatessen + exact-fn-of-means |
| continuous trt, shift MTP | `apply_intervention()` writes `A + delta`; `X^*` extrapolates | delicatessen + marginaleffects |
| categorical (k>2) treatment | design matrix only | delicatessen |
| K = 4 (and general K) | gradient loops over all K classes | delicatessen (schema) |
| ATE / ATT / by-strata / subset | `target_idx` + Channel-1 `(n / n_t)` normalisation | treated/stratum/subset delicatessen or marginaleffects |
| >= 3 interventions | per-intervention IF columns | sandwich-vs-bootstrap parity |
| spline / `I(.)` confounders | `iv_design_matrix()` term expansion | delicatessen |
| sandwich vs bootstrap SE agreement | both variance paths on one DGP | internal MC cross-check (large n) |
| scalar (binomial/gaussian) outcome | early branch never taken | byte-identical guard (unchanged) |

**23a-2b — survey / external weights (SHIPPED).** Weighted score `w_i s_i`,
weighted information `H_w = sum_i w_i (...) kron X_i X_i^T`; Channel 1 and the
marginal-mean gradient weight by `w_i` and normalise by `sum_w`. `weights = NULL`
collapses to 23a-2a. Oracle: weighted `delicatessen`-style stack
(`multinom_gcomp_weighted_sandwich.py`) + weighted sandwich-vs-bootstrap parity.

| Scenario | Mechanism | Test / oracle |
|---|---|---|
| binary trt, static, weighted, means + diff/ratio/OR | weighted bread/score/Ch1/grad | weighted delicatessen stack ~1e-7 |
| weighted sandwich vs weighted bootstrap | both paths, one DGP | SE ratio in (0.9, 1.1) |
| `weights = NULL` ≡ ones | weighted formulas collapse | per-class vcov identical to 1e-12 |

**23a-2c — IPCW (missing Y).** Stacked multinomial censoring cross-term added to
the per-class IF; reuses the `make_ipcw_weight_fn()` / `numDeriv::jacobian`
plumbing of `compute_ipcw_if_correction()` but over the stacked multinomial
score (`phi_bar_cens(gamma)` returns a `(K-1)*p` vector). Oracle: Monte-Carlo
coverage of the MAR-debiased estimand (documented; no external stack covers
multinomial-outcome + IPCW M-estimation).

### Still owned by other chunks (not part of 23a-2)

- **Stochastic interventions** under multinomial — 23a-4b (MC-marginalisation
  over the matrix predictions; the sandwich would follow the gcomp
  stochastic-gradient MC path of `build_point_channel_pieces()`).
- **Transport (`target =`)** — 23a-8.
- **IPW / AIPW / matching / longitudinal / SNM** — their own 23a / 23b chunks;
  `causat()` still rejects them for a categorical outcome.

### Transitional classed gates (removed as each slice lands)

- `causatr_categorical_outcome_sandwich` — raised at `contrast()` when
  `ci_method == "sandwich"` **and** `isTRUE(fit$details$ipcw)` (until 23a-2c).
  Complete-case (23a-2a) and survey/external-weighted (23a-2b) sandwiches no
  longer raise it. The gate keys on the explicit `ipcw` flag rather than on
  weight presence because IPCW also populates `fit$details$weights` (the
  combined survey × IPCW weight). Removed entirely when 23a-2c ships.
