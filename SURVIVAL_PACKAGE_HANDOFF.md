# Handoff: Causal Survival Analysis as a Standalone etverse Package

> **Purpose.** This document is the full hand-off to start a new etverse
> package dedicated to causal survival analysis. It replaces Phase 7 of
> `causatr`. Survival analysis is being removed from `causatr` (scope
> was "causal effect estimation for methodological triangulation on
> scalar / longitudinal outcomes") and moved into its own package that
> **uses causatr as the engine** for outcome-modelling, IPW, variance,
> interventions, and diagnostics.
>
> **Audience.** A fresh Claude Code chat in the new survival-package
> repo. Treat this as the single source of truth for scope, design, and
> the parts of causatr you can depend on.

## 1. Why a separate package

Causal survival analysis differs from scalar-outcome causal inference
in ways that cleanly delimit a package boundary:

1. **Data shape.** Time-to-event data is almost always censored. The
   natural runtime is the person-period (long) table with an event
   indicator at each interval, not a single row per individual.
2. **Estimand shape.** The target is a **curve** (survival, cumulative
   incidence, RMST up to $t^*$) or a contrast of curves — not a
   scalar. The `causatr_result` S3 is scalar-shaped; survival results
   need a time-indexed shape.
3. **Competing risks.** Cause-specific hazards + cumulative incidence
   decomposition is a first-class requirement that only exists here.
4. **Pooled logistic / hazard link.** The outcome model is a
   discrete-time hazard (pooled logistic, following Hernán & Robins
   Ch. 17), not a regression on a continuous outcome.
5. **Variance.** Sandwich variance for a survival curve needs a
   cross-time delta-method chain through the cumulative product — a
   different aggregation than anything in `causatr::variance_if()`.
6. **Matching is rejected.** `MatchIt` weights + pooled logistic on
   person-period data is awkward and duplicative of what
   `survival::coxph(..., weights=, cluster=)` already does. Out of
   scope.

Keeping all of this inside `causatr` would either bloat the API
(survival branches inside `contrast()`, curve-shaped results inside a
scalar-shaped class, new S3 methods everywhere) or stall progress on
scalar estimators that are better served by a focused package.

## 2. What already exists in causatr (pre-removal baseline)

Before removal, causatr shipped a **scaffolded** Phase 7 with:

- `causat_survival(data, outcome, treatment, confounders, id, time, censoring, competing, time_formula, weights, ...)`
  — fit-only. Converts to person-period if needed, builds the risk set
  (drops rows at/after first event per id, drops rows at/after first
  censor per id), fits a pooled logistic regression
  $\text{logit}\,h(t \mid A, L) = \alpha(t) + \beta_A A + \beta_L L$
  on at-risk rows. Returns a `causatr_fit` with `type = "survival"`.
- `to_person_period(data, id, time_varying, time_invariant, time_name)`
  — wide→long reshape (kept in causatr, general-purpose).
- `contrast()` on a survival fit aborts with "not implemented".
- `causat_survival(competing = ...)` aborts: the argument was reserved
  but would have silently fit a biased cause-deleted hazard.
- Risk-set bookkeeping used `.causatr_prev_event` / `.causatr_prev_cens`
  within-id lagged cumsum columns.
- `T6`: weighted fits switch to `quasibinomial()` (same score equations,
  free dispersion, drops the "non-integer #successes" warning).
- `T9`: internal columns stripped from `fit$data` before return.
- Smoke tests pinned fit-time behaviour.

**That entire surface is being removed from causatr.** Section 6
summarises what ports over to the new package and what needs to be
written fresh.

## 3. Target scope for the new package

### Two tracks + competing risks

1. **Track A — Point survival.** Baseline treatment, time-to-event
   outcome. Pooled logistic hazard on person-period data. This is the
   Ch. 17 workflow in Hernán & Robins.
2. **Track B — Longitudinal survival (ICE hazards).** Time-varying
   treatment, time-to-event outcome. Iterated conditional hazards
   (Zivich et al. 2024 extended to the hazard link): per-step logistic
   fit at the final period, backward iteration with the **survival
   tail** $\tilde{Y}_k = 1 - \hat{S}^d_{k:K}$ as pseudo-outcome, fit
   with `quasibinomial` at earlier steps.
3. **Competing risks.** Cause-specific hazards + cumulative incidence
   function (CIF) across both tracks. Fine–Gray / subdistribution
   hazards are **out of scope** (different estimand, different data
   structure — document the choice).

### Matching + survival: out of scope

Matching on survival outcomes is covered well by
`survival::coxph(..., weights = match_weights, cluster = subclass)`
directly on the `MatchIt` output. The new package should **hard-abort**
any matching + survival call with a clear error pointing users to
`gcomp` / `ipw` (or to `coxph` directly).

## 4. Inheritance principle: survival is an outcome-model swap

**Design thesis.** Survival is *not* a new treatment-side axis. It is a
change on the **outcome side** (pooled logistic hazard instead of a
scalar GLM on Y). So nearly every treatment-side combination that
causatr already supports inherits into survival for free:

| Feature axis | Point survival (A) | Longitudinal survival (B) | Source in causatr |
|---|---|---|---|
| Binary / continuous / categorical / count treatment (gcomp) | inherits | inherits | `R/gcomp.R`, `R/ice.R` |
| Binary / continuous / categorical / count treatment (IPW) | inherits | via longitudinal IPW composition | `R/ipw.R` |
| Multivariate treatment (gcomp) | inherits | inherits from ICE multivariate | `R/gcomp.R` |
| Static / shift / scale_by / threshold / dynamic | as per causatr rules | as per ICE rules | `R/interventions.R` |
| IPSI (`ipsi()`) | IPW only, inherits | via longitudinal IPW composition | `R/ipw_weights.R` |
| Stochastic interventions | via composition | via composition | `PHASE_12_STOCHASTIC.md` in causatr |
| ATE / ATT / ATC | inherits (ATT/ATC: gcomp any, IPW static binary only) | ICE: ATE only | existing estimand gating |
| `by`-stratification | inherits | inherits | `R/contrast.R` |
| Effect modification (`A:modifier`) | inherits | inherits (lag-expanded) | `R/effect_modification.R` |
| External weights (survey / IPCW) | inherits | inherits | `check_weights()` |
| Bootstrap variance | inherits (resample individuals) | inherits | `R/variance_bootstrap.R` |
| Numeric Tier 1/2 fallback | inherits | inherits | `variance_if_numeric()` |
| `censoring =` row filter | inherits | inherits | `get_fit_rows()` |

**Survival-specific work** — the three items that do NOT compose away
and must be written in the new package:

1. **Pooled-logistic outcome model** on person-period data (fit and
   contrast paths).
2. **Survival-curve estimand shape** — `contrast()` returns a
   time-indexed `data.table`, not a scalar.
3. **Cross-time variance aggregation** — delta-method chain on the
   cumulative product $\hat{S}^a_i(t) = \prod_{k \leq t}(1 -
   \hat{h}^a_{i,k})$; per-individual IF becomes a
   $n \times |t\text{-grid}|$ matrix. Pointwise bands only (no
   simultaneous bands in v1).

Everything else is wiring, type gating, and tests.

## 5. What causatr provides that the new package should call into

causatr is the **engine**. The new package should depend on causatr
(`Imports:` line in DESCRIPTION) and reuse:

### Public causatr API

- **`causat(data, outcome, treatment, confounders, estimator, ...)`**
  — fits a gcomp / IPW / matching model on scalar outcomes. Use this
  when the outcome model is a **standard** GLM (e.g., for the treatment
  side of an IPW pipeline where the "outcome" is the treatment). For
  hazard models, the new package fits its own pooled-logistic model
  but should match the `model_fn` convention (see below).
- **`contrast(fit, interventions, ...)`** — scalar contrast path. Not
  used directly for survival curves, but is the template for the
  survival-curve contrast path.
- **`diagnose(fit)`** — positivity, balance, weight summaries. Extend
  with survival-aware panels (per-period hazard positivity, per-period
  balance at each time step, competing-risks decomposition).
- **Interventions**: `static()`, `shift()`, `scale_by()`,
  `threshold()`, `dynamic()`, `ipsi()`. Fully reusable on the
  treatment column(s) of person-period data.
- **`to_person_period()`** — wide→long reshape. Remains in causatr.
- **`nhefs`** dataset — remains in causatr. The new package can
  augment with a `death` + `yrdth` column subset of the full NHEFS
  extract; the stripped-down bundled version in causatr does **not**
  have those columns.

### Design conventions worth preserving in the new package

- **`model_fn` parameter** — user passes the fitting function
  (`stats::glm`, `mgcv::gam`, ...) instead of hardcoding. Survival's
  hazard model should accept this same contract (default
  `stats::glm` with `family = binomial()` / `quasibinomial()` for
  weighted fits).
- **`estimator`, not `method`** — never use `method =` as a top-level
  argument because `MatchIt::matchit()` and other delegated functions
  use `method =` themselves.
- **Two-step API** — fit once, contrast many interventions. Keep this
  for survival (fit the hazard model once, then compute curves under
  each intervention).
- **Single influence-function engine** — one dispatcher, per-method
  branches, share primitives (`prepare_model_if()`,
  `apply_model_correction()`, `vcov_from_if()`). The survival variance
  engine should follow the same pattern, with the cross-time delta
  chain as a new aggregation layer.
- **`rlang::abort()` / `rlang::warn()` / `rlang::inform()`** — no
  `stop()` / `warning()` / `message()`.
- **`data.table` internally**, return `data.table` from user-facing
  functions.
- **Roxygen on every function, `@noRd` for internals.**
- **Reserved-column guard** — causatr's `.causatr_prev_event` /
  `.causatr_prev_cens` / `.pseudo_y` convention. The survival package
  should add its own (e.g. `.cs_prev_event`, `.cs_prev_cens`) and
  guard against user-data collisions upfront.
- **`rlang::check_weights()`-style validation** — external weights are
  validated at the boundary (NA, Inf, NaN, negative, non-numeric,
  mis-sized all rejected). Zero weights allowed as row exclusion.
- **Reject `na.action = na.exclude`** — under `na.exclude`,
  `residuals(model, "working")` is padded with NAs while
  `model.matrix()` drops them; the IF engine's
  `d_fit * r_score` then silently recycles. causatr's
  `check_dots_na_action()` is the template.

### Internals that can be copied / adapted verbatim

The following causatr internals are the **direct ancestors** of what
the survival package needs. The handoff includes a note at each site
on whether to copy, wrap, or re-derive:

- `R/causat_survival.R` (in the pre-removal baseline) — pooled
  logistic fit + risk-set construction + weighted fit + T6 family
  choice + T9 internal column strip. **Copy verbatim** into the new
  package as the starting point for Track A's fit path.
- `R/utils.R :: is_uncensored()` — `NA | 0 ⇒ TRUE` convention for
  censoring indicators. **Copy.**
- `R/checks.R :: check_weights()` — external weights validation.
  **Copy.**
- `R/checks.R :: check_dots_na_action()` — rejection path for
  `na.exclude`. **Copy.**
- `R/utils.R :: check_reserved_cols()` pattern — column-name guard.
  **Copy, adapt reserved list.**
- `R/variance_if.R :: prepare_model_if()` / `apply_model_correction()`
  / `vcov_from_if()` — the IF primitives. **Do not copy.** Call them
  via `causatr:::` (or re-export from causatr as `@export` internals)
  and layer the cross-time delta aggregation on top.
- `R/variance_bootstrap.R` — individual-level bootstrap. **Pattern
  transferable**; the survival package should resample individuals
  (all person-period rows together) and refit the hazard model per
  replicate.

## 6. The pre-removal Phase 7 design in full

The rest of this doc is the Phase 7 design that was in causatr,
extracted here so the new package has the full plan without needing
to look back at causatr's git history. Chunks 7a–7h map to
recommended implementation milestones; section numbering matches the
original design doc.

### Track A — Point survival via pooled logistic

#### Outcome model

On person-period data:
$$
\text{logit}\, h(t \mid A, L) = \alpha(t) + \beta_A A + \beta_L L
$$

with $\alpha(t)$ a flexible function of time (dummies, splines, or
`s(t)` via `mgcv::gam`). The observed response is the event indicator
$Y_{i,k} = \mathbb{1}\{T_i = k, C_i \geq k\}$ at each row; censored
rows contribute $Y_{i,k} = 0$ up to their last-observed period.

#### Intervention and standardization

`apply_intervention()` acts on the treatment column(s) of the
person-period data (broadcasting baseline `A` to every period for
point treatments). The contrast path predicts $\hat{h}^a_{i,k}$ per
person-period row, cumulates within individual to
$\hat{S}^a_i(t) = \prod_{k \leq t}(1 - \hat{h}^a_{i,k})$, and averages
across individuals to $\hat{S}^a(t) = (1/n) \sum_i \hat{S}^a_i(t)$.

#### Estimands available at contrast time

- **Survival at $t$:** $\hat{S}^a(t)$
- **Risk at $t$:** $1 - \hat{S}^a(t)$
- **Risk difference at $t$:**
  $(1 - \hat{S}^{a_1}(t)) - (1 - \hat{S}^{a_0}(t))$
- **Risk ratio at $t$:**
  $(1 - \hat{S}^{a_1}(t)) / (1 - \hat{S}^{a_0}(t))$
- **RMST difference:**
  $\int_0^{t^*} [\hat{S}^{a_1}(u) - \hat{S}^{a_0}(u)]\, du$
  (trapezoidal over the time grid)
- **Hazard ratio (averaged)** — optional; the book discourages it.

The contrast path returns a `data.table` indexed by `time`.

#### Variance

Sandwich variance for $\hat{S}^a(t)$ propagates through the cumulative
product via the delta method:
$$
d\hat{S}^a_i(t) = -\hat{S}^a_i(t) \sum_{k \leq t}
  \frac{d\hat{h}^a_{i,k}}{1 - \hat{h}^a_{i,k}}
$$
The outcome-model IF (via causatr's `prepare_model_if()` /
`apply_model_correction()`) produces per-row IFs on $\hat{h}^a_{i,k}$.
Aggregate across rows within individual via the delta chain to get
per-individual IFs on $\hat{S}^a_i(t)$; average to an IF on
$\hat{S}^a(t)$. Cross-time covariance is the cross-product of the
per-individual IF vectors evaluated at each $t$.

Bootstrap resamples individuals and refits.

### Track B — Longitudinal survival via ICE hazards

Zivich et al. (2024) iterated conditional expectations extended to the
hazard link:

1. At the final time $K$: fit
   $\hat{h}_K(a, \bar{l}_K, \bar{a}_{K-1}) = P(Y_K = 1 \mid \text{at risk}, A_K = a, \bar{L}_K, \bar{A}_{K-1})$
   as a logistic model on rows at risk at $K$.
2. Step backward to $k < K$: form the pseudo-outcome as the predicted
   conditional survival tail through the remaining periods under the
   intervention, $\tilde{Y}_k = 1 - \hat{S}^d_{k:K}$, and fit
   $\hat{h}_k$ targeting this pseudo-outcome (quasibinomial).
3. Iterate to $k = 0$.
4. Counterfactual risk:
   $\hat{R}^d(t) = (1/n) \sum_i [1 - \hat{S}^d_{0:t, i}]$.

This is exactly causatr's Phase 5 ICE with the **hazard indicator /
survival tail** replacing the scalar outcome. The forward sensitivity
recursion in `variance_if_ice()` is **agnostic to what the per-step
response is**, as long as the per-step model exposes `family$mu.eta`
and `family$variance` — both `binomial` and `quasibinomial` do. So the
ICE variance engine reuses directly; survival-Track-B only adds the
cumulative-product survival-curve shape on top.

#### What survives from causatr's ICE as-is

- `ice_iterate()` structure (backward loop).
- `variance_if_ice_one()` (forward sensitivity; block-triangular
  bread).
- Treatment-lag EM auto-expansion (`expand_em_lag_terms()`).
- External weights propagation.
- Bootstrap at individual level.

#### What needs a thin survival-aware wrapper

- Per-step target construction (hazard indicator at $K$,
  survival-tail pseudo-outcome at $k < K$).
- Per-step link forcing (binomial at $K$, quasibinomial at $k < K$).
- Cumulative-product aggregation inside the contrast path.

### Competing risks — cause-specific hazards + CIF

For $J$ competing event types, fit $J$ parallel cause-specific hazard
models $h^{(j)}(t \mid A, L)$ on person-period data, each treating
rows with event type $j' \neq j$ as censored at their event time. The
cause-specific cumulative incidence under intervention $a$ is
$$
F^{(j),a}(t) = \sum_{k=1}^{t} \hat{S}^a(k-1) \cdot \hat{h}^{(j), a}(k)
$$
where $\hat{S}^a$ uses the all-cause hazard $\sum_j \hat{h}^{(j), a}$.
This is the cause-specific decomposition from Hernán & Robins Ch. 17.
Subdistribution-hazard (Fine–Gray) models are **not** targeted.
Document the choice explicitly.

### Survival-curve estimand shape

The contrast API returns:

```r
result$estimates     # intervention | time | s_hat | se_s | ...
result$contrasts     # contrast | time | estimate | se | ...
result$time_grid     # numeric vector of time points
```

S3 methods (`print`, `plot`, `tidy`) dispatch on the survival shape to
render curves, melt over time, etc. `forrest`-style forest plots are
still available at a user-chosen reference time $t^*$. **This is the
largest user-visible API change vs. scalar causatr.**

### Cross-time variance aggregation

Sandwich variance for a survival curve is a
$|t\text{-grid}| \times |t\text{-grid}|$ covariance. Per-individual IF
is stored as an $n \times |t\text{-grid}|$ matrix; aggregate to the
full time-covariance via `crossprod(IF_mat) / n^2`. Pointwise SEs are
the diagonal. RMST SEs use quadratic forms $a^\top V a$ with $a$ the
trapezoidal weights. Simultaneous bands (Hall–Wellner, etc.) are
**out of scope for v1**.

## 7. Cross-phase composition (planned)

The new package should hold its own design notes on composition with
causatr's pending phases. Short summary:

- **Multivariate IPW + point survival.** Joint density
  $f(A_1, A_2 \mid L)$ fit once on the **original-row** data (baseline
  properties), product density-ratio weight **broadcast** onto every
  person-period row. Hazard MSM absorbs the weight per row.
- **Longitudinal IPW + survival.** Per-period treatment density models
  produce cumulative density-ratio weights indexed by $(i, k)$.
  Weighted pooled-logistic hazard MSM on uncensored person-period
  rows; survival curve from cumulative product. This is the IPW
  analogue of ICE-survival.
- **Stochastic + survival.** MC draws are taken **at the individual
  cumulative-product level**, not at the hazard level:
  $\hat{S}^g_i(t) = (1/M) \sum_m \prod_{k \leq t} (1 - \hat{h}(k \mid A_{i,m}, L_{i,k}))$.
  Averaging hazards before cumulating is biased (Jensen's inequality:
  the cumulative product is nonlinear).
- **Built-in IPCW + survival.** Survival is the *motivating* use case
  for IPCW, not a supported afterthought. Per-period cumulative IPCW
  weights multiply into the hazard MSM; stacked EE extends with
  censoring model blocks.
- **AIPW + survival.** Well-studied (Bai et al. 2013; Zhang &
  Schaubel 2012). Composes cleanly with Track A:
  $\hat{S}^g_{\mathrm{AIPW}}(t) = (1/n) \sum_i [\hat{S}^g_i(t) + W_i(g) \cdot (\mathbb{1}\{T_i > t, C_i > t\} - \hat{S}^{A_i}_i(t))]$.
- **Transportability + survival.** Sampling weight broadcast onto
  person-period rows; weighted hazard MSM on study-sample
  person-period rows; survival curve via cumulative product.
  Cross-check against `transport` / `transported` point-outcome, no
  direct lmtp equivalent for survival transport.

## 8. Oracle / reference implementations to triangulate against

- **`lmtp::lmtp_tmle(outcome_type = "survival", cens = <col>, trt = <list>, shift = ...)`**
  — direct oracle for Tracks A and B under static, shift, IPSI, and
  stochastic interventions. EIF-based SE is not directly comparable to
  the M-estimation sandwich, so it is a **point-estimate oracle** for
  variance.
- **`gfoRmula::gformula_survival()`** — forward-simulation
  g-formula reference for Track B.
- **`survival::survfit()`** — unadjusted Kaplan–Meier for sanity
  checks.
- **`survival::coxph()`** — continuous-time Cox partial likelihood.
  Pooled logistic with small per-interval hazard (< 0.1) is a close
  approximation (Hernán & Robins Technical Point 17.1; D'Agostino et
  al. 1990).
- **`causatr::causat()`** with multivariate gcomp on the person-period
  + pooled logistic route — multivariate survival reference.

## 9. NHEFS Ch. 17 replication targets

- 120-month survival: ≈ 80.7% under treatment, ≈ 80.5% under no
  treatment.
- Risk difference: ≈ 0.2% (95% CI: −4.1% to 3.7%) — essentially null.

These are the acceptance targets for the Track A truth-based test on
NHEFS.

## 10. Implementation chunks (proposed)

| Chunk | Scope | Depends on |
|---|---|---|
| 1 | Package skeleton: DESCRIPTION, NAMESPACE, lint, CI. Copy + adapt `causat_survival()` fit path. Copy `is_uncensored()`, `check_weights()`, `check_dots_na_action()`, reserved-col guard. | — |
| 2 | Track A contrast path: per-individual hazards → survival curve → risk/RMST contrasts, **no variance yet**. Time-indexed `data.table` result shape. | chunk 1 |
| 3 | Track A sandwich variance: delta-method cross-time IF aggregation. Depends on `causatr::prepare_model_if()` / `apply_model_correction()` — import or re-export as `@keywords internal`. | chunk 2 |
| 4 | Track A bootstrap + S3 methods (`print` / `plot` / `tidy` for survival curves). | chunk 2 |
| 5 | Track A under IPW: baseline density-ratio weights from `causatr::fit_ipw()`-style treatment model, **broadcast** onto person-period rows, weighted hazard MSM. | chunk 2, causatr IPW |
| 6 | Track B (ICE-hazards): per-step hazard target + survival-tail pseudo-outcome, **reuse** causatr's `ice_iterate()` and `variance_if_ice()` via internal imports. | chunk 3, causatr ICE |
| 7 | Competing risks: parallel cause-specific hazards + CIF contrast + sandwich via stacked EE across cause-specific models. | chunks 2, 3 |
| 8 | Matching rejection path + classed error. | — |
| 9 | NHEFS Ch. 17 replication test + `survival` vignette. | chunks 2–7 |
| 10 | Survival-aware `diagnose()` (per-period hazard positivity, cross-time balance, competing-risks decomposition). | chunks 2, 7 |

## 11. Package naming / placement

- **Repo:** already created under `etverse/` (user-confirmed — name TBD
  by the user on first commit).
- **Suggested package name:** something short + mnemonic in the
  etverse style (e.g. `survcausatr`, `causurvtr`, `causatr.surv`).
  Not prescribed here.
- **Dependencies (`Imports`):** `causatr`, `data.table`, `rlang`,
  `stats`, `sandwich`, `numDeriv`, `boot`.
- **Dependencies (`Suggests`):** `survival`, `lmtp`, `gfoRmula`,
  `MatchIt` (for the rejection path test), `mgcv`, `splines`,
  `testthat`, `quarto` (vignettes).

## 12. What was removed from causatr (for reference)

The following files and symbols were removed from causatr in the
same commit that introduced this handoff:

- `R/causat_survival.R` (entire file)
- `CAUSATR_SURVIVAL_INTERNAL_COLS` constant in `R/utils.R`
- `.causatr_prev_event` / `.causatr_prev_cens` from
  `CAUSATR_RESERVED_COLS`
- Survival branch in `R/contrast.R::compute_contrast()`
- Export of `causat_survival` from `NAMESPACE`
- `man/causat_survival.Rd`
- `vignettes/survival.qmd` and built artifacts
- `PHASE_7_SURVIVAL.md`
- Survival row / sections in `FEATURE_COVERAGE_MATRIX.md`
- Phase 7 line in `CLAUDE.md` and `DESCRIPTION`
- Survival-composition subsections in `PHASE_8`, `PHASE_10`,
  `PHASE_11`, `PHASE_12`, `PHASE_14`, `PHASE_16`, `PHASE_17` design
  docs
- `survival` from `Suggests` in `DESCRIPTION` (kept only if another
  non-survival use remains; audited at removal time)
- All survival-focused tests (`test-s3-methods.R` survival blocks,
  `test-simulation.R` survival blocks, `test-causat.R` competing
  block, `test-weights-edge-cases.R` survival block,
  `test-critical-review-2026-04.R` B5 block) and corresponding
  snapshots

`to_person_period()` and `is_uncensored()` **stay in causatr** because
they serve non-survival longitudinal use cases as well.

## 13. References

- Hernán MA, Robins JM (2025). *Causal Inference: What If*. Chapman &
  Hall/CRC. **Chapter 17** (survival analysis, IP weighting and
  standardization).
- Zivich PN, Ross RK, Shook-Sa BE, Cole SR, Edwards JK (2024).
  Empirical sandwich variance estimator for iterated conditional
  expectation g-computation. *Stat Med* 43:5562–5572.
- Young JG, Tchetgen Tchetgen EJ (2014). Simulation from a known
  cause-specific cumulative incidence function. *Stat Med*
  33:1098–1114.
- D'Agostino RB, Lee ML, Belanger AJ, Cupples LA, Anderson K,
  Kannel WB (1990). Relation of pooled logistic regression to time
  dependent Cox regression analysis: the Framingham Heart Study.
  *Statistics in Medicine* 9:1501–1515.
- Cole SR, Hernán MA (2004). Adjusted survival curves with inverse
  probability weights. *Computer Methods and Programs in Biomedicine*
  75:45–49.
- Fine JP, Gray RJ (1999). A proportional hazards model for the
  subdistribution of a competing risk. *JASA* 94:496–509.
  *(out of scope for v1)*
- Bai X, Tsiatis AA, O'Brien SM (2013). Doubly robust estimators of
  treatment-specific survival distributions. *Biometrics* 69:830–839.
- Zhang M, Schaubel DE (2012). Contrasting treatment-specific survival
  using double-robust estimators. *Stat Med* 31:4255–4268.
