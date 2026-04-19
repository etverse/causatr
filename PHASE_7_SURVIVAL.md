# Phase 7 — Causal Survival Analysis

> **Status: SCAFFOLDED** (`causat_survival()` fits pooled logistic; `to_person_period()` done; survival-curve contrast path pending)
>
> **Book chapter:** 17
>
> **Depends on:** Phase 2 (point gcomp), Phase 4 (self-contained IPW), Phase 5 (ICE infrastructure)
>
> **Composes with (planned):** Phase 8 (multivariate IPW), Phase 10 (longitudinal IPW), Phase 12 (stochastic), Phase 14 (IPCW)

## Scope

Two tracks, one cross-cutting concern:

1. **Track A — Point survival:** baseline treatment, time-to-event outcome, pooled logistic hazard model on person-period data (Hernán & Robins Ch. 17).
2. **Track B — Longitudinal survival:** time-varying treatment, time-to-event outcome, ICE with iterated conditional hazards (Zivich et al. 2024 extended to the hazard link).
3. **Competing risks** (cause-specific hazards + cumulative incidence) across both tracks.

Matching + survival is explicitly out of scope (§ "Matching + survival").

## Guiding principle: inheritance, not reinvention

Survival in causatr is **an outcome-model structure, not a new dimension of the design matrix**. Switching to pooled logistic on person-period data is a change on the outcome side; it is orthogonal to how the treatment is modelled, how weights are computed, how interventions are applied, which estimand is targeted, and how effect modification is parsed. The design consequence:

- **Track A reduces to:** point gcomp / IPW with the outcome model swapped to `glm(event ~ A + L + f(t), family = binomial, data = person_period)` and the contrast swapped to `S^a(t) = (1/n) Σᵢ ∏ₖ₌₁ᵗ (1 − ĥ^a_{i,k})`.
- **Track B reduces to:** Phase 5 ICE with the per-step outcome being the hazard indicator and the per-step link being logistic. `R/ice.R`, `ice_iterate()`, and `variance_if_ice()` generalize — they do not duplicate.

The corollary is that **almost every combination supported by the gcomp/IPW/ICE tables inherits into survival for free**. Phase 7's real work is confined to three genuinely survival-specific items (§ "Survival-specific design"); everything else is wiring, type gating, and tests.

## Track A — Point survival via pooled logistic

### Outcome model

On person-period data from `to_person_period()`:
$$
\text{logit}\, h(t \mid A, L) = \alpha(t) + \beta_A A + \beta_L L
$$
with $\alpha(t)$ a flexible function of time (dummies, splines, or `s(t)` via `mgcv::gam`). The observed response is the event indicator $Y_{i,k} = \mathbb{1}\{T_i = k, C_i \geq k\}$ at each row; censored rows contribute $Y_{i,k} = 0$ up to their last-observed period.

### Intervention and standardization

`apply_intervention()` acts on the treatment column(s) of the person-period data (broadcasting baseline `A` to every period for point treatments). `contrast()` then predicts $\hat{h}^a_{i,k}$ per person-period row, cumulates within individual to $\hat{S}^a_i(t) = \prod_{k \leq t}(1 - \hat{h}^a_{i,k})$, and averages across individuals to $\hat{S}^a(t) = (1/n) \sum_i \hat{S}^a_i(t)$. This is the standardized survival curve under intervention $a$.

### Estimands available at contrast time

- **Survival at $t$:** $\hat{S}^a(t)$
- **Risk at $t$:** $1 - \hat{S}^a(t)$
- **Risk difference at $t$:** $(1 - \hat{S}^{a_1}(t)) - (1 - \hat{S}^{a_0}(t))$
- **Risk ratio at $t$:** $(1 - \hat{S}^{a_1}(t)) / (1 - \hat{S}^{a_0}(t))$
- **Hazard ratio (averaged):** ratio of pooled hazards under intervention (requested only if interpretable; book discourages it)
- **RMST difference:** $\int_0^{t^*} [\hat{S}^{a_1}(u) - \hat{S}^{a_0}(u)]\, du$ (trapezoidal over the time grid)

The contrast path returns a `data.table` of estimates indexed by `time`, not a scalar. This is the biggest structural change to `contrast()`.

### Variance

Sandwich variance for $\hat{S}^a(t)$ propagates through the cumulative product via the delta method: $d\hat{S}^a_i(t) = -\hat{S}^a_i(t) \sum_k d\hat{h}^a_{i,k} / (1 - \hat{h}^a_{i,k})$. The outcome-model IF machinery (`prepare_model_if` / `apply_model_correction`) produces per-row IFs on $\hat{h}^a_{i,k}$; these are aggregated across rows within individual via the delta-method chain to yield per-individual IFs on $\hat{S}^a_i(t)$, then averaged to an IF on $\hat{S}^a(t)$. Covariance across the time grid is the cross-product of the per-individual IF vectors evaluated at each $t$.

Bootstrap resamples individuals (all their person-period rows together), refits, and recomputes the full survival curve per replicate — exactly like ICE.

## Track B — Longitudinal survival via ICE hazards

### Iterated conditional hazards (Zivich et al. 2024)

For time-varying treatment and time-to-event outcome under ICE:

1. At the final time $K$: fit $\hat{h}_K(a, \bar{l}_K, \bar{a}_{K-1}) = P(Y_K = 1 \mid \text{at risk}, A_K = a, \bar{L}_K, \bar{A}_{K-1})$ as a logistic model on rows at risk at $K$.
2. Step backward to $k < K$: form the pseudo-outcome as the predicted **conditional survival tail** through the remaining periods under the intervention, $\tilde{Y}_k = 1 - \hat{S}^d_{k:K}$, and fit $\hat{h}_k$ targeting this pseudo-outcome (fractional-logistic / quasibinomial when $\tilde{Y}_k \in [0, 1]$).
3. Iterate to $k = 0$.
4. Counterfactual risk: $\hat{R}^d(t) = (1/n) \sum_i [1 - \hat{S}^d_{0:t, i}]$.

This is exactly Phase 5 ICE with the **hazard indicator / survival tail** replacing the scalar outcome. The forward sensitivity recursion in `variance_if_ice()` carries over unchanged — it is agnostic to what the per-step response is, as long as the per-step model exposes `family$mu.eta` and `family$variance` (both `binomial` and `quasibinomial` do).

### What survives from Phase 5 as-is

- `ice_iterate()` structure (backward loop over per-step models).
- `variance_if_ice_one()` (forward sensitivity recursion; same block-triangular bread).
- Treatment-lag auto-expansion (effect modification; Chunk 6d).
- External weights propagation (including manual IPCW today).
- Bootstrap at the individual level.

### What needs a thin survival-aware wrapper

- Per-step target construction (hazard indicator at $K$, survival-tail pseudo-outcome at $k < K$).
- Per-step link forcing (binomial at $K$, quasibinomial at $k < K$).
- Cumulative-product aggregation inside the contrast path (same delta-method chain as Track A, but operating on the ICE per-step hazards rather than pooled-logistic predictions).

## Survival-specific design

These are the items that are **genuinely survival** and do not compose away:

### 1. Competing risks — cause-specific hazards + CIF

For $J$ competing event types, fit $J$ parallel cause-specific hazard models $h^{(j)}(t \mid A, L)$ on person-period data, each treating rows with event type $j'\neq j$ as censored at their event time. The cause-specific cumulative incidence under intervention $a$ is
$$
F^{(j),a}(t) = \sum_{k=1}^{t} \hat{S}^a(k-1) \cdot \hat{h}^{(j), a}(k)
$$
where $\hat{S}^a$ uses the all-cause hazard $\sum_j \hat{h}^{(j), a}$. This is the cause-specific decomposition from Hernán & Robins Ch. 17. Subdistribution-hazard (Fine-Gray) models are **not** targeted in Phase 7 — they require a different data structure and produce a different estimand. Document the choice explicitly.

### 2. Survival-curve estimand shape

`contrast()` currently returns one row per intervention per contrast (scalar estimand). Survival returns one row per intervention per contrast **per time point**. The API shape change:

```r
# Current (scalar)
result$estimates     # intervention | mu_hat | se | ...
result$contrasts     # contrast | estimate | se | ...

# Survival
result$estimates     # intervention | time | s_hat | se_s | ...
result$contrasts     # contrast | time | estimate | se | ...
result$time_grid     # numeric vector of time points
```

S3 methods (`print`, `plot`, `tidy`) dispatch on the `fit_type == "survival"` attribute to render curves (plot) or melt over time (tidy). `forrest`-style forest plots are still available at a user-chosen reference time $t^*$. This is the largest user-visible API change in Phase 7.

### 3. Cross-time variance aggregation

Sandwich variance for a survival curve is a $|t\text{-grid}| \times |t\text{-grid}|$ covariance, not a scalar. Internally, the per-individual IF is stored as an `n × |t-grid|` matrix, and `vcov_from_if()` aggregates to the full time-covariance via `crossprod(IF_mat) / n^2`. Pointwise SEs are the diagonal; confidence bands use the full covariance. RMST SEs use quadratic forms $a^\top V a$ with $a = \Delta t$ trapezoid weights. Confidence bands (simultaneous, Hall-Wellner, etc.) are out of scope for Phase 7 — pointwise bands only.

## Inheritance from other phases

| Feature axis                  | Point survival (Track A) | Longitudinal survival (Track B) | Source                     |
|-------------------------------|--------------------------|---------------------------------|----------------------------|
| Binary / continuous / categorical / count treatment (gcomp) | ✓ inherit from point gcomp | ✓ inherit from Phase 5 ICE  | `R/gcomp.R`, `R/ice.R`     |
| Binary / continuous / categorical / count treatment (IPW)   | ✓ inherit from Phase 4     | via Phase 10 composition         | `R/ipw.R`                  |
| Multivariate treatment (gcomp) | ✓ inherit directly         | ✓ inherit from ICE multivariate  | `R/gcomp.R`                |
| Multivariate treatment (IPW)   | via Phase 8 composition    | via Phases 8+10 composition      | `PHASE_8_MULTIVARIATE_IPW.md` § Survival |
| Static / shift / scale_by / threshold / dynamic interventions | ✓ as per gcomp rules | ✓ as per ICE rules              | `R/interventions.R`        |
| `ipsi()` intervention          | ✓ IPW only, inherit        | via Phase 10 composition         | `R/ipw_weights.R`          |
| Stochastic interventions       | via Phase 12 composition   | via Phase 12 composition         | `PHASE_12_STOCHASTIC.md` § Survival |
| ATE / ATT / ATC                | ✓ inherit (ATT/ATC: gcomp any, IPW static binary only) | ICE: ATE only (unchanged)    | existing estimand gating   |
| `by`-stratification            | ✓ inherit                 | ✓ inherit                        | `R/contrast.R`             |
| Effect modification (`A:modifier`) | ✓ inherit from gcomp/IPW  | ✓ inherit from ICE lag-expansion | `R/effect_modification.R`  |
| External weights               | ✓ inherit                 | ✓ inherit                        | `check_weights()`          |
| Sandwich variance              | ✓ extended with cross-time delta | ✓ extended with cross-time delta | § "Cross-time variance"  |
| Bootstrap variance             | ✓ inherit (resample individuals) | ✓ inherit (resample individuals) | `R/variance_bootstrap.R` |
| Numeric Tier 1/2 fallback       | ✓ inherit                 | ✓ inherit                        | `variance_if_numeric()`    |
| `censoring =` row filter        | ✓ already supported        | ✓ already supported              | `get_fit_rows()`           |
| Built-in IPCW                   | via Phase 14 composition   | via Phase 14 composition         | `PHASE_14_IPCW.md` § Survival |

## Composition with pending phases

Each pending phase that interacts with survival carries a dedicated "Survival composition" subsection. The design patterns are **decided now**, so the phase's eventual implementation cannot silently skip survival. Summary:

- **Phase 8 (multivariate IPW + survival).** Joint density $f(A_1, A_2 \mid L)$ is fit on the **original-row** data (one row per individual). The product density-ratio weight is computed once per individual and **broadcast** onto every person-period row for that individual. The hazard MSM absorbs the weight per row. See `PHASE_8_MULTIVARIATE_IPW.md` § "Survival composition".

- **Phase 10 (longitudinal IPW + survival).** Per-period treatment density models produce cumulative density-ratio weights indexed by $(i, k)$. The weighted pooled-logistic hazard MSM is fit on the person-period data; the survival curve is recovered by the standard cumulative product. This is the IPW analogue of ICE-survival. See `PHASE_10_LONGITUDINAL_IPW.md` § "Survival composition".

- **Phase 12 (stochastic + survival).** MC draws are taken **at the individual cumulative-product level**, not at the hazard level: $\hat{S}^g_i(t) = (1/M) \sum_m \prod_{k \leq t} (1 - \hat{h}(k \mid A_{i,m}, L_{i,k}))$. Averaging hazards before cumulating is incorrect (nonlinearity of the cumulative product). See `PHASE_12_STOCHASTIC.md` § "Survival composition".

- **Phase 14 (IPCW + survival).** Survival is the **motivating use case**, not a supported afterthought. Censoring is the rule, not the exception, in time-to-event data. The Phase 14 design treats survival + IPCW as the primary test target, with point non-survival outcomes as the secondary case. See `PHASE_14_IPCW.md` — survival promoted to §1.

## Matching + survival: out of scope

Matching on survival outcomes is genuinely awkward: `MatchIt` produces match weights on original rows, but survival estimation requires person-period rows with subclass identifiers that survive the expansion. The standard practice is `survival::coxph(..., weights = match_weights, cluster = subclass)` — which is (a) a fundamentally different estimator (Cox partial likelihood, not pooled logistic), (b) a different variance story (robust sandwich on the partial-likelihood score, not IF on the pooled-logistic hazards), and (c) duplicative with what `survival::coxph` already does well outside causatr.

Decision: `fit_matching()` + `type = "survival"` will hard-abort with a clear error pointing to `estimator = "gcomp"` / `"ipw"`, exactly like matching today rejects non-binary treatment and longitudinal data. No Phase 7 effort is spent on matching.

## What exists today

- `causat_survival(data, outcome, treatment, time, event, censoring, competing, ...)` fits a pooled logistic hazard on person-period data. Stores `fit_type = "survival"`.
- `to_person_period(data, time, event, censoring)` converts wide → long.
- `contrast()` on a survival fit currently aborts with "survival contrast not yet implemented" (tested ✅).
- `causat_survival()` rejects `competing != NULL` at fit time (tested ✅).
- Vignette scaffold: `vignettes/survival.qmd` exists with basic structure.

## What remains

### Track A — Point survival

- [ ] Survival-aware contrast pathway in `compute_contrast()`: per-individual hazards → per-individual survival → standardized survival.
- [ ] Time-indexed estimates/contrasts data.tables.
- [ ] Risk difference / risk ratio / RMST contrasts at user-supplied times or full grid.
- [ ] Sandwich variance with cross-time covariance via delta-method chain on the cumulative product.
- [ ] Bootstrap individuals with cumulative-product recomputation per replicate.
- [ ] S3 methods (`print.causatr_result`, `plot.causatr_result`, `tidy.causatr_result`) dispatching on `fit_type = "survival"`.
- [ ] Integration with `estimator = "ipw"` (Track A under IPW): user supplies a baseline propensity model, density-ratio weights are computed once per individual and broadcast onto person-period rows; the hazard MSM is weighted.

### Track B — Longitudinal survival

- [ ] ICE-hazards path in `ice_iterate()`: per-step hazard indicator at $K$, per-step survival-tail pseudo-outcome at $k < K$.
- [ ] `causat_survival()` accepts `type = "longitudinal"` with time-varying treatment.
- [ ] `variance_if_ice()` generalization to cumulative-product survival via the same cross-time delta-method aggregation as Track A.
- [ ] Bootstrap inheritance (no structural change — resample individuals).

### Competing risks

- [ ] Parallel cause-specific hazard models fit inside `causat_survival(competing = "event_type")`.
- [ ] Cumulative incidence function under intervention per cause.
- [ ] Sandwich variance via stacked EE across cause-specific models (shared $S^a$ tail; one IF block per cause).
- [ ] Contrasts on CIF: difference, ratio.

### Tests

- [ ] Truth-based point survival: analytical cumulative hazard → $S(t)$ on a linear-Gaussian/exponential DGP.
- [ ] External reference: `lmtp::lmtp_tmle(outcome_type = "survival")` on point and longitudinal DGPs.
- [ ] NHEFS replication (Ch. 17): 120-month survival, risk difference.
- [ ] Competing risks: two-cause DGP with known CIFs.
- [ ] Cross-method agreement: gcomp-survival vs IPW-survival on a DGP where both are consistent.
- [ ] Matching rejection path.

### Documentation

- [ ] Complete `vignettes/survival.qmd` (NHEFS walkthrough + competing risks example).
- [ ] `FEATURE_COVERAGE_MATRIX.md` survival rows populated from the "What remains" checklist.

## NHEFS replication targets (Ch. 17)

- 120-month survival: ≈ 80.7% under treatment, ≈ 80.5% under no treatment.
- Risk difference: ≈ 0.2% (95% CI: −4.1% to 3.7%) — essentially null.

These are the acceptance targets for the Track A truth-based test on NHEFS.

## Implementation chunks

| Chunk | Scope | Depends on |
|---|---|---|
| 7a | Track A contrast path: per-individual hazards → survival curve → risk/RMST contrasts, no variance yet | scaffold |
| 7b | Track A sandwich variance: delta-method cross-time IF aggregation | 7a |
| 7c | Track A bootstrap + S3 method dispatch (`print` / `plot` / `tidy` for survival curves) | 7a |
| 7d | Track A under IPW: baseline density-ratio weights broadcast onto person-period, weighted hazard MSM | 7a |
| 7e | Track B (ICE-hazards): per-step hazard target + survival-tail pseudo-outcome, reuse `ice_iterate()` / `variance_if_ice()` | Phase 5, 7a, 7b |
| 7f | Competing risks: parallel cause-specific hazards + CIF contrast | 7a, 7b |
| 7g | Matching rejection path + error class | — |
| 7h | NHEFS Ch. 17 replication test + vignette | 7a–7f |

Cross-phase composition chunks (owned by each respective phase doc, not Phase 7):

- Phase 8 ships "Survival composition" chunk once Phase 8 core is done.
- Phase 10 ships "Survival composition" chunk once Phase 10 core is done.
- Phase 12 ships "Survival composition" chunk once Phase 12 core is done.
- Phase 14 ships survival-first from the start (survival promoted to motivating case).

## References

- Hernán MA, Robins JM (2025). *Causal Inference: What If*. Chapter 17 (survival analysis, IP weighting and standardization).
- Zivich PN, Ross RK, Shook-Sa BE, Cole SR, Edwards JK (2024). Empirical sandwich variance estimator for iterated conditional expectation g-computation. *Stat Med* 43:5562–5572.
- Young JG, Tchetgen Tchetgen EJ (2014). Simulation from a known cause-specific cumulative incidence function. *Stat Med* 33:1098–1114.
- Fine JP, Gray RJ (1999). A proportional hazards model for the subdistribution of a competing risk. *JASA* 94:496–509. *(out of scope for Phase 7)*
