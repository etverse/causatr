# Phase 17 — Transportability and Generalizability

> **Status: IN PROGRESS** — 17a–17l shipped; 17m pending.
>
> **Depends on:** Phase 2 (point gcomp), Phase 4 (self-contained IPW)
>
> **Composes with:** Phase 10 (longitudinal, chunk 17i), Phase 14 (IPCW, chunk 17k), Phase 16 (AIPW, chunk 17e ✅).

## Motivation

Every causal effect estimated in causatr today is an effect in the **study population** — the sample `data` supplied to `causat()`. In applied work this is almost never what is scientifically relevant. The question is usually "what is the effect of treatment in *my* target population (clinic, country, cohort)?" — a different population from the study. Moving the estimand from study to target is the job of **transportability and generalizability weights**, and it composes naturally with the self-contained density-ratio IPW engine because the sampling model is just another density to reweight by.

**Generalizability** (Cole & Stuart 2010; Hernán & VanderWeele 2017): the study sample is drawn *from* the target population, possibly with non-random selection. Upweight under-represented subgroups.

**Transportability** (Pearl & Bareinboim 2011; Dahabreh et al. 2020): the study sample is *external* to the target (different cohort, clinic, country). Reweight based on covariate-distribution differences.

The distinction is a matter of data structure and interpretation, not a matter of methods — both use a selection indicator S and weights proportional to $f_{\text{target}}(L) / f_{\text{study}}(L)$. Phase 17 implements them under a single API.

## Scope

1. **Target-ATE estimation** in gcomp, IPW, and (post-Phase-16) AIPW under a mixed-data input: the user supplies one `data` containing both study rows (S = 1, full A / L / Y) and target rows (S = 0 or a separate target-population sample, L only).
2. **Sampling model** `P(S = 1 \mid L)` fit via a user-supplied `sampling_model_fn` (default `stats::glm(..., family = binomial)`), with the Phase 4-style treatment-family generalisation later if continuous-selection indicators become useful (they almost never are — S is nearly always binary).
3. **Three estimators** with the same two-step API as today: gcomp transport (standardize on target covariates), IPW transport (sampling × treatment density-ratio weights), AIPW transport (compose Phase 16 + Phase 17 for double robustness across both the confounding model and the selection model).
4. **Sandwich variance** via a stacked estimating-equation system that extends the existing engines with the sampling-model block.
5. **Bootstrap variance** (natural — resample individuals, refit everything per replicate).
6. **Static binary interventions** are the primary target. MTP interventions (shift / scale\_by / ipsi) require observed A on target rows; for transportability (S = 0 has A = NA) these need MC marginalization (chunk 17l). MTP + generalizability works where all target rows have observed A.
7. **External cross-check** against `TransportHealth` R package (Dahabreh group) on shared DGPs where feasible.

## Non-scope

- **Matching + transportability.** MatchIt's match-weight path does not compose cleanly with external sampling weights (weighted matching with two weight sources is awkward and rarely published). Hard reject.
- **Effect modification by time-varying covariates under IPW transport.** Same MSM-level restriction as Phase 6 under IPW: only baseline modifiers are safe. Gcomp transport handles time-varying EM correctly because it uses the full outcome model.
- **Transportability under sampling-on-unobservables.** Phase 17 assumes exchangeability of potential outcomes across populations given L — i.e., a correctly-specified sampling model recovers the target estimand. Sensitivity analysis for violations of this assumption (Dahabreh & Hernán 2019) is separate and out of scope.
- **Target population defined by a stochastic rule** on L (e.g. "target = individuals with L > 0 weighted by f(L)"). Out of scope; requires a different identification story.
- **Data fusion beyond two populations.** Multi-site transportability (more than one study; e.g. Dahabreh et al. 2023) is a natural extension but adds complexity; Phase 17 ships the two-population case and defers multi-site.

## Design

### Terminology and data shape

Let S be a binary indicator with S = 1 for study-population rows and S = 0 for target-population rows. The mixed-data input is a single `data.table` whose rows are tagged by S; study rows carry full (A, L, Y), target rows carry L only (A and Y may be NA or absent — they are not used). This matches the `transport` package's convention and the Dahabreh 2020 paper's data layout.

Generalizability vs transportability: same data shape, same estimators, same weights. The only difference is whether the target population includes S = 1. For **generalizability**, the target is the union S = 0 ∪ S = 1 (study is a biased subsample of the target). For **transportability**, the target is S = 0 only (study is external). Phase 17 exposes both via a `target_subset` argument:

```r
# Generalizability: target = whole population
fit <- causat(data, ..., target = "S",
              sampling_model_fn = stats::glm,
              target_subset = "all")     # S = 0 ∪ S = 1

# Transportability: target = external population only
fit <- causat(data, ..., target = "S",
              sampling_model_fn = stats::glm,
              target_subset = "target")  # S = 0 only
```

Default `target_subset = "target"` (the more common case in applied work).

### Identification assumptions

For the target ATE ($\psi_{\text{target}} = E_{\text{target}}[Y^{a_1} - Y^{a_0}]$), Dahabreh et al. 2020 give four conditions:

1. **Conditional exchangeability in study:** $Y^a \perp A \mid L, S = 1$.
2. **Treatment positivity in study:** $0 < P(A = a \mid L, S = 1) < 1$ on the study support.
3. **Exchangeability over populations (transportability):** $E[Y^a \mid L, S = 1] = E[Y^a \mid L, S = 0]$ — the conditional potential-outcome regression is the same in study and target.
4. **Sampling positivity:** $P(S = 1 \mid L) > 0$ on the target support.

Assumption 3 is the transportability assumption and the identifying leap Phase 17 makes; it is testable only indirectly (via negative-control outcomes or sensitivity bounds). Assumption 4 is the analog of treatment positivity and is diagnosable from the sampling-model fit (similar to weight-distribution diagnostics).

### Three estimators

**Gcomp transport.**
$$
\hat\psi_{\text{target}}(a) = \frac{1}{n_{\text{target}}} \sum_{i: \, \text{target}} \hat{m}(a, L_i), \qquad \hat{m} \text{ fit on } S = 1.
$$
The outcome model is fit on study rows only; standardization averages over target covariates. This is the minimal-assumptions pathway — it does not use the sampling model at all. `contrast()` needs only a subset-filter change to predict on target rows.

**IPW transport.**
$$
\hat\psi_{\text{target}}(a) = \frac{\sum_{i: S = 1} w_i^S \cdot w_i^A(a) \cdot Y_i}{\sum_{i: S = 1} w_i^S \cdot w_i^A(a)},
$$
where $w_i^A(a) = \mathbb{1}\{A_i = a\} / \hat P(A = a \mid L_i, S = 1)$ is the treatment density-ratio weight (the Phase 4 weight) and $w_i^S = [1 - \hat P(S = 1 \mid L_i)] / \hat P(S = 1 \mid L_i)$ is the sampling odds weight (for `target_subset = "target"`; for `"all"` it is $1 / \hat P(S = 1 \mid L_i)$). The two weights multiply pointwise; the Hájek (self-normalised) form in the numerator and denominator is the default, matching the Phase 4 `Y ~ 1` MSM pathway. **General interventions** (shift, scale_by, binary dynamic, IPSI) slot in by replacing $w_i^A(a)$ with the corresponding density-ratio weight from `make_weight_fn()`, no other changes.

**AIPW transport.** Post-Phase-16 composition:
$$
\hat\psi_{\text{target}}(a) = \frac{1}{n_{\text{target}}} \sum_{i: \, \text{target}} \hat{m}(a, L_i) + \frac{1}{n_{\text{target}}} \sum_{i: S = 1} w_i^S \cdot w_i^A(a) \cdot (Y_i - \hat{m}(a, L_i)),
$$
where both terms use $1/n_{\text{target}}$ as denominator (Dahabreh et al. 2020, Section 4.2). The augmentation term sums over study rows only but is normalized by the target population size, matching the Hájek form used elsewhere. Consistent if **any two** of (outcome model, treatment model, sampling model) are correctly specified (triple robustness under some conventions; this is the "2-out-of-3 DR" structure of Dahabreh et al. 2020 Section 4.2).

### Sampling model

Fit on **all** rows (S = 1 and S = 0):
```
sampling_model <- sampling_model_fn(S ~ L, data = data, family = binomial)
```
`sampling_model_fn` defaults to `stats::glm` with logit link; users can pass `mgcv::gam` for spline-based selection models. Phase 17 does not attempt to generalise to continuous S (contrast to the propensity-family dispatch in Phase 4): selection indicators are nearly always binary in applied work, and a continuous S would need a separate identification story that is out of scope here.

The sampling model is stored in `fit$details$sampling_model` alongside the existing `$propensity_model` (for IPW fits) or `$model` (for gcomp fits).

### Stacked sandwich variance

The estimating-equation system grows by one block (the sampling model):

For gcomp transport:
$$
\omega = \begin{pmatrix} s_\beta^{\text{outcome}}(L, A, Y; \text{fit on } S = 1) \\ s_\gamma^{\text{sampling}}(L, S; \text{fit on all}) \\ \omega_\psi(L, S, A, Y; \beta, \gamma, \psi) \end{pmatrix}.
$$
The bread is block-triangular; the plug-in row has cross-derivatives with respect to $\beta$ (outcome) and $\gamma$ (sampling).

For IPW transport: add the propensity-model block to the above stack. Three score blocks + plug-in.

For AIPW transport: the stack adds the outcome-model block to the IPW-transport stack, giving four score blocks + plug-in. This is the "triple-robust" stacked EE.

The sandwich primitives (`prepare_model_if`, `compute_ipw_if_self_contained_one`, `vcov_from_if`) handle the block extension with the same pattern established in Phase 4 and carried into Phase 16 — no fundamentally new variance machinery.

### Diagnostics

`diagnose()` gains a new sampling-model panel, mirroring the propensity-score panel:

- Sampling-score histogram facet by S (check overlap of $\hat P(S = 1 \mid L)$ across populations).
- Extreme-sampling-weight flagging (positivity violations for assumption 4).
- Covariate balance between S = 1 and S = 0 after reweighting.

Phase 11 (diagnose rewrite) will fold this in; Phase 17 ships a minimal shim analogous to the Phase 4 `diagnose()` shim for IPW.

## Composition with other phases

Each composition is tracked as a concrete chunk below (17i–17l). Summary:

- **Phase 8 (multivariate IPW) × Phase 17** → chunk 17j.
- **Phase 10 (longitudinal IPW) × Phase 17** → chunk 17i.
- **Phase 12 (stochastic) × Phase 17.** MC integration on the outcome-model augmentation term; sampling weight is deterministic in L. Stochastic AIPW is tracked in **Phase 24** (`PHASE_24_STOCHASTIC_AIPW.md`). Also subsumes MTP + transportability (see 17l).
- **Phase 14 (IPCW) × Phase 17** → chunk 17k.
- **Phase 16 (AIPW) × Phase 17** → chunk 17e (✅ done).
- **MTP (shift/scale_by/ipsi) × transportability** → chunk 17l. For static interventions, $\hat{m}(a, L)$ on target rows (S=0) is well-defined. For MTPs where $d(A, L)$ depends on observed $A$, target rows have $A = \mathrm{NA}$. The solution (Díaz & Hejazi 2020; Hejazi et al. 2024) is to marginalize $E_{A|L,S=1}[\hat{m}(d(A,L), L) \mid L]$ over the study treatment distribution via MC integration or sequential regression. Until then, MTP + transportability is unsupported; MTP + generalizability works where all target rows have observed $A$.
- **Survival composition.** Owned by `etverse/survatr`; imports Phase 17 for the scalar sampling-weight primitive. Not a causatr chunk.

## Chunks

| Chunk | Scope | Depends on | Status |
|---|---|---|---|
| 17a | `fit_sampling_model()`: fit $P(S = 1 \mid L)$ via `sampling_model_fn`; validate S is binary and present in data; reject if S has NAs; store in `fit$details$sampling_model` | — | ✅ done |
| 17b | Gcomp transport: target-subset filter in `compute_contrast()`; outcome model fit on S = 1, standardization over target rows; sandwich with sampling-model cross-derivative | 17a, Phase 2 | ✅ done |
| 17c | IPW transport: sampling × treatment weight product in the IPW weight pipeline (`ipw_weights.R` / `make_weight_fn()`); weighted MSM on study rows; stacked sandwich | 17a, Phase 4 | ✅ done |
| 17d | Bootstrap (refit sampling + propensity + outcome per replicate) | 17a–17c | ✅ done |
| 17e | AIPW transport: compose Phase 16 + Phase 17; 2-out-of-3 DR test (deliberately misspecify any one of outcome / treatment / sampling — verify consistency) | 17a–17c, Phase 16 | ✅ done |
| 17f | `diagnose()` shim: sampling-score panel + extreme-sampling-weight flags | 17a | ✅ done |
| 17g | Sampling model predictor-set validation: emit `rlang::warn()` when confounders formula RHS is a strict subset of outcome/treatment model predictors | 17a | ✅ done |
| 17h | External cross-check against `TransportHealth` R package on shared DGPs | 17b, 17c | ✅ done |
| 17i | Longitudinal transport (Phase 10 × Phase 17): broadcast sampling weight onto person-period rows; multiply into per-period treatment weight; weighted longitudinal MSM | Phase 10, 17c | ✅ done |
| 17j | Multivariate IPW × transport (Phase 8 composition): joint treatment density × sampling density; stacked sandwich with multivariate propensity + sampling blocks | Phase 8, 17c | ✅ done |
| 17k | IPCW × transport (Phase 14 composition): four-way stacked EE (outcome/propensity + censoring + sampling + plug-in); triply-weighted estimator | Phase 14, 17c | ✅ done |
| 17l | MTP + transportability: MC marginalization $E_{A|L,S=1}[\hat{m}(d(A,L), L) \mid L]$ for shift/scale\_by/ipsi on target rows where $A = \mathrm{NA}$; sequential regression or Monte Carlo integration (Díaz & Hejazi 2020; Hejazi et al. 2024) | 17e | ✅ done |
| 17m | Documentation, vignette (`transportability.qmd`), `FEATURE_COVERAGE_MATRIX.md` rows, `CLAUDE.md` update | 17a–17l | pending |

## Invariants

- `target = "S"` requires S to be a binary 0/1 column with no NAs. The sampling model is fit on **all rows**, not just study rows.
- When `target = NULL` (default), Phase 17's pathway is inactive and `causat()` behaves exactly as pre-Phase-17 — the study estimand is returned as today. Non-breaking change.
- Under `target_subset = "target"`, target-ATE rows MUST have no treatment (A) or outcome (Y) dependency — the estimator only uses target-row L. If A or Y are present on target rows, they are ignored with a silent `rlang::inform()` note (not an error — users sometimes leave observed A in target rows for convenience).
- The sampling model's predictor set should be a superset of the outcome model's and treatment model's predictor sets — otherwise the transportability assumption may fail silently. Tracked as chunk 17g.
- Bootstrap MUST refit the sampling model per replicate (it has estimated parameters); a replicate that fails (degenerate S distribution in the bootstrap sample) returns NA and is excluded, same convention as Phase 2 bootstrap.

## DGP for truth-based tests

### Generalizability (chunk 17b–17c)

Simplified example (constant treatment effect — study and target ATEs coincide, so not useful for truth-based tests but illustrates the data structure):
```
L ~ N(0, 1)                                (population)
P(S = 1 | L) = expit(-0.5 + 1.0 · L)      (sampling: under-represents L < 0)
A | L, S = 1 ~ Bernoulli(expit(0.2 + 0.3 · L))
Y | A, L ~ N(2 + 3 · A + 1.5 · L, 1)

Target ATE = E[Y^1 - Y^0] = 3              (over full population N(0,1))
Study ATE = E[Y^1 - Y^0 | S = 1] = 3       (coincidentally 3 here because
                                            treatment effect is constant in L)
```
To make study-target divergence visible (this is the DGP used in `simulate_transport()` and truth-based tests), add an interaction:
```
Y | A, L ~ N(2 + 3 · A + 1.5 · L + 1.0 · A · L, 1)
Target ATE = 3 + 1.0 · E[L] = 3
Study ATE = 3 + 1.0 · E[L | S = 1] = 3 + 1.0 · E[L | S = 1]
```
which differs from 3 by the L-conditional expectation in the biased sample. Truth-based test: Phase 17 estimator must recover 3; unadjusted study estimator must not.

### Transportability (chunk 17b–17c)

Same as above but target is S = 0 only and the interaction makes study $\neq$ target.

## References

- Cole SR, Stuart EA (2010). Generalizing evidence from randomized clinical trials to target populations. *Am J Epidemiol* 172:107–115.
- Pearl J, Bareinboim E (2011). Transportability of causal and statistical relations. *AAAI*.
- Hernán MA, VanderWeele TJ (2017). Compound treatments and transportability of causal inference. *Epidemiology* 22:368–377.
- Dahabreh IJ, Robertson SE, Tchetgen EJ, Stuart EA, Hernán MA (2019). Generalizing causal inferences from individuals in randomized trials to all trial-eligible individuals. *Biometrics* 75:685–694.
- Dahabreh IJ, Robertson SE, Steingrimsson JA, Stuart EA, Hernán MA (2020). Extending inferences from a randomized trial to a new target population. *Stat Med* 39:1999–2014.
- Westreich D, Edwards JK, Lesko CR, Stuart EA, Cole SR (2017). Transportability of trial results using inverse odds of sampling weights. *Am J Epidemiol* 186:1010–1014.
- `TransportHealth` R package (Dahabreh group).
