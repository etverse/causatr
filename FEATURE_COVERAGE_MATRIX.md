# Feature Coverage Matrix

Single source of truth for what works, what's tested, and at what fidelity.
**Every PR that changes a feature MUST update this file.**

## Legend

| Symbol | Meaning |
|---|---|
| ✅ | Truth-based: estimate/SE/CI checked against analytical or external reference |
| 🟡 | Smoke: runs without error, finite output, target not pinned |
| ❌ | No test |
| ⛔ | Rejected by design (rejection path tested) |

References: lmtp (`lmtp_tmle`, `lmtp_sdr`), Hernán & Robins book values, closed-form analytical truth, WeightIt oracles, `delicatessen`.

---

## G-comp (point)

| Trt | Outcome | Model | Intervention | Estimand | Contrast | Variance | Wt | Status | Test |
|---|---|---|---|---|---|---|---|---|---|
| bin | gauss | GLM | static | ATE | diff | sandwich | — | ✅ | test-gcomp.R, test-simulation.R |
| bin | gauss | GLM | static | ATE | diff | boot | — | ✅ | test-gcomp.R, test-simulation.R |
| bin | gauss | GLM | static | ATT | diff | sandwich | — | ✅ | test-gcomp.R, test-simulation.R |
| bin | gauss | GLM | static | ATC | diff | sandwich | — | ✅ | test-simulation.R |
| bin | gauss | GLM | static | ATE | diff | sandwich | survey | ✅ | test-simulation.R |
| bin | gauss | GLM | static | ATE | diff | boot | survey | ✅ | test-simulation.R |
| bin | binom | GLM | static | ATE | diff/ratio/OR | sandwich | — | ✅ | test-simulation.R |
| bin | binom | GLM | static | ATE | diff | boot | — | ✅ | test-simulation.R |
| bin | poisson | GLM | static | ATE | ratio | sandwich | — | ✅ | test-gcomp.R |
| bin | quasibinom | GLM | static | ATE | diff | sandwich | — | ✅ | test-simulation.R |
| bin | gamma | GLM | static | ATE | ratio | sandwich | — | ✅ | test-simulation.R |
| bin | gauss | GAM | static | ATE | diff | sandwich | — | ✅ | test-complex-dgp.R |
| bin | gauss | GAM | static | ATE | diff | boot | — | ✅ | test-complex-dgp.R |
| bin | gauss | GLM+splines | static | ATE | diff | sandwich | — | ✅ | test-complex-dgp.R |
| bin | gauss | GLM | static | by(L) | diff | sandwich | — | ✅ | test-by-estimand.R |
| bin | gauss | GLM | static | by(L) | diff | boot | — | ✅ | test-by-estimand.R |
| cont | gauss | GLM | shift | ATE | diff | sandwich | — | ✅ | test-gcomp.R, test-simulation.R |
| cont | gauss | GLM | shift | ATE | diff | boot | — | ✅ | test-simulation.R |
| cont | gauss | GLM | scale_by | ATE | diff | sandwich | — | ✅ | test-simulation.R |
| cont | gauss | GLM | threshold | ATE | diff | sandwich | — | ✅ | test-simulation.R |
| bin | gauss | GLM | dynamic | ATE | diff | sandwich | — | ✅ | test-simulation.R |
| multi | gauss | GLM | static | ATE | diff | sandwich | — | ✅ | test-multivariate.R |
| multi | gauss | GLM | shift | ATE | diff | sandwich | — | ✅ | test-multivariate.R |
| multi | gauss | GLM | scale_by | ATE | diff | sandwich | — | ✅ | test-multivariate.R |
| multi | gauss | GLM | threshold | ATE | diff | sandwich | — | ✅ | test-multivariate.R |
| multi | gauss | GLM | dynamic | ATE | diff | sandwich | — | ✅ | test-multivariate.R |
| multi | binom | GLM | static | ATE | diff/ratio/OR | sandwich | — | ✅ | test-multivariate.R |
| multi | gauss | GLM | static | subset | diff | sandwich | — | ✅ | test-multivariate.R |
| multi | gauss | GLM | static | by(L) | diff | sandwich | — | ✅ | test-multivariate.R |
| multi | gauss | GLM | static | ATE | diff | boot | — | ✅ | test-multivariate.R |
| bin | gauss | GLM | static | ATE | diff | sandwich | cens | ✅ | test-simulation.R |
| bin | gauss | custom | static | ATE | diff | numeric T1 | — | ✅ | test-variance-if.R |
| bin | gauss | custom | static | ATE | diff | numeric T2 | — | ✅ | test-variance-if.R |
| cat | gauss | GLM | static | ATE | diff | sandwich | — | ✅ | test-simulation.R |

Rejections: invalid family string ✅, missing outcome/treatment col ✅ (test-gcomp.R, test-causat.R).

---

## IPW (point)

| Trt | Outcome | Intervention | Estimand | Contrast | Variance | Wt | Status | Test |
|---|---|---|---|---|---|---|---|---|
| bin | gauss | static | ATE | diff | sandwich | — | ✅ +oracle | test-simulation.R, test-ipw.R, test-ipw-weights.R, test-ipw-weightit-oracle.R, test-treatment-model.R |
| bin | gauss | static | ATE | diff | boot | — | ✅ | test-simulation.R |
| bin | gauss | static | ATT | diff | sandwich | — | ✅ +oracle | test-simulation.R, test-ipw-weightit-oracle.R |
| bin | gauss | static | ATC | diff | sandwich | — | ✅ +oracle | test-simulation.R, test-ipw-weightit-oracle.R |
| bin | gauss | static | ATE | diff | sandwich (GAM PS) | — | ✅ +oracle | test-ipw.R, test-ipw-weightit-oracle.R |
| bin | gauss | static | ATE | diff | sandwich | survey | ✅ | test-simulation.R |
| bin | binom | static | ATE | diff/ratio/OR | sandwich | — | ✅ | test-simulation.R |
| bin | binom | static | ATE | diff | boot | — | ✅ | test-simulation.R |
| bin | gauss | dynamic | ATE | diff | sandwich | — | ✅ | test-simulation.R, test-ipw-weights.R |
| bin | gauss | ipsi($\delta$) | ATE | diff | sandwich | — | ✅ +oracle | test-simulation.R, test-ipw-weights.R, test-ipw-lmtp-oracle.R |
| cont | gauss | shift | ATE | diff | sandwich | — | ✅ +oracle | test-simulation.R, test-ipw-weights.R, test-ipw-lmtp-oracle.R |
| cont | gauss | scale_by | ATE | diff | sandwich | — | ✅ | test-simulation.R, test-ipw-weights.R |
| cat | gauss | static | ATE | diff | sandwich | — | ✅ | test-simulation.R |
| cat | gauss | dynamic | ATE | diff | sandwich | — | 🟡 | test-simulation.R |
| cont | gauss | static | ATE | diff | sandwich | — | ✅ | test-simulation.R |
| count (pois) | gauss | shift | ATE | diff | sandwich | — | ✅ | test-ipw-count.R |
| count (negbin) | gauss | shift | ATE | diff | sandwich | — | ✅ | test-ipw-count.R |
| bin | gauss | static | ATE + EM | diff | sandwich | — | ✅ | test-effect-modification.R |
| bin | gauss | static | ATE + EM | diff | boot | — | ✅ | test-effect-modification.R |

Variance internals: self-contained IF ✅ (hand-derived cross-derivative + end-to-end stacked-sandwich; test-ipw-branch-b.R, test-ipw-cross-derivative.R, test-variance-if.R). Non-static variance regression ✅ (shift ~5-8% SE reduction, IPSI ~90% off-diagonal covariance; test-ipw-variance-regression.R). Bootstrap parity ✅ (within 30% MC tolerance; test-ipw-variance-regression.R).

Rejections (all ✅ tested):
- `static`/`threshold`/`dynamic` on continuous $\to$ test-ipw.R
- `static`/`threshold`/`dynamic`/`ipsi` on count $\to$ test-ipw-count.R
- non-integer `shift()` on count $\to$ test-ipw-count.R
- non-integer-preserving `scale_by()` on count $\to$ test-ipw-count.R
- shift/scale_by/dynamic/ipsi + ATT/ATC $\to$ test-estimand-intervention-compat.R
- multivariate $\to$ test-s3-methods.R, test-multivariate.R
- longitudinal $\to$ test-ipw.R

---

## Matching (point)

| Trt | Outcome | Estimand | Contrast | Variance | Wt | Status | Test |
|---|---|---|---|---|---|---|---|
| bin | gauss | ATT | diff | sandwich | — | ✅ | test-simulation.R |
| bin | gauss | ATT | diff | boot | — | ✅ | test-simulation.R |
| bin | gauss | ATE (full) | diff | sandwich | — | ✅ | test-simulation.R |
| bin | gauss | ATC | diff | sandwich | — | ✅ | test-simulation.R |
| bin | gauss | ATT | diff | sandwich | survey | ✅ | test-matching.R |
| bin | gauss | ATT | diff | boot | survey | ✅ | test-matching.R |
| bin | binom | ATT | diff | sandwich | — | ✅ | test-simulation.R |
| bin | binom | ATT | ratio/OR | sandwich | — | ✅ | test-simulation.R |
| bin | quasibinom | ATT | diff | sandwich | — | ✅ | test-simulation.R |
| bin | gauss | ATE (full) + EM | diff | sandwich | — | ✅ | test-effect-modification.R |
| bin | gauss | ATE (full) + EM | diff | boot | — | ✅ | test-effect-modification.R |

Rejections (all ✅ tested):
- categorical treatment $\to$ test-matching.R
- continuous treatment $\to$ test-matching.R
- non-static interventions $\to$ test-contrast.R
- multivariate $\to$ test-s3-methods.R
- longitudinal $\to$ test-matching.R

---

## ICE (longitudinal g-comp)

| Trt | Outcome | Intervention | Periods | Variance | Wt | Status | Test |
|---|---|---|---|---|---|---|---|
| bin | gauss | static | 2 | sandwich | — | ✅ (Table 20.1) | test-ice.R |
| bin | gauss | static | 2 | boot | — | ✅ | test-ice.R |
| bin | gauss | dynamic | 2 | sandwich | — | ✅ | test-ice.R, test-simulation.R |
| bin | gauss | static | 2 | sandwich | survey | ✅ | test-simulation.R |
| bin | gauss | static | 2 | sandwich | survey | ✅ vs unwt | test-ice.R, test-simulation.R |
| bin | binom | static | 3+ | sandwich | — | ✅ | test-ice.R, test-simulation.R |
| bin | binom | dynamic | 3+ | sandwich | — | ✅ | test-ice.R |
| bin | binom | static | 3+ | boot | — | ✅ | test-ice.R |
| bin | binom | static | 3+ | boot (parallel) | — | ✅ | test-ice.R |
| cont | gauss | shift | 2 | sandwich | — | ✅ (lmtp) | test-simulation.R |
| cont | gauss | scale_by | 2 | sandwich | — | ✅ (lmtp) | test-simulation.R |
| cont | gauss | threshold | 2 | sandwich | — | ✅ | test-simulation.R |
| cont | gauss | shift | 2 | boot | — | ✅ | test-simulation.R |
| bin | gauss | static | 2 | sandwich | cens | ✅ | test-ice.R |
| multi | gauss | static | 2 | sandwich | — | ✅ | test-multivariate.R |
| multi | gauss | shift | 2 | sandwich | — | ✅ | test-multivariate.R |
| bin | gauss | static | 2 | boot | survey | 🟡 | test-simulation.R |
| bin | gauss | static + EM | 2 | sandwich | — | ✅ | test-effect-modification.R |
| bin | gauss | static + EM | 3 | sandwich | — | ✅ | test-effect-modification.R |
| bin | gauss | static + EM | 2 | boot | — | ✅ | test-effect-modification.R |
| bin | gauss | static + multi-EM | 2 | sandwich | — | 🟡 | test-effect-modification.R |
| cont | gauss | shift + EM | 2 | sandwich | — | ✅ | test-effect-modification.R |
| cont | gauss | scale_by + EM | 2 | sandwich | — | 🟡 | test-effect-modification.R |
| cont | gauss | threshold + EM | 2 | sandwich | — | 🟡 | test-effect-modification.R |
| bin | gauss | dynamic + EM | 2 | sandwich | — | 🟡 | test-effect-modification.R |
| bin | binom | static + EM | 2 | sandwich | — | ✅ | test-effect-modification.R |
| multi | gauss | static + EM | 2 | sandwich | — | 🟡 | test-effect-modification.R |

Rejections (all ✅ tested): ipsi $\to$ test-ice.R, ATT/ATC $\to$ test-ice.R.

---

## Survival

Phase 7 defines two tracks: **Track A** = point survival via pooled logistic; **Track B** = longitudinal survival via ICE-hazards (iterated conditional hazards; Zivich et al. 2024 extended to the hazard link). See `PHASE_7_SURVIVAL.md` for the inheritance principle — survival is an outcome-model swap, not a new treatment-side axis, so most treatment-side features inherit as-is.

### Fit path (`causat_survival()`, pooled logistic hazard)

| Trt | Censoring | Competing | Status | Test |
|---|---|---|---|---|
| bin | none | none | 🟡 | test-s3-methods.R, test-simulation.R |
| bin | present (row filter) | none | 🟡 | test-simulation.R |

### Track A — Point survival (contrast path)

All ❌ pending Phase 7 chunks 7a–7d. `contrast()` on a survival fit aborts today ✅ (test-contrast.R).

| Estimator | Trt | Intervention | Estimand | Contrast | Variance | Status | Blocking chunk |
|---|---|---|---|---|---|---|---|
| gcomp | bin | static | ATE | S(t) / risk / RD / RR / RMST | sandwich | ❌ | 7a, 7b |
| gcomp | bin | static | ATE | S(t) / risk / RD | boot | ❌ | 7a, 7c |
| gcomp | bin | static | ATT | RD / RR | sandwich | ❌ | 7a, 7b |
| gcomp | bin | static | ATC | RD / RR | sandwich | ❌ | 7a, 7b |
| gcomp | cont | shift | ATE | S(t) / RD | sandwich | ❌ | 7a, 7b |
| gcomp | cont | scale_by | ATE | S(t) / RD | sandwich | ❌ | 7a, 7b |
| gcomp | cont | threshold | ATE | S(t) / RD | sandwich | ❌ | 7a, 7b |
| gcomp | bin | dynamic | ATE | S(t) / RD | sandwich | ❌ | 7a, 7b |
| gcomp | multi | static | ATE | S(t) / RD | sandwich | ❌ | 7a, 7b |
| gcomp | cat | static | ATE | S(t) / RD | sandwich | ❌ | 7a, 7b |
| gcomp | bin | static + by(L) | ATE | RD per stratum | sandwich | ❌ | 7a, 7b |
| gcomp | bin | static + EM | ATE | RD | sandwich | ❌ | 7a, 7b |
| gcomp | bin | static | ATE (competing risks) | CIF per cause | sandwich | ❌ | 7f |
| IPW | bin | static | ATE | S(t) / RD | sandwich | ❌ | 7d |
| IPW | bin | static | ATT/ATC | RD | sandwich | ❌ | 7d |
| IPW | cont | shift | ATE | S(t) / RD | sandwich | ❌ | 7d |
| IPW | bin | ipsi($\delta$) | ATE | S(t) / RD | sandwich | ❌ | 7d |

### Track B — Longitudinal survival via ICE hazards (contrast path)

All ❌ pending Phase 7 chunk 7e.

| Estimator | Trt | Intervention | Periods | Variance | Status | Blocking chunk |
|---|---|---|---|---|---|---|
| gcomp (ICE) | bin | static | 2 | sandwich | ❌ | 7e |
| gcomp (ICE) | bin | static | 3+ | sandwich | ❌ | 7e |
| gcomp (ICE) | bin | dynamic | 2 | sandwich | ❌ | 7e |
| gcomp (ICE) | cont | shift | 2 | sandwich | ❌ | 7e |
| gcomp (ICE) | cont | scale_by | 2 | sandwich | ❌ | 7e |
| gcomp (ICE) | bin | static | 2 | boot | ❌ | 7e, 7c |
| gcomp (ICE) | bin | static + EM | 2 | sandwich | ❌ | 7e |
| gcomp (ICE) | bin | static (competing risks) | 2+ | sandwich | ❌ | 7e, 7f |

### Cross-phase composition (survival × pending phases)

All ❌. Pattern fixed in the "Survival composition" subsections of the respective phase docs.

| Composition | Pattern | Status | Owning phase |
|---|---|---|---|
| Multivariate IPW + point survival | Joint density on original rows; product weight broadcast onto person-period; weighted hazard MSM | ❌ | Phase 8 |
| Longitudinal IPW + survival | Cumulative $(i, k)$-indexed weights; weighted hazard MSM on uncensored rows | ❌ | Phase 10 |
| Stochastic + point survival | MC draws at the individual cumulative-product level (not the hazard level) | ❌ | Phase 12 |
| Stochastic + longitudinal survival | Per-step MC with survival-tail pseudo-outcome (existing Phase 12 ICE architecture) | ❌ | Phase 12 |
| Built-in IPCW + point survival | Per-period stabilized censoring weights on person-period grid; weighted pooled logistic | ❌ | Phase 14 (primary use case) |
| Built-in IPCW + longitudinal survival | Per-period IPCW through ICE backward loop; cascade-gradient extension | ❌ | Phase 14 (primary use case) |
| Built-in IPCW + longitudinal IPW + survival | Product of cumulative treatment and censoring weights per $(i, k)$; stacked EE with $K+K+1$ blocks | ❌ | Phase 14 |
| Multivariate + longitudinal + survival | Three-way composition (Phase 8 + Phase 10 + Phase 7); deferred past the individual phases | ❌ | Phases 8+10, post-ship |

### Rejections

| Condition | Status | Test |
|---|---|---|
| `contrast()` on survival fit (pre-Phase 7) | ✅ | test-contrast.R |
| `competing != NULL` in `causat_survival()` | ✅ | test-causat.R |
| Matching + survival (architectural; see `PHASE_7_SURVIVAL.md` § "Matching + survival") | ❌ reject path pending | 7g |

---

## Missing data

### Rejections (all ✅)

| Condition | Test |
|---|---|
| Treatment NAs without `censoring=` | test-missing-data.R, test-causat.R |
| Weight NAs | test-weights-edge-cases.R, test-causat.R |
| `na.action = na.exclude` | test-causat.R |
| IPW: covariate NAs diverge from outcome mask | test-missing-data.R |

### G-comp + missing data

| Trt | Outcome | NA in | Intervention | Variance | Status | Test |
|---|---|---|---|---|---|---|
| bin | gauss | Y (15%) | static | sandwich | ✅ | test-missing-data.R |
| bin | binom | Y (15%) | static | sandwich | ✅ | test-missing-data.R |
| bin | gauss | L (10%) | static | sandwich | ✅ | test-missing-data.R |
| bin | gauss | Y+L (10%) | static | sandwich | ✅ | test-missing-data.R |
| bin | gauss | Y (10%) | static | boot | ✅ | test-missing-data.R |
| cont | gauss | L (10%) | shift | sandwich | ✅ | test-missing-data.R |
| bin | gauss | Y (10%) | static ATT | sandwich | ✅ | test-missing-data.R |
| bin | gauss | Y (10%) | static ATC | sandwich | ✅ | test-missing-data.R |
| bin | gauss | Y (10%) | static+by | sandwich | ✅ | test-missing-data.R |
| multi | gauss | Y (10%) | static | sandwich | ✅ | test-missing-data.R |
| bin | binom | Y (10%) | static | sandwich (ratio/OR) | 🟡 | test-missing-data.R |
| cat | gauss | Y (10%) | static | sandwich | ✅ | test-missing-data.R |
| bin | gauss | Y (10%) | static | sandwich | survey | ✅ | test-missing-data.R |
| bin | gauss | Y (10%) | static (GAM) | sandwich | ✅ | test-missing-data.R |
| bin | poisson | Y (10%) | static | sandwich (ratio) | ✅ | test-missing-data.R |

### IPW + missing data

| Trt | NA in | Intervention | Status | Test |
|---|---|---|---|---|
| bin | Y (15%) | static | ✅ | test-missing-data.R |
| bin | Y+L (same rows) | static | ✅ | test-missing-data.R |
| cont | Y (10%) | shift | ✅ | test-missing-data.R |
| cat | Y (10%) | static | ✅ | test-missing-data.R |

### Matching + missing data

| NA in | Estimand | Status | Test |
|---|---|---|---|
| Y (10%) | ATT | ✅ | test-missing-data.R |
| L (5%) | ATT | 🟡 | test-missing-data.R |

### ICE + missing data

| NA in | Intervention | Variance | Status | Test |
|---|---|---|---|---|
| Y final (10%) | static | sandwich | ✅ | test-missing-data.R |
| time-varying L (8%) | static | boot | ✅ | test-missing-data.R |
| Y final (10%) | static | boot | ✅ | test-missing-data.R |
| time-varying L | static | sandwich | ❌ (cascade gradient alignment) | — |

### MAR outcome / IPCW

| Scenario | Method | Status | Test |
|---|---|---|---|
| MAR, correct outcome model | test-gcomp.R (complete-case) | ✅ | test-missing-data.R |
| MAR, manual IPCW weights | test-gcomp.R (weighted) | ✅ | test-missing-data.R |
| MAR longitudinal, manual IPCW | ICE (weighted) | ✅ | test-missing-data.R |
| Built-in IPCW (Phase 14) | any | ❌ | — |

---

## Cross-cutting

| Concern | Status | Test |
|---|---|---|
| `to_person_period()` round-trip | ✅ | test-simulation.R |
| `to_person_period()` rejects dup ids / mismatched lengths | ✅ | test-simulation.R |
| `causat()` type auto-detection | ✅ | test-causat.R |
| `causat()` rejects bad inputs | ✅ | test-causat.R |
| `contrast()` rejects bad estimand/reference/intervention | ✅ | test-contrast.R |
| `diagnose()` gcomp/matching | ✅ 🟡 | test-diagnose.R |
| `diagnose()` IPW (binary static ATE shim) | 🟡 | test-diagnose.R |
| S3 methods (print/summary/plot/coef/vcov/confint/tidy/glance) | 🟡 | test-s3-methods.R |
| `confint()` respects `level` | ✅ | test-s3-methods.R |
| Weight validation (NA/Inf/neg/mis-sized) | ✅ | test-causat.R |
| Intervention constructors | ✅ | test-interventions.R |
| External reference cross-checks (stdReg2, delicatessen) | ✅ | test-variance-reference.R |
| Coverage matrix ↔ tests audit | ✅ | test-coverage-matrix.R |
| Numeric variance Tier 1 / Tier 2 | ✅ | test-variance-if.R |
| Cluster-robust test-matching.R variance | ✅ | test-variance-if.R, test-simulation.R |
| `causat_mice()` stub | 🟡 | test-causat-mice.R |
| Cross-method EM triangulation (gcomp + IPW + matching) | ✅ | test-effect-modification.R |

### Critical-review regression tests (2026-04-15)

All in `test-critical-review-2026-04.R`:

| Fix | Failure mode |
|---|---|
| B1 | `subset = quote()` could not resolve session vars |
| B2 | Bootstrap refit dropped user `...` |
| B5 | `causat_survival()` censoring dropped only current row |
| B6 | External weights post-multiplied instead of entering via `s.weights` |
| B7 | ICE Ch1 IF weighted/unweighted formula drift |
| B8 | `by` aborted on empty strata |
| R6 | OR validation aborted on NA `mu_hat` |
| R12 | Reserved column collision checked in only one place |

Supplementary: `test-weights-edge-cases.R` (external weights edge cases), `test-replay-fit.R` (`replay_fit()` helper).

---

## Planned (future phases)

### Phase 6 — Effect modification
Unified `A:modifier` API across gcomp / IPW / matching / ICE. IPW MSM expansion ✅ (chunk 6b). Matching MSM expansion ✅ (chunk 6c). ICE lag auto-expansion ✅ (chunk 6d). Cross-method triangulation ✅ (chunk 6e). **Phase complete.** **Known limitation** (documented in `PHASE_6_INTERACTIONS.md` § "Known limitation: modifier must be **baseline**"): under IPW / matching / longitudinal IPW, the modifier must be baseline because the MSM conditions on post-treatment variables otherwise (silent bias; Robins 2000). Time-varying effect modification is deferred to Phase 18 (SNMs). Currently doc-level only; a `check_em_baseline_only()` runtime guard via an explicit `baseline_cols =` contract is a follow-up.

### Phase 7 — Survival analysis
Two tracks: Track A (point survival via pooled logistic), Track B (longitudinal survival via ICE hazards — Zivich et al. 2024 extended to the hazard link). Plus competing risks (cause-specific hazards + CIF), survival-curve contrast shape ($S(t)$ / risk / RD / RR / RMST as time-indexed `data.table` rows), cross-time variance via delta-method on the cumulative product, matching rejection path, NHEFS Ch. 17 replication. Most treatment-side features (binary / continuous / categorical / count / multivariate, all IPW-legal interventions, ATE / ATT / ATC where the estimator allows, effect modification, external weights, `by`-stratification) inherit from gcomp/IPW/ICE as-is. All ❌.

### Phase 8 — Multivariate treatment IPW
Joint density via sequential factorisation, product density-ratio weights, multi-model propensity sandwich. Plus survival composition: baseline joint weight broadcast onto person-period grid, weighted pooled-logistic hazard MSM, cross-time delta IF aggregation (depends on Phase 7 Track A). All ❌.

### Phase 9 — Inference infrastructure
Survey design integration, general cluster-robust sandwich, `future` backend. Mix of ❌ and partial.

### Phase 10 — Longitudinal IPW
Sequential density-ratio weights, cumulative product weights, stabilized weights, time-varying MSM. Plus survival composition: cumulative $(i, k)$-indexed weights, weighted hazard MSM on uncensored person-period rows, cumulative-product survival curve, stacked sandwich with cross-time delta (depends on Phase 7 chunks 7a–7b and 7e). All ❌.

### Phase 11 — diagnose() rewrite
Intervention-aware, treatment-type-aware, estimand-aware, longitudinal-aware diagnostics. All ❌.

### Phase 12 — Stochastic interventions
`stochastic()` under gcomp (point + ICE), MC g-formula, MC-averaged IFs. IPW/matching rejected. Plus survival composition: MC draws at the individual cumulative-product level for Track A (not at the hazard level — nonlinear cumulative product requires MC-after-cumulating); per-step MC with survival-tail pseudo-outcomes for Track B (existing per-step architecture). All ❌.

### Phase 13 — Outcome types
Negative binomial tests, beta regression (`resolve_family("beta")`), multinomial/ordinal outcomes. All ❌.

### Phase 14 — Built-in IPCW
**Survival is the motivating use case** (censoring is the rule in time-to-event data). Survival-first chunks: internal per-period censoring model, cumulative stabilized IPCW weights on person-period grid, weighted hazard MSM, stacked EE sandwich extension for censoring model blocks (Track A + Track B), Phase 10 composition (cumulative treatment × cumulative censoring weights). Non-survival MAR outcomes drop in as a one-period special case. All ❌.

### Phase 15 — Polish and documentation
Continuous treatment vignette, `target_trial()` helper, misc release-prep. All ❌.

### Phase 16 — AIPW / doubly-robust estimator
Composes Phase 2 gcomp + Phase 4 IPW into the classical analytical doubly-robust estimator: $\hat\psi_{\mathrm{AIPW}}(a) = E[\hat{m}(a, L)] + E[W \cdot (Y - \hat{m}(A, L))]$. Stacked EE sandwich (outcome block + propensity block + plug-in); distinct from `lmtp` (TMLE/SDR with ML + cross-fitting). Double-robustness + efficiency tests; survival composition (Phase 7 Track A under AIPW); ICE-AIPW longitudinal extension (Bang & Robins 2005). All ❌.

### Phase 17 — Transportability / Generalizability
Sampling model `P(S=1 | L)` + sampling-odds weights multiply into IPW Hájek MSM; gcomp / IPW / AIPW transport paths; stacked EE extends with sampling-model block. References: Dahabreh et al. 2020; Westreich et al. 2017. Cross-check against `transport` / `transported` R packages. Survival + transport composition subsection. All ❌.

### Phase 18 — G-estimation of Structural Nested Mean Models
Third leg of the Robins triangle. Motivating use case: **correct handling of time-varying effect modification** — SNMs parameterise the per-stage blip $\gamma_k(a_k, \bar{l}_k, \bar{a}_{k-1}; \psi)$ directly and identify it via a moment condition that uses the treatment model as instrument, so time-varying modifiers are supported by design (closes the Phase 6 limitation under MSM-based estimators). Scope: linear-blip additive SNMMs for point + longitudinal, stacked EE sandwich ($K$ treatment blocks + blip block), bootstrap, `gesttools` cross-check. Survival SNMs (SNFTMs/SNCFTMs) out of scope. All ❌.

### `causat_mice()` — Multiple imputation
Pool across `mice` imputations via Rubin's rules. All ❌.

---

## Maintenance rules

1. Update this matrix when adding/removing/changing a feature.
2. Add truth-based tests when feasible; use lmtp/delicatessen as external reference.
3. Test rejection paths with `expect_snapshot(error = TRUE)`.
4. Divergence between this matrix and test files is a bug.
