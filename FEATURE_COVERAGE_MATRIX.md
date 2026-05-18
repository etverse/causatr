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
| multi (bin × bin) | gauss | static | ATE | diff | sandwich | — | ✅ +gcomp | test-multivariate-ipw.R |
| multi (bin × bin) | gauss | static | ATE | diff | boot | — | ✅ | test-multivariate-ipw.R |
| multi (bin × cont) | gauss | static + shift | ATE | diff | sandwich | — | ✅ | test-multivariate-ipw.R |
| multi (cont × cont) | gauss | shift + shift | ATE | diff | sandwich | — | ✅ seq-MTP | test-multivariate-ipw.R |
| multi (bin × bin) | binom | static | ATE | diff/ratio/OR | sandwich | — | ✅ | test-multivariate-ipw.R |
| multi (bin × bin) | gauss | static | by(L) | diff | sandwich | — | ✅ | test-multivariate-ipw.R |
| multi (bin × bin) | gauss | static | subset | diff | sandwich | — | ✅ | test-multivariate-ipw.R |
| multi (bin × bin) | gauss | dynamic + static | ATE | diff | sandwich | — | ✅ | test-multivariate-ipw.R |
| multi (bin³) | gauss | static | ATE | diff | sandwich | — | ✅ | test-multivariate-ipw.R |
| multi (bin × bin) | gauss | static + EM (`A1:sex`) | by(sex) | diff | sandwich | — | ✅ +gcomp | test-multivariate-ipw.R |
| multi (bin × cat) | gauss | static | ATE | diff | sandwich | — | ✅ +gcomp | test-multivariate-ipw.R |
| multi (cat × bin) | gauss | static | ATE | diff | sandwich | — | ✅ | test-multivariate-ipw.R |
| multi (cat × cat) | gauss | static | ATE | diff | sandwich | — | ✅ | test-multivariate-ipw.R |
| multi (bin × poisson) | gauss | static + shift | ATE | diff | sandwich | — | ✅ | test-multivariate-ipw.R |
| multi (bin × negbin) | gauss | static + shift | ATE | diff | sandwich | — | ✅ | test-multivariate-ipw.R |
| multi (poisson × cont) | gauss | shift + natural | ATE | diff | sandwich | — | ✅ seq-MTP | test-multivariate-ipw.R |
| multi (bin × bin) | gauss | static + stabilize="marginal" | ATE | diff | sandwich | — | ✅ | test-multivariate-ipw.R |
| multi (cont × cont) | gauss | shift + shift + stabilize="marginal" | ATE | diff | sandwich | — | ✅ | test-multivariate-ipw.R |
| multi (bin × bin) | gauss | static + stabilize="marginal" | ATE | diff | bootstrap | — | ✅ | test-multivariate-ipw.R |
| multi (bin × cat) | gauss | static + stabilize="marginal" | ATE | diff | sandwich | — | ✅ | test-multivariate-ipw.R |

**Estimand note.** Multivariate IPW implements the *sequential MTP* estimand (Díaz et al. 2023, `lmtp`). Multivariate gcomp implements deterministic joint transformation. For static-only interventions the two coincide; for non-static interventions on non-final components they can differ by the upstream $\to$ downstream cross-dependence (e.g. shift $+$ shift on a cont $\times$ cont DGP with $A_1 \to A_2$ coefficient $0.3$: gcomp gives $-0.9$, IPW gives $-1.02$).

Variance internals: self-contained IF ✅ (hand-derived cross-derivative + end-to-end stacked-sandwich; test-ipw-branch-b.R, test-ipw-cross-derivative.R, test-variance-if.R). Non-static variance regression ✅ (shift ~5-8% SE reduction, IPSI ~90% off-diagonal covariance; test-ipw-variance-regression.R). Bootstrap parity ✅ (within 30% MC tolerance; test-ipw-variance-regression.R).

Rejections (all ✅ tested):
- `static`/`threshold`/`dynamic` on continuous $\to$ test-ipw.R
- `static`/`threshold`/`dynamic`/`ipsi` on count $\to$ test-ipw-count.R
- non-integer `shift()` on count $\to$ test-ipw-count.R
- non-integer-preserving `scale_by()` on count $\to$ test-ipw-count.R
- shift/scale_by/dynamic/ipsi + ATT/ATC $\to$ test-estimand-intervention-compat.R
- multivariate + ATT/ATC $\to$ test-multivariate.R, test-multivariate-ipw.R
- multivariate + ipsi() in any component $\to$ test-multivariate-ipw.R
- multivariate + bare treatment in confounders (`~ L + A1`) $\to$ test-multivariate-ipw.R
- multivariate + `propensity_family` invalid shape $\to$ test-multivariate-ipw.R
- multivariate count + `static()` / `threshold()` / `dynamic()` / non-integer `shift()` $\to$ test-multivariate-ipw.R
- univariate IPW + `stabilize = "marginal"` $\to$ test-multivariate-ipw.R (`causatr_stabilize_univariate`)

---

## Longitudinal IPW (Phase 10)

Per-period treatment density chain $f(A_k \mid \bar A_{k-1}, \bar L_k)$ + cumulative product weights + final-period intercept-only Hájek MSM + stacked sandwich (block-diagonal propensity bread across periods) + id-clustered bootstrap.

| Trt | Outcome | Intervention | Periods | Variance | Wt | Status | Test |
|---|---|---|---|---|---|---|---|
| cont | gauss | shift | 2 | sandwich | — | ✅ vs ICE | test-longitudinal-ipw.R |
| cont | gauss | shift | 2 | bootstrap | — | ✅ vs sandwich | test-longitudinal-ipw.R |
| cont | gauss | shift | 2 | sandwich | — | ✅ vs lmtp::lmtp_sdr | test-longitudinal-ipw.R |
| bin | gauss | static (always vs never) | 2 | sandwich | — | ✅ vs ICE | test-longitudinal-ipw.R |
| bin | gauss | natural course (NULL) | 2 | sandwich | — | ✅ recovers observed mean | test-longitudinal-ipw.R |
| bin | gauss | static (always vs never) + stabilize="marginal" | 2 | sandwich | — | ✅ identical to unstabilized | test-longitudinal-ipw.R |
| bin | gauss | static (always vs never) + stabilize="marginal" | 2 | bootstrap | — | ✅ vs sandwich | test-longitudinal-ipw.R |
| cont | gauss | shift + stabilize="marginal" + numerator=~L0 | 2 | sandwich | — | ✅ vs ICE | test-longitudinal-ipw.R |
| (numerator structure) | — | default `~ 1` / `~ lag1_A` per period | 2 | — | — | ✅ formula check | test-longitudinal-ipw.R |
| (numerator structure) | — | custom `numerator = ~ L0` keeps treatment lags | 2 | — | — | ✅ formula check | test-longitudinal-ipw.R |
| bin | gauss | static (always vs never) + EM (`A:sex`) | 2 | sandwich | — | ✅ ATE\|sex=0 = 5, ATE\|sex=1 = 8 | test-longitudinal-ipw.R |
| bin | gauss | static (always vs never) + EM (`A:sex`) | 2 | sandwich | — | ✅ vs ICE EM | test-longitudinal-ipw.R |
| (EM structure) | — | per-period propensity strips `A:sex`; MSM expands to `Y ~ 1 + sex` | 2 | — | — | ✅ formula check | test-longitudinal-ipw.R |

**Known limitation (Phase 6, Robins 2000):** under longitudinal IPW the modifier MUST be a **baseline** covariate. A time-varying modifier in an MSM conditions on a post-treatment variable, biasing the estimand via mediator + collider paths. Not enforced at runtime (time-varying status isn't inferable from data); doc-level constraint only. The scientifically correct tool for time-varying effect modification is a structural nested model (Phase 18).

Sequential positivity warning (`causatr_longitudinal_seq_positivity`): fires when any per-period weight max exceeds 100 (default threshold); silent below threshold. Tested directly on `warn_seq_positivity()`.

Rejections (all ✅ tested, "_pending" classed errors deferred to follow-up chunks):
- `ipsi()` $\to$ test-longitudinal-ipw.R (`causatr_longitudinal_ipsi_pending`; per-period IPSI extension deferred)
- `numerator =` without `stabilize = "marginal"` $\to$ test-longitudinal-ipw.R (`causatr_longitudinal_numerator_without_stabilize`)
- bare treatment in confounders (`~ L + A`) $\to$ test-longitudinal-ipw.R (`causatr_bare_treatment_in_confounders`)
- multivariate longitudinal IPW $\to$ test-longitudinal-ipw.R (`causatr_longitudinal_multivariate_pending`)
- ATT / ATC under longitudinal IPW $\to$ test-longitudinal-ipw.R (inherits `check_estimand_trt_compat`)
- single-period data labelled `type = "longitudinal"` $\to$ test-longitudinal-ipw.R (`causatr_longitudinal_too_few_times`)

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

Causal survival analysis has been **moved out of causatr** into a
separate etverse package (`etverse/survatr`). No survival rows
remain in this matrix.

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
| Built-in IPCW: censoring model primitive | fit + weights (14a) | ✅ | test-censoring.R |
| Built-in IPCW: point estimators | gcomp/IPW/matching (14b) | ✅ | test-ipcw.R |
| Built-in IPCW: lmtp cross-check | point + longitudinal (14d) | ✅ | test-ipcw-lmtp-oracle.R |
| Built-in IPCW: ICE longitudinal | ICE + long IPW (14c) | ✅ | test-ipcw.R, test-ipcw-lmtp-oracle.R |
| Built-in IPCW: variance regression | sandwich + bootstrap (14e) | ✅ | test-ipcw-variance.R |
| Built-in IPCW: diagnose integration | point + longitudinal (14f) | ✅ | test-diagnose.R |

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
| `diagnose()` IPW binary ATE — default `obs` panel | ✅ | test-diagnose.R |
| `diagnose()` IPW binary ATE — per-intervention panels (`interventions = list(...)`) | ✅ | test-diagnose.R |
| `diagnose()` IPW continuous ATE — density-range positivity + overall weight summary | ✅ | test-diagnose.R |
| `diagnose()` IPW continuous ATE — shift() weight reconstruction | ✅ | test-diagnose.R |
| `diagnose()` IPW categorical ATE — per-level positivity + overall weight summary | ✅ | test-diagnose.R |
| `diagnose()` IPW count ATE (Poisson) — density-range positivity + weight summary | ✅ | test-diagnose.R |
| `diagnose()` IPW count ATE (negbin) — density-range positivity + weight summary | ✅ | test-diagnose.R |
| `diagnose()` IPW multivariate ATE — per-component positivity + product weight | ✅ | test-diagnose.R |
| `diagnose()` IPW multivariate ATE — per-component interventions | ✅ | test-diagnose.R |
| `diagnose()` longitudinal ICE — per-period balance, no weights/positivity | ✅ | test-diagnose.R |
| `diagnose()` longitudinal IPW binary — per-period positivity + weights + balance | ✅ | test-diagnose.R |
| `diagnose()` longitudinal IPW continuous — shift intervention per-period diagnostics | ✅ | test-diagnose.R |
| `diagnose()` longitudinal IPW — cumulative weight matches manual product | ✅ | test-diagnose.R |
| `diagnose()` IPW ATT — treated weight = 1, control = p/(1-p) | ✅ | test-diagnose.R |
| `diagnose()` IPW ATC — control weight = 1, treated = (1-p)/p | ✅ | test-diagnose.R |
| `diagnose()` IPW ATT ESS sanity: ESS_treated = n_treated | ✅ | test-diagnose.R |
| `diagnose()` `by =` stratified balance via cobalt cluster | ✅ | test-diagnose.R |
| `diagnose()` `by =` validation (missing column, non-scalar) | ✅ | test-diagnose.R |
| S3 methods (print/summary/plot/coef/vcov/confint/tidy/glance) | 🟡 | test-s3-methods.R |
| `knit_print.causatr_result` + `.onLoad` registration | ✅ | test-knit-print.R |
| Input-validation helpers (`check_string`, `check_formula`, `check_col_exists`, `check_pkg`) | ✅ | test-checks.R |
| `prepare_data()` + `warn_confounder_variation()` | ✅ | test-prepare-data.R |
| Bootstrap result processor + `refit_model()` dispatch | ✅ | test-variance-bootstrap.R |
| `apply_intervention_to_values()` (IPW HT branch helper) | ✅ | test-ipw-weights.R |
| `make_weight_fn()` defensive guards | ✅ | test-ipw-weights.R |
| `apply_intervention()` multivariate dispatch | ✅ | test-interventions.R |
| `extract_sigma()` GAM + residual fallbacks | ✅ | test-treatment-model.R |
| `evaluate_density()` + `evaluate_categorical_density()` defensive guards | ✅ | test-treatment-model.R |
| `format_family()` display helper | ✅ | test-s3-methods.R |
| `apply_single_intervention()` dynamic-rule type guards | ✅ | test-interventions.R |
| `check_intervention_list()` structural validation | ✅ | test-checks.R |
| `check_causat_inputs()` (`treatment` shape, outcome==treatment, `history`) | ✅ | test-checks.R |
| `check_intervention_family_compat()` rejection branches | ✅ | test-ipw-weights.R |
| `compute_density_ratio_weights()` defensive guards | ✅ | test-ipw-weights.R |
| `ht_bayes_numerator()` per-estimand (ATE/ATT/ATC) + error guards | ✅ | test-ipw-weights.R |
| `fit_treatment_model()` NB-theta + aliased-column aborts | ✅ | test-treatment-model.R |
| `bread_inv()` GAM-missing-`$Vp` abort | ✅ | test-variance-if.R |
| `combine_match_and_external_weights()` row-name invariant | ✅ | test-matching.R |
| `compute_positivity_binary()` short-circuits (multivariate / non-binary / non-bernoulli IPW) | ✅ | test-diagnose.R |
| `compute_contrast()` ratio/OR validation aborts + empty target | ✅ | test-contrast.R |
| `confint()` respects `level` | ✅ | test-s3-methods.R |
| Weight validation (NA/Inf/neg/mis-sized) | ✅ | test-causat.R |
| Intervention constructors | ✅ | test-interventions.R |
| External reference cross-checks (stdReg2, delicatessen) | ✅ | test-variance-reference.R |
| `target_trial()` constructor + S3 print | ✅ | test-target-trial.R |
| `diagnose()` feasibility (pct_intervened) — static/shift/threshold/ipsi/NULL | ✅ | test-diagnose.R |
| `diagnose()` feasibility — longitudinal per-period + overall | ✅ | test-diagnose.R |
| `ice_build_formula()` transformed `confounders_tv` (`poly()`, `log()`) | ✅ | test-ice.R |
| `substitute_vars_in_term()` + `is_bare_term()` + `term_vars()` helpers | ✅ | test-ice.R |
| Coverage matrix ↔ tests audit | ✅ | test-coverage-matrix.R |
| Numeric variance Tier 1 / Tier 2 | ✅ | test-variance-if.R |
| Cluster-robust test-matching.R variance | ✅ | test-variance-if.R, test-simulation.R |
| General cluster-robust sandwich (`contrast(cluster = ...)`) — gcomp + sandwich::vcovCL oracle | ✅ | test-cluster-sandwich.R |
| General cluster-robust sandwich — IPW, ICE, by-stratified | ✅ | test-cluster-sandwich.R |
| `causat(cluster = ...)` preservation through `prepare_data()` | ✅ | test-cluster-sandwich.R |
| Cluster rejections (matching + fit/contrast, unknown col, wrong length, NA) | ✅ | test-cluster-sandwich.R |
| `causat(weights = survey::svydesign(...))` — weights + PSU auto-extract | ✅ | test-survey-design.R |
| svydesign + matching rejected, row-count mismatch rejected | ✅ | test-survey-design.R |
| `contrast(parallel = "future")` — gcomp + ICE via `future::plan()` | ✅ | test-future-backend.R |
| `dispatch_boot()` — future vs boot MC equivalence | ✅ | test-future-backend.R |
| `causat_mice()` stub | 🟡 | test-causat-mice.R |
| Cross-method EM triangulation (gcomp + IPW + matching) | ✅ | test-effect-modification.R |

### Critical-review regression tests (2026-04-15)

All in `test-critical-review-2026-04.R`:

| Fix | Failure mode |
|---|---|
| B1 | `subset = quote()` could not resolve session vars |
| B2 | Bootstrap refit dropped user `...` |
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

### Phase 8 — Multivariate treatment IPW
Joint density via sequential factorisation, product density-ratio weights under sequential-MTP semantics (Díaz et al. 2023), multi-model propensity sandwich (block-diagonal bread). **Done** (chunks 8a, 8a-fix, 8b, 8c, 8d, 8e). Supports binary / continuous / categorical / count components, mixed types, K = 2..K, effect modification (`A_k:modifier`), and optional stabilized weights (`stabilize = "marginal"` with per-component numerator models $g_k(A_k \mid A_{1..k-1})$). Rejects IPSI. See multi rows in the IPW table above.

### Phase 9 — Inference infrastructure
Survey design integration, general cluster-robust sandwich, `future` backend. **All shipped.** General cluster-robust sandwich ✅ (`cluster =` on `causat()` and `contrast()`; `sandwich::vcovCL` oracle). `svydesign` integration ✅ (`causat(weights = svydesign_obj)`). Parallel bootstrap ✅. `future` backend ✅ (`parallel = "future"`).

### Phase 10 — Longitudinal IPW
Sequential density-ratio weights, cumulative product weights, stabilized weights, time-varying MSM. All ❌.

### Phase 11 — diagnose() rewrite
Intervention-aware, treatment-type-aware, estimand-aware, longitudinal-aware diagnostics. **All chunks shipped (11a–11e).** 11a (foundation): nested per-intervention `causatr_diag` shape, `interventions =` argument, binary IPW per-intervention weight summary. 11b (treatment-type dispatch): density-range positivity for continuous/count, per-level positivity for categorical, per-component positivity for multivariate, overall weight summaries for all non-binary types, truth-based weight reconstruction test. 11c (longitudinal dispatch): per-period positivity/balance/weight diagnostics for longitudinal IPW, per-period balance for ICE, cumulative weight summary. 11d (estimand + EM): IPW ATT/ATC observed-treatment weights, `by =` stratified balance via cobalt `cluster`. 11e (plot overhaul + vignette): `plot.causatr_diag()` with `which = c("balance", "weights", "positivity")`, `diagnostics.qmd` vignette.

| Feature | Status | Test |
|---|---|---|
| plot(diag, which="balance") — Love plot (IPW) | 🟡 | test-diagnose.R |
| plot(diag, which="balance") — Love plot (matching) | 🟡 | test-diagnose.R |
| plot(diag, which="weights") — weight histogram | 🟡 | test-diagnose.R |
| plot(diag, which="weights", log_scale=TRUE) | 🟡 | test-diagnose.R |
| plot(diag, which="positivity") — binary PS density | 🟡 | test-diagnose.R |
| plot(diag, which="positivity") — continuous density | 🟡 | test-diagnose.R |
| plot(diag, which="balance") errors for gcomp | ⛔ | test-diagnose.R |
| plot(diag, which="weights") warns for gcomp (no wts) | ⛔ | test-diagnose.R |

### Phase 12 — Stochastic interventions
`stochastic()` under gcomp (point + ICE), MC g-formula, MC-averaged IFs. IPW/matching rejected.

| Feature | Status | Test file |
|---|---|---|
| `stochastic()` constructor + validation | ✅ | test-stochastic.R |
| `stochastic()` n_mc < 10 warning | ✅ | test-stochastic.R |
| `stochastic()` n_mc coercion | ✅ | test-stochastic.R |
| `apply_single_intervention()` stochastic (numeric) | ✅ | test-stochastic.R |
| `apply_single_intervention()` stochastic (factor) | ✅ | test-stochastic.R |
| Stochastic rejected under IPW | ✅ | test-stochastic.R |
| Stochastic rejected under matching | ✅ | test-stochastic.R |
| Point gcomp × binary trt × gaussian × ATE (analytical truth) | ✅ | test-stochastic.R |
| Point gcomp × binary trt × gaussian × ATE (bootstrap) | ✅ | test-stochastic.R |
| Point gcomp × continuous trt × gaussian × ATE (analytical truth) | ✅ | test-stochastic.R |
| Point gcomp × continuous trt × gaussian × ATE (txshift oracle) | ✅ | test-stochastic.R |
| Point gcomp × binary trt × binomial × ATE (analytical truth) | ✅ | test-stochastic.R |
| Point gcomp × categorical trt × gaussian × ATE (analytical truth) | ✅ | test-stochastic.R |
| Point gcomp × multivariate trt × gaussian × ATE (analytical truth) | ✅ | test-stochastic.R |
| Point gcomp × binomial × ratio contrast | ✅ | test-stochastic.R |
| Point gcomp × binomial × OR contrast | ✅ | test-stochastic.R |
| Point gcomp × ATT estimand | ✅ | test-stochastic.R |
| Point gcomp × ATC estimand | ✅ | test-stochastic.R |
| Point gcomp × by-stratified estimand | ✅ | test-stochastic.R |
| Point gcomp × subset estimand | ✅ | test-stochastic.R |
| Point gcomp × GAM model | ✅ | test-stochastic.R |
| Point gcomp × Poisson × ratio contrast | ✅ | test-stochastic.R |
| Point gcomp × n_mc = 1 degenerate draw | ✅ | test-stochastic.R |
| Point gcomp × set.seed() reproducibility | ✅ | test-stochastic.R |
| Point gcomp × mixed stochastic + static interventions | ✅ | test-stochastic.R |
| MC-averaged gradient in sandwich variance (point) | ✅ | test-stochastic.R |
| ICE × binary trt × gaussian × 2 periods (simulation truth) | ✅ | test-stochastic.R |
| ICE × binary trt × gaussian × 2 periods (bootstrap) | ✅ | test-stochastic.R |
| ICE × continuous trt × gaussian × 2 periods (simulation truth) | ✅ | test-stochastic.R |
| ICE × set.seed() reproducibility | ✅ | test-stochastic.R |
| ICE × mixed stochastic + static interventions | ✅ | test-stochastic.R |
| MC-averaged gradient in ICE sandwich variance | ✅ | test-stochastic.R |
| Sandwich vs bootstrap SE agreement (point, binary, gaussian) | ✅ | test-stochastic.R |
| Sandwich vs bootstrap SE agreement (point, binary, binomial) | ✅ | test-stochastic.R |
| Sandwich vs bootstrap SE agreement (point, continuous) | ✅ | test-stochastic.R |
| Sandwich vs bootstrap SE agreement (point, categorical) | ✅ | test-stochastic.R |
| Sandwich vs bootstrap SE agreement (point, multivariate) | ✅ | test-stochastic.R |
| Sandwich vs bootstrap SE agreement (ICE) | ✅ | test-stochastic.R |
| lmtp_sdr cross-check (point treatment) | ✅ | test-stochastic.R |
| lmtp_sdr cross-check (longitudinal) | ✅ | test-stochastic.R |

### Phase 13 — Extended outcome type coverage

Comprehensive truth-based + cross-estimator tests for Poisson, Gamma, quasibinomial, negative binomial, and beta regression outcomes.

#### Poisson outcome (expanded)

| Estimator | Trt | Intervention | Estimand | Contrast | Variance | Status | Test |
|---|---|---|---|---|---|---|---|
| gcomp | bin | static | ATE | diff | sandwich | ✅ | test-simulation.R |
| gcomp | bin | static | ATE | ratio | sandwich | ✅ | test-gcomp.R (existing) |
| gcomp | bin | static | ATE | OR | sandwich | ⛔ rejection | test-simulation.R |
| gcomp | bin | static | ATE | diff | bootstrap | ✅ SE agreement | test-simulation.R |
| gcomp | bin | static | ATT | diff | sandwich | ✅ | test-simulation.R |
| gcomp | cont | shift | ATE | diff | sandwich | ✅ | test-simulation.R |
| gcomp | cont | scale_by | ATE | diff | sandwich | 🟡 | test-simulation.R |
| ipw | bin | static | ATE | diff | sandwich | ✅ vs gcomp | test-simulation.R |
| ipw | bin | static | ATE | ratio | sandwich | ✅ vs gcomp | test-simulation.R |
| matching | bin | — | ATT | diff | sandwich | ✅ | test-simulation.R (existing) |

#### Gamma outcome (expanded)

| Estimator | Trt | Intervention | Estimand | Contrast | Variance | Status | Test |
|---|---|---|---|---|---|---|---|
| gcomp | bin | static | ATE | ratio | sandwich | ✅ | test-simulation.R (existing) |
| gcomp | bin | static | ATE | diff | sandwich | ✅ | test-simulation.R |
| gcomp | bin | static | ATE | OR | sandwich | ⛔ rejection | test-simulation.R |
| gcomp | bin | static | ATE | diff | bootstrap | ✅ SE agreement | test-simulation.R |
| gcomp | bin | static | ATT | diff | sandwich | ✅ | test-simulation.R |
| gcomp | cont | shift | ATE | diff | sandwich | ✅ | test-simulation.R |
| gcomp | cont | scale_by | ATE | ratio | sandwich | 🟡 | test-simulation.R |
| ipw | bin | static | ATE | diff | sandwich | ✅ vs gcomp | test-simulation.R |
| ipw | bin | static | ATE | ratio | sandwich | ✅ vs gcomp | test-simulation.R |
| matching | bin | — | ATT | diff | sandwich | 🟡 | test-simulation.R |

#### Quasibinomial outcome (expanded)

| Estimator | Trt | Intervention | Estimand | Contrast | Variance | Status | Test |
|---|---|---|---|---|---|---|---|
| gcomp | bin | static | ATE | diff | sandwich | ✅ | test-simulation.R (existing) |
| gcomp | bin | static | ATE | ratio | sandwich | ✅ | test-simulation.R |
| gcomp | bin | static | ATE | OR | sandwich | ✅ | test-simulation.R |
| gcomp | bin | static | ATE | diff | bootstrap | ✅ SE agreement | test-simulation.R |
| gcomp | bin | static | ATT | diff | sandwich | ✅ | test-simulation.R |
| gcomp | cont | shift | ATE | diff | sandwich | ✅ | test-simulation.R |
| gcomp | cont | scale_by | ATE | diff | sandwich | 🟡 | test-simulation.R |
| ipw | bin | static | ATE | diff | sandwich | ✅ vs gcomp | test-simulation.R |
| ipw | bin | static | ATE | OR | sandwich | ✅ vs gcomp | test-simulation.R |
| matching | bin | — | ATT | diff | sandwich | ✅ | test-simulation.R (existing) |

#### Negative binomial outcome (new)

| Estimator | Trt | Intervention | Estimand | Contrast | Variance | Status | Test |
|---|---|---|---|---|---|---|---|
| gcomp | bin | static | ATE | diff | sandwich | ✅ truth | test-simulation.R |
| gcomp | bin | static | ATE | ratio | sandwich | ✅ truth | test-simulation.R |
| gcomp | bin | static | ATE | OR | sandwich | ⛔ rejection | test-simulation.R |
| gcomp | bin | static | ATE | diff | bootstrap | ✅ SE agreement | test-simulation.R |
| gcomp | bin | static | ATT | diff | sandwich | ✅ | test-simulation.R |
| gcomp | cont | shift | ATE | diff | sandwich | ✅ | test-simulation.R |
| gcomp | cont | scale_by | ATE | diff | sandwich | 🟡 | test-simulation.R |
| ipw | bin | static | ATE | diff | sandwich | ✅ vs gcomp | test-simulation.R |
| ipw | bin | static | ATE | ratio | sandwich | ✅ vs gcomp | test-simulation.R |
| ipw | bin | static | ATT | diff | sandwich | 🟡 | test-simulation.R |
| matching | bin | — | ATT | diff | sandwich | 🟡 | test-simulation.R |

#### Beta regression outcome (new)

| Estimator | Trt | Intervention | Estimand | Contrast | Variance | Status | Test |
|---|---|---|---|---|---|---|---|
| gcomp | bin | static | ATE | diff | sandwich | ✅ truth | test-simulation.R |
| gcomp | bin | static | ATE | ratio | sandwich | ✅ truth | test-simulation.R |
| gcomp | bin | static | ATE | OR | sandwich | ✅ truth | test-simulation.R |
| gcomp | bin | static | ATE | diff | bootstrap | ✅ SE agreement | test-simulation.R |
| gcomp | bin | static | ATT | diff | sandwich | ✅ | test-simulation.R |
| gcomp | bin | static | ATE | diff | sandwich | ✅ family="beta" | test-simulation.R |
| gcomp | cont | shift | ATE | diff | sandwich | ✅ | test-simulation.R |
| gcomp | cont | scale_by | ATE | diff | sandwich | 🟡 | test-simulation.R |
| ipw | bin | static | ATE | diff | sandwich | ✅ vs gcomp | test-simulation.R |
| ipw | bin | static | ATE | ratio | sandwich | ✅ vs gcomp | test-simulation.R |
| ipw | bin | static | ATT | diff | sandwich | 🟡 | test-simulation.R |
| matching | bin | — | ATT | diff | sandwich | 🟡 | test-simulation.R |

#### Infrastructure

| Feature | Status | Test |
|---|---|---|
| `fn_accepts_family()` helper | ✅ | test-checks.R |
| `resolve_family("beta")` extension | ✅ | test-checks.R |
| betareg Jacobian fix in `variance_if_numeric()` | ✅ | test-simulation.R (beta SE tests) |

#### Not in scope (Phase 13)
Multinomial/ordinal outcomes — requires structural extensions beyond GLM families.

### Phase 14 — Built-in IPCW
Scalar-outcome IPCW for MAR censoring: internal censoring model, stabilized IPCW weights, weighted fit, stacked EE sandwich extension for censoring model blocks. Point + ICE scalar final outcome. Survival-specific IPCW (per-period cumulative weights + hazard MSM) is owned by the separate survival package. All ❌.

### Phase 15 — Polish and documentation

| Feature | Status | Test |
|---|---|---|
| `target_trial()` constructor + print method | ✅ | test-target-trial.R |
| Confounder selection warning (colliders/LASSO/dagitty) | ✅ | vignettes/introduction.qmd (callout) |
| ML g-formula debiasing warning | ✅ | vignettes/gcomp.qmd (callout) |
| Feasibility diagnostic (`compute_pct_intervened()`) — point | ✅ | test-diagnose.R |
| Feasibility diagnostic — longitudinal | ✅ | test-diagnose.R |
| Feasibility diagnostic — ipsi/NULL returns NULL | ✅ | test-diagnose.R |
| ICE formula builder: transformed `confounders_tv` (e.g. `poly(L, 2)`) | ✅ | test-ice.R |
| ICE formula builder: `substitute_vars_in_term()` helper | ✅ | test-ice.R |
| ICE self-consistency: bare vs poly() vs ns() vs I() on linear DGP | ✅ | test-ice.R |
| ICE self-consistency: log() transform on linear DGP | ✅ | test-ice.R |
| ICE nonlinear DGP: ns() handles sin() confounding | ✅ | test-ice.R |
| ICE lmtp cross-check: ns() on nonlinear DGP vs lmtp_sdr | ✅ | test-ice.R |
| ICE lmtp cross-check: poly() on linear DGP vs lmtp_sdr | ✅ | test-ice.R |
| Grace period / visit-process interventions | → Phase 22 | — |
| Stratified ICE | → Phase 22 | — |
| Multinomial outcomes | → Phase 23 | — |
| Ordinal outcomes | → Phase 23 | — |

### Phase 16 — AIPW / doubly-robust estimator
Composes Phase 2 gcomp + Phase 4 IPW into the classical analytical doubly-robust estimator: $\hat\psi_{\mathrm{AIPW}}(a) = E[\hat{m}(a, L)] + E[W \cdot (Y - \hat{m}(A, L))]$. Stacked EE sandwich (outcome block + propensity block + plug-in); distinct from `lmtp` (TMLE/SDR with ML + cross-fitting). Double-robustness + efficiency tests; ICE-AIPW longitudinal extension (Bang & Robins 2005).

| Treatment | Outcome | Intervention | Estimand | Contrast | Variance | Status | Test file |
|---|---|---|---|---|---|---|---|
| binary | gaussian | static | ATE | diff | sandwich | ✅ | test-aipw.R |
| binary | gaussian | static | ATE | diff | bootstrap | ✅ | test-aipw.R |
| binary | gaussian | static | ATT | diff | sandwich | ✅ | test-aipw.R |
| binary | gaussian | static | ATC | diff | sandwich | ✅ | test-aipw.R |
| binary | gaussian | static | ATE + EM (by=sex) | diff | sandwich | ✅ | test-aipw.R |
| binary | binomial | static | ATE | diff/ratio/OR | sandwich | ✅ | test-aipw.R |
| binary | gaussian | dynamic | ATE | diff | sandwich | ✅ | test-aipw.R |
| binary | gaussian | ipsi(δ) | ATE | diff | sandwich | ✅ | test-aipw.R |
| continuous | gaussian | shift | ATE | diff | sandwich | ✅ | test-aipw.R |
| continuous | gaussian | scale_by | ATE | diff | sandwich | ✅ | test-aipw.R |
| categorical | gaussian | static | ATE | diff | sandwich | ✅ | test-aipw.R |
| count (pois) | gaussian | shift | ATE | diff | sandwich | ✅ | test-aipw.R |
| binary | gaussian | DR: wrong outcome | ATE | diff | sandwich | ✅ | test-aipw.R |
| binary | gaussian | DR: wrong propensity | ATE | diff | sandwich | ✅ | test-aipw.R |
| binary | gaussian | efficiency: SE ≤ gcomp & IPW | ATE | diff | sandwich | ✅ | test-aipw.R |
| binary | gaussian | delicatessen cross-check (stacked EE) | ATE | diff | sandwich | ✅ point + SE match | test-aipw.R |
| continuous | gaussian | delicatessen cross-check (shift δ=1) | ATE | diff | sandwich | ✅ point + SE match | test-aipw.R |
| binary | gaussian | longitudinal static (always vs never) | ATE | diff | sandwich | ✅ ATE ≈ 5 | test-aipw.R |
| binary | gaussian | longitudinal static (always vs never) | ATE | diff | bootstrap | ✅ ATE ≈ 5 | test-aipw.R |
| continuous | gaussian | longitudinal shift | ATE | diff | sandwich | ✅ vs ICE + long-IPW | test-aipw.R |
| continuous | gaussian | longitudinal shift | ATE | diff | bootstrap | ✅ | test-aipw.R |
| binary | gaussian | longitudinal dynamic | ATE | diff | sandwich | ✅ | test-aipw.R |
| continuous | gaussian | longitudinal scale_by | ATE | diff | sandwich | ✅ | test-aipw.R |
| binary | gaussian | longitudinal DR: wrong outcome | ATE | diff | bootstrap | ✅ ATE ≈ 5 | test-aipw.R |
| binary | gaussian | longitudinal DR: wrong propensity | ATE | diff | bootstrap | ✅ ATE ≈ 5 | test-aipw.R |
| binary | gaussian | longitudinal static + EM (by=sex) | ATE | diff | sandwich | ✅ per-sex ATE ≈ 5,8 | test-aipw.R |
| binary | gaussian | longitudinal cross-method ICE vs AIPW vs long-IPW | ATE | diff | — | ✅ | test-aipw.R |
| binary | gaussian | longitudinal 3-period AIPW vs IPW | ATE | diff | sandwich | ✅ | test-aipw.R |
| — | — | longitudinal sandwich vs bootstrap SE agreement | ATE | — | — | ✅ ratio ∈ (0.5, 2) | test-aipw.R |
| binary | gaussian | longitudinal lmtp cross-check (binary static) | ATE | diff | — | ✅ vs lmtp_sdr | test-aipw.R |
| continuous | gaussian | longitudinal lmtp cross-check (continuous shift) | ATE | diff | — | ✅ vs lmtp_sdr | test-aipw.R |
| binary | binomial | longitudinal static (always vs never) | ATE | diff | sandwich | ✅ | test-aipw.R |
| binary | binomial | longitudinal lmtp cross-check (binary outcome) | ATE | diff | — | ✅ vs lmtp_sdr | test-aipw.R |

**Chunk 16m — Multivariate point AIPW**

| Feature | Status | Test |
|---|---|---|
| mv AIPW: binary × binary static recovers truth (ATE ≈ 2.5) | ✅ | test-aipw.R |
| mv AIPW: cross-checks gcomp and IPW (all agree within 0.15) | ✅ | test-aipw.R |
| mv AIPW: bootstrap parity with sandwich (within 30%) | ✅ | test-aipw.R |
| mv AIPW: continuous × continuous shift (truth ≈ −0.9) | ✅ | test-aipw.R |
| mv AIPW: DR — wrong outcome, correct propensity (ATE ≈ 2.5) | ✅ | test-aipw.R |
| mv AIPW: DR — correct outcome, wrong propensity (ATE ≈ 2.5) | ✅ | test-aipw.R |
| mv AIPW: binary × binary binomial outcome (finite RD) | ✅ | test-aipw.R |
| mv AIPW: stabilize = "marginal" matches unstabilized | ✅ | test-aipw.R |
| AIPW rejects stabilize for univariate treatment | ✅ | test-aipw.R |
| AIPW warns when propensity_model_fn not specified | ✅ | test-aipw.R (implicit) |

**Rejections (point, same as IPW):** static/threshold/dynamic on continuous → Dirac rejection ✅; stochastic (without `density`) → rejected ✅; stochastic (with `density`) → Phase 24.

**Rejections (longitudinal):** multivariate → deferred ✅; ATT/ATC → rejected ✅.

**Chunk 16o — AIPW + IPCW (triply-weighted DR)**

| Feature | Status | Test |
|---|---|---|
| AIPW + IPCW recovers ATE on DGP-M2b | ✅ | test-aipw-ipcw.R |
| AIPW + IPCW: sandwich SE finite, CI covers truth | ✅ | test-aipw-ipcw.R |
| AIPW + IPCW reduces bias vs naive complete-case AIPW | ✅ | test-aipw-ipcw.R |
| AIPW + IPCW: sandwich vs bootstrap SE ratio in (0.5, 2.0) | ✅ | test-aipw-ipcw.R |
| AIPW + IPCW DR: wrong outcome model still recovers truth | ✅ | test-aipw-ipcw.R |
| AIPW + IPCW + external weights composes correctly | ✅ | test-aipw-ipcw.R |
| AIPW + IPCW stashes correct details | ✅ | test-aipw-ipcw.R |

**Chunk 16p — AIPW test coverage parity**

| Feature | Status | Test |
|---|---|---|
| AIPW × gamma(log) × diff | ✅ | test-aipw.R |
| AIPW × gamma(log) × ratio | ✅ | test-aipw.R |
| AIPW × quasibinomial × diff | ✅ | test-aipw.R |
| AIPW × quasibinomial × ratio | ✅ | test-aipw.R |
| AIPW × negbin × diff | ✅ | test-aipw.R |
| AIPW × negbin × ratio | ✅ | test-aipw.R |
| AIPW × beta (betareg) × diff × bootstrap | ✅ | test-aipw.R |
| AIPW × beta (betareg) × ratio × bootstrap | ✅ | test-aipw.R |
| AIPW × external weights × sandwich | ✅ | test-aipw.R |
| AIPW × external weights × bootstrap SE agreement | ✅ | test-aipw.R |
| AIPW × svydesign auto-extract | ✅ | test-aipw.R |
| AIPW × svydesign equivalence to manual | ✅ | test-aipw.R |
| AIPW × cluster-robust SE inflation | ✅ | test-aipw.R |
| AIPW × cluster × truth recovery | ✅ | test-aipw.R |
| AIPW × cluster × by(sex) | ✅ | test-aipw.R |
| AIPW × MCAR outcome NAs × censoring × sandwich | ✅ | test-aipw.R |
| AIPW × MCAR outcome NAs × bootstrap SE agreement | ✅ | test-aipw.R |
| AIPW × covariate NAs → error (require pre-processing) | ✅ | test-aipw.R |
| AIPW × external weights + MCAR outcome NAs | ✅ | test-aipw.R |
| AIPW × subset(L > 0) restricts population | ✅ | test-aipw.R |
| AIPW × subset + by composition | ✅ | test-aipw.R |
| AIPW × subset + ATT | ✅ | test-aipw.R |
| AIPW × GAM outcome × nonlinear DGP × sandwich | ✅ | test-aipw.R |
| AIPW × GAM propensity DR (wrong outcome, correct GAM) | ✅ | test-aipw.R |
| AIPW × GAM outcome+propensity × sandwich vs bootstrap SE | ✅ | test-aipw.R |

**Chunk 16q — model_fn / propensity_model_fn default warnings audit**

| Feature | Status | Test |
|---|---|---|
| `model_fn` default warns for gcomp | ✅ | test-checks.R |
| `model_fn` default warns for IPW | ✅ | test-checks.R |
| `propensity_model_fn` default warns for IPW | ✅ | test-checks.R |
| `propensity_model_fn` default warns for AIPW | ✅ | test-checks.R |
| No `model_fn` warning for matching | ✅ | test-checks.R |
| Explicit `model_fn` suppresses warning | ✅ | test-checks.R |

### Phase 17 — Transportability / Generalizability
Sampling model `P(S=1 | L)` + sampling-odds weights multiply into IPW Hájek MSM; gcomp / IPW / AIPW transport paths; stacked EE extends with sampling-model block. References: Dahabreh et al. 2020; Westreich et al. 2017. Cross-check against `transport` / `transported` R packages.

**Chunk 17a — sampling model primitive**

| Feature | Status | Test |
|---|---|---|
| `fit_sampling_model()`: fit P(S=1\|L) via logistic GLM | ✅ | test-sampling-model.R |
| `causat(target = "S")`: stores sampling model in `fit$details` | ✅ | test-sampling-model.R |
| `causat(target_subset = "all")`: generalizability mode | ✅ | test-sampling-model.R |
| Validation: NA in S → error | ✅ | test-sampling-model.R |
| Validation: non-binary S → error | ✅ | test-sampling-model.R |
| Validation: degenerate S (all 0 or all 1) → error | ✅ | test-sampling-model.R |
| Validation: matching + target → error | ✅ | test-sampling-model.R |
| Validation: extreme selection rate → warning | ✅ | test-sampling-model.R |
| Target column preserved in `fit$data` | ✅ | test-sampling-model.R |

**Chunk 17b — gcomp transport**

| Feature | Status | Test |
|---|---|---|
| `causat(target = "S")` with gcomp restricts `fit_rows` to S=1 | ✅ | test-gcomp-transport.R |
| `target_subset = "target"`: standardize over S=0 rows (transportability) | ✅ | test-gcomp-transport.R |
| `target_subset = "all"`: standardize over all rows (generalizability) | ✅ | test-gcomp-transport.R |
| Truth-based test (transportability): ATE ≈ 3 + E[L\|S=0] | ✅ | test-gcomp-transport.R |
| Truth-based test (generalizability): ATE ≈ 3 + E[L] ≈ 3 | ✅ | test-gcomp-transport.R |
| Transport corrects study-population bias | ✅ | test-gcomp-transport.R |
| Sandwich SE plausible (ratio to bootstrap in (0.5, 2)) | ✅ | test-gcomp-transport.R |
| Target rows with NA outcome/treatment handled | ✅ | test-gcomp-transport.R |
| Validation: ATT/ATC with `target` → error | ✅ | test-gcomp-transport.R |
| Bootstrap refits outcome model on S=1 rows per replicate | ✅ | test-gcomp-transport.R |
| `target = NULL` (non-transport): behaviour unchanged | ✅ | test-gcomp-transport.R |
| Treatment-covariate interactions stripped from sampling model formula | ✅ | test-gcomp-transport.R |

**Chunk 17c — IPW transport**

| Feature | Status | Test file |
|---|---|---|
| `causat(target = "S", estimator = "ipw")` restricts `fit_rows` to S=1 | ✅ | test-ipw-transport.R |
| `target_subset` stored on fit | ✅ | test-ipw-transport.R |
| Sampling × treatment weight product in MSM refit | ✅ | test-ipw-transport.R |
| Truth-based test (transportability): ATE ≈ 3 + E[L\|S=0] | ✅ | test-ipw-transport.R |
| Truth-based test (generalizability): ATE ≈ 3 + E[L] | ✅ | test-ipw-transport.R |
| Transport corrects study-population bias vs. naive IPW | ✅ | test-ipw-transport.R |
| Sandwich SE plausible (ratio to bootstrap in (0.6, 1.8)) | ✅ | test-ipw-transport.R |
| Target rows with NA outcome/treatment handled | ✅ | test-ipw-transport.R |
| Bootstrap refits sampling + propensity models per replicate | ✅ | test-ipw-transport.R |
| `target = NULL` (non-transport): behaviour unchanged | ✅ | test-ipw-transport.R |

**Chunk 17d — Bootstrap (standalone)**

| Feature | Status | Test |
|---|---|---|
| gcomp bootstrap (transportability): point estimate near truth | ✅ | test-transport-bootstrap.R |
| gcomp bootstrap (generalizability): point estimate near truth | ✅ | test-transport-bootstrap.R |
| gcomp bootstrap CI brackets truth | ✅ | test-transport-bootstrap.R |
| gcomp bootstrap point estimate equals sandwich estimate | ✅ | test-transport-bootstrap.R |
| IPW bootstrap (transportability): point estimate near truth | ✅ | test-transport-bootstrap.R |
| IPW bootstrap (generalizability): point estimate near truth | ✅ | test-transport-bootstrap.R |
| IPW bootstrap CI brackets truth | ✅ | test-transport-bootstrap.R |
| IPW bootstrap: sampling model refitted per replicate (`gamma_hat` differs) | ✅ | test-transport-bootstrap.R |
| Cross-estimator: gcomp vs IPW bootstrap transport estimates agree | ✅ | test-transport-bootstrap.R |

**Chunk 17e — AIPW transport**

| Feature | Status | Test |
|---|---|---|
| AIPW transport: fit_rows restricted to study (S=1) rows | ✅ | test-aipw-transport.R |
| AIPW transport: target_subset stored on fit | ✅ | test-aipw-transport.R |
| AIPW transport (transportability): recovers target ATE | ✅ | test-aipw-transport.R |
| AIPW transport (generalizability): recovers marginal ATE | ✅ | test-aipw-transport.R |
| AIPW transport corrects study bias vs naive | ✅ | test-aipw-transport.R |
| AIPW transport 2-of-3 DR: wrong outcome model | ✅ | test-aipw-transport.R |
| AIPW transport 2-of-3 DR: wrong propensity model | ✅ | test-aipw-transport.R |
| AIPW transport 2-of-3 DR: wrong sampling model | ✅ | test-aipw-transport.R |
| AIPW transport: sandwich SE plausible (ratio to bootstrap) | ✅ | test-aipw-transport.R |
| AIPW transport: bootstrap point estimate near truth | ✅ | test-aipw-transport.R |
| AIPW transport: continuous treatment (binary, static) | ✅ | test-aipw-transport.R |
| AIPW transport: dynamic intervention | ✅ | test-aipw-transport.R |
| AIPW transport: ipsi intervention (static fallback) | ✅ | test-aipw-transport.R |
| AIPW transport: binary treatment stability (different seeds) | ✅ | test-aipw-transport.R |
| AIPW transport: binomial outcome (diff/ratio/OR) | ✅ | test-aipw-transport.R |
| AIPW transport: categorical treatment | ✅ | test-aipw-transport.R |
| AIPW transport: count (Poisson) treatment (fit only) | ✅ | test-aipw-transport.R |
| AIPW transport: effect modification (by=sex) | ✅ | test-aipw-transport.R |
| AIPW transport: cross-estimator agreement (gcomp, IPW, AIPW) | ✅ | test-aipw-transport.R |
| AIPW transport: efficiency (SE_AIPW ≤ SE_gcomp ∧ SE_IPW) | ✅ | test-aipw-transport.R |
| AIPW transport: target rows with NA outcome/treatment handled | ✅ | test-aipw-transport.R |
| AIPW without transport unaffected (target=NULL) | ✅ | test-aipw-transport.R |

**Chunk 17f — Sampling-model diagnostics**

| Feature | Status | Test |
|---|---|---|
| `diagnose()` sampling panel for gcomp transport | ✅ | test-diagnose.R |
| `diagnose()` sampling panel for IPW transport | ✅ | test-diagnose.R |
| `diagnose()` sampling panel for AIPW transport | ✅ | test-diagnose.R |
| `diagnose()` sampling panel NULL when target = NULL | ✅ | test-diagnose.R |
| Sampling weight plausibility (ESS, quantiles, extreme count) | ✅ | test-diagnose.R |
| Per-intervention panels carry shared sampling panel | ✅ | test-diagnose.R |
| `print(diagnose(...))` renders sampling panel | ✅ | test-diagnose.R |

**Chunk 17g — Sampling model predictor-set validation**

| Feature | Status | Test |
|---|---|---|
| `fit_sampling_model()` warns when baseline covariate only in treatment interaction | ✅ | test-sampling-model.R |
| No warning when covariate has a main effect alongside interaction | ✅ | test-sampling-model.R |
| `causat()` propagates warning end-to-end | ✅ | test-sampling-model.R |

**Chunk 17h — External cross-check (TransportHealth oracle)**

| Feature | Status | Test |
|---|---|---|
| gcomp transport TATE matches `TransportHealth::transportGC()` (tolerance 1e-4) | ✅ | test-gcomp-transport.R |
| IPW transport TATE matches `TransportHealth::transportIP()` (tolerance 1e-4) | ✅ | test-ipw-transport.R |

**Chunk 17i — Longitudinal IPW × transport**

| Feature | Status | Test |
|---|---|---|
| `causat(target="S", estimator="ipw", id=, time=)`: sampling model fit on first-period rows only | ✅ | test-longitudinal-transport.R |
| Sampling odds weight broadcast to all person-period rows for each individual | ✅ | test-longitudinal-transport.R |
| Sampling weight × per-period treatment DR weight in longitudinal MSM | ✅ | test-longitudinal-transport.R |
| `target_subset="target"` (transportability): truth-based test ATE ≈ 4 + E[L\|S=0] ≈ 3.67 | ✅ | test-longitudinal-transport.R |
| `target_subset="all"` (generalizability): truth-based test ATE = 4 | ✅ | test-longitudinal-transport.R |
| Transport corrects study-population bias (naive longitudinal IPW disagrees by >0.25) | ✅ | test-longitudinal-transport.R |
| Long IPW transport agrees with long gcomp transport (linear DGP; diff < 0.15) | ✅ | test-longitudinal-transport.R |
| Sandwich SE plausible (ratio to bootstrap SE in (0.4, 2.5)) | ✅ | test-longitudinal-transport.R |
| Bootstrap refits sampling model per replicate; original fit unchanged | ✅ | test-longitudinal-transport.R |
| `stabilize="marginal"` + transport: finite ATE and SE | ✅ | test-longitudinal-transport.R |
| `target=NULL` longitudinal IPW (non-transport): behaviour unchanged | ✅ | test-longitudinal-transport.R |
| ICE gcomp + longitudinal + transport: NA lag columns filled for target individuals under static interventions | ✅ | test-longitudinal-transport.R |

**Chunk 17j — Multivariate IPW × transport**

| Feature | Status | Test |
|---|---|---|
| `fit_rows` restricted to S=1 study rows | ✅ | test-multivariate-ipw-transport.R |
| mv IPW transport (transportability): recovers target ATE ≈ 3.5 + E[L\|S=0] | ✅ | test-multivariate-ipw-transport.R |
| Transport corrects study-population bias (naive disagrees by >0.1) | ✅ | test-multivariate-ipw-transport.R |
| mv IPW transport (generalizability): recovers ATE = 3.5 | ✅ | test-multivariate-ipw-transport.R |
| Sandwich SE plausible vs bootstrap (ratio in (0.4, 2.5)) | ✅ | test-multivariate-ipw-transport.R |
| Bootstrap point estimate near truth (within 0.25) | ✅ | test-multivariate-ipw-transport.R |
| `target=NULL`: behaviour unchanged | ✅ | test-multivariate-ipw-transport.R |
| `stabilize="marginal"` + transport: recovers target ATE, finite SE | ✅ | test-multivariate-ipw-transport.R |

**Implementation note.** The only source change was removing the early `return()` in the multivariate branch of `variance_if_ipw()` and applying `compute_ipw_sampling_correction()` when `fit$details$transport == TRUE`. Weight multiplication and bootstrap already handled multivariate transport correctly; only the sandwich gamma-block correction was missing.

Remaining chunks (17k–17l): IPCW × transport, MTP + transportability. Chunk 17m (documentation) pending. All ❌.

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
