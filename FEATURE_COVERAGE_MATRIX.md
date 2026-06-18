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

## Test-file organization

A few large test files are split into thematic siblings purely for CI test-file
parallelism; the estimand coverage and the rows below are unchanged.

- `test-aipw.R` (point AIPW) + `test-aipw-longitudinal.R` (ICE-AIPW + multivariate longitudinal) + `test-aipw-lmtp-oracle.R` (`lmtp` SDR cross-checks) + `test-aipw-point-extra.R` (non-Gaussian families, external/survey/cluster weights, missing data, subset, GAM)
- `test-longitudinal-ipw.R` (univariate) + `test-longitudinal-ipw-mv.R` (multivariate)
- `test-ipcw-variance.R` (sandwich-vs-bootstrap + jacobian) + `test-ipcw-variance-coverage.R` (95% CI coverage)

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

### Categorical (multinomial) outcome — point g-comp (Phase 23a-1 / 23a-2a / 23a-2b)

Single K-level factor outcome via `model_fn = nnet::multinom`. Estimand is the
K-vector `P(Y = k | do(A = a))` per intervention; `estimates` / `contrasts` gain
a `class` column. Variance is **bootstrap** (23a-1) or the analytic per-class IF
**sandwich** (23a-2a complete-case + 23a-2b survey/external-weighted;
`variance_if_gcomp_multinom()`). Oracles:
large-n softmax g-computation truth + exact `marginaleffects::avg_predictions()`
on causatr's own multinom fit (point parity ~1e-15); the sandwich SE is pinned to
a Python M-estimation stack (`fixtures/python/multinom_gcomp_sandwich.py` and the
weighted `multinom_gcomp_weighted_sandwich.py`, ~1e-7 / ~1e-10). All rows
`test-gcomp-categorical-outcome.R`.

| Trt | Outcome | Intervention | Estimand | Contrast | Variance | Extra | Status |
|---|---|---|---|---|---|---|---|
| bin | 3-class | static | ATE | diff | boot | — | ✅ (truth + marginaleffects) |
| bin | 3-class | static | ATE | ratio | boot | — | ✅ (exact fn of means) |
| bin | 3-class | static | ATE | or | boot | — | ✅ (exact fn of means) |
| cont | 3-class | shift | ATE | diff | boot | — | ✅ (truth + marginaleffects) |
| cat | 3-class | static | ATE | diff | boot | — | ✅ (marginaleffects) |
| bin | 4-class | static | ATE | diff | boot | — | ✅ (schema + marginaleffects) |
| bin | 3-class | static | ATT | diff | boot | — | ✅ (treated-standardised oracle) |
| bin | 3-class | static | by(G) | diff | boot | — | ✅ (per-stratum oracle) |
| bin | 3-class | static | subset | diff | boot | — | ✅ (subset oracle) |
| bin | 3-class | static | ATE | diff | boot | survey wts | ✅ (weighted oracle) |
| bin | 3-class | static | ATE | diff | boot | IPCW | ✅ (MAR de-biasing → truth) |
| cont | 3-class | shift(×3) | ATE | diff | boot | ≥3 interventions | ✅ |
| bin | 3-class | static | ATE | diff | boot | spline confounders | ✅ |

S3 layer (coef / confint / tidy / glance / print / plot per-class) ✅.
Scalar byte-identical guard ✅. Bootstrap reproducibility + SE-vs-delta sanity ✅.

#### Analytic sandwich (Phase 23a-2a, complete-case)

| Trt | Outcome | Intervention | Estimand | Contrast | Extra | Status |
|---|---|---|---|---|---|---|
| bin | 3-class | static | ATE | diff / ratio / or | — | ✅ delicatessen M-est stack ~1e-7 / ~1e-10 |
| bin | 3-class | static | ATE | diff | — | ✅ sandwich vs bootstrap (≤10%) |
| cont | 3-class | shift | ATE | diff | — | ✅ marginaleffects point + boot SE parity |
| cat | 3-class | static | ATE | diff | — | ✅ sandwich vs bootstrap |
| bin | 4-class | static | ATE | diff | — | ✅ schema + boot SE parity |
| bin | 3-class | static | ATT | diff | — | ✅ treated-standardised point + boot SE |
| bin | 3-class | static | by(G) / subset | diff | spline confounders | ✅ per-stratum / subset sandwich + boot SE |
| bin | 3-class | static | ATE | diff | S3 (boot_t NULL) | ✅ print / summary / tidy / coef / confint |

Sandwich SE is the full M-estimation IF (Channel 1 empirical-distribution +
Channel 2 nuisance); `marginaleffects` omits Channel 1, so the tight SE oracle is
the Python `delicatessen`-style stack.

#### Analytic sandwich — survey / external weights (Phase 23a-2b)

The outcome `nnet::multinom` is a weighted MLE, so the sandwich generalises by
weighting the stacked bread `H = Σ wᵢ (diag(pᵢ) − pᵢpᵢᵀ) ⊗ XᵢXᵢᵀ`, the score
residual (`r_score = wᵢ`), Channel 1 (`(n/Σw) wᵢ (p_{k,i} − μ)`), and the softmax
marginal-mean gradient (`(1/Σw) Σ wᵢ p_{k,i}(δ − p_{l+1,i}) Xᵢ*`). With `w ≡ 1`
all four collapse to the 23a-2a complete-case formulas byte-for-byte.

| Trt | Outcome | Intervention | Estimand | Contrast | Extra | Status |
|---|---|---|---|---|---|---|
| bin | 3-class | static | ATE | diff / ratio / or | survey wts | ✅ weighted delicatessen M-est stack ~1e-7 |
| bin | 3-class | static | ATE | diff | survey wts | ✅ weighted sandwich vs weighted bootstrap (≤10%) |
| bin | 3-class | static | ATE | diff | `weights = NULL` ≡ ones | ✅ byte-identical to complete-case vcov (1e-12) |

**IPCW** sandwich stays on the bootstrap (23a-2c).

Rejections (all ✅, `test-gcomp-categorical-outcome.R`): `snm` ⛔
(`causatr_snm_categorical_outcome`); `ipw` / `aipw` / `matching` / longitudinal /
transport ⛔ (`causatr_categorical_outcome_unsupported`); `ci_method = "sandwich"`
with IPCW ⛔ (`causatr_categorical_outcome_sandwich`, until 23a-2c); stochastic
intervention ⛔ (`causatr_categorical_outcome_unsupported`, until 23a-4b).

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
| **MV binary** | gauss | static (both1 vs both0) | 2 | sandwich | — | ✅ SE ratio ~1 | test-longitudinal-ipw.R |
| **MV binary** | gauss | static (both1 vs both0) | 2 | bootstrap | — | ✅ vs sandwich | test-longitudinal-ipw.R |
| **MV cont** | gauss | shift + NULL (natural course) | 2 | sandwich | — | ✅ SE ratio ~1 | test-longitudinal-ipw.R |
| **MV cont** | gauss | shift + NULL (natural course) | 2 | bootstrap | — | ✅ vs sandwich | test-longitudinal-ipw.R |
| **MV binary** | gauss | static vs NULL (natural course) | 2 | sandwich | — | ✅ finite | test-longitudinal-ipw.R |
| **MV** (numerator structure) | — | `stabilize="marginal"` per-period × per-component chain numerators (L dropped, lags + prior components kept) | 2 | — | — | ✅ formula check | test-longitudinal-ipw.R |
| **MV binary** | gauss | static (both1 vs both0) + `stabilize="marginal"` | 2 | sandwich | — | ✅ identical to unstabilized | test-longitudinal-ipw.R |
| **MV cont** | gauss | shift + NULL + `stabilize="marginal"` | 2 | sandwich + bootstrap | — | ✅ SE ratio ~1 | test-longitudinal-ipw.R |
| **MV binary** | gauss | static (both1 vs both0) + EM (`A1:sex` + `A2:sex`), by(sex) | 2 | sandwich | — | ✅ ATE\|sex=0 = 9, ATE\|sex=1 = 15 | test-longitudinal-ipw.R |
| **MV binary** | gauss | static + EM, by(sex) | 2 | sandwich | — | ✅ vs ICE g-comp | test-longitudinal-ipw.R |
| **MV** (EM structure) | — | each per-period × per-component propensity strips `A1:sex`/`A2:sex` (keeps `sex` as confounder); MSM expands to `Y ~ 1 + sex` | 2 | — | — | ✅ formula check | test-longitudinal-ipw.R |
| **MV binary** | gauss | static + EM + `stabilize="marginal"`, by(sex) | 2 | sandwich | — | ✅ identical to unstabilized | test-longitudinal-ipw.R |
| bin | gauss | `ipsi($\delta$=2)` (per-period Kennedy product) | 2 | sandwich | — | ✅ +oracle (manual per-period Kennedy product) | test-longitudinal-ipw.R |
| bin | gauss | `ipsi($\delta$=2)` | 2 | bootstrap | — | ✅ vs sandwich (within 15%) | test-longitudinal-ipw.R |
| bin | gauss | `ipsi($\delta$=2)` vs `ipsi($\delta$=1)` (difference) | 2 | sandwich | — | ✅ +oracle | test-longitudinal-ipw.R |

**Per-period IPSI (Phase 20):** univariate `ipsi(delta)` extends Kennedy's (2019) closed-form weight to a per-period product $W_i = \prod_t (\delta A_{t} + (1 - A_{t})) / (\delta p_{t} + (1 - p_{t}))$. No new density-evaluation path — each period reuses `compute_density_ratio_weights()`'s IPSI branch; the stacked sandwich reuses `make_weight_fn()`'s IPSI sub-closure per period.

**Known limitation (Phase 6, Robins 2000):** under longitudinal IPW the modifier MUST be a **baseline** covariate. A time-varying modifier in an MSM conditions on a post-treatment variable, biasing the estimand via mediator + collider paths. Not enforced at runtime (time-varying status isn't inferable from data); doc-level constraint only. The scientifically correct tool for time-varying effect modification is a structural nested model (Phase 18).

Sequential positivity warning (`causatr_longitudinal_seq_positivity`): fires when any per-period weight max exceeds 100 (default threshold); silent below threshold. Tested directly on `warn_seq_positivity()`.

Rejections (all ✅ tested):
- multivariate `ipsi()` (any component) $\to$ test-longitudinal-ipw.R (`causatr_longitudinal_ipsi_multivariate`; Kennedy's closed form is binary-univariate)
- `ipsi()` + `stabilize = "marginal"` $\to$ test-longitudinal-ipw.R (`causatr_longitudinal_ipsi_stabilize`; IPSI weight has no separate marginal numerator)
- `numerator =` without `stabilize = "marginal"` $\to$ test-longitudinal-ipw.R (`causatr_longitudinal_numerator_without_stabilize`)
- bare treatment in confounders (`~ L + A`) $\to$ test-longitudinal-ipw.R (`causatr_bare_treatment_in_confounders`)
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
| bin | poisson | static | 2 | sandwich | — | ✅ vs analytical truth + Python M-estimation | test-ice-outcome-types.R |
| cont | gamma | shift | 2 | sandwich | — | ✅ finite + lmtp cross-check | test-ice-outcome-types.R |
| bin | glm.nb | static | 2 | sandwich | — | ✅ vs analytical truth + lmtp cross-check | test-ice-outcome-types.R |
| bin | betareg | static | 2 | boot | — | 🟡 smoke (betareg coef structure not sandwich-compatible) | test-ice-outcome-types.R |

### Stratified ICE (`stratified = "G"`, Phase 22a)

Per-baseline-stratum outcome models at each backward step. **Both** the ID-cluster
bootstrap and the analytic per-stratum × per-time stacked-EE sandwich are supported
(`R/variance_if_ice_stratified.R`; derivation in `PHASE_22_ICE_ENHANCEMENTS.md` §1.6).

| Trt | Outcome | Intervention | Periods | Variance | Wt | Status | Test |
|---|---|---|---|---|---|---|---|
| bin | gauss | static/shift | 2 | (point) | — | ✅ exact vs per-stratum pooled | test-ice-stratified.R |
| bin | gauss | static | 2 | sandwich | — | ✅ exact vs pooled sandwich at S=1 | test-ice-stratified.R |
| bin | gauss | static | 2 | sandwich | — | ✅ vs delicatessen ~1e-7 (mean/SE) | test-ice-stratified.R |
| bin | binom | static | 2 | sandwich | — | ✅ vs delicatessen ~1e-7 (mean/SE) | test-ice-stratified.R |
| cont | gauss | shift | 2 | sandwich | — | ✅ vs ID-cluster bootstrap | test-ice-stratified.R |
| bin | gauss | static | 2 | sandwich | ext | ✅ vs ID-cluster bootstrap | test-ice-stratified.R |
| bin | binom | static | 2 | boot | — | ✅ truth + lmtp; boot vs sandwich/delicatessen | test-ice-stratified.R |
| cont | gauss | shift | 2 | boot | — | ✅ | test-ice-stratified.R |
| bin | gauss | dynamic (cov) | 2 | (point) | — | ✅ | test-ice-stratified.R |
| bin | gauss | static | 2 | (point) | ext / cens | ✅ | test-ice-stratified.R |
| factor | gauss | static | 2 | sandwich | — | ✅ byte-identical to integer coding | test-ice-stratified.R |

Rejections (all ✅ tested): ipsi $\to$ test-ice.R, ATT/ATC $\to$ test-ice.R.
Stratified rejections (all ✅ tested, `test-ice-stratified.R`): non-ICE estimator / point treatment
$\to$ `causatr_stratified_not_ice`; missing column $\to$ `causatr_stratified_not_found`;
time-varying column $\to$ `causatr_stratified_not_baseline`; continuous column $\to$
`causatr_stratified_too_many`.

### Flexible-treatment ICE term (`treatment_form = ~ factor(A)` / `~ ns(A)`, Phase 22b-5)

By default the treatment enters every per-step ICE outcome model as a bare numeric main effect, so a
non-monotone or curved counterfactual dose-response is misspecified. `treatment_form` lets the treatment
enter via a transformed term while the intervention still sets the *numeric* treatment column; lag terms
are expanded automatically (`factor(A)` $\to$ `factor(lag1_A)`). **Longitudinal g-computation (ICE) only**,
including the natural-history MTPs above. Validated by exact analytic-contrast truths (DGPs with exogenous
covariates, so the contrast cancels the model-independent baseline term).

| Trt | Outcome | `treatment_form` | Intervention | Periods | Variance | Status | Test |
|---|---|---|---|---|---|---|---|
| categorical {0,1,2} | gauss | `~ factor(A)` | static(1) vs static(0) | 2 | (point) | ✅ recovers non-monotone truth (= 6); bare ~4.3 off | test-ice-treatment-form.R |
| cont | gauss | `~ splines::ns(A, df)` | shift(0.5) vs shift(0) | 2 | (point) | ✅ recovers curved truth (= −1.5); linear ~0.5 off | test-ice-treatment-form.R |
| categorical {0,1,2} | gauss | `~ factor(A)` | static | 2 | sandwich | ✅ vs ID-cluster bootstrap (~5%) | test-ice-treatment-form.R |
| any | — | `NULL` / `~ A` | — | 2 | (point) | ✅ byte-identical to bare numeric (1e-12) | test-ice-treatment-form.R |

Rejections (all ✅ tested, `test-ice-treatment-form.R`): non-ICE estimator / point treatment $\to$
`causatr_treatment_form_not_ice`; non-treatment variable / two-sided formula / non-formula $\to$
`causatr_treatment_form_bad`; intervention leaving the observed `factor(A)` support $\to$ native
`predict()` "new level" error (documented constraint; `ns(A)` extrapolates).

### Natural-history MTPs (G-LMTPs, `grace_period()` / `carry_forward()` / `cap_escalation()`, Phase 22b)

Modified treatment policies whose intervened value at time $t$ depends on the
**natural-value history of treatment** (delays, grace periods, carry-forward). The
standard ICE recursion is provably wrong here (it conditions on the *observed* lag);
the **augmented-data sequential regression** of Diaz, Williams, Morzywołek & Rudolph
(2026, arXiv:2605.24167) carries every natural-history sequence $\bar{s}_t$ as a label
through the backward recursion (`R/glmtp.R`, `R/glmtp_augment.R`,
`R/glmtp_interventions.R`). **Longitudinal g-computation only, discrete treatment.**
Variance: ID-cluster **bootstrap** and analytic **M-estimation sandwich**
(`R/variance_if_glmtp.R`, 22b-4).

| Trt | Outcome | Intervention | Periods | Variance | Wt | Status | Test |
|---|---|---|---|---|---|---|---|
| bin (absorbing) | gauss | `grace_period(1)` | 4 | (point) | — | ✅ vs forward-MC truth (Prop. 1) | test-glmtp.R |
| bin (absorbing) | binom | `grace_period(1)` | 5 | (point) | — | ✅ replicates Diaz et al. (2026) §6 truth | test-glmtp.R |
| bin | gauss | `grace_period(1)` vs naive `dynamic(lag1_A)` | 4 | (point) | — | ✅ engine necessity (naive off by ~56%) | test-glmtp.R |
| bin | gauss | `carry_forward()` | 4 | (point) | — | ✅ exact vs standard-ICE baseline regime (1e-8) | test-glmtp.R |
| bin | gauss | `grace_period(0)` | 4 | (point) | — | ✅ exact vs natural course (1e-9) | test-glmtp.R |
| bin | gauss | `grace_period(1)` | 4 | boot | — | ✅ CI covers forward-MC truth | test-glmtp.R |
| bin | gauss | `grace_period(1)` | 4 | boot | ext / cens | ✅ composes (finite point + SE) | test-glmtp.R |
| bin | gauss | `grace_period(1)` | 4 | sandwich | — | ✅ bootstrap parity (0.7% diff at n=1500) | test-glmtp.R |
| bin | gauss | `grace_period(1)` | 2 | sandwich | — | ✅ +oracle Python M-estimation (SE ~1e-4) | test-glmtp.R |
| bin | gauss | `grace_period(1)` vs `NULL` multi-arm | 4 | sandwich | — | ✅ finite contrasts + per-arm SEs | test-glmtp.R |
| bin | gauss | `grace_period(1)` | 4 | sandwich | ext | ✅ finite SE with external weights | test-glmtp.R |
| ordinal {0,1,2} | gauss | `carry_forward()` | 3 | (point) | — | ✅ exact vs baseline regime (1e-14) — `|A|^t` machinery | test-glmtp.R |
| ordinal {0,1,2} | gauss | `cap_escalation` | 3 | (point) | — | ✅ engine == independent hand-coded recursion (1e-9) | test-glmtp.R |
| ordinal {0,1,2} | gauss | `cap_escalation` + `treatment_form = ~ factor(A)` | 3 | (point) | — | ✅ vs forward-MC truth: comparative n=40000 (gap 0.0015 vs 0.034 bare) | test-glmtp.R |
| ordinal {0,1,2} | gauss | `cap_escalation(1)` + `~ factor(A)` | 3 | (point) | — | ✅ tight vs forward-MC truth, n=80000 (gap <0.0035; Tier-2) | test-glmtp.R |
| ordinal {0,1,2} | gauss | `cap_escalation(1)` + `~ factor(A)` | 3 | boot | — | ✅ CI covers forward-MC truth; `contrast()` == `glmtp_iterate()` | test-glmtp.R |
| ordinal {0,1,2} | gauss | `cap_escalation(1)` + `~ factor(A)` | 3 | sandwich | — | ✅ bootstrap parity (4.4% diff at n=1500) | test-glmtp.R |
| ordinal {0,1,2} | gauss | `cap_escalation()` ctor + rejections | 3 | (point) | — | ✅ arg validation + mixed/continuous/non-ICE classed | test-glmtp.R |
| bin | gauss | `grace_period(1)` + `~ factor(A)` | 4 | boot | — | ✅ composes (finite CI; == bare for binary) | test-glmtp.R |
| bin (absorbing) | poisson | `grace_period(1)` | 4 | (point) | — | ✅ vs forward-MC truth | test-glmtp.R |
| bin (absorbing) | poisson | `grace_period(1)` | 4 | sandwich | — | ✅ bootstrap parity | test-glmtp.R |
| bin (absorbing) | gamma | `grace_period(1)` | 4 | (point) | — | ✅ vs forward-MC truth | test-glmtp.R |

The engine handles **any single discrete treatment** ($|\mathcal{A}| \ge 2$); the ordinal rows pin the
`|A|^t` label machinery. Public policy constructors: `grace_period` (binary delay), `carry_forward`
(any-discrete LOCF), and `cap_escalation` (ordered-dose increase cap, 22b-6). The **flexible-treatment
ICE term** (`treatment_form =`, Phase 22b-5) lets the dose enter the per-step model as `factor(A)` /
`splines::ns(A)`, removing the bare-numeric misspecification of a kinked capped-dose response (the
cap-policy gap closes from $\approx0.034$ to $\approx0.0015$ vs the forward-MC truth) — which is what
makes `cap_escalation()` a consistent public estimator. The **augmented-data sandwich** (22b-4,
`R/variance_if_glmtp.R`) stacks per-(step, label) GLM scores + the estimand EE in the same
block-triangular M-estimation structure as the ICE chain; validated against a Python M-estimation
oracle (SE agreement ~1e-4) and against ID-cluster bootstrap parity (~0.7% at n=1500, ~4.4% for
`cap_escalation + factor(A)`). **Rejected by design** (`PHASE_22` §3.7):
**multivariate** (vector-treatment) G-LMTP (22b-7) — the paper develops the augmented-data
enumeration for a scalar discrete exposure only, and the joint-support blow-up
$|\mathcal{A}^{(1)}\times\cdots\times\mathcal{A}^{(K)}|^{\tau-1}$ is intractable.

Augment helpers (`glmtp_support`, `glmtp_enumerate_labels`, `glmtp_label_key`,
`glmtp_check_tractable`) — all ✅ unit-tested (`test-glmtp.R`).
Rejections (all ✅ tested, `test-glmtp.R`): non-ICE estimator / point treatment / transport /
stratified $\to$ `causatr_glmtp_not_ice`; mixing with a standard intervention $\to$
`causatr_glmtp_mixed`; continuous / factor treatment $\to$
`causatr_glmtp_continuous_trt`; multivariate (vector) treatment $\to$
`causatr_glmtp_multivariate`; history blow-up $\to$ `causatr_glmtp_too_many`.

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
| time-varying L (8%) | static | sandwich | ✅ (cascade restricted to model-complete rows; matches bootstrap/jackknife ~1.0) | test-missing-data.R |
| time-varying L, 3-period intermediate NA | static | sandwich | ✅ (two-step cascade) | test-missing-data.R |

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
| `boot_ci = "percentile"` (default) — means = replicate quantiles | ✅ | test-bootstrap-ci.R |
| `boot_ci = "percentile"` — diff/ratio/OR = per-replicate-contrast quantiles | ✅ | test-bootstrap-ci.R |
| `boot_ci = "normal"` — Wald from bootstrap SE (= legacy) | ✅ | test-bootstrap-ci.R |
| `boot_ci` point/SE/vcov invariant; cross-check vs `boot::boot.ci` | ✅ | test-bootstrap-ci.R |
| `boot_ci` honoured by `confint()` (+ override) and `tidy()` | ✅ | test-bootstrap-ci.R |
| `boot_ci` composes with IPW / AIPW / SNM (blip Path A) bootstrap | ✅ | test-bootstrap-ci.R |
| `boot_ci` invalid value rejected | ⛔ | test-bootstrap-ci.R |
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
| delicatessen extended: categorical AIPW + binary IPW + longitudinal IPW | ✅ | test-delicatessen-extended.R |
| delicatessen NHEFS: g-comp/IPW/AIPW/IPCW/EM SE match + Ross binary + Cole transport | ✅ | test-delicatessen-nhefs.R |
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
| `causat_mice()` multiple-imputation pooling (see MI section below) | ✅ | test-causat-mice.R |
| Cross-method EM triangulation (gcomp + IPW + matching) | ✅ | test-effect-modification.R |
| Per-component confounders (`confounders_outcome`, `confounders_treatment`, etc.) — routing, backward compat, validation, ground truth | ✅ | test-separate-confounders.R |
| Per-component confounders — Kang-Schafer S2/S3 DR with split formulas | ✅ | test-kang-schafer.R |
| Per-component confounders — delicatessen cross-check (synthetic + NHEFS) | ✅ | test-variance-reference.R |
| Per-component confounders — stdReg2 cross-check (synthetic DGPs) | ✅ | test-variance-reference.R |

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
Sequential density-ratio weights, cumulative product weights, stabilized weights, time-varying MSM. **All shipped.**

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
Scalar-outcome IPCW for MAR censoring: internal censoring model, stabilized IPCW weights, weighted fit, stacked EE sandwich extension for censoring model blocks. Point + ICE scalar final outcome. Survival-specific IPCW (per-period cumulative weights + hazard MSM) is owned by the separate survival package. **All shipped** (chunks 14a–14f).

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
| Grace period / delay / carry-forward (natural-history MTPs) — boot + sandwich | ✅ | test-glmtp.R |
| Stratified ICE | ✅ | test-ice-stratified.R |
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
| — | — | longitudinal sandwich vs bootstrap SE agreement | ATE | — | — | ✅ ratio ≈ 1.0 (±0.12) | test-aipw.R |
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

**Rejections (longitudinal AIPW):** multivariate longitudinal AIPW → supported (Phase 25) ✅ (see the MV longitudinal AIPW block below); ATT/ATC → rejected ✅; `stabilize = "marginal"` → rejected ✅ (`causatr_stabilize_longitudinal_aipw`; use `estimator = "ipw"` for stabilized longitudinal weights).

**Longitudinal AIPW sandwich — full stacked-EE M-estimation sandwich (balanced + UNBALANCED panels):** `ci_method = "sandwich"` builds `V = B⁻¹ M B⁻ᵀ / n` where `ψ(θ)` stacks the per-period propensity scores, the per-step ICE outcome / pseudo-outcome scores, and the marginal-mean equation; the numerical bread (`numDeriv`) captures every block-triangular cross-term (including the baseline-pseudo-regression block), so it is consistent under monotone dropout / censoring row-filter. This replaced the forward-cascade assembly that underestimated the SE by ~50% on unbalanced panels and aborted (`causatr_longitudinal_aipw_unbalanced_sandwich`, removed). The EE is faithful by construction (`~1e-11` summed-score vanishing); `causatr_longitudinal_aipw_sandwich_unfaithful` steers to bootstrap if it does not.

**Sandwich nuisance-model coverage (explicit):** GLM families (gaussian / binomial / poisson / Gamma / quasi-* / inverse-gaussian) + `MASS::glm.nb` + multinomial (`nnet::multinom`, categorical treatment) are scored analytically. **Penalised / non-likelihood fitters `mgcv::gam` and `betareg::betareg` are an explicit limitation** — they abort with `causatr_longitudinal_aipw_sandwich_model` and route to the bootstrap (their bread is not the score Jacobian; mirrors the longitudinal-ICE betareg bootstrap-only path).

| Treatment | Outcome | Intervention | Panel | Variance | Status | Test |
|---|---|---|---|---|---|---|
| binary | gaussian | static | balanced + unbalanced | sandwich | ✅ vs bootstrap + jackknife (ratio ≈ 1.0) | test-aipw-longitudinal.R |
| binary | gaussian | static | balanced + unbalanced | sandwich | ✅ vs `delicatessen` MEstimator (SE ~1e-7, point ~1e-13) | test-aipw-longitudinal.R |
| binary | gaussian | static | unbalanced, T=3 | sandwich | ✅ vs bootstrap | test-aipw-longitudinal.R |
| binary | binomial | static | unbalanced | sandwich | ✅ vs bootstrap | test-aipw-longitudinal.R |
| binary | gaussian | static + EM (by=sex) | unbalanced | sandwich | ✅ vs bootstrap | test-aipw-longitudinal.R |
| mv (A1,A2) | gaussian | static | unbalanced | sandwich | ✅ vs bootstrap + jackknife | test-aipw-longitudinal.R |
| binary | gaussian | static + external weights | unbalanced | sandwich | ✅ vs bootstrap | test-aipw-longitudinal.R |
| categorical (multinom) | gaussian | static | balanced | sandwich | ✅ vs bootstrap (softmax score) | test-aipw-longitudinal.R |
| binary | gaussian | static, `mgcv::gam` nuisance | any | sandwich | ✅ explicit limitation → bootstrap (`causatr_longitudinal_aipw_sandwich_model`) | test-aipw-longitudinal.R |

**Chunk 25 — Multivariate longitudinal AIPW**

Joint time-varying treatments `treatment = c("A1", "A2")` with `id =`/`time =`
under the doubly-robust ICE-AIPW estimator. The outcome side reuses the ICE
backward recursion (already loops over vector treatment); the propensity side
reuses the Phase 19 per-period × per-component density chain
(`fit_treatment_models()`, `compute_density_ratio_weights_mv()`,
`make_weight_fn_mv()`). The sandwich stacks T×K block-diagonal propensity
corrections in `variance_if_aipw_long_one()` Channel 2b.

| Feature | Status | Test |
|---|---|---|
| MV longitudinal AIPW: static recovers analytical truth by sex (9 / 15) | ✅ | test-aipw.R (`T-long-mv-aipw1`) |
| MV longitudinal AIPW: sandwich vs bootstrap SE parity (within 30%) | ✅ | test-aipw.R (`T-long-mv-aipw2`) |
| MV longitudinal AIPW: cross-method agreement with MV ICE g-comp + MV long-IPW | ✅ | test-aipw.R (`T-long-mv-aipw3`) |
| MV longitudinal AIPW: double-robustness (one-sided misspecification) | ✅ | test-aipw.R (`T-long-mv-aipw4`) |
| MV longitudinal AIPW: continuous shift finite + SE parity | ✅ | test-aipw.R (`T-long-mv-aipw5`) |
| `stabilize = "marginal"` rejected for longitudinal AIPW (univ + MV) | ✅ | test-aipw.R (`R-long-mv-aipw1`) |

**Rejections (MV longitudinal AIPW):** `ipsi()` components → rejected ✅ (`causatr_longitudinal_ipsi_pending`, inherited from the univariate path); transport (`target =`) → owned by Phase 26 (design doc).

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

**Chunk 17k — IPCW × transport**

| Feature | Status | Test |
|---|---|---|
| gcomp + IPCW + transport (transportability): recovers target ATE | ✅ | test-ipcw-transport.R |
| gcomp + IPCW + transport (generalizability): recovers marginal ATE | ✅ | test-ipcw-transport.R |
| gcomp + IPCW + transport: sandwich SE finite, CI covers truth | ✅ | test-ipcw-transport.R |
| IPW + IPCW + transport (transportability): recovers target ATE | ✅ | test-ipcw-transport.R |
| IPW + IPCW + transport (generalizability): recovers marginal ATE | ✅ | test-ipcw-transport.R |
| IPW + IPCW + transport: sandwich vs bootstrap SE ratio in (0.4, 2.5) | ✅ | test-ipcw-transport.R |
| AIPW + IPCW + transport (transportability): recovers target ATE | ✅ | test-ipcw-transport.R |
| AIPW + IPCW + transport: sandwich SE finite, CI covers truth | ✅ | test-ipcw-transport.R |
| AIPW + IPCW + transport: 2-of-3 DR (wrong outcome model) | ✅ | test-ipcw-transport.R |
| AIPW + IPCW + transport: sandwich vs bootstrap SE ratio in (0.5, 2.0) | ✅ | test-ipcw-transport.R |
| Cross-estimator agreement: gcomp ≈ IPW ≈ AIPW under correct specification | ✅ | test-ipcw-transport.R |
| IPCW + transport corrects bias vs naive complete-case transport | ✅ | test-ipcw-transport.R |
| Bootstrap refits sampling + censoring + propensity per replicate | ✅ | test-ipcw-transport.R |
| Details stashing: transport + ipcw flags and models present | ✅ | test-ipcw-transport.R |
| TransportHealth cross-check: gcomp + manual IPCW weights | ✅ | test-ipcw-transport.R |
| delicatessen stacked-EE cross-check: IPW point + SE | ✅ | test-ipcw-transport.R |
| Longitudinal IPW + IPCW + transport: smoke test (finite output) | ✅ | test-ipcw-transport.R |

**Implementation note.** No source code changes were required. The runtime path in `causat.R` already composes transport (sampling model) and IPCW (censoring model) independently and sequentially; the sandwich variance engines (`variance_if_ipw`, `variance_if_gcomp`, `variance_if_aipw`) already apply both the IPCW and sampling-model corrections when the respective flags are active. This chunk is purely test coverage validating the four-block composition.

**Chunk 17l — MTP + transportability (MC marginalization)**

For MTP interventions (shift, scale_by, threshold, dynamic) combined with transportability, target rows (S=0) have no observed treatment. MC marginalization integrates $E_{A|L,S=1}[\hat{m}(d(A,L), L)]$ over the study treatment distribution: exact enumeration for binary treatment, Monte Carlo draws for continuous. Sandwich variance is rejected (classed error `causatr_mtp_transport_sandwich`) for gcomp and AIPW because the MC marginalization introduces treatment-model dependence not captured by the current IF chain; bootstrap is supported. IPW + MTP + transport does not require MC marginalization (density-ratio weights operate on study rows only) so sandwich works.

| Feature | Status | Test |
|---|---|---|
| gcomp + shift + transport: recovers analytical truth | ✅ | test-mtp-transport.R |
| gcomp + scale_by + transport: recovers truth | ✅ | test-mtp-transport.R |
| gcomp + threshold + transport: runs without error | 🟡 | test-mtp-transport.R |
| gcomp + dynamic + transport: recovers truth | ✅ | test-mtp-transport.R |
| gcomp + shift + generalizability: recovers truth | ✅ | test-mtp-transport.R |
| IPW + shift + transport: recovers truth | ✅ | test-mtp-transport.R |
| IPW + ipsi + transport: runs and produces finite estimates (binary DGP) | 🟡 | test-mtp-transport.R |
| AIPW + shift + transport: recovers truth | ✅ | test-mtp-transport.R |
| AIPW + shift + transport: DR under wrong outcome model | ✅ | test-mtp-transport.R |
| AIPW + shift + transport: smoke test | 🟡 | test-mtp-transport.R |
| Cross-estimator agreement: gcomp shift ≈ IPW shift | ✅ | test-mtp-transport.R |
| gcomp + shift + transport: bootstrap SE finite and reasonable | ✅ | test-mtp-transport.R |
| gcomp + shift + transport: sandwich is rejected | ✅ | test-mtp-transport.R |
| AIPW + shift + transport: sandwich is rejected | ✅ | test-mtp-transport.R |
| IPW + shift + transport: sandwich works (no MC marginalization) | ✅ | test-mtp-transport.R |
| gcomp + static + transport: sandwich still works (no MC needed) | ✅ | test-mtp-transport.R |
| gcomp + dynamic + binary treatment + transport: exact marginalization | ✅ | test-mtp-transport.R |

**Phase 17 complete.** Chunk 17m (documentation, vignette, coverage matrix) shipped 2026-05-20.

---

## Cross-validation & stress tests

### NHEFS published-result replication (Hernán & Robins 2025)

| Estimator | Estimand | Published | Status | Test |
|---|---|---|---|---|
| G-comp (point) | ATE (wt82_71) | ≈ 3.5 kg | ✅ | test-nhefs-replication.R |
| IPW (Hájek) | ATE | ≈ 3.44 kg | ✅ | test-nhefs-replication.R |
| AIPW | ATE | ≈ 3.46 kg | ✅ | test-nhefs-replication.R |
| Matching | ATT | ballpark 3–4 kg | ✅ | test-nhefs-replication.R |
| Cross-method agreement | ATE | pairwise < 0.5 kg | ✅ | test-nhefs-replication.R |

### LaLonde / Dehejia-Wahba replication

| Estimator | Estimand | Published | Status | Test |
|---|---|---|---|---|
| Matching | ATT (re78) | ≈ $1,794 | ✅ | test-external-replication.R |
| G-comp + IPW + matching | ATT agreement | pairwise < $2000 | ✅ | test-external-replication.R |
| AIPW | ATT | agrees with gcomp < $1500 | ✅ | test-external-replication.R |

### Kang-Schafer adversarial DR test (2007)

| Scenario | Models | Expected | Status | Test |
|---|---|---|---|---|
| S1 | Both correct (Z) | All recover 210 | ✅ | test-kang-schafer.R |
| S2 | Outcome wrong (X), PS correct (Z) | IPW recovers; gcomp biased | ✅ | test-kang-schafer.R |
| S3 | PS wrong (X), outcome correct (Z) | gcomp recovers; IPW biased | ✅ | test-kang-schafer.R |
| S4 | Both wrong (X) | All biased (neg. control) | ✅ | test-kang-schafer.R |
| DR property | AIPW(Z) << AIPW(X) | Correct < misspecified bias | ✅ | test-kang-schafer.R |

### Cross-estimator triangulation

| Scenario | Estimators | Status | Test |
|---|---|---|---|
| Het. effect DGP — ATE | gcomp, IPW, AIPW, matching | ✅ | test-triangulation.R |
| Het. effect DGP — ATT | gcomp, IPW, matching | ✅ | test-triangulation.R |
| Binary outcome — RD | gcomp, IPW, AIPW | ✅ | test-triangulation.R |
| Continuous trt — shift | gcomp, IPW, AIPW | ✅ | test-triangulation.R |
| Categorical trt — static | gcomp, IPW | ✅ | test-triangulation.R |

### Naimi longitudinal replication (2017-inspired)

| Estimator | Truth (MC) | Status | Test |
|---|---|---|---|
| ICE g-comp | ATE ≈ MC truth | ✅ | test-naimi-replication.R |
| Longitudinal IPW | ATE ≈ MC truth | ✅ | test-naimi-replication.R |
| Longitudinal AIPW | ATE ≈ MC truth | ✅ | test-naimi-replication.R |
| Cross-method | 3-way < 100 | ✅ | test-naimi-replication.R |
| Naive regression | biased > 30 | ✅ | test-naimi-replication.R |

### Combinatorial stress test

Systematic grid (52 tests): estimator × treatment type × outcome family × intervention × variance method. All valid cells produce finite estimates, positive SEs, and match truth within tolerance.

Test file: test-combinatorial-stress.R

---

### Phase 18 — G-estimation of Structural Nested Mean Models
Third leg of the Robins triangle. Motivating use case: **correct handling of time-varying effect modification** — SNMs parameterise the per-stage blip $\gamma_k(a_k, \bar{l}_k, \bar{a}_{k-1}; \psi)$ directly and identify it via a moment condition that uses the treatment model as instrument, so time-varying modifiers are supported by design (closes the Phase 6 limitation under MSM-based estimators). Scope: linear-blip additive SNMMs for point + longitudinal, stacked EE sandwich ($K$ treatment blocks + blip block), bootstrap, `DTRreg` + `delicatessen` cross-checks. Survival SNMs (SNFTMs/SNCFTMs) out of scope.

**Chunk 18a — Routing and validation (shipped)**

| Feature | Status | Test |
|---|---|---|
| `causat(estimator = "snm")` routing to `fit_snm()` | ✅ | test-snm.R |
| Treatment model fit via `fit_treatment_model()` | ✅ | test-snm.R |
| Blip spec from EM terms (`A:modifier`) | ✅ | test-snm.R |
| Per-component confounders (`confounders_outcome`, `confounders_treatment`) | ✅ | test-snm.R |
| No-EM inform (blip = constant ATE) | ✅ | test-snm.R |
| Continuous treatment support | ✅ | test-snm.R |
| Multivariate treatment → rejection (`causatr_snm_multivariate`) | ✅ | test-snm.R |
| Longitudinal data accepted (`type = "longitudinal"`) | ✅ | test-snm.R |
| `contrast(interventions=)` → rejection (`causatr_snm_no_interventions`) | ✅ | test-snm.R |
| Missing treatment confounders → error | ✅ | test-snm.R |
| `propensity_model_fn` default warning | ✅ | test-snm.R |
| No `model_fn` default warning for SNM | ✅ | test-snm.R |
| `build_blip_spec()` — no EM / single / multiple modifiers | ✅ | test-snm.R |

**Chunk 18b — Point estimation, contrast, and sandwich variance (shipped)**

| Feature | Status | Test |
|---|---|---|
| `compute_snm_blip_point()` — closed-form g-estimation | ✅ | test-snm.R |
| `contrast(fit)` — blip parameter table (Path A) | ✅ | test-snm.R |
| `contrast(fit, treatment_values = c(0, 1))` — averaged blip (Path B) | ✅ | test-snm.R |
| `variance_if_snm()` — stacked EE sandwich | ✅ | test-snm.R |
| Continuous trt + EM (design doc DGP) — truth | ✅ | test-snm.R |
| Continuous trt + no EM — truth | ✅ | test-snm.R |
| Binary trt + EM — truth | ✅ | test-snm.R |
| Multiple modifiers (2 EM terms) — truth | ✅ | test-snm.R |
| Large-sample convergence (n=20000) | ✅ | test-snm.R |
| Sandwich SE vs manual bootstrap consistency | ✅ | test-snm.R |
| delicatessen cross-check — continuous trt, EM | ✅ | test-snm.R |
| delicatessen cross-check — binary trt, EM | ✅ | test-snm.R |
| delicatessen cross-check — continuous trt, 2 modifiers | ✅ | test-snm.R |
| DTRreg cross-check — binary trt, no EM (point + SE) | ✅ | test-snm.R |
| `contrast(interventions=)` → rejection | ✅ | test-snm.R |
| `ci_method = "bootstrap"` accepted (point + longitudinal) | ✅ | test-snm-bootstrap.R (see Chunk 18i) |
| `treatment_values` on non-SNM → rejection | ✅ | test-snm.R |
| `treatment_values` length ≠ 2 → rejection | ✅ | test-snm.R |
| Vcov dimensions and PSD | ✅ | test-snm.R |

**Chunk 18b½ — Treatment-free outcome model for efficiency augmentation (shipped)**

Joint estimation of (β, ψ) following Vansteelandt & Joffe (2014) and DTRreg's `tf.mod`. The treatment-free model E[Y | L] absorbs the L→Y variance, reducing SEs by 30–45% with unchanged point estimates.

| Feature | Status | Test file |
|---|---|---|
| `treatment_free = ~ L` argument in `causat()` | ✅ | test-snm.R |
| Joint (β, ψ) estimation in `compute_snm_blip_point()` | ✅ | test-snm.R |
| 3-block sandwich: (α_trt, θ_joint = (β, ψ)) | ✅ | test-snm.R |
| Continuous trt + EM + TF — truth | ✅ | test-snm.R |
| No EM + TF — truth | ✅ | test-snm.R |
| Binary trt + EM + TF — truth | ✅ | test-snm.R |
| `treatment_values` + TF — averaged blip | ✅ | test-snm.R |
| SE reduction across all DGPs (TF < no-TF) | ✅ | test-snm.R |
| TF sandwich SE vs manual bootstrap consistency | ✅ | test-snm.R |
| delicatessen cross-check — continuous trt, EM, TF | ✅ | test-snm.R |
| delicatessen cross-check — binary trt, EM, TF | ✅ | test-snm.R |
| delicatessen cross-check — continuous trt, 2 mods, TF | ✅ | test-snm.R |
| DTRreg cross-check — binary trt, EM, TF (point + SE) | ✅ | test-snm.R |
| DTRreg cross-check — binary trt, no EM, TF (point + SE) | ✅ | test-snm.R |
| `treatment_free` on non-SNM → rejection (`causatr_treatment_free_not_snm`) | ✅ | test-snm.R |
| `treatment_free` non-formula → rejection (`causatr_snm_bad_treatment_free`) | ✅ | test-snm.R |
| Vcov PSD with TF model | ✅ | test-snm.R |

**Chunk 18c — Time-varying (post-treatment) effect modification (shipped)**

SNMs identify the blip under treatment-model correctness alone, so modifiers that depend on treatment (post-treatment) are valid — unlike IPW-MSM, which conditions on a descendant of A (collider bias; Robins 2000). This chunk adds truth-based tests, delicatessen cross-checks, DTRreg cross-checks, and an IPW bias demonstration on DGPs with genuinely post-treatment modifiers.

| Feature | Status | Test file |
|---|---|---|
| Continuous trt + TV modifier + TF — truth (ψ₀ = 3, ψ_M = 2) | ✅ | test-snm.R |
| Binary trt + TV modifier + TF — truth | ✅ | test-snm.R |
| TV modifier without TF — runs, different blip quantity | ✅ | test-snm.R |
| TF model reduces SEs on TV modifier DGP | ✅ | test-snm.R |
| delicatessen cross-check — continuous trt, TV modifier, no TF | ✅ | test-snm.R |
| delicatessen cross-check — continuous trt, TV modifier, TF | ✅ | test-snm.R |
| delicatessen cross-check — binary trt, TV modifier, no TF | ✅ | test-snm.R |
| delicatessen cross-check — binary trt, TV modifier, TF | ✅ | test-snm.R |
| DTRreg cross-check — binary trt, TV modifier, TF | ✅ | test-snm.R |
| IPW biased with post-treatment modifier (SNM is not) | ✅ | test-snm.R |
| Vcov PSD for TV modifier (no TF) | ✅ | test-snm.R |
| Vcov PSD for TV modifier (with TF) | ✅ | test-snm.R |

**Chunk 18d — Longitudinal SNMM with per-stage blip estimation (shipped)**

Backward sequential g-estimation (Robins 1994): per-stage blip parameters estimated by backward induction from stage K to stage 0. Cluster-robust sandwich with cross-stage derivatives. DTRreg cross-check validates point estimates.

| Feature | Status | Test file |
|---|---|---|
| `compute_snm_blip_longitudinal()` — backward sequential estimation | ✅ | test-snm.R |
| Per-stage blip parameters (`stage0_psi_intercept`, ...) | ✅ | test-snm.R |
| Binary treatment, 2-period DGP — truth (ψ₀ = 3.15, ψ₁ = 3) | ✅ | test-snm.R |
| Continuous treatment, 2-period DGP — truth | ✅ | test-snm.R |
| `variance_if_snm_longitudinal()` — cluster-robust sandwich | ✅ | test-snm.R |
| Sandwich SE vs cluster bootstrap consistency | ✅ | test-snm.R |
| Vcov dimensions (K×p_psi × K×p_psi) and PSD | ✅ | test-snm.R |
| DTRreg cross-check — binary trt, 2-period | ✅ | test-snm.R |
| Treatment-free model accepted (point estimates consistent) | ✅ | test-snm.R |
| TF joint (β,ψ) variance — psi marginal extraction | ✅ | test-snm.R |
| TF variance efficiency gain when E[R·Z] ≠ 0 | ✅ | test-snm.R |
| `treatment_values` rejected for longitudinal | ✅ | test-snm.R |
| Bootstrap accepted (smoke test) | ✅ | test-snm.R |
| `< 2 time points` → rejection | ✅ | test-snm.R |
| `n_obs` and metadata in result | ✅ | test-snm.R |

**delicatessen reference CSV:** `tests/testthat/fixtures/snm_delicatessen_reference.csv` (10 scenarios, generated by `data-raw/snm_reference.py`).

**Chunk 18e — Longitudinal time-varying EM truth-based test (shipped)**

Headline scientific demonstration: time-varying modifier $M_1 = 1\{L_1 > 0\}$ is post-treatment ($L_1$ depends on $A_0$). SNMs identify blip parameters correctly; IPW-MSM conditioning on $M_1$ is biased (collider bias). DGP: `simulate_snm_longitudinal_tv_em()` with truth $\psi_{00} = 1.15$, $\psi_{0M} = 2$, $\psi_{10} = 2$, $\psi_{1M} = 2$.

| Feature | Status | Test |
|---|---|---|
| Per-stage blip with TV modifier — truth (4 params) | ✅ | test-snm.R |
| TF model + TV-EM: same truth, tighter SEs | ✅ | test-snm.R |
| DTRreg cross-check — TV-EM, `history=0` for exact PS match | ✅ | test-snm.R |
| IPW biased with post-treatment modifier; SNM is not | ✅ | test-snm.R |
| Vcov PSD and correct size (4×4) for TV-EM | ✅ | test-snm.R |
| Sandwich SE vs cluster bootstrap for TV-EM | ✅ | test-snm.R |

**Chunk 18f — Triangulation tests + delicatessen cross-check + `history = 0` (shipped)**

Validates the Robins triangle invariant: under correct specification and no TV-EM, SNM blip sum = IPW-MSM ATE = ICE g-comp ATE for binary static interventions. Also adds `delicatessen` cross-language reference values and `history = 0` (no-lag) support for all longitudinal estimators.

| Feature | Status | Test |
|---|---|---|
| SNM blip sum ≈ IPW ATE (binary, no EM, n=5000) | ✅ | test-snm-triangulation.R |
| SNM blip sum ≈ ICE ATE (binary, no EM, n=5000) | ✅ | test-snm-triangulation.R |
| 3-way triangulation: SNM ≈ IPW ≈ ICE (n=8000, tighter tol) | ✅ | test-snm-triangulation.R |
| SNM correct but IPW biased with TV-EM (negative control) | ✅ | test-snm-triangulation.R |
| delicatessen cross-check — base (no EM, no TF) point + SE | ✅ | test-snm-triangulation.R |
| delicatessen cross-check — TF (no EM, treatment-free ~L) point + SE | ✅ | test-snm-triangulation.R |
| delicatessen cross-check — TV-EM + TF (A:M modifier, treatment-free ~L) point + SE | ✅ | test-snm-triangulation.R |
| `history = 0` — longitudinal IPW (no lags, ATE recovery) | ✅ | test-snm-triangulation.R |
| `history = 0` — longitudinal ICE (no lags, finite output) | ✅ | test-snm-triangulation.R |
| `history = 0` — longitudinal SNM (no lags, blip recovery) | ✅ | test-snm.R |
| `history` validator: non-negative integer or Inf | ✅ | test-checks.R, test-causat.R |

**delicatessen fixtures:** `tests/testthat/fixtures/snm_longitudinal_delicatessen.csv` (base), `snm_longitudinal_tf_delicatessen.csv` (TF), `snm_longitudinal_tv_em_delicatessen.csv` (TV-EM + TF). Generated by `data-raw/snm_longitudinal_reference.py` from shared fixtures in `data-raw/`.

**Chunk 18i — Bootstrap variance for SNM (shipped)**

Bootstrap variance for point and longitudinal SNM estimators. Point bootstrap resamples rows and re-solves the g-estimating equation (with or without treatment-free model). Longitudinal bootstrap uses ID-clustered resampling + full `fit_snm()` refit. Forwards per-component confounders (`confounders_outcome`, `confounders_tv_*`) to ensure correct blip spec on bootstrap replicates. `treatment_values` path bootstraps the averaged blip effect directly (scalar statistic, no delta method needed).

| Test | Status | File |
|---|---|---|
| Point bootstrap SE consistent with sandwich (continuous, EM) | ✅ | test-snm-bootstrap.R |
| Point bootstrap no EM (single psi) | ✅ | test-snm-bootstrap.R |
| Point bootstrap + treatment-free: SE consistent | ✅ | test-snm-bootstrap.R |
| Point bootstrap with `treatment_values` | ✅ | test-snm-bootstrap.R |
| Point bootstrap: binary treatment + EM | ✅ | test-snm-bootstrap.R |
| Longitudinal bootstrap SE consistent with sandwich | ✅ | test-snm-bootstrap.R |
| Longitudinal bootstrap + treatment-free | ✅ | test-snm-bootstrap.R |
| Longitudinal bootstrap with TV-EM (`confounders_outcome`) | ✅ | test-snm-bootstrap.R |
| Bootstrap accepted (previously rejected) smoke test | ✅ | test-snm-bootstrap.R |

**Chunk 18h — Categorical and count treatment extensions (shipped)**

Extends SNM g-estimation to categorical (k>2, multinomial residualisation via `nnet::multinom`) and count (Poisson/NB, scalar residual R = A − λ̂) treatment types. `snm_treatment_design()` abstracts treatment design matrices (AM, RM) uniformly across all treatment types. For categorical, psi expands to (K−1) blocks with per-level indicators and multinomial residuals. `treatment_values` is rejected for categorical SNM (per-level blip parameters are the estimand).

| Test | Status | File |
|---|---|---|
| Poisson count trt, no EM — truth (ψ = 0.5) | ✅ | test-snm.R |
| Poisson count trt + EM — truth (ψ₀ = 0.5, ψ_M = 0.3) | ✅ | test-snm.R |
| Poisson sandwich SE vs bootstrap consistency | ✅ | test-snm.R |
| Poisson + TF model — truth | ✅ | test-snm.R |
| Poisson + `treatment_values` — averaged blip | ✅ | test-snm.R |
| Longitudinal Poisson — per-stage truth | ✅ | test-snm.R |
| 3-level categorical, no EM — truth (ψ₁ = 0.8, ψ₂ = 1.5) | ✅ | test-snm.R |
| Categorical + EM — 4 blip params (per-level × per-modifier) | ✅ | test-snm.R |
| Categorical sandwich SE vs bootstrap consistency | ✅ | test-snm.R |
| Categorical + TF model — truth | ✅ | test-snm.R |
| Categorical vcov PSD and correct dimensions | ✅ | test-snm.R |
| `treatment_values` rejected for categorical (`causatr_snm_cat_no_treatment_values`) | ✅ | test-snm.R |
| Longitudinal categorical — per-stage per-level truth | ✅ | test-snm.R |

**Chunk 18j — S3 dispatch for SNM results (shipped)**

| Feature | Status | Test |
|---|---|---|
| `print.causatr_result()` — "Blip parameters" header for SNM Path A | ✅ | test-s3-methods.R |
| `print.causatr_result()` — "Averaged blip effect" for SNM Path B | ✅ | test-s3-methods.R |
| `print.causatr_result()` — "Per-stage blip parameters" for longitudinal SNM | ✅ | test-s3-methods.R |
| `summary.causatr_result()` — skips intervention details for SNM + blip vcov label | ✅ | test-s3-methods.R |
| `coef.causatr_result()` — returns named blip parameter vector | ✅ | test-s3-methods.R |
| `tidy.causatr_result()` — handles `parameter` column, type = "parameter" | ✅ | test-s3-methods.R |
| `tidy.causatr_result()` — SNM Path B (averaged blip) | ✅ | test-s3-methods.R |
| `glance.causatr_result()` — works for SNM | ✅ | test-s3-methods.R |
| `plot.causatr_result()` — forest plot with "Parameter" header + ref_line = 0 | ✅ | test-s3-methods.R |
| `plot.causatr_result()` — longitudinal SNM title variation | ✅ | test-s3-methods.R |

**Chunk 18b-ext — Combination matrix follow-ups (shipped 2026-05-27)**

Groups 1–6 from the Phase 18 combination matrix audit. Addresses all 🟡 items deferred from the initial shipping.

| Feature | Status | Test |
|---|---|---|
| Binomial outcome truth test (linear probability DGP, n=3000) | ✅ | test-snm.R |
| Quasibinomial outcome truth test | ✅ | test-snm.R |
| Poisson outcome truth test | ✅ | test-snm.R |
| Gamma outcome truth test | ✅ | test-snm.R |
| Negative binomial outcome truth test | ✅ | test-snm.R |
| Betareg outcome truth test | ✅ | test-snm.R |
| GAM propensity smoke test (`mgcv::gam`, finite estimates+SEs) | ✅ | test-snm.R |
| `predict.gam()` 1-D array bug fix (`as.numeric()` wrap) | ✅ | test-snm.R |
| Cluster-robust SE: `cluster =` threading via `vcov_from_if()` | ✅ | test-snm.R |
| Cluster-robust SE: singleton cluster = i.i.d. exactly | ✅ | test-snm.R |
| Cluster-robust SE: fit-time propagation (`causat(cluster=)` → `contrast()`) | ✅ | test-snm.R |
| Cluster-robust SE: dimension-safe subsetting to fit rows | ✅ | test-snm.R |
| `by =`-stratified averaged blip: Path C in `compute_snm_contrast()` | ✅ | test-snm.R |
| `by =` per-stratum truth (M=0 → 3, M=1 → 5 for simulate_snm_point DGP) | ✅ | test-snm.R |
| `by =` weighted average of strata = pooled exactly | ✅ | test-snm.R |
| `by =` rejection without `treatment_values` (`causatr_snm_by_needs_treatment_values`) | ✅ | test-snm.R |
| `by =` rejection: column not in data (`causatr_snm_by_not_found`) | ✅ | test-snm.R |
| `weights =` propagation to blip EE (RM_w = RM * w) | ✅ | test-snm.R |
| `weights =` propagation to sandwich bread and score | ✅ | test-snm.R |
| `weights =` manual WLS formula match (`sum(w*R*Y)/sum(w*R*A)`) | ✅ | test-snm.R |
| `weights = rep(1, n)` recovers unweighted exactly | ✅ | test-snm.R |
| Sandwich fallback: `causatr_snm_no_estfun` for unsupported model class | ✅ | test-snm.R |
| DTRreg cross-check — binary × no EM: point estimate to 1e-6 | ✅ | test-snm-dtrreg-oracle.R |
| DTRreg cross-check — binary × EM + TF: point + SE to 1e-6 / 1e-2 | ✅ | test-snm-dtrreg-oracle.R |
| Continuous truth (no DTRreg — linear blip unsupported by DTRreg) | ✅ | test-snm-dtrreg-oracle.R |
| MCAR outcomes: causatr CC agrees with DTRreg CC to 1e-6 | ✅ | test-snm-dtrreg-oracle.R |
| Observation weights: WLS oracle `sum(w*R*Y)/sum(w*R*A)` match | ✅ | test-snm-dtrreg-oracle.R |
| Cluster singleton: `cluster_vec = 1:n` = i.i.d. exactly | ✅ | test-snm-dtrreg-oracle.R |
| Cluster propagation: fit-time cluster = explicit `cluster =` argument | ✅ | test-snm-dtrreg-oracle.R |
| `by =` per-stratum: manual delta-method formula match to 1e-6 | ✅ | test-snm-dtrreg-oracle.R |
| Longitudinal × binary × no EM: causatr vs DTRreg within 0.05 | ✅ | test-snm-dtrreg-oracle.R |
| Longitudinal × TV-EM: blip params within 0.20 of truth | ✅ | test-snm-dtrreg-oracle.R |
| Sandwich fallback: `variance_if_snm()` throws `causatr_snm_no_estfun` | ✅ | test-snm-dtrreg-oracle.R |

**Phase 18 COMPLETE.** All chunks shipped (18a–18f, 18h–18k, 18b-ext). 18g dropped (gesttools archived).

### Phase 19-trim — Cross-cutting weight truncation (shipped)

`trim =` argument on `contrast()`. Winsorizes per-component density ratios at the `trim`-th quantile before products (Cole & Hernan 2008; lmtp default 0.999).

| Feature | Status | Test |
|---|---|---|
| `truncate_weights()` helper: trim=1 no-op | ✅ | test-ipw-weights.R |
| `truncate_weights()` helper: clips at quantile | ✅ | test-ipw-weights.R |
| Point IPW + trim: max weight reduced | ✅ | test-ipw-weights.R |
| Continuous IPW shift + trim | ✅ | test-ipw-weights.R |
| IPSI + trim | ✅ | test-ipw-weights.R |
| `make_weight_fn` closure agrees at alpha_hat with trim | ✅ | test-ipw-weights.R |
| Multivariate IPW + trim | ✅ | test-ipw-weights.R |
| Sandwich SE agrees with bootstrap under trim | ✅ | test-ipw-weights.R |
| Monotonic max-weight reduction (0.999 > 0.99 > 0.95) | ✅ | test-ipw-weights.R |
| Near-positivity DGP: trim stabilizes extreme weights | ✅ | test-ipw-weights.R |
| AIPW + trim: doubly-robust property preserved | ✅ | test-ipw-weights.R |
| Longitudinal IPW + trim | ✅ | test-longitudinal-ipw.R |
| Longitudinal IPW + trim + bootstrap | ✅ | test-longitudinal-ipw.R |
| MV IPW + trim: lmtp_sdr cross-check | ✅ | test-ipw-weights.R |
| Long IPW + trim: lmtp_sdr cross-check | ✅ | test-longitudinal-ipw.R |

### `causat_mice()` — Multiple imputation
Pool `causat()` + `contrast()` across `mice` imputations. `pool_method = "rubin"`
(default, Rubin's rules + Barnard-Rubin df) and `pool_method = "boot_mi"` (von
Hippel two-stage bootstrap-then-impute).

| Feature | Status | Test file |
|---|---|---|
| Rubin pooling math matches `mice::pool.scalar()` (1e-8) | ✅ | test-pool-rubin.R |
| Barnard-Rubin df matches `mice:::barnard.rubin()` | ✅ | test-pool-rubin.R |
| Boot MI von Hippel decomposition (unit, exact) | ✅ | test-pool-boot-mi.R |
| gcomp + MI: MAR covariate, recovers ATE | ✅ | test-causat-mice.R |
| ipw / aipw / matching + MI: recover ATE | ✅ | test-causat-mice.R |
| snm + MI: per-blip-parameter pooling | ✅ | test-causat-mice.R |
| MAR treatment: impute A from P(A\|L,Y) | ✅ | test-causat-mice.R |
| Continuous / categorical / count treatment + MI | ✅ | test-causat-mice.R |
| Multivariate treatment (gcomp + IPW, incl. stabilized) + MI | ✅ | test-causat-mice.R |
| Outcome families: binomial (ratio/OR), poisson, Gamma | ✅ | test-causat-mice.R |
| Interventions: shift, scale_by, dynamic, threshold, ipsi | ✅ | test-causat-mice.R |
| ATT + MI; `by`-stratified pooling (per-row) | ✅ | test-causat-mice.R |
| Longitudinal (ICE + longitudinal IPW) + MI | ✅ | test-causat-mice.R |
| IPCW + MI: missing L imputed, censored Y reweighted | ✅ | test-causat-mice.R |
| Boot MI recovers ATE; tighter SE than Rubin under misspecified imputation | ✅ | test-pool-boot-mi.R |
| MC coverage: Rubin ~nominal when estimator consistent (DGP-MI1) | ✅ | test-mi-coverage.R |
| MC coverage: Boot MI ~nominal, SE no wider than Rubin (consistent regime) | ✅ | test-mi-coverage.R |
| MC coverage: omitting Y biases estimate, collapses coverage for both (DGP-MI5) | ✅ | test-mi-coverage.R |
| S3: confint / tidy honor Barnard-Rubin df | ✅ | test-causat-mice.R |
| `parallel = "future"` Boot MI | ✅ | test-pool-boot-mi.R |
| Edge: m=1 degenerate (warn); all-fit-failure (abort) | ✅ | test-causat-mice.R |
| Warn when analysis var absent/unused in imputation model | ✅ | test-causat-mice.R |
| Reject non-`mids` input; reject unknown `pool_method` | ✅ | test-causat-mice.R |
| Reject Boot MI with `B < 2` or `M < 2` | ✅ | test-pool-boot-mi.R |
| Boot MI warns on floored (non-positive) variance component | ✅ | test-pool-boot-mi.R |

Transport (`target =`) + MI is not yet validated (owned by the pending Phase
17i / 26 transport work).

---

## Maintenance rules

1. Update this matrix when adding/removing/changing a feature.
2. Add truth-based tests when feasible; use lmtp/delicatessen as external reference.
3. Test rejection paths with `expect_snapshot(error = TRUE)`.
4. Divergence between this matrix and test files is a bug.
