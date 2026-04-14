# Feature × Test Coverage Matrix

This file enumerates every supported (and planned) combination of
features in causatr and tracks the test-coverage status of each.
**Every PR that adds, removes, or changes a feature MUST update this
table and the corresponding tests.** Coverage is graded as:

| Symbol | Meaning |
|---|---|
| ✅ **truth** | Truth-based simulation: estimate / SE / CI all checked against an analytical or external reference (lmtp, delicatessen, Hernán & Robins published numbers). |
| 🟡 **smoke** | Runs without error, returns finite estimates / SEs, but the numerical target is not pinned. |
| ❌ **none** | No test exists. |
| ⛔ **n/a** | Not supported by this method (and the rejection path is also tracked separately). |

External references used:
- **lmtp**: gold standard for longitudinal stochastic interventions
  (`lmtp::lmtp_tmle`, `lmtp::lmtp_sdr`).
- **Hernán & Robins** book values (NHEFS examples, Table 20.1).
- **Closed-form analytical truth** for linear-Gaussian / saturated
  binary DGPs that admit a hand derivation.

## How to read the matrix

Each row pins down ONE combination across every dimension. The matrix
is split by **method** (gcomp, IPW, matching, ICE, survival) because
the supported intervention / variance shapes differ across methods.

## Dimensions

| Dimension | Values |
|---|---|
| **Method** | gcomp, IPW, matching, ICE, survival |
| **Time structure** | point, longitudinal |
| **Treatment type** | binary, categorical (k > 2), continuous, multivariate |
| **Outcome family** | gaussian, binomial, poisson, quasibinomial, gamma |
| **Model class** | GLM (`stats::glm`), GAM (`mgcv::gam`), other `model_fn` |
| **Intervention** | NULL (natural), `static`, `shift`, `scale_by`, `threshold`, `dynamic`, `ipsi` |
| **Estimand** | ATE, ATT, ATC, by-stratified |
| **Contrast type** | difference, ratio, OR |
| **Variance** | sandwich (analytic), bootstrap, numeric Tier 1, numeric Tier 2 |
| **External weights** | none, survey/IPCW |
| **Censoring** | none, present (filtered out) |
| **Competing risks** | n/a (point), placeholder (survival, Phase 6) |

---

## Method 1 — Point-treatment g-computation (`method = "gcomp"`)

| Treatment | Outcome | Model | Intervention | Estimand | Contrast | Variance | Weights | Status | Test file |
|---|---|---|---|---|---|---|---|---|---|
| binary | gaussian | GLM | static | ATE | difference | sandwich | none | ✅ truth (analytical + NHEFS) | test-gcomp.R, test-simulation.R |
| binary | gaussian | GLM | static | ATE | difference | bootstrap | none | ✅ truth | test-gcomp.R, test-simulation.R |
| binary | gaussian | GLM | static | ATT | difference | sandwich | none | ✅ truth | test-gcomp.R, test-simulation.R |
| binary | gaussian | GLM | static | ATC | difference | sandwich | none | ✅ truth | test-simulation.R |
| binary | gaussian | GLM | static | ATE | difference | sandwich | survey | ✅ truth | test-simulation.R |
| binary | gaussian | GLM | static | ATE | difference | bootstrap | survey | 🟡 smoke | test-simulation.R |
| binary | binomial | GLM | static | ATE | difference (RD) | sandwich | none | ✅ truth | test-simulation.R |
| binary | binomial | GLM | static | ATE | ratio (RR) | sandwich | none | ✅ truth | test-simulation.R |
| binary | binomial | GLM | static | ATE | OR | sandwich | none | ✅ truth | test-simulation.R |
| binary | binomial | GLM | static | ATE | difference | bootstrap | none | ✅ truth | test-simulation.R |
| binary | poisson | GLM | static | ATE | ratio | sandwich | none | ✅ truth | test-gcomp.R |
| binary | quasibinomial | GLM | static | ATE | difference | sandwich | none | ❌ none | — |
| binary | gamma | GLM | static | ATE | ratio | sandwich | none | ❌ none | — |
| binary | gaussian | GAM | static | ATE | difference | sandwich | none | ✅ truth | test-complex-dgp.R |
| binary | gaussian | GAM | static | ATE | difference | bootstrap | none | ✅ truth | test-complex-dgp.R |
| binary | gaussian | GLM (splines) | static | ATE | difference | sandwich | none | ✅ truth | test-complex-dgp.R |
| binary | gaussian | GLM | static | by(L) | difference | sandwich | none | ✅ truth | test-by-estimand.R |
| binary | gaussian | GLM | static | by(L) | difference | bootstrap | none | ✅ truth | test-by-estimand.R |
| continuous | gaussian | GLM | shift | ATE | difference | sandwich | none | ✅ truth | test-gcomp.R, test-simulation.R |
| continuous | gaussian | GLM | shift | ATE | difference | bootstrap | none | ✅ truth | test-simulation.R |
| continuous | gaussian | GLM | scale_by | ATE | difference | sandwich | none | ✅ truth | test-simulation.R |
| continuous | gaussian | GLM | threshold | ATE | difference | sandwich | none | ✅ truth | test-simulation.R |
| binary | gaussian | GLM | dynamic | ATE | difference | sandwich | none | ✅ truth | test-simulation.R |
| multivariate | gaussian | GLM | static | ATE | difference | sandwich | none | ✅ truth | test-multivariate.R |
| multivariate | gaussian | GLM | shift | ATE | difference | sandwich | none | ✅ truth | test-multivariate.R |
| binary | gaussian | GLM | static | ATE | difference | sandwich | none + censoring | ✅ truth | test-simulation.R |
| binary | gaussian | unsupported `model_fn` | static | ATE | difference | numeric Tier 1 | none | 🟡 smoke | test-variance-if.R |
| binary | gaussian | unsupported `model_fn` | static | ATE | difference | numeric Tier 2 | none | ✅ truth (vs main path, to ~1%) | test-variance-if.R |
| categorical (k>2) | gaussian | GLM | static | ATE | difference | sandwich | none | ✅ truth | test-simulation.R |

**Error / rejection paths** (must also be tested):
- gcomp + invalid family string → ✅ test-gcomp.R
- gcomp + missing outcome col → ✅ test-causat.R
- gcomp + missing treatment col → ✅ test-causat.R

---

## Method 2 — Inverse probability weighting (`method = "ipw"`)

| Treatment | Outcome | Intervention | Estimand | Contrast | Variance | Weights | Status | Test file |
|---|---|---|---|---|---|---|---|---|
| binary | gaussian | static | ATE | difference | sandwich | none | ✅ truth | test-simulation.R |
| binary | gaussian | static | ATE | difference | bootstrap | none | ✅ truth | test-simulation.R |
| binary | gaussian | static | ATT | difference | sandwich | none | ✅ truth | test-simulation.R |
| binary | gaussian | static | ATC | difference | sandwich | none | ✅ truth | test-simulation.R |
| binary | gaussian | static | ATE | difference | sandwich | survey | 🟡 smoke | test-simulation.R |
| binary | binomial | static | ATE | difference (RD) | sandwich | none | ✅ truth | test-simulation.R |
| binary | binomial | static | ATE | ratio | sandwich | none | ✅ truth | test-simulation.R |
| binary | binomial | static | ATE | OR | sandwich | none | ✅ truth | test-simulation.R |
| binary | binomial | static | ATE | difference | bootstrap | none | ✅ truth | test-simulation.R |
| categorical (k>2) | gaussian | static | ATE | difference | sandwich | none | 🟡 smoke | test-ipw.R |
| continuous | gaussian | static | ATE | difference | sandwich | none | 🟡 smoke | test-ipw.R |
| binary | gaussian | shift / scale / threshold / dynamic / ipsi | — | — | — | — | ⛔ **rejected** (Phase 4) | test-ipw.R, test-contrast.R |
| binary | gaussian | static | ATE | difference | sandwich (non-Mparts WeightIt method) | none | ✅ e2e smoke (finite understated SE + bootstrap fallback) | test-ipw.R |
| multivariate | any | any | any | any | any | any | ⛔ **rejected** (Phase 4) | test-s3-methods.R, test-multivariate.R |
| longitudinal (any shape) | any | any | any | any | any | any | ⛔ **rejected** (Phase 4) | test-ipw.R |

---

## Method 3 — Matching (`method = "matching"`)

| Treatment | Outcome | Intervention | Estimand | Contrast | Variance | Weights | Status | Test file |
|---|---|---|---|---|---|---|---|---|
| binary | gaussian | static | ATT | difference | sandwich (cluster-robust) | none | ✅ truth | test-simulation.R |
| binary | gaussian | static | ATT | difference | bootstrap | none | ✅ truth | test-simulation.R |
| binary | gaussian | static | ATE (full matching auto) | difference | sandwich | none | ✅ truth | test-simulation.R |
| binary | gaussian | static | ATC | difference | sandwich | none | ✅ truth | test-simulation.R |
| binary | gaussian | static | ATT | difference | sandwich | survey | ✅ truth | test-matching.R |
| binary | gaussian | static | ATT | difference | bootstrap | survey | ✅ truth | test-matching.R |
| binary | binomial | static | ATT | difference | sandwich | none | ✅ truth | test-simulation.R |
| binary | binomial | static | ATT | ratio / OR | sandwich | none | ✅ truth | test-simulation.R |
| binary | quasibinomial | static | ATT | difference | sandwich | none | ✅ truth | test-simulation.R |
| categorical (k>2) | any | any | any | any | any | any | ⛔ **rejected** (MatchIt is binary-only; verified via docs) | test-matching.R |
| continuous | any | any | any | any | any | any | ⛔ **rejected** (no continuous matching in MatchIt) | test-matching.R |
| binary | gaussian | shift / scale / threshold / dynamic / ipsi | — | — | — | — | ⛔ **rejected** (Phase 4) | test-contrast.R |
| multivariate | any | any | any | any | any | any | ⛔ **rejected** (Phase 7) | test-s3-methods.R |
| longitudinal | any | any | any | any | any | any | ⛔ **rejected** | test-matching.R |

---

## Method 4 — ICE longitudinal g-computation (`method = "gcomp"`, `type = "longitudinal"`)

| Treatment | Outcome | Intervention | n_periods | Variance | Weights | Status | Test file |
|---|---|---|---|---|---|---|---|
| binary | gaussian | static (always vs never) | 2 | sandwich | none | ✅ truth (Table 20.1) | test-ice.R |
| binary | gaussian | static (always vs never) | 2 | bootstrap | none | ✅ truth | test-ice.R |
| binary | gaussian | dynamic | 2 | sandwich | none | ✅ truth | test-ice.R, test-simulation.R |
| binary | gaussian | static | 2 | sandwich | survey | ✅ truth | test-simulation.R |
| binary | gaussian | static | 2 | sandwich | survey | ✅ vs unweighted regression test | test-ice.R, test-simulation.R |
| binary | binomial | static | 3+ | sandwich | none | ✅ truth | test-ice.R, test-simulation.R |
| binary | binomial | dynamic | 3+ | sandwich | none | ✅ truth | test-ice.R |
| binary | binomial | static | 3+ | bootstrap | none | ✅ truth | test-ice.R |
| binary | binomial | static | 3+ | bootstrap (parallel) | none | ✅ truth | test-ice.R |
| continuous | gaussian | shift | 2 | sandwich | none | ✅ truth (vs lmtp) | test-simulation.R |
| continuous | gaussian | scale_by | 2 | sandwich | none | ✅ truth (vs lmtp) | test-simulation.R |
| continuous | gaussian | threshold | 2 | sandwich | none | 🟡 smoke | test-ice.R |
| continuous | gaussian | shift | 2 | bootstrap | none | ❌ none | **GAP** |
| binary | gaussian | static | 2 | sandwich | none + censoring | ✅ truth | test-ice.R |
| multivariate | gaussian | static | 2 | sandwich | none | ✅ truth | test-multivariate.R |
| multivariate | gaussian | shift | 2 | sandwich | none | ✅ truth | test-multivariate.R |
| binary | gaussian | ipsi | any | — | — | ⛔ **rejected** (Phase 4) | test-ice.R |
| binary | gaussian | static | 2 | bootstrap | survey | ❌ none | **GAP** |
| binary | gaussian | static | 2 | sandwich | none + ATT | ⛔ **rejected** (only ATE for ICE) | test-ice.R |

---

## Method 5 — Survival (`causat_survival()`)

| Treatment | Outcome | Hazard model | Intervention | Variance | Censoring | Competing | Status | Test file |
|---|---|---|---|---|---|---|---|---|
| binary | survival | pooled logistic + ns(time) | (no contrast yet) | — | none | none | ✅ fit-only smoke | test-s3-methods.R, test-simulation.R |
| binary | survival | pooled logistic + factor(time) | (no contrast yet) | — | none | none | ✅ fit-only smoke | test-s3-methods.R |
| binary | survival | pooled logistic | (no contrast yet) | — | present | none | ❌ none | **GAP** |
| any | survival | any | any | any | any | non-NULL | ⛔ **rejected** (Phase 6) | test-causat.R |
| binary | survival | pooled logistic | static | analytic | none | none | ❌ contrast for survival not implemented | **GAP / Phase 6** |

---

## Cross-cutting concerns

| Concern | Status | Notes |
|---|---|---|
| `to_person_period()` wide → long round-trip | ✅ truth | test-simulation.R |
| `to_person_period()` rejects duplicate ids | ✅ snapshot | test-simulation.R |
| `to_person_period()` rejects mismatched lengths | ✅ snapshot | test-simulation.R |
| `causat()` type auto-detection (NULL → point/longitudinal) | ✅ | test-causat.R |
| `causat()` rejects bad inputs (missing cols, bad estimand) | ✅ snapshots | test-causat.R |
| `contrast()` rejects bad estimand × method combos | ✅ snapshots | test-contrast.R |
| `contrast()` rejects bad reference / intervention shape | ✅ snapshots | test-contrast.R |
| `diagnose()` for gcomp / IPW / matching | ✅ smoke + snapshots | test-diagnose.R |
| `diagnose()` for longitudinal | ⛔ **rejected** (per-period diagnostics deferred to a future phase) | test-diagnose.R |
| `diagnose()` aborts on missing WeightIt treat.type | 🟡 unit | manual check; needs test |
| S3 methods: print, summary, plot, coef, vcov, confint, tidy, glance | ✅ smoke | test-s3-methods.R |
| `confint()` respects `level` arg (sandwich path) | ✅ smoke | test-s3-methods.R |
| `confint()` consistency between sandwich and bootstrap on `level` | ✅ unit + monotonicity | test-s3-methods.R |
| `causat_mice()` (multiple imputation wrapper) | 🟡 stubbed | test-causat-mice.R |
| Numeric variance Tier 1 (estfun-based fallback) | 🟡 unit | test-variance-if.R |
| Numeric variance Tier 2 (delta shortcut + warn) | ✅ e2e via custom `model_fn` | test-variance-if.R |
| Cluster-robust matching variance (vs an IF-level reference) | 🟡 e2e | test-simulation.R; **GAP** for IF-level test |

---

## Maintenance rules

> **EVERY PR that adds, removes, or changes a feature MUST:**
> 1. Update this matrix to reflect the new state (add a row / change
>    the status column).
> 2. Add or update the corresponding test(s). Truth-based when feasible.
> 3. If the feature is currently unsupported, add an explicit error
>    snapshot to lock in the rejection path.
> 4. If the feature is supported but the closed-form truth is hard,
>    use `lmtp::lmtp_tmle()` (R) or `delicatessen` (Python) as the
>    external reference and pin the expected value with the validation
>    in a comment block above the test.

This matrix is the single source of truth for "what works". Any
divergence between the matrix and the test files is a bug.
