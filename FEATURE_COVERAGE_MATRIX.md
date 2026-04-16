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

## 2026-04-15 critical-review regression coverage

The 2026-04-15 critical-review sweep (B1-B8, R3-R12) landed a set of
silent-correctness fixes. Each blocking fix has a paired regression
test in `tests/testthat/test-critical-review-2026-04.R`, targeted at
the specific failure mode rather than the full feature combination:

| Fix | Failure mode | Test |
|---|---|---|
| B1 | `subset = quote(age > cutoff)` could not resolve session vars inside `compute_contrast()` / boot workers. | `B1: subset expression resolves session-scoped variables`, `B1: subset works inside bootstrap workers too`, `B1: subset length mismatch aborts with a clear message` |
| B2 | Bootstrap refit dropped user `...`, so non-default IPW / matching / gcomp bootstrap SEs corresponded to a different estimator. | `B2: bootstrap refit replays user's ... (gcomp quasipoisson)`, `B2: IPW bootstrap replays stashed WeightIt dots` |
| B5 | `causat_survival()` censoring dropped only the current row, not subsequent rows. | `B5: causat_survival drops all rows at/after first censor` |
| B6 | External weights were post-multiplied onto `w$weights` after WeightIt ran, so Mparts IF under-corrected. | `B6: external weights enter WeightIt via s.weights` |
| B7 | ICE Ch1 IF weighted/unweighted formulas drifted unless `sum(w) == n_target`. | `B7: ICE Ch1 IF uniform-weighted matches unweighted (unification)` |
| B8 | `by` enumerated over full data and aborted on empty strata. | `B8: \`by\` skips empty strata instead of aborting` |
| R6 | OR validation aborted on NA `mu_hat` via `if (NA)`. | `R6: OR validation does not abort on NA mu_hat` |
| R12 | Reserved column collision (`.pseudo_y`, `.causatr_prev_event`) was only checked in one place. | `R12: reserved column name is rejected up front` |

Non-blocking fixes (R3, R7, R8, R9, R10, R11, S1, S5) are exercised
indirectly by the existing truth-based tests in
`test-variance-if.R` / `test-variance-reference.R` and the simulation
harness.

### Supplementary test files from the 2026-04-15 review rounds

| File | Scope |
|---|---|
| `test-weights-edge-cases.R` | Edge cases for the external-weights path across every method: `check_weights()` gate (NA, Inf, negative, non-numeric, wrong length); zero-weight pass-through; uniform weights recovering the unweighted result at machine precision (gcomp + ICE); heterogeneous weights producing non-trivial SE shifts; IPW survey-weight `s.weights` plumbing; matching weight composition; `causat_survival()` weight flow into the hazard fit. |
| `test-replay-fit.R` | Unit + end-to-end tests for the central `replay_fit()` helper (`R/utils.R`) introduced by the third-round dots audit: base-args precedence, positional-dot drop, `reserved`-key blocking, NULL dots, backward-compat with fit objects missing `$dots`, and real gcomp (quasipoisson) / IPW (cbps) bootstrap replay through `variance_bootstrap()`. |

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

## Method 1 — Point-treatment g-computation (`estimator = "gcomp"`)

| Treatment | Outcome | Model | Intervention | Estimand | Contrast | Variance | Weights | Status | Test file |
|---|---|---|---|---|---|---|---|---|---|
| binary | gaussian | GLM | static | ATE | difference | sandwich | none | ✅ truth (analytical + NHEFS) | test-gcomp.R, test-simulation.R |
| binary | gaussian | GLM | static | ATE | difference | bootstrap | none | ✅ truth | test-gcomp.R, test-simulation.R |
| binary | gaussian | GLM | static | ATT | difference | sandwich | none | ✅ truth | test-gcomp.R, test-simulation.R |
| binary | gaussian | GLM | static | ATC | difference | sandwich | none | ✅ truth | test-simulation.R |
| binary | gaussian | GLM | static | ATE | difference | sandwich | survey | ✅ truth | test-simulation.R |
| binary | gaussian | GLM | static | ATE | difference | bootstrap | survey | ✅ truth (ATE = 3) | test-simulation.R |
| binary | binomial | GLM | static | ATE | difference (RD) | sandwich | none | ✅ truth | test-simulation.R |
| binary | binomial | GLM | static | ATE | ratio (RR) | sandwich | none | ✅ truth | test-simulation.R |
| binary | binomial | GLM | static | ATE | OR | sandwich | none | ✅ truth | test-simulation.R |
| binary | binomial | GLM | static | ATE | difference | bootstrap | none | ✅ truth | test-simulation.R |
| binary | poisson | GLM | static | ATE | ratio | sandwich | none | ✅ truth | test-gcomp.R |
| binary | quasibinomial | GLM | static | ATE | difference | sandwich | none | ✅ truth (vs empirical E_L) | test-simulation.R |
| binary | gamma (log) | GLM | static | ATE | ratio | sandwich | none | ✅ truth (exp(β_A)) | test-simulation.R |
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
| multivariate | gaussian | GLM | scale_by | ATE | difference | sandwich | none | ✅ truth | test-multivariate.R |
| multivariate | gaussian | GLM | threshold | ATE | difference | sandwich | none | ✅ truth | test-multivariate.R |
| multivariate | gaussian | GLM | dynamic | ATE | difference | sandwich | none | ✅ truth | test-multivariate.R |
| multivariate | binomial | GLM | static | ATE | difference (RD) | sandwich | none | ✅ truth | test-multivariate.R |
| multivariate | binomial | GLM | static | ATE | ratio (RR) | sandwich | none | ✅ truth | test-multivariate.R |
| multivariate | binomial | GLM | static | ATE | OR | sandwich | none | ✅ truth | test-multivariate.R |
| multivariate | gaussian | GLM | static | subset | difference | sandwich | none | ✅ truth | test-multivariate.R |
| multivariate | gaussian | GLM | static | by(L) | difference | sandwich | none | ✅ truth | test-multivariate.R |
| multivariate | gaussian | GLM | static | ATE | difference | bootstrap | none | ✅ ratio vs sandwich | test-multivariate.R |
| binary | gaussian | GLM | static | ATE | difference | sandwich | none + censoring | ✅ truth | test-simulation.R |
| binary | gaussian | unsupported `model_fn` | static | ATE | difference | numeric Tier 1 | none | ✅ unit (vs analytic IF on logistic GLM; no e2e through `causat()` because the analytic gcomp branch wins for any standard fitter, so Tier 1 only fires in synthetic unit setup) | test-variance-if.R |
| binary | gaussian | unsupported `model_fn` | static | ATE | difference | numeric Tier 2 | none | ✅ truth (vs main path, to ~1%) | test-variance-if.R |
| categorical (k>2) | gaussian | GLM | static | ATE | difference | sandwich | none | ✅ truth | test-simulation.R |

**Error / rejection paths** (must also be tested):
- gcomp + invalid family string → ✅ test-gcomp.R
- gcomp + missing outcome col → ✅ test-causat.R
- gcomp + missing treatment col → ✅ test-causat.R

---

## Method 2 — Inverse probability weighting (`estimator = "ipw"`)

The Phase 4 chunk 3c rewrite landed the self-contained IPW engine: a
single density-ratio + weighted-MSM path (no WeightIt runtime) covering
static (HT indicator), shift / scale_by (continuous pushforward with
Jacobian), dynamic (binary HT indicator), and ipsi (closed-form weight).
Variance goes through one `variance_if_ipw()` loop calling
`compute_ipw_if_self_contained_one()` — no Branch A / Branch B split.
Chunk 3d added `WeightIt::glm_weightit()` contrast-level oracle tests
for the static-binary ATE / ATT / ATC cases (T-oracle1..3) and a
`mgcv::gam` propensity smoke test (T-oracle4); `WeightIt` moved from
`Imports:` to `Suggests:` in the same chunk. The oracle forced a
correctness fix to the ATT / ATC density-ratio weights: the HT indicator
branch now carries the Bayes-rule numerator `f*(L) = f(A* | L)` so ATT
uses `p(L)` and ATC uses `1 - p(L)`, making the per-arm Hájek sums
equal to the textbook conditional-mean weighted sums (Imbens 2004;
Hernán & Robins Ch. 12) and matching `glm_weightit` on the same PS
model to ~1e-6 on point estimates and SEs.

| Treatment | Outcome | Intervention | Estimand | Contrast | Variance | Weights | Status | Test file |
|---|---|---|---|---|---|---|---|---|
| binary | gaussian | static | ATE | difference | sandwich (HT + self-contained IF) | none | ✅ truth + WeightIt oracle (T-oracle1) | test-simulation.R, test-ipw.R, test-ipw-weights.R, test-ipw-weightit-oracle.R, test-treatment-model.R |
| binary | gaussian | static | ATE | difference | bootstrap | none | ✅ truth | test-simulation.R |
| binary | gaussian | static | ATT | difference | sandwich (HT, ATT Bayes numerator) | none | ✅ truth + WeightIt oracle (T-oracle2) | test-simulation.R, test-ipw-weightit-oracle.R |
| binary | gaussian | static | ATC | difference | sandwich (HT, ATC Bayes numerator) | none | ✅ truth + WeightIt oracle (T-oracle3) | test-simulation.R, test-ipw-weightit-oracle.R |
| binary | gaussian | static | ATE | difference | sandwich (mgcv::gam propensity) | none | ✅ sanity + WeightIt GLM-PS soft match (T-oracle4) | test-ipw.R, test-ipw-weightit-oracle.R |
| binary | gaussian | static | ATE | difference | sandwich | survey | ✅ truth | test-simulation.R |
| binary | binomial | static | ATE | difference (RD) | sandwich | none | ✅ truth | test-simulation.R |
| binary | binomial | static | ATE | ratio | sandwich | none | ✅ truth | test-simulation.R |
| binary | binomial | static | ATE | OR | sandwich | none | ✅ truth | test-simulation.R |
| binary | binomial | static | ATE | difference | bootstrap | none | ✅ truth | test-simulation.R |
| binary | gaussian | dynamic | ATE | difference | sandwich (HT indicator) | none | ✅ truth | test-simulation.R, test-ipw-weights.R |
| binary | gaussian | ipsi(δ) | ATE | difference | sandwich (closed-form weight) | none | ✅ truth | test-simulation.R, test-ipw-weights.R |
| continuous | gaussian | shift | ATE | difference | sandwich (pushforward ratio) | none | ✅ truth | test-simulation.R, test-ipw-weights.R |
| continuous | gaussian | scale_by | ATE | difference | sandwich (pushforward + Jacobian) | none | ✅ truth | test-simulation.R, test-ipw-weights.R |
| categorical (k>2) | gaussian | static | ATE | difference | sandwich | none | ✅ truth | test-simulation.R |
| continuous | gaussian | static (specific levels) | ATE | difference | sandwich | none | ✅ truth | test-simulation.R |
| Self-contained IPW sandwich variance (single path) | — | — | — | sandwich | ✅ truth (T-A_β_α hand-derived cross-derivative + T-end-to-end stacked-sandwich by hand) | test-ipw-branch-b.R, test-ipw-cross-derivative.R, test-variance-if.R |
| continuous | gaussian | static | — | — | — | — | ⛔ **rejected** (HT indicator degenerate on continuous treatment) | test-ipw.R |
| continuous | gaussian | threshold | — | — | — | — | ⛔ **rejected** (pushforward under clamp is a mixed measure; use gcomp) | test-ipw.R |
| continuous | gaussian | dynamic | — | — | — | — | ⛔ **rejected** (deterministic rule on continuous treatment is a Dirac pushforward; use gcomp) | test-ipw.R |
| any | any | shift / scale_by / dynamic / ipsi | ATT or ATC | any | any | any | ⛔ **rejected** — `check_estimand_intervention_compat()` aborts with `causatr_bad_estimand_intervention` | test-estimand-intervention-compat.R |
| any | any | any — with `A:modifier` in `confounders` | any | any | any | any | ⛔ **rejected** (Phase 8; saturated MSM cannot represent effect modification — use gcomp) | test-ipw.R |
| multivariate | any | any | any | any | any | any | ⛔ **rejected** (Phase 7) | test-s3-methods.R, test-multivariate.R |
| longitudinal (any shape) | any | any | any | any | any | any | ⛔ **rejected** (deferred — pooled-over-time MSM not yet implemented) | test-ipw.R |

---

## Method 3 — Matching (`estimator = "matching"`)

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
| any | any | any — with `A:modifier` in `confounders` | any | any | any | any | ⛔ **rejected** (Phase 8; saturated MSM cannot represent effect modification — use gcomp) | test-matching.R |
| multivariate | any | any | any | any | any | any | ⛔ **rejected** (Phase 7) | test-s3-methods.R |
| longitudinal | any | any | any | any | any | any | ⛔ **rejected** | test-matching.R |

---

## Method 4 — ICE longitudinal g-computation (`estimator = "gcomp"`, `type = "longitudinal"`)

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
| continuous | gaussian | threshold | 2 | sandwich | none | ✅ truth (2*E[max(A,0)]) | test-simulation.R |
| continuous | gaussian | shift | 2 | bootstrap | none | ✅ truth (2*δ) | test-simulation.R |
| binary | gaussian | static | 2 | sandwich | none + censoring | ✅ truth | test-ice.R |
| multivariate | gaussian | static | 2 | sandwich | none | ✅ truth | test-multivariate.R |
| multivariate | gaussian | shift | 2 | sandwich | none | ✅ truth | test-multivariate.R |
| binary | gaussian | ipsi | any | — | — | ⛔ **rejected** (Phase 4) | test-ice.R |
| binary | gaussian | static | 2 | bootstrap | survey | ✅ smoke (ratio vs sandwich) | test-simulation.R |
| binary | gaussian | static | 2 | sandwich | none + ATT | ⛔ **rejected** (only ATE for ICE) | test-ice.R |

---

## Method 5 — Survival (`causat_survival()`)

### Fit path (implemented today)

| Treatment | Outcome | Hazard model | Censoring | Competing | Status | Test file |
|---|---|---|---|---|---|---|
| binary | survival | pooled logistic + `ns(time)` | none | none | ✅ fit-only smoke | test-s3-methods.R, test-simulation.R |
| binary | survival | pooled logistic + `factor(time)` | none | none | ✅ fit-only smoke | test-s3-methods.R |
| binary | survival | pooled logistic | present (IPCW-style row filter) | none | ✅ fit-only smoke | test-simulation.R |

### Contrast path (Phase 6, not yet implemented)

| Quantity | Intervention | Variance | Status | Notes |
|---|---|---|---|---|
| Survival curve S^a(t) | static | sandwich | ❌ Phase 6 | `contrast(fit_surv, ...)` aborts today |
| Risk at time t (1 − S^a(t)) | static | sandwich | ❌ Phase 6 | |
| Risk difference | static | sandwich | ❌ Phase 6 | |
| Risk ratio | static | sandwich + log-scale delta | ❌ Phase 6 | |
| Any survival quantity | dynamic | sandwich | ❌ Phase 6 | |
| NHEFS 120-month RD replication (Ch. 17) | static | sandwich | ❌ Phase 6 | target ≈ 0.2%, 95% CI ≈ (−4.1%, 3.7%) |

### Rejection paths (tested today)

| Condition | Behaviour | Test file |
|---|---|---|
| `competing` argument non-`NULL` | ⛔ rejected (Phase 6) | test-causat.R |
| `contrast()` called on a `causatr_fit` with `type = "survival"` | ⛔ rejected (Phase 6) | test-contrast.R |

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
| `diagnose()` for gcomp / matching | ✅ smoke + snapshots | test-diagnose.R |
| `diagnose()` for IPW (self-contained runtime via `diagnose_ipw_point()` shim on binary static ATE) | 🟡 smoke + snapshots; full rewrite deferred to **Phase 9** (`PHASE_9_DIAGNOSE.md`) | test-diagnose.R |
| `diagnose()` intervention-aware (`intervention =` arg, per-intervention panels) | ❌ planned (**Phase 9**) | — |
| `diagnose()` for longitudinal | ⛔ **rejected** today; full longitudinal dispatch planned for **Phase 9** | test-diagnose.R |
| `diagnose()` treatment-type-aware dispatch via `detect_treatment_family()` (binary / categorical / continuous) | 🟡 shim in Phase 4; full dispatch in **Phase 9** | — |
| `diagnose()` estimand-aware balance (ATE vs ATT vs ATC reference populations) | ❌ planned (**Phase 9**) | — |
| `diagnose()` effect-modification-aware stratification (`by = "modifier"`) | ❌ planned (**Phase 9**, depends on Phase 8) | — |
| S3 methods: print, summary, plot, coef, vcov, confint, tidy, glance | ✅ smoke | test-s3-methods.R |
| `confint()` respects `level` arg (sandwich path) | ✅ smoke | test-s3-methods.R |
| `confint()` consistency between sandwich and bootstrap on `level` | ✅ unit + monotonicity | test-s3-methods.R |
| `causat()` rejects NA / Inf / negative / mis-sized / non-numeric weights | ✅ snapshot | test-causat.R |
| `contrast()` rejects duplicated intervention names | ✅ snapshot | test-contrast.R |
| `contrast()` rejects empty intervention list | ✅ snapshot | test-contrast.R |
| `causat_mice()` (multiple imputation wrapper) | 🟡 stubbed | test-causat-mice.R |
| Numeric variance Tier 1 (estfun-based fallback) | ✅ unit (analytic IF to ~1e-6) | test-variance-if.R |
| Numeric variance Tier 2 (delta shortcut + warn) | ✅ e2e via custom `model_fn` | test-variance-if.R |
| Cluster-robust matching variance (vs an IF-level reference) | ✅ unit (hand-computed sum-then-square) + e2e | test-variance-if.R, test-simulation.R |
| Intervention constructors (`static`, `shift`, `scale_by`, `threshold`, `dynamic`, `ipsi`) | ✅ unit (class + slot shape, reject paths) | test-interventions.R |
| External-reference cross-checks (`stdReg2`, hand-computed SE, `delicatessen` ICE values) | ✅ truth (tight tolerance against independent implementations) | test-variance-reference.R |
| FEATURE_COVERAGE_MATRIX.md ↔ tests/ round-trip audit | ✅ | test-coverage-matrix.R |

---

## Planned coverage (future phases)

These rows capture the target coverage for features the scaffold and
`PHASE_*.md` files commit to, but that are not yet implemented. They
are listed here so (a) the matrix reflects the full design intent, not
just what's shipped, and (b) a future PR that lands the implementation
has an explicit checklist to fill in.

### Phase 4 — Self-contained IPW (remaining work)

The Phase 4 chunk 3c runtime rewrite has shipped — the shipped rows now
live in the main IPW method table above. The rows below are the Phase 4
items still outstanding after 3c (categorical sub-chunk 3e, contrast-
level oracles 3d / 3f, non-static regression tests 3g, bootstrap parity,
and multivariate / longitudinal IPW).

| Treatment | Outcome | Intervention | Estimand | Variance | Status | Notes |
|---|---|---|---|---|---|---|
| categorical (k>2) | gaussian | static | ATE | sandwich | ❌ planned (chunk 3e) | multinomial density via `nnet::multinom` through `propensity_model_fn`; `fit_treatment_model()` currently aborts with `causatr_phase4_categorical_pending` |
| categorical (k>2) | gaussian | dynamic | ATE | sandwich | ❌ planned (chunk 3e) | HT indicator weight on deterministic rule |
| binary | gaussian | static | ATE × {ATE/ATT/ATC} | sandwich | ❌ planned (chunk 3d) | T-oracle1..3: contrast-level parity vs `WeightIt::glm_weightit()` |
| binary | gaussian | static | ATE | sandwich | ❌ planned (chunk 3d) | T-oracle4: `propensity_model_fn = mgcv::gam` |
| continuous | gaussian | shift | ATE | sandwich | ❌ planned (chunk 3f) | T-oracle5: `lmtp::lmtp_ipw()` point-estimate parity |
| binary | gaussian | ipsi(δ) | ATE | sandwich | ❌ planned (chunk 3f) | T-oracle6: manual IPSI weight parity |
| non-static IPW | — | shift / ipsi | ATE | sandwich | ❌ planned (chunk 3g) | T-non-static: full IF variance must materially exceed the `J V_β Jᵀ` delta-only shortcut |
| non-static IPW | — | shift / ipsi | ATE | bootstrap | ❌ planned (chunk 3g) | bootstrap parity: `ci_method = "bootstrap"` agrees with `ci_method = "sandwich"` within Monte Carlo error |
| multivariate | any | any | any | any | ❌ deferred (Phase 7) | IPW multivariate needs density-ratio engine with joint treatment density |
| longitudinal | any | any | any | any | ❌ deferred | extends density-ratio machinery to pooled-over-time MSM; explicitly out of Phase 4 scope |

### Phase 6 — Survival contrasts and competing risks

| Treatment | Outcome | Intervention | Quantity | Variance | Status | Notes |
|---|---|---|---|---|---|---|
| binary | survival | static | survival curve S^a(t) | sandwich | ❌ planned | per-individual hazard → cumulative product |
| binary | survival | static | risk at time t | sandwich | ❌ planned | 1 − S^a(t) |
| binary | survival | static | risk difference | sandwich | ❌ planned | (1 − S^{a1}) − (1 − S^{a0}) |
| binary | survival | static | risk ratio | sandwich | ❌ planned | delta method on cumulative product |
| binary | survival | static | any | bootstrap | ❌ planned | refit → cumulative product per replicate |
| any | survival | any | cause-specific CIF | sandwich | ❌ planned | competing risks |
| binary | survival | dynamic | survival curve | sandwich | ❌ planned | Ch. 17 dynamic strategies |
| longitudinal + ICE survival | any | static/dynamic | survival curve | sandwich | ❌ planned | backward iteration with hazard models |
| NHEFS 120-month survival replication (Ch. 17) | — | — | risk difference | sandwich | ❌ planned | target ≈ 0.2%, 95% CI ≈ (−4.1%, 3.7%) |

### Phase 7 — Advanced

| Feature | Status | Notes |
|---|---|---|
| Survey weights pass-through + adjusted sandwich | 🟡 partial | Current `weights` arg supports basic survey weights; needs explicit `survey` design integration |
| Cluster-robust sandwich beyond matching subclass | ❌ planned | `sandwich::vcovCL()` integration for cluster-resampled designs |
| Parallel bootstrap via `boot::boot(parallel=, ncpus=)` | ✅ done | tests in `test-ice.R`, `test-s3-methods.R` |
| Parallel bootstrap via `future` backend | ❌ planned | alternative to `boot::boot()` built-in parallelism |
| Multivariate treatment (g-comp, point + longitudinal) | ✅ done | test-multivariate.R |
| Multivariate treatment (IPW) | ❌ planned | Phase 4 — needs self-contained density-ratio engine |
| Multivariate treatment (matching) | ❌ planned | Phase 7 — `MatchIt` does not natively support joint treatments |
| Continuous treatment vignette | ❌ planned | companion to `vignettes/gcomp.qmd` |
| `target_trial()` metadata / specification helper | ❌ planned | Ch. 22 |
| Documentation: collider / LASSO confounder selection → `dagitty` | ❌ planned | Ch. 18 |
| Documentation: ML in g-formula → `lmtp` | ❌ planned | Ch. 18 |
| Sequential positivity warnings (longitudinal) | ❌ planned | Phase 7 (deferred from Phase 5) |
| Stratified ICE option (`stratified = TRUE`) | ❌ planned | Phase 7 (deferred from Phase 5) |
| Multinomial outcomes across methods | ❌ planned | `nnet::multinom()` or `VGAM::vglm(multinomial())` |
| Ordinal outcomes (cumulative probability contrasts) | ❌ planned | `MASS::polr()` / `ordinal::clm()` |
| Grace period / visit process interventions | ❌ planned | cf. `gfoRmula` visitprocess |
| "% intervened on" feasibility diagnostic | ❌ planned | cf. `gfoRmula` |
| ICE: auto-expand treatment × baseline-modifier interaction to all lag periods | ❌ planned | current behavior handles current-period A × modifier only; lagged treatments do not get auto-expanded modifier interactions. See the `~ A:sex` discussion in the longitudinal section — the fix is a formula-template expansion in `ice_build_formula()` |

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
