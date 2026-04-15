# causatr (development version)

## 2026-04-15 — Documentation overhaul: vignettes, math rendering, and `knit_print`

A pass over every package vignette plus a small package-level addition
so that `causatr_result` objects render as readable HTML tables inside
Quarto / knitr documents.

### New package feature: knitr-aware `causatr_result` rendering

- **`R/knit_print.R`** defines a `knit_print.causatr_result()` S3 method
  that renders a `causat() |> contrast()` result as a metadata header
  plus two tinytable HTML tables (intervention means + contrasts) when
  the result object is evaluated inside a knitr chunk. The method
  degrades gracefully to the existing `print.causatr_result()` if
  tinytable or knitr is unavailable.
- **`R/zzz.R`** conditionally registers the method via `.onLoad()`
  (`requireNamespace("knitr", quietly = TRUE)`), so the dependency on
  knitr stays Suggests-only.
- **`DESCRIPTION`** adds `tinytable` and `litedown` to `Suggests` (the
  latter is needed by `tinytable::format_tt(markdown = TRUE)`, which the
  variance-theory vignette uses for math-in-cells rendering).

### Vignette fixes (no behavioural changes to the engine)

- **`longitudinal.qmd` — Table 20.1 null-effect section rewritten.**
  The previous explanation claimed "individual effects of A_0 and A_1
  are both zero" and that naive regression gives non-zero estimates for
  both. In fact (i) the true *marginal counterfactual* means are all
  equal to 60, not the conditional effects, (ii) `lm(Y ~ A0 + A1 + L1)`
  gives `A1 = 0` exactly (in the DGP Y does not depend on A1), and
  (iii) the bias shows up in `A0` only, as the conditional-within-L1
  effect of −8. The section now shows *both* naive strategies
  (unadjusted and L1-adjusted) failing, demonstrating the null paradox
  as stated in Hernán & Robins Chapter 20.

- **`by = "sex"` examples fixed across three vignettes.**
  - `gcomp.qmd` now refits with an explicit `qsmk:sex` term in
    `confounders`, so the stratified estimates actually move (Male
    ≈ 3.61, Female ≈ 3.44).
  - `ipw.qmd` and `matching.qmd` replace the `by = "sex"` code chunks
    with a Phase 8 callout explaining that IPW and matching wrap a
    saturated `Y ~ A` MSM and so `by = "sex"` is a no-op until the
    unified effect-modification API lands
    (`PHASE_8_INTERACTIONS.md`).

- **Unified "Summary of covered combinations" tables.** gcomp, ipw,
  matching, and longitudinal now share one 8-column schema
  (Treatment | Outcome | Intervention | Estimand | Contrast | Inference
  | Weights | Status) with a common ✅ / 🟡 / ⛔ legend and a pointer
  to `FEATURE_COVERAGE_MATRIX.md`.

- **`code-fold: show` + `code-tools: true` in every vignette YAML.**
  Combined with the new `knit_print` hook, readers can fold the R code
  of any example chunk and still see a nicely formatted result table.

### `variance-theory.qmd` — roadmap + verified numerics

- New **Section 0 "Roadmap"** walks through Sections 1 → 6 as a single
  narrative and introduces the running example (a small logistic GLM)
  used by every verification block below.
- Eight new **Math / Code / Result** tabset panels, each of which
  instantiates a key equation on the running example and checks it
  against either `sandwich::sandwich()`, `sandwich::vcovCL()`, or
  `causatr::variance_if()`:
  1. **Sec 1.4** — Sandwich variance of β̂ vs `sandwich::sandwich(fit)`
     (agreement ≲ 1e-8).
  2. **Sec 2.3** — `(1/n²) Σ IFᵢ IFᵢᵀ` equals the sandwich to machine
     precision.
  3. **Sec 2.4** — IF for a sample mean reproduces the textbook
     $s/\sqrt n$ up to the $(n-1)/n$ d.f. correction.
  4. **Sec 3.2** — Two-channel IF for g-comp $\hat\mu_1$ matches
     causatr's sandwich SE to $\sim$1e-9.
  5. **Sec 3.3** — Delta-method decomposition: shows numerically that
     `J V_β Jᵀ` equals `Ch.2²`, and that dropping `Ch.1²` + cross term
     underestimates the SE by about 3% on this example.
  6. **Sec 4.2** — Horvitz-Thompson IPW $\hat\mu_1$ equals causatr's
     IPW point estimate, with SE from the propensity-corrected
     adjusted score.
  7. **Sec 4.3** — Cluster-robust sandwich: per-cluster IFs summed
     before squaring reproduce `sandwich::vcovCL(..., cadjust = FALSE)`
     on a toy subclass assignment.
  8. **Sec 5.5** — ICE chained IF on a 2-period linear-Gaussian DGP
     recovers the structural truth $0.8 + 0.6 = 1.4$ within the
     sandwich CI.

### `validation.qmd` — tinytable comparisons

- `compare_row()` now returns a tinytable with a third row of
  **absolute differences** between causatr and each reference
  implementation (stdReg2, WeightIt, marginaleffects, lmtp). Agreement
  — or the lack of it — is readable without having to compare 6
  decimal places by eye.

## Breaking: `causat()` / `causat_survival()` argument `method` renamed to `estimator`

`causat(method = ...)` shadowed `WeightIt::weightit(method = ...)` and
`MatchIt::matchit(method = ...)`, making it impossible to forward those
packages' own `method` selector through `...`. The causatr argument has
been renamed to `estimator`; extra `...` arguments now flow cleanly into
WeightIt / MatchIt. This is a hard rename with no deprecation shim —
update call sites from `method = "gcomp" | "ipw" | "matching"` to
`estimator = "..."`, and from `fit$method` / `result$method` /
`diag$method` to `fit$estimator` / `result$estimator` / `diag$estimator`.
The `glance.causatr_result()` output column is renamed in the same way.

After the rename, forwarding e.g. a CBPS propensity model is a one-liner:

```r
fit <- causat(data, outcome = "Y", treatment = "A", confounders = ~ L1 + L2,
              estimator = "ipw", method = "cbps")
```

See `vignette("ipw")` for a worked example. `ci_method` is a separate
concept and is unchanged.

## 2026-04-14 — Critical review sweep (R1–R9 + S1–S9)

Systematic audit and hardening of the whole package in a single pass.
All fixes below are covered by existing or newly-added tests; the
FEATURE_COVERAGE_MATRIX.md now carries a CI-enforced round-trip check
so future drift between the matrix and the test files fails loudly.

### Required changes (correctness / silent-failure fixes)

- **`contrast()` ratio / odds-ratio guard (R2).** The log-scale delta
  method was silently returning `NaN` CIs when `mu_hat` had mixed signs
  (legal for Gaussian outcomes, fatal for `log(R)`). `contrast(type =
  "ratio")` now aborts upfront if any intervention mean is non-positive,
  and `type = "or"` additionally aborts if any mean is ≥ 1. Error
  messages point to `type = "difference"` or a binomial / poisson /
  gamma refit.
- **`variance_bootstrap()` warning suppression narrowed (R3).** The
  refit closure previously wrapped the whole pipeline in
  `suppressWarnings()`, hiding GLM non-convergence, bread-singular
  warnings, and everything else. Replaced with a
  `withCallingHandlers()` that demotes only the three expected noisy
  warnings (`fitted probabilities numerically 0 or 1`, `X'WX is
  singular`, `Fewer (control|treated) units`). Genuine failures now
  surface through `tryCatch`'s error path and show up as failed
  replicates in `boot_info`.
- **Shared match-weight alignment helper (R4).** Extracted
  `combine_match_and_external_weights()` (in `R/matching.R`) so both
  `fit_matching()` and `refit_matching()` share the rowname-invariant
  check. A future change to `as.data.frame.data.table`'s rowname
  behavior now fails loudly in exactly one place instead of silently
  producing misaligned weights in the bootstrap.
- **`dynamic()` docstring corrected (R5).** The previous wording
  claimed the rule received data "subset to the current time point";
  in reality it is called once per intervention with the full
  counterfactual panel. Updated the roxygen prose and the example to
  show how to branch on the time column or a lag column.
- **ICE all-NA column drop is no longer silent (R6).**
  `ice_build_formula()` still drops `lag*_` columns at step 0 (a
  structural requirement), but now emits a `rlang::inform()` when a
  user-supplied column is dropped at a step because it happened to be
  all NA (e.g. a time-varying confounder first measured at time > 0).
  Gated via `.frequency = "regularly"` so the message does not spam
  long pipelines.
- **`boot_info$n_fail_by_int` (R7).** `process_boot_results()` now
  records per-intervention NA counts in addition to the whole-row
  totals. This lets users see whether bootstrap failures cluster on a
  single intervention (e.g. a rare `static()` level triggering
  separation) versus being spread across the whole replicate. The
  field is surfaced in `boot_info` and tested in `test-s3-methods.R`.
- **`causat_survival()` wide-format check (R8).** The previous
  `max(rows_per_id$N) == 1L` test only caught uniformly wide data; a
  mixed frame with some single-row ids and some multi-row ids would
  slip through. Now uses `any(rows_per_id$N == 1L)` and additionally
  asserts that `time` varies within each id (no duplicated rows at
  the same time value).
- **`FEATURE_COVERAGE_MATRIX.md` ↔ `tests/` round-trip check (R9).**
  New `tests/testthat/test-coverage-matrix.R` parses the matrix,
  extracts every `test-*.R` reference, and asserts both that the
  referenced file exists and that every test file on disk is
  referenced somewhere in the matrix. Guards against the exact kind
  of drift the rest of this sweep was cleaning up.

### Suggestions (hardening / clarity)

- **`fit_gcomp_point()` resolves `family` via `resolve_family()` (S1)**
  so custom `model_fn` implementations that don't auto-coerce a bare
  string see a fully-evaluated family object. Matches the behavior
  already in `fit_ipw()` and `fit_matching()`.
- **Multivariate sub-intervention name validation (S2).**
  `apply_intervention()` now asserts `setequal(names(iv), treatment)`
  before iterating the sub-interventions. A typo like
  `list(X = static(1), A2 = static(0))` for `treatment = c("A1",
  "A2")` now aborts with a clear message instead of silently creating
  a spurious `X` column via data.table.
- **`ipsi()` constructor informs about its Phase-4 status (S3).** The
  constructor succeeds (so users can build their interventions list)
  but surfaces a one-per-session `rlang::inform()` pointing to the
  known `contrast()` abort, so users hit the wall at construction
  rather than deep in a pipeline. Use `.frequency_id =
  "causatr_ipsi_dead_end"` to silence or grep for it.
- **`check_estimand_trt_compat()` error message (S4)** now explicitly
  says "binary point treatments coded as 0/1" and suggests recoding.
  Factor- or `1/2`-coded treatments previously got the same message
  without any hint at the coding issue.
- **Tier 2 numeric variance warning is now classed (S5).**
  `variance_if_numeric()`'s delta-shortcut warning carries
  `class = "causatr_tier2_fallback"` so batch pipelines and CI can
  grep for it via `withCallingHandlers()`.
- **`compute_weight_summary()` binary split pulls levels from
  WeightIt (S6)** instead of hard-coding `trt == 0` / `trt == 1`.
  Factor-coded binary treatments and integer codings like `c(1, 2)`
  now get correctly-populated "treated" / "control" rows.
- **`contrast()` subset evaluation preserves lexical scope (S7).** The
  three `eval(subset, envir = as.list(data))` call sites were copying
  every column of `data` into a fresh list and losing the caller's
  enclosing environment. Switched to
  `eval(subset, envir = data, enclos = parent.frame())`, which skips
  the copy and lets `subset` expressions reference session variables
  like `quote(age > cutoff)`.
- **FEATURE_COVERAGE_MATRIX survival section split (S8).** Replaced
  the mixed "fit-only smoke (Phase 6 contrast pending)" rows with an
  explicit three-table layout (fit path / contrast path / rejection
  path) so the Phase 6 gap is visible without status annotations that
  mean different things on different rows.
- **`.pseudo_y` column collision guard in ICE (S9).** `fit_ice()` now
  aborts upfront if `data` already has a column named `.pseudo_y`,
  which is the internal bookkeeping column ICE writes predicted
  pseudo-outcomes into during the backward iteration.

### Doc drift — `A:modifier` behaviour (R1)

`CLAUDE.md`, `NEWS.md` (this file's older entries), `PHASE_8_INTERACTIONS.md`,
and `CAUSATR_SCAFFOLD.md` all still described IPW and matching as
"silently ignoring" `A:modifier` interaction terms in `confounders`.
That code path was hardened months ago: `fit_ipw()` and `fit_matching()`
now call `check_confounders_no_treatment()` and **abort** upfront with
a Phase-8 pointer. Every doc that mis-described the old behavior has
been updated to match the current behavior.

### Vignette expansion

Every method vignette (excluding `variance-theory.qmd`) now covers
every combination supported by `FEATURE_COVERAGE_MATRIX.md`:

- **`gcomp.qmd`** — added worked examples for categorical (k > 2)
  treatment, multivariate (joint) treatment, external survey weights,
  and a Poisson count outcome with a log-link ratio contrast. Summary
  table now records the family / weights columns.
- **`ipw.qmd`** — added categorical and continuous (GPS) treatment
  sections, an external-survey-weights section, a parallel-bootstrap
  section, and a non-Mparts WeightIt fallback example. Updated the
  Channel-2 technical note to reference the current
  `prepare_propensity_if()` / `apply_propensity_correction()`
  architecture (was still referring to the renamed
  `correct_propensity()`).
- **`matching.qmd`** — added an ATT + external survey weights section
  demonstrating `combine_match_and_external_weights()`, and extended
  the summary table with a categorical / continuous rejection row.
- **`longitudinal.qmd`** — added continuous-treatment + MTP section
  (shift / scale_by / threshold), multivariate time-varying
  treatment, external weights (with a note pointing to the ICE
  cascade-gradient correction in §5.4 of the theory vignette), and a
  censoring / IPCW-style example.
- **`survival.qmd`** — replaced the 3-line TODO stub with a full
  narrative of the current fit-only Phase 6 status, the intended
  end-to-end contrast algorithm, and the tracking pointers.
- **`triangulation.qmd`** — introduction.qmd and the comparison table
  both corrected for the drift above (matching is binary-only;
  continuous / categorical are not "Limited").

## 2026-04-14 — Sweep and audit

- **ICE sandwich fix under non-uniform external weights.** The
  cascade gradient at step k > 0 was dropping the `w_{k-1}` factor
  from `A_{k-1,k} = E[∂s_{k-1}/∂β_k]`, silently underestimating the
  SE by ~2x on multi-period DGPs. Verified against `lmtp::lmtp_tmle`
  and a 300-replicate Monte-Carlo. (Followed by a separate fix for
  the earlier regression where the ICE backward loop dropped
  external weights entirely.)
- **ICE double-counting fix for `shift`/`scale_by`/`threshold`.**
  `ice_apply_intervention_long()` was recomputing lag columns from
  the intervened treatment, double-counting the shift through both
  the lag path and the current-A prediction path. Now matches
  `lmtp` on the canonical linear-Gaussian DGP.
- **New `FEATURE_COVERAGE_MATRIX.md`.** Single source of truth for
  what's tested and at what fidelity (truth / smoke / none /
  rejected). New test-writing rule in `CLAUDE.md`: every feature
  change MUST update the matrix in the same PR.
- **Coverage sweep** closing most residual gaps:
  - gcomp × quasibinomial, Gamma(log), probit, cloglog, Poisson ×
    bootstrap, categorical (3-level) treatment
  - IPW × categorical (3-level), continuous, survey weights, ATT
    bootstrap, non-Mparts WeightIt method at contrast time
  - matching × Poisson, Gamma(log), ATC bootstrap, categorical
    rejection (MatchIt is binary-only; verified against upstream
    docs)
  - ICE × `scale_by`, `threshold`, `shift` bootstrap, bootstrap ×
    survey weights, binomial × ratio/OR, `by(sex)` stratification,
    external-weights truth test
  - Survival × censored fit (smoke; contrasts still Phase 6)
  - End-to-end Tier 2 numeric variance via a custom `model_fn`
  - `vcov_from_if(cluster = ...)` unit test vs hand-computed
    reference
  - `confint()` level monotonicity on both sandwich and bootstrap
- **Defensive guards across `causat()`, `causat_survival()`, IPW,
  matching, ICE, diagnose.** Front-door validation for NA / Inf /
  negative / mis-sized / non-numeric `weights`; rejection of
  duplicated or empty intervention names; explicit longitudinal
  rejection in `diagnose()` (previously crashed inside cobalt);
  categorical treatment rejection in `matching` with a clearer
  causatr-side message.
- **`~ sex + A:sex` no longer triggers a false-positive** "Baseline
  confounder 'A' varies within some individuals" warning. Treatment
  columns are now excluded from `warn_confounder_variation()`.
- **altdoc sidebar fix.** Single-item sections ("Getting started",
  "Theory") were silently dropped from the deployed site because
  `yaml::yaml.load()` collapses a one-element scalar sequence into a
  scalar that Quarto's sidebar can't parse. Both sections now use
  the explicit `text:`/`file:` form, making the `variance-theory`
  vignette visible on the website.
- **`vignettes/variance-theory.qmd` Section 5.4** now derives the
  ICE cascade gradient with the explicit `w_{k-1}` factor and the
  pseudo-code shows the `prior.weights` lookup.
- **`PHASE_8_INTERACTIONS.md` design doc.** Plans a unified
  effect-modification API across gcomp / IPW / matching / ICE. IPW
  and matching currently hardcode a saturated `Y ~ A` MSM and
  **hard-abort** on any `A:modifier` term in `confounders` via
  `check_confounders_no_treatment()` (rather than silently dropping
  it and returning a pooled ATE); ICE handles current-period A ×
  modifier but not lagged A × modifier. Phase 8 will replace both
  the abort and the lag-gap with a proper MSM builder.

## 2026-04-13 — Variance engine unification

- **Single influence-function engine.** Replaced four method-specific
  variance paths with one `variance_if()` dispatcher that routes to
  `variance_if_gcomp()` / `variance_if_ipw()` /
  `variance_if_matching()` / `variance_if_ice()`. All four share
  `prepare_model_if()` / `apply_model_correction()` primitives and
  `vcov_from_if()` aggregator.
- **Prep / apply split.** `prepare_model_if()` computes the
  gradient-independent half of Channel 2 once per model;
  `apply_model_correction()` finishes it per intervention. Saves an
  O(k) bread solve when `contrast()` is called with many
  interventions.
- **Matching uses `cluster = subclass`** in the aggregator
  (cluster-robust IF-level sum-then-square).

## 2026-04-12 — Large refactor + feature pass

- **Point-treatment sandwich now uses the full influence function**
  (previous version was effectively a delta-method shortcut and
  could underestimate the SE on non-saturated MSMs).
- **`variance_if()` as single entry point.** Introduces
  `correct_model()` and `correct_propensity()` as the two per-model
  correction kernels. New `vignettes/variance-theory.qmd` derives
  the two-channel IF from first principles.
- **Log-scale CIs for `ratio` and `or` contrasts.** Respects the
  (0, ∞) support and gives better coverage than Wald bounds.
- **`static()` accepts character values** for categorical
  treatments (e.g. `static("never")`).
- **Analytical Jacobian for GLMs** in sandwich variance (previously
  used numerical differencing).
- **External weights moved to `fit$details$weights`** (never a
  column in `data`) so every method consumes them uniformly.
- **New `test-complex-dgp.R`** covering nonlinear confounding,
  GAMs, splines, and multi-confounder designs.
- **ICE IF scaling fix** for subset / by-stratified targets, plus a
  `confint()` ordering fix and a vectorized IF loop for speed.

## 2026-04-11 — Bootstrap percentile CIs

- **Bootstrap path stores raw replicates** in `object$boot_t` and
  `confint()` recomputes percentile CIs from them per `level`. The
  sandwich path continues to reconstruct Wald CIs from the stored
  SE.
- External weights now propagate through every method (gcomp, IPW,
  matching, ICE); previously ICE dropped them in the backward loop.

## 2026-04-10 — Polish and hardening

- **`by`-stratification `vcov` fix** — stratum-specific vcov
  matrices now land in `res$vcov[[level]]` as a list instead of
  being collapsed to the pooled matrix.
- **`family` now threads through IPW and matching** (previously
  silently forced Gaussian on the MSM regardless of the user's
  `family` argument).
- **`scale()` → `scale_by()` rename** to avoid colliding with the
  base R `scale()`.
- **`by + estimand` bug fix** — previously the ATT/ATC restriction
  and the by-level restriction fought over the target population
  and returned the wrong per-stratum estimate.
- **`type` parameter on `causat()`** for explicit point /
  longitudinal selection; autodetects from `id` + `time` otherwise.

## 2026-04-09 — Effect modification, S3 polish, broom integration

- **New `by` argument on `contrast()`** for stratified (by-level)
  effect estimates. Computes separate target means and vcov per
  stratum; `confint()` and `tidy()` respect the stratification.
- **broom integration** via `generics` — `tidy.causatr_result()`
  and `glance.causatr_result()` return the expected columns.
- **Forest plots** via the `forrest` package for `causatr_result`
  objects.
- **S3 method polish** across `print`, `summary`, `plot`, `coef`,
  `vcov`, `confint`, `tidy`, `glance`.

## 2026-04-02 / 2026-04-03 — Longitudinal scaffolding + R CMD check

- ICE sandwich variance with scaling fix; first batch of
  simulation-based ICE tests.
- Clean R CMD check: man pages regenerated, MASS / knitr declared,
  Windows vignette skip to work around a toolchain issue.

## 2026-04-01 — Phase 2–3 implementation

- **`causat(estimator = "gcomp" / "ipw" / "matching")`** fully
  working, with:
  - parametric g-formula via a user-supplied `model_fn`
  - IPW via `WeightIt::weightit()` and `WeightIt::glm_weightit()`
  - matching via `MatchIt::matchit()` + weighted MSM
- **`diagnose()`** for all three methods: positivity, balance (via
  `cobalt`), weight summaries, match quality, Love plots.
- **First vignettes**: `gcomp.qmd`, `ipw.qmd`, `matching.qmd`,
  `triangulation.qmd`.

## 2026-03-30 / 2026-03-31 — Phase 1 scaffolding

- Initial package structure: `DESCRIPTION`, `NAMESPACE`, `R/`
  layout, testthat setup, dev tooling, GitHub Actions CI.
- Core API stubs: `causat()`, `contrast()`, `diagnose()`,
  intervention constructors (`static`, `shift`, `scale_by`,
  `threshold`, `dynamic`, `ipsi`), S3 classes.
- NHEFS dataset bundled; input-validation helpers in `checks.R`.
