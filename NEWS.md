# causatr (development version)

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
  silently ignore `A:modifier` terms; ICE handles current-period
  A × modifier but not lagged A × modifier. Phase 8 will fix all
  four in one pass.

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

- **`causat(method = "gcomp" / "ipw" / "matching")`** fully
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
