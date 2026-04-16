# Phase 11 — Full `diagnose()` rewrite

> **Status: PENDING** (design doc)
> Dependencies: Phase 4 (done), Phase 5 (done), Phase 6 (effect modification), Phase 8 (multivariate IPW), Phase 10 (longitudinal IPW)

## Why this deserves its own phase

`diagnose()` is the user-facing diagnostic entry point for every causatr fit. It currently (Phase 3) handles positivity, covariate balance, and weight-distribution summaries for the three cross-sectional estimators (gcomp, IPW via WeightIt, matching via MatchIt). It does **not** currently handle:

- Longitudinal fits (ICE) — `diagnose()` on a longitudinal fit aborts with "deferred to a future phase".
- Per-intervention diagnostics — the current implementation assumes a fit-time weight vector exists (true for Phase 3 IPW via WeightIt) and shows that one vector regardless of which intervention the user will eventually hand to `contrast()`. Under Phase 4's deferred-MSM architecture, weights are per-intervention, so the "which weights to show" question becomes real.
- Treatment-type-aware diagnostics — the current implementation reads `fit$weights_obj$treat` and `attr(..., "treat.type")` from WeightIt to dispatch between binary / categorical / continuous summaries. Under Phase 4 this dispatch lives in `R/treatment_model.R:detect_treatment_family()` and has to be re-wired.
- Estimand-aware diagnostics — balance diagnostics under ATT vs ATE are fundamentally different (standardized mean difference within the treated vs across the full pseudo-population). The current implementation partially handles this via `cobalt::bal.tab()`'s `estimand = ` argument on the `weightit` object, which we'll lose when WeightIt moves to Suggests.
- Effect-modification-aware diagnostics (Phase 6 interaction) — when the outcome model carries `A:modifier` terms, balance should be reported within each modifier stratum, not just overall.

All of these are real user-facing needs. Lumping them into the Phase 4 IPW engine rewrite risks under-scoping the diagnostic story: we'd either ship a Phase 4 `diagnose()` that silently works differently than Phase 3 users expect, or bolt on a dozen conditional branches for each new feature as it lands. Cleaner to call out the full rewrite as its own phase and commit to a comprehensive redesign.

## What Phase 4 DOES ship for `diagnose()`

Phase 4 (Chunk 3c) delivers a **minimal shim** so `diagnose()` keeps working on the common cases while the full rewrite is planned:

- For `estimator = "ipw"` fits under the self-contained engine:
  - Reads `fit$details$propensity_model` (the treatment density model from `fit_treatment_model()`) instead of `fit$weights_obj$ps`.
  - Builds a default weight vector at diagnose time by calling `compute_density_ratio_weights(tm, data, static(1))` for binary-treatment fits, or `compute_density_ratio_weights(tm, data, NULL)` (natural course weights = 1) for non-binary fits.
  - Dispatches between binary / categorical / continuous via `detect_treatment_family()` rather than WeightIt's `treat.type`.
  - Accepts an optional `intervention = ` argument so users can request diagnostics under a specific intervention; defaults to the binary/natural-course choice above.
- For `estimator = "gcomp"` and `estimator = "matching"` fits: no changes (still uses `fit$model` / `fit$match_obj` / `fit$weights_obj` for matching).
- Balance diagnostics still go through `cobalt::bal.tab()` but the weight vector passed in comes from the new pipeline, not from a `weightit` object.

The shim is explicitly marked as incomplete in the `diagnose()` docstring. Tests in `test-diagnose.R` are updated for the new slot layout but their coverage stays the same — binary static ATE on a cross-sectional fit, the bread-and-butter case.

## Phase 11 full rewrite scope

### 1. Intervention-aware diagnostics

Currently `diagnose()` takes `(fit)`. Under Phase 11 it takes `(fit, interventions = NULL)`, mirroring `contrast()`'s signature:

```r
diag <- diagnose(fit, interventions = list(a1 = static(1), a0 = static(0)))
```

Each intervention gets its own diagnostic panel: propensity distribution under that intervention's reweighting, weight summary (min / max / 99th percentile / effective sample size), covariate balance (standardized mean differences before and after weighting), positivity check (fraction of rows with extreme weights).

For natural-course (`NULL` or not supplied), the diagnostic reports the "raw" fit-time weights (1 for IPW, match weights for matching, unweighted for gcomp).

### 2. Cross-sectional vs longitudinal

For longitudinal (ICE) fits, per-time-point diagnostics:
- Positivity check per time step (per-period weight distribution).
- Covariate balance at each time step after weighting.
- Censoring diagnostics (fraction lost by each step).
- Visit-process diagnostics (grace-period violations if applicable).

The longitudinal diagnostic is inherently bigger than a single plot or table — probably a gt/tinytable report with a collapsible section per time step.

### 3. Treatment-type dispatch

Binary:
- Positivity: histogram of `p_i = Pr[A=1|L_i]` on [0, 1] with tails highlighted.
- Balance: standardized mean differences of confounders across treatment groups, before and after weighting. `cobalt` is the right backend.
- Weight distribution: histogram of `w_i = f(d(A_i,L_i)|L_i) / f(A_i|L_i)` (the density ratio), with log-scale option and extreme-weight count.

Categorical (k > 2):
- Positivity: per-level `Pr[A=k|L]` distributions.
- Balance: pairwise standardized mean differences or Love plots faceted by level.
- Weight distribution: per-level.

Continuous:
- Positivity: range of `f(A_i|L_i)` with low-density tail identification.
- Balance: correlation between treatment and confounders (Spearman or Pearson), before and after weighting.
- Weight distribution: histogram of density-ratio weights with truncation diagnostics.

Each dispatch should share infrastructure where possible (e.g., the effective-sample-size formula `ESS = (sum w)^2 / sum(w^2)` is common to all three).

### 4. Estimand-aware diagnostics

- **ATE**: balance across the full pseudo-population.
- **ATT**: balance within the treated, with controls reweighted to match the treated distribution. The "balance before" panel should report the observational SMDs; the "balance after" panel should report the reweighted SMDs relative to the treated marginals, not the pooled marginals.
- **ATC**: symmetric statement for the untreated.

The rejection gates from Phase 4 `check_estimand_intervention_compat()` carry over: ATT/ATC are only defined for static interventions under IPW/matching, and `diagnose()` should respect that gating (refuse to produce ATT diagnostics for a shift intervention).

### 5. Effect-modification awareness

When the outcome model carries interaction terms (Phase 6: `confounders = ~ L + sex + A:sex`), `diagnose()` should optionally report balance **within each modifier stratum**. Surface-level invocation:

```r
diagnose(fit, interventions = list(a1 = static(1)), by = "sex")
```

Under the hood: `cobalt::bal.tab(..., cluster = "sex")` or equivalent.

### 6. Output shape

Phase 3's `diagnose()` returns a `causatr_diag` object with `$positivity`, `$balance`, `$weights` slots. Phase 11 generalizes this to a nested structure:

```
causatr_diag
├─ fit_info: treatment type, estimator, estimand, interventions
├─ positivity
│    ├─ intervention_1: { density summary, extreme counts }
│    ├─ intervention_2: ...
├─ balance
│    ├─ intervention_1: { smd_before, smd_after, [stratified by modifier] }
│    ├─ intervention_2: ...
├─ weights
│    ├─ intervention_1: { min, max, quartiles, ess, extreme_mask }
│    ├─ intervention_2: ...
└─ longitudinal (optional)
     ├─ time_0: ...
     ├─ time_1: ...
```

Print / summary / plot methods dispatch appropriately.

### 7. Plot methods

Phase 3's `plot.causatr_diag()` produces a Love plot via `cobalt` when balance data is present. Phase 11 extends this to:

- Propensity-score histogram (facet by intervention).
- Weight-distribution histogram (log-scale option).
- Love plot (facet by intervention and optionally by modifier stratum).
- Longitudinal: per-time-step panels for each of the above.

### 8. Backwards compatibility

The Phase 4 shim keeps `diagnose(fit)` working on the common binary static ATE case. Phase 11's `diagnose(fit)` (no `interventions = ` argument) should produce the same output it does today — the new signature is additive, not a breaking change. Tests in `test-diagnose.R` from Phase 4 survive into Phase 11 as regression anchors.

## Items (to be landed in Phase 11)

**Design:**
- [ ] Decide on the final shape of `causatr_diag` (probably nested list of named sub-diagnostics per intervention, with a print method that dispatches on presence/absence of longitudinal / modifier panels).
- [ ] Decide on the `intervention = ` argument default: `list(NULL = NULL)` (natural course) or `list(obs = NULL)` (named natural course). The latter is better for the print output.

**Core rewrite:**
- [ ] `R/diagnose.R` — full rewrite. Remove all `fit$weights_obj` accesses; route everything through `fit$details$propensity_model`, `fit$model`, `fit$match_obj`, `fit$type`, `fit$family`, `fit$estimand`.
- [ ] Positivity + balance + weight helpers, one per treatment type (bernoulli / gaussian / categorical), each intervention-aware.
- [ ] Longitudinal dispatch path (`variance_if_ice`-style loop over time points).

**Tests:**
- [ ] Extend `test-diagnose.R` with one test per (estimator × treatment type × intervention type × estimand) combination that's supported.
- [ ] `cobalt` oracle: when `cobalt` is installed, balance numbers should match `cobalt::bal.tab()` output to high precision on the cross-sectional static binary case (this is what Phase 3 already does, carried forward).

**Vignettes:**
- [ ] `vignettes/diagnostics.qmd` — user-facing tour of all diagnose outputs across estimators / treatments / estimands.

**FEATURE_COVERAGE_MATRIX:**
- [ ] New section "diagnose() coverage" enumerating which (fit shape × intervention × estimand) combinations are covered, matching the main feature matrix's granularity.

## Open questions

- Should `diagnose()` be pipe-friendly? (`fit |> diagnose() |> plot()`.) Current behavior is pipe-friendly for the simple case; confirm it survives the rewrite.
- Should per-intervention diagnostics be computed lazily or eagerly? Eager is simpler; lazy saves compute when the user only looks at one of many intervention panels. Probably eager for Phase 11 with an optional `keep = c("positivity", "balance", "weights")` arg for users who want to disable some.
- How does Phase 11 interact with the survey-weight / clustered-data work in Phase 9? Balance diagnostics under survey weights need the survey weights to enter the SMD computation; clustered diagnostics need cluster-robust balance. Both are natural extensions but should be scoped inside Phase 9, not Phase 11.
