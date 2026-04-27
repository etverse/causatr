# Phase 11 — Full `diagnose()` rewrite

> **Status: IN PROGRESS** — chunk 11a shipped 2026-04-26.
> Dependencies: Phase 4 (done), Phase 5 (done), Phase 6 (effect modification), Phase 8 (multivariate IPW), Phase 10 (longitudinal IPW), Phase 14 (IPCW)

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

**Design (resolved in chunk 11a):**
- [x] Final shape of `causatr_diag`: nested `per_intervention` list, top-level
  `match_quality` and `fit_info` slots, backward-compat flat slots
  (`positivity`, `balance`, `weights`) shadow the first panel's contents.
- [x] `intervention =` argument default: `interventions = NULL` auto-injects a
  named `list(obs = ...)` so the print output stays readable.

**Core rewrite (in progress):**
- [x] `R/diagnose.R` chunk 11a — `interventions =` arg, per-intervention
  dispatch, seven `causatr_diag_*_pending` rejection gates, helpers split
  into `compute_weight_summary_observed()` / `compute_weight_summary_intervention()`
  / `summarise_weights_by_arm()`.
- [ ] Positivity + balance + weight helpers, one per treatment type
  (chunk 11b — `gaussian` / `categorical` / `count` / multivariate).
- [ ] Longitudinal dispatch path (chunk 11c).
- [ ] Estimand-aware ATT / ATC balance + `by =` stratification (chunk 11d).

**Tests:**
- [x] `test-diagnose.R` chunk 11a — per-intervention dispatch on binary IPW;
  manual `1/p` ESS reconstruction; one `expect_error(class = "causatr_diag_*_pending")`
  per pending gate; multi-panel `print()` output check.
- [ ] One test per (estimator × treatment type × intervention type × estimand)
  combination as each follow-up chunk lifts its rejection.

**Vignettes:**
- [ ] `vignettes/diagnostics.qmd` (chunk 11e).

**FEATURE_COVERAGE_MATRIX:**
- [x] Phase 11 row updated with chunk 11a status.
- [ ] Per-treatment-type / per-estimand cells fill in across chunks 11b–11d.

## Chunks

The rewrite is staged into five focused chunks so each lands with its own
plan / test / commit cycle. The chunk pattern mirrors phases 4 (chunks
4a–4j), 6 (6a–6d), 8 (8a–8e), and 10 (10a–10c).

### Chunk 11a — Foundation (DONE, 2026-04-26)

What shipped:

- New nested `causatr_diag` shape — top-level slots `per_intervention`
  (named list, key → `{positivity, balance, weights}`), `interventions`
  (character vector of keys), `match_quality` (top-level, intervention-
  agnostic), `fit_info` (`treatment_type`, `estimand`, `type`, `has_em`).
  Backward-compat top-level `positivity` / `balance` / `weights` slots
  point to the first panel so the pre-chunk-11a flat shape continues to
  work for every existing caller and test.
- `interventions =` argument on `diagnose()`, parallel to `contrast()`'s
  signature. Default `NULL` auto-injects a single `obs` panel using the
  observed-treatment Horvitz-Thompson view (`1/p` on treated, `1/(1-p)`
  on controls) — preserves the legacy default behaviour. Each user-named
  entry spawns its own per-intervention density-ratio weight summary via
  `compute_density_ratio_weights()`.
- Hard-pending classed errors for every sub-feature deferred to a later
  chunk: `causatr_diag_continuous_pending`,
  `causatr_diag_categorical_pending`, `causatr_diag_count_pending`,
  `causatr_diag_multivariate_pending`,
  `causatr_diag_longitudinal_pending`,
  `causatr_diag_estimand_pending` (IPW + ATT/ATC only — gcomp/matching
  pass through as before), and `causatr_diag_em_pending` (the `by =`
  argument is reserved at the signature level so future chunks lift the
  rejection without a signature break).
- S3 method updates: `print.causatr_diag()` dispatches on number of
  panels (single-panel layout matches pre-chunk-11a output verbatim;
  multi-panel layout sections by intervention key);
  `plot.causatr_diag()` and `summary.causatr_diag()` carry through the
  nested shape unchanged because the cobalt object is intervention-
  independent in chunk 11a.
- Test suite expanded from 12 to 24 test blocks, including per-
  intervention dispatch checks, manual `1/p` reconstruction sanity, and
  one `expect_error(class = "causatr_diag_*_pending")` per pending
  rejection gate.

What deferred to later chunks (with the rejection class that gates each
combination today):

| Sub-feature | Chunk | Pending class |
|---|---|---|
| Continuous IPW per-intervention diagnostics | 11b | `causatr_diag_continuous_pending` |
| Categorical IPW per-intervention diagnostics | 11b | `causatr_diag_categorical_pending` |
| Count (poisson / negbin) IPW | 11b | `causatr_diag_count_pending` |
| Multivariate IPW | 11b | `causatr_diag_multivariate_pending` |
| Longitudinal (ICE) per-period diagnostics | 11c | `causatr_diag_longitudinal_pending` |
| ATT / ATC balance for IPW | 11d | `causatr_diag_estimand_pending` |
| Effect modification (`by = "modifier"`) | 11d | `causatr_diag_em_pending` |
| Propensity histograms + weight-distribution plots | 11e | (no rejection — current Love plot still works) |
| `vignettes/diagnostics.qmd` | 11e | — |

Note: chunk 11a's balance view is the unadjusted SMD across treatment
groups — the same view the cobalt formula interface produces today. Per-
intervention post-weighting balance (which couples tightly to estimand
semantics) lands in chunk 11d alongside the ATT / ATC rewrite.

### Chunk 11b — Treatment-type dispatch (PENDING)

Lift `causatr_diag_continuous_pending`, `causatr_diag_categorical_pending`,
`causatr_diag_count_pending`, and `causatr_diag_multivariate_pending`.
Add per-treatment-type positivity and weight-distribution helpers:
density-range checks for gaussian, per-level distributions for
categorical, discrete-density quantiles for poisson / negbin, and per-
component summaries for multivariate. Reuse `detect_treatment_family()`
and the existing `compute_density_ratio_weights*` family.

### Chunk 11c — Longitudinal (ICE) dispatch (PENDING)

Lift `causatr_diag_longitudinal_pending`. Add per-period positivity,
balance, censoring-attrition, and weight-distribution diagnostics on the
ICE / longitudinal-IPW fits. Loop over `time_points` similar to
`variance_if_ice()`'s per-period iteration. Surface results under a
`time_points` nested list in the `causatr_diag` output.

### Chunk 11d — Estimand and effect-modification awareness (PENDING)

Lift `causatr_diag_estimand_pending` and `causatr_diag_em_pending`. Per-
intervention post-weighting balance using the standard ATE-IPW combined
weight `A/p + (1-A)/(1-p)` for static interventions, ATT-flavoured
within-treated balance using `A + (1-A) * p/(1-p)`, and the symmetric
ATC view. EM stratification via `cobalt::bal.tab(..., cluster = by)` or
manual subset loop. Both gate on baseline-modifier status (Phase 6's
known limitation under MSMs).

### Chunk 11e — Plot overhaul and vignette (PENDING)

Faceted propensity-score histograms (one panel per intervention),
faceted weight-distribution histograms with log-scale option, and a Love
plot grid faceted by intervention × modifier stratum. Write
`vignettes/diagnostics.qmd` as a user-facing tour of every supported
combination.

## Open questions

- Should `diagnose()` be pipe-friendly? (`fit |> diagnose() |> plot()`.) Current behavior is pipe-friendly for the simple case; confirm it survives the rewrite.
- Should per-intervention diagnostics be computed lazily or eagerly? Eager is simpler; lazy saves compute when the user only looks at one of many intervention panels. Probably eager for Phase 11 with an optional `keep = c("positivity", "balance", "weights")` arg for users who want to disable some.
- How does Phase 11 interact with the survey-weight / clustered-data work in Phase 9? Balance diagnostics under survey weights need the survey weights to enter the SMD computation; clustered diagnostics need cluster-robust balance. Both are natural extensions but should be scoped inside Phase 9, not Phase 11.
