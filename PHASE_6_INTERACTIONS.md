# Phase 6 — Unified Effect-Modification API

> **Status: PENDING (design doc)**
> Book chapters: 13 (standardization + effect modification), 15 (matching-based effect modification), 21 (ICE effect modification across time)

## Scope

Unify how causatr handles effect modification (EM) — i.e. interactions between the treatment and a baseline or time-varying variable — across all estimation methods. Today each method does something different; two of them (IPW and matching) don't handle EM at all, and ICE only handles EM at the current-period treatment slot.

This phase does **not** introduce new interventions, methods, or outcome families. It's entirely about getting the existing methods to correctly surface effect modification via a consistent user-facing API.

## Today's state (pre-Phase 6)

The test matrix in `FEATURE_COVERAGE_MATRIX.md` currently marks `by(L)` as "✅ truth" only for point gcomp. The other cells either hard-abort (IPW / matching) or return partially correct numbers (ICE).

| Method | How EM is currently supported | What happens when a user writes `~ L + sex + A:sex` |
|---|---|---|
| **Point gcomp** | ✅ Outcome model includes the full `confounders` RHS. `A:sex` feeds the outcome GLM directly. `by(sex)` averages the predictions per stratum and recovers the correct effect. | Works. MC-verified: sex-specific contrasts recovered to ~1% on a linear-Gaussian DGP. |
| **Point IPW** | ⛔ `fit_ipw()` builds `msm_formula <- stats::reformulate(treatment, response = outcome)` — a hardcoded saturated `Y ~ A`. `confounders` is used for the propensity model only. `fit_ipw()` calls `check_confounders_no_treatment()` and **aborts** before any weights are estimated. | Hard error with a Phase-8 pointer. Preferred to silently returning a pooled ATE. |
| **Point matching** | ⛔ Same MSM shape as IPW: `msm_formula <- stats::reformulate(treatment, response = outcome)`. `fit_matching()` also calls `check_confounders_no_treatment()` and aborts upfront. | Hard error with a Phase-8 pointer. |
| **Longitudinal ICE** | ⚠️ `ice_build_formula()` injects `baseline_terms` verbatim at every period. `A:sex` resolves to `A_k:sex` — the **current-period** treatment × modifier. Lagged treatments (`lag1_A`, `lag2_A`, ...) do NOT get a corresponding `lag1_A:sex` term. | Partial. Current-period effect modification is captured; earlier-period contributions are collapsed. MC evidence: 2-period DGP with persistent `(1 + 1.5·sex)·A_k` effect returns contrasts 2.81 / 4.43 vs structural truth 2.10 / 5.10 — roughly half the intended heterogeneity. |

Phase 6 will replace the hard-abort path with a proper MSM builder that honors `A:modifier` terms for IPW and matching, and will auto-expand `lag{k}_A:modifier` for ICE.

The root cause in each case is different, which is part of why this needs a dedicated phase rather than four scattered fixes.

## Design goals

1. **One formula convention.** Users write `~ L + sex + A:sex` and all four methods do the right thing. The convention should be: a term involving the treatment (e.g. `A:sex`, `A:L`, `A:I(age > 65)`) is interpreted as effect modification of A on the outcome scale.
2. **No silent failures.** If a user writes an EM term and the chosen method genuinely can't support it (e.g. IPW + multivariate treatment + continuous modifier in Phase 4), abort with a specific error that names the method, the term, and the workaround.
3. **Leave (B) and (C) alone.** Interactions *among confounders* (e.g. `L1:L2`) and *within the current period* (e.g. `A:L` where L is TV) already work across methods. This phase only touches the `A × baseline` (and its time-varying analogs) case.
4. **No breaking changes to existing tests.** Every truth-test currently in the matrix must still pass after the refactor. Correct behavior under additive models collapses to the current behavior by construction.

## Plan

### 1. Introduce a shared helper for parsing effect-modification terms

In `R/checks.R` (or a new `R/effect_modification.R`), add:

```r
parse_effect_mod <- function(confounders, treatment) {
  # Walk the term.labels of `confounders`. For each term that
  # involves any variable in `treatment`, return a list with
  #   $term: the full term label (e.g. "A:sex")
  #   $treatment_var: the matching treatment variable
  #   $modifier_vars: the other variables in the interaction
  # Terms without any treatment variable are left alone.
}
```

This is the canonical detector every method will consult. Anchoring it here (rather than ad-hoc regex in each fitter) keeps the convention consistent and makes future extensions (e.g. three-way interactions) mechanical.

### 2. Fix point IPW and matching

The MSM formula in `fit_ipw()` and `fit_matching()` currently ignores everything except the treatment. Change both so:

- If `parse_effect_mod(confounders, treatment)` returns any EM terms, build the MSM formula as `reformulate(c(treatment, EM terms, modifier main effects), response = outcome)` — e.g. `Y ~ A + sex + A:sex`.
- Otherwise, keep the saturated `Y ~ A` (no regression in shape; this is what every existing test expects).

The weighted outcome fit (`WeightIt::glm_weightit()` for IPW, plain weighted `stats::glm()` for matching) already supports arbitrary MSM formulas; the only change is the formula construction itself.

Tests to add:

- IPW × binary × binary-modifier × gcomp-matched truth test (shares DGP with the existing point gcomp × by test).
- Matching × binary × binary-modifier × gcomp-matched truth test.
- IPW × saturated MSM sanity check — a DGP without EM terms must still give the exact result the current tests pin.

### 3. Fix ICE longitudinal

Extend `ice_build_formula()` in `R/ice.R` so that when a baseline term of the form `treatment:modifier` appears, it auto-expands to include the same interaction with every currently-available lag:

```r
# At time_idx = 2, max_lag = 2, treatment = "A", term = "A:sex":
#   emit "A:sex", "lag1_A:sex", "lag2_A:sex"
```

The expansion is per-period (later periods have more lags) and defaults to "this interaction applies uniformly across time". This handles the **time-invariant effect modifier** semantics. The rarer **time-varying effect modifier** semantics (different functional form per period) remains out of scope until a per-period formula DSL lands — document the limitation and point users at wide-format + point gcomp in the meantime.

Tests to add:

- ICE × 2-period DGP × `(1 + γ·sex)·A_k` × `by(sex)` — must recover the stratum-specific contrast to ~5% of MC truth (vs the current ~30% compression).
- ICE × DGP without EM — additive-model case must give identical numbers before and after the refactor (regression guard).
- ICE × `A:sex + A:age` — multiple EM terms stacked; auto-expansion must handle both.

### 4. Unify `by()` stratification across methods

Once (2) and (3) are in place, the `by` branch in `compute_contrast()` just works — it already averages predictions per stratum, and after the MSM fix the predictions actually depend on the stratum. No changes to `by` itself.

Add a cross-method truth test:

- gcomp, IPW, matching, ICE all run on the same EM DGP with the same formula; their stratum-specific contrasts agree to within the usual cross-method tolerance.

### 5. Documentation

- New vignette section in `vignettes/gcomp.qmd`, `vignettes/ipw.qmd`, `vignettes/matching.qmd`, and `vignettes/longitudinal.qmd` showing how to specify effect modification for each method.
- Expand `vignettes/triangulation.qmd` with an EM example.
- Rewrite the "Known limitations" note in `NEWS.md` to reflect the fix.
- Update `FEATURE_COVERAGE_MATRIX.md`: replace the three ❌/⚠️ EM rows with ✅ truth and add a unified cross-method EM row.

## Out of scope for Phase 6

| Topic | Reason | Deferred to |
|---|---|---|
| Per-period formula DSL (different EM shape per time period) | Requires a new formula syntax and a model-list parameter | Future |
| Effect modification by a **continuous** modifier interpreted as a smooth function | Requires GAM interior-smooth terms or tensor product bases | Works via GAM already in point gcomp; no API change needed |
| EM for multivariate treatment (joint A1:L, A2:L) | Needs the multivariate treatment expansion in Phase 8 | Phase 8 |
| EM for self-contained IPW with non-static interventions | Requires density ratio weights from Phase 4 (done) | — |
| EM for survival contrasts | Survival contrasts themselves are Phase 7 | Phase 7 |

## Test matrix rows added by Phase 6

| Method | Treatment | Modifier | Variance | Status (target) |
|---|---|---|---|---|
| gcomp (point) | binary | binary baseline | sandwich | ✅ truth |
| gcomp (point) | binary | continuous baseline (main effect only) | sandwich | ✅ truth |
| gcomp (point) | binary | GAM interior-smooth modifier | sandwich | ✅ truth |
| IPW | binary | binary baseline | sandwich | ✅ truth (new in this phase) |
| IPW | binary | binary baseline | bootstrap | ✅ truth (new) |
| matching | binary | binary baseline | sandwich | ✅ truth (new) |
| matching | binary | binary baseline | bootstrap | ✅ truth (new) |
| ICE | binary | binary baseline | sandwich | ✅ truth (upgrade; current cell compresses heterogeneity) |
| ICE | binary | binary baseline (3+ periods) | sandwich | ✅ truth (new) |
| ICE | binary | binary baseline | bootstrap | ✅ truth (new) |
| cross-method | binary | binary baseline | sandwich | ✅ truth (triangulation-style) |

## Items

- [ ] `R/effect_modification.R` — shared `parse_effect_mod()` helper
- [ ] `R/ipw.R` — MSM formula expansion
- [ ] `R/matching.R` — MSM formula expansion
- [ ] `R/ice.R` — `ice_build_formula()` auto-expansion of `A:modifier` across lags
- [ ] `tests/testthat/test-effect-modification.R` — new file with the rows above
- [ ] `FEATURE_COVERAGE_MATRIX.md` — upgrade the three EM cells + add the cross-method row
- [ ] Vignette updates per section 5 above
- [ ] Regression guards: every DGP in `FEATURE_COVERAGE_MATRIX.md` that does NOT use EM must give exactly the same numbers pre- and post-refactor (a non-EM fit must produce the saturated Y~A as before)
