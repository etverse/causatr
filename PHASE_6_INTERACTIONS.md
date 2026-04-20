# Phase 6 — Unified Effect-Modification API

> **Status: DONE**
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
| **Longitudinal ICE** | ⚠️ `ice_build_formula()` injects `baseline_terms` verbatim at every period. `A:sex` resolves to `A_k:sex` — the **current-period** treatment × modifier. Lagged treatments (`lag1_A`, `lag2_A`, ...) do NOT get a corresponding `lag1_A:sex` term. | Partial. Current-period effect modification is captured; earlier-period contributions are collapsed. MC evidence: 2-period DGP with persistent $(1 + 1.5 \cdot \text{sex}) \cdot A_k$ effect returns contrasts 2.81 / 4.43 vs structural truth 2.10 / 5.10 — roughly half the intended heterogeneity. |

Phase 6 will replace the hard-abort path with a proper MSM builder that honors `A:modifier` terms for IPW and matching, and will auto-expand `lag{k}_A:modifier` for ICE.

The root cause in each case is different, which is part of why this needs a dedicated phase rather than four scattered fixes.

## Design goals

1. **One formula convention.** Users write `~ L + sex + A:sex` and all four methods do the right thing. The convention should be: a term involving the treatment (e.g. `A:sex`, `A:L`, `A:I(age > 65)`) is interpreted as effect modification of A on the outcome scale.
2. **No silent failures.** If a user writes an EM term and the chosen method genuinely can't support it (e.g. IPW + multivariate treatment + continuous modifier in Phase 4), abort with a specific error that names the method, the term, and the workaround.
3. **Leave (B) and (C) alone.** Interactions *among confounders* (e.g. `L1:L2`) and *within the current period* (e.g. `A:L` where L is TV) already work across methods. This phase only touches the `A` $\times$ `baseline` (and its time-varying analogs) case.
4. **No breaking changes to existing tests.** Every truth-test currently in the matrix must still pass after the refactor. Correct behavior under additive models collapses to the current behavior by construction.

## Design correction: IPW MSM is `Y ~ 1`, not `Y ~ A`

The original draft described the IPW MSM as `Y ~ A`. In practice, the self-contained
IPW engine refits an **intercept-only** `Y ~ 1` per intervention (Hájek mean),
because the density-ratio weights collapse every nonzero row onto the same
counterfactual treatment value. A saturated `Y ~ A` would be rank-deficient under
HT weights on a binary treatment.

This means the EM expansion for IPW is:

- **Without EM:** `Y ~ 1` per intervention (unchanged).
- **With EM:** `Y ~ 1 + modifier_main_effects` per intervention (e.g. `Y ~ 1 + sex`).
  The modifier main effect interacts with the intercept, so `predict()` on the
  target population returns stratum-specific counterfactual means. No treatment
  term goes into the MSM — the weights absorb it.

For matching, the MSM is already `Y ~ A`, so the expansion is `Y ~ A + sex + A:sex`.

## Known limitation: modifier must be **baseline** under IPW / matching MSM

**Under IPW and matching, the modifier in `A:modifier` must be baseline — not affected by earlier treatment.** Under ICE this extends to: the modifier must be measurable **before** the treatment whose effect it modifies (so at stage $k$, any component of $\bar{L}_{k-1}$ is fine, but $L_k$ itself is not, because it can be affected by $A_{k-1}$).

### Why

MSMs parameterise the marginal counterfactual mean:

$$
E[Y^{\bar{a}} \mid V] = g^{-1}\!\big(\beta_0 + \beta \cdot \text{cum\_a} + \gamma \cdot V\big),
$$

where $V$ is **baseline** covariates. The density-ratio weights explicitly integrate out time-varying covariates $\bar{L}$ by design — that is what "marginal" means in *marginal* structural model. Putting a time-varying $L_k$ into the MSM as an effect modifier is equivalent to **conditioning on a post-treatment variable**: $L_k$ can be affected by earlier $A_j$ ($j < k$), so the resulting stratum-specific estimand mixes (i) genuine effect modification with (ii) selection bias through mediators and (iii) collider bias on treatment–outcome paths (Robins 2000; Hernán & Robins Ch. 12.6). The estimate has no clean causal interpretation.

This restriction is a property of MSMs, not of causatr's implementation — the same warning applies to every MSM-based estimator in any package. It is not enforced by the current `check_em_compat()` gate; a Phase 6 follow-up should add an opt-in `check_em_baseline_only()` check that fires on `estimator = "ipw"` / `"matching"` when the modifier is not demonstrably baseline (hard to do automatically — time-varying status is not encoded in the data structure — so the check will likely rely on a user-supplied `baseline_cols = ...` contract).

### What each estimator actually supports

| Estimator | Baseline modifier | Time-varying modifier |
|---|---|---|
| **gcomp (point)** | ✓ (outcome model standard) | ✓ — the outcome model conditions on the full $(A, L)$ and standardizes over $L$; no post-treatment selection |
| **IPW (point)** | ✓ | ⛔ — silent bias; MSM conditions on post-treatment variable |
| **matching (point)** | ✓ | ⛔ — same issue as IPW |
| **ICE (longitudinal gcomp)** | ✓ | ⚠️ with care — the current-period modifier $L_k$ is post-treatment relative to $A_{k-1}$, so $A_{k-1}:L_k$ has the same interpretive problem. Modifiers should be from $\bar{L}_{k-1}$ (pre-period history). |
| **IPW (longitudinal, Phase 10)** | ✓ | ⛔ — same marginal-vs-conditional issue, repeated per period |

### Scientific workaround: use an SNM

The correct tool for **genuine time-varying effect modification** is g-estimation of a structural nested mean model (SNM). SNMs parameterise the per-stage blip $\gamma(a_k, \bar{l}_k, \bar{a}_{k-1}; \psi)$ directly — so a stage-$k$ modifier drawn from $\bar{L}_k$ is supported by design, not by accident. causatr's SNM pathway is scoped in `PHASE_18_SNMS.md`. Users who write `A:L_k` with time-varying $L_k$ should be pointed to `estimator = "snm"` (once Phase 18 ships).

### Runtime behaviour until Phase 18

Pre-Phase-18, causatr has no way to detect time-varying modifiers from data alone, so it cannot raise a targeted error. The doc-level mitigation is:

- This section in `PHASE_6_INTERACTIONS.md`.
- A parallel warning in the IPW / matching / longitudinal vignettes.
- A sentence in the `?causat` roxygen `@details` naming the restriction.
- A single-line note in CLAUDE.md's Phase 6 description.

A runtime `check_em_baseline_only()` via an explicit `baseline_cols = character()` argument is a natural follow-up once users start tripping over the restriction in practice. Flagged as an open item, not as Phase 6 scope creep.

## Plan — chunk sequence

### Chunk 6a — `parse_effect_mod()` helper + gate refactoring

**Scope:** parser infrastructure only; no estimation changes.

In a new `R/effect_modification.R`, add:

```r
parse_effect_mod <- function(confounders, treatment) {
  # Walk the term.labels of `confounders`. For each term that
  # involves any variable in `treatment`, return a list with
  #   $term: the full term label (e.g. "A:sex")
  #   $treatment_var: the matching treatment variable
  #   $modifier_vars: the other variables in the interaction
  # Terms without any treatment variable are left alone.
  # Returns a list of class "causatr_em_info":
  #   $em_terms — list of parsed EM terms (each with $term, $treatment_var, $modifier_vars)
  #   $confounder_terms — character vector of non-EM term labels
  #   $modifier_vars — unique modifier variable names across all EM terms
  #   $has_em — logical scalar
}
```

This is the canonical detector every method will consult. Anchoring it here (rather
than ad-hoc regex in each fitter) keeps the convention consistent and makes future
extensions (e.g. three-way interactions) mechanical.

Refactor `check_confounders_no_treatment()` in `R/utils.R`: instead of aborting
unconditionally for IPW/matching, call `parse_effect_mod()` and only abort for
*genuinely unsupported* EM patterns (EM + multivariate treatment under IPW,
EM + non-binary treatment under matching, EM + non-static intervention under IPW).
The function becomes a targeted guard rather than a blanket ban.

**Tests:**

- Unit tests for `parse_effect_mod()` on various formula shapes: `A:sex`, `sex:A`,
  `A:I(age>65)`, `A:sex + A:race`, `L1:L2` (no treatment), multivariate `c("A1","A2")`.
- The existing hard-abort tests in `test-ipw.R` and `test-matching.R` must be updated
  to reflect the new (narrower) rejection conditions.
- Regression: all existing tests must still pass — no estimation paths change.

**Files touched:** `R/effect_modification.R` (new), `R/utils.R`, `tests/testthat/test-effect-modification.R` (new).

---

### Chunk 6b — IPW MSM expansion for effect modification

**Scope:** IPW estimation path only.

In `compute_ipw_contrast_point()`, when `parse_effect_mod()` detects EM terms,
expand the per-intervention MSM from `Y ~ 1` to `Y ~ 1 + modifier_main_effects`
(e.g. `Y ~ 1 + sex`). Predict on the full target population (modifier-aware), then
average per stratum under `by` or overall.

Wire `parse_effect_mod()` into `fit_ipw()` so the detected modifiers are stored in
`fit$details$em_info` for downstream use by `compute_ipw_contrast_point()` and
`variance_if_ipw()`.

**Variance:** `variance_if_ipw()` needs the expanded MSM's score/bread. Since it
already works with arbitrary `glm` MSMs via `sandwich::estfun` / `sandwich::bread`,
this should flow through — verify.

**Tests:**

- Truth-based: binary treatment × binary modifier × IPW sandwich — DGP with known
  stratum-specific ATE, cross-checked against gcomp.
- Truth-based: IPW bootstrap on same DGP.
- Regression guard: IPW without EM terms produces identical results.

**Files touched:** `R/ipw.R`, `R/variance_if.R` (if needed), `tests/testthat/test-effect-modification.R`.

---

### Chunk 6c — Matching MSM expansion for effect modification

**Scope:** matching estimation path only.

In `fit_matching()`, when `parse_effect_mod()` detects EM terms, expand the outcome
MSM from `Y ~ A` to `Y ~ A + modifier + A:modifier`.

Store the EM metadata in `fit$details$em_info`.

**Variance:** matching already uses `cluster = subclass` with `prepare_model_if()` on
the weighted GLM — the expanded formula should flow through.

**Tests:**

- Truth-based: binary treatment × binary modifier × matching sandwich — same DGP as 6b.
- Truth-based: matching bootstrap on same DGP.
- Regression guard: matching without EM terms gives identical results.

**Files touched:** `R/matching.R`, `tests/testthat/test-effect-modification.R`.

---

### Chunk 6d — ICE lag auto-expansion for `A:modifier` terms

**Scope:** ICE formula builder only.

Extend `ice_build_formula()` in `R/ice.R` so that when a baseline term of the form
`treatment:modifier` appears, it auto-expands to include the same interaction with
every currently-available lag:

```r
# At time_idx = 2, max_lag = 2, treatment = "A", term = "A:sex":
#   emit "A:sex", "lag1_A:sex", "lag2_A:sex"
```

The expansion is per-period (later periods have more lags) and defaults to "this
interaction applies uniformly across time". This handles the **time-invariant effect
modifier** semantics. The rarer **time-varying effect modifier** semantics (different
functional form per period) remains out of scope until a per-period formula DSL
lands — document the limitation and point users at wide-format + point gcomp in the
meantime.

**Tests:**

- Truth-based: ICE × 2-period DGP × $(1 + \gamma \cdot \text{sex}) \cdot A_k$ × `by(sex)` — must recover
  stratum-specific contrast to ~5% of MC truth (vs current ~30% compression).
- Truth-based: ICE × 3-period DGP for deeper lag coverage.
- Regression guard: ICE without EM terms gives identical numbers pre/post refactor.
- ICE bootstrap on the EM DGP.
- ICE × `A:sex + A:age` — multiple EM terms; auto-expansion handles both.

**Files touched:** `R/ice.R`, `tests/testthat/test-effect-modification.R`.

---

### Chunk 6e — Cross-method triangulation test + docs + matrix update

**Scope:** integration testing and documentation.

Once chunks 6b–6d are done, the `by` branch in `compute_contrast()` just works — it
already averages predictions per stratum, and the MSM fix ensures predictions depend
on the stratum. No changes to `by` itself.

**Tests:**

- Cross-method truth test: gcomp, IPW, matching, ICE all run on the same EM DGP with
  the same formula; stratum-specific contrasts agree within cross-method tolerance.

**Documentation:**

- Update `FEATURE_COVERAGE_MATRIX.md`: upgrade the three EM cells + add cross-method row.
- Update `NEWS.md` to reflect the fix.
- Update `CLAUDE.md` architecture notes (Phase 6 status, EM design notes).
- Mark all items in this Phase 6 doc as done.
- New vignette section in `vignettes/gcomp.qmd`, `vignettes/ipw.qmd`,
  `vignettes/matching.qmd`, and `vignettes/longitudinal.qmd` showing how to specify
  effect modification for each method.
- Expand `vignettes/triangulation.qmd` with an EM example.

**Files touched:** `tests/testthat/test-effect-modification.R`, `FEATURE_COVERAGE_MATRIX.md`,
`NEWS.md`, `CLAUDE.md`, `PHASE_6_INTERACTIONS.md`, vignettes.

## Out of scope for Phase 6

| Topic | Reason | Deferred to |
|---|---|---|
| Per-period formula DSL (different EM shape per time period) | Requires a new formula syntax and a model-list parameter | Future |
| Effect modification by a **continuous** modifier interpreted as a smooth function | Requires GAM interior-smooth terms or tensor product bases | Works via GAM already in point gcomp; no API change needed |
| EM for multivariate treatment (joint A1:L, A2:L) | Needs the multivariate treatment expansion in Phase 8 | Phase 8 |
| EM for self-contained IPW with non-static interventions | Requires density ratio weights from Phase 4 (done) | — |
| EM for survival contrasts | Survival analysis is now a separate etverse package | Separate package (see `SURVIVAL_PACKAGE_HANDOFF.md`) |
| **EM by a time-varying modifier under IPW / matching** | MSM conditioning on post-treatment variables is a **structural limitation of MSM identification**, not an implementation gap. See § "Known limitation: modifier must be **baseline**". Correct tool is a structural nested model. | Phase 18 (`PHASE_18_SNMS.md`) |
| **Runtime `check_em_baseline_only()` guard** | Time-varying status is not inferable from the data structure; user must declare via an explicit `baseline_cols =` contract | Follow-up once users trip over the silent-bias risk in practice |

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

### Chunk 6a — parser + gate refactoring
- [x] `R/effect_modification.R` — `parse_effect_mod()` + `build_ipw_msm_formula()` + `build_matching_msm_formula()` + `check_em_compat()` + `em_confounder_terms()`
- [x] `R/utils.R` — replaced `check_confounders_no_treatment()` with `check_confounders_treatment()` + `build_ps_formula()` now strips EM terms
- [x] `tests/testthat/test-effect-modification.R` — 65 parser unit tests + updated rejection tests
- [x] All existing tests pass (no estimation changes)

### Chunk 6b — IPW MSM expansion
- [x] `R/ipw.R` — `compute_ipw_contrast_point()` MSM expansion via `build_ipw_msm_formula()` + `fit$details$em_info`
- [x] `R/variance_if.R` — verified: expanded MSM flows through `variance_if_ipw()` unchanged (J, X_star, phi_bar all generalize from p_beta=1 to p_beta>1)
- [x] `tests/testthat/test-effect-modification.R` — IPW truth (DGP 4, sandwich), bootstrap, gcomp cross-check, regression guard
- [x] `tests/testthat/test-ipw.R` — updated rejection test $\to$ acceptance test

### Chunk 6c — Matching MSM expansion
- [x] `R/matching.R` — `fit_matching()` MSM expansion via `build_matching_msm_formula()` + `fit$details$em_info` + formula env fix
- [x] `R/variance_bootstrap.R` — `refit_matching()` replays EM-expanded MSM + formula env fix
- [x] `tests/testthat/test-effect-modification.R` — matching truth (DGP 4, sandwich), bootstrap, gcomp cross-check, regression guard
- [x] `tests/testthat/test-matching.R` — updated rejection test → acceptance test

### Chunk 6d — ICE lag auto-expansion
- [x] `R/effect_modification.R` — `expand_em_lag_terms()` helper: splits interaction on `:`, substitutes treatment with `lag{k}_{trt}`, returns expanded terms
- [x] `R/ice.R` — `fit_ice()` stores `em_info` from `parse_effect_mod()`; `ice_build_formula()` accepts `em_info` and appends lag-expanded EM terms with all-NA column validity check
- [x] `tests/testthat/helper-dgp.R` — `make_em_ice_scm()` DGP with sex-specific treatment effects
- [x] `tests/testthat/test-effect-modification.R` — ICE truth (2-period ATE|sex=0=5/sex=1=8, 3-period ATE|sex=0=8/sex=1=12.5) + multiple EM terms smoke + bootstrap + regression guard + `expand_em_lag_terms()` unit tests

### Chunk 6e — Cross-method triangulation + docs + matrix
- [x] `tests/testthat/test-effect-modification.R` — cross-method triangulation test (gcomp + IPW + matching on DGP 4, truth + pairwise agreement)
- [x] `FEATURE_COVERAGE_MATRIX.md` — added cross-method EM triangulation row, marked Phase 6 complete
- [x] `NEWS.md` — Phase 6 completion entry
- [x] `CLAUDE.md` — updated Phase 6 status to done
- [x] `PHASE_6_INTERACTIONS.md` — all items marked done
- [x] Vignette updates: ipw.qmd (worked example), matching.qmd (worked example), longitudinal.qmd (ICE EM section), triangulation.qmd (cross-method EM example), gcomp.qmd (updated callout), introduction.qmd (updated EM description)
