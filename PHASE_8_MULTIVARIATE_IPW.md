# Phase 8 — Multivariate Treatment IPW

> **Status: DONE (2026-04-20)**

## Scope

Extend the self-contained IPW density-ratio engine to support multivariate (joint) treatments `treatment = c("A1", "A2")`. Point g-comp already handles multivariate treatments; this phase closes the IPW gap.

## Current state

- **G-comp (point):** `fit_gcomp()` accepts `treatment = c("A1", "A2")` — the outcome model `Y ~ A1 + A2 + L` handles it natively. `apply_intervention()` loops over each treatment variable. Sandwich and bootstrap work.
- **IPW:** `fit_ipw()` hard-aborts at `length(treatment) > 1L` (R/ipw.R) with "Multivariate treatments are not supported under `estimator = 'ipw'`."
- **Matching:** `fit_matching()` also blocks multivariate — MatchIt itself only handles binary treatments.

## Design

### Joint treatment density

The key question: how to model $f(A_1, A_2 \mid L)$.

**Option A: Product of conditionals (sequential factorisation)**
$$
f(A_1, A_2 \mid L) = f_1(A_1 \mid L)\, f_2(A_2 \mid A_1, L)
$$
Each factor is a standard univariate density model that `fit_treatment_model()` already handles. $A_1$ enters the second model as an additional covariate. This is the simplest approach and mirrors how longitudinal IPW chains per-period models.

**Option B: Multivariate normal (continuous case)**
$$
f(A_1, A_2 \mid L) = \mathcal{N}\bigl(\mu(L), \Sigma(L)\bigr)
$$
Requires a multivariate regression model. More statistically efficient but harder to implement and less flexible (can't mix binary + continuous treatments).

**Recommendation:** Option A. It reuses the existing `fit_treatment_model()` machinery, supports mixed treatment types (e.g. one binary, one continuous), and the product factorisation is valid under the full model.

### Joint density-ratio weights (sequential MTP semantics)

The implementation follows the **sequential MTP** semantics of Díaz, Williams, Hoffman, Schenck (2023), the standard in the causal-inference MTP literature and what `lmtp` implements. For a joint intervention $d = (d_1, d_2, \ldots)$ applied stage-sequentially:

$$
w_i = \prod_{k=1}^{K} \frac{f_k\bigl(d_k^{-1}(A_{k,i}) \,\big|\, A_{1..k-1,i}^{\mathrm{obs}}, L_i\bigr) \cdot |\mathrm{Jac}\,d_k^{-1}|}{f_k\bigl(A_{k,i} \,\big|\, A_{1..k-1,i}^{\mathrm{obs}}, L_i\bigr)}.
$$

Critically, **both numerator and denominator condition on OBSERVED upstream values** — only the $k$-th argument changes ($A_k$ vs $d_k^{-1}(A_k)$). This is because the sequential MTP estimand is
$$
E[Y^d] = \int E[Y \mid A, L] \prod_k f_k\bigl(d_k^{-1}(a_k) \,\big|\, a_{1..k-1}, L\bigr)\, |\mathrm{Jac}\,d_k^{-1}|\, da\, dP(L),
$$
where the conditioning variables $a_{1..k-1}$ at the evaluation point are the observed values.

Under natural course on component $k$ (no intervention), the ratio is identically $1$ — downstream components do not accumulate an "upstream conditioning shift" from natural-course components.

**Note on mv gcomp:** multivariate gcomp implements a *different* estimand (deterministic joint transformation: simultaneous per-individual substitution of each counterfactual column). For static-only interventions the two estimators coincide; for non-static interventions on non-final components they differ by the upstream $\to$ downstream cross-dependence. See `CLAUDE.md` "Architecture notes" for the detailed semantics contract.

### Sandwich variance

The stacked M-estimation system includes parameters from both density models + the MSM. The variance engine extends `compute_ipw_if_self_contained_one()` to chain through multiple propensity models — each contributes a block to the bread matrix and a column block to the cross-derivative.

### Interventions

All intervention types that work on univariate IPW (`static`, `shift`, `scale_by`, `dynamic`, `ipsi`) should work per-component. The user specifies a list of interventions, one per treatment variable:
```r
contrast(fit,
  interventions = list(
    both_up = list(shift(1), shift(0.5)),
    control = list(static(0), static(0))
  ),
  reference = "control"
)
```

## Items (chunks 8a – 8c)

### 8a — core (commits `62a25f8` + `8a-fix` sequential-MTP semantics correction)

- [x] Remove `length(treatment) > 1L` gate in `fit_ipw()` and `check_causat_inputs()`.
- [x] Add `fit_treatment_models()` (plural) in `R/treatment_model.R` to fit the sequential factorisation $f_k(A_k \mid A_{1..k-1}, L)$ per component.
- [x] Add `compute_density_ratio_weights_mv()` in `R/ipw_weights.R` implementing the sequential-MTP per-component ratio. **Both numerator and denominator condition on observed upstream treatments** — no intervened-newdata substitution (Diaz et al. 2023). Under natural course on any component the ratio is identically $1$.
- [x] Add `make_weight_fn_mv()` building a stacked-alpha closure across $K$ models for the variance engine; per-component sub-closures via `mv_ht_closure()` (static / dynamic on discrete) / `mv_pushforward_closure()` (shift / scale on continuous or count). Natural-course components have constant-$1$ closures. `force()` every captured arg to avoid the R for-loop promise gotcha.
- [x] Add `compute_ipw_if_self_contained_mv_one()` in `R/variance_if.R`. Computes the stacked cross-derivative $[A_{\beta\alpha_1}, \ldots, A_{\beta\alpha_K}]$ via `numDeriv::jacobian` on the product-weight closure, then sums $K$ block-diagonal propensity corrections (one `apply_model_correction()` per propensity model).
- [x] Wire dispatch in `fit_ipw()`, `compute_ipw_contrast_point()`, `variance_if_ipw()` on `length(treatment) > 1L` (stored as `fit$details$is_multivariate`).
- [x] Bootstrap path (`refit_ipw()`) replays the same `treatment` slot — multivariate just works through the existing `replay_fit()` plumbing.
- [x] Reject IPSI under multivariate IPW (`causatr_multivariate_ipsi`). Matching stays binary-only.

### 8b — effect modification (commit `7f641e6`)

- [x] Lift the `causatr_multivariate_em` rejection.
- [x] Strip all treatment-touching terms from per-component propensity formulas via `parse_effect_mod()`'s `confounder_terms` slot in `fit_treatment_models()`. This covers `A_k:modifier` (EM interaction) and `A_j:A_k` (treatment-treatment interaction) uniformly.
- [x] Reuse `build_ipw_msm_formula()` to expand the per-intervention MSM from `Y ~ 1` to `Y ~ 1 + modifier_main_effects`.
- [x] Carry the Phase 6 baseline-modifier constraint (Robins 2000) as a doc-level note.

### 8c — categorical + count components (commit `fb8f04e`)

- [x] Lift the `causatr_multivariate_categorical` rejection. Lift the `propensity_family` rejection for multivariate.
- [x] Per-component fitter dispatch in `fit_treatment_models()` — categorical components get `nnet::multinom`, negbin gets `MASS::glm.nb`, rest use `model_fn`. New `propensity_model_fn` argument allows a single user-chosen fitter for all components.
- [x] `propensity_family` accepts `NULL`, length $1$ (broadcast), or length $K$ (per-component opt-in).
- [x] Categorical branch in `mv_ht_closure()` using the multinomial softmax (flattened $(K-1) \times p$ alpha reshape).
- [x] Multinomial propensity bread dispatch in `compute_ipw_if_self_contained_mv_one()` via `inherits(model, "multinom")` routing to `prepare_model_if_multinom()`.

### 8d — count components (`5f87acb`)

- [x] Add truth-based tests for count components. Bin $\times$ Poisson and bin $\times$ negbin tests use an $A_2 \perp A_1$ DGP so the truth collapses to the marginal $A_1$ effect ($1.5$). Poisson $\times$ continuous tests use the chain-rule DGP where shift($+1$) on $A_1$ propagates through $f_2(A_2 \mid A_1, L)$ via the $A_1 \to A_2$ coefficient; sequential-MTP truth is $0.3 + 0.4 \cdot 0.2 = 0.38$.
- [x] Add rejection tests for count-component interventions: `static()`, `threshold()`, `dynamic()`, non-integer `shift()`. All already dispatched to by `check_intervention_family_compat()` per component.
- Note: most of 8d's plumbing landed in 8c (per-component `propensity_family` in `fit_treatment_models()`, count branches in `mv_pushforward_closure()`). 8d is therefore test-only — no implementation changes.

### 8e — stabilized weights (`TBD`)

- [x] Add `stabilize = c("none", "marginal")` argument to `causat()`, thread through `fit_ipw()`, stash in `fit$details$stabilize`.
- [x] Extend `fit_treatment_models()` to fit per-component numerator models $g_k(A_k \mid A_{1..k-1})$ under `stabilize = "marginal"`. Stored as `attr(treatment_models, "numerator_models")`.
- [x] `compute_density_ratio_weights_mv()` swaps the numerator density to $g_k$ when stabilized; denominator stays at the full-$L$ $f_k$.
- [x] `make_weight_fn_mv()` precomputes a fixed-gamma numerator vector $f^{\mathrm{num}}_{\mathrm{fixed}}$ and routes to new `mv_stabilized_closure()` helper. The closure only varies with $\alpha$ (denominator parameters) — numerator $\gamma$ is held fixed under numDeriv perturbation (nuisance-fixed convention, matching $\sigma$ for Gaussian and $\theta$ for negbin).
- [x] `refit_ipw()` replays `stabilize` so bootstrap captures the full (including $\gamma$) uncertainty.
- [x] Reject stabilize under univariate IPW (`causatr_stabilize_univariate`).
- [x] Tests: numerator-model structure check, static+static agreement with unstabilized, shift+shift recovery of sequential-MTP truth, bootstrap refit, univariate rejection, categorical component compatibility.

### Tests

- [x] `tests/testthat/test-multivariate-ipw.R` covers: bin × bin (static, binomial outcome with diff/ratio/OR, by, subset, dynamic, bootstrap parity), bin × cont (static + shift), cont × cont (shift + shift — sequential-MTP truth), K = 3 binary, gcomp cross-check for static, gcomp ESTIMAND DIVERGENCE pin for shift + shift, cross-method $A_1{:}\mathrm{sex}$ EM, bin × cat / cat × bin / cat × cat, bin × Poisson / bin × negbin / Poisson × cont count components, plus rejection tests (IPSI, bare treatment in confounders, invalid `propensity_family` shape, count-component `static` / `threshold` / `dynamic` / non-integer shift).

## Dependencies

Phase 4 (self-contained IPW engine), Phase 6 (EM infrastructure used by 8b).

## Out of scope (deferred to other phases)

- Multivariate matching (MatchIt limitation).
- Longitudinal multivariate IPW (Phase 10).
- Sandwich variance that includes numerator-model uncertainty under `stabilize = "marginal"`. Currently the numerator parameters $\gamma$ are held fixed in the sandwich (nuisance-fixed convention); bootstrap gives the full variance. Follow-up if users need analytic full-variance SEs.
- `stabilize = "baseline"` mode (keep baseline modifiers from EM, drop other covariates from numerator). Deferred — the current `"marginal"` drops all covariates.
