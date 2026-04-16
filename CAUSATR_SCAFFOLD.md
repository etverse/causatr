# CAUSATR_SCAFFOLD.md — Package Architecture for `causatr`

> Master reference for Claude Code. Defines the complete architecture of `causatr`,
> an R package for causal effect estimation via outcome-model-based g-computation.
> Part of the [etverse](https://github.com/etverse) ecosystem.
>
> Companion files (one per phase):
> - `PHASE_1_FOUNDATION.md` through `PHASE_9_DIAGNOSE.md` — per-phase implementation plans and status
> - `chapters/` — Individual chapter PDFs for implementation-level detail
> - Implementation guides for each phase are created from claude.ai chapter summaries when that phase begins

---

## 1. What causatr IS and IS NOT

### IS

- A **causal effect estimation package** supporting multiple methods for methodological triangulation
- **Primary engine:** outcome-model-based g-computation — the parametric g-formula (point treatments) and iterated conditional expectation (ICE) g-computation (longitudinal treatments)
- **Also supports IPW and matching** for triangulation — owns the IPW pipeline end-to-end (treatment density model, density-ratio weights, per-intervention weighted MSM) and delegates only matching to `MatchIt`. `WeightIt` is invoked only by tests as a contrast-level reference for the static binary case.
- Provides a **unified two-step API** (`causat()` → `contrast()`) across all methods
- Supports multiple **inference approaches**: sandwich (default), bootstrap
- Includes **diagnostics**: positivity checks, covariate balance summaries (via `cobalt` integration)
- Supports a comprehensive range of treatment types, outcome types, and intervention types
- Uses `data.table` internally for performance
- Part of the `etverse` ecosystem (follows `negatr` conventions)

### IS NOT

- Not reimplementing the WeightIt propensity-method zoo — we DO own the IPW pipeline (treatment density model fit via the user's `propensity_model_fn`, density-ratio weights, weighted MSM), but stick to GLM/GAM-class fitters via `model_fn`. WeightIt is a test-only contrast-level reference for the static binary case.
- Not reimplementing matching algorithms — `MatchIt` does matching, we estimate effects on matched data
- Not a TMLE/debiased ML package — `lmtp` exists for this; totally out of scope
- Not a mediation package — out of scope
- Not a sensitivity analysis package — separate `etverse` package planned
- Not a heterogeneous treatment effects package — `grf` / `causal_forest` exist; out of scope
- Not a forward-simulation g-formula — we use ICE (backward iteration), not `gfoRmula`-style Monte Carlo

### Methodological triangulation

`causatr` enables comparing results from multiple estimation strategies on the same data:

```r
# G-computation (our core)
fit_gcomp <- causat(data, outcome = "Y", treatment = "A", confounders = ~ L1 + L2,
                    estimator = "gcomp")

# IPW (density-ratio weights from a self-contained treatment density model)
fit_ipw <- causat(data, outcome = "Y", treatment = "A", confounders = ~ L1 + L2,
                  estimator = "ipw")

# Matching (matched data from MatchIt)
fit_match <- causat(data, outcome = "Y", treatment = "A", confounders = ~ L1 + L2,
                    estimator = "matching")

# Same contrast() call works for all
contrast(fit_gcomp, interventions = list(a1 = static(1), a0 = static(0)))
contrast(fit_ipw,   interventions = list(a1 = static(1), a0 = static(0)))
contrast(fit_match, interventions = list(a1 = static(1), a0 = static(0)))
```

If all three give similar answers — more confidence in the result. If they diverge — investigate model specification, positivity, etc.

---

## 2. Why ICE g-computation, not forward simulation

The book (Hernán & Robins, Ch. 21) describes two g-computation approaches:

**Forward simulation (gfoRmula-style, Ch. 21.6):**
- Requires models for *every time-varying covariate* + the outcome
- Simulates covariate trajectories forward in time under intervention
- Highly susceptible to model misspecification (need p+1 correct models)
- Bootstrap-only inference (slow)

**ICE g-computation (Zivich et al., 2024; Stat in Med 43:5562–5572):**
- Requires models for the *outcome only* (at each time point)
- Iterates backward: fit E[Y_τ | Ā, L̄], predict under intervention, use predictions as pseudo-outcomes for E[Ŷ_{τ-1} | Ā_{τ-2}, L̄_{τ-2}], etc.
- Far fewer models to specify and get right
- Can be expressed as **stacked estimating equations** — empirical sandwich variance estimator (no bootstrap needed for standard inference)
- Substantially faster than bootstrap, even when bootstrap is parallelized

**causatr implements ICE g-computation.**

For point (time-fixed) treatments, ICE reduces to standard g-computation: fit E[Y | A, L], predict under A=a for all individuals, average. This is what the book calls "standardization" or the "parametric g-formula" (Ch. 13).

### On ML and the g-formula

Using black-box ML (random forests, neural nets) directly in g-computation does NOT automatically yield valid causal estimates. The theoretical problem: ML learners optimize prediction, not the bias properties needed for causal inference. Naïve plug-in g-computation with ML can be √n-inconsistent due to overfitting bias unless combined with debiasing (cross-fitting, TMLE, etc. — which are out of scope).

**causatr's approach:** Support flexible-but-parametric models — GLMs with splines, GAMs (`mgcv`), fractional polynomials — that are expressive enough to capture nonlinearity while remaining within the parametric g-formula framework where the sandwich variance estimator is valid. Users who want ML-based debiased estimation should use `lmtp`.

---

## 3. R Landscape — What Exists, What We Integrate

| Need | Package | Relationship to causatr |
|---|---|---|
| Weight estimation | `WeightIt` | Test oracle only — `causat(estimator = "ipw")` runs a self-contained density-ratio engine using `model_fn` / `propensity_model_fn` (default `stats::glm`). WeightIt is invoked only by tests as a contrast-level reference for the static binary case. |
| Matching | `MatchIt` | **Imports**: `causat(estimator = "matching")` calls `MatchIt::matchit()` internally |
| Balance diagnostics | `cobalt` | **Suggests**: `diagnose()` calls `cobalt::bal.tab()` for balance summaries |
| Semiparametric/TMLE | `lmtp` | Out of scope. Complementary package, not a competitor. |
| Forward-simulation g-formula | `gfoRmula` | Out of scope. We supersede with ICE for longitudinal. |
| GAMs | `mgcv` | **Suggests**: used as model engine via `model_fn = mgcv::gam`. |
| Sandwich variance | `sandwich` | **Imports**: used for inference. |
| Numerical derivatives | `numDeriv` | **Imports**: Jacobian for J V_β Jᵀ propagation. |
| Bootstrap | `boot` | **Imports**: used for inference. |
| Survival | `survival` | **Suggests**: for `Surv` objects and baseline hazard utilities. |

---

## 4. Package Metadata

```
Package: causatr
Title: Unified Causal Effect Estimation for Methodological Triangulation
Version: 0.0.0.9000
Imports:
    boot,
    data.table (>= 1.14.0),
    generics,
    MatchIt,
    numDeriv,
    rlang (>= 1.0.0),
    sandwich,
    stats
Suggests:
    cobalt,
    forrest,
    knitr,
    lmtp,
    MASS,
    mgcv,
    mice,
    optmatch,
    quarto,
    survival,
    tinyplot,
    WeightIt,
    testthat (>= 3.0.0)
```

---

## 5. Feature Matrix

| Feature | Supported | Notes |
|---|---|---|
| **Treatment timing** | | |
| Point treatment | Yes | Standard g-computation (Ch. 13) |
| Longitudinal treatment | Yes (Phase 5) | ICE g-computation (Zivich et al.) |
| **Intervention types** | | |
| Static intervention | Yes (gcomp / IPW / matching) | Set A=a for all (e.g., "always treat") |
| Dynamic intervention | Yes (gcomp; IPW for binary treatment via HT indicator) | Deterministic per-individual rule `A_k = d(L̄_k)` (e.g., "treat if CD4 < 200"). Continuous treatment under IPW is rejected — see the dynamic vs MTP note below. |
| Modified treatment policy (`shift`, `scale_by`) | Yes (gcomp / IPW) | Smooth diffeomorphic transformations of the observed treatment (e.g., `shift(δ)` = "add δ to every dose"). IPW handles them via the pushforward density ratio with the correct Jacobian. |
| Threshold intervention | Yes (gcomp only) | `threshold(lo, hi)` clamps continuous treatment into a range. Rejected under IPW because the pushforward of a continuous density under a clamp is a mixed measure (point masses at `lo` / `hi` plus a truncated density), so the density ratio w.r.t. Lebesgue is not well-defined. |
| Incremental propensity score (IPSI) | Yes (IPW only, binary) | Kennedy 2019 closed-form weight `w_i = (δ·A_i + (1−A_i)) / (δ·p_i + (1−p_i))`. Multiplies the treatment odds by δ. |
| **Treatment types** | | |
| Binary treatment | Yes | All methods |
| Continuous treatment | Yes (g-comp, IPW via GPS) | Matching rejects continuous with an explicit error (MatchIt is binary-only) |
| Categorical treatment (k > 2) | Yes (g-comp, IPW via multinomial PS) | Matching rejects categorical with an explicit error (MatchIt is binary-only) |
| Multivariate treatment | Yes (g-comp); IPW and matching rejected | Multivariate IPW and multivariate matching are Phase 7 work. |
| **Outcome types** | | |
| Continuous outcome | Yes | family = "gaussian" |
| Binary outcome | Yes | family = "binomial" (logit / probit / cloglog all tested) |
| Count outcome | Yes | family = "poisson" or `MASS::glm.nb` via `model_fn` |
| Fractional outcome in [0, 1] | Yes | family = "quasibinomial" |
| Positive continuous outcome | Yes | family = `Gamma(link = "log")` |
| Censored outcome (survival) | Scaffolded (Phase 6) | Pooled logistic fit done; survival curves in contrast() pending |
| Competing risks | Planned (Phase 6) | Cause-specific hazard models; parameter scaffolded in causat_survival() |
| Multinomial outcome | Planned (Phase 7) | Multi-category outcomes via `nnet::multinom()` or `VGAM`; all methods |
| Ordinal outcome | Planned (Phase 7) | Ordered categories via `MASS::polr()` or `ordinal::clm()` |
| **Estimation methods** | | |
| G-computation | Yes | Core engine; supports all intervention types |
| IPW (self-contained density-ratio engine) | Yes — `static`, `shift`, `scale_by`, `dynamic` (binary), `ipsi`. `threshold` and continuous `static` / `dynamic` rejected by design. Categorical treatment pending the multinomial sub-chunk. | One engine for every supported intervention type: treatment density model via `fit_treatment_model()` + density-ratio weight builder + per-intervention `Y ~ 1` Hájek MSM refit. WeightIt is a contrast-level test oracle, not on the runtime path. |
| Matching (via MatchIt) | Yes | Static interventions only |
| **Effect modification** | | |
| Effect modification (A × baseline modifier) | Partial | Point gcomp: yes; IPW / matching: **hard-abort** on any `A:modifier` term in `confounders` via `check_confounders_no_treatment()` (MSM hardcoded to `Y ~ A` cannot carry it); ICE: current-period A only, lags miss the interaction. Unified API planned in **Phase 8** — see `PHASE_8_INTERACTIONS.md`. |
| **Inference** | | |
| Sandwich variance | Yes | J V_β Jᵀ via sandwich + numDeriv |
| Bootstrap | Yes | Full-pipeline resampling via boot |
| **Infrastructure** | | |
| Survey weights | Planned (Phase 7) | Pass through to outcome model fitting |
| Clustered data | Planned (Phase 7) | Cluster-robust sandwich variance |
| Parallel processing | Partial | `boot::boot(parallel=, ncpus=)` done; optional `future` backend in Phase 7 |

### Dynamic regimes vs modified treatment policies — terminology gotcha

The causal-inference literature uses "dynamic regime" loosely for two distinct objects:

1. **Deterministic rules** `d(L_i)` that pick out a single counterfactual treatment value per individual (Robins; Murphy; Hernán & Robins Ch. 21). On a binary treatment the output is 0 or 1; on a continuous treatment it is a single specific real number per individual.
2. **Modified treatment policies (MTPs)** — smooth / stochastic transformations of the observed exposure such as `A → A + δ`, `A → δ·A`, or incremental propensity-score shifts (Díaz & van der Laan; Haneuse & Rotnitzky; Kennedy 2019).

`causatr`'s `dynamic()` constructor is reserved for sense (1), which is why it rejects continuous treatment under IPW — a deterministic rule on a continuous A is a Dirac per individual with no Lebesgue density, so the density ratio in the IPW weight is degenerate. MTPs in sense (2) **are** fully supported under IPW on continuous treatment via `shift()`, `scale_by()`, and `ipsi()` — they go through the smooth pushforward branch in `R/ipw_weights.R` with the correct Jacobian and serve the same scientific use cases the literature's "continuous dynamic regimes" usually mean. Users coming from an MTP background should reach for `shift()` / `scale_by()` / `ipsi()`, not `dynamic()`.

---

## 6. Inference Methods

`causatr` supports two variance estimation / CI approaches. The user selects via `ci_method` in `contrast()`.

### 6a. Sandwich variance estimator (default)

**When:** GLM-based outcome models (including pooled logistic for survival).

**How it works:** Compute the robust (Huber–White) variance of the model coefficients V_β, then propagate to the marginal means via the multivariate delta method: Var(μ̂) = J V_β Jᵀ, where J is the Jacobian computed numerically via `numDeriv::jacobian()`.

**Method-specific V_β:**
- **G-comp**: `sandwich::sandwich(model)` — standard Huber–White HC0
- **IPW**: per-individual influence function on the stacked `(α, β)` M-estimation system, built by `compute_ipw_if_self_contained_one()` — Channel 1 from the propensity score via `sandwich::estfun` + `sandwich::bread` on the propensity model, Channel 2 from the cross-derivative `A_{βα}` of the per-intervention weight closure via `numDeriv::jacobian`. The dispatcher routes every IPW fit through a single straight loop `variance_if_ipw()` that aggregates IFs via `vcov_from_if()`.
- **Matching**: `sandwich::vcovCL(model, cluster = subclass)` — cluster-robust SE on matched-pair subclass

**Contrast SE:** For risk differences, the delta method gradient is trivial ([1, -1]). For risk ratios and odds ratios, the gradient involves the ratio/odds formula — applied automatically in `compute_contrast()`.

**ICE (longitudinal):** Stacked estimating equations sandwich (Zivich et al. 2024) implemented natively via manual influence function computation — no `geex` dependency. The J V_β Jᵀ approach used for point treatments does NOT work for ICE (ignores upstream model uncertainty).

### 6b. Nonparametric bootstrap

**When:** GAMs, complex models, or as a robustness check against sandwich.

**How it works:** Resample the data B times (default 500), refit the full pipeline (model + standardization + contrast) each time, collect the distribution of estimates.

**Implementation:** Uses `boot::boot()`. Dispatches to method-specific refit functions (`refit_gcomp`, `refit_ipw`, `refit_matching`).

### API for inference

```r
# Sandwich (default) — fast, valid for GLMs
result <- contrast(fit, interventions, ci_method = "sandwich")

# Bootstrap — universal, slower
result <- contrast(fit, interventions, ci_method = "bootstrap", n_boot = 500)
```

---

## 7. Implementation Priority

### Phase 1: Foundation — DONE
1. Package scaffolding (DESCRIPTION, NAMESPACE, CLAUDE.md, Makefile, air.toml)
2. NHEFS dataset (`data/nhefs.rda`, `R/data.R`)
3. Input validation (`R/checks.R`, `R/prepare_data.R`)
4. Intervention constructors (`R/interventions.R`): `static()`, `shift()`, `scale_by()`, `threshold()`, `dynamic()`, `ipsi()`
5. S3 class definitions + print/summary methods

### Phase 2: Point Treatment G-Computation + Inference — DONE
6. `causat(estimator = "gcomp")` for point treatments (`R/causat.R`, `R/gcomp.R`)
7. `contrast()` for point treatments (`R/contrast.R`)
8. Sandwich variance via J V_β Jᵀ (`R/variance_sandwich.R`)
9. Bootstrap variance (`R/variance_bootstrap.R`)
10. Support: binary + continuous outcomes, binary + continuous treatments
11. `model_fn` parameter for pluggable fitting functions

### Phase 3: IPW + Matching (Triangulation) — DONE
12. `causat(estimator = "ipw")` — original Phase 3 implementation delegated to WeightIt (`R/ipw.R`); Phase 4 chunks 3c.i–3c.iii replaced the runtime with a self-contained density-ratio engine and chunk 3d moved `WeightIt` to `Suggests:`. ✓
13. `causat(estimator = "matching")` — wraps MatchIt (`R/matching.R`) ✓
14. Vignettes: `ipw.qmd`, `matching.qmd`, `triangulation.qmd` ✓
15. Comprehensive simulation-based tests (binary/continuous outcome, all contrast types, all estimands) ✓
16. `diagnose()` — positivity, balance (cobalt), weight distribution, match quality, Love plots (`R/diagnose.R`) ✓

### Phase 4: Intervention Types + Self-Contained IPW — RUNTIME SHIPPED
16. Self-contained IPW engine — **single** density-ratio + per-intervention weighted-MSM path covering `static` + `shift` + `scale_by` + `dynamic` (binary) + `ipsi`. WeightIt is no longer on the runtime path; it remains as a contrast-level test oracle for the static binary case. ✓ (chunks 3c.i–3c.iii)
17. `R/treatment_model.R` + `R/ipw_weights.R` — treatment density model via the user's `propensity_model_fn` (default = `model_fn`) + density-ratio weight builder with HT-indicator / smooth-pushforward / IPSI-closed-form branches. ✓ (foundation chunk; categorical multinomial branch pending in chunk 3e)
18. Per-intervention `Y ~ 1` Hájek MSM refit in `compute_contrast()` — recovers each counterfactual marginal mean from a weighted intercept; uniform across HT and pushforward weights. ✓ (chunk 3c.i)
19. IPSI (incremental propensity score interventions) via the closed-form weight, consumed by the unified `compute_ipw_contrast_point()` path. ✓ (chunk 3c.i)
20. Categorical treatment support across checks + IPW path. **Pending — chunk 3e.**
21. Single-path IPW sandwich variance: `variance_if_ipw()` straight loop calling `compute_ipw_if_self_contained_one()` (Channel 1 + Channel 2 IF on the stacked `(α, β)` system). ✓ (chunk 3b helper, chunks 3c.i–3c.ii dispatcher unification, ≈213 lines of Branch A scaffolding deleted from `R/variance_if.R`)
22. `WeightIt::glm_weightit()` contrast-level oracle tests (T-oracle1..4) for the static-binary ATE / ATT / ATC + GAM propensity cases, `WeightIt` moved from `Imports:` to `Suggests:`, and a latent correctness fix: the HT weight branch now carries a per-estimand Bayes-rule numerator `f*(L) = f(A* \| L)` so ATT / ATC produce the textbook conditional-mean weighted sums instead of silently collapsing to the ATE. Sandwich SEs agree with `glm_weightit` on the same propensity model to ~1e-6. ✓ (chunk 3d)
23. Vignettes: `interventions.qmd` (intervention-type tour), `ipw-variance-theory.qmd` (density-ratio derivation + cross-derivative). **Pending — chunks 3h / 3i.**

### Phase 5: Longitudinal / ICE — DONE
21. ICE g-computation engine (`R/ice.R`) ✓
22. `causat()` longitudinal path (detects `id` + `time`) ✓
23. Sandwich variance for ICE (stacked EE, manual influence functions) ✓
24. Censoring handling (within ICE, restrict to uncensored at each step) ✓
25. Bootstrap for ICE (resample individuals, parallel support) ✓
26. Vignette: `longitudinal.qmd` ✓

### Phase 6: Survival — SCAFFOLDED
26. `causat_survival()` — pooled logistic hazard models ✓ (basic fit)
27. `to_person_period()` — wide → long conversion ✓
28. Survival curves under intervention in `contrast()` — pending
29. Competing risks (cause-specific hazards) — pending
30. Vignette: `survival.qmd` ✓

### Phase 7: Advanced Features — PARTIAL
31. Survey weights (pass through to model fitting + adjust sandwich) — pending (basic `weights` pass-through done; explicit `survey` design integration pending)
32. Clustered data (cluster-robust sandwich via `sandwich::vcovCL()`) — pending (matching uses subclass cluster-robust today; general designs pending)
33. Parallel processing for bootstrap — `boot::boot(parallel=, ncpus=)` ✓; `future` backend pending
34. Multivariate treatment — g-comp (point + longitudinal) ✓; IPW pending (Phase 4); matching pending (Phase 7)
35. Continuous treatment vignette — pending

### Phase 8: Unified Effect-Modification API — PENDING (design doc)
36. `parse_effect_mod(confounders, treatment)` helper that detects
    treatment × modifier terms (e.g. `A:sex`) in the `confounders`
    formula — shared across every method so the convention is
    consistent.
37. IPW — `fit_ipw()` currently hardcodes a saturated `Y ~ A` MSM
    and aborts upfront on any `A:modifier` term via
    `check_confounders_no_treatment()`. Extend the MSM builder to
    include detected EM terms and drop the abort.
38. Matching — same abort path as IPW (`fit_matching()` also
    hardcodes `Y ~ A` and calls the same guard). Same fix.
39. ICE — `ice_build_formula()` resolves `A:sex` only for the
    current-period treatment slot; lagged treatments do not get
    auto-expanded modifier interactions, compressing heterogeneity
    in multi-period DGPs. Auto-expand `A:modifier` to include the
    same interaction with every currently-available lag per period.
40. Cross-method truth test: gcomp / IPW / matching / ICE all run
    on the same EM DGP with the same formula and agree within the
    usual triangulation tolerance.
41. Regression guards: every existing DGP without EM terms must
    give identical numbers pre- and post-refactor — a non-EM fit
    must still produce the saturated `Y ~ A` MSM it does today.
42. Vignette updates across `gcomp.qmd`, `ipw.qmd`, `matching.qmd`,
    `longitudinal.qmd`, and `triangulation.qmd` documenting the
    convention. Full plan: `PHASE_8_INTERACTIONS.md`.

### Phase 9: Full `diagnose()` rewrite for Phase 4 architecture — PENDING (design doc)
43. Intervention-aware `diagnose()` with `intervention =` argument,
    mirroring `contrast()`. One diagnostic panel per intervention
    the user passes in.
44. Longitudinal (ICE) dispatch path with per-time-step positivity,
    balance, and censoring diagnostics.
45. Treatment-type-aware dispatch (binary / categorical / continuous),
    routed through `detect_treatment_family()` rather than WeightIt's
    `treat.type` attribute.
46. Estimand-aware balance: SMDs under ATE vs ATT vs ATC use
    different reference populations, and the rewrite makes that
    explicit.
47. Effect-modification-aware stratification (Phase 8 dependency):
    optional `by = "modifier"` argument to report balance within
    modifier strata.
48. Redesigned `causatr_diag` object (nested per-intervention) and
    corresponding print / summary / plot methods.
49. `vignettes/diagnostics.qmd` — user-facing tour.
50. Phase 4 ships a **minimal shim** so `diagnose()` keeps working
    on the common binary static ATE cross-sectional case; the full
    rewrite is scoped in `PHASE_9_DIAGNOSE.md` and depends on
    Phases 4, 5, and 8 being in place.

---

## 8. Key Design Decisions

| Decision | Rationale |
|---|---|
| ICE over forward simulation | Fewer models, sandwich variance, faster |
| GLMs/GAMs over ML | Valid sandwich inference without debiasing; users wanting ML should use `lmtp` |
| Three estimation methods | Triangulation: gcomp (outcome model), IPW (treatment model), matching (both). Agreement → confidence. |
| MatchIt delegated; WeightIt as test oracle | Matching is genuinely hard and we delegate. IPW we own end-to-end via a single density-ratio engine that handles `static` + `shift` + `scale_by` + `dynamic` (binary) + `ipsi` uniformly — reusing `model_fn` / `propensity_model_fn` for the treatment density model. `threshold` and continuous `dynamic` / `static` are rejected by design (mixed measure / Dirac pushforward). WeightIt is an independent reference for the static binary case in tests and is not on the runtime path. |
| `model_fn` parameter | User passes fitting function (glm, gam, etc.) — no hardcoded model detection |
| Sandwich via J V_β Jᵀ | `sandwich` package for V_β, `numDeriv` for J. Asymptotically equivalent to stacked EE (Zivich et al. 2024). |
| ci_method: sandwich + bootstrap only | Delta method applied internally for ratio/OR contrasts; not a separate ci_method. |
| `data.table` internally | Performance for large longitudinal datasets |
| Intervention functions (not formulas) | Maximum flexibility for modified treatment policies, dynamic rules, IPSI |
| `contrast()` as a separate step | Fit once, contrast many interventions; clean separation of concerns |
| `diagnose()` integrates cobalt | Balance diagnostics are essential but cobalt already does it perfectly. Phase 4 ships a minimal shim for the new IPW architecture (binary static ATE works); full intervention-/treatment-/estimand-aware rewrite is `PHASE_9_DIAGNOSE.md`. |
| No TMLE/DML/mediation/HTE | Totally out of scope — these are separate problems with dedicated packages |
