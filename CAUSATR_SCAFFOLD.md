# CAUSATR_SCAFFOLD.md — Package Architecture for `causatr`

> Master reference for Claude Code. Defines the complete architecture of `causatr`,
> an R package for causal effect estimation via outcome-model-based g-computation.
> Part of the [etverse](https://github.com/etverse) ecosystem.
>
> Companion files (one per phase):
> - `PHASE_1_FOUNDATION.md` through `PHASE_7_ADVANCED.md` — per-phase implementation plans and status
> - `chapters/` — Individual chapter PDFs for implementation-level detail
> - Implementation guides for each phase are created from claude.ai chapter summaries when that phase begins

---

## 1. What causatr IS and IS NOT

### IS

- A **causal effect estimation package** supporting multiple methods for methodological triangulation
- **Primary engine:** outcome-model-based g-computation — the parametric g-formula (point treatments) and iterated conditional expectation (ICE) g-computation (longitudinal treatments)
- **Also supports IPW and matching** for triangulation — delegates the hard work (weight estimation, matching algorithm) to `WeightIt` / `MatchIt` but performs the causal effect estimation and inference itself
- Provides a **unified two-step API** (`causat()` → `contrast()`) across all methods
- Supports multiple **inference approaches**: sandwich (default), bootstrap
- Includes **diagnostics**: positivity checks, covariate balance summaries (via `cobalt` integration)
- Supports a comprehensive range of treatment types, outcome types, and intervention types
- Uses `data.table` internally for performance
- Part of the `etverse` ecosystem (follows `negatr` conventions)

### IS NOT

- Not reimplementing weight estimation — `WeightIt` computes weights, we consume them (for static interventions; self-contained IPW for dynamic/MTP planned for Phase 4)
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
                    method = "gcomp")

# IPW (weights from WeightIt)
fit_ipw <- causat(data, outcome = "Y", treatment = "A", confounders = ~ L1 + L2,
                  method = "ipw")

# Matching (matched data from MatchIt)
fit_match <- causat(data, outcome = "Y", treatment = "A", confounders = ~ L1 + L2,
                    method = "matching")

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
| Weight estimation | `WeightIt` | **Imports**: `causat(method = "ipw")` calls `WeightIt::weightit()` internally |
| Matching | `MatchIt` | **Imports**: `causat(method = "matching")` calls `MatchIt::matchit()` internally |
| Balance diagnostics | `cobalt` | **Suggests**: `diagnose()` calls `cobalt::bal.tab()` for balance summaries |
| Semiparametric/TMLE | `lmtp` | Out of scope. Complementary package, not a competitor. |
| Forward-simulation g-formula | `gfoRmula` | Out of scope. We supersede with ICE for longitudinal. |
| GAMs | `mgcv` | **Suggests**: used as model engine via `model_fn = mgcv::gam`. |
| Sandwich variance | `sandwich` | **Imports**: used for inference. |
| Numerical derivatives | `numDeriv` | **Imports**: Jacobian for J V_β Jᵀ propagation. |
| Bootstrap | `boot` | **Imports**: used for inference. |
| DAGs | `dagitty` | Optional (Suggests). For adjustment set identification. |
| Survival | `survival` | For `Surv` objects and baseline hazard utilities. |

---

## 4. Package Metadata

```
Package: causatr
Title: Causal Effect Estimation via G-Computation
Version: 0.0.0.9000
Imports:
    boot,
    data.table (>= 1.14.0),
    MatchIt,
    numDeriv,
    rlang (>= 1.0.0),
    sandwich,
    stats,
    WeightIt
Suggests:
    cobalt,
    mgcv,
    mice,
    quarto,
    survival,
    tinyplot,
    testthat (>= 3.0.0)
```

---

## 5. Feature Matrix

| Feature | Supported | Notes |
|---|---|---|
| **Treatment timing** | | |
| Point treatment | Yes | Standard g-computation (Ch. 13) |
| Longitudinal treatment | Planned (Phase 5) | ICE g-computation (Zivich et al.) |
| **Intervention types** | | |
| Static intervention | Yes | Set A=a for all (e.g., "always treat") |
| Dynamic intervention | Yes (gcomp only) | A_k = g(L̄_k) (e.g., "treat if CD4 < 200") |
| Modified treatment policy | Yes (gcomp only) | d(a, l) shifts observed treatment (e.g., "reduce by 10%") |
| Incremental propensity score | Planned (Phase 4) | Multiply treatment odds by δ |
| **Treatment types** | | |
| Binary treatment | Yes | |
| Continuous treatment | Yes | |
| Categorical treatment | Planned (Phase 4) | |
| Multivariate treatment | Planned (Phase 7) | Multiple simultaneous treatments |
| **Outcome types** | | |
| Continuous outcome | Yes | family = "gaussian" |
| Binary outcome | Yes | family = "binomial" |
| Censored outcome (survival) | Scaffolded (Phase 6) | Pooled logistic / discrete hazard models |
| Competing risks | Planned (Phase 6) | Cause-specific hazard models |
| **Estimation methods** | | |
| G-computation | Yes | Core engine; supports all intervention types |
| IPW (via WeightIt) | Yes | Static interventions only; self-contained IPW for dynamic/MTP in Phase 4 |
| Matching (via MatchIt) | Yes | Static interventions only |
| **Inference** | | |
| Sandwich variance | Yes | J V_β Jᵀ via sandwich + numDeriv |
| Bootstrap | Yes | Full-pipeline resampling via boot |
| **Infrastructure** | | |
| Survey weights | Planned (Phase 7) | Pass through to outcome model fitting |
| Clustered data | Planned (Phase 7) | Cluster-robust sandwich variance |
| Parallel processing | Planned (Phase 7) | `future` backend for bootstrap |

---

## 6. Inference Methods

`causatr` supports two variance estimation / CI approaches. The user selects via `ci_method` in `contrast()`.

### 6a. Sandwich variance estimator (default)

**When:** GLM-based outcome models (including pooled logistic for survival).

**How it works:** Compute the robust (Huber–White) variance of the model coefficients V_β, then propagate to the marginal means via the multivariate delta method: Var(μ̂) = J V_β Jᵀ, where J is the Jacobian computed numerically via `numDeriv::jacobian()`.

**Method-specific V_β:**
- **G-comp**: `sandwich::sandwich(model)` — standard Huber–White HC0
- **IPW**: `stats::vcov(glm_weightit_model)` — M-estimation vcov from WeightIt that accounts for weight-estimation uncertainty
- **Matching**: `sandwich::vcovCL(model, cluster = subclass)` — cluster-robust SE on matched-pair subclass

**Contrast SE:** For risk differences, the delta method gradient is trivial ([1, -1]). For risk ratios and odds ratios, the gradient involves the ratio/odds formula — applied automatically in `compute_contrast()`.

**Future:** Stacked estimating equations via `geex` package (Zivich et al. 2024) for ICE sandwich in Phase 5.

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
4. Intervention constructors (`R/interventions.R`): `static()`, `shift()`, `scale()`, `threshold()`, `dynamic()`, `ipsi()`
5. S3 class definitions + print/summary methods

### Phase 2: Point Treatment G-Computation + Inference — DONE
6. `causat(method = "gcomp")` for point treatments (`R/causat.R`, `R/gcomp.R`)
7. `contrast()` for point treatments (`R/contrast.R`)
8. Sandwich variance via J V_β Jᵀ (`R/variance_sandwich.R`)
9. Bootstrap variance (`R/variance_bootstrap.R`)
10. Support: binary + continuous outcomes, binary + continuous treatments
11. `model_fn` parameter for pluggable fitting functions

### Phase 3: IPW + Matching (Triangulation) — DONE
12. `causat(method = "ipw")` — wraps WeightIt (`R/ipw.R`) ✓
13. `causat(method = "matching")` — wraps MatchIt (`R/matching.R`) ✓
14. Vignettes: `ipw.qmd`, `matching.qmd`, `triangulation.qmd` ✓
15. Comprehensive simulation-based tests (binary/continuous outcome, all contrast types, all estimands) ✓
16. `diagnose()` — positivity, balance (cobalt), weight distribution, match quality, Love plots (`R/diagnose.R`) ✓

### Phase 4: Intervention Types + Self-Contained IPW — PENDING
16. Self-contained IPW engine (density ratio weights for dynamic/MTP interventions)
17. Modified treatment policies (shift-based interventions) via IPW
18. IPSI (incremental propensity score interventions)
19. Categorical treatment support
20. Vignette: `interventions.Rmd`

### Phase 5: Longitudinal / ICE — PENDING
21. ICE g-computation engine (`R/ice.R`)
22. `causat()` longitudinal path (detects `id` + `time`)
23. Sandwich variance for ICE (stacked estimating equations via `geex`)
24. Censoring handling (within ICE or via IPCW)
25. Vignette: `longitudinal.Rmd`

### Phase 6: Survival — SCAFFOLDED
26. `causat_survival()` — pooled logistic hazard models ✓ (basic fit)
27. `to_person_period()` — wide → long conversion ✓
28. Survival curves under intervention in `contrast()` — pending
29. Competing risks (cause-specific hazards) — pending
30. Vignette: `survival.Rmd`

### Phase 7: Advanced Features — PENDING
31. Survey weights (pass through to model fitting + adjust sandwich)
32. Clustered data (cluster-robust sandwich via `sandwich::vcovCL()`)
33. Parallel processing (`future` for bootstrap)
34. Multivariate treatment
35. Continuous treatment vignette

---

## 8. Key Design Decisions

| Decision | Rationale |
|---|---|
| ICE over forward simulation | Fewer models, sandwich variance, faster |
| GLMs/GAMs over ML | Valid sandwich inference without debiasing; users wanting ML should use `lmtp` |
| Three estimation methods | Triangulation: gcomp (outcome model), IPW (treatment model), matching (both). Agreement → confidence. |
| WeightIt/MatchIt as Imports | Don't reinvent the wheel for static interventions. Plan self-contained IPW for dynamic/MTP in Phase 4. |
| `model_fn` parameter | User passes fitting function (glm, gam, etc.) — no hardcoded model detection |
| Sandwich via J V_β Jᵀ | `sandwich` package for V_β, `numDeriv` for J. Asymptotically equivalent to stacked EE (Zivich et al. 2024). |
| ci_method: sandwich + bootstrap only | Delta method applied internally for ratio/OR contrasts; not a separate ci_method. |
| `data.table` internally | Performance for large longitudinal datasets |
| Intervention functions (not formulas) | Maximum flexibility for modified treatment policies, dynamic rules, IPSI |
| `contrast()` as a separate step | Fit once, contrast many interventions; clean separation of concerns |
| `diagnose()` integrates cobalt | Balance diagnostics are essential but cobalt already does it perfectly |
| No TMLE/DML/mediation/HTE | Totally out of scope — these are separate problems with dedicated packages |
