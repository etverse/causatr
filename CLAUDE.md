# causatr

Unified causal effect estimation for methodological triangulation: g-computation
(parametric g-formula + ICE), IPW (self-contained density-ratio engine), and matching (via MatchIt).
Part of the [etverse](https://github.com/etverse) ecosystem.

## Guide files

- `FEATURE_COVERAGE_MATRIX.md` — **single source of truth for "what works".** Every PR that changes a feature MUST update this file.
- `PHASE_*.md` — per-phase implementation guides in the project root. Completed: 2–6, 8–10. In progress: 11 (diagnose). Pending: 12–20 (design docs).

## Project structure

This is an R package: `R/` (source), `tests/testthat/` (tests, `test-foo.R` mirrors `R/foo.R`), `man/` (generated — do not edit), `NAMESPACE` (generated — do not edit), `vignettes/` (long-form docs).

### R/ file layout

**Core API:** `causat.R` (main fitting), `contrast.R` (causal contrasts), `diagnose.R` (diagnostics).
**Interventions:** `interventions.R` — `static()`, `shift()`, `scale_by()`, `threshold()`, `dynamic()`, `ipsi()`.
**Estimation:** `gcomp.R`, `ice.R`, `ipw.R`, `longitudinal_ipw.R`, `matching.R`.
**Inference:** `variance_if.R` (IF sandwich engine), `variance_bootstrap.R`.
**Data:** `to_person_period.R`, `prepare_data.R`.
**S3:** `print.R`, `summary.R`, `plot.R`, `coef.R`, `confint.R`, `tidy.R`, `knit_print.R`.
**Support:** `effect_modification.R`, `ipw_weights.R`, `treatment_model.R`, `utils.R`, `checks.R`, `zzz.R`.

### S3 classes

- `causatr_fit` (from `causat()`), `causatr_result` (from `contrast()`), `causatr_diag` (from `diagnose()`), `causatr_intervention` (from intervention constructors).

## Two-step API

```r
fit <- causat(data, outcome = "Y", treatment = "A", confounders = ~ L1 + L2,
              estimator = "gcomp", model_fn = stats::glm)
result <- contrast(fit,
  interventions = list(a1 = static(1), a0 = static(0)),
  type = "difference", ci_method = "sandwich")
```

## Development commands

```r
devtools::load_all()     # load for dev
devtools::test()         # run tests
devtools::check()        # R CMD check
devtools::document()     # regenerate roxygen
```

Shell: `air format .` (format all R files).

## Code style

- `pkg::fun()` for external calls (no bare `library()`)
- `rlang::abort()` / `rlang::warn()` / `rlang::inform()` (not `stop()` / `warning()` / `message()`)
- `data.table` internally; return `data.table` from user-facing functions
- Roxygen on every function, including `@noRd` helpers
- Generous inline comments for math, design rationale, subtle invariants
- Do NOT remove existing comments unless the related code is also removed
- Exported functions at top of files, internal helpers below

## Testing rules

- `expect_snapshot(error = TRUE)` for error conditions
- NEVER delete/mock failing tests — fix the source
- Truth-based simulation tests mandatory for new features (see `FEATURE_COVERAGE_MATRIX.md`)
- External reference cross-checks against `lmtp::lmtp_tmle()` or `delicatessen` when analytical truth is hard to derive
- Update `FEATURE_COVERAGE_MATRIX.md` in the same PR as test changes

## Constraints

- Run `devtools::test()` before committing
- Do not modify `man/` or `NAMESPACE` directly
- Run `devtools::document()` after changing roxygen comments

## Scope

causatr owns g-comp (parametric g-formula + ICE), a self-contained IPW density-ratio engine, and delegates matching to MatchIt. NOT: TMLE/debiased-ML (`lmtp`), mediation, sensitivity analysis, HTE (`grf`), or forward-simulation (`gfoRmula`). Uses GLMs/GAMs (not ML) because the sandwich variance estimator requires parametric models.

## R ecosystem integration

| Need | Package | Relationship |
|---|---|---|
| Matching | `MatchIt` | **Imports** (delegated) |
| Balance diagnostics | `cobalt` | **Suggests** |
| Sandwich variance | `sandwich` | **Imports** |
| Numerical derivatives | `numDeriv` | **Imports** |
| Bootstrap | `boot` | **Imports** |
| Weight estimation | `WeightIt` | **Suggests** (test oracle only) |
| GAMs | `mgcv` | **Suggests** |

## Supported features

| Dimension | Values |
|---|---|
| **Treatment timing** | point, longitudinal (ICE + longitudinal IPW) |
| **Treatment type** | binary, continuous, categorical k>2, count (IPW: Poisson/NB), multivariate (gcomp + IPW) |
| **Outcome family** | gaussian, binomial, quasibinomial, poisson, Gamma, any GLM family, `MASS::glm.nb` |
| **Interventions** | `static`, `shift`, `scale_by`, `threshold` (gcomp only), `dynamic`, `ipsi` (IPW only) |
| **Estimand** | ATE, ATT, ATC, `by`-stratified |
| **Contrast** | difference, ratio, OR |
| **Variance** | sandwich (analytic IF), bootstrap, numeric Tier 1/2 fallback |

## Key design decisions

- **`estimator`, not `method`** — avoids shadowing `MatchIt::matchit(method = ...)`.
- **IPW MSM is `Y ~ 1` (Hájek intercept)** per intervention, not `Y ~ A`. With effect modification, expands to `Y ~ 1 + modifier`. Treatment absorbed by density-ratio weights.
- **Matching MSM is `Y ~ A`** (or `Y ~ A + modifier + A:modifier` with EM).
- **`dynamic()` = deterministic rules**, not MTPs. MTPs use `shift()` / `scale_by()` / `ipsi()`.
- **Multivariate IPW = sequential MTP** (Díaz et al. 2023); multivariate gcomp = deterministic joint transformation. They coincide for static interventions, diverge otherwise by design.
- **ICE applies intervention to current-time treatment only** — lag columns hold observed values. Recomputing lags double-counts interventions.
- **Single IF engine** — `variance_if()` in `R/variance_if.R` serves all four methods via Channel 1 (sampling) + Channel 2 (nuisance correction).
- **ICE defers model fitting to `contrast()`** — sequential outcome models are intervention-dependent.
- **`censoring =` is a row filter, not IPCW.** No censoring model is fit. Built-in IPCW is Phase 14.
- **`na.action = na.exclude` is rejected** — causes silent IF corruption via residual padding mismatch.
- **ATT/ATC only for static interventions on binary 0/1 treatment.**
- **Effect modifier must be baseline** under IPW/matching/longitudinal IPW (doc-level constraint, not enforced at runtime).
- **Stabilized weights** (`stabilize = "marginal"`) supported for multivariate IPW (Phase 8) and longitudinal IPW (Phase 10). Numerator parameters held fixed in sandwich; bootstrap refits both.
