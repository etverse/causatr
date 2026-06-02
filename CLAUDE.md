# causatr

Unified causal effect estimation for methodological triangulation: g-computation
(parametric g-formula + ICE), IPW (self-contained density-ratio engine), AIPW (doubly-robust), and matching (via MatchIt).
Part of the [etverse](https://github.com/etverse) ecosystem.

## Guide files

- `FEATURE_COVERAGE_MATRIX.md` — **single source of truth for "what works".** Every PR that changes a feature MUST update this file.
- `PHASE_*.md` — per-phase implementation guides in the project root. Completed: 2–6, 8–21 (21a–21e + the longitudinal/IPCW combinations), 22a (stratified ICE), 25. Pending: 22b (natural-history MTPs — designed in PHASE_22, not implemented), 23–24, 26 (design docs); Phase 21 vignette (21i) still to write.

## Project structure

This is an R package: `R/` (source), `tests/testthat/` (tests, `test-foo.R` mirrors `R/foo.R`), `man/` (generated — do not edit), `NAMESPACE` (generated — do not edit), `vignettes/` (long-form docs).

### R/ file layout

**Core API:** `causat.R` (main fitting), `contrast.R` (causal contrasts), `diagnose.R` (main dispatch + panel helpers) + `diagnose_longitudinal.R` + `diagnose_positivity.R` + `diagnose_balance.R` + `diagnose_weights.R` + `diagnose_censoring.R` + `diagnose_intervention.R` + `diagnose_sampling.R` (sampling-model diagnostics for transport).
**Interventions:** `interventions.R` — `static()`, `shift()`, `scale_by()`, `threshold()`, `dynamic()`, `ipsi()`, `stochastic()`.
**Estimation:** `gcomp.R`, `ice.R` (+ `ice_stratified.R` — per-stratum fit/predict helpers for `stratified = "G"`), `ipw.R`, `aipw.R` (point AIPW), `aipw_longitudinal.R` (longitudinal AIPW), `longitudinal_ipw.R` (longitudinal IPW fitter) + `longitudinal_ipw_formula.R` (per-period propensity/numerator formula builders) + `longitudinal_ipw_contrast.R` (per-intervention contrast bundles), `matching.R`, `snm.R` (SNM fitter + blip spec), `snm_point.R` (point g-estimation + blip design matrix), `snm_contrast.R` (SNM contrast dispatch + averaged blip), `snm_longitudinal.R` (longitudinal backward-sequential g-estimation), `censoring.R` (IPCW model fitting + weights).
**Inference (IF sandwich):** `variance_if_core.R` (model-correction primitives + vcov aggregation), `variance_if.R` (main dispatcher + numeric fallback + point channel + gcomp/matching IF), `variance_if_ice.R`, `variance_if_ipw.R` (point + mv IPW IF) + `variance_if_ipw_longitudinal.R` (longitudinal IPW IF) + `variance_if_ipw_sampling.R` (transport sampling-model corrections), `variance_if_aipw.R` (point AIPW IF) + `variance_if_aipw_longitudinal.R` (longitudinal AIPW IF), `variance_if_snm.R` (point SNM sandwich), `variance_if_snm_longitudinal.R` (longitudinal SNM cluster sandwich).
**Inference (bootstrap):** `variance_bootstrap.R` (core + refitters), `variance_bootstrap_longitudinal.R` (longitudinal IPW/ICE/AIPW bootstrap), `variance_bootstrap_snm.R` (point + longitudinal SNM bootstrap).
**Multiple imputation:** `causat_mice.R` (exported `causat_mice()` + per-imputation loop) + `pool_rubin.R` (Rubin's rules + Barnard-Rubin df) + `pool_boot_mi.R` (von Hippel two-stage bootstrap-then-impute).
**Data:** `to_person_period.R`, `prepare_data.R`, `data.R` (dataset documentation).
**S3:** `print.R`, `summary.R`, `plot.R`, `coef.R`, `confint.R`, `tidy.R`, `knit_print.R`.
**Support:** `effect_modification.R`, `ipw_weights.R` (point weights), `ipw_weights_mv.R` (mv weights) + `ipw_weights_mv_closure.R` (mv HT/pushforward/stabilized alpha-closures), `ipw_weights_shared.R` (cross-cutting helpers: `apply_intervention_to_values`, `ipsi_weight_formula`, `ht_bayes_numerator`, `check_intervention_family_compat`), `ipw_weights_longitudinal.R` (longitudinal weights) + `ipw_weights_longitudinal_closure.R` (stacked-alpha replay closures), `treatment_model.R`, `sampling_model.R` (sampling model for transport), `constructors.R` (`new_causatr_*` S3 constructors), `family.R` (GLM family helpers), `utils.R` (misc helpers), `checks.R`, `target_trial.R` (target trial protocol), `causat_mice.R` (multiple-imputation pooling — see Multiple imputation above), `causatr-package.R` (package-level doc), `zzz.R`.

### S3 classes

- `causatr_fit` (from `causat()`), `causatr_result` (from `contrast()`), `causatr_diag` (from `diagnose()`), `causatr_intervention` (from intervention constructors).

## Two-step API

```r
# Simple: same confounders for all models (backward-compatible)
fit <- causat(data, outcome = "Y", treatment = "A", confounders = ~ L1 + L2,
              estimator = "gcomp", model_fn = stats::glm)

# Explicit: separate confounders for outcome and treatment models
fit <- causat(data, outcome = "Y", treatment = "A",
              confounders_outcome = ~ L1 + L2 + I(L1^2),
              confounders_treatment = ~ L1 + L2,
              estimator = "aipw", propensity_model_fn = stats::glm)

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

## Cost discipline

- **Targeted tests**: Use `devtools::test(filter = "foo")` during development. Only run the full `devtools::test()` suite before committing. A hook blocks unfiltered test runs.
- **Foreground tests**: Run all test and check commands in **foreground** with `timeout: 600000` (10 min). Never use `run_in_background` for `devtools::test()`, `testthat::test_file()`, `devtools::check()`, or `R CMD check`. A hook enforces this. Output comes back directly — no polling needed.
- **Batch R scripts**: Combine multiple diagnostic/validation checks into a single `Rscript -e '...'` call instead of running them one at a time. Each R process startup costs 10-30 seconds of idle context.
- **Model awareness**: For routine work (formatting, simple edits, running tests, git operations), Sonnet is 5x cheaper than Opus. Suggest `/model sonnet` to the user when entering a routine-work phase, and `/model opus` when hard reasoning is needed (debugging subtle bugs, designing new features, variance derivations).

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
| Transport cross-check | `TransportHealth` | **Suggests** (test oracle only) |
| SNM cross-check | `DTRreg` | **Suggests** (test oracle only) |

## Supported features

| Dimension | Values |
|---|---|
| **Estimator** | gcomp, ipw, aipw, matching, snm |
| **Treatment timing** | point, longitudinal (ICE + longitudinal IPW + longitudinal AIPW + longitudinal SNM), transportability (`target =`) |
| **Treatment type** | binary, continuous, categorical k>2, count (IPW: Poisson/NB), multivariate (gcomp + IPW + AIPW + longitudinal IPW + longitudinal AIPW) |
| **Outcome family** | gaussian, binomial, quasibinomial, poisson, Gamma, any GLM family, `MASS::glm.nb`, `betareg::betareg` (beta regression) |
| **Interventions** | `static`, `shift`, `scale_by`, `threshold` (gcomp only), `dynamic`, `ipsi` (IPW only — point + univariate longitudinal), `stochastic` (gcomp only; IPW/AIPW when `density` supplied — Phase 24) |
| **Estimand** | ATE, ATT, ATC, `by`-stratified |
| **Contrast** | difference, ratio, OR |
| **Variance** | sandwich (analytic IF), bootstrap, numeric Tier 1/2 fallback |
| **Missing data** | complete-case (default), IPCW (missing Y, `ipcw = TRUE`), multiple imputation (missing L/A, `causat_mice()`) |

## Key design decisions

- **`estimator`, not `method`** — avoids shadowing `MatchIt::matchit(method = ...)`.
- **IPW MSM is `Y ~ 1` (Hájek intercept)** per intervention, not `Y ~ A`. With effect modification, expands to `Y ~ 1 + modifier`. Treatment absorbed by density-ratio weights.
- **Matching MSM is `Y ~ A`** (or `Y ~ A + modifier + A:modifier` with EM).
- **`dynamic()` = deterministic rules**, not MTPs. MTPs use `shift()` / `scale_by()` / `ipsi()`.
- **Longitudinal IPSI (Phase 20)** — univariate `ipsi(delta)` under longitudinal IPW extends Kennedy's (2019) closed form to a per-period product $W_i = \prod_t (\delta A_t + (1 - A_t)) / (\delta p_t + (1 - p_t))$. No new machinery: each period reuses `compute_density_ratio_weights()`'s IPSI branch and the stacked sandwich reuses `make_weight_fn()`'s IPSI sub-closure. **Rejected**: multivariate IPSI (`causatr_longitudinal_ipsi_multivariate` — closed form is binary-univariate) and IPSI + `stabilize = "marginal"` (`causatr_longitudinal_ipsi_stabilize` — Kennedy's weight has no separate marginal numerator). Longitudinal **AIPW** IPSI stays rejected (no counterfactual treatment stream for the ICE recursion).
- **Multivariate IPW = sequential MTP** (Díaz et al. 2023); multivariate gcomp = deterministic joint transformation. They coincide for static interventions, diverge otherwise by design.
- **ICE applies intervention to current-time treatment only** — lag columns hold observed values. Recomputing lags double-counts interventions.
- **Single IF engine** — `variance_if()` in `R/variance_if.R` serves all five estimators via Channel 1 (sampling) + Channel 2 (nuisance correction).
- **ICE defers model fitting to `contrast()`** — sequential outcome models are intervention-dependent.
- **`censoring =` is a row filter by default.** With `ipcw = TRUE`, an internal censoring model is fit and IPCW weights are computed (Phase 14, shipped).
- **`na.action = na.exclude` is rejected** — causes silent IF corruption via residual padding mismatch.
- **ATT/ATC only for static interventions on binary 0/1 treatment.**
- **Effect modifier must be baseline** under IPW/matching/longitudinal IPW (doc-level constraint, not enforced at runtime). Baseline EM composes with multivariate longitudinal IPW (Phase 19c): each per-period × per-component propensity strips its `A:modifier` interactions (keeping the modifier as a confounder main effect) and the single final-period MSM expands to `Y ~ 1 + modifier`. **SNM lifts the baseline restriction** — time-varying effect modification is the headline use case for `estimator = "snm"` (Phase 18).
- **SNM estimates blip parameters, not counterfactual means** — `contrast()` with `interventions =` is rejected for SNM fits. The blip ψ is the estimand directly; no intervention specification is needed.
- **Stabilized weights** (`stabilize = "marginal"`) supported for multivariate IPW (Phase 8), multivariate AIPW (Phase 16m), longitudinal IPW (Phase 10), and multivariate longitudinal IPW (Phase 19b — per-period × per-component chain numerator that drops time-varying L). Numerator parameters held fixed in sandwich; bootstrap refits both. **Rejected for longitudinal AIPW** (univariate + multivariate, Phase 25): `fit_aipw()` never threaded `stabilize` into `fit_aipw_longitudinal()`, and the ICE-AIPW recursion uses Hájek-normalized single-period weights — a classed error (`causatr_stabilize_longitudinal_aipw`) points users to `estimator = "ipw"` for stabilized longitudinal weights.
- **Multivariate longitudinal AIPW** (Phase 25) composes the ICE-AIPW backward recursion (Bang & Robins 2005; outcome side already loops over vector `treatment`) with Phase 19's per-period × per-component density chain (`fit_treatment_models()`, `compute_density_ratio_weights_mv()`). The longitudinal AIPW sandwich (Channel 2b in `variance_if_aipw_long_one()`) stacks $T \times K$ block-diagonal propensity corrections. Transport (`target =`) for longitudinal AIPW is **not** supported — owned by the pending Phase 26 design doc.
- **Transport uses a sampling model** \eqn{P(S=1 \mid L)} fit on combined study+target data; weights are \eqn{(1-\hat{p})/\hat{p}} sampling-odds weights. `target_subset = "target"` (transportability) vs `"all"` (generalizability).
- **MTP + transport uses MC marginalization** over \eqn{P(A \mid L, S=1)} because target rows lack observed treatment. Exact enumeration for binary, Monte Carlo integration for continuous. Sandwich variance is not supported for this combination (bootstrap only).
- **Matching + transport is rejected** — matching estimands are fixed at fitting time and cannot incorporate sampling-odds reweighting.
- **Per-component confounders** — `confounders_outcome`, `confounders_treatment`, `confounders_censoring`, `confounders_sampling` (+ TV variants `confounders_tv_outcome`, `confounders_tv_treatment`) allow separate covariate specifications per model component. The old `confounders` / `confounders_tv` are soft-deprecated but still work as a convenience shorthand. No cross-defaults between new args — each model component resolves independently via `%||%` fallback to the deprecated arg.
- **Multiple imputation (`causat_mice()`, Phase 21)** — the analysis step of an MI workflow for MAR covariates/treatment (orthogonal to IPCW, which owns missing Y). It is **estimator-agnostic**: it loops `causat()` + `contrast()` over the completed datasets of a `mice` `mids` object and pools, so every `causat()` configuration (all five estimators, every treatment type / family / intervention / contrast scale, `by`/ATT, stabilized weights, longitudinal ICE + longitudinal IPW, IPCW) works for free. Pooling is **per row independently** (the design contract): means/difference on the natural scale, ratio/OR on the log scale; SNM pools per blip parameter. `pool_method = "rubin"` (default, Barnard-Rubin df, mirrors `mice::pool.scalar()`) is exactly valid under congeniality; under uncongeniality its variance can be biased in *either* direction (Meng 1994; Bartlett & Hughes 2020) — conservative only for certain kinds of uncongeniality, not in general. `pool_method = "boot_mi"` (von Hippel 2020 two-stage) gives a resampling variance that attains nominal coverage **provided the point estimator stays consistent**: a bootstrap fixes variance, not bias, so a misspecified imputation that makes the causal estimate inconsistent (e.g. dropping Y from the predictor matrix, which biases the imputed L) defeats both pooling rules. Per-row pooling diagnostics live in the `"mi_details"` attribute, which `confint()`/`tidy()` read to build `t`-based intervals. **`causat_mice()` does not impute** (the user runs `mice()` upstream) and does not impute Y (Y should be a *predictor* in the imputation model, never imputed — use IPCW for missing Y). Transport (`target =`) + MI is not yet validated.
- **Treatment-free model** (`treatment_free = ~ L`) — SNM efficiency augmentation following Vansteelandt & Joffe (2014). Joint estimation of (β, ψ) absorbs L→Y variance, reducing SEs by 30–45%. Does not change point estimates under correct blip specification.
- **Stratified ICE** (`stratified = "G"`, Phase 22a) — fits one per-step outcome/pseudo-outcome model per level of a **baseline** column `G` in the ICE backward recursion (`ice_fit_step()`/`ice_predict_step()` in `R/ice_stratified.R`). Because `G` is time-invariant and individuals never cross strata, this is **exactly** pooled ICE run per stratum subset — the estimand is unchanged; it is a nuisance-model choice. The within-stratum formula drops the constant `G` term (`strip_stratum_terms()`). **Longitudinal gcomp only**; rejected otherwise (`causatr_stratified_not_ice`), as are time-varying/continuous/NA stratifiers. **Variance is bootstrap-only** — the analytic stacked-EE sandwich (per-stratum × per-time block-diagonal bread) is derived in `PHASE_22_ICE_ENHANCEMENTS.md` but not yet implemented; `ci_method = "sandwich"` aborts (`causatr_stratified_ice_sandwich`).
- **Natural-history MTPs are NOT thin `dynamic()` wrappers** (Phase 22b, designed not shipped) — `grace_period()`/`carry_forward()` depend on the *natural-value history of treatment*, for which the standard ICE recursion is provably wrong (it would set a past treatment to two values at once). They require an augmented-data sequential regression (G-LMTP, arXiv:2605.24167); `dynamic(\(d,trt) d$lag1_A)` runs but silently targets the wrong estimand (observed, not counterfactual, lag). Design + chunk plan in `PHASE_22_ICE_ENHANCEMENTS.md`.
