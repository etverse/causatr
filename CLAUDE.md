# causatr

Unified causal effect estimation for methodological triangulation: g-computation
(parametric g-formula + ICE), IPW (self-contained density-ratio engine), AIPW (doubly-robust), and matching (via MatchIt).
Part of the [etverse](https://github.com/etverse) ecosystem.

## Guide files

- `FEATURE_COVERAGE_MATRIX.md` — **single source of truth for "what works".** Every PR that changes a feature MUST update this file.
- `PHASE_*.md` — per-phase implementation guides in the project root. Completed: 2–6, 8–21 (21a–21e + the longitudinal/IPCW combinations), 22a (stratified ICE), 22b-core (natural-history MTPs / G-LMTPs — `grace_period()` / `carry_forward()`, bootstrap variance), 22b-4 (G-LMTP augmented-data sandwich — `R/variance_if_glmtp.R`), 22b-5 (flexible-treatment ICE term — `treatment_form = ~ factor(A)` / `~ ns(A)`), 22b-6 (`cap_escalation()` public release — ordered-dose increase cap, consistent under `treatment_form`), 25. Pending: **22b-7** (multivariate G-LMTP) — designed with a chunk plan in `PHASE_22` §3.7; 23–24, 26 (design docs); Phase 21 vignette (21i) still to write.

## Project structure

This is an R package: `R/` (source), `tests/testthat/` (tests, `test-foo.R` mirrors `R/foo.R`), `man/` (generated — do not edit), `NAMESPACE` (generated — do not edit), `vignettes/` (long-form docs).

### R/ file layout

**Core API:** `causat.R` (main fitting), `contrast.R` (causal contrasts), `diagnose.R` (main dispatch + panel helpers) + `diagnose_longitudinal.R` + `diagnose_positivity.R` + `diagnose_balance.R` + `diagnose_weights.R` + `diagnose_censoring.R` + `diagnose_intervention.R` + `diagnose_sampling.R` (sampling-model diagnostics for transport).
**Interventions:** `interventions.R` — `static()`, `shift()`, `scale_by()`, `threshold()`, `dynamic()`, `ipsi()`, `stochastic()`; `glmtp_interventions.R` — `grace_period()`, `carry_forward()`, `cap_escalation()` (natural-history MTPs / G-LMTPs).
**Estimation:** `gcomp.R`, `ice.R` (+ `ice_stratified.R` — per-stratum fit/predict helpers for `stratified = "G"`), `glmtp.R` (+ `glmtp_augment.R` — augmented-data sequential regression for natural-history MTPs; support enumeration + tractability guard), `ipw.R`, `aipw.R` (point AIPW), `aipw_longitudinal.R` (longitudinal AIPW), `longitudinal_ipw.R` (longitudinal IPW fitter) + `longitudinal_ipw_formula.R` (per-period propensity/numerator formula builders) + `longitudinal_ipw_contrast.R` (per-intervention contrast bundles), `matching.R`, `snm.R` (SNM fitter + blip spec), `snm_point.R` (point g-estimation + blip design matrix), `snm_contrast.R` (SNM contrast dispatch + averaged blip), `snm_longitudinal.R` (longitudinal backward-sequential g-estimation), `censoring.R` (IPCW model fitting + weights).
**Inference (IF sandwich):** `variance_if_core.R` (model-correction primitives + vcov aggregation), `variance_if.R` (main dispatcher + numeric fallback + point channel + gcomp/matching IF), `variance_if_ice.R` (pooled ICE IF: `ice_if_setup()` Channel-1 + shared state, `variance_if_ice_chain()` per-step Channel-2 cascade) + `variance_if_ice_stratified.R` (per-stratum block-diagonal sandwich), `variance_if_ipw.R` (point + mv IPW IF) + `variance_if_ipw_longitudinal.R` (longitudinal IPW IF) + `variance_if_ipw_sampling.R` (transport sampling-model corrections), `variance_if_aipw.R` (point AIPW IF) + `variance_if_aipw_longitudinal.R` (longitudinal AIPW IF), `variance_if_snm.R` (point SNM sandwich), `variance_if_snm_longitudinal.R` (longitudinal SNM cluster sandwich).
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
| **Interventions** | `static`, `shift`, `scale_by`, `threshold` (gcomp only), `dynamic`, `ipsi` (IPW only — point + univariate longitudinal), `stochastic` (gcomp only; IPW/AIPW when `density` supplied — Phase 24), `grace_period` / `carry_forward` / `cap_escalation` (natural-history MTPs / G-LMTPs — longitudinal gcomp only, discrete treatment, bootstrap + sandwich variance) |
| **Treatment term** | bare numeric (default); flexible per-step ICE term via `treatment_form = ~ factor(A)` / `~ splines::ns(A, df)` (longitudinal gcomp only, incl. G-LMTPs) |
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
- **Stratified ICE** (`stratified = "G"`, Phase 22a) — fits one per-step outcome/pseudo-outcome model per level of a **baseline** column `G` in the ICE backward recursion (`ice_fit_step()`/`ice_predict_step()` in `R/ice_stratified.R`). Because `G` is time-invariant and individuals never cross strata, this is **exactly** pooled ICE run per stratum subset — the estimand is unchanged; it is a nuisance-model choice. The within-stratum formula drops the constant `G` term (`strip_stratum_terms()`). **Longitudinal gcomp only**; rejected otherwise (`causatr_stratified_not_ice`), as are time-varying/continuous/NA stratifiers. **Both variance paths** are supported: ID-cluster bootstrap and the analytic stacked-EE sandwich (per-stratum × per-time block-diagonal bread, `R/variance_if_ice_stratified.R`). The sandwich reuses the pooled engine — `ice_if_setup()` builds the global Channel-1 term + shared state and `variance_if_ice_chain()` runs the per-step correction cascade once per stratum (block-diagonal across strata; Channel 1 stays global because the estimand is marginal). Validated against `delicatessen`'s `MEstimator` to ~1e-7 (the same plug-in M-estimation sandwich) for gaussian and binomial.
- **Natural-history MTPs are NOT thin `dynamic()` wrappers** (Phase 22b, shipped) — `grace_period()`/`carry_forward()` depend on the *natural-value history of treatment*, for which the standard ICE recursion is provably wrong (it conditions on the observed lag, not the counterfactual natural value). They use the **augmented-data sequential regression** of Díaz, Williams, Morzywołek & Rudolph (2026, arXiv:2605.24167): `glmtp_iterate()` (`R/glmtp.R`) decouples the conditioning treatment value from the policy-input value by carrying each natural-treatment-history sequence $\bar{s}_t$ as a label through the backward recursion — the base $m_\tau = E[Y \mid A_\tau, H_\tau]$ is fit once, then each earlier step fits one per-label model ($|\mathcal{A}|^t$ at step $t$, guarded by `glmtp_check_tractable()`) regressing the label-specific pseudo-response on $(A_t, H_t)$ and predicting at the policy-shifted treatment. **Longitudinal gcomp + discrete treatment only**; rejected otherwise (`causatr_glmtp_not_ice` / `_mixed` / `_continuous_trt` / `_too_many`). `grace_period(window)` = the paper's §6 delay-initiation policy; `carry_forward()` = degenerate LOCF (equals a standard-ICE baseline regime, the limiting-case validator); `cap_escalation(delta)` = the paper's dose-escalation cap for an ordered/count dose (public, 22b-6; consistent under `treatment_form = ~ factor(A)`). The three share the augmented path with no subtype dispatch — `glmtp_iterate()` calls each policy closure generically. Reuses `fit_ice` metadata / `ice_build_formula()` / `ice_fit_step()` (with `is_pseudo = TRUE` to muffle count-family integrality warnings at pseudo-steps) / quasibinomial or quasipoisson pseudo family (resolved by `family_pseudo` switch in `fit_ice()`) / ID-cluster bootstrap. Validated against the **forward-MC truth** (the paper's Proposition 1; `lmtp`'s one-shot shift is the *wrong* oracle here — it computes the standard LMTP), replicating Díaz et al. §6. **Variance**: ID-cluster bootstrap and analytic M-estimation sandwich (`R/variance_if_glmtp.R`, 22b-4) — both supported. The sandwich stacks per-(step, label) GLM scores + estimand EE in the same block-triangular structure as the ICE chain; the per-label sensitivity dictionary `D[s_t]` generalises the single `d_vec` of `variance_if_ice_chain()`. Validated against a Python M-estimation oracle (SE ~1e-4) and bootstrap parity (~0.7% at n=1500). The naive `dynamic(\(d,trt) d$lag1_A)` runs but silently targets the wrong estimand (observed, not counterfactual, lag — ~56% off in tests).
- **Flexible-treatment ICE term** (`treatment_form =`, Phase 22b-5) — by default the treatment enters every per-step ICE / G-LMTP outcome model as a **bare numeric** main effect (`ice_build_formula()`), so a nonlinear or kinked counterfactual dose-response (the canonical case: a capped dose) is misspecified. `treatment_form = ~ factor(A)` / `~ splines::ns(A, df)` makes the treatment enter via a transformed term while the intervention/policy still writes the **numeric** treatment column — only the model's design term changes. Lag terms expand by the same parse-tree substitution as transformed TV confounders (`factor(A)` → `factor(lag1_A)`); `treatment_form = NULL` is byte-identical to the old bare-numeric behaviour (the substitution `A` → `lag1_A` reproduces the historical string). Resolved to `treatment_terms` in `fit_ice()` (stored in `details`, raw form re-passed by the bootstrap refit) and consumed by both `ice_iterate()` and `glmtp_iterate()`. **Longitudinal gcomp only**; rejected otherwise via `check_treatment_form()` (`causatr_treatment_form_not_ice`; non-treatment term / non-formula → `causatr_treatment_form_bad`). Under `factor(A)` an intervention must keep the treatment within observed levels (a new level makes `predict()` error); `ns(A)` extrapolates. This is the enabler that closed the `cap_escalation()` bias (gap 0.034 → 0.0015 vs the forward-MC truth), making the now-public `cap_escalation()` (22b-6) a consistent estimator, and unblocks multivariate G-LMTP (22b-7). The analytic ICE sandwich handles the wider design matrix transparently; the G-LMTP sandwich (`variance_if_glmtp()`) handles `treatment_form` transparently via `iv_design_matrix()` / `model.matrix()`.
