# causatr hard rules

Project-specific rules that override / extend the etverse-wide rules at
`~/Documents/personal/software/etverse/.claude/skills/*/SKILL.md`. Read by the
`implement-feature` and `critical-review-loop` skills before they do anything.

## Project conventions

- **Design-doc pattern.** `PHASE_*.md` at the repo root. When a skill says "read
  the design doc", read the numbered phase doc that covers the feature.
- **Feature coverage file.** `FEATURE_COVERAGE_MATRIX.md`. Every PR that adds,
  removes, or changes a feature MUST update this file.
- **Error-class prefix.** `causatr_*` for all `rlang::abort()` calls.
- **Repro-script prefix.** `/tmp/causatr_repro_<slug>.R`.
- **Test-log paths.** `/tmp/causatr-test-<file>.txt` for per-file runs,
  `/tmp/causatr-test-results.txt` for the full suite.

## Supported dimensions (for combination audits)

| Dimension | Values |
|---|---|
| **Estimator** | gcomp, ipw, matching, ice |
| **Treatment timing** | point, longitudinal (ICE) |
| **Treatment type** | binary, continuous, categorical (k>2), count (Poisson/NB, IPW only), multivariate (gcomp + IPW) |
| **Outcome family** | gaussian, binomial, quasibinomial, poisson, Gamma, quasipoisson, inverse.gaussian, MASS::glm.nb, betareg::betar() |
| **Model class** | GLM, GAM, custom model_fn |
| **Intervention** | static, shift, scale_by, threshold (gcomp only), dynamic, ipsi (IPW only), stochastic (gcomp only, pending) |
| **Estimand** | ATE, ATT, ATC, by-stratified, subset |
| **Contrast type** | difference, ratio, OR |
| **Variance method** | sandwich (analytic IF), bootstrap, numeric Tier 1, numeric Tier 2 |
| **Weights** | none, survey/external, censoring row-filter |
| **Missing data** | complete, Y-missing, L-missing |

## Hard rules (appended to the skill's generic rules)

### Architecture invariants — DO NOT flag these as bugs without a numerical reproducer

- **IPW MSM is `Y ~ 1` per intervention (Hájek intercept)**, not `Y ~ A`. With
  effect modification (`A:modifier` in confounders), it expands to
  `Y ~ 1 + modifier` via `build_ipw_msm_formula()`. The treatment is absorbed
  by the density-ratio weights. Don't add treatment terms to the IPW MSM.
- **Matching MSM is `Y ~ A`** (or `Y ~ A + modifier + A:modifier` with EM via
  `build_matching_msm_formula()`). Bare treatment in confounders (`~ L + A`)
  is rejected by `check_em_compat()`.
- **ICE applies intervention to current-time treatment only.** Lag columns
  hold OBSERVED `A_{k-1}, A_{k-2}, ...` at every backward step. Don't recompute
  lag columns from the intervened treatment — it double-counts shift/scale/
  threshold/dynamic interventions.
- **ICE sandwich under non-uniform weights.** `variance_if_ice_one()` builds
  the step > 0 cascade gradient with the `w_{k-1,j}` factor from the prior
  weighted fit. Dropping this factor underestimates SE by ~2x. The fit-row
  bread × target-row gradient is correct by the delta method; don't flag it
  as a "scope mismatch" without running sandwich vs bootstrap numerically.
- **Matching is binary-only.** MatchIt rejects non-binary; `fit_matching()`
  intercepts upstream with a clear error pointing to `gcomp` / `ipw`.
- **WeightIt is test-only (Suggests).** Never on the runtime path.
- **`dynamic()` is deterministic per-individual, not MTP-as-stochastic.**
  Rejects continuous treatment under IPW by design (no Lebesgue density for
  a Dirac). MTPs in the Díaz/Kennedy sense are `shift()` / `scale_by()` /
  `ipsi()`.
- **Multivariate IPW implements sequential MTP** (Díaz et al. 2023), NOT the
  joint deterministic transformation that multivariate gcomp implements. They
  coincide for static-only interventions and diverge otherwise by design; the
  cross-method divergence test in `test-multivariate-ipw.R` pins this as
  intentional.

### Invariants enforced by code — tests must exercise these, not flag them

- **`na.action = na.exclude` is REJECTED** by `check_dots_na_action()` (error
  class `causatr_bad_na_action`). Only `na.omit` and `na.fail` are accepted.
  Do not try to "harden the pipeline" to accept `na.exclude` — the rejection
  is correct (residuals padding vs model.matrix dropping causes silent IF
  corruption).
- **`censoring =` is a row filter, not IPCW.** No censoring model is fit, no
  IPCW weights are computed. Built-in IPCW is Phase 14.
- **ATT/ATC are only defined for static interventions on binary 0/1 treatment.**
  `check_estimand_intervention_compat()` enforces this with error class
  `causatr_bad_estimand_intervention`.
- **Effect modifier in IPW / matching must be a baseline covariate.** Not
  enforced at runtime (time-varying status isn't inferable from data) — doc-
  level constraint only. Do not flag this as a "missing check".

### Implementation conventions

- **Weights live in `fit$details$weights`**, never as a column in data.
- **`estimator`, not `method`** — `causat()` uses `estimator = c("gcomp",
  "ipw", "matching")`. `method` is reserved for `...`-forwarding to
  `MatchIt::matchit(method = ...)`. Never rename.
- **Bootstrap must replay `...` via `replay_fit()`.** Any new fit path that
  supports bootstrap must wire through the central refit helper.
- **External oracle cross-checks** use contrast-level (not IF-level)
  comparisons. WeightIt + lmtp as test-only references.

### Review-time heuristics

- **Saturated MSM neutralizes some variance concerns** — always verify against
  a non-saturated MSM before flagging an IPW/matching variance issue.
- **ICE delta-method derivations.** Before flagging a variance bug in
  `R/ice.R` / `variance_if_ice_one()`, derive the IF decomposition on paper
  and run sandwich vs bootstrap numerically.
- **The matching variance engine uses `cluster = subclass`** for within-
  subclass aggregation. Don't flag this as an unnecessary cluster adjustment.
