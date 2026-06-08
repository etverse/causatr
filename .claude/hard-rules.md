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

- **Natural-history MTPs (`grace_period()` / `carry_forward()`) are an
  INTERVENTION TYPE, not an estimator.** They return class `causatr_glmtp`
  (inheriting `causatr_intervention`); the estimator stays `gcomp`. There is no
  `estimator = "glmtp"` (considered and rejected, PHASE_22 §3.4 — it would
  duplicate ICE plumbing). When `contrast()` sees a `causatr_glmtp` intervention
  it routes to the augmented-data sequential regression `glmtp_iterate()`
  (`R/glmtp.R`) instead of `ice_iterate()`. Don't flag the absence of a glmtp
  estimator.
- **`glmtp_iterate()` keeps OBSERVED lag columns for conditioning and uses the
  natural-history LABEL `s̄_t` (not the lags) as the policy input** — that
  decoupling is the whole point of the augmented recursion. Don't flag "the
  lags should be counterfactual"; the tower property integrates over the
  observed distribution and the intervention enters only via the current-period
  policy value. Verified by the forward-MC truth test (`test-glmtp.R`).
- **`lmtp` is NOT a valid oracle for genuinely history-dependent G-LMTP
  policies** — its one-shot shift computes the standard LMTP, a different
  estimand. The truth oracle is forward Monte-Carlo of the natural-history
  regime (Díaz et al. 2026 Proposition 1; `helper-glmtp-dgp.R`). Don't add an
  lmtp cross-check for `grace_period()`/`carry_forward()`.
- **G-LMTP sandwich is deferred by design.** `ci_method = "sandwich"` aborts
  with `causatr_glmtp_sandwich`; bootstrap is the only variance path. Don't flag
  the missing sandwich as a bug.
- **`list[[""]]` returns `NULL` in R.** `glmtp_iterate()`'s final `q1 <- Q[[1L]]`
  reads the single empty-label element positionally (its name is `""`, which
  `[[` cannot match) behind a `length(Q) == 1L` guard. Do NOT "simplify" it to
  an empty-string-key access — that silently returns `NULL` (the 2026-06-03
  review caught exactly this).

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
- **Longitudinal AIPW sandwich SE underestimates on unbalanced panels by
  design.** When the panel is unbalanced (monotone dropout / censoring
  row-filter), `variance_if_aipw_longitudinal()` emits a classed warning
  (`causatr_longitudinal_aipw_unbalanced_sandwich`, `.frequency = "once"`)
  and the SE is ~15% low (Monte-Carlo verified, 300 reps). This is a known
  limitation of the row-filtering IF: dropped ids contribute zero to
  later-period channels instead of their unobserved counterfactual, and a
  constant `n / n_period_k` rescale cannot repair it. The bootstrap is
  correct. Do NOT flag the low sandwich SE as a fixable bug — the contract
  is the warning + bootstrap fallback, not a re-derived IF.

### Implementation conventions

- **Weights live in `fit$details$weights`**, never as a column in data.
- **`estimator`, not `method`** — `causat()` uses `estimator = c("gcomp",
  "ipw", "matching")`. `method` is reserved for `...`-forwarding to
  `MatchIt::matchit(method = ...)`. Never rename.
- **Bootstrap must replay `...` via `replay_fit()`.** Any new fit path that
  supports bootstrap must wire through the central refit helper.
- **External oracle cross-checks** use contrast-level (not IF-level)
  comparisons. WeightIt + lmtp as test-only references.

- **SNM `treatment_free` + longitudinal uses the joint (β,ψ) variance.**
  `variance_if_snm_longitudinal_tf()` builds the full per-stage
  `(β_k, ψ_k)` joint system with cross-stage derivatives, then extracts
  the ψ marginal. The psi-only residual `H_k - γ_k` is WRONG when
  `has_tf = TRUE` — the correct residual is `H_k - Z_k β_k - γ_k`. The
  error cancels when `E[R·Z] ≈ 0` (orthogonality), so tests with a
  correctly specified treatment model won't catch it. Test under broken
  orthogonality (TF formula includes covariates absent from PS model).

### Review-time heuristics

- **Saturated MSM neutralizes some variance concerns** — always verify against
  a non-saturated MSM before flagging an IPW/matching variance issue.
- **ICE delta-method derivations.** Before flagging a variance bug in
  `R/ice.R` / `variance_if_ice_one()`, derive the IF decomposition on paper
  and run sandwich vs bootstrap numerically.
- **The matching variance engine uses `cluster = subclass`** for within-
  subclass aggregation. Don't flag this as an unnecessary cluster adjustment.

### Test-DGP conventions (learned the hard way)

- **Person-period outcome placement.** In hand-built longitudinal DGPs with
  interleaved rows (`id = rep(1:n, each = 2)`, `time = rep(1:2, n)`), the
  end-of-follow-up outcome MUST be interleaved onto the final-period rows:
  `Y = c(rbind(rep(NA_real_, n), Y_obs))` (NA at time 1, Y at time 2),
  matching the `A = c(rbind(A_1, A_2))` pattern. The shorthand
  `Y = c(rep(NA_real_, n), Y_obs)` is WRONG — it blocks all NAs into the
  first n rows and all outcomes into the last n rows, scrambling Y against
  the covariates. Tests with this bug still "pass" if their only assertions
  are loose (`expect_gt(est, -0.3)`, sandwich-vs-bootstrap parity), because
  those hold on a meaningless estimate. When reviewing a longitudinal test,
  check the Y construction first.
- **lmtp is the MV oracle too.** `lmtp::lmtp_sdr()` accepts multivariate
  treatment via `trt = list(c("A1_t", "A2_t"), ...)` and a shift function
  over those columns. Use it to pin multivariate longitudinal MTP contrasts
  to truth, not just sandwich-vs-bootstrap parity. causatr's MV longitudinal
  IPW shift matches lmtp to ~3% — it is correct; do not "fix" it.
- **G-LMTP `cap_escalation` + `~ factor(A)`: assert the bare-vs-factor
  *comparison* at n = 40000, not a tight absolute band.** The factor(A) plug-in
  is consistent for the dose-cap policy, but its finite-sample gap to the
  forward-MC truth is ~0.0015–0.010 at n = 40000 (seed-dependent; settles to
  ~0.001–0.0025 only by n = 80000). The bare-numeric bias is seed-stable at
  ~0.03–0.04 (the cap kink). So the headline test (`test-glmtp.R`) keeps the
  affordable n = 40000 and asserts the seed-robust facts (`gap_factor < 0.015`,
  `gap_bare > 0.02`, `gap_factor < gap_bare/2`) rather than a tight 0.005 band
  (which is seed-fragile at n = 40000 and needs n >= 80000 to be robust — too
  slow for CI). Do NOT re-tighten to an absolute 0.005 band at n = 40000, and do
  NOT read a transient n = 40000 factor gap as estimator bias — verify at
  n >= 80000 first.
- **The `cap_escalation` *bootstrap* test (n = 3000) is plumbing-only; the
  consistency oracle is the n = 80000 `~ factor(A)` tight-truth test.** At
  n = 3000 the bootstrap CI half-width (~0.05) dwarfs a consistent estimator's
  gap to the truth, so coverage there cannot distinguish the consistent
  `factor(A)` plug-in from the biased bare-numeric one (the 2026-06-08 review
  verified the *biased* estimator also "covers" and also satisfies a
  `cap < natural` check). The load-bearing assertion in that test is the `1e-8`
  agreement between `contrast()`'s per-arm mean and the direct `glmtp_iterate()`
  engine call — do NOT flag it as redundant. Do NOT add coverage / one-sided
  `cap < natural` assertions back as "consistency" checks (they are vacuous at
  this n), and do NOT move the tight absolute truth band onto the bootstrap test.
