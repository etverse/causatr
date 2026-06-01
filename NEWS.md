# causatr (development version)

## 2026-06-01 — Phase 21: Multiple imputation via `causat_mice()`

`causat_mice()` is now implemented and exported. It is the analysis step of a
multiple-imputation workflow: the user imputes missing covariates and/or
treatment upstream with `mice::mice()`, and `causat_mice()` fits `causat()` +
`contrast()` on every completed dataset and pools the per-imputation estimates
into a single `causatr_result`. This is the supported way to handle
missing-at-random covariates (L) or treatment (A) — orthogonal to IPCW, which
handles missing outcomes (Y).

* **Two pooling engines.** `pool_method = "rubin"` (default) applies Rubin's
  (1987) rules with Barnard-Rubin (1999) degrees of freedom, cross-checked to
  machine precision against `mice::pool.scalar()`. `pool_method = "boot_mi"`
  implements von Hippel's (2020) two-stage bootstrap-then-impute variance,
  which attains nominal coverage without relying on Rubin's congeniality
  assumption — **provided the point estimator stays consistent** (a bootstrap
  corrects the variance, not bias; Bartlett & Hughes 2020). Under
  uncongeniality Rubin's variance can be biased in either direction (Meng 1994),
  so it is not uniformly conservative. Means and difference contrasts pool on
  the natural scale; ratio / odds-ratio contrasts pool on the log scale.
* **Estimator-agnostic.** Every `causat()` configuration flows through
  unchanged: all five estimators (gcomp, IPW, AIPW, matching, SNM), every
  treatment type and outcome family, all interventions and contrast scales,
  `by`-stratification and ATT, stabilized weights, longitudinal (ICE +
  longitudinal IPW), and IPCW (missing L imputed while censored Y is reweighted)
  all pool. SNM pools per blip parameter.
* **S3 integration.** `confint()` and `tidy()` reconstruct `t`-based intervals
  from the per-row pooling degrees of freedom (stored in the `"mi_details"`
  attribute) so the user's confidence level is honored. `parallel = "future"`
  forwards the Boot MI bootstrap to `future.apply`.
* **Guards.** `causat_mice()` warns when an analysis variable is absent from or
  unused by the imputation model (a uncongeniality risk), drops individual
  imputations that fail to fit (aborting only if none survive), and degrades
  gracefully to a single complete-data analysis when `m = 1`. Boot MI warns
  (`causatr_mi_boot_floor`) when a variance component goes non-positive and is
  floored — a signal that `B`/`M` are too small for the random-effects
  decomposition rather than a genuinely tiny standard error.

* **Coverage study + theory correction (21h).** A Monte-Carlo study
  (`test-mi-coverage.R`) confirms Rubin and Boot MI both attain ~nominal
  coverage when the imputation model keeps the causal estimator consistent
  (DGP-MI1), with Rubin's intervals marginally wider. It also corrects an
  earlier claim in the design doc: excluding the outcome from the predictor
  matrix (DGP-MI5) is a *misspecification that biases the estimate*, so neither
  pooling rule recovers coverage — a bootstrap fixes variance, not bias
  (Bartlett & Hughes 2020). The MI workflow guide is the `missing-data.qmd`
  vignette (21i).

New internal engines live in `R/pool_rubin.R` and `R/pool_boot_mi.R`. The
`check_treatment_nas()` hint now advertises MI as a third way to handle missing
treatment alongside IPCW and complete-case analysis.

## 2026-05-31 — Phase 20: Per-period IPSI under longitudinal IPW

Univariate `ipsi(delta)` interventions now work under longitudinal IPW
(`estimator = "ipw"` with `id =` and `time =`). The previously emitted
`causatr_longitudinal_ipsi_pending` rejection in
`compute_ipw_contrast_longitudinal()` has been lifted for the univariate case.

* **No new estimation math**: Kennedy's (2019) closed-form incremental
  propensity-score weight extends to a per-period product
  `W_i = prod_t (delta * A_t + (1 - A_t)) / (delta * p_t + (1 - p_t))`. Each
  period already routes through `compute_density_ratio_weights()`'s IPSI
  branch on the period-`t` subset, and `compute_longitudinal_weights()`
  multiplies the per-period factors within an id exactly as for
  `static` / `shift` / `scale_by`.
* **Stacked sandwich is automatic**: `make_weight_fn_longitudinal()` builds
  per-period sub-closures via `make_weight_fn()`, whose IPSI branch is smooth
  in the propensity coefficients through the `delta * p_t + (1 - p_t)`
  denominator, so `numDeriv::jacobian()` on the stacked closure produces the
  block-diagonal cross-derivative without modification. Sandwich and bootstrap
  SEs agree within Monte-Carlo error.
* **Two combinations stay rejected** with dedicated classes: multivariate IPSI
  (`causatr_longitudinal_ipsi_multivariate` — Kennedy's closed form is
  binary-univariate, the per-component density chain has no IPSI shortcut) and
  IPSI + `stabilize = "marginal"` (`causatr_longitudinal_ipsi_stabilize` —
  the IPSI weight is already a bounded propensity-odds ratio with no separate
  marginal numerator). Longitudinal **AIPW** IPSI remains unsupported (the ICE
  recursion needs a counterfactual treatment stream IPSI does not provide).

## 2026-05-30 — Phase 25: Multivariate longitudinal AIPW

Doubly-robust estimation of joint time-varying interventions now works:
`estimator = "aipw"` with `treatment = c("A1", "A2")`, `id =`, and `time =`.
The previously emitted `causatr_longitudinal_multivariate_pending` rejection
in `fit_aipw_longitudinal()` has been removed.

* **A low-risk composition, no new estimation math**: the outcome (ICE) side
  of the Bang & Robins (2005) backward recursion already loops over vector
  `treatment`, and Phase 19 already built the multivariate per-period density
  chain (`fit_treatment_models()`, `compute_density_ratio_weights_mv()`,
  `make_weight_fn_mv()`). `fit_aipw_longitudinal()` now dispatches on
  `length(treatment)` to fit a per-period K-component propensity chain, and
  `ice_aipw_iterate()` uses the multivariate per-period density-ratio weight.
* **Stacked sandwich**: Channel 2b of `variance_if_aipw_long_one()` stacks
  `T × K` block-diagonal propensity corrections (block-diagonal across time
  via disjoint period subsets, block-diagonal across components via the chain
  rule) instead of `T`. Bootstrap is unchanged — id-clustered resampling
  preserves both axes.
* **Stabilization rejected**: `stabilize = "marginal"` is now an explicit
  classed error (`causatr_stabilize_longitudinal_aipw`) for longitudinal AIPW
  (univariate and multivariate). `fit_aipw()` never threaded `stabilize` into
  the longitudinal path, so a stabilized request was previously honoured
  nowhere; users wanting stabilized longitudinal weights should use
  `estimator = "ipw"`.
* **Truth-based + cross-method tests**: on the static binary
  `make_em_mv_long_scm()` DGP (analytical truths 9 / 15 by sex), MV AIPW
  recovers the stratified truths, agrees with MV ICE g-computation and MV
  longitudinal IPW cross-method, and satisfies one-sided double-robustness
  (wrong outcome + correct propensity, and vice versa). A continuous-shift DGP
  checks finiteness and sandwich-vs-bootstrap SE parity. No external
  single-call oracle exists for longitudinal AIPW (lmtp's SDR uses a different
  nuisance-fixing convention), so triangulation + DR + truth recovery stand in.

### Unbalanced-panel sandwich caveat

The longitudinal AIPW sandwich SE now emits a classed warning
(`causatr_longitudinal_aipw_unbalanced_sandwich`) when the panel is
unbalanced (some individuals are not observed at every period, e.g. monotone
dropout or a censoring row-filter). A Monte-Carlo study (300 reps) showed the
per-period rescaled sandwich underestimates the true SE by ~15% under
informative dropout, because dropped ids contribute exactly zero to
later-period channels rather than their unobserved counterfactual
contribution — a limitation of the row-filtering influence function that a
constant rescale cannot repair. The bootstrap reproduces the truth, so the
warning steers users to `ci_method = "bootstrap"`. Affects both univariate
and multivariate longitudinal AIPW.

## 2026-05-29 — Phase 19c: Effect modification for multivariate longitudinal IPW

Baseline effect modification (`A:modifier` terms in `confounders`) now works
for joint time-varying treatments under `estimator = "ipw"`,
`type = "longitudinal"`. The previously emitted
`causatr_longitudinal_mv_em_pending` rejection has been removed.

* **No multivariate-specific EM code needed**: lifting the rejection lets the
  existing wirings compose. Each per-period, per-component propensity in
  `fit_treatment_models()` strips its own treatment-touching interaction
  terms (Phase 8b), the per-period formula strips them via
  `em_info$confounder_terms` (Phase 10c), and the single final-period MSM
  expands from `Y ~ 1` to `Y ~ 1 + modifier` via `build_ipw_msm_formula()`.
  The modifier stays as a confounder main effect in every propensity; only
  the `A:modifier` interactions are dropped.
* **Baseline-only constraint carries over** (Robins 2000): the modifier must be
  a pre-treatment covariate. A time-varying modifier in the MSM conditions on a
  post-treatment variable. Not enforced at runtime; use `estimator = "snm"` for
  time-varying effect modification.
* **Truth-based test**: `make_em_mv_long_scm()` is a 2-period × 2-component
  binary DGP with sex effect modification whose static contrast has analytical
  truths 9 (sex = 0) and 15 (sex = 1). Tests confirm IPW recovers these, agrees
  with ICE g-computation cross-method, strips the interactions per component,
  and is invariant to stabilization under a static intervention.

## 2026-05-29 — Phase 19b: Stabilized multivariate longitudinal IPW

Stabilized weights (`stabilize = "marginal"`) now compose with joint
time-varying treatments under `estimator = "ipw"`, `type = "longitudinal"`.
The marginal numerator factorises over both axes, mirroring the T × K
denominator:

    g(Ā1, Ā2 | V) = ∏_t ∏_k g_{t,k}(A_{t,k} | A_{t,1..k-1}, Ā_{t-1}, V)

* **Per-period × per-component numerator sweep**: when `stabilize = "marginal"`
  and `length(treatment) > 1L`, each period fits a K-component numerator via
  the denominator `fit_treatment_models()` path, but with confounders set to
  `remove_response(build_longitudinal_numerator_ps_formula(...))` — this drops
  the time-varying L while retaining the treatment lags and baseline numerator
  covariates `V`. The per-component chain rule then layers the prior components
  `A_{1..k-1}` on automatically. Numerator models are stashed as the
  `numerator_models` / `stabilize` attributes on each period's denominator
  `causatr_treatment_models` list.
* **Weight + variance composition**: the MV weight and IF closures
  (`compute_density_ratio_weights_mv()`, `make_weight_fn_mv()`) read the
  numerator attributes and apply the stabilized per-component ratio. Numerator
  parameters (gamma) are held fixed under the sandwich (nuisance-fixed
  convention, as for sigma/theta); bootstrap refits both numerator and
  denominator.
* **Invariants**: under static interventions on discrete treatments the
  stabilized estimate and sandwich SE are *identical* to unstabilized (the
  Hájek-mean invariant carries over from Phase 8e / 10b). The
  `causatr_longitudinal_mv_stabilize_pending` rejection is removed.

## 2026-05-24 — Phase 19a: Multivariate longitudinal IPW

New support for joint time-varying treatments (`treatment = c("A1", "A2")`)
under `estimator = "ipw"` with `type = "longitudinal"`. Composes Phase 8
(multivariate point IPW, K components) with Phase 10 (longitudinal IPW, T
periods) to produce the T × K propensity factorisation:

    W_i = ∏_t ∏_k w_{t,k,i}

* **Per-period multivariate fitting**: each time period fits K per-component
  propensity models via `fit_treatment_models()` with sequential conditioning
  (A_k ~ A_{1..k-1} + confounders + lags). The resulting per-period models are
  `causatr_treatment_models` objects.
* **Weight computation**: `compute_longitudinal_weights()` detects multivariate
  models and dispatches to `compute_density_ratio_weights_mv()` per period,
  then takes the row-product across periods.
* **Sandwich variance**: `make_weight_fn_longitudinal()` delegates to
  `make_weight_fn_mv()` per period, producing a T×K block-diagonal stacked
  alpha. The variance engine (`compute_ipw_if_self_contained_long_one()`)
  iterates over periods, then over components within each period, for the
  per-model corrections.
* **Bootstrap**: works automatically via `refit_ipw()` → `fit_longitudinal_ipw()`.
* **Deferred**: stabilized weights (`stabilize = "marginal"`) and effect
  modification (`A:modifier`) with multivariate longitudinal IPW are rejected
  with classed errors (`causatr_longitudinal_mv_stabilize_pending`,
  `causatr_longitudinal_mv_em_pending`).

## 2026-05-23 — Phase 19-trim: Cross-cutting weight truncation

New `trim` argument on `contrast()` for density-ratio weight truncation
(winsorization). Clips per-component density ratios at the `trim`-th
quantile before they enter any product (multivariate cross-component or
longitudinal cross-period), following the approach of Cole & Hernan
(2008) and the lmtp default of 0.999. Recommended values: `0.999`
(mild) or `0.99` (moderate, per Spreafico et al. 2025).

* `trim = 1` (default): no truncation, all existing behavior unchanged.
* Applied to all weight-computing paths: point IPW, multivariate IPW,
  longitudinal IPW, point AIPW, longitudinal AIPW.
* Sandwich variance: truncation threshold is fixed at the original data's
  quantile and held constant under numDeriv perturbation, so the SE
  reflects the truncated weight surface at a fixed cutoff.
* Bootstrap: each replicate recomputes its own quantile threshold on the
  resampled data, capturing the full variance including truncation.
* NOT applied to sampling weights (transport) or IPCW weights (censoring),
  which define the target population rather than causal reweighting.
* Cross-validated against `lmtp::lmtp_sdr()` with matching `.trim` values
  for both multivariate point IPW and longitudinal IPW.

## 2026-05-24 — Fix: MV IPW sandwich crash with natural-course intervention

* **numDeriv::jacobian crash**: multivariate continuous IPW sandwich
  variance crashed with "subscript out of bounds" when any intervention
  was `NULL` (natural course). Root cause: `make_weight_fn_mv()` returns
  `alpha_hat = numeric(0)` for natural course (no propensity parameters
  to optimize), and `numDeriv::jacobian(fn, x = numeric(0))` cannot
  handle zero-length input. Fixed by adding an early return in
  `compute_ipw_if_self_contained_mv_one()` that skips the
  cross-derivative computation when no propensity parameters exist —
  the propensity correction is zero in this case since the weight
  function is constant (all ones).

## 2026-05-24 — Phase 19-trim critical review fix

* **Sandwich/bootstrap truncation semantics mismatch**: the variance-engine
  closures (`make_weight_fn_longitudinal`, `make_weight_fn_mv`) applied
  truncation to the post-product cumulative weight, while the
  weight-computation path (`compute_longitudinal_weights`,
  `compute_density_ratio_weights_mv`) applied per-component/per-period
  truncation before the product. This caused sandwich and bootstrap SEs to
  disagree by ~10% under aggressive trim values. Fixed by distributing
  per-component/per-period truncation thresholds to each sub-closure,
  matching the per-component semantics throughout.

## 2026-05-22 — Phase 18 critical review fixes

3rd-round critical review of Phase 18 (SNM implementation):

* **Bootstrap metadata fix**: longitudinal SNM bootstrap now correctly
  forwards `confounders_outcome_raw` to replicated fits, ensuring
  metadata consistency on the bootstrap `causatr_fit` objects.
* **Predict consistency**: `snm_treatment_design()` now uses
  `newdata = data` in `predict()` for scalar treatments, matching the
  categorical branch and guarding against future data-vs-model mismatch.
* **File split**: extracted `compute_snm_contrast()` from `snm_point.R`
  (488 lines) into `snm_contrast.R` (191 lines), bringing both files
  under the 300-line guideline.
* **Roxygen**: added missing `@return` tags to
  `variance_if_snm_longitudinal_notf()` and
  `variance_if_snm_longitudinal_tf()`.
* **delicatessen cross-checks**: added two new cross-language oracle
  tests for longitudinal SNM — treatment-free model (TF, no EM) and
  time-varying effect modification with TF (TV-EM + TF). Both use
  Python delicatessen stacked M-estimation as the reference, with
  tight point-estimate tolerances (1e-4) and SE tolerances (0.3–1%).

## 2026-05-22 — Phase 18 complete: Structural Nested Mean Models

Phase 18 is complete. The `snm.qmd` vignette now covers all shipped
features: categorical and count treatment extensions (18h), bootstrap
confidence intervals (18i), and S3 method output (18j), in addition to
the existing coverage of point and longitudinal SNMs with time-varying
effect modification, treatment-free efficiency augmentation, and
sandwich inference.

## 2026-05-22 — SNM-aware S3 dispatch (Phase 18j)

All user-facing S3 methods now recognize SNM results and display
blip-parameter-specific output instead of the intervention-mean layout
used by gcomp/IPW/AIPW/matching.

`print.causatr_result()` shows "Blip parameters:" (point), "Per-stage
blip parameters:" (longitudinal), or "Averaged blip effect:" (Path B
with `treatment_values`), and suppresses the empty contrasts table for
Path A. `summary.causatr_result()` skips the "Intervention details"
section for SNM and labels the vcov as "blip parameters".
`plot.causatr_result()` produces a forest plot with parameter names as
row labels, a zero reference line, and an SNM-specific title.
`coef()` and `tidy()` correctly use the `parameter` column (instead of
`intervention`) present in SNM result objects, and `tidy()` tags rows
with `type = "parameter"`.

## 2026-05-22 — Categorical and count treatment extensions for SNM (Phase 18h)

Extends SNM g-estimation to categorical (k>2) and count (Poisson/NB)
treatment types. The new `snm_treatment_design()` helper abstracts the
construction of treatment-action (AM) and treatment-residual (RM)
matrices so the g-estimating equation solve, variance engine, and
bootstrap are uniform across all treatment types.

**Categorical treatments** use multinomial residualisation via
`nnet::multinom`. The per-level indicators D_j = 1{A=j} and residuals
R_j = D_j - P(A=j|L) produce (K−1) × p_mod blip parameters, one set
per non-reference treatment level. The sandwich variance uses
`prepare_model_if_multinom()` for the multinomial treatment model IF
and a softmax-based cross-derivative closure. `treatment_values` is
rejected for categorical SNM — the per-level blip parameters are the
estimand directly.

**Count treatments** work with zero code changes: the existing scalar
residual path R = A − E[A|L] handles Poisson and negative binomial
treatment models via `propensity_family = "poisson"` or `"negbin"`.

13 new truth-based tests covering point estimation, sandwich variance,
bootstrap consistency, treatment-free model, and longitudinal backward
sequential estimation for both categorical and count treatments.

## 2026-05-22 — Bootstrap variance for SNM (Phase 18i)

Adds bootstrap variance estimation for both point and longitudinal SNM
estimators, completing the dual inference path (sandwich + bootstrap)
that all other estimators already support.

**Point bootstrap** resamples rows and re-solves the g-estimating equation
on each replicate, handling all three paths: standard (no effect
modification), with `A:modifier` interactions, and with `treatment_free`
joint (β, ψ) estimation. When `treatment_values` is supplied, the
bootstrap statistic is the scalar averaged blip effect directly — no
delta-method transformation needed on the bootstrap vcov.

**Longitudinal bootstrap** uses ID-clustered resampling (resample entire
individual trajectories with clone-and-reassign for duplicates) and
re-runs the full backward sequential g-estimation on each replicate.
Forwards per-component confounders (`confounders_outcome`,
`confounders_tv_outcome`, etc.) to ensure the blip specification is
correctly reconstructed on bootstrap replicates even when outcome and
treatment confounders are specified separately.

All bootstrap SE estimates agree with sandwich SEs within a factor of 2
across 9 test configurations (point ± EM ± treatment-free ±
treatment_values ± binary/continuous treatment; longitudinal ± TF ±
TV-EM).

## 2026-05-21 — Triangulation tests, delicatessen cross-check, history=0 support (Phase 18f)

Adds the Robins triangle invariant tests: under correct specification and
no time-varying effect modification, the SNM blip sum, longitudinal
IPW-MSM ATE, and ICE g-comp ATE must agree on binary static interventions.
Three pairwise comparisons (SNM vs IPW, SNM vs ICE, 3-way) plus a negative
control where IPW is biased by collider conditioning on a post-treatment
modifier but SNM recovers the correct blip.

Cross-language validation against Python's delicatessen package (Zivich
2022) using stacked M-estimation on shared fixture data
(data-raw/snm_longitudinal_fixture.csv). Both R and Python operate on
the exact same dataset, enabling tight tolerances: 1e-4 for point
estimates, 0.3% for sandwich SEs.

The `history` parameter now accepts 0 (no-lag Markov model) for all
longitudinal estimators (SNM, IPW, ICE). Previously history < 1 was
rejected. The validator error message updated from "positive integer"
to "non-negative integer." Tests added for all three estimators with
history=0.

Removed 18 suppressMessages/suppressWarnings wrappers around causatr
calls in the test suite (test-snm.R, test-missing-data.R,
test-longitudinal-ipw.R). All underlying messages were either eliminated
at the source or are expected informational output.

## 2026-05-21 — Longitudinal TV-EM truth test + formula dedup fix (Phase 18e)

Adds the headline truth-based test for longitudinal SNMs with time-varying
effect modification. The 2-period DGP has a post-treatment modifier
M_1 = 1{L_1 > 0} (where L_1 depends on A_0) at both stages, with known
blip parameters (psi_00 = 1.15, psi_0M = 2, psi_10 = 2, psi_1M = 2).
Validates that the SNM correctly recovers all four parameters, while
point-treatment IPW conditioning on M_1 is detectably biased due to
collider bias (A_0 -> L_1 -> M_1).

DTRreg cross-check uses `history = 0` for exact treatment-model match.
Cluster bootstrap consistency check and TF-model efficiency test included.

Also fixes `build_longitudinal_ps_formula()` to deduplicate baseline terms
that overlap with time-varying terms. R formulas silently dropped the
duplicate, so the fitted model was unaffected, but the formula object was
misleading to anyone who inspected it.

## 2026-05-21 — Fix longitudinal SNM variance with treatment-free model

The sandwich variance for longitudinal `estimator = "snm"` with
`treatment_free` now correctly accounts for the joint (beta, psi) system
at each stage. Previously, the variance engine used the psi-only
estimating equation scores and bread, ignoring the treatment-free nuisance
parameters. When E[R * Z] ≈ 0 (correctly specified treatment model, TF
covariates subset of PS covariates), the error cancelled exactly. When
E[R * Z] ≠ 0 (TF formula includes covariates absent from the treatment
model), the sandwich SEs were inflated by the unabsorbed TF prediction.
The fix implements the full joint system with per-stage (beta_k, psi_k)
blocks, cross-stage derivatives, and psi-marginal extraction — mirroring
the point-treatment TF variance that was already correct.

## 2026-05-21 — Longitudinal SNMM with per-stage blip estimation (Phase 18d)

`estimator = "snm"` now supports `type = "longitudinal"` for multi-period
panel data. The implementation uses backward sequential g-estimation (Robins
1994): starting at the final stage K, each stage's blip parameters are
estimated using the backward-transformed outcome H_k, yielding per-stage
parameters (e.g. `stage0_psi_intercept`, `stage1_psi_intercept`). Each stage
gets its own treatment model fit via the existing longitudinal propensity
infrastructure.

The cluster-robust sandwich variance accounts for cross-stage dependencies:
earlier-stage blip parameters depend on later-stage estimates through the
backward-transformed outcome, creating off-diagonal blocks in the bread
matrix. Treatment-model corrections are applied per-stage via the standard
numDeriv::jacobian approach.

Testing: truth-based tests on binary and continuous treatment DGPs, DTRreg
cross-checks (binary treatment), vcov PSD checks, rejection tests for
treatment_values and bootstrap. DGP truth values account for the mediated
A_0 -> L_1 -> Y path correctly (validated against DTRreg at n = 500k).

## 2026-05-21 — Time-varying effect modification for SNMs (Phase 18c)

SNMs now explicitly support **time-varying (post-treatment) modifiers** — the
headline feature motivating the addition of SNMMs to causatr. The modifier M
can depend on treatment A, unlike IPW-MSM where this opens collider paths
(Robins 2000). With `treatment_free = ~ L`, the blip parameters recover the
structural truth; without it, the moment-condition estimates capture a
different (still valid) quantity that absorbs the A → M → Y indirect path.

Documentation: updated `parse_effect_mod()` to clarify the SNM exception
to the baseline-only modifier constraint. New vignette `snm.qmd` with a
worked example demonstrating IPW bias with post-treatment modifiers and
correct SNM identification.

Testing: truth-based tests on DGPs where M = f(A, L), delicatessen
cross-checks (4 scenarios: continuous/binary × no-TF/TF), DTRreg cross-check,
IPW bias demonstration, vcov PSD checks.

## 2026-05-20 — Treatment-free outcome model for SNM efficiency (Phase 18b½)

`causat(..., treatment_free = ~ L)` enables a treatment-free outcome model that
reduces blip parameter standard errors by 30–45% without changing point
estimates. The approach follows Vansteelandt & Joffe (2014) and matches DTRreg's
`tf.mod` argument: the treatment-free model parameters β and blip parameters ψ
are solved jointly from a linear system, rather than fitting E[Y | L] separately.

The sandwich variance engine extends to a stacked system
(α_treatment, θ_joint = (β, ψ)) with the joint bread computed analytically and
the cross-derivative A_{θ,α} via `numDeriv::jacobian()`. Validated against
delicatessen on three DGPs (continuous + EM, binary + EM, two modifiers) and
against DTRreg on the EM case — point estimates and SEs match to 0.01
tolerance.

Breaking change: none. The `treatment_free` parameter defaults to `NULL`,
preserving the existing behavior. Non-SNM estimators reject `treatment_free`
with a `causatr_treatment_free_not_snm` classed error.

## 2026-05-20 — SNM point estimation and sandwich variance (Phase 18b)

Point-treatment SNM g-estimation is now fully operational. `contrast(fit)`
returns the blip parameter table (ψ₀, ψ_M, ...) with stacked-EE sandwich
standard errors and confidence intervals. `contrast(fit, treatment_values =
c(0, 1))` computes the population-averaged blip effect with delta-method
variance.

The sandwich variance engine (`variance_if_snm()`) implements the full
stacked M-estimation system: treatment model score block + blip estimating
equation block + cross-derivative A_{ψ,α} via `numDeriv::jacobian()`.
Validated against `delicatessen` (Zivich et al. 2024) on three DGPs
(continuous treatment + single modifier, binary treatment + single modifier,
continuous treatment + two modifiers) with point estimate and SE agreement
to 0.01 tolerance.

New files: `R/variance_if_snm.R` (split from snm.R following the existing
variance file pattern), `data-raw/snm_reference.py` (delicatessen reference
script), fixture CSVs.

## 2026-05-20 — SNM routing and validation (Phase 18a)

Added `estimator = "snm"` to `causat()` for structural nested mean model
(SNMM) g-estimation — the third leg of the Robins triangle alongside
g-computation and IPW-MSM. Chunk 18a ships the routing infrastructure:
`fit_snm()` fits the treatment model via the existing `fit_treatment_model()`
machinery, parses effect-modification terms (`A:modifier` in confounders)
into a linear blip specification, and stores the result for downstream
g-estimation. Supports binary and continuous point treatments.

Rejection paths: multivariate treatment (`causatr_snm_multivariate`),
longitudinal data (`causatr_snm_longitudinal_pending`), and
`contrast(interventions = ...)` (`causatr_snm_no_interventions`) — SNM
estimates blip parameters directly, not counterfactual means.

Point estimation, sandwich variance, and the contrast pathway are deferred
to chunk 18b. `DTRreg` added to Suggests as the primary R oracle for
cross-checks.

## 2026-05-20 — Phase 17 complete: transportability documentation (Phase 17m)

The transportability vignette (`transportability.qmd`) is expanded with
sections for longitudinal transport (IPW and ICE with `target =` and
`id =` / `time =`), multivariate treatment transport (joint-treatment
IPW with sampling weights), separate confounder specifications for the
sampling model (`confounders_sampling =`), and a summary-of-combinations
table. All worked examples include full `causat()` + `contrast()` code
and rendered output.

This completes Phase 17: transportability and generalizability. All
chunks (17a–17m) are shipped. The `CLAUDE.md` guide files section and
`PHASE_17_TRANSPORTABILITY.md` status are updated accordingly.

## 2026-05-18 — MTP + transportability via MC marginalization (Phase 17l)

Modified treatment policy (MTP) interventions (`shift()`, `scale_by()`,
`threshold()`, `dynamic()`, and natural-course `NULL`) can now be combined with
`target = "S"` under g-computation and AIPW. Previously, target-population rows
(S=0) had no observed treatment, so MTP interventions that depend on observed A
produced NA predictions and silently dropped all target rows.

The fix uses Monte Carlo marginalization: for each target row with covariates
L_i, causatr draws treatment values from P(A|L_i, S=1) — the study treatment
distribution — applies the intervention d(A,L), predicts outcomes from the
fitted model, and averages. For binary treatments, exact enumeration replaces
MC sampling. A simple GLM treatment model is fitted internally at contrast time
for g-computation; AIPW reuses its existing propensity model.

IPW + MTP + transport requires no MC marginalization (density-ratio weights
operate on study rows only) and continues to support sandwich variance. Sandwich
variance for gcomp and AIPW + MTP + transport is not yet supported (the MC
marginalization introduces treatment-model dependence not captured by the
current IF chain); `ci_method = "sandwich"` aborts with class
`causatr_mtp_transport_sandwich`, pointing users to bootstrap.

## 2026-05-18 — IPCW × transportability (Phase 17k)

`causat()` now supports `target = "S"` together with `ipcw = TRUE` for all
three estimators (gcomp, IPW, AIPW). The composition produces a
triply/quadruply-weighted estimator: treatment density-ratio weights, stabilized
IPCW weights, and sampling-odds weights all multiply pointwise into the MSM
prior weights (IPW) or the outcome-model fitting weights (gcomp/AIPW). The
sandwich variance engine already handled the IPCW and sampling-model score
blocks independently; this chunk verifies the full four-block composition
(propensity + censoring + sampling + plug-in) works end-to-end.

Truth-based simulation tests validate transportability and generalizability ATE
recovery for all three estimators, sandwich-vs-bootstrap SE agreement, AIPW
double-robustness under outcome misspecification, and cross-estimator agreement
(gcomp ≈ IPW ≈ AIPW). External cross-checks against delicatessen (stacked
M-estimation sandwich) and TransportHealth (gcomp transport with manual IPCW
weights) confirm the composition. A longitudinal IPW + IPCW + transport smoke
test verifies the person-period pathway produces finite estimates.

## 2026-05-18 — Multivariate IPW × transportability (Phase 17j)

`causat()` now supports `target = "S"` together with multivariate treatment
(`treatment = c("A1", "A2", ...)`, IPW estimator). Previously, the weight
computation (joint density-ratio × sampling-odds) was already correct, and the
bootstrap refitted the sampling model per replicate. The only missing piece was
the sandwich variance: the multivariate branch of `variance_if_ipw()` exited
early via `return()` before the sampling-model score block could be appended,
so the gamma uncertainty was absent from the SE. The fix is a single structural
change: store the multivariate IF, conditionally subtract
`compute_ipw_sampling_correction()`, and return the corrected value. The
correction is independent of the number of propensity components.

Truth-based simulation tests verify transportability (target ATE ≈ 3.5 + E[L|S=0])
and generalizability (ATE = 3.5) recovery under binary × binary static
interventions, with SE plausibility against bootstrap. Stabilized weights
(`stabilize = "marginal"`) also produce a finite, positive sandwich SE.

## 2026-05-18 — Longitudinal IPW × transportability (Phase 17i)

`causat()` now supports `target = "S"` together with `id =` / `time =`
(longitudinal IPW). Passing a `sampling_model_fn` alongside the longitudinal
arguments fits P(S=1 | L) on the first period of each individual, broadcasts
the resulting sampling-odds weight `(1 − p̂) / p̂` to all person-period rows for
that individual, and multiplies it into the per-period DR treatment weight
inside the longitudinal Hájek MSM. The result is a transportable (or
generalizable, via `target_subset = "all"`) longitudinal causal estimate.

Sandwich variance extends the Phase 10 stacked EE to include a sampling-model
score block; bootstrap refits the sampling model per replicate. ICE
g-computation with longitudinal transport is also supported: NA lag treatment
columns for target individuals (who have no observed treatment history) are
filled with the static intervention value before the backward ICE iteration, so
counterfactuals can be predicted over the target population. Stabilized weights
(`stabilize = "marginal"`) compose without modification.

Truth-based simulation tests verify transportability (ATE ≈ 4 + E[L|S=0] ≈ 3.67)
and generalizability (ATE = 4) recovery, cross-estimator agreement between
longitudinal IPW and ICE gcomp transport, SE plausibility against bootstrap, and
correct bootstrap refitting of the sampling model.

## 2026-05-18 — External cross-check against TransportHealth (Phase 17h)

Two cross-check tests confirm that causatr's transportability estimates agree
to numerical precision with the independent `TransportHealth` package
(CoreClinicalSciences/TransportHealth, GitHub). `TransportHealth::transportGC()`
implements the same Dahabreh 2020 standardization — fit outcome model on study
rows, average counterfactual predictions over target rows — and agrees with
causatr gcomp to within `1e-4`. `TransportHealth::transportIP()` fits the
same P(S=1|L) sampling model and P(A=1|L) propensity model, combines them into
a weighted MSM, and agrees with causatr IPW to within `1e-4`. Both tests skip
on CRAN; `TransportHealth` is added to `Suggests`. The DGP used is the
standard observational transport DGP with confounded treatment and
treatment-effect heterogeneity (A:L interaction).

## 2026-05-18 — Sampling model predictor-set validation (Phase 17g)

`fit_sampling_model()` now emits a `causatr_transport_predictor_subset` warning
when a baseline covariate appears in the `confounders` formula only via an
interaction with treatment (e.g. `~ L1 + L2:A`). In this situation the
stripping step that removes treatment terms from the sampling formula also
removes the covariate, so the sampling model silently omits a variable that the
outcome and treatment models condition on. The transportability identification
assumption requires the sampling model to control for all baseline confounders;
a misspecified sampling model breaks the exchangeability-over-populations
condition (Dahabreh et al. 2020, assumption 3). The fix is to add the missing
variable as a main effect in `confounders`. The warning fires through
`causat()` as well as when calling `fit_sampling_model()` directly.

## 2026-05-17 — Sampling-model diagnostics for transportability (Phase 17f)

`diagnose()` now includes a **sampling-model panel** when the fit uses
transportability or generalizability (`target =`). The panel reports:
the sampling-score distribution (quantiles of P(S=1|L) overall and per
population), the sampling-weight distribution on study rows (mean, SD,
min, max, ESS), and a count of extreme sampling weights that may signal
positivity violations. The panel appears in every per-intervention
diagnostic and is shared across interventions (it depends only on the
sampling model, not the treatment intervention). `print.causatr_diag()`
renders it automatically.

## 2026-05-17 — AIPW test coverage parity + model_fn default warnings (Phase 16p, 16q)

AIPW now has full test coverage parity with gcomp and IPW for
non-standard outcome families: Gamma, quasibinomial, negative binomial,
and beta regression (betareg). Also covers external weights, survey
design auto-extract, cluster-robust variance, MCAR outcome NAs, subset,
GAM outcome and propensity models, and all compositions thereof.

`causat()` now emits a one-per-session warning when `model_fn` or
`propensity_model_fn` is left at the implicit `stats::glm` default.
The warning encourages deliberate model specification and is silenced
once the user supplies an explicit function. Matching is excluded (it
delegates model choice to MatchIt).

Fixed a cluster-robust variance bug where the cluster vector was not
subset to `fit_rows` before being passed to the AIPW IF engine,
producing misaligned cluster assignments when outcome NAs caused
row-dropping.

## 2026-05-17 — AIPW + IPCW + transport sandwich fix

Fixed an error in `variance_if_aipw()` where the IPCW censoring
correction (Channel 2d) referenced the non-transport variable
`target_fit` unconditionally. When `ipcw = TRUE` and `target =` were
both supplied, `contrast(..., ci_method = "sandwich")` failed with
"object 'target_fit' not found". The censoring correction closure now
branches on transport mode: the transport path varies IPCW weights in
the study-row augmentation term with a fixed target-population
denominator, while the non-transport path retains the original Hájek
ratio.

## 2026-05-07 — Multivariate point AIPW + AIPW-IPCW composition (Phase 16m, 16o)

**Multivariate AIPW (chunk 16m).** `estimator = "aipw"` now supports
multivariate (joint) treatments `treatment = c("A1", "A2", ...)`. The
doubly-robust augmentation term uses per-component density-ratio
weights from the sequential-MTP factorisation. Binary × binary,
continuous × continuous, and mixed-type combinations are supported.
`stabilize = "marginal"` composes correctly (stabilized weights in the
augmentation term match unstabilized point estimates). Rejected for
univariate AIPW (univariate already benefits from Hájek
normalization).

**AIPW + IPCW (chunk 16o).** AIPW now composes with built-in IPCW
(`ipcw = TRUE`) for triply-weighted doubly-robust estimation. The
censoring weights multiply into the augmentation term alongside the
treatment and (optionally) sampling weights. The stacked estimating
equation adds the censoring-model block. Consistent if any two of the
three models (outcome, treatment, censoring) are correctly specified.

## 2026-05-06 — Transportability and generalizability core (Phase 17a–17e)

`causat(target = "S")` enables transportability and generalizability
estimation. The sampling indicator S (1 = study, 0 = target) tags rows
in a mixed-data input; the sampling model P(S = 1 | L) is fit via
`sampling_model_fn` (default logistic GLM) on all rows.

Three estimators transport the estimand to the target population:

- **Gcomp transport**: outcome model fit on S = 1 rows, standardized
  over target-population covariates.
- **IPW transport**: sampling × treatment density-ratio weights multiply
  pointwise in the Hájek MSM on study rows.
- **AIPW transport**: 2-out-of-3 doubly-robust composition (Dahabreh
  et al. 2020) — consistent if any two of outcome, treatment, and
  sampling models are correctly specified.

`target_subset = "target"` (default) targets S = 0 only
(transportability); `target_subset = "all"` targets S = 0 ∪ S = 1
(generalizability). Sandwich variance extends the stacked estimating
equation with the sampling-model block. Bootstrap refits the sampling
model per replicate.

Truth-based simulation tests verify that the transport estimators
recover the target-population ATE on a DGP where study and target ATEs
diverge via an A × L interaction. Cross-estimator agreement (gcomp,
IPW, AIPW) and AIPW efficiency (SE ≤ gcomp ∧ IPW) are pinned.

## 2026-05-05 — Point AIPW estimator (Phase 16a–16j)

New `estimator = "aipw"` implements the classical doubly-robust
estimator: the augmented IPW functional

  ψ(g) = E[m(g(A,L), L)] + E[W(g) · (Y − m(A, L))]

where m is the outcome model and W is the density-ratio weight from
the self-contained IPW engine. Consistent if either the outcome model
or the treatment model is correctly specified.

Supported: binary / continuous / categorical / count (Poisson)
treatments, all contrast types, static / dynamic / shift / scale_by /
ipsi interventions, ATT / ATC / ATE estimands, effect modification
(by =), sandwich and bootstrap variance. Stacked estimating equation
adds the outcome-model and propensity-model score blocks plus the
plug-in row, computed via the existing `prepare_model_if()` and
`apply_model_correction()` primitives.

Cross-validated against `delicatessen` (Python) for both binary static
and continuous shift DGPs — point estimates and SEs match to 3
significant figures.

## 2026-05-05 — Longitudinal AIPW (Phase 16i)

`estimator = "aipw"` with `id =` and `time =` implements the
ICE-AIPW sequential doubly-robust estimator (Bang & Robins 2005). At
each backward step k, the augmentation term corrects the ICE pseudo-
outcome using the period-k density-ratio weight and the period-k
outcome model residual. The cumulative product of per-period
augmentation terms gives the longitudinal DR correction.

Supported: binary / continuous treatments, static / shift / scale_by /
dynamic interventions, gaussian / binomial outcomes, effect
modification (by =), sandwich and bootstrap variance. Double robustness
verified: deliberately misspecifying either the outcome model or the
propensity model still recovers the correct ATE.

Cross-validated against `lmtp::lmtp_sdr()` for binary static and
continuous shift settings.

## 2026-05-05 — Polish: target_trial(), ICE formula builder, feasibility diagnostic (Phase 15)

- `target_trial()` constructor + S3 print method for documenting
  target trial emulation parameters (eligibility, treatment strategy,
  follow-up, outcome, causal contrast, analysis plan).
- `ice_build_formula()` now handles transformed time-varying
  confounders (`poly()`, `ns()`, `log()`, `I()`) via the new
  `substitute_vars_in_term()` helper. Previously, ICE silently dropped
  any term that wasn't a bare column name when building per-step
  formulas.
- `compute_pct_intervened()` reports what fraction of observations are
  affected by each intervention (feasibility diagnostic). Appears in
  the `pct_intervened` slot of per-intervention diagnostic panels.
  Returns NULL for natural-course and ipsi interventions. Longitudinal
  variant reports per-period and overall feasibility.

## 2026-05-04 — Built-in IPCW for MAR outcome censoring (Phase 14)

`causat(censoring = "C", ipcw = TRUE)` now fits an internal censoring
model P(C = 0 | A, L) and applies inverse-probability-of-censoring
weights to the estimator, replacing the previous row-filter-only
`censoring =` behavior when `ipcw = TRUE` is specified.

**Point estimators (14b).** Gcomp, IPW, and matching all receive IPCW
weights. For gcomp and matching, the IPCW weights enter the outcome
model via the `weights` argument. For IPW, the IPCW weight multiplies
into the treatment density-ratio weight. The sandwich variance engine
extends with the censoring-model score block.

**Longitudinal IPCW (14c).** ICE and longitudinal IPW fit per-period
censoring models and apply cumulative IPCW weights. The per-period
censoring probability enters the per-step pseudo-outcome weighting.

**Variance regression (14e).** Sandwich SEs for IPCW fits agree with
bootstrap SEs within the standard MC tolerance (ratio in 0.5–2.0).

**Diagnostics (14f).** `diagnose()` includes a censoring panel with
IPCW weight distribution, censoring prevalence, and censoring-model
positivity quantiles. Longitudinal IPCW shows per-period and
cumulative censoring diagnostics.

Cross-validated against `lmtp::lmtp_sdr()` for both point and
longitudinal IPCW settings.

## 2026-05-04 — Extended outcome type coverage (Phase 13)

`causat(model_fn = MASS::glm.nb)` and `causat(model_fn = betareg::betareg,
family = "beta")` now work end-to-end with sandwich and bootstrap
variance. Previously both errored because `fit_gcomp_point()` always
passed a `family =` argument that neither function accepts.

**Infrastructure changes:**

- New `fn_accepts_family()` helper conditionally strips `family` from
  model args when the fitting function doesn't accept it.
- `resolve_family("beta")` returns a sentinel family object after
  checking that `betareg` is installed.
- `variance_if_numeric()` Jacobian now handles betareg's list-structured
  coefficients (`$mean` + `$precision`).
- `betareg` added to Suggests.

**Expanded test coverage** for all non-standard outcome families
(Poisson, Gamma, quasibinomial, negative binomial, beta regression):
binary and continuous treatments, multiple contrast types
(difference/ratio/OR with rejection tests where invalid), sandwich and
bootstrap SE agreement, IPW cross-checks against g-computation, and
matching smoke tests. ~50 new truth-based simulation tests.

## 2026-05-03 — Stochastic interventions under g-computation (Phase 12)

New `stochastic()` intervention constructor accepts a user-supplied
sampler function and implements Monte Carlo g-formula integration for
both point and longitudinal (ICE) estimators.

```r
sampler <- function(a, data) rnorm(length(a), mean = a - 1, sd = 0.5)
contrast(fit,
  interventions = list(stoch = stochastic(sampler, n_mc = 100)),
  reference = "observed")
```

- Sandwich variance uses MC-averaged influence functions; bootstrap
  draws fresh MC samples per replicate.
- Rejection paths for IPW and matching with informative errors
  (`causatr_stochastic_ipw`, `causatr_stochastic_matching`).
- Cross-validated against `lmtp::lmtp_sdr()` for point and
  longitudinal settings.
- 57 new tests covering binary, continuous, categorical, and
  multivariate treatments across gaussian, binomial, Poisson, and
  Gamma outcome families.

## 2026-05-03 — Suppress non-integer binomial and NB iteration warnings

- IPW outcome MSMs fit with non-integer density-ratio weights now use
  `quasibinomial()` instead of `binomial()` via the new `msm_family()`
  helper (identical coefficients and SEs, no "non-integer #successes"
  warning). Applied in both point IPW and longitudinal IPW.
- Propensity models under external survey weights use `quasibinomial()`
  for the same reason.

## 2026-05-03 — Plot overhaul and diagnostics vignette (Phase 11 chunk 11e)

`plot.causatr_diag()` now supports three plot types via
`which = c("balance", "weights", "positivity")`. The balance plot
(default) remains the cobalt Love plot. Weight-distribution histograms
(`which = "weights"`) use tinyplot, facet by arm for binary treatments,
and accept a `log_scale` argument. Positivity plots
(`which = "positivity"`) show the propensity-score density (binary),
conditional density histogram (continuous/count), or per-level
probability histogram (categorical).

New `vignettes/diagnostics.qmd` provides a user-facing tour of every
supported diagnostic scenario: binary/continuous/categorical IPW,
matching, gcomp, ATT/ATC, `by =` stratification, longitudinal, and
all three plot types.

This completes Phase 11. All `causatr_diag_*_pending` gates from
chunk 11a are lifted, and `diagnose()` now handles every supported
combination of estimator, treatment type, estimand, and timing.

## 2026-05-02 — Treatment-type and longitudinal dispatch for diagnose() (Phase 11 chunks 11b, 11c)

**Chunk 11b — treatment-type dispatch.** `diagnose()` now handles
continuous, categorical, count (Poisson / NB), and multivariate IPW
fits. Continuous IPW reports density-range positivity (flagging
observations in the low-density tails), categorical IPW reports
per-level probability positivity, count IPW reports density-range
positivity, and multivariate IPW reports per-component positivity plus
the combined product-weight distribution. Previously all non-binary
IPW fits hit a `_pending` classed error.

**Chunk 11c — longitudinal dispatch.** Longitudinal IPW and ICE fits
now produce per-period diagnostics: per-period positivity and weight
distribution for longitudinal IPW, per-period covariate balance for
both. The cumulative (product) weight summary flags compounding
instability across periods. Previously longitudinal fits hit a
`_pending` classed error.

## 2026-04-26 — `diagnose()` rewrite foundation (Phase 11 chunk 11a)

`diagnose()` now accepts an `interventions =` argument that mirrors
`contrast()`'s signature, and the returned `causatr_diag` object exposes
a nested `per_intervention` named list (one panel per intervention key
with `positivity`, `balance`, and `weights` slots). The default call
`diagnose(fit)` is unchanged: a single panel keyed `obs` is auto-
injected, populated with the observed-treatment Horvitz-Thompson view
(`1/p` on treated, `1/(1-p)` on controls). Top-level `positivity` /
`balance` / `weights` slots continue to point at the first panel so every
existing flat-shape consumer keeps working without modification.

```r
diag <- diagnose(fit,
                 interventions = list(a1 = static(1), a0 = static(0)))
diag$per_intervention$a1$weights  # weight summary under static(1)
diag$per_intervention$a0$weights  # weight summary under static(0)
```

User-named interventions route through `compute_density_ratio_weights()`
to build the per-intervention HT-style weight; the resulting
treated / control / overall summary uses the same Kish effective sample
size formula as the legacy view, with a divide-by-zero guard so the
zero-weight off-arm under a static-arm intervention reports `ess = 0`
cleanly.

Seven `causatr_diag_*_pending` classed errors gate every sub-feature
that the rewrite has not yet absorbed:
`causatr_diag_continuous_pending` (continuous IPW),
`causatr_diag_categorical_pending` (categorical IPW),
`causatr_diag_count_pending` (Poisson / negbin IPW),
`causatr_diag_multivariate_pending` (multivariate IPW),
`causatr_diag_longitudinal_pending` (replaces the pre-chunk-11a flat
abort on ICE fits), `causatr_diag_estimand_pending` (IPW under ATT /
ATC — gcomp and matching are unaffected because their diagnostic paths
are estimand-agnostic), and `causatr_diag_em_pending` (the `by =`
argument is reserved at the signature level so future chunks can lift
the rejection without a signature break). The chunk 11a balance view
remains the unadjusted SMD across treatment groups; per-intervention
post-weighting balance lands in chunk 11d alongside the ATT / ATC
rewrite.

The 12 pre-existing test blocks in `test-diagnose.R` continue to pass
unchanged as regression anchors; 12 new test blocks cover the per-
intervention dispatch path (with manual `1/p` ESS reconstruction
sanity), the multi-panel `print()` output, and one classed-error check
per `_pending` gate.

## 2026-04-26 — Effect modification under longitudinal IPW (Phase 10 chunk 10c)

`A:modifier` interactions in `confounders =` are now supported under
longitudinal IPW. The per-period propensity formulas strip every
treatment-touching term via `parse_effect_mod()` (modifier main
effects are kept; only the `A:sex` interaction is removed because A
is the response of the propensity model), and the final-period MSM
expands from `Y ~ 1` to `Y ~ 1 + modifier` via the existing
`build_ipw_msm_formula()` so `predict()` returns stratum-specific
counterfactual means.

```r
fit <- causat(d, outcome = "Y", treatment = "A",
              confounders = ~ L0 + sex + A:sex,
              confounders_tv = ~ L,
              id = "id", time = "time", estimator = "ipw")
contrast(fit, list(a = static(1), z = static(0)),
         reference = "z", by = "sex")
```

**Known limitation.** The modifier MUST be a baseline covariate
(Robins 2000). A time-varying modifier in an MSM conditions on a
post-treatment variable, biasing the estimand. Not enforced at
runtime (time-varying status isn't inferable from data); doc-level
constraint only. Use a structural nested model (Phase 18) for
time-varying effect modification.

The chunk 10a `causatr_longitudinal_em_pending` rejection is removed.
Bare treatment in confounders (`~ L + A`) continues to be rejected
upstream by `check_em_compat()` with `causatr_bare_treatment_in_confounders`.

## 2026-04-26 — Stabilized longitudinal IPW (Phase 10 chunk 10b)

`causat(estimator = "ipw", id, time, stabilize = "marginal")` now
fits per-period numerator models $g_k(A_k \mid \bar A_{k-1}, V)$
that drop time-varying confounders from the conditioning set. The
per-period weight becomes $g_k(d(A_k) \mid \bar A_{k-1}, V) \cdot
\lvert \mathrm{Jac} \rvert / f_k(A_k \mid \bar A_{k-1}, \bar L_k)$,
dampening the multiplicative L-dependence across periods that
inflates the cumulative product.

```r
fit <- causat(d, outcome = "Y", treatment = "A",
              confounders = ~ L0, confounders_tv = ~ L,
              id = "id", time = "time", estimator = "ipw",
              stabilize = "marginal", numerator = ~ L0)
```

- Default `numerator = NULL` keeps treatment lags only (drops L from
  every period). Custom `numerator = ~ V` adds the user-supplied
  baseline conditioning set on top of the treatment lags so the
  chain-rule factorisation of $g_k$ stays valid.
- Sandwich variance treats numerator parameters $\gamma$ as fixed
  (same nuisance-fixed convention as Phase 8e for multivariate IPW
  and $\sigma$ / $\theta$ for Gaussian / negbin densities). Bootstrap
  refits both numerator and denominator per replicate and captures
  the full uncertainty.
- Stabilized + static binary recovers identical point estimate **and**
  SE as unstabilized (HT indicator collapses surviving rows so
  numerator and denominator densities coincide).
- Stabilized + `numerator =` without `stabilize = "marginal"` is
  rejected with a `causatr_longitudinal_numerator_without_stabilize`
  classed error (no semantically valid use case for a custom
  numerator without stabilization).

## 2026-04-26 — Longitudinal IPW core (Phase 10 chunk 10a)

`causat(estimator = "ipw", id = ..., time = ...)` now ships:
fits a per-period treatment density chain $f(A_k \mid \bar A_{k-1},
\bar L_k)$ across `time_points`, builds the cumulative density-ratio
weight per individual under sequential-MTP semantics (Robins, Hernán,
Brumback 2000; Díaz et al. 2023), refits an intercept-only Hájek MSM
on final-period outcomes, and reads $\hat\mu_a$ off the intercept.

```r
fit <- causat(d, outcome = "Y", treatment = "A",
              confounders = ~ L0, confounders_tv = ~ L,
              id = "id", time = "time", estimator = "ipw")
contrast(fit, list(s = shift(0.5), n = NULL),
         reference = "n", ci_method = "sandwich")
```

- Supported interventions: `static()` / `shift()` / `scale_by()` /
  `dynamic()` (per-period dispatch through the existing
  `check_intervention_family_compat()` gates).
- Variance: stacked sandwich (`variance_if_ipw_longitudinal()`) with
  block-diagonal propensity bread across periods, plus id-clustered
  nonparametric bootstrap that resamples entire individual
  trajectories.
- Sequential positivity: per-period weight tails are scanned and a
  `causatr_longitudinal_seq_positivity` warning surfaces the offending
  period(s) (default threshold 100).
- Rejected (deferred to follow-up chunks): effect-modification
  (`A:modifier`), `stabilize = "marginal"` / `numerator =`,
  multivariate longitudinal treatment, `ipsi()` per period. Each
  rejection ships with a `_pending` classed error pointing at the
  alternative.
- The longitudinal-IPW rejection in `check_causat_inputs()` is
  removed; matching is the only estimator still gated for longitudinal
  data.

## 2026-04-22 — `future` backend for bootstrap (Phase 9c)

`contrast(parallel = "future")` now dispatches bootstrap replicates
through `future.apply::future_lapply()`, so any active
`future::plan()` is honoured. Users with remote / cluster /
`future.batchtools` workflows no longer have to set up socket
clusters by hand via `parallel = "snow"`.

```r
library(future)
plan(multisession, workers = 4)
contrast(fit, list(a1 = static(1), a0 = static(0)),
         ci_method = "bootstrap", n_boot = 500,
         parallel = "future")
```

- The `"no" / "multicore" / "snow"` paths remain unchanged and still
  go through `boot::boot()`. `"future"` bypasses `boot::boot()` to
  let the `future` plan own the scheduling.
- `ncpus` is ignored under `"future"` — the plan's worker count is
  the single source of truth.
- Reproducibility: `future.seed = TRUE` gives each worker a
  deterministic L'Ecuyer stream, so `set.seed()` at the dispatch site
  governs every replicate.
- `future` and `future.apply` are Suggests-only; `check_pkg()` gates
  the dispatch.

## 2026-04-22 — `survey::svydesign` integration (Phase 9b)

`causat(weights = )` now accepts a `survey::svydesign` object directly.
The sampling weights are extracted via `stats::weights(design)` and the
design's first-stage cluster (PSU) is auto-propagated into the fit's
`cluster` slot so the contrast-time sandwich is cluster-robust by
default.

```r
library(survey)
des <- svydesign(ids = ~psu, weights = ~pw, data = d)
fit <- causat(d, outcome = "Y", treatment = "A",
              confounders = ~ L, weights = des)
# weights and cluster both applied automatically
contrast(fit, list(a1 = static(1), a0 = static(0)))
```

Behaviour:

- Equivalent to `causat(weights = weights(des), cluster = names(des$cluster)[1])`
  — identical point estimate and SE (pinned in `test-survey-design.R`).
- Users can override the auto-adopted PSU by passing `cluster =`
  explicitly; the override wins.
- Single-PSU designs (`svydesign(ids = ~1, ...)`) do not auto-adopt a
  cluster — only the weights flow through.
- Row counts must match: a design built on a larger / different data
  frame aborts with `causatr_bad_svydesign`.
- `survey` is Suggests-only; `check_pkg("survey")` gates the unpack
  path.

## 2026-04-22 — General cluster-robust sandwich (Phase 9a)

`causat()` and `contrast()` accept a new `cluster` argument naming a
column for cluster-robust sandwich variance aggregation. Before this
change the only cluster-robust path was matching's internal subclass
aggregation; users with design clusters (site, household, PSU, repeated
measures) had to hand-roll their own SE.

The implementation threads through the existing influence-function
engine: `variance_if()` now accepts a `cluster_vec` argument that
arrives at `vcov_from_if(cluster = ...)`, which has supported the
sum-within-cluster-then-square Liang & Zeger (1986) aggregation since
the matching branch was wired. Numerically equivalent to
`sandwich::vcovCL(type = "HC0")` applied to the final
predict-then-average step (agreement to within a `G / (G - 1)`
small-cluster correction, pinned in `test-cluster-sandwich.R`).

### Scope

- **Supported**: `estimator = "gcomp"`, `"ipw"`, `"ice"` (including `by`
  stratification and longitudinal ICE; for ICE the cluster vector is
  read from the first-time-point rows).
- **Rejected at fit + contrast**: `estimator = "matching"` already
  clusters on its own `subclass`; layering a user cluster on top would
  double-count within-subclass rows. Classed as `causatr_bad_cluster`.
- **Tier 2 numerical fallback** (no `sandwich::estfun` method and no
  analytic bread): only the Channel-1 block is cluster-adjusted; the
  `J V_beta J^T` block uses the model's standard vcov. Emits a
  classed warning `causatr_tier2_cluster_partial`.

### API

```r
fit <- causat(data, outcome = "Y", treatment = "A",
              confounders = ~ L, cluster = "site")
# cluster auto-propagates from fit to contrast
contrast(fit, list(a1 = static(1), a0 = static(0)))

# or specify at contrast time only
fit <- causat(data, outcome = "Y", treatment = "A", confounders = ~ L)
contrast(fit, list(a1 = static(1), a0 = static(0)),
         cluster = "site")     # column name
contrast(fit, list(a1 = static(1), a0 = static(0)),
         cluster = data$site)  # direct vector
```

### Truth testing

- **gcomp oracle**: direct match against
  `sandwich::vcovCL(fit$model, cluster = d$cl, type = "HC0")` on a
  saturated `Y ~ A + L` model. Agreement to machine precision after
  the `G / (G - 1)` correction.
- **IPW + ICE**: cluster-randomised design (cluster-level treatment
  assignment + cluster-level outcome shock) produces SE inflation
  vs independent aggregation; regression-tested at 1.05x-1.1x lower
  bounds.
- **Singleton clusters**: one cluster per row reduces to the standard
  independent aggregation (tested by equality to the `cluster = NULL`
  path).

---

## 2026-04-22 — Multivariate treatment IPW (Phase 8)

`estimator = "ipw"` now accepts multivariate (joint) treatments
`treatment = c("A1", "A2", ...)`. Previously the path hard-aborted;
users needed `estimator = "gcomp"` for joint interventions. The new
implementation factorises the joint conditional density via the chain
rule

$$
f(A_1, \ldots, A_K \mid L) = \prod_{k=1}^{K} f_k(A_k \mid A_{1..k-1}, L),
$$

fits one univariate density model per component (sequential
factorisation), and forms the product density-ratio weight per
intervention.

### Estimand semantics: sequential MTP (Díaz et al. 2023)

Multivariate IPW implements the **sequential modified treatment policy**
estimand (Díaz, Williams, Hoffman, Schenck 2023, *JASA*), the standard
in the causal-inference MTP literature and what `lmtp` computes. Under
this interpretation, the counterfactual at stage $k$ is $A_k^d =
d_k(A_k^n)$ where $A_k^n$ is drawn from the conditional density of
$A_k$ given the COUNTERFACTUAL history up to stage $k$. At
identification time this plugs observed upstream treatments into the
$k$-th factor's conditioning:

$$
w_k(\alpha_k) = \frac{f_k(d_k^{-1}(A_k) \mid A_{1..k-1}^{\mathrm{obs}}, L; \alpha_k) \cdot |\mathrm{Jac}\,d_k^{-1}|}{f_k(A_k \mid A_{1..k-1}^{\mathrm{obs}}, L; \alpha_k)}.
$$

Under natural course on component $k$ (no intervention), the ratio is
identically $1$ — downstream components do not accumulate an "upstream
conditioning shift" from natural-course components.

**Multivariate gcomp implements a different estimand** — deterministic
joint transformation, $E[Y^d] = E\{E[Y \mid A_1 = d_1(A_1^{\mathrm{obs}}),\, A_2 = d_2(A_2^{\mathrm{obs}}),\, L]\}$
via simultaneous per-individual substitution of each counterfactual
column. For static-only interventions the two estimators coincide (the
indicator on the HT branch collapses the conditioning); for non-static
interventions on non-final components they can differ by the upstream
→ downstream cross-dependence (e.g. shift + shift on a
continuous × continuous DGP with $A_1 \to A_2$ coefficient $0.3$:
gcomp gives $-0.9$, IPW gives $-1.02$). This difference is intentional
and pinned by a regression test.

### Supported combinations

- Binary × binary (and K = 2, 3 binary), all-static, gauss / binom
  outcome, RD / RR / OR contrasts.
- Binary × continuous mixed type, static + shift / scale_by /
  dynamic.
- Continuous × continuous, shift + shift, scale + scale, dynamic
  components.
- Effect modification via `A_k:modifier` (per-component propensity
  formulas strip treatment-touching terms; MSM expands to `Y ~ 1 +
  modifier_main_effects`). Baseline-modifier constraint from Phase 6
  carries over.
- Categorical components (via `nnet::multinom` per component;
  multinomial softmax closure in `mv_ht_closure`). Cat × bin, bin ×
  cat, and cat × cat all supported.
- Count components (Poisson / negative binomial) via per-component
  `propensity_family` opt-in (`NULL` / length-1 / length-K).
  Truth-based tests cover bin $\times$ Poisson, bin $\times$ negbin,
  and Poisson $\times$ continuous (the last exercises the
  upstream $\to$ downstream chain-rule structure under sequential MTP:
  shift($+1$) on $A_1$ propagates to $E[A_2 \mid A_1^d, L]$ through
  the fitted $A_1 \to A_2$ coefficient, giving truth
  $0.3 + 0.4 \cdot 0.2 = 0.38$ on the canonical DGP).
- ATE only (matches the existing single-treatment ATT/ATC gate that
  rejects multivariate at fit time).
- Sandwich + bootstrap variance, `by =`, `subset =`, survey weights.

### Internal bits

`fit$details$treatment_models` is the list of per-component
`causatr_treatment_model` objects (length K, named by treatment
column). The legacy `treatment_model` slot stays populated for
univariate fits. `fit_treatment_models()` in `R/treatment_model.R`
does per-component fitter auto-selection. The mv sandwich
(`compute_ipw_if_self_contained_mv_one()`) exploits that the K
propensity models are fit independently, so the stacked propensity
bread is block-diagonal — the propensity correction sums K
single-model `apply_model_correction()` results.

### Rejected with classed errors

- `ipsi()` in any component (`causatr_multivariate_ipsi`) — Kennedy's
  closed form has no joint analogue.
- Multivariate matching still rejected (MatchIt is binary-only).

### Stabilized weights (Phase 8e)

Multivariate IPW now accepts a `stabilize = c("none", "marginal")`
argument. Under `"marginal"`, the per-component numerator factor is
swapped for a separately-fit density $g_k(A_k \mid A_{1..k-1})$ that
drops the full covariate vector $L$ from the conditioning (Robins,
Hernán, Brumback 2000, extended to multivariate). The per-component
weight becomes

$$
w_k^{s}(\alpha_k) = \frac{g_k\bigl(d_k^{-1}(A_k) \mid A_{1..k-1}\bigr) \cdot |\mathrm{Jac}|}{f_k(A_k \mid A_{1..k-1}, L; \alpha_k)}.
$$

The numerator model's parameters are held fixed in the sandwich
variance (nuisance-fixed, same convention as $\sigma$ for Gaussian
and $\theta$ for negative binomial). Bootstrap refits both numerator
and denominator models on every replicate and captures the full
variance correctly. Implemented for all treatment families the mv
engine already supports (binary, continuous, categorical, Poisson,
negative binomial) and every intervention shape (static, shift,
scale, dynamic, natural course). Rejected with
`causatr_stabilize_univariate` under univariate IPW (univariate
already benefits from the Hájek normalization on the MSM).

Whether stabilization reduces or inflates the sandwich SE depends on
the DGP: when the marginal density of $A_k$ has smaller variance than
the full-$L$ conditional, weights are dampened and the SE shrinks;
when the full-$L$ conditional is tighter than the marginal,
stabilization can inflate SE. This is a known property of MV MTP
weights (Díaz et al. 2023) and a DGP-dependent tradeoff, not a bug.

## 2026-04-20 — Breaking change: survival analysis removed

Causal survival analysis has been extracted from causatr into its own
etverse package. This removes the scaffolded Phase 7 surface — the
`causat_survival()` export, `type = "survival"` branch in `contrast()`,
`PHASE_7_SURVIVAL.md`, the `survival.qmd` vignette, and all survival
tests and snapshots. `to_person_period()` stays in causatr because it
also serves general longitudinal data.

If you were using `causat_survival()`, switch to the new survival
package (`etverse/survatr`). Until then, manual pooled-logistic
regression on person-period data (via `causat(..., estimator = "gcomp")`
plus custom intervention + cumulative-product code) remains possible.

## 2026-04-20 — Bug fixes

- `bread_inv()` now aborts with `causatr_gam_missing_vp` when a
  GAM-classed fit object lacks `$Vp`. Previously it fell through to
  the GLM-style `X'WX` bread on `model.matrix(model)`, silently
  miscomputing the sandwich variance for any caller that built a
  GAM-subclassed object without the Bayesian posterior covariance.
  Properly fitted `mgcv::gam()` objects always carry `$Vp`, so this
  is a defensive guard rather than a regression for standard users.
  Sixth-round critical review Issue L2.

- ICE g-computation now emits a rate-limited, classed
  `causatr_ice_boundary_saturation` warning when any binary-outcome
  prediction lands exactly at 0 or 1. Such rows enter the next
  backward pseudo-outcome model's quasibinomial fit with zero IWLS
  working weight `mu * (1 - mu)` and are silently excluded from that
  fit, biasing the counterfactual chain. The warning is emitted on
  both the final-time model's predictions and each subsequent backward
  step, and is rate-limited via `rlang::warn(.frequency = "once")` so
  long chains and bootstrap loops don't flood stderr.
  Sixth-round critical review Issue L3.

- Tier-2 sandwich variance fallback in `variance_if_numeric()` now tags
  the returned vcov with `attr(., "tier2_approximate") = TRUE`, and
  `print.causatr_result()` surfaces a one-line note when the flag is
  present. Previously the approximate-variance path was only announced
  via a classed `rlang::warn()` that scrolled off quickly and was not
  detectable post-hoc without re-running the variance calculation.
  Users and batch pipelines can now query
  `attr(result$vcov, "tier2_approximate")` to decide whether to prefer
  `ci_method = "bootstrap"`. Sixth-round critical review Issue L1.

- `expand_em_lag_terms()` now correctly finds all treatment components in
  multivariate treatment EM terms (e.g. `A1:A2:sex` with
  `treatment = c("A1", "A2")`). Previously, `which(components == trt_var)`
  relied on R's `==` vector recycling, which silently returned zero matches
  when `treatment_var` was a vector of length > 1 whose length did not
  evenly divide the component count. The function returned `character(0)`,
  dropping all lag EM terms from ICE formulas. Fixed by switching to `%in%`.
  Fifth-round critical review Issue #B1.

- IPW `compute_ipw_contrast_point()` now predicts and averages on
  `fit_data` (outcome-complete rows) instead of the full data. When the
  MSM includes modifier terms (e.g. `Y ~ 1 + sex` from effect
  modification), predictions vary across rows. Previously, the marginal
  mean included NA-outcome rows that the variance engine excluded,
  creating a centering mismatch in Channel 1 of the sandwich IF. Under
  modifier-dependent outcome missingness the sandwich SE was biased
  (9% underestimate at 70% missingness in a repro). For `Y ~ 1` (no
  EM) the bug was invisible because all predictions were constant.
  Fifth-round critical review Issue #B2.

- Removed unused `em_confounder_terms()` helper from
  `R/effect_modification.R`. Added a diagnostic `rlang::warn()` in
  `expand_em_lag_terms()` when no treatment components are found in an
  EM term (previously silent). Fifth-round critical review Issues #S1, #S2.

## 2026-04-17 — Phase 6 complete: unified effect-modification API

Phase 6 (chunks 6a--6e) delivers a unified effect-modification API across
all four estimation methods. Users write `confounders = ~ L + sex + A:sex`
and every method does the right thing:

- **G-comp** feeds the interaction to the outcome model directly.
- **IPW** expands the per-intervention MSM from `Y ~ 1` to `Y ~ 1 + modifier`
  via `build_ipw_msm_formula()`.
- **Matching** expands the outcome MSM from `Y ~ A` to
  `Y ~ A + modifier + A:modifier` via `build_matching_msm_formula()`.
- **ICE** auto-expands the interaction across treatment lags
  (`lag1_A:modifier`, `lag2_A:modifier`, ...) via `expand_em_lag_terms()`.

Cross-method triangulation test confirms that gcomp, IPW, and matching
agree on stratum-specific ATEs within cross-method tolerance on the same
DGP. Vignettes updated with worked examples for all four methods.

## 2026-04-17 — Phase 6 chunk 6d: ICE lag auto-expansion for effect modification

ICE g-computation now auto-expands `A:modifier` interaction terms across
treatment lags. When `confounders = ~ L0 + sex + A:sex`, the ICE outcome
model at each backward step includes not only `A:sex` (current period) but
also `lag1_A:sex`, `lag2_A:sex`, etc. for all available lags. Without this
expansion, ICE captured effect modification only at the current treatment
period, collapsing roughly half the intended heterogeneity.

- New `expand_em_lag_terms()` helper in `R/effect_modification.R` splits
  the interaction term on `:`, substitutes the treatment component with
  `lag{k}_{treatment}`, and returns the expanded terms.
- `fit_ice()` now parses and stores `em_info` from `parse_effect_mod()`.
- `ice_build_formula()` accepts `em_info` and appends lag-expanded EM
  terms to the formula RHS, subject to the same all-NA column validity
  check as other dynamic terms.
- The sandwich variance engine and bootstrap path require no changes:
  both operate generically on the model matrices, which simply grow to
  accommodate the additional interaction columns.
- Tests: 2-period truth (ATE|sex=0 = 5, ATE|sex=1 = 8), 3-period truth
  (ATE|sex=0 = 8, ATE|sex=1 = 12.5), multiple EM terms smoke test,
  bootstrap coverage, regression guard.

## 2026-04-17 — Phase 6 chunks 6a--6c: effect modification infrastructure + IPW/matching MSM expansion

**Phase 6 chunk 6a — `parse_effect_mod()` infrastructure + gate refactoring.**

- New `R/effect_modification.R` with the shared effect-modification parser
  that all four estimation methods consult. `parse_effect_mod()` classifies
  each term in `confounders` as either a pure confounder or an
  effect-modification interaction (`A:modifier`), returning a structured
  `causatr_em_info` object.
- Companion helpers: `build_ipw_msm_formula()`,
  `build_matching_msm_formula()`, and `check_em_compat()` provide
  method-specific MSM construction and validation logic.
- The old blanket-reject `check_confounders_no_treatment()` is refactored
  into `check_confounders_treatment()`, which delegates to the new parser.
  Bare treatment terms (`~ L + A`) are still rejected; true EM terms
  (`A:sex`) now pass through the parser for method-specific dispatch.

**Phase 6 chunk 6b — IPW MSM expansion for effect modification.**

- The IPW MSM expands from intercept-only (`Y ~ 1`) to
  `Y ~ A + modifier + A:modifier` when `parse_effect_mod()` detects an
  interaction term. `build_ipw_msm_formula()` constructs the expanded
  formula; `compute_ipw_contrast_point()` passes it to the per-intervention
  weighted GLM refit.
- The variance engine (`variance_if_ipw`) generalises unchanged: `J`,
  `X_star`, and `phi_bar` all extend from `p_beta = 1` to `p_beta > 1`.
- Tests: DGP-4 truth-based simulation (binary treatment x binary modifier,
  sandwich SE), bootstrap parity, gcomp cross-check.

**Phase 6 chunk 6c — Matching MSM expansion for effect modification.**

- The matching MSM expands from `Y ~ A` to `Y ~ A + modifier + A:modifier`
  when `parse_effect_mod()` detects an interaction term in `confounders`.
  `build_matching_msm_formula()` constructs the saturated MSM; the expanded
  model recovers stratum-specific treatment effects via the matched-sample
  outcome regression.
- `refit_matching()` (bootstrap path) replays the EM-expanded MSM formula
  so bootstrap SEs correctly track the effect-modification structure.
- Formula environment fix: `stats::glm()` evaluates `weights` in the
  formula's environment, but formulas returned from helper functions carry
  the helper's frame. Both `fit_matching()` and `refit_matching()` now
  reset the formula environment to the local frame.
- The variance engine (`variance_if_matching`) generalises unchanged:
  `prepare_model_if()` and `apply_model_correction()` work on the
  expanded GLM's coefficients, score, and bread without modification.
- Tests: DGP-4 truth-based simulation (binary treatment x binary modifier,
  sandwich SE), bootstrap parity, gcomp cross-check, regression guard.

## 2026-04-17 — CRAN compliance: non-ASCII character cleanup

`R CMD check --as-cran` flags non-ASCII characters in `.R` files. Replaced
all occurrences across 22 `R/` files: em/en dashes, arrows, multiplication
signs to ASCII equivalents; accented names to unaccented; Unicode math
symbols to `\eqn{}` LaTeX in roxygen and plain ASCII in comments.

## 2026-04-17 — Documentation overhaul

**Data-generating mechanism descriptions in vignettes.** All vignettes that
use simulated data now include explicit DGM descriptions before each
simulation code block, providing:

- A plain-language description of the causal structure.
- Structural equations (SEM) in display math.
- The known true value of the estimand when analytically available.

Vignettes using the NHEFS dataset now include a description of the assumed
causal structure (treatment, outcome, confounders, DAG) at the start.
Affected vignettes: introduction, gcomp, ipw, matching, interventions,
longitudinal, triangulation, missing-data, validation.

**Vignette coverage expansion.** New sections in existing vignettes to cover
feature combinations that had test coverage but no user-facing
documentation: Gamma outcome and quasibinomial outcome (gcomp, matching),
GLM+splines model specification, IPSI dose-response curve (ipw), 4-period
ICE example (longitudinal), missing-data vignette (new).

**Unicode-to-LaTeX conversion in docs.** Converted Unicode math symbols in
vignettes and `CLAUDE.md` to proper LaTeX for correct MathJax rendering.

**Other documentation changes:**

- `gcomp.qmd` gains a callout note explaining that GAM sandwich SEs use
  the penalised Bayesian covariance (`model$Vp`) as the bread inverse.
- `longitudinal.qmd` gains a truth reference line at `E[Y] = 60` on the
  marginal-means comparison plot for the Table 20.1 null-effect example.
- Variance-theory vignettes enable `tinytable_html_mathjax` for correct
  LaTeX rendering in table captions.
- Numerical tables in `variance-theory.qmd` reduced from 6 to 4
  significant digits for readability.

## 2026-04-17 — `print.causatr_result()` uses `digits = 3`

Console output for `causatr_result` objects now prints with `digits = 3`
(down from 4), matching the `knit_print` method for more compact display.

## 2026-04-16 — Phase 4 chunk 3j: count treatment IPW (Poisson / negative binomial)

`causat()` gains a `propensity_family` parameter (`"poisson"` or
`"negbin"`) that enables IPW estimation on integer-valued count
treatments. Auto-detection never infers count (non-negative integers
like age are legitimately Gaussian), so this is an explicit opt-in.

**Runtime:**
- `fit_count_density()` fits either `stats::glm(family = poisson())` or
  `MASS::glm.nb()` (auto-selected for `"negbin"` when
  `propensity_model_fn = NULL`).
- `evaluate_density()` gains Poisson (`dpois`) and NB (`dnbinom`)
  branches.
- `make_weight_fn()` gains a count closure using the log link
  `lambda = exp(X %*% alpha)`.
- `check_intervention_family_compat()` allows `shift()` (integer delta
  only) and `scale_by()` (inverse must preserve integer support) on
  count treatments; rejects `static`, `dynamic`, `threshold`, `ipsi`.

**Variance:** no changes needed — Poisson and NB GLMs inherit from
`glm`, so `bread_inv()` and `prepare_model_if()` work out of the box.
NB `theta` is treated as fixed (same convention as Gaussian `sigma`).

**Tests:** truth-based Poisson DGP (`shift(1)` recovers ATE = 1.5),
NB parity (NB nests Poisson on the same DGP), 6 rejection snapshot
tests for disallowed intervention types.

## 2026-04-16 — Fifth-round critical review S2: near-zero intervened density warning

`compute_density_ratio_weights()` now warns (class
`causatr_near_zero_intervened_density`) when > 80% of observations have
near-zero density at the intervened treatment value under `shift()` or
`scale_by()`. This catches the case where the intervention pushes
treatment far outside the fitted distribution's support, producing
density-ratio weights that are effectively zero for the entire sample.
The Hájek estimator still runs (the math is valid — zero weights just
mean no contribution), but the resulting estimate is degenerate and
should not be interpreted. Previously this situation was silent.

## 2026-04-16 — Fifth-round critical review R1: `diagnose()` censoring alignment

`diagnose()` now correctly excludes censored rows when computing
positivity, balance, and simple-balance diagnostics for gcomp fits that
use a `censoring` column. Previously, `compute_positivity()`,
`compute_balance()`, and `compute_balance_simple()` called
`get_fit_rows(data, outcome)` without passing `fit$censoring`, so
diagnostics were computed on the full non-NA-outcome sample (including
censored individuals) while the main gcomp pipeline excluded them. The
diagnostics thus described a different sample than the one used for
estimation. All three sites now call `get_fit_rows(data, outcome,
fit$censoring)`. Regression test in `tests/testthat/test-diagnose.R`.

## 2026-04-15 — Fourth-round critical review: reject `na.action = na.exclude`

`causat()` and `causat_survival()` now abort via `check_dots_na_action()`
when `na.action = na.exclude` is forwarded through `...`. Under
`na.exclude`, `stats::residuals(model, "working")` is padded with NAs to
the original data length while `model.matrix(model)` drops NA rows — the
length mismatch propagated into `prepare_model_if()`'s `r_score` and
then silently corrupted the Channel-2 correction via R's recycling (only
a "longer object length" warning, no abort). Sandwich SEs were
mathematically wrong. Only `na.omit` (default) and `na.fail` are
accepted now; users with NAs should drop them explicitly before calling
`causat()`. New error class `causatr_bad_na_action`. Regression test in
`tests/testthat/test-causat.R`.

## 2026-04-15 — Third-round critical review: dots audit + `replay_fit()` helper + T6/T7/T9

A full audit of the `fit$details$dots` plumbing (capture sites in
`fit_gcomp_point()` / `fit_ipw()` / `fit_matching()` / `fit_ice()`,
replay sites in `refit_gcomp()` / `refit_ipw()` / `refit_matching()` /
`ice_iterate()`) flagged three systemic risks:

- **R1 (T2)**: positional dots with empty names slipped through the
  hand-written `setdiff(names(dots), ...)` strips.
- **R2 (C5)**: four duplicate-strip blocks had to stay in sync by hand.
- **R4**: the strip lists were drifting — adding a new reserved key to
  one site but forgetting the others was a silent regression.

**Fix — Alternative E in the audit: one central helper.** New
`replay_fit(fn, base_args, dots, reserved)` in `R/utils.R` owns the
`c(base_args, dots) + do.call(fn, ...)` composition. All four refit
sites (plus `fit_ipw()` itself) now call it. Behavior:

- Unnamed (positional) entries in `dots` are dropped. Fixes T2.
- Keys already present in `base_args` are stripped from `dots` so the
  caller's explicit value always wins. Covers C5 centrally.
- A `reserved` vector lets callers block additional keys the target
  function uses implicitly (e.g. `WeightIt::weightit`'s `s.weights`
  when external weights are not supplied).

Other non-dots fixes in the same pass:

- **T6** — `causat_survival()` now uses `stats::quasibinomial()` when
  external weights are supplied, suppressing the spurious
  "non-integer #successes in a binomial glm" warning without
  changing coefficients or SEs (same score equations, free dispersion).
- **T7** — `bread_inv()`'s singular-bread warning now carries a
  `causatr_singular_bread` class. The bootstrap `withCallingHandlers`
  in `variance_bootstrap()` matches on class, not on the English
  substring, so a future wording change cannot silently break the
  demotion.
- **T9** — `causat_survival()`'s internal-column strip now sources the
  list from a new `CAUSATR_SURVIVAL_INTERNAL_COLS` constant in
  `R/utils.R` (derived as `setdiff(CAUSATR_RESERVED_COLS, ".pseudo_y")`),
  so adding a new reserved column propagates automatically.

Audit items R3 (unknown dots silently ignored by target function), R5
(backward-compat with old `fit` objects missing `$dots`), R6 (replay
opacity), and R7 (arg aliasing across target functions) were verified
as non-bugs via direct R scripts — base R catches unknown args at fit
time, `%||% list()` handles missing `$dots`, `fit$details$dots` is
inspectable, and explicit `base_args` locks override dots.

**Tests.** New `tests/testthat/test-replay-fit.R` adds 14 unit tests
pinning `replay_fit()` edge cases (base-wins precedence, positional
dots dropping, reserved-key blocking, NULL dots, `glm.control` catching
unknown args, backward-compat with missing `$dots`, end-to-end gcomp
and IPW bootstrap replay). All 44 earlier critical-review +
weight-edge-case tests still pass unchanged.

## 2026-04-15 — Second-round critical-review follow-ups (C3, C5, C12) + B7 gradient fix + weight edge-case tests

A second-round self-review of the 2026-04-15 critical-review sweep
flagged twelve concerns (C1–C12). C1, C2, C4, C6, C7, C8, C9, C10, C11
were either false alarms (misread R semantics) or acceptable as-is.
C3, C5, and C12 landed as real fixes; running the regression tests
also surfaced an incomplete B7 fix. A new
`tests/testthat/test-weights-edge-cases.R` adds 24 weight-specific
tests covering every branch the critical review flagged as
weight-sensitive.

- **C3 — Classed abort for empty target population.**
  `build_point_channel_pieces()` now emits
  `rlang::abort(class = "causatr_empty_target", ...)`, and the `by`
  skip path in `compute_contrast()` matches on class, not on the
  English message text. Any future wording drift can no longer
  silently turn skipped strata into real aborts.

- **C5 — Duplicate key stripping in `c(args, dots)` compositions.**
  `fit_ipw()`, `refit_ipw()`, `refit_gcomp()`, and
  `ice_iterate()` now strip keys already present in `args` from the
  stashed `dots` before the `c()`. `do.call()` always takes the first
  named duplicate, so the previous behavior was correct but fragile —
  explicit stripping makes the precedence unambiguous.

- **C12 — `.causatr_prev_event` / `.causatr_prev_cens` stripped
  from `fit$data`.** The within-id risk-set bookkeeping columns in
  `causat_survival()` are now dropped before the `causatr_fit` is
  returned, so they no longer leak into user-visible state.

- **B7 gradient follow-up — unify ICE step_i == 1 gradient on the
  weighted form.** The B7 fix in the first sweep only touched the
  Channel-1 IF vector, not the corresponding sensitivity gradient
  computed at step 1 of the backward iteration, which still used
  `n_target` (now removed from the enclosing scope). Running the new
  regression test surfaced a `n_target not found` abort. The gradient
  is now computed via the same `w_t * mu_eta / sum_w_target` form
  and reduces to the original unweighted expression when `w_t = 1`.

- **Weight edge-case coverage.** New test file exercises: NA / Inf /
  NaN / negative / non-numeric / wrong-length weights rejected by
  `check_weights()`; zero weights as row exclusion; uniform weights
  producing identical SEs to the unweighted path (gcomp + ICE);
  heterogeneous weights changing SEs non-trivially; IPW survey
  weights attaching as `s.weights`; matched outcome model carrying
  `match_weight * external_weight`; `causat_survival()` weight
  plumbing to the hazard fit.

## 2026-04-15 — Critical-review sweep (B1-B8, R3-R12)

A second full-codebase critical review flagged ~twenty correctness and
robustness issues across the variance engine, bootstrap refit path,
`contrast()` dispatch, and `causat_survival()`. Every blocking fix
ships with a regression test in
`tests/testthat/test-critical-review-2026-04.R`.

### Blocking correctness fixes

- **B1 — `subset` expressions now resolve caller-scoped variables.**
  `contrast(fit, subset = quote(age > cutoff))` previously evaluated
  `subset` with `parent.frame()` deep inside internal dispatch frames,
  so `cutoff` (defined in the user's session) failed to resolve. The
  caller's frame is now captured once at `contrast()`'s top level and
  threaded through `compute_contrast()`, `get_target_idx()`,
  `variance_bootstrap()`, and `ice_variance_bootstrap()` as an explicit
  `subset_env` argument. `get_target_idx()` also validates that the
  evaluated expression has length `nrow(data)` so a scalar `TRUE`
  cannot recycle silently.

- **B2 — Bootstrap replicates replay the user's `...` verbatim.**
  `refit_gcomp()`, `refit_ipw()`, and `refit_matching()` previously
  dropped the user's extras (`method = "cbps"`, `stabilize = TRUE`,
  `caliper`, `ratio`, `gamma`, GAM smoothing, etc.), so bootstrap SEs
  for any non-default estimator silently corresponded to a different
  estimator than the point estimate. `fit_ipw()` / `fit_matching()` /
  `fit_gcomp_point()` now stash `list(...)` onto `fit$details$dots`
  and each refit function threads those back through `do.call()`.

- **B5 — `causat_survival()` implements the documented censoring rule.**
  The roxygen promised "subsequent rows for that individual are also
  dropped", but the implementation only excluded the current row. A
  `.causatr_prev_cens` within-id cumsum + lag column now drops all
  rows at and after the first censor per id.

- **B6 — External weights enter WeightIt as `s.weights`.** `fit_ipw()`
  post-multiplied external weights onto `w$weights` *after*
  `WeightIt::weightit()` had finished, so the Mparts sandwich
  correction silently under-corrected survey-weighted IPW. External
  weights now flow in as `s.weights`, and `refit_ipw()` mirrors the
  change at bootstrap time.

- **B7 — ICE Channel-1 IF unified on a single formula.**
  `variance_if_ice_one()` previously used `(n / n_target) * (Y - mu)`
  unweighted and `n * (w / sum_w) * (Y - mu)` weighted. These agree
  only when `sum(w) == n_target`, so heterogeneous weights mis-scaled
  the Channel-1/Channel-2 cross-term. The unweighted branch now sets
  `w = 1` and shares the weighted formula. A regression test pins
  uniform-weighted and unweighted SEs to machine precision.

- **B8 — `by` skips empty strata instead of aborting the whole call.**
  `by_vals` is now enumerated from `fit$details$fit_rows` (not the
  full data), and any per-level call that still ends up with an empty
  target is caught by a `tryCatch`, collected into a single
  `rlang::warn()` at the end, and dropped from the combined result.

### Required robustness fixes

- **R3 — Branch-A propensity IF reconciles via `na.action`.** When the
  MSM and propensity fits drop different rows,
  `prepare_propensity_if_weightit()` now intersects `fit_idx` with the
  MSM's `na.action` instead of aborting.

- **R6 — OR validation no longer aborts on NA `mu_hat`.** The
  `mu_hat >= 1` check now pre-filters NAs consistently with the
  `<= 0` check.

- **R7 — `bread_inv()` prefers `weights(model, type = "working")`.**
  Keeps GLM subclasses (`glm_weightit`, glmnet's glm wrapper)
  correctly aligned with the sandwich convention.

- **R8 — Matched cluster aggregation uses a first-seen factor.**
  `vcov_from_if(..., cluster = subclass)` no longer sorts the cluster
  levels alphabetically / numerically, which had desynchronised
  `rowsum(..., reorder = FALSE)` from `IF_mat`'s row order.

- **R9 — `causat_survival()` sorts once up front via `setkeyv()`.**
  Previously two `by = id` aggregations relied on data.table's
  group-order staying consistent between calls.

- **R10 — Robust Mparts detection.** Handles both list-shaped and
  flag-shaped `attr(w, "Mparts")` across WeightIt versions.

- **R11 — Dead `is.integer` branch removed in `interventions.R`.**
  `is.integer(x)` already implies `is.numeric(x)`.

- **R12 — Reserved column names centralised.** New constant
  `CAUSATR_RESERVED_COLS` and helper `check_reserved_cols()` in
  `R/utils.R`, called from `prepare_data()`, `fit_ice()`, and
  `causat_survival()`. The previously-unprefixed `prev_event`
  internal column was renamed to `.causatr_prev_event`.

### Suggestions addressed

- **S1** — `subset` length validated in `get_target_idx()`.
- **S5** — `Reduce(`&`, ...)` in `compute_contrast()` gets an explicit
  `init = rep(TRUE, nrow(data))` so an empty `preds_list` falls
  through cleanly.

## 2026-04-15 — Documentation overhaul: vignettes, math rendering, and `knit_print`

A pass over every package vignette plus a small package-level addition
so that `causatr_result` objects render as readable HTML tables inside
Quarto / knitr documents.

### New package feature: knitr-aware `causatr_result` rendering

- **`R/knit_print.R`** defines a `knit_print.causatr_result()` S3 method
  that renders a `causat() |> contrast()` result as a metadata header
  plus two tinytable HTML tables (intervention means + contrasts) when
  the result object is evaluated inside a knitr chunk. The method
  degrades gracefully to the existing `print.causatr_result()` if
  tinytable or knitr is unavailable.
- **`R/zzz.R`** conditionally registers the method via `.onLoad()`
  (`requireNamespace("knitr", quietly = TRUE)`), so the dependency on
  knitr stays Suggests-only.
- **`DESCRIPTION`** adds `tinytable` and `litedown` to `Suggests` (the
  latter is needed by `tinytable::format_tt(markdown = TRUE)`, which the
  variance-theory vignette uses for math-in-cells rendering).

### Vignette fixes (no behavioural changes to the engine)

- **`longitudinal.qmd` — Table 20.1 null-effect section rewritten.**
  The previous explanation claimed "individual effects of A_0 and A_1
  are both zero" and that naive regression gives non-zero estimates for
  both. In fact (i) the true *marginal counterfactual* means are all
  equal to 60, not the conditional effects, (ii) `lm(Y ~ A0 + A1 + L1)`
  gives `A1 = 0` exactly (in the DGP Y does not depend on A1), and
  (iii) the bias shows up in `A0` only, as the conditional-within-L1
  effect of $-8$. The section now shows *both* naive strategies
  (unadjusted and L1-adjusted) failing, demonstrating the null paradox
  as stated in Hernán & Robins Chapter 20.

- **`by = "sex"` examples fixed across three vignettes.**
  - `gcomp.qmd` now refits with an explicit `qsmk:sex` term in
    `confounders`, so the stratified estimates actually move (Male
    $\approx$ 3.61, Female $\approx$ 3.44).
  - `ipw.qmd` and `matching.qmd` replace the `by = "sex"` code chunks
    with a Phase 8 callout explaining that IPW and matching wrap a
    saturated `Y ~ A` MSM and so `by = "sex"` is a no-op until the
    unified effect-modification API lands
    (`PHASE_6_INTERACTIONS.md`).

- **Unified "Summary of covered combinations" tables.** gcomp, ipw,
  matching, and longitudinal now share one 8-column schema
  (Treatment | Outcome | Intervention | Estimand | Contrast | Inference
  | Weights | Status) with a common ✅ / 🟡 / ⛔ legend and a pointer
  to `FEATURE_COVERAGE_MATRIX.md`.

- **`code-fold: show` + `code-tools: true` in every vignette YAML.**
  Combined with the new `knit_print` hook, readers can fold the R code
  of any example chunk and still see a nicely formatted result table.

### `variance-theory.qmd` — roadmap + verified numerics

- New **Section 0 "Roadmap"** walks through Sections 1 $\to$ 6 as a single
  narrative and introduces the running example (a small logistic GLM)
  used by every verification block below.
- Eight new **Math / Code / Result** tabset panels, each of which
  instantiates a key equation on the running example and checks it
  against either `sandwich::sandwich()`, `sandwich::vcovCL()`, or
  `causatr::variance_if()`:
  1. **Sec 1.4** — Sandwich variance of $\hat{\beta}$ vs `sandwich::sandwich(fit)`
     (agreement ≲ 1e-8).
  2. **Sec 2.3** — $(1/n^2) \sum IF_i IF_i^\top$ equals the sandwich to machine
     precision.
  3. **Sec 2.4** — IF for a sample mean reproduces the textbook
     $s/\sqrt n$ up to the $(n-1)/n$ d.f. correction.
  4. **Sec 3.2** — Two-channel IF for g-comp $\hat\mu_1$ matches
     causatr's sandwich SE to $\sim$1e-9.
  5. **Sec 3.3** — Delta-method decomposition: shows numerically that
     $J V_\beta J^\top$ equals $\text{Ch.2}^2$, and that dropping $\text{Ch.1}^2$ + cross term
     underestimates the SE by about 3% on this example.
  6. **Sec 4.2** — Horvitz-Thompson IPW $\hat\mu_1$ equals causatr's
     IPW point estimate, with SE from the propensity-corrected
     adjusted score.
  7. **Sec 4.3** — Cluster-robust sandwich: per-cluster IFs summed
     before squaring reproduce `sandwich::vcovCL(..., cadjust = FALSE)`
     on a toy subclass assignment.
  8. **Sec 5.5** — ICE chained IF on a 2-period linear-Gaussian DGP
     recovers the structural truth $0.8 + 0.6 = 1.4$ within the
     sandwich CI.

### `validation.qmd` — tinytable comparisons

- `compare_row()` now returns a tinytable with a third row of
  **absolute differences** between causatr and each reference
  implementation (stdReg2, WeightIt, marginaleffects, lmtp). Agreement
  — or the lack of it — is readable without having to compare 6
  decimal places by eye.

## 2026-04-17 — Breaking: `causat()` / `causat_survival()` argument `method` renamed to `estimator`

`causat(method = ...)` shadowed `WeightIt::weightit(method = ...)` and
`MatchIt::matchit(method = ...)`, making it impossible to forward those
packages' own `method` selector through `...`. The causatr argument has
been renamed to `estimator`; extra `...` arguments now flow cleanly into
WeightIt / MatchIt. This is a hard rename with no deprecation shim —
update call sites from `method = "gcomp" | "ipw" | "matching"` to
`estimator = "..."`, and from `fit$method` / `result$method` /
`diag$method` to `fit$estimator` / `result$estimator` / `diag$estimator`.
The `glance.causatr_result()` output column is renamed in the same way.

After the rename, forwarding e.g. a CBPS propensity model is a one-liner:

```r
fit <- causat(data, outcome = "Y", treatment = "A", confounders = ~ L1 + L2,
              estimator = "ipw", method = "cbps")
```

See `vignette("ipw")` for a worked example. `ci_method` is a separate
concept and is unchanged.

## 2026-04-14 — Critical review sweep (R1–R9 + S1–S9)

Systematic audit and hardening of the whole package in a single pass.
All fixes below are covered by existing or newly-added tests; the
FEATURE_COVERAGE_MATRIX.md now carries a CI-enforced round-trip check
so future drift between the matrix and the test files fails loudly.

### Required changes (correctness / silent-failure fixes)

- **`contrast()` ratio / odds-ratio guard (R2).** The log-scale delta
  method was silently returning `NaN` CIs when `mu_hat` had mixed signs
  (legal for Gaussian outcomes, fatal for `log(R)`). `contrast(type =
  "ratio")` now aborts upfront if any intervention mean is non-positive,
  and `type = "or"` additionally aborts if any mean is $\geq 1$. Error
  messages point to `type = "difference"` or a binomial / poisson /
  gamma refit.
- **`variance_bootstrap()` warning suppression narrowed (R3).** The
  refit closure previously wrapped the whole pipeline in
  `suppressWarnings()`, hiding GLM non-convergence, bread-singular
  warnings, and everything else. Replaced with a
  `withCallingHandlers()` that demotes only the three expected noisy
  warnings (`fitted probabilities numerically 0 or 1`, `X'WX is
  singular`, `Fewer (control|treated) units`). Genuine failures now
  surface through `tryCatch`'s error path and show up as failed
  replicates in `boot_info`.
- **Shared match-weight alignment helper (R4).** Extracted
  `combine_match_and_external_weights()` (in `R/matching.R`) so both
  `fit_matching()` and `refit_matching()` share the rowname-invariant
  check. A future change to `as.data.frame.data.table`'s rowname
  behavior now fails loudly in exactly one place instead of silently
  producing misaligned weights in the bootstrap.
- **`dynamic()` docstring corrected (R5).** The previous wording
  claimed the rule received data "subset to the current time point";
  in reality it is called once per intervention with the full
  counterfactual panel. Updated the roxygen prose and the example to
  show how to branch on the time column or a lag column.
- **ICE all-NA column drop is no longer silent (R6).**
  `ice_build_formula()` still drops `lag*_` columns at step 0 (a
  structural requirement), but now emits a `rlang::inform()` when a
  user-supplied column is dropped at a step because it happened to be
  all NA (e.g. a time-varying confounder first measured at time > 0).
  Gated via `.frequency = "regularly"` so the message does not spam
  long pipelines.
- **`boot_info$n_fail_by_int` (R7).** `process_boot_results()` now
  records per-intervention NA counts in addition to the whole-row
  totals. This lets users see whether bootstrap failures cluster on a
  single intervention (e.g. a rare `static()` level triggering
  separation) versus being spread across the whole replicate. The
  field is surfaced in `boot_info` and tested in `test-s3-methods.R`.
- **`causat_survival()` wide-format check (R8).** The previous
  `max(rows_per_id$N) == 1L` test only caught uniformly wide data; a
  mixed frame with some single-row ids and some multi-row ids would
  slip through. Now uses `any(rows_per_id$N == 1L)` and additionally
  asserts that `time` varies within each id (no duplicated rows at
  the same time value).
- **`FEATURE_COVERAGE_MATRIX.md` ↔ `tests/` round-trip check (R9).**
  New `tests/testthat/test-coverage-matrix.R` parses the matrix,
  extracts every `test-*.R` reference, and asserts both that the
  referenced file exists and that every test file on disk is
  referenced somewhere in the matrix. Guards against the exact kind
  of drift the rest of this sweep was cleaning up.

### Suggestions (hardening / clarity)

- **`fit_gcomp_point()` resolves `family` via `resolve_family()` (S1)**
  so custom `model_fn` implementations that don't auto-coerce a bare
  string see a fully-evaluated family object. Matches the behavior
  already in `fit_ipw()` and `fit_matching()`.
- **Multivariate sub-intervention name validation (S2).**
  `apply_intervention()` now asserts `setequal(names(iv), treatment)`
  before iterating the sub-interventions. A typo like
  `list(X = static(1), A2 = static(0))` for `treatment = c("A1",
  "A2")` now aborts with a clear message instead of silently creating
  a spurious `X` column via data.table.
- **`ipsi()` constructor informs about its Phase-4 status (S3).** The
  constructor succeeds (so users can build their interventions list)
  but surfaces a one-per-session `rlang::inform()` pointing to the
  known `contrast()` abort, so users hit the wall at construction
  rather than deep in a pipeline. Use `.frequency_id =
  "causatr_ipsi_dead_end"` to silence or grep for it.
- **`check_estimand_trt_compat()` error message (S4)** now explicitly
  says "binary point treatments coded as 0/1" and suggests recoding.
  Factor- or `1/2`-coded treatments previously got the same message
  without any hint at the coding issue.
- **Tier 2 numeric variance warning is now classed (S5).**
  `variance_if_numeric()`'s delta-shortcut warning carries
  `class = "causatr_tier2_fallback"` so batch pipelines and CI can
  grep for it via `withCallingHandlers()`.
- **`compute_weight_summary()` binary split pulls levels from
  WeightIt (S6)** instead of hard-coding `trt == 0` / `trt == 1`.
  Factor-coded binary treatments and integer codings like `c(1, 2)`
  now get correctly-populated "treated" / "control" rows.
- **`contrast()` subset evaluation preserves lexical scope (S7).** The
  three `eval(subset, envir = as.list(data))` call sites were copying
  every column of `data` into a fresh list and losing the caller's
  enclosing environment. Switched to
  `eval(subset, envir = data, enclos = parent.frame())`, which skips
  the copy and lets `subset` expressions reference session variables
  like `quote(age > cutoff)`.
- **FEATURE_COVERAGE_MATRIX survival section split (S8).** Replaced
  the mixed "fit-only smoke (Phase 7 contrast pending)" rows with an
  explicit three-table layout (fit path / contrast path / rejection
  path) so the Phase 7 gap is visible without status annotations that
  mean different things on different rows.
- **`.pseudo_y` column collision guard in ICE (S9).** `fit_ice()` now
  aborts upfront if `data` already has a column named `.pseudo_y`,
  which is the internal bookkeeping column ICE writes predicted
  pseudo-outcomes into during the backward iteration.

### Doc drift — `A:modifier` behaviour (R1)

`CLAUDE.md`, `NEWS.md` (this file's older entries), `PHASE_6_INTERACTIONS.md`,
and `CAUSATR_SCAFFOLD.md` all still described IPW and matching as
"silently ignoring" `A:modifier` interaction terms in `confounders`.
That code path was hardened months ago: `fit_ipw()` and `fit_matching()`
now call `check_confounders_no_treatment()` and **abort** upfront with
a Phase-8 pointer. Every doc that mis-described the old behavior has
been updated to match the current behavior.

### Vignette expansion

Every method vignette (excluding `variance-theory.qmd`) now covers
every combination supported by `FEATURE_COVERAGE_MATRIX.md`:

- **`gcomp.qmd`** — added worked examples for categorical (k > 2)
  treatment, multivariate (joint) treatment, external survey weights,
  and a Poisson count outcome with a log-link ratio contrast. Summary
  table now records the family / weights columns.
- **`ipw.qmd`** — added categorical and continuous (GPS) treatment
  sections, an external-survey-weights section, a parallel-bootstrap
  section, and a non-Mparts WeightIt fallback example. Updated the
  Channel-2 technical note to reference the current
  `prepare_propensity_if()` / `apply_propensity_correction()`
  architecture (was still referring to the renamed
  `correct_propensity()`).
- **`matching.qmd`** — added an ATT + external survey weights section
  demonstrating `combine_match_and_external_weights()`, and extended
  the summary table with a categorical / continuous rejection row.
- **`longitudinal.qmd`** — added continuous-treatment + MTP section
  (shift / scale_by / threshold), multivariate time-varying
  treatment, external weights (with a note pointing to the ICE
  cascade-gradient correction in §5.4 of the theory vignette), and a
  censoring / IPCW-style example.
- **`survival.qmd`** — replaced the 3-line TODO stub with a full
  narrative of the current fit-only Phase 7 status, the intended
  end-to-end contrast algorithm, and the tracking pointers.
- **`triangulation.qmd`** — introduction.qmd and the comparison table
  both corrected for the drift above (matching is binary-only;
  continuous / categorical are not "Limited").

## 2026-04-14 — Sweep and audit

- **ICE sandwich fix under non-uniform external weights.** The
  cascade gradient at step k > 0 was dropping the `w_{k-1}` factor
  from $A_{k-1,k} = E[\partial s_{k-1} / \partial \beta_k]$, silently underestimating the
  SE by ~2x on multi-period DGPs. Verified against `lmtp::lmtp_tmle`
  and a 300-replicate Monte-Carlo. (Followed by a separate fix for
  the earlier regression where the ICE backward loop dropped
  external weights entirely.)
- **ICE double-counting fix for `shift`/`scale_by`/`threshold`.**
  `ice_apply_intervention_long()` was recomputing lag columns from
  the intervened treatment, double-counting the shift through both
  the lag path and the current-A prediction path. Now matches
  `lmtp` on the canonical linear-Gaussian DGP.
- **New `FEATURE_COVERAGE_MATRIX.md`.** Single source of truth for
  what's tested and at what fidelity (truth / smoke / none /
  rejected). New test-writing rule in `CLAUDE.md`: every feature
  change MUST update the matrix in the same PR.
- **Coverage sweep** closing most residual gaps:
  - gcomp × quasibinomial, Gamma(log), probit, cloglog, Poisson ×
    bootstrap, categorical (3-level) treatment
  - IPW × categorical (3-level), continuous, survey weights, ATT
    bootstrap, non-Mparts WeightIt method at contrast time
  - matching × Poisson, Gamma(log), ATC bootstrap, categorical
    rejection (MatchIt is binary-only; verified against upstream
    docs)
  - ICE × `scale_by`, `threshold`, `shift` bootstrap, bootstrap ×
    survey weights, binomial × ratio/OR, `by(sex)` stratification,
    external-weights truth test
  - Survival × censored fit (smoke; contrasts still Phase 7)
  - End-to-end Tier 2 numeric variance via a custom `model_fn`
  - `vcov_from_if(cluster = ...)` unit test vs hand-computed
    reference
  - `confint()` level monotonicity on both sandwich and bootstrap
- **Defensive guards across `causat()`, `causat_survival()`, IPW,
  matching, ICE, diagnose.** Front-door validation for NA / Inf /
  negative / mis-sized / non-numeric `weights`; rejection of
  duplicated or empty intervention names; explicit longitudinal
  rejection in `diagnose()` (previously crashed inside cobalt);
  categorical treatment rejection in `matching` with a clearer
  causatr-side message.
- **`~ sex + A:sex` no longer triggers a false-positive** "Baseline
  confounder 'A' varies within some individuals" warning. Treatment
  columns are now excluded from `warn_confounder_variation()`.
- **altdoc sidebar fix.** Single-item sections ("Getting started",
  "Theory") were silently dropped from the deployed site because
  `yaml::yaml.load()` collapses a one-element scalar sequence into a
  scalar that Quarto's sidebar can't parse. Both sections now use
  the explicit `text:`/`file:` form, making the `variance-theory`
  vignette visible on the website.
- **`vignettes/variance-theory.qmd` Section 5.4** now derives the
  ICE cascade gradient with the explicit `w_{k-1}` factor and the
  pseudo-code shows the `prior.weights` lookup.
- **`PHASE_6_INTERACTIONS.md` design doc.** Plans a unified
  effect-modification API across gcomp / IPW / matching / ICE. IPW
  and matching currently hardcode a saturated `Y ~ A` MSM and
  **hard-abort** on any `A:modifier` term in `confounders` via
  `check_confounders_no_treatment()` (rather than silently dropping
  it and returning a pooled ATE); ICE handles current-period A ×
  modifier but not lagged A × modifier. Phase 8 will replace both
  the abort and the lag-gap with a proper MSM builder.

## 2026-04-13 — Variance engine unification

- **Single influence-function engine.** Replaced four method-specific
  variance paths with one `variance_if()` dispatcher that routes to
  `variance_if_gcomp()` / `variance_if_ipw()` /
  `variance_if_matching()` / `variance_if_ice()`. All four share
  `prepare_model_if()` / `apply_model_correction()` primitives and
  `vcov_from_if()` aggregator.
- **Prep / apply split.** `prepare_model_if()` computes the
  gradient-independent half of Channel 2 once per model;
  `apply_model_correction()` finishes it per intervention. Saves an
  O(k) bread solve when `contrast()` is called with many
  interventions.
- **Matching uses `cluster = subclass`** in the aggregator
  (cluster-robust IF-level sum-then-square).

## 2026-04-12 — Large refactor + feature pass

- **Point-treatment sandwich now uses the full influence function**
  (previous version was effectively a delta-method shortcut and
  could underestimate the SE on non-saturated MSMs).
- **`variance_if()` as single entry point.** Introduces
  `correct_model()` and `correct_propensity()` as the two per-model
  correction kernels. New `vignettes/variance-theory.qmd` derives
  the two-channel IF from first principles.
- **Log-scale CIs for `ratio` and `or` contrasts.** Respects the
  (0, ∞) support and gives better coverage than Wald bounds.
- **`static()` accepts character values** for categorical
  treatments (e.g. `static("never")`).
- **Analytical Jacobian for GLMs** in sandwich variance (previously
  used numerical differencing).
- **External weights moved to `fit$details$weights`** (never a
  column in `data`) so every method consumes them uniformly.
- **New `test-complex-dgp.R`** covering nonlinear confounding,
  GAMs, splines, and multi-confounder designs.
- **ICE IF scaling fix** for subset / by-stratified targets, plus a
  `confint()` ordering fix and a vectorized IF loop for speed.

## 2026-04-11 — Bootstrap percentile CIs

- **Bootstrap path stores raw replicates** in `object$boot_t` and
  `confint()` recomputes percentile CIs from them per `level`. The
  sandwich path continues to reconstruct Wald CIs from the stored
  SE.
- External weights now propagate through every method (gcomp, IPW,
  matching, ICE); previously ICE dropped them in the backward loop.

## 2026-04-10 — Polish and hardening

- **`by`-stratification `vcov` fix** — stratum-specific vcov
  matrices now land in `res$vcov[[level]]` as a list instead of
  being collapsed to the pooled matrix.
- **`family` now threads through IPW and matching** (previously
  silently forced Gaussian on the MSM regardless of the user's
  `family` argument).
- **`scale()` $\to$ `scale_by()` rename** to avoid colliding with the
  base R `scale()`.
- **`by + estimand` bug fix** — previously the ATT/ATC restriction
  and the by-level restriction fought over the target population
  and returned the wrong per-stratum estimate.
- **`type` parameter on `causat()`** for explicit point /
  longitudinal selection; autodetects from `id` + `time` otherwise.

## 2026-04-09 — Effect modification, S3 polish, broom integration

- **New `by` argument on `contrast()`** for stratified (by-level)
  effect estimates. Computes separate target means and vcov per
  stratum; `confint()` and `tidy()` respect the stratification.
- **broom integration** via `generics` — `tidy.causatr_result()`
  and `glance.causatr_result()` return the expected columns.
- **Forest plots** via the `forrest` package for `causatr_result`
  objects.
- **S3 method polish** across `print`, `summary`, `plot`, `coef`,
  `vcov`, `confint`, `tidy`, `glance`.

## 2026-04-02 / 2026-04-03 — Longitudinal scaffolding + R CMD check

- ICE sandwich variance with scaling fix; first batch of
  simulation-based ICE tests.
- Clean R CMD check: man pages regenerated, MASS / knitr declared,
  Windows vignette skip to work around a toolchain issue.

## 2026-04-01 — Phase 2–3 implementation

- **`causat(estimator = "gcomp" / "ipw" / "matching")`** fully
  working, with:
  - parametric g-formula via a user-supplied `model_fn`
  - IPW via `WeightIt::weightit()` and `WeightIt::glm_weightit()`
  - matching via `MatchIt::matchit()` + weighted MSM
- **`diagnose()`** for all three methods: positivity, balance (via
  `cobalt`), weight summaries, match quality, Love plots.
- **First vignettes**: `gcomp.qmd`, `ipw.qmd`, `matching.qmd`,
  `triangulation.qmd`.

## 2026-03-30 / 2026-03-31 — Phase 1 scaffolding

- Initial package structure: `DESCRIPTION`, `NAMESPACE`, `R/`
  layout, testthat setup, dev tooling, GitHub Actions CI.
- Core API stubs: `causat()`, `contrast()`, `diagnose()`,
  intervention constructors (`static`, `shift`, `scale_by`,
  `threshold`, `dynamic`, `ipsi`), S3 classes.
- NHEFS dataset bundled; input-validation helpers in `checks.R`.
