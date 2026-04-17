# causatr (development version)

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

## Breaking: `causat()` / `causat_survival()` argument `method` renamed to `estimator`

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
  the mixed "fit-only smoke (Phase 6 contrast pending)" rows with an
  explicit three-table layout (fit path / contrast path / rejection
  path) so the Phase 6 gap is visible without status annotations that
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
  narrative of the current fit-only Phase 6 status, the intended
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
  - Survival × censored fit (smoke; contrasts still Phase 6)
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
