# Phase 3 — IPW + Matching + Diagnostics

> **Status: DONE**
> Book chapters: 12, 15

## Scope

IPW estimation via WeightIt, matching via MatchIt, diagnostics (`diagnose()`), and triangulation vignette.

## What was implemented

1. `causat(method = "ipw")` → `fit_ipw()` in `R/ipw.R`
   - Builds treatment model formula from `confounders` term labels
   - Calls `WeightIt::weightit()` for propensity-score weights
   - Calls `WeightIt::glm_weightit()` for the weighted MSM (Y ~ A)
   - Supports external weights multiplication (survey weights)
   - Supports ATE, ATT, ATC estimands (fixed at fit time)
   - Rejects longitudinal data (future: `WeightIt::weightitMSM()`)

2. `causat(method = "matching")` → `fit_matching()` in `R/matching.R`
   - Calls `MatchIt::matchit()` for matched sets
   - Auto-selects `method = "full"` for ATE (nearest-neighbor only supports ATT/ATC); user can override via `...`
   - Extracts matched data + weights via `MatchIt::match.data()`
   - Fits `glm(Y ~ A, data = matched_data, weights = match_weights)`
   - Supports ATE, ATT, ATC estimands (fixed at fit time)
   - Stores `matched_data` and `match_obj` in `causatr_fit`

3. Sandwich variance (method-specific V_β):
   - IPW: `stats::vcov(glm_weightit_model)` — M-estimation vcov accounting for weight uncertainty
   - Matching: `sandwich::vcovCL(model, cluster = subclass)` — cluster-robust on matched pairs
   - Both propagated via `compute_vcov_marginal()` (J V_β Jᵀ)

4. Bootstrap variance (method-specific refit):
   - `refit_ipw()`: resample → re-estimate weights → refit MSM
   - `refit_matching()`: resample → re-match → refit on matched data

5. WeightIt and MatchIt moved from Suggests to Imports

6. `diagnose()` in `R/diagnose.R` — method-specific diagnostics:
   - **All methods**: positivity checks (propensity score distribution, flags near-violations)
   - **IPW**: covariate balance via `cobalt::bal.tab()` (SMD, variance ratios, before/after weighting), weight distribution summary (mean, SD, min, max, ESS by group)
   - **Matching**: covariate balance via `cobalt::bal.tab()` (SMD, variance ratios, before/after matching), match quality summary (n matched, n discarded, % retained)
   - **G-comp**: unadjusted covariate balance via cobalt formula interface
   - Graceful fallback to simple SMD table when cobalt is not installed

7. `plot.causatr_diag()` in `R/plot.R` — Love plot via `cobalt::love.plot()`:
   - Shows absolute SMDs before and after adjustment
   - Threshold line at SMD = 0.1 (customisable)
   - Works for IPW and matching fits

8. Comprehensive test coverage:
   - IPW × binary treatment × binary/continuous outcome × difference/ratio contrasts × sandwich/bootstrap
   - Matching × binary treatment × binary/continuous outcome × difference/ratio contrasts × sandwich/bootstrap
   - IPW/matching × ATE/ATT/ATC estimands
   - Triangulation tests (continuous + binary outcome)
   - Diagnose × gcomp/ipw/matching structure, positivity, balance (cobalt), weights, match quality, print/summary, plot, custom parameters

9. Vignettes:
   - `vignettes/ipw.qmd` — IPW with NHEFS examples + diagnostics section
   - `vignettes/matching.qmd` — matching with NHEFS examples + diagnostics section
   - `vignettes/triangulation.qmd` — method comparison with forest plots + diagnostics

## Key decisions

- **WeightIt/MatchIt as Imports** (not Suggests) — they're core to the triangulation mission
- **cobalt in Suggests** — optional for diagnostics, graceful fallback to simple SMD table
- **IPW currently static-only** — WeightIt doesn't support density ratio weights for dynamic/MTP interventions. Self-contained IPW planned for Phase 4.
- **Predict-then-average** for all methods — `compute_contrast()` is unified across gcomp/ipw/matching. For IPW/matching with saturated MSMs, the Jacobian is trivial so J V_β Jᵀ gives the same result as reading off β₁.
- **diagnose() stores fit reference** — `causatr_diag$fit` stores the original `causatr_fit` so `plot.causatr_diag()` can pass the weightit/matchit object directly to `cobalt::love.plot()`.

## Testing targets

| Test                                                 | Expected                      | Status  |
| ---------------------------------------------------- | ----------------------------- | ------- |
| IPW ATE on simulated DGP (continuous outcome)        | ≈ 3                           | Passing |
| IPW RD on simulated DGP (binary outcome)             | ≈ 0.33                        | Passing |
| IPW risk ratio on simulated DGP                      | > 1                           | Passing |
| IPW ATT estimand                                     | ≈ 3                           | Passing |
| IPW sandwich vs bootstrap SE agreement               | ratio ∈ (0.3, 3.0)            | Passing |
| Matching ATT on simulated DGP (continuous outcome)   | ≈ 3                           | Passing |
| Matching RD on simulated DGP (binary outcome)        | ≈ 0.33                        | Passing |
| Matching risk ratio on simulated DGP                 | > 1                           | Passing |
| Matching sandwich vs bootstrap SE agreement          | ratio ∈ (0.3, 3.0)            | Passing |
| Triangulation: gcomp ≈ ipw ≈ matching (continuous)   | Within ±0.5                   | Passing |
| Triangulation: gcomp ≈ ipw ≈ matching (binary)       | Within ±0.15                  | Passing |
| diagnose() returns causatr_diag for all methods      | correct class + slots         | Passing |
| diagnose() positivity table for binary treatment     | ps summary + violations       | Passing |
| diagnose() balance via cobalt for ipw/matching/gcomp | bal.tab object                | Passing |
| diagnose() weight summary for IPW                    | treated/control/overall + ESS | Passing |
| diagnose() match quality for matching                | n_matched, pct_retained       | Passing |
| plot.causatr_diag Love plot for ipw/matching         | ggplot object                 | Passing |
