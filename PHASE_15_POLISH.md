# Phase 15 — Polish and Documentation

> **Status: DONE (2026-05-04)**

## Scope

Documentation, convenience features, and release-prep items.

## Items

- [x] ~~Continuous treatment vignette~~ — covered by `vignettes/interventions.qmd` (shift/scale with g-comp & IPW)
- [x] ~~`target_trial()` metadata/specification object with print method~~ — `R/target_trial.R`
- [x] ~~Documentation: warn about colliders, stepwise/LASSO for confounder selection, point to `dagitty`~~ — callout in `vignettes/introduction.qmd`
- [x] ~~Documentation: ML in g-formula requires debiasing, point to `lmtp`~~ — callout in `vignettes/gcomp.qmd`
- [x] ~~Percent intervened on diagnostic (feasibility metric for each intervention)~~ — `compute_pct_intervened()` in `R/diagnose.R`
- [x] ~~ICE formula builder: support function-transformed TV confounders~~ — `substitute_vars_in_term()` in `R/ice.R`; `confounders_tv = ~ ns(L, 3)` now produces `ns(lag1_L, 3)` at lag expansion
- ~~Grace period / visit-process interventions~~ → deferred to [Phase 22](PHASE_22_ICE_ENHANCEMENTS.md)
- ~~Stratified ICE option~~ → deferred to [Phase 22](PHASE_22_ICE_ENHANCEMENTS.md)
- ~~Multinomial outcomes~~ → deferred to [Phase 23](PHASE_23_CATEGORICAL_OUTCOMES.md)
- ~~Ordinal outcomes~~ → deferred to [Phase 23](PHASE_23_CATEGORICAL_OUTCOMES.md)

## Dependencies

None. Can run at any time. Items can be cherry-picked independently.
Phases 2–6, 8–14 are all complete.

## Out of scope (confirmed across all guides)

| Topic | Reason | Alternative |
|---|---|---|
| G-estimation of SNMMs (Ch. 14) | Rarely used in practice | — |
| Instrumental variables (Ch. 16) | Different identification strategy | `ivreg` |
| TMLE / AIPW / Debiased ML (Ch. 18, 21) | Separate problem class | `lmtp` |
| Forward-simulation g-formula (Ch. 21.6) | ICE is superior | `gfoRmula` |
| Doubly robust sequential (Ch. 21.3) | Out of scope | `lmtp` |
| Causal mediation (Ch. 23) | Different estimands | `mediation`, `medflex` |
| Heterogeneous treatment effects | Different problem | `grf` |
| Sensitivity analysis | Planned as separate etverse package | — |
| Regression discontinuity (RDD) | Different identification strategy | `rdrobust`, `rdd` |
| Difference-in-differences (DiD) | Different identification strategy | `did`, `fixest` |
| Synthetic control | Different identification strategy | `Synth`, `gsynth` |
| Interference / spillover effects | Different assumption set (SUTVA violation) | `interferference` |
| Dynamic treatment regimes (DTR) | Optimal policy learning, not effect estimation | `DTRreg`, `DynTxRegime` |
| Causal discovery / DAG learning | Structure learning, not effect estimation | `pcalg`, `dagitty` |
| Mendelian randomization | Genetic IV, different identification | `MendelianRandomization` |
