# Phase 15 — Polish and Documentation

> **Status: PENDING**

## Scope

Documentation, convenience features, and release-prep items.

## Items

- [ ] Continuous treatment vignette (shift/scale MTP examples with both g-comp and IPW)
- [ ] `target_trial()` metadata/specification object with print method (Ch. 22)
- [ ] Documentation: warn about colliders, stepwise/LASSO for confounder selection, point to `dagitty` (Ch. 18)
- [ ] Documentation: ML in g-formula requires debiasing, point to `lmtp` (Ch. 18)
- [ ] Percent intervened on diagnostic (feasibility metric for each intervention)
- [ ] Grace period / visit-process interventions for longitudinal data (carry forward, censor after N missed)
- [ ] Stratified ICE option (`causat(..., stratified = TRUE)`)
- [ ] ICE formula builder: support function-transformed treatments and TV confounders (e.g. `ns(A, 3)`, `log(L)`) — currently `ice_build_formula()` handles bare column names and `A:modifier` lag expansion, but not `ns(A, 3)` → `ns(lag1_A, 3)`. Requires term-level variable substitution inside parsed expressions. Affects `confounders` with treatment transforms, `confounders_tv` with confounder transforms, and EM interactions with transformed treatments.
- [ ] Multinomial outcomes: multi-category via `nnet::multinom()` or `VGAM::vglm(multinomial())`
- [ ] Ordinal outcomes: ordered categorical via `MASS::polr()` or `ordinal::clm()`

## Dependencies

None. Can run at any time. Items can be cherry-picked independently.

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
