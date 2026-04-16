# Phase 7 — Advanced Features

> **Status: PENDING**

## Scope

Survey weights, clustered data, parallel bootstrap, multivariate treatment, documentation items, and the `target_trial()` metadata helper.

## Items

- [ ] Survey weights: pass through to model fitting + adjust sandwich
- [ ] Clustered data: cluster-robust sandwich via `sandwich::vcovCL()` (beyond matching subclass)
- [x] Parallel processing: `boot::boot(parallel=, ncpus=)` for all methods (gcomp, IPW, matching, ICE). `contrast()` accepts `parallel` and `ncpus` params.
- [ ] Parallel processing (optional): `future` backend as alternative to `boot::boot()` built-in parallelism
- [ ] Multivariate treatment: `treatment = c("A1", "A2")` with joint interventions
- [ ] Continuous treatment vignette
- [ ] `target_trial()` metadata/specification object with print method (Ch. 22)
- [ ] Documentation: warn about colliders, stepwise/LASSO for confounder selection, point to `dagitty` (Ch. 18)
- [ ] Documentation: ML in g-formula requires debiasing, point to `lmtp` (Ch. 18)
- [ ] Sequential positivity warnings in `causat()` for longitudinal data (deferred from Phase 5)
- [ ] Stratified ICE option (`causat(..., stratified = TRUE)`) (deferred from Phase 5)
- [ ] Multinomial outcomes: support multi-category outcomes across all methods (g-comp, IPW, matching) via `nnet::multinom()` or `VGAM::vglm(multinomial())`. Requires generalising `compute_contrast()` to handle vector-valued predictions, per-category marginal means, and adapted sandwich/bootstrap variance. Category-specific and pairwise contrasts (risk differences, relative risks) should be supported.
- [ ] Ordinal outcomes: support ordered categorical outcomes via `MASS::polr()` or `ordinal::clm()`. Cumulative probability contrasts and category-specific marginal means.
- [ ] Negative binomial outcomes: already works via `model_fn = MASS::glm.nb` (analytic sandwich path), needs truth-based test coverage. See `PHASE_10_OUTCOME_TYPES.md` for the full plan.
- [ ] Beta regression outcomes: bounded-continuous outcomes in (0,1) via `betareg::betar()` family or `betareg::betareg()` as `model_fn`. Extend `resolve_family()` to accept `"beta"` string, add `betareg` to Suggests. See `PHASE_10_OUTCOME_TYPES.md` for the full plan.
- [ ] Grace period / visit process interventions: for longitudinal data, support interventions that allow missed visits (carry forward last treatment) and censor after exceeding a maximum number of consecutive missed visits (cf. gfoRmula `visitprocess`).
- [ ] Percent intervened on: diagnostic tracking what fraction of the population has their treatment modified under each intervention (feasibility metric, cf. gfoRmula).

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
| Interference / spillover effects | Different assumption set (SUTVA violation) | `inferference` |
| Dynamic treatment regimes (DTR) | Optimal policy learning, not effect estimation | `DTRreg`, `DynTxRegime` |
| Causal discovery / DAG learning | Structure learning, not effect estimation | `pcalg`, `dagitty` |
| Mendelian randomization | Genetic IV, different identification | `MendelianRandomization` |
