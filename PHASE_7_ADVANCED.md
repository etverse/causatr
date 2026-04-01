# Phase 7 — Advanced Features

> **Status: PENDING**

## Scope

Survey weights, clustered data, parallel bootstrap, multivariate treatment, documentation items, and the `target_trial()` metadata helper.

## Items

- [ ] Survey weights: pass through to model fitting + adjust sandwich
- [ ] Clustered data: cluster-robust sandwich via `sandwich::vcovCL()` (beyond matching subclass)
- [ ] Parallel processing: `future` backend for bootstrap
- [ ] Multivariate treatment: `treatment = c("A1", "A2")` with joint interventions
- [ ] Continuous treatment vignette
- [ ] `target_trial()` metadata/specification object with print method (Ch. 22)
- [ ] Documentation: warn about colliders, stepwise/LASSO for confounder selection, point to `dagitty` (Ch. 18)
- [ ] Documentation: ML in g-formula requires debiasing, point to `lmtp` (Ch. 18)

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
