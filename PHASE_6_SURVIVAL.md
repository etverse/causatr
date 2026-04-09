# Phase 6 — Causal Survival Analysis

> **Status: SCAFFOLDED** (`causat_survival()` fits pooled logistic; `to_person_period()` done; survival curves in `contrast()` pending)
> Book chapter: 17

## Scope

Complete `causat_survival()`, implement survival curves under intervention in `contrast()`, competing risks, ICE integration for longitudinal survival.

## What exists

- `causat_survival()` fits a pooled logistic regression on person-period data (with censoring + competing risks parameters)
- `to_person_period()` converts wide → long format
- `contrast()` can be called on survival fits but treats it as a standard binary outcome (incorrect for survival)
- Vignette: `survival.qmd` exists (basic structure)

## What remains

- [ ] Survival-specific contrast pathway in `compute_contrast()`:
  - Per-individual hazards: ĥ_{i,k}^a = predict(model, newdata = data_a_at_time_k)
  - Per-individual survival: S_i^a(k) = Π_{j≤k} (1 − ĥ_{i,j}^a)
  - Standardised survival: S^a(k) = (1/n) Σᵢ S_i^a(k)
  - Risk at time t: 1 − S^a(t)
  - Risk difference: (1 − S^{a1}(t)) − (1 − S^{a0}(t))
- [ ] Return time-specific estimates (survival curve as a data.table, not a single point)
- [ ] Competing risks: cause-specific hazard models (one pooled logistic per event type), cumulative incidence under intervention
- [ ] ICE + survival integration (Phase 5 dependency): backward iteration with hazard models
- [ ] Variance: sandwich for survival curves (propagate through cumulative product), bootstrap
- [ ] Vignette: `survival.qmd` (file exists; needs full content)

## NHEFS survival replication targets (Ch. 17)

- 120-month survival: ≈ 80.7% under treatment, ≈ 80.5% under no treatment
- Risk difference: ≈ 0.2% (95% CI: −4.1% to 3.7%) — essentially null

## Implementation guide

> To be created from claude.ai summarizing Ch. 17 survival details when this phase begins.
