# Phase 22 — ICE Enhancements: Stratified ICE + Grace Period Interventions

> **Status: PENDING** (design doc)
>
> **Depends on:** Phase 5 (longitudinal ICE), Phase 15 (ICE formula builder fix for transformed TV confounders)

## Scope

Two ICE-specific enhancements that modify the backward iteration loop in
`ice_iterate()` and add longitudinal-aware intervention constructors.

### 22a — Stratified ICE (`causat(..., stratified = TRUE)`)

Fit separate outcome models per stratum of a baseline variable at each
backward step. Motivation: when the outcome–treatment relationship
differs structurally across subgroups (e.g., different functional forms
by sex), a single pooled model is misspecified. Stratified ICE fits
S models per time step instead of one.

**Key changes:**
- Add `stratified` parameter to `causat()` (character: column name, or `NULL`)
- `ice_iterate()` splits fitting data by stratum, fits per-stratum models, merges predictions
- Models stored as `list(time = list(stratum = model))` instead of `list(time = model)`
- Variance: per-stratum IF contributions assembled into the stacked EE sandwich
- Bootstrap: refit per-stratum models at each replicate

### 22b — Grace period / visit-process interventions

New intervention constructors for longitudinal settings:

- `grace_period(inner, n)` — apply `inner` intervention, but tolerate up to
  `n` consecutive periods of deviation from target before enforcing. Wraps
  `dynamic()` with deviation-counting logic using lag columns.
- `carry_forward()` — LOCF intervention: `A_t = A_{t-1}`. Wraps `dynamic()`
  with `rule = \(data, trt) data[["lag1_<trt>"]]`.

These are factory functions returning `dynamic()` interventions. They
require the lag columns materialized by `prepare_data()`.

## Out of scope

- Forward-simulation g-formula (owned by `gfoRmula`)
- Optimal dynamic treatment regimes (owned by `DTRreg`)
