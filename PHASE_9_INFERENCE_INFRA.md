# Phase 9 — Inference Infrastructure

> **Status: DONE** (parallel bootstrap, general cluster-robust sandwich, survey `svydesign` integration, and `future` backend all shipped)

## Scope

Variance and inference enhancements that harden the existing estimators. Survey weight integration, general cluster-robust sandwich, and optional `future` backend.

## Items

- [x] Survey weights: full integration with `survey::svydesign()` — `causat(weights = svydesign_obj)` extracts sampling weights via `stats::weights()` and auto-adopts the first-stage PSU as the contrast-time cluster. Explicit `cluster =` overrides. Single-PSU designs (`~1`) propagate only weights. See `test-survey-design.R`.
- [x] Clustered data: cluster-robust sandwich via `vcov_from_if(cluster = ...)` (Liang & Zeger 1986 sum-within-cluster-then-square), equivalent to `sandwich::vcovCL()` on the predict-then-average step. `contrast(cluster = ...)` and `causat(cluster = ...)` accept a column name or vector. Supported for gcomp, IPW, ICE; matching rejects (subclass cluster is structural). See `test-cluster-sandwich.R`.
- [x] Parallel processing: `boot::boot(parallel=, ncpus=)` for all methods (gcomp, IPW, matching, ICE). `contrast()` accepts `parallel` and `ncpus` params.
- [x] Parallel processing (optional): `future` backend — `contrast(parallel = "future")` dispatches bootstrap replicates through `future.apply::future_lapply()` so any active `future::plan()` is honoured. Reproducible under `set.seed()` via `future.seed = TRUE`. See `test-future-backend.R`.

## Chunk plan

- **9a (done)** — General cluster-robust sandwich. `cluster` parameter on `causat()` (preserves the column through `prepare_data()`) and `contrast()` (threads into `variance_if()` → `vcov_from_if(cluster = ...)`). Oracle: `sandwich::vcovCL(fit$model, cluster = ..., type = "HC0")` on a saturated outcome model. Matching rejects user-provided cluster via `causatr_bad_cluster` at fit and contrast time.
- **9b (done)** — `survey::svydesign` object as `weights =`: extract weights via `stats::weights()` and auto-propagate PSU into the fit's cluster slot via `unpack_svydesign()` in `R/causat.R`. User `cluster =` override wins.
- **9c (done)** — `future` backend: `parallel = "future"` on `contrast()` dispatches bootstrap replicates through `future.apply::future_lapply()` via `dispatch_boot()` in `R/variance_bootstrap.R`. Honours any active `future::plan()`; `future.seed = TRUE` keeps the seed-aware contract of `boot::boot()`.

## Dependencies

Phases 1–5 only. Fully independent of Phases 6–8 and 10. Can run in parallel with any other phase.
