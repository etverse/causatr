# Phase 9 — Inference Infrastructure

> **Status: PARTIAL** (parallel bootstrap done; survey weights, clustered SE, `future` backend pending)

## Scope

Variance and inference enhancements that harden the existing estimators. Survey weight integration, general cluster-robust sandwich, and optional `future` backend.

## Items

- [ ] Survey weights: full integration with `survey::svydesign()` (pass through to model fitting + adjust sandwich)
- [ ] Clustered data: cluster-robust sandwich via `sandwich::vcovCL()` (beyond matching subclass — general cluster variable)
- [x] Parallel processing: `boot::boot(parallel=, ncpus=)` for all methods (gcomp, IPW, matching, ICE). `contrast()` accepts `parallel` and `ncpus` params.
- [ ] Parallel processing (optional): `future` backend as alternative to `boot::boot()` built-in parallelism

## Dependencies

Phases 1–5 only. Fully independent of Phases 6–8 and 10. Can run in parallel with any other phase.
