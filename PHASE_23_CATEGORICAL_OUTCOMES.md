# Phase 23 — Multinomial and Ordinal Outcome Models

> **Status: PENDING** (design doc)
>
> **Depends on:** Phase 2 (point gcomp), Phase 13 (extended outcome types)

## Scope

Extend the outcome model layer to support multi-category and ordered
categorical outcomes. Phase 13 expanded outcome *families* within the
GLM framework (Poisson, Gamma, NB, beta). Phase 23 goes beyond GLM to
model classes that return matrix-valued predictions.

### 23a — Multinomial outcomes

Multi-category outcomes via `nnet::multinom()` or
`VGAM::vglm(multinomial())`.

**Key challenges:**
- `predict()` returns an n x K probability matrix, not a vector.
  `causatr_result` currently assumes scalar E[Y^a].
- `contrast()` must be adapted: the estimand becomes a K-vector of
  class probabilities under intervention. Contrasts (difference, ratio)
  apply per-class.
- Sandwich variance requires the multinomial score for outcome models.
  The treatment-model multinomial machinery in `variance_if.R`
  (`prepare_model_if_multinom()`) can be reused as a template.
- `coef()`, `tidy()`, `print()` must handle K-vector estimates.

### 23b — Ordinal outcomes

Ordered categorical outcomes via `MASS::polr()` or `ordinal::clm()`.

**Key challenges:**
- Proportional-odds model returns cumulative probabilities P(Y <= k).
  The natural estimand is P(Y >= k | do(A=a)) for each cutpoint.
- Same vector-valued result structure as multinomial.
- `MASS::polr()` has a specific Hessian structure for sandwich variance.
- `ordinal::clm()` would be a Suggests dependency.

### Design note

23a and 23b share the same architectural requirement: generalizing
`causatr_result` from scalar to vector-valued marginal means. They
should be designed and implemented together to avoid rework.

## Out of scope

- Survival / time-to-event outcomes (owned by `survatr`)
- Zero-inflated models (would need custom predict + variance)
