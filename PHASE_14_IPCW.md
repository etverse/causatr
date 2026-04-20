# Phase 14 — Built-in IPCW for MAR outcome censoring

## Motivation

causatr's `censoring =` parameter today is a **row filter** only. Phase 14 adds a **built-in IPCW** path so users do not have to compute censoring weights externally for MAR outcome censoring.

Survival-specific IPCW (per-period cumulative censoring weights over a person-period grid, weighted pooled-logistic hazard, survival-curve variance) is owned by the separate survival package (see `SURVIVAL_PACKAGE_HANDOFF.md`). Phase 14 in causatr targets **scalar outcomes** only — point outcomes with MAR missingness and longitudinal ICE with a final scalar outcome. The survival package imports Phase 14 for the censoring-model fitting primitive and layers cross-time aggregation on top.

### How `censoring =` works today

causatr's `censoring =` parameter currently performs **row filtering only**:
`get_fit_rows()` drops rows where `C != 0` before model fitting, and
`contrast()` predicts on all rows regardless of censoring. This is the
approach used by Zivich et al. (2024) in the ICE g-computation paper —
censoring enters the estimating equations as indicator functions
`I(C_{i,k} = 0)` that restrict each model to uncensored individuals.

This row-restriction approach is **correct under censoring
exchangeability**:

> $Y_k \perp C_k \mid \bar{A}_{k-1}, \bar{L}_{k-1}, C_{k-1} = 0$

When this holds and the outcome model E[Y | A, L] is correctly
specified, g-computation on the uncensored sample is consistent — the
regression surface is the same for censored and uncensored individuals.
No IPCW weights are needed.

### Why IPCW adds value beyond the row filter

1. **Doubly robust protection.** If the outcome model is misspecified,
   IPCW weights reduce bias. IPCW + g-comp gives a "doubly robust"
   estimator (Bang & Robins, 2005).

2. **IPW estimator.** The IPW path does not fit an outcome model on
   complete cases — it relies entirely on treatment weights. Under MAR
   censoring, the IPW weights must be combined with IPCW weights (the
   "doubly weighted" MSM; Hernán & Robins Ch. 12).

3. **Efficiency.** Under MAR censoring, IPCW-weighted estimation can be
   more efficient than unweighted complete-case analysis because the
   weights stabilize the contribution of each observation.

## What Zivich et al. (2024) say about censoring

The paper (Statistics in Medicine 43:5562–5572) handles censoring through
**indicator functions** restricting estimating equations to uncensored rows.
No IPCW weights are used. Their simulation induces informative
monotonic dropout via logistic models dependent on treatment history:

```
C_{i,1} ~ Bernoulli(expit(-2.5 - 0.5 * A_{i,0}))
C_{i,2} ~ 1 if C_{i,1} = 1; else Bernoulli(expit(-2.5 - 0.5 * A_{i,1}))
```

Identification requires:
- **Censoring exchangeability**: loss-to-follow-up at time k is
  non-informative conditional on L̅_k, A̅_k, Y̅_k
- **Censoring positivity**: P(C_k = 0 | a̅*_{k-1}, l̅_{k-1}, C_{k-1} = 0) > 0

The paper does **not** discuss combining ICE with IPCW weights.
causatr already supports external IPCW weights via `weights =`
(validated and propagated through every ICE backward step), so the
doubly-robust combination is available today as a manual step.

## Current state: what `censoring =` does today

```
                    censoring = "C"
                         │
                         ▼
              ┌─────────────────────┐
              │   get_fit_rows()    │
              │  is_uncensored(C)   │
              │  & !is.na(Y)        │
              └────────┬────────────┘
                       │
              fit_rows = TRUE for rows
              where C == 0 and Y is observed
                       │
          ┌────────────┼────────────┐
          ▼            ▼            ▼
      fit_gcomp    fit_ice      fit_ipw / fit_matching
      (point)      (longitudinal)  (censoring = NULL;
       │            │               ignored at fit time)
       │            │
       │            ▼
       │        Per-step backward:
       │        fit on uncensored rows,
       │        predict on all rows
       ▼
    Outcome model fit on
    uncensored + non-NA-Y rows.
    contrast() predicts on ALL rows.
```

**Key fact**: No censoring model is fit. No IPCW weights are computed.
The `censoring =` parameter is a row filter, not IPCW.

## Proposed: built-in IPCW

### User API

```r
# Option A: built-in censoring model (new)
fit <- causat(
  data,
  outcome = "Y", treatment = "A", confounders = ~ L1 + L2,
  censoring = "C",
  ipcw = TRUE,                    # NEW: fit censoring model internally
  censoring_model_fn = stats::glm # NEW: fitter for P(C=0 | A, L)
)

# Option B: manual IPCW (works today)
# User computes w_i = 1/P(C=0 | A_i, L_i) externally
fit <- causat(data, ..., weights = ipcw_weights)
```

### Censoring model

For **point treatments**:
- Fit: P(C = 0 | A, L) using `censoring_model_fn` (default: `stats::glm`)
- Formula: `C ~ A + <confounders>` (same RHS as the outcome model, but
  with the censoring indicator as response)
- Weights: stabilized IPCW = P(C = 0) / P(C = 0 | A_i, L_i)

For **longitudinal treatments** (ICE):
- At each time step k, fit: $P(C_k = 0 \mid \bar{A}_k, \bar{L}_k, C_{k-1} = 0)$
- Cumulative weight: $w^C_i = \prod_k 1/P(C_k = 0 \mid \ldots)$
- Stabilized: $w^C_i = \prod_k P(C_k = 0 \mid \bar{A}_{k-1}, \bar{L}_{k-1}) / P(C_k = 0 \mid \bar{A}_k, \bar{L}_k)$
- These weights enter multiplicatively with any existing external weights
  and are propagated through every ICE backward step (architecture
  already supports this via `weights =`)

### Sandwich variance extension

When `ipcw = TRUE`, the censoring model parameters enter the stacked
estimating equation system. The IF becomes:

```
IF_i = IF_outcome_i + IF_censoring_correction_i
```

where the censoring correction accounts for the fact that the IPCW
weights are estimated (not known). This follows the standard
M-estimation result for estimated weights (Stefanski & Boos, 2002).

For ICE, the stacked system grows from K+1 outcome models to
K+1 outcome models + K censoring models (one per time step). The
cascade gradient must be extended to include the censoring model
parameters at each step.

### Integration with estimators (scalar outcomes)

| Estimator | IPCW role | Implementation |
|---|---|---|
| **gcomp (point)** | IPCW-weighted outcome model: `glm(Y ~ A + L, weights = ipcw)` | Multiply IPCW into `model_weights` inside `fit_gcomp_point()` |
| **gcomp (ICE, scalar final outcome)** | IPCW-weighted pseudo-outcome models at each step | Multiply IPCW into step-level weights in `ice_iterate()` |
| **IPW** | IPCW × IPW "doubly weighted" MSM: `glm(Y ~ 1, weights = ipw * ipcw)` | Multiply IPCW into `w_final` inside `compute_ipw_contrast_point()` |
| **matching** | IPCW-weighted outcome model on matched data | Multiply IPCW into `matched_weights` |

### Diagnostics (Phase 11 integration)

`diagnose()` should report:
- Censoring model balance (SMD of covariates weighted by IPCW)
- Weight distribution summary (mean, sd, min, max, ESS)
- Positivity warnings (extreme IPCW weights)

## Missing data: complete scenario map

```
┌──────────────────────────────────────────────────────────────────┐
│                    MISSING DATA IN causatr                       │
├─────────────┬────────────────┬───────────────────────────────────┤
│  Variable   │  Mechanism     │  Handling                         │
├─────────────┼────────────────┼───────────────────────────────────┤
│             │  MCAR          │  ✅ Complete-case: drop via       │
│  Y          │                │     get_fit_rows(). Unbiased.     │
│  (outcome)  ├────────────────┼───────────────────────────────────┤
│             │  MAR on (A,L)  │  ✅ G-comp: correct if outcome    │
│             │                │     model correctly specified.    │
│             │                │  ❌ IPW: needs IPCW weights.      │
│             │                │  🔧 Manual: weights= today.      │
│             │                │  ❌ Built-in: Phase 14.           │
│             ├────────────────┼───────────────────────────────────┤
│             │  MNAR          │  ⛔ Out of scope (sensitivity     │
│             │                │     analysis).                    │
├─────────────┼────────────────┼───────────────────────────────────┤
│             │  MCAR          │  ✅ User drops rows before        │
│  A          │                │     causat(). Unbiased.           │
│  (treatment)├────────────────┼───────────────────────────────────┤
│             │  MAR on L      │  ❌ MI: impute A from P(A|L)      │
│             │                │     via causat_mice(). Planned.   │
│             ├────────────────┼───────────────────────────────────┤
│             │  MNAR          │  ⛔ Out of scope.                 │
├─────────────┼────────────────┼───────────────────────────────────┤
│             │  MCAR          │  ✅ G-comp: glm na.omit.          │
│  L          │                │     IPW: abort if masks diverge.  │
│  (covariates)                │     Matching: MatchIt handles.    │
│             ├────────────────┼───────────────────────────────────┤
│             │  MAR           │  ❌ MI: impute L from P(L|A,Y,L') │
│             │                │     via causat_mice(). Planned.   │
│             ├────────────────┼───────────────────────────────────┤
│             │  MNAR          │  ⛔ Out of scope.                 │
└─────────────┴────────────────┴───────────────────────────────────┘
```

## Why IPCW conditions on (A, L), not just L

The censoring exchangeability assumption is:

> $Y^a \perp C \mid A, L$

This says potential outcomes are independent of censoring **given both
treatment and covariates**. The censoring model must therefore condition
on A and L:

P(C = 0 | A, L)

Why A matters: censoring can depend on treatment received (e.g., side
effects cause dropout in the treated arm). If we only condition on L,
we'd miss this source of informative censoring and the IPCW weights
would be wrong.

For longitudinal data, the sequential version is:

P(C_k = 0 | A̅_k, L̅_k, C_{k-1} = 0)

## Implementation chunks

| Chunk | Scope | Depends on |
|---|---|---|
| 14a | `fit_censoring_model()`: fit $P(C = 0 \mid A, L)$ (point) and $P(C_k = 0 \mid \bar{A}_k, \bar{L}_k, C_{k-1} = 0)$ (longitudinal) using `censoring_model_fn`; return `causatr_censoring_model` with per-row stabilized weights | — |
| 14b | **Point scalar outcome + gcomp/IPW/matching:** IPCW-weighted fit, sandwich extension | 14a |
| 14c | **ICE scalar final outcome:** IPCW in the ICE backward loop; cascade-gradient extension for censoring model blocks | 14a, 14b |
| 14d | lmtp cross-checks: truth-based tests against `lmtp::lmtp_tmle(cens = "C")` on point and longitudinal DGPs with informative censoring | 14b, 14c |
| 14e | Sandwich extension: formalise the stacked EE system, document the block structure, regression tests | 14a–14c |
| 14f | `diagnose()` integration: censoring balance, weight distribution summary, positivity warnings | 14b, Phase 11 |
| 14g | Vignette: missing-data handling guide with built-in IPCW examples | 14b, 14c |

## References

- Bang H, Robins JM (2005). Doubly robust estimation in missing data
  and causal inference models. *Biometrics* 61:962–973.
- Cole SR, Hernán MA (2004). Adjusted survival curves with inverse
  probability weights. *Comput Methods Programs Biomed* 75:45–49.
- Hernán MA, Robins JM (2025). *Causal Inference: What If*. Ch. 12, 17.
- Howe CJ, Cole SR, Lau B, Napravnik S, Eron JJ (2016). Selection bias
  due to loss to follow up in cohort studies. *Epidemiology* 27:91–97.
- Stefanski LA, Boos DD (2002). The calculus of M-estimation. *Am Stat*
  56:29–38.
- Zivich PN, Ross RK, Shook-Sa BE, Cole SR, Edwards JK (2024).
  Empirical sandwich variance estimator for iterated conditional
  expectation g-computation. *Stat Med* 43:5562–5572.
