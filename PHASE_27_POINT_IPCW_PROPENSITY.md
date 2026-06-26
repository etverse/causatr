# Phase 27 — Point IPW / AIPW IPCW estimator standardization

> **Status: PLANNED** — surfaced by the 2026-06-25 package-wide IPCW audit (PR #16,
> the longitudinal IPCW sandwich). Run via the `implement-feature` workflow.

## Why this phase exists

Built-in IPCW (`ipcw = TRUE`) currently folds the censoring weights into the
master `weights` vector inside `causat()`, and that vector feeds **every nuisance
fit — including the treatment-density (propensity) models**. The textbook
IPW+IPCW construction estimates the treatment and censoring models *separately*
and multiplies their weights, applying the censoring weight to the
outcome/MSM side only (Hernán & Robins 2025, *Causal Inference: What If*,
Ch. 12.6 & 17: `W = 1/P(A|H) × 1/P(C=0|H)`). causatr's folding makes the
propensity IPCW-weighted on the uncensored rows — a non-standard estimator.

The **longitudinal** IPW path was fixed in PR #16 (a `propensity_weights`
argument keeps the propensity unweighted; the sandwich gained the γ→β cross-term).
This phase finishes the job for the **point** estimators.

## The confirmed bug (delicatessen-verified)

**Point IPW + IPCW analytic sandwich SE is wrong by 3–7% per arm.** Because the
propensity is IPCW-weighted it depends on γ, so the point sandwich
(`compute_ipw_ipcw_correction()` carries γ→β only) omits the **γ→α** cross-term.

- delicatessen exact own-estimator sandwich: `se_ate = 0.04351`
- causatr reports: `se_ate = 0.04519`  (per-arm `se_mu1` ratio 0.965, `se_mu0` 0.927)
- The **point estimate is consistent** (matches delicatessen to ~1e-12) — this is
  a variance bug, not bias. The bootstrap (ratio ~1.01 at 800 reps) masked it
  under Monte-Carlo noise; the delicatessen M-estimation oracle exposed it.

**Point AIPW + IPCW**: propensity also IPCW-weighted, but doubly-robust
orthogonality makes the γ→α effect negligible (sandwich/bootstrap 1.00–1.02) —
the SE is acceptable; only the estimator *convention* is non-standard.

**SNM longitudinal + IPCW (chunk 18g)**: runs and produces a censoring-block
sandwich; not independently validated against an oracle in the audit.

## Scope

In scope:
1. **Point IPW + IPCW** — standardize the propensity to **unweighted on all
   rows** (the textbook estimator, matching the standardized longitudinal path),
   and make the sandwich exactly correct.
2. **Point AIPW + IPCW** — same propensity standardization (for consistency; its
   SE is already orthogonal-correct, so this is convention alignment + a guard).
3. **SNM longitudinal + IPCW** — validate the censoring-block sandwich against the
   bootstrap / an oracle, or reject the combination with a classed error.

Out of scope: the longitudinal IPW/AIPW/ICE IPCW sandwiches (done / valid — see
`PHASE_14_IPCW.md` "Post-shipment audit"); changing the IPCW *weight* semantics
(per-period vs cumulative — a separate established design choice).

## Key design decision — why this is deeper than the IF

Standardizing the point estimator is **not** an IF-only change. `make_weight_fn()`
(`R/ipw_weights.R`) and the point contrast assembly are hard-coupled to the
propensity being fit on the **same rows** as the MSM — they read
`treatment_model$fit_rows` and `treatment_model$X_prop`. Fitting the propensity
on *all* rows (the standard estimator) breaks that coupling in the point-estimate
path, not just the variance path.

**Chosen approach (user-approved):** gate the new behaviour on `ipcw` so the
heavily-validated **non-IPCW point IPW path stays byte-identical**, and build a
parallel point-IPW-IPCW estimation + variance path:

- Fit the propensity **unweighted on all rows** (external weights only, via the
  existing `propensity_weights` thread).
- Compute the per-intervention Hájek marginal means directly (propensity-on-all-
  rows density ratio × IPCW weight on the uncensored MSM rows).
- Variance via a **stacked-EE numerical-bread sandwich** (θ = α + γ + per-arm μ),
  `V = B⁻¹ M B⁻ᵀ / n` with `B` the `numDeriv` Jacobian of the summed score —
  the same construction `variance_if_aipw_longitudinal_ee.R` uses, which matches
  `delicatessen.MEstimator` by construction.

With the propensity IPCW-unweighted, γ reaches μ only through the MSM (γ→β), so
no γ→α term is needed — the bug disappears with the standardization.

## Implementation chunks

| Chunk | Scope | Depends on |
|---|---|---|
| 27a | `causat()` / `fit_ipw()` point path: under `ipcw`, fit the propensity unweighted on all rows (reuse `propensity_weights`); keep the MSM on uncensored rows | PR #16 infra |
| 27b | Gated stacked-EE sandwich for point IPW + IPCW (α all-rows + γ + Hájek μ); route `variance_if()` to it when `ipcw && estimator == "ipw" && type == "point"` | 27a |
| 27c | Point AIPW + IPCW: same propensity standardization; verify SE stays orthogonal-correct | 27a |
| 27d | SNM longitudinal + IPCW: validate the 18g censoring-block sandwich vs bootstrap/oracle, or reject | — |
| 27e | Tighten the loose NHEFS IPCW test (`test-delicatessen-nhefs.R`, tol 0.1 → tight vs delicatessen); add a point-IPW-IPCW delicatessen fixture | 27a–27c |

## Oracle strategy

- **delicatessen** M-estimation stack for point IPW + IPCW (propensity unweighted
  on all rows + censoring γ + IPCW-weighted Hájek μ) — the primary cross-language
  oracle. The audit already prototyped it (`se_ate = 0.04351` on the test DGP);
  promote it to a committed fixture under `tests/testthat/fixtures/python/`.
- Monte-Carlo SE-vs-empirical-SD calibration; bootstrap parity.
- Tighten the existing NHEFS IPCW oracle test once the estimator matches.

## References

- Hernán MA, Robins JM (2025). *Causal Inference: What If*. Ch. 12.6, 17.
- Robins JM, Rotnitzky A, Zhao LP (1994). Estimation of regression coefficients
  when some regressors are not always observed. *JASA* 89:846–866.
- Stefanski LA, Boos DD (2002). The calculus of M-estimation. *Am Stat* 56:29–38.
