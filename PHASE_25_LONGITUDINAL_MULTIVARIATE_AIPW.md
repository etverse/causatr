# Phase 25 — Multivariate Longitudinal AIPW

> **Status: COMPLETE — 25a, 25b, 25c all shipped.**
>
> **Depends on:** Phase 16i (univariate longitudinal AIPW), Phase 16m
> (multivariate point AIPW), Phase 19 (multivariate longitudinal IPW).

## Scope

Combine the doubly-robust ICE-AIPW machinery (Phase 16i) with the
multivariate per-period density chain (Phase 19) to support **joint
time-varying treatments under AIPW**: `treatment = c("A1", "A2")` combined
with `id =`, `time =`, and `estimator = "aipw"`. This fills the last
unfilled cell in the "multivariate × longitudinal" matrix — MV gcomp/ICE
and MV longitudinal IPW already shipped; this phase adds the doubly-robust
analogue.

The estimand is the joint time-varying intervention mean under sequential
MTP semantics (Díaz et al. 2023). For a **static** contrast, MV longitudinal
IPW = MV ICE g-computation = MV longitudinal AIPW coincide, so the binary
static analytical truths from `make_em_mv_long_scm()` (9 at sex = 0, 15 at
sex = 1) are valid AIPW oracles.

## Why this is a small, low-risk composition

Both halves already support multivariate treatment independently:

- **Outcome (ICE) side already loops over vector `treatment`.**
  `ice_build_formula()` and `ice_apply_intervention_long()` iterate over the
  treatment columns. The Phase 19c cross-method test confirmed MV
  longitudinal ICE g-computation runs correctly.
- **Propensity (IPW) side has a ready-made MV chain.** Phase 19 built
  `fit_treatment_models()` (per-period K-component chain),
  `compute_density_ratio_weights_mv()`, and `make_weight_fn_mv()`.

So the work is swapping the univariate propensity calls for their existing
multivariate counterparts and stacking $T \times K$ (not $T$) propensity
blocks in the sandwich. No new estimation math.

## Design

### ICE-AIPW recursion (unchanged)

The Bang & Robins (2005) backward recursion is

$$
\tilde Q_k = m_k(d_k) + W_k \cdot \bigl(\tilde Q_{k+1} - m_k(A_k)\bigr),
$$

where $m_k$ is the period-$k$ sequential outcome model, $d_k$ the
intervened treatment, and $W_k$ the **single-period** density-ratio weight.
The recursion uses the single-period weight; only the per-period weight
computation changes for the multivariate case.

### Multivariate per-period weight

For `length(treatment) > 1L`, each period's weight is the product of the
$K$ within-period component ratios via `compute_density_ratio_weights_mv()`:

$$
W_{t,k=\text{period}} = \prod_{k=1}^{K}
\frac{f_{t,k}\bigl(d_k(A_{t,k}) \mid A_{t,1..k-1}^{\mathrm{obs}}, \bar A_{t-1}^{\mathrm{obs}}, \bar L_t\bigr) \cdot \lvert \mathrm{Jac}\,d_k^{-1} \rvert}{f_{t,k}\bigl(A_{t,k} \mid A_{t,1..k-1}^{\mathrm{obs}}, \bar A_{t-1}^{\mathrm{obs}}, \bar L_t\bigr)}.
$$

The single per-period scalar this returns feeds the recursion exactly as the
univariate single-period weight did, so the backward iteration downstream is
unchanged.

### Per-period propensity chain

In `fit_aipw_longitudinal()` the per-period loop dispatches on
`length(treatment)`: the multivariate path fits the per-period **chain** with
`fit_treatment_models()` (returns a `causatr_treatment_models` list of $K$
component models) instead of the single `fit_treatment_model()`. Each period
stores a $K$-model list in `treatment_models_by_time`. EM stripping and
baseline-term wiring happen inside `fit_treatment_models()`, exactly as in the
IPW MV path.

### Stacked sandwich (Channel 2b, $T \times K$)

The longitudinal AIPW IF has three channels in
`variance_if_aipw_long_one()`:

1. **Channel 1** — direct (the pseudo-outcome IF), unchanged.
2. **Channel 2a** — outcome-model cascade, unchanged (the ICE outcome
   machinery already handles vector treatment).
3. **Channel 2b** — propensity correction. For the multivariate case the
   per-period sub-closures are built with `make_weight_fn_mv()` and the
   block-diagonal bread stacks $T \times K$ α-blocks (the within-period
   component parameters concatenated) instead of $T$. This mirrors how
   `make_weight_fn_longitudinal()` delegates to `make_weight_fn_mv()` in the
   IPW MV path. The `aug_mean(alpha)` numDeriv closure reconstructs
   per-period weights from perturbed α via the MV sub-closure; the rest of
   the block-diagonal assembly is unchanged.

The propensity blocks are fit on **disjoint** time subsets (block-diagonal
across time) and are chain-rule-independent within each period
(block-diagonal across components), so the full $\sum_{t,k} p_{t,k}$ bread is
block-diagonal and the engine sums $T \times K$ single-model
`apply_model_correction()` calls.

### Bootstrap (unchanged)

`aipw_longitudinal_variance_bootstrap()` refits via `fit_aipw_longitudinal()`
and re-iterates `ice_aipw_iterate()`. Once the point engine dispatches on
`length(treatment)`, id-clustered resampling already preserves both axes
(same as the Phase 19 IPW bootstrap). Verified by `T-long-mv-aipw2/5`.

## Items

### 25a — core engine ✅

- [x] Remove the `causatr_longitudinal_multivariate_pending` rejection in
  `fit_aipw_longitudinal()`.
- [x] Dispatch the per-period loop on `length(treatment)`: multivariate path
  fits `fit_treatment_models()` (K-component chain) per period; univariate
  path unchanged. NULL-safe baseline-term derivation and the multivariate
  `prop_model_fn` fallback mirror `fit_longitudinal_ipw()`.
- [x] `ice_aipw_iterate()` forward per-period weight precompute calls
  `compute_density_ratio_weights_mv()` for the multivariate case (reset
  `fit_rows` + class on a local copy). Backward outcome iteration unchanged.
- [x] Channel 2b in `variance_if_aipw_long_one()` stacks $T \times K$
  α-blocks via `make_weight_fn_mv()`; the correction loop slices each
  period's stacked Jacobian into per-component blocks and applies
  `apply_model_correction()` per component.

### 25b — stabilization ✅ (rejected, not implemented)

- [x] Stabilization (`stabilize = "marginal"`) does **not** compose for free
  for longitudinal AIPW: `fit_aipw()` silently dropped `stabilize` before
  reaching `fit_aipw_longitudinal()`, so a stabilized request would have been
  honoured nowhere. Rather than implement a per-period × per-component
  numerator sweep into the AIPW recursion, this combination is **explicitly
  rejected** with a classed error (`causatr_stabilize_longitudinal_aipw`),
  pointing the user to `estimator = "ipw"` for stabilized longitudinal
  weights. The rejection fires for both univariate and multivariate
  longitudinal AIPW.

### 25c — docs + vignettes + dependency cross-refs ✅

- [x] This design doc.
- [x] `FEATURE_COVERAGE_MATRIX.md` — chunk-16m deferral updated to "supported
  (Phase 25)"; MV longitudinal AIPW block (6-row table) added; rejection
  notes for `stabilize` and `ipsi()`.
- [x] `PHASE_19` downstream note marked addressed.
- [x] `PHASE_16_AIPW.md`, `PHASE_17_TRANSPORTABILITY.md` downstream notes.
- [x] `CLAUDE.md` phase index, supported-features row, stabilize
  design-decision bullet.
- [x] `NEWS.md` dated entry.
- [x] Vignettes (`longitudinal.qmd`, `aipw.qmd`) combination tables + worked
  code chunk.

## Tests

All in `tests/testthat/test-aipw.R` (colocated with the existing longitudinal
AIPW tests).

- `T-long-mv-aipw1` — binary static, sandwich: recovers the stratified
  analytical truths 9 (sex = 0) / 15 (sex = 1) via `by = "sex"` on
  `make_em_mv_long_scm()`.
- `T-long-mv-aipw2` — sandwich vs bootstrap SE parity (ratio in 0.7–1.4).
- `T-long-mv-aipw3` — cross-method: MV AIPW ≈ MV ICE g-comp ≈ MV longitudinal
  IPW on the static binary DGP (per-stratum abs diff < 0.5).
- `T-long-mv-aipw4` — double-robustness: (a) wrong outcome model + correct
  propensity recovers truth; (b) correct outcome + wrong propensity recovers
  truth.
- `T-long-mv-aipw5` — continuous shift (inline 2-period × 2-component DGP),
  sandwich + bootstrap finite and SE parity (no closed-form truth under
  shift; parity + finiteness only).
- `R-long-mv-aipw1` — `stabilize = "marginal"` rejected for longitudinal AIPW
  (class `causatr_stabilize_longitudinal_aipw`).

**Oracle note:** longitudinal AIPW has no external single-call oracle in this
suite (lmtp's SDR differs in the nuisance-fixing convention). The established
pattern is cross-method triangulation (`T-long-mv-aipw3`) + double-robustness
(`T-long-mv-aipw4`) + truth recovery on a static binary DGP
(`T-long-mv-aipw1`), documented in the test-file header comment.

## Out of scope

- **Transport × longitudinal AIPW** (univariate + multivariate) — owned by
  **Phase 26** (`PHASE_26_LONGITUDINAL_AIPW_TRANSPORT.md`, PENDING design
  doc). `variance_if_aipw_long_one()` has no sampling-model correction
  channel and there is no longitudinal AIPW sampling-correction helper;
  Phase 26 designs the bootstrap-only path.
- **Longitudinal IPSI under AIPW** — stays rejected
  (`causatr_longitudinal_ipsi_pending`, inherited from the univariate path;
  Phase 20 territory).
- **Stabilization** — explicitly rejected (see 25b); use `estimator = "ipw"`.

## References

Bang H, Robins JM (2005). Doubly robust estimation in missing data and
causal inference models. *Biometrics* 61(4):962–973.

Díaz I, Williams N, Hoffman KL, Schenck EJ (2023). Non-parametric causal
effects based on longitudinal modified treatment policies. *JASA*
118:846–857.
