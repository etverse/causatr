# Phase 26 — Transport × Longitudinal AIPW

> **Status: PENDING** (design doc)
>
> **Depends on:** Phase 16i (univariate longitudinal AIPW — done), Phase 17
> (transportability: chunks 17e/17i/17j/17l — done), Phase 25 (multivariate
> longitudinal AIPW — done).

## Motivation

Phase 17 built transport/generalizability (`target =`) across the point
estimators (gcomp, IPW, AIPW) and the **longitudinal IPW** path (chunk 17i),
plus the MTP + transport MC-marginalization machinery (chunk 17l). It
deliberately left **one cell open**: transport × **longitudinal AIPW**. This
is the single remaining unfilled combination in the transport matrix.

The gap is concrete. The point AIPW IF handles transport via the
`is_transport` branch in `variance_if_aipw.R`, but
`variance_if_aipw_long_one()` has **no sampling-model correction channel**,
and there is no longitudinal-AIPW analogue of the IPW sampling-correction
helper (`variance_if_ipw_sampling.R`). So a longitudinal AIPW fit with
`target =` cannot currently incorporate sampling-odds reweighting.

This phase designs that composition for **both** univariate and multivariate
longitudinal AIPW (the latter composing with Phase 25's per-period density
chain).

## Scope

1. **Transportability / generalizability** (`target =`) for the longitudinal
   AIPW path — univariate and multivariate.
2. Both `target_subset` modes: `"target"` (transportability, S = 0 rows) and
   `"all"` (generalizability), matching the rest of Phase 17.
3. **Bootstrap-only variance** — the sandwich is rejected for this
   combination (see "Variance" below), following the chunk-17l precedent for
   transport AIPW cases needing MC marginalization.

## Non-scope

- **Sandwich variance.** Rejected by a classed error (see below). The MC
  marginalization over $P(A \mid L, S = 1)$ on target rows makes the
  influence function non-differentiable in closed form, the same reason
  chunk 17l rejects the sandwich for MTP + transport.
- **Longitudinal IPSI under AIPW** — inherits the
  `causatr_longitudinal_ipsi_pending` rejection from the Phase 25 path.
- **Stabilization** — inherits the `causatr_stabilize_longitudinal_aipw`
  rejection from Phase 25.

## Design

### Point estimate — sampling-odds broadcast + ICE-AIPW recursion

Reuse the chunk-17i broadcast: fit the first-period sampling model
$P(S = 1 \mid L)$ on combined study + target baseline rows, form the
sampling-odds weight

$$
\omega_i = \frac{1 - \hat p(L_i)}{\hat p(L_i)}, \qquad
\hat p(L_i) = \widehat{P}(S = 1 \mid L_i),
$$

and **broadcast** $\omega_i$ from the baseline row to all person-period rows
of individual $i$ (sampling status is time-invariant). Fold $\omega_i$ into
the ICE-AIPW backward recursion's target average — i.e. the final
standardization over the target population is the $\omega$-weighted mean of
the period-1 pseudo-outcome.

The wrinkle: under `target_subset = "target"`, the S = 0 (target) rows lack
observed treatment **and** outcome, so the sequential outcome models cannot
produce in-sample predictions for them. The intervened outcome predictions
for target rows therefore need the **MC marginalization** machinery already
built for chunk 17l — `mc_marginalize_preds()` /
`mc_marginalize_binary()` / `mc_marginalize_continuous()` in `contrast.R` —
to integrate $\hat m_k(d_k(A), L)$ over $P(A \mid L, S = 1)$ at each period.

### Multivariate composition (Phase 25)

For `length(treatment) > 1L`, the per-period propensity chain is Phase 25's
`fit_treatment_models()` and the per-period density-ratio weight is
`compute_density_ratio_weights_mv()`. The sampling-odds broadcast is
orthogonal to the treatment-component axis (sampling is per-individual, not
per-component), so the multivariate composition is: MV per-period density
chain × single broadcast sampling weight, folded into the same
$\omega$-weighted target average. No new MV-specific transport math.

### Variance — bootstrap only

Reject the sandwich with a classed error (either extend
`causatr_mtp_transport_sandwich` or add a dedicated
`causatr_longitudinal_aipw_transport_sandwich`). Support **id-clustered
bootstrap** that refits the sampling model per replicate: the
`aipw_longitudinal_variance_bootstrap()` path already refits via
`fit_aipw_longitudinal()`, so the extension is to also refit the sampling
model on each resampled study + target dataset, exactly as the longitudinal
IPW transport bootstrap (17i) does. Replicates with a degenerate S
distribution return NA and are excluded (Phase 2 convention).

## Reuse map

| Need | Existing machinery | Source |
|---|---|---|
| Sampling model $P(S=1\mid L)$ | `fit_sampling_model()` | `R/sampling_model.R` |
| Sampling-odds broadcast onto person-periods | chunk 17i broadcast | `R/variance_if_ipw_sampling.R` + longitudinal IPW path |
| Target-row outcome predictions (no observed A) | `mc_marginalize_preds()` / `_binary()` / `_continuous()` | `R/contrast.R` (chunk 17l) |
| MV per-period propensity chain | `fit_treatment_models()`, `compute_density_ratio_weights_mv()` | Phase 25 / Phase 19 |
| ICE-AIPW backward recursion | `ice_aipw_iterate()` | `R/aipw_longitudinal.R` |
| Bootstrap refit (extend to refit sampling model) | `aipw_longitudinal_variance_bootstrap()` | `R/variance_bootstrap_longitudinal.R` |

## Chunk plan

1. **26a — univariate longitudinal AIPW × transport**: sampling-odds
   broadcast into `ice_aipw_iterate()`; MC marginalization for target-row
   outcome predictions; sandwich-rejection classed error; id-clustered
   bootstrap refitting the sampling model. Truth-based transportability /
   generalizability ATE on a study + target longitudinal DGP.
2. **26b — multivariate longitudinal AIPW × transport**: compose 26a with the
   Phase 25 per-period density chain. Truth-based MV transport ATE.
3. **26c — docs / vignettes**: `transportability.qmd` + `longitudinal.qmd`
   combination rows, `FEATURE_COVERAGE_MATRIX.md`, `CLAUDE.md`, NEWS.

## Tests (planned)

- Truth-based transportability (S = 0 target) and generalizability (all)
  ATE on a study + target longitudinal DGP — univariate and multivariate.
- Bias-correction vs naive (study-only) estimate: transport estimate should
  recover the target-population truth where the naive estimate is biased.
- Cross-method: longitudinal AIPW transport ≈ longitudinal IPW transport
  (17i) ≈ ICE g-computation transport on a static contrast (where the
  estimands coincide).
- Sandwich-rejection classed-error test.
- Bootstrap finiteness + SE sanity.

**Oracle note:** no external single-call oracle covers longitudinal AIPW
transport. Validation is triangulation against longitudinal IPW transport
(17i) and ICE g-comp transport cross-checks, plus truth recovery on a
constructed study + target DGP — the same pattern Phase 25 uses for
longitudinal AIPW.

## References

Bang H, Robins JM (2005). Doubly robust estimation in missing data and
causal inference models. *Biometrics* 61(4):962–973.

Díaz I, Williams N, Hoffman KL, Schenck EJ (2023). Non-parametric causal
effects based on longitudinal modified treatment policies. *JASA*
118:846–857.

Dahabreh IJ, Robertson SE, Steingrimsson JA, Stuart EA, Hernán MA (2020).
Extending inferences from a randomized trial to a new target population.
*Statistics in Medicine* 39(14):1999–2014.
