# Phase 20 — Per-period IPSI under Longitudinal IPW

> **Status: COMPLETE (univariate).**
>
> **Depends on:** Phase 4 (point IPSI via Kennedy 2019 closed form), Phase 10 (longitudinal IPW core).

## Scope

Extend the chunk 10a longitudinal IPW engine to support `ipsi(delta)`
interventions applied at every period. Kennedy's (2019) closed-form
incremental propensity score weight extends naturally to a per-period
product:

$$
W_i = \prod_{t=1}^{T} \frac{\delta \cdot A_{t,i} + (1 - A_{t,i})}{\delta \cdot p_{t,i} + (1 - p_{t,i})},
$$

where $p_{t,i} = P(A_t = 1 \mid \bar A_{t-1}, \bar L_t)$ is the
per-period propensity. The intervention acts on the propensity at
each period rather than substituting a counterfactual treatment value,
exactly as in point IPSI.

## Current state (implemented)

`compute_ipw_contrast_longitudinal()` accepts univariate `ipsi(delta)`:
no new estimation code was required because the point IPSI machinery
(`ipsi_weight_formula()` in `R/ipw_weights.R`, `make_weight_fn()`'s IPSI
branch) was already per-period-ready --
`compute_density_ratio_weights()` routes IPSI to the closed form on
Bernoulli treatments without touching the time index, so
`compute_longitudinal_weights()` multiplies the per-period factors and
`make_weight_fn_longitudinal()` stacks the per-period IPSI sub-closures
for the sandwich.

Two combinations remain rejected with dedicated classes:
multivariate IPSI (`causatr_longitudinal_ipsi_multivariate`) and IPSI +
`stabilize = "marginal"` (`causatr_longitudinal_ipsi_stabilize`).

## Design

### Per-period IPSI weight

At each period $t$, reuse `compute_density_ratio_weights(tm_t, data_t,
ipsi(delta))`: the closed form already returns the Kennedy weight on
the period-$t$ subset. Multiply across periods within an id, identical
to the longitudinal product structure for static / shift / scale_by.

No new density-evaluation path is needed; the IPSI shortcut sidesteps
`evaluate_density()` entirely.

### Stacked sandwich

`make_weight_fn_longitudinal()` already builds per-period sub-closures
via `make_weight_fn()`, which has an IPSI branch. The cross-derivative
`A_{beta, alpha_t}` for IPSI is well-defined (the closed form is
smooth in alpha through the $\delta p_t + (1 - p_t)$ denominator),
and `numDeriv::jacobian` on the stacked closure handles it
transparently.

### Estimand

Sequential IPSI is a stochastic intervention: at each period the
treatment odds are multiplied by $\delta$. The longitudinal target
parameter is

$$
\psi(\delta) = \mathbb{E}\!\left[\,Y \cdot \prod_{t=1}^{T}
  \frac{\delta \cdot A_t + (1 - A_t)}{\delta \cdot p_t + (1 - p_t)}\,\right],
$$

which interpolates between observed treatment ($\delta = 1$ gives unit
weight, recovering the natural-course mean) and the always-treated
counterfactual ($\delta \to \infty$). Same interpretation as point
IPSI extended across periods.

## Items

- [x] Lift the univariate `ipsi()` rejection in `compute_ipw_contrast_longitudinal()`. Multivariate IPSI now rejects with `causatr_longitudinal_ipsi_multivariate`; IPSI + `stabilize = "marginal"` rejects with `causatr_longitudinal_ipsi_stabilize`.
- [x] Verify per-period dispatch: `compute_longitudinal_weights()` calls `compute_density_ratio_weights()` per period, which routes IPSI to the closed form — confirmed, no change needed.
- [x] Verify `make_weight_fn_longitudinal()` builds correct IPSI sub-closures via `make_weight_fn()` per period — confirmed, no change needed.
- [x] Confirm the variance engine handles IPSI without modification (the cross-derivative is well-defined; `numDeriv::jacobian` on the stacked closure produces the right block structure). Sandwich vs bootstrap parity within 15%.
- [x] Truth-based test (`R-long-ipw5`): 2-period binary DGP, IPSI($\delta = 2$) matches the per-period product-of-Kennedy-weights oracle to 1e-6. New oracle `manual_long_ipsi_oracle()` in `helper-ipw-lmtp-oracle.R`. Difference contrast (`R-long-ipw5c`) also oracle-checked.
- [x] Bootstrap parity test (`R-long-ipw5b`).
- [x] Sequential positivity: default Phase 10a behaviour applies — IPSI per-period weights are bounded by construction (denominator $\delta p_t + (1 - p_t) \in [\min(1, \delta), \max(1, \delta)]$); `warn_seq_positivity()` is intervention-agnostic and already tested.
- [x] Updated `FEATURE_COVERAGE_MATRIX.md`, `NEWS.md`, `CLAUDE.md`, and the longitudinal vignette.

## Out of scope

- IPSI on continuous treatment: same out-of-scope as point IPSI (Kennedy's closed form is binary-only).
- IPSI under multivariate longitudinal: combine with Phase 19 if both ship.

## References

Kennedy EH (2019). Nonparametric causal effects based on incremental
propensity score interventions. *JASA* 114:645–656.
