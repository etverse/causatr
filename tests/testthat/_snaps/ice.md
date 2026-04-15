# ipsi() is explicitly rejected by ICE/contrast() (Phase 4 placeholder)

    Code
      contrast(fit, interventions = list(boost = ipsi(2), nat = ipsi(1)), ci_method = "sandwich")
    Condition
      Error in `apply_single_intervention()`:
      ! `ipsi()` interventions are only supported under `estimator = 'ipw'`.
      i The intervention shifts the propensity, not the treatment value, so there is no counterfactual treatment to predict at under g-computation or matching.
      i Use `causat(..., estimator = 'ipw')` with an IPSI intervention, or rewrite the intervention as a `shift()` / `scale_by()` / `static()` for g-comp / matching.

