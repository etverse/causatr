# ipsi() is explicitly rejected by ICE/contrast() (Phase 4 placeholder)

    Code
      contrast(fit, interventions = list(boost = ipsi(2), nat = ipsi(1)), ci_method = "sandwich")
    Condition
      Error in `apply_single_intervention()`:
      ! IPSI interventions require a propensity model and are not yet supported in contrast().

