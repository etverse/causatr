# ipsi() is explicitly rejected by ICE/contrast() (Phase 4 placeholder)

    Code
      contrast(fit, interventions = list(boost = ipsi(2), nat = ipsi(1)), ci_method = "sandwich")
    Message
      `ipsi()` is currently dead-ended: the constructor succeeds but `contrast()` aborts.
      i IPSI requires a fitted propensity model that is not wired through any estimation engine yet (planned for Phase 4).
      i Use `shift()`, `scale_by()`, or `static()` with `estimator = 'gcomp'` in the meantime.
      This message is displayed once per session.
    Condition
      Error in `apply_single_intervention()`:
      ! IPSI interventions require a propensity model and are not yet supported in contrast().

