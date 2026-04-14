# IPW rejects shift intervention until Phase 4

    Code
      contrast(fit, interventions = list(s = shift(1), c = static(0)), ci_method = "sandwich")
    Condition
      Error in `contrast()`:
      ! Non-static interventions (shift, dynamic, scale, threshold, ipsi) are not supported for method = 'ipw'. The weights/matched sets were estimated under the original treatment regime and are not valid under a different intervention. Use method = 'gcomp' instead.

