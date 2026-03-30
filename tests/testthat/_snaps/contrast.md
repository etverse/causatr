# contrast() rejects non-causatr_fit input

    Code
      contrast(list(), list(a = static(1)))
    Condition
      Error in `contrast()`:
      ! `fit` must be a `causatr_fit` object returned by `causat()`.

# contrast() rejects unnamed intervention list

    Code
      contrast(fit, list(static(1), static(0)))
    Condition
      Error in `contrast()`:
      ! All elements of `interventions` must be named.

# contrast() rejects non-static interventions for IPW

    Code
      contrast(fit, list(a1 = shift(1), a0 = static(0)))
    Condition
      Error in `contrast()`:
      ! Non-static interventions (shift, dynamic, scale, threshold, ipsi) are not supported for method = 'ipw'. The weights/matched sets were estimated under the original treatment regime and are not valid under a different intervention. Use method = 'gcomp' instead.

# contrast() rejects non-static interventions for matching

    Code
      contrast(fit, list(a1 = dynamic(function(d, a) 1), a0 = static(0)))
    Condition
      Error in `contrast()`:
      ! Non-static interventions (shift, dynamic, scale, threshold, ipsi) are not supported for method = 'matching'. The weights/matched sets were estimated under the original treatment regime and are not valid under a different intervention. Use method = 'gcomp' instead.

