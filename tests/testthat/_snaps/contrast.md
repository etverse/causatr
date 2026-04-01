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

# contrast() rejects estimand and subset together

    Code
      contrast(fit, list(a1 = static(1), a0 = static(0)), estimand = "ATT", subset = quote(
        A == 1))
    Condition
      Error in `contrast()`:
      ! Specify either 'estimand' or 'subset', not both.

# contrast() aborts when IPW estimand is changed

    Code
      contrast(fit, list(a1 = static(1), a0 = static(0)), estimand = "ATT")
    Condition
      Error in `contrast()`:
      ! For method = 'ipw', the estimand is fixed at fitting time because it determines the weights. Refit with causat(estimand = 'ATT').

# contrast() aborts when matching estimand is changed

    Code
      contrast(fit, list(a1 = static(1), a0 = static(0)), estimand = "ATE")
    Condition
      Error in `contrast()`:
      ! For method = 'matching', the estimand is fixed at fitting time because it determines the weights. Refit with causat(estimand = 'ATE').

# contrast() rejects ATT for longitudinal fit

    Code
      contrast(fit, list(a1 = static(1), a0 = static(0)), estimand = "ATT")
    Condition
      Error in `contrast()`:
      ! estimand = 'ATT' is only defined for binary point treatments. Use estimand = 'ATE' or subset = quote(...) for subgroup effects.

# contrast() rejects reference not in interventions

    Code
      contrast(fit, list(a1 = static(1), a0 = static(0)), reference = "a2")
    Condition
      Error in `contrast()`:
      ! `reference` ('a2') must be the name of one of the interventions.

# contrast() rejects multivariate intervention, missing trt var

    Code
      contrast(fit, list(a1 = list(static(1), static(0))))
    Condition
      Error in `contrast()`:
      ! `interventions$a1` is a list but not all elements are named. For multivariate treatment, supply a named list with one entry per treatment variable.

