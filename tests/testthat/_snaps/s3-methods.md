# by parameter rejects missing variable

    Code
      contrast(fit, list(a1 = static(1), a0 = static(0)), by = "nonexistent")
    Condition
      Error in `contrast()`:
      ! `by` variable 'nonexistent' not found in fitted data.

# multivariate treatment blocked for IPW

    Code
      causat(data.frame(Y = 1:10, A1 = 1:10, A2 = 1:10, L = 1:10), outcome = "Y",
      treatment = c("A1", "A2"), confounders = ~L, estimator = "ipw")
    Condition
      Error in `causat()`:
      ! Multivariate treatments are not supported for estimator = 'ipw'. Use estimator = 'gcomp' for joint interventions on multiple treatments, or fit separate models for each treatment.

# multivariate treatment blocked for matching

    Code
      causat(data.frame(Y = 1:10, A1 = 1:10, A2 = 1:10, L = 1:10), outcome = "Y",
      treatment = c("A1", "A2"), confounders = ~L, estimator = "matching")
    Condition
      Error in `causat()`:
      ! Multivariate treatments are not supported for estimator = 'matching'. Use estimator = 'gcomp' for joint interventions on multiple treatments, or fit separate models for each treatment.

