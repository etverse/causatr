# by parameter rejects missing variable

    Code
      contrast(fit, list(a1 = static(1), a0 = static(0)), by = "nonexistent")
    Condition
      Error in `contrast()`:
      ! `by` variable 'nonexistent' not found in fitted data.

# survival type aborts in contrast()

    Code
      contrast(fit, list(a1 = static(1), a0 = static(0)))
    Condition
      Error in `compute_contrast()`:
      ! Survival curve estimation via contrast() is not yet implemented. causat_survival() currently fits a pooled logistic model; survival curve contrasts are planned for a future release.

# multivariate treatment blocked for IPW

    Code
      causat(data.frame(Y = 1:10, A1 = 1:10, A2 = 1:10, L = 1:10), outcome = "Y",
      treatment = c("A1", "A2"), confounders = ~L, estimator = "ipw")
    Condition
      Error in `causat()`:
      ! Multivariate treatments are not yet supported for estimator = 'ipw'. Use estimator = 'gcomp' for joint interventions on multiple treatments, or fit separate models for each treatment.

# multivariate treatment blocked for matching

    Code
      causat(data.frame(Y = 1:10, A1 = 1:10, A2 = 1:10, L = 1:10), outcome = "Y",
      treatment = c("A1", "A2"), confounders = ~L, estimator = "matching")
    Condition
      Error in `causat()`:
      ! Multivariate treatments are not yet supported for estimator = 'matching'. Use estimator = 'gcomp' for joint interventions on multiple treatments, or fit separate models for each treatment.

