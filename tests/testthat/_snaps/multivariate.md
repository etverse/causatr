# mv point: ATT/ATC rejected

    Code
      causat(df, "Y", c("A1", "A2"), ~L, estimand = "ATT")
    Condition
      Error in `causat()`:
      ! estimand = 'ATT' is only defined for binary point treatments coded as 0/1. Use estimand = 'ATE' or subset = quote(...) for subgroup effects (or recode the treatment as integer 0/1 if it already has two levels).

# mv point: IPW rejected

    Code
      causat(df, "Y", c("A1", "A2"), ~L, method = "ipw")
    Condition
      Error in `causat()`:
      ! Multivariate treatments are not yet supported for method = 'ipw'. Use method = 'gcomp' for joint interventions on multiple treatments, or fit separate models for each treatment.

# mv point: matching rejected

    Code
      causat(df, "Y", c("A1", "A2"), ~L, method = "matching")
    Condition
      Error in `causat()`:
      ! Multivariate treatments are not yet supported for method = 'matching'. Use method = 'gcomp' for joint interventions on multiple treatments, or fit separate models for each treatment.

