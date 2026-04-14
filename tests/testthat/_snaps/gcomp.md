# g-comp rejects unknown family string

    Code
      causat(df, outcome = "Y", treatment = "A", confounders = ~L, family = "not_a_family")
    Condition
      Error in `value[[3L]]()`:
      ! Unknown family: 'not_a_family'.

