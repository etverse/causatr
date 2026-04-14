# g-comp rejects unknown family string

    Code
      causat(df, outcome = "Y", treatment = "A", confounders = ~L, family = "not_a_family")
    Condition
      Error in `get()`:
      ! object 'not_a_family' of mode 'function' was not found

