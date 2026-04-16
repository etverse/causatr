# matching aborts on categorical (k > 2) treatment with a clear error

    Code
      causat(df, outcome = "Y", treatment = "A", confounders = ~L, estimator = "matching")
    Condition
      Error in `fit_matching()`:
      ! Matching supports only binary treatments, but `A` has 3 levels. Use `estimator = "gcomp"` or `estimator = "ipw"` for categorical treatments.

