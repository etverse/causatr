# matching aborts on categorical (k > 2) treatment with a clear error

    Code
      causat(df, outcome = "Y", treatment = "A", confounders = ~L, method = "matching")
    Condition
      Error in `fit_matching()`:
      ! Matching supports only binary treatments, but `A` has 3 levels. Use `method = "gcomp"` or `method = "ipw"` for categorical treatments (Phase 4 / Phase 7 will revisit multi-category matching).

