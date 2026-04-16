# treatment NAs without censoring → abort (gcomp)

    Code
      causat(d, outcome = "Y", treatment = "A", confounders = ~L)
    Condition
      Error in `causat()`:
      ! Treatment variable 'A' has 3 missing values.
      i Use `censoring = '...'` for inverse probability of censoring weights.
      i Or remove incomplete cases before calling `causat()`.

# treatment NAs without censoring → abort (ipw)

    Code
      causat(d, outcome = "Y", treatment = "A", confounders = ~L, estimator = "ipw")
    Condition
      Error in `causat()`:
      ! Treatment variable 'A' has 3 missing values.
      i Use `censoring = '...'` for inverse probability of censoring weights.
      i Or remove incomplete cases before calling `causat()`.

# treatment NAs without censoring → abort (matching)

    Code
      causat(d, outcome = "Y", treatment = "A", confounders = ~L, estimator = "matching")
    Condition
      Error in `causat()`:
      ! Treatment variable 'A' has 3 missing values.
      i Use `censoring = '...'` for inverse probability of censoring weights.
      i Or remove incomplete cases before calling `causat()`.

