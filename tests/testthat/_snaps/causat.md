# causat() rejects missing outcome column

    Code
      causat(df, outcome = "Y", treatment = "A", confounders = ~L)
    Condition
      Error in `causat()`:
      ! Column `Y` (outcome) not found in `data`.

# causat() rejects missing treatment column

    Code
      causat(df, outcome = "Y", treatment = "A", confounders = ~L)
    Condition
      Error in `causat()`:
      ! Column `A` (treatment) not found in `data`.

# causat() rejects missing confounder column

    Code
      causat(df, outcome = "Y", treatment = "A", confounders = ~L)
    Condition
      Error in `causat()`:
      ! Confounder variable(s) not found in `data`: L

# causat() rejects id without time

    Code
      causat(df, outcome = "Y", treatment = "A", confounders = ~L, id = "id")
    Condition
      Error in `causat()`:
      ! Both `id` and `time` must be provided together for longitudinal data.

# causat() rejects time without id

    Code
      causat(df, outcome = "Y", treatment = "A", confounders = ~L, time = "t")
    Condition
      Error in `causat()`:
      ! Both `id` and `time` must be provided together for longitudinal data.

# causat() rejects ATT for continuous treatment

    Code
      causat(df, outcome = "Y", treatment = "A", confounders = ~L, estimand = "ATT")
    Condition
      Error in `causat()`:
      ! estimand = 'ATT' is only defined for binary point treatments. Use estimand = 'ATE' or subset = quote(...) for subgroup effects.

# causat() rejects ATT for multivariate treatment

    Code
      causat(df, outcome = "Y", treatment = c("A1", "A2"), confounders = ~L,
      estimand = "ATT")
    Condition
      Error in `causat()`:
      ! estimand = 'ATT' is only defined for binary point treatments. Use estimand = 'ATE' or subset = quote(...) for subgroup effects.

# causat() aborts when treatment has NAs and no censoring

    Code
      causat(df, outcome = "Y", treatment = "A", confounders = ~L)
    Condition
      Error in `causat()`:
      ! Treatment variable 'A' has 1 missing value.
      i Use `censoring = '...'` for inverse probability of censoring weights.
      i Use `causat_mice()` with a mice `mids` object for multiple imputation.
      i Or remove incomplete cases before calling `causat()`.

# causat() rejects missing confounders_tv column

    Code
      causat(df, outcome = "Y", treatment = "A", confounders = ~L, confounders_tv = ~
        CD4, id = "id", time = "time")
    Condition
      Error in `causat()`:
      ! Time-varying confounder variable(s) not found in `data`: CD4

# causat() rejects invalid history value

    Code
      causat(df, outcome = "Y", treatment = "A", confounders = ~L, id = "id", time = "time",
        history = 0)
    Condition
      Error in `causat()`:
      ! `history` must be a positive integer or `Inf`.

