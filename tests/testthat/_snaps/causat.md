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

# causat() errors on unimplemented gcomp

    Code
      causat(df, outcome = "Y", treatment = "A", confounders = ~L)
    Condition
      Error in `fit_gcomp_point()`:
      ! G-computation for point treatments is not yet implemented.

