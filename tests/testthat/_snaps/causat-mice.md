# causat_mice() aborts when mice is not installed gracefully

    Code
      causat_mice(list(), outcome = "Y", treatment = "A", confounders = ~L,
      interventions = list(a1 = static(1)))
    Condition
      Error in `check_pkg()`:
      ! Package 'mice' is required but not installed. Install it with: install.packages('mice')

