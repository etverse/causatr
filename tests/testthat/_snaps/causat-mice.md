# causat_mice() rejects non-mids input

    Code
      causat_mice(list(), outcome = "Y", treatment = "A", confounders = ~L,
      interventions = list(a1 = static(1)))
    Condition
      Error in `causat_mice()`:
      ! `imp` must be a `mids` object returned by `mice::mice()`. Got an object of class: list

# causat_mice() errors on unimplemented stub

    Code
      causat_mice(imp, outcome = "Y", treatment = "A", confounders = ~L,
        interventions = list(a1 = static(1)))
    Condition
      Error in `causat_mice()`:
      ! causat_mice() is not yet implemented.

