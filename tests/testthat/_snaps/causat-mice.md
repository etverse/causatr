# causat_mice rejects non-mids input

    Code
      causat_mice(list(), outcome = "Y", treatment = "A", confounders = ~L,
      interventions = list(a1 = static(1)))
    Condition
      Error in `causat_mice()`:
      ! `imp` must be a `mids` object returned by `mice::mice()`.
      x Got an object of class: list.
      i Run `mice::mice(data, ...)` first, then pass the result here.

