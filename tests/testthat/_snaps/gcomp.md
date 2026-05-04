# g-comp rejects unknown family string

    Code
      causat(df, outcome = "Y", treatment = "A", confounders = ~L, family = "not_a_family")
    Condition
      Error in `resolve_family()`:
      ! Unknown family: 'not_a_family'. Supported: gaussian, binomial, poisson, quasibinomial, quasipoisson, Gamma, inverse.gaussian, beta. Or pass a family object directly.

