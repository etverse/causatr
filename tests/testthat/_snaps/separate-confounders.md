# Error on non-formula confounders_outcome

    Code
      causat(d, outcome = "Y", treatment = "A", confounders_outcome = "L1")
    Condition
      Error in `causat()`:
      ! `confounders_outcome` must be a formula (e.g. `~ L1 + L2`).

# Error on missing column in confounders_outcome

    Code
      causat(d, outcome = "Y", treatment = "A", confounders_outcome = ~nonexistent)
    Condition
      Error in `causat()`:
      ! confounders_outcome variable(s) not found in `data`: nonexistent

# Error when AIPW missing required confounders_treatment

    Code
      causat(d, outcome = "Y", treatment = "A", confounders_outcome = ~L1, estimator = "aipw",
        propensity_model_fn = stats::glm)
    Condition
      Error in `causat()`:
      ! `estimator = "aipw"` requires treatment-model confounders. Supply `confounders_treatment` or `confounders`.

# Error when IPW has no treatment confounders at all

    Code
      causat(d, outcome = "Y", treatment = "A", estimator = "ipw",
        propensity_model_fn = stats::glm)
    Condition
      Error in `causat()`:
      ! `estimator = "ipw"` requires treatment-model confounders. Supply `confounders_treatment` or `confounders`.

