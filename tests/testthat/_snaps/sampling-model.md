# check_transport_inputs: rejects matching estimator

    Code
      check_transport_inputs("S", d$S, "target", "matching")
    Condition
      Error:
      ! `target` is not supported with `estimator = "matching"`.
      i Transportability does not compose cleanly with matching. Use `estimator = "gcomp"` or `"ipw"` instead.

# check_transport_inputs: rejects NA in target column

    Code
      check_transport_inputs("S", d$S, "target", "gcomp")
    Condition
      Error:
      ! Column `S` contains NA values.
      i The sampling indicator S must be complete (no NAs) for transportability.

# check_transport_inputs: rejects non-binary target column

    Code
      check_transport_inputs("S", d$S, "target", "gcomp")
    Condition
      Error:
      ! `S` must be binary (0 = target, 1 = study). Found values: 0, 1, 2.
      i Recode the sampling indicator to 0/1 before calling `causat()`.

# check_transport_inputs: rejects degenerate target (all S=1)

    Code
      check_transport_inputs("S", d$S, "target", "gcomp")
    Condition
      Error:
      ! `S` must have both study (S = 1) and target (S = 0) rows. Found only: 1.
      i The sampling model requires rows from both populations.

# causat: matching + target aborts

    Code
      causat(d_study, outcome = "Y", treatment = "A", confounders = ~L, estimator = "matching",
        target = "S", method = "nearest")
    Condition
      Error in `causat()`:
      ! `target` is not supported with `estimator = "matching"`.
      i Transportability does not compose cleanly with matching. Use `estimator = "gcomp"` or `"ipw"` instead.

