# diagnose() aborts when the WeightIt object is missing treat.type

    Code
      diagnose(fit)
    Condition
      Error in `compute_weight_summary()`:
      ! WeightIt object is missing the `treat.type` attribute. This indicates a non-standard or serialized WeightIt fit. Refit the model with `causat(..., method = 'ipw')` so causatr can label weight summaries correctly.

# diagnose() rejects longitudinal fits with a clear error

    Code
      diagnose(fit)
    Condition
      Error in `diagnose()`:
      ! `diagnose()` is not yet supported for longitudinal fits.
      i Per-period positivity and per-period balance tables will land in a future phase. For now, run `diagnose()` on a point-treatment subset of the data (e.g. baseline with `time == min(time)`).

