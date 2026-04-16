# shift(0.5) on count treatment aborts

    Code
      contrast(fit, interventions = list(s = shift(0.5)), ci_method = "sandwich")
    Condition
      Error in `check_intervention_family_compat()`:
      ! `shift(0.5)` is not integer-valued.
      i Count treatments (poisson) require integer shift deltas because `dpois()` / `dnbinom()` return 0 at non-integer arguments.
      i Use `estimator = 'gcomp'` for fractional shifts on count data.

# scale_by(2) on count treatment aborts (inverse not integer)

    Code
      contrast(fit, interventions = list(s = scale_by(2)), ci_method = "sandwich")
    Condition
      Error in `check_intervention_family_compat()`:
      ! `scale_by(2)` does not produce integer inverse values (A / 2) for all observed treatment values.
      i Count treatments (poisson) require that A / factor is a non-negative integer for every observation, because the density ratio evaluates the pmf at A / factor.
      i Use `estimator = 'gcomp'` for non-integer-preserving scales on count data.

# static() on count treatment aborts

    Code
      contrast(fit, interventions = list(s = static(2)), ci_method = "sandwich")
    Condition
      Error in `check_intervention_family_compat()`:
      ! `static(v)` on a count treatment is degenerate for IPW.
      i The Horvitz-Thompson indicator weight is zero for almost all observations.
      i Use `shift()` for integer shifts, or switch to `estimator = 'gcomp'`.

# dynamic() on count treatment aborts

    Code
      contrast(fit, interventions = list(s = dynamic(function(data, a) pmax(a - 1L,
      0L))), ci_method = "sandwich")
    Condition
      Error in `check_intervention_family_compat()`:
      ! `dynamic()` rules on count treatments are not supported by the IPW engine.
      i Use `shift()` for integer shifts, or switch to `estimator = 'gcomp'` for deterministic rules.

# threshold() on count treatment aborts

    Code
      contrast(fit, interventions = list(s = threshold(0, 5)), ci_method = "sandwich")
    Condition
      Error in `check_intervention_family_compat()`:
      ! `threshold()` on a count treatment is not supported by the IPW engine.
      i The pushforward of a count density under a boundary clamp is a mixed measure.
      i Use `estimator = 'gcomp'`.

# ipsi() on count treatment aborts

    Code
      contrast(fit, interventions = list(s = ipsi(2)), ci_method = "sandwich")
    Condition
      Error in `check_intervention_family_compat()`:
      ! `ipsi()` interventions are only defined for binary (0/1) treatments.
      i The treatment column is classified as 'poisson'. Use `shift()` or `scale_by()` for a continuous MTP instead.

