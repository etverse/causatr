# IPW rejects static(1) on a continuous treatment

    Code
      contrast(fit, interventions = list(a = static(1), b = static(0)), ci_method = "sandwich")
    Condition
      Error in `check_intervention_family_compat()`:
      ! `static(v)` on a continuous treatment is degenerate for IPW.
      i No observations lie exactly at `v`, so the Horvitz-Thompson weight is zero almost surely.
      i Use `shift()` or `scale_by()` to move the whole treatment distribution, or switch to `estimator = 'gcomp'`.

