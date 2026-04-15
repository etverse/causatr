# IPW rejects static(1) on a continuous treatment

    Code
      contrast(fit, interventions = list(a = static(1), b = static(0)), ci_method = "sandwich")
    Condition
      Error in `check_intervention_family_compat()`:
      ! `static(v)` on a continuous treatment is degenerate for IPW.
      i No observations lie exactly at `v`, so the Horvitz-Thompson weight is zero almost surely.
      i Use `shift()` or `scale_by()` to move the whole treatment distribution, or switch to `estimator = 'gcomp'`.

# IPW rejects A:modifier interaction terms in confounders

    Code
      causat(d, outcome = "Y", treatment = "A", confounders = ~ L + sex + A:sex,
      estimator = "ipw")
    Condition
      Error in `check_confounders_no_treatment()`:
      ! `confounders` contains term(s) involving the treatment, which are not supported for `estimator = "ipw"`.
      x Offending term(s): sex:A.
      i IPW and matching wrap a saturated MSM `Y ~ A` around the propensity/match model, so treatment-by-modifier interactions cannot be estimated here.
      i Use `estimator = "gcomp"` for heterogeneous treatment effects, or `by = "modifier"` in `contrast()` for stratum-specific summaries of a homogeneous effect.

# IPW rejects bare treatment in confounders

    Code
      causat(d, outcome = "Y", treatment = "A", confounders = ~ L + A, estimator = "ipw")
    Condition
      Error in `check_confounders_no_treatment()`:
      ! `confounders` contains term(s) involving the treatment, which are not supported for `estimator = "ipw"`.
      x Offending term(s): A.
      i IPW and matching wrap a saturated MSM `Y ~ A` around the propensity/match model, so treatment-by-modifier interactions cannot be estimated here.
      i Use `estimator = "gcomp"` for heterogeneous treatment effects, or `by = "modifier"` in `contrast()` for stratum-specific summaries of a homogeneous effect.

