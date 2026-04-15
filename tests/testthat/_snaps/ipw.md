# IPW rejects shift intervention until Phase 4

    Code
      contrast(fit, interventions = list(s = shift(1), c = static(0)), ci_method = "sandwich")
    Condition
      Error in `contrast()`:
      ! Non-static interventions (shift, dynamic, scale, threshold, ipsi) are not supported for estimator = 'ipw'. The weights/matched sets were estimated under the original treatment regime and are not valid under a different intervention. Use estimator = 'gcomp' instead.

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

