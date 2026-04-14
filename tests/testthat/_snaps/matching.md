# matching aborts on categorical (k > 2) treatment with a clear error

    Code
      causat(df, outcome = "Y", treatment = "A", confounders = ~L, method = "matching")
    Condition
      Error in `fit_matching()`:
      ! Matching supports only binary treatments, but `A` has 3 levels. Use `method = "gcomp"` or `method = "ipw"` for categorical treatments (Phase 4 / Phase 7 will revisit multi-category matching).

# matching rejects A:modifier interaction terms in confounders

    Code
      causat(d, outcome = "Y", treatment = "A", confounders = ~ L + sex + A:sex,
      method = "matching")
    Condition
      Error in `check_confounders_no_treatment()`:
      ! `confounders` contains term(s) involving the treatment, which are not supported for `method = "matching"`.
      x Offending term(s): sex:A.
      i IPW and matching wrap a saturated MSM `Y ~ A` around the propensity/match model, so treatment-by-modifier interactions cannot be estimated here.
      i Use `method = "gcomp"` for heterogeneous treatment effects, or `by = "modifier"` in `contrast()` for stratum-specific summaries of a homogeneous effect.
      i See `PHASE_8_INTERACTIONS.md` for the planned unified API.

