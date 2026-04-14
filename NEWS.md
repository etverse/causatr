# causatr (development version)

## Bug fixes

* ICE sandwich variance now propagates external weights correctly
  through the cascade gradient. The previous implementation dropped
  the `w_{k-1}` factor when building `g_k^eff = A_{k-1,k}^T h_{k-1}`,
  silently underestimating the SE by ~2x under non-uniform external
  weights on multi-period DGPs. Verified against a Monte-Carlo truth
  and `lmtp::lmtp_tmle()`.
* Fixed a false-positive "Baseline confounder 'A' varies within some
  individuals" warning when users passed treatment-by-modifier
  interactions like `confounders = ~ sex + A:sex`. Treatment columns
  are now excluded from the baseline-variation check.
* `diagnose()` on a longitudinal fit used to silently fall through
  to the point-treatment path and crash inside cobalt's printer.
  It now rejects up front with a clear message pointing users at
  running `diagnose()` on a baseline-only subset.
* `matching` rejects categorical (k > 2 levels) treatments upstream
  with a clear error pointing users to `method = "gcomp"` or
  `method = "ipw"`, instead of propagating `MatchIt`'s internal
  "the treatment must be a binary variable" error.
* `causat()` and `causat_survival()` now validate external weights
  up front and reject NA, Inf / NaN, negative, non-numeric, or
  mis-sized weight vectors with a specific error.
* `contrast()` rejects duplicated intervention names (previously
  silently accepted and produced stale rows in the contrasts table).

## New tests

* Truth-based coverage for previously-unpinned combinations:
  * gcomp × quasibinomial and Gamma(log) outcomes
  * gcomp × binomial with probit and cloglog links
  * gcomp × external weights × bootstrap (upgraded smoke → truth)
  * gcomp × Poisson × ratio × bootstrap
  * gcomp × categorical (3-level) treatment
  * IPW × categorical (3-level) treatment (shares DGP with gcomp)
  * IPW × continuous treatment × static at specific levels
  * IPW × survey weights × sandwich
  * IPW × ATT × bootstrap
  * IPW × non-Mparts WeightIt method (via stubbed weightit) at
    contrast time — fit-time warning already existed
  * matching × Poisson and Gamma(log) outcomes
  * matching × ATC × bootstrap
  * ICE × continuous TV treatment × `scale_by()` and `threshold()`
    against analytic truth (`2 * (c - 1)` and `2 * dnorm(0)`)
  * ICE × continuous TV treatment × `shift()` × bootstrap
  * ICE × binary TV × static × bootstrap × survey weights
  * ICE × binomial outcome × ratio and OR contrasts
  * ICE × by(sex) stratification
  * Survival × censoring (fit-only smoke; Phase 6 contrast still
    pending)
* End-to-end Tier 2 numeric variance test via a custom `model_fn`
  with no `sandwich::estfun` method.
* Direct unit test for `vcov_from_if(cluster = ...)` against a
  hand-computed sum-then-square reference.
* `confint()` level monotonicity on both the sandwich and bootstrap
  paths.
* Snapshot rejection tests for every new input validation branch
  (weights, intervention names, longitudinal diagnose, categorical
  matching, missing WeightIt `treat.type`).

## Documentation

* Rewrote the single-item `Getting started` and `Theory` sections in
  `altdoc/quarto_website.yml` to use the explicit `text: / file:`
  form. The previous bare `- vignettes/foo.qmd` form was silently
  collapsed by `yaml::yaml.load()` -> `yaml::as.yaml()` into a
  scalar, breaking the sidebar link and hiding the
  `variance-theory.qmd` vignette from the deployed site.
* `vignettes/variance-theory.qmd` Section 5.4 now derives the
  cascade gradient with the explicit `w_{k-1}` factor that the
  engine needs under external weights, and the pseudo-code shows
  the `prior.weights` lookup.
* `causat()` and `vignette("matching")` now explicitly state the
  binary-only restriction for matching (verified against the
  upstream MatchIt documentation).

## Coverage matrix

* `FEATURE_COVERAGE_MATRIX.md` is now the single source of truth
  for "what works" across every supported feature combination.
  Every cell has been updated to reflect the tests above, and a
  new "Planned coverage" section enumerates the target rows for
  Phase 4 (self-contained IPW), Phase 6 (survival contrasts and
  competing risks, including the NHEFS Ch. 17 replication target),
  and Phase 7 (survey integration, future parallel backend,
  multinomial / ordinal outcomes, ICE lag × modifier auto-expansion).

## Known limitations

* `ice_build_formula()` resolves the `A:modifier` interaction only
  for the current-period treatment slot. Lagged treatments
  (`lag1_A`, ...) do not get the auto-expanded modifier term, so a
  2-period DGP with a persistent `A × sex` interaction will have
  its heterogeneity compressed across strata. Workaround: use
  wide-format + point gcomp with an explicit `A1:sex` interaction
  until the ICE formula-template expansion lands.
