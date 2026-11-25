---
paths:
  - "R/causat.R"
  - "R/contrast.R"
  - "R/checks.R"
  - "R/prepare_data.R"
  - "R/interventions.R"
  - "R/effect_modification.R"
  - "R/variance_bootstrap*.R"
  - "tests/testthat/test-causat.R"
  - "tests/testthat/test-checks.R"
  - "tests/testthat/test-contrast.R"
  - "tests/testthat/test-replay-fit.R"
---

# Fitting and input checks

These are deliberate. Do not report them as bugs without a numerical reproducer.

- `estimator`, not `method`: `method` is passed through `...` to `MatchIt::matchit(method = )`. Never rename.
- Weights live in `fit$details$weights`, never as a column of the data.
- Every fit path that supports the bootstrap refits through `replay_fit()`, which replays `...`; a new fit path wires into it.
- `na.action = na.exclude` is rejected by `check_dots_na_action()` (class `causatr_bad_na_action`); only `na.omit` and `na.fail` are accepted. `na.exclude` pads residuals while `model.matrix()` drops rows, which silently corrupts the influence function, so do not "harden" the pipeline to accept it.
- `censoring =` alone is a row filter: no censoring model is fitted. `ipcw = TRUE` fits the censoring model and computes IPCW weights.
- ATT and ATC are defined only for static interventions on a binary 0/1 treatment; `check_estimand_intervention_compat()` enforces it (class `causatr_bad_estimand_intervention`).
- An effect modifier under IPW, matching and longitudinal IPW must be a baseline covariate, because a marginal structural model that conditions on a post-treatment variable mixes effect modification with selection through mediators and collider bias (Robins 2000). Under ICE the modifier at step k must come from the history before the treatment it modifies, L̄ₖ₋₁. This is documented, not checked, because time-varying status cannot be inferred from the data; it is not a missing check. SNM lifts the restriction: time-varying effect modification is what `estimator = "snm"` is for.
- An explicit `cluster =` overrides the primary sampling unit that `causat()` adopts from a `survey::svydesign` passed as `weights`. Matching rejects a user `cluster` (`causatr_bad_cluster`), because its matched subclasses are the cluster.
- `dynamic()` is a deterministic rule per individual, not a stochastic modified treatment policy. Under IPW it rejects a continuous treatment, because a point mass has no Lebesgue density. Modified treatment policies in the Díaz and Kennedy sense are `shift()`, `scale_by()` and `ipsi()`.
- Bootstrap intervals: `boot_ci = "percentile"` and `"normal"` both come from the same stored replicate matrix; the point estimate, SE and vcov do not depend on the choice. For `ci_method = "sandwich"` the argument is validated and unused.
