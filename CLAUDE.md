# causatr

causatr estimates causal effects of point and time-varying treatments by g-computation (including iterated conditional expectation, ICE), inverse probability weighting (IPW), augmented IPW (AIPW), matching (through MatchIt) and structural nested models (SNM), with sandwich and bootstrap variance.

It is part of the etverse, whose design (https://etverse.github.io/design.html) fixes what causatr owns and the interface every package shares. causatr does not follow that design yet; until it is brought into line, this file and the package documentation describe what the code does today.

## API

```r
fit <- causat(data, outcome = "Y", treatment = "A", confounders = ~ L1 + L2,
              estimator = "gcomp")   # or "ipw", "aipw", "matching", "snm"

res <- contrast(fit,
  interventions = list(a1 = static(1), a0 = static(0)),
  type = "difference", ci_method = "sandwich")

diagnose(fit)
```

Longitudinal fits add `id` and `time`; transport adds `target`; per-model covariates use `confounders_outcome`, `confounders_treatment`, `confounders_censoring` and `confounders_sampling`. `causat_mice()` runs the same analysis over multiply imputed data.

S3 classes: `causatr_fit` (from `causat()`), `causatr_result` (from `contrast()`), `causatr_diag` (from `diagnose()`) and `causatr_intervention` (from the intervention constructors). Constructors are in `R/constructors.R`.

## Code layout

Files in `R/` are named by what they hold; `tests/testthat/test-foo.R` tests `R/foo.R`.

- Fit and contrast: `causat.R`, `contrast.R`, `prepare_data.R`, `checks.R`, `to_person_period.R`
- Interventions: `interventions.R`, `glmtp_interventions.R`
- Estimators: `gcomp.R`, `ice*.R`, `glmtp*.R`, `ipw.R`, `longitudinal_ipw*.R`, `aipw*.R`, `matching.R`, `snm*.R`
- Nuisance models and weights: `treatment_model.R`, `censoring.R`, `sampling_model.R`, `ipw_weights*.R`, `effect_modification.R`, `family.R`
- Sandwich variance: `variance_if_core.R` (shared primitives), `variance_if.R` (dispatcher and numeric fallback), then `variance_if_<estimator>*.R`
- Bootstrap: `variance_bootstrap*.R`
- Multiple imputation: `causat_mice.R`, `pool_rubin.R`, `pool_boot_mi.R`
- Diagnostics: `diagnose*.R`
- S3 methods: `print.R`, `summary.R`, `plot.R`, `coef.R`, `confint.R`, `tidy.R`, `knit_print.R`
- Other: `target_trial.R`, `utils.R`, `data.R`, `causatr-package.R`, `zzz.R`

Invariants that look like bugs and are not, and the axes a combination audit runs over, are in `.claude/rules/`; each file loads when you open the code it covers.

## Commands

The `Makefile` has `make document`, `make test`, `make check`, `make format` (Air), `make lint` and `make site` (altdoc). While iterating, run one file:

```bash
Rscript -e 'devtools::load_all(quiet = TRUE); testthat::test_file("tests/testthat/test-ipw.R")'
```

`CAUSATR_TEST_TIER=fast` skips the blocks guarded by `skip_if_fast()`; the default tier is `full`, which takes over an hour on four processes. Python oracle fixtures are regenerated with the virtual environments in `data-raw/` (`deli_venv` for `delicatessen`, `zepid_venv` for zEpid).

## Planned work

Planned work and known bugs are GitHub issues (`gh issue list`, `gh issue view <n>`), grouped in `ROADMAP.md`; the reason behind shipped code is in the roxygen of the function it explains. `FEATURE_COVERAGE_MATRIX.md` records which combinations are tested.

## Dependencies

Imports: boot, data.table, generics, MatchIt, nnet, numDeriv, rlang, sandwich. Test oracles in Suggests include WeightIt, lmtp, TransportHealth, DTRreg and marginaleffects. survatr and matchatr call causatr's code, so a change to an exported or internal function they use can break them.
