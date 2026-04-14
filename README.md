
<!-- README.md is generated from README.Rmd. Please edit that file -->

# causatr

<!-- badges: start -->

[![R-CMD-check](https://github.com/etverse/causatr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/etverse/causatr/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/etverse/causatr/graph/badge.svg)](https://app.codecov.io/gh/etverse/causatr)
<!-- badges: end -->

**causatr** provides a unified interface for causal effect estimation
via three complementary methods: g-computation (parametric g-formula),
inverse probability weighting (IPW), and propensity score matching. When
all three methods agree, you can be more confident in your findings —
this is called **methodological triangulation**.

The package implements the methods described in Hernán & Robins (2025)
*Causal Inference: What If* with a simple two-step API:

1.  **Fit** the causal model with `causat()`
2.  **Contrast** interventions with `contrast()`

## Installation

Install the development version from GitHub:

``` r
# install.packages("pak")
pak::pak("etverse/causatr")
```

## Quick example

Estimate the average causal effect of quitting smoking on weight gain
using the NHEFS dataset from Hernán & Robins (2025):

``` r
library(causatr)
data("nhefs")

# Step 1: Fit the outcome model via g-computation
fit <- causat(
  nhefs,
  outcome = "wt82_71",
  treatment = "qsmk",
  confounders = ~ sex + age + I(age^2) + race + factor(education) +
    smokeintensity + I(smokeintensity^2) + smokeyrs + I(smokeyrs^2) +
    factor(exercise) + factor(active) + wt71 + I(wt71^2) +
    qsmk:smokeintensity,
  censoring = "censored"
)

# Step 2: Contrast interventions
result <- contrast(
  fit,
  interventions = list(quit = static(1), continue = static(0)),
  reference = "continue"
)
#> 117 row(s) with NA predictions excluded from the target population.
result
#> <causatr_result>
#>  Method:    G-computation
#>  Estimand:  ATE
#>  Contrast:  Difference
#>  CI method: sandwich
#>  N:         1629
#> 
#> Intervention means:
#>    intervention estimate     se ci_lower ci_upper
#>          <char>    <num>  <num>    <num>    <num>
#> 1:         quit    5.176 0.4367     4.32    6.032
#> 2:     continue    1.660 0.2195     1.23    2.090
#> 
#> Contrasts:
#>          comparison estimate     se ci_lower ci_upper
#>              <char>    <num>  <num>    <num>    <num>
#> 1: quit vs continue    3.516 0.4791    2.577    4.455
```

## Methodological triangulation

Compare g-computation, IPW, and matching on the same data:

``` r
# G-computation (outcome model)
fit_gc <- causat(nhefs, outcome = "wt82_71", treatment = "qsmk",
  confounders = ~ sex + age + race + smokeintensity + smokeyrs +
    factor(exercise) + factor(active) + wt71,
  censoring = "censored")

# IPW (treatment model)
fit_ipw <- causat(nhefs, outcome = "wt82_71", treatment = "qsmk",
  confounders = ~ sex + age + race + smokeintensity + smokeyrs +
    factor(exercise) + factor(active) + wt71,
  estimator = "ipw")

# Matching (propensity score)
fit_m <- causat(nhefs, outcome = "wt82_71", treatment = "qsmk",
  confounders = ~ sex + age + race + smokeintensity + smokeyrs +
    factor(exercise) + factor(active) + wt71,
  estimator = "matching", estimand = "ATT")

# All three estimates
rbind(
  data.frame(estimator = "gcomp", contrast(fit_gc,
    list(quit = static(1), cont = static(0)), reference = "cont")$contrasts),
  data.frame(estimator = "ipw", contrast(fit_ipw,
    list(quit = static(1), cont = static(0)), reference = "cont")$contrasts),
  data.frame(estimator = "matching", contrast(fit_m,
    list(quit = static(1), cont = static(0)), reference = "cont")$contrasts)
)
#>     method   comparison estimate        se ci_lower ci_upper
#> 1    gcomp quit vs cont 3.155727 0.4487520 2.276190 4.035265
#> 2      ipw quit vs cont 3.205240 0.4513455 2.320619 4.089861
#> 3 matching quit vs cont 2.984411 0.5091996 1.986398 3.982424
```

## Intervention types

Beyond static interventions, g-computation supports modified treatment
policies:

``` r
fit_cont <- causat(nhefs, outcome = "wt82_71",
  treatment = "smokeintensity",
  confounders = ~ sex + age + race + wt71,
  censoring = "censored")

contrast(fit_cont,
  interventions = list(
    reduce10 = shift(-10),
    halved = scale_by(0.5),
    cap20 = threshold(0, 20),
    observed = NULL
  ),
  reference = "observed"
)
#> <causatr_result>
#>  Method:    G-computation
#>  Estimand:  ATE
#>  Contrast:  Difference
#>  CI method: sandwich
#>  N:         1746
#> 
#> Intervention means:
#>    intervention estimate     se ci_lower ci_upper
#>          <char>    <num>  <num>    <num>    <num>
#> 1:     reduce10    2.621 0.2403    2.150    3.092
#> 2:       halved    2.619 0.2432    2.142    3.095
#> 3:        cap20    2.664 0.1987    2.275    3.053
#> 4:     observed    2.699 0.1914    2.324    3.074
#> 
#> Contrasts:
#>              comparison estimate      se ci_lower ci_upper
#>                  <char>    <num>   <num>    <num>    <num>
#> 1: reduce10 vs observed -0.07803 0.16323  -0.3980   0.2419
#> 2:   halved vs observed -0.08028 0.16792  -0.4094   0.2488
#> 3:    cap20 vs observed -0.03488 0.07294  -0.1778   0.1081
```

## Diagnostics

Check covariate balance and positivity after fitting:

``` r
diag <- diagnose(fit_ipw)
diag          # positivity + balance summary
plot(diag)    # Love plot (requires cobalt)
```

## Features

- **Three estimation methods**: g-computation, IPW (via
  [WeightIt](https://ngreifer.github.io/WeightIt/)), matching (via
  [MatchIt](https://kosukeimai.github.io/MatchIt/)). Matching is
  binary-only; continuous / categorical treatments use gcomp or IPW.
- **Longitudinal support**: ICE g-computation (Zivich et al. 2024) for
  time-varying treatments with sandwich variance via stacked estimating
  equations, plus parallel bootstrap via `boot::boot()`.
- **Flexible interventions**: `static()`, `shift()`, `scale_by()`,
  `threshold()`, `dynamic()` for modified treatment policies. IPSI is
  scaffolded for Phase 4.
- **Any outcome family**: gaussian, binomial (logit / probit / cloglog),
  Poisson, quasibinomial (fractional responses), Gamma with log link,
  plus any family you pass through `model_fn`.
- **Pluggable models**: `stats::glm`, `mgcv::gam`, splines via `ns()` /
  `bs()`, or any fit function with signature
  `(formula, data, family, weights, ...)`. A two-tier numeric-variance
  fallback handles model classes without a `sandwich::estfun` method.
- **Robust inference**: analytic sandwich SE (default, via a unified
  influence-function engine) or nonparametric bootstrap with percentile
  CIs. `confint()` respects the `level` argument on both paths.
- **External weights**: survey / IPCW weights pass through every method
  and propagate to the sandwich variance.
- **Contrast types**: risk difference, risk ratio, odds ratio — ratio
  and OR use log-scale CIs.
- **Estimands**: ATE, ATT, ATC, or custom subgroups via `subset=` /
  `by=`.
- **Built-in diagnostics**: positivity checks, covariate balance via
  [cobalt](https://ngreifer.github.io/cobalt/), weight summaries, Love
  plots.
- **Tidy integration**: `tidy()` / `glance()` / `confint()` / `coef()` /
  `vcov()` / `plot()` (forest plot via `forrest`) / broom-compatible
  output.

## References

Hernán MA, Robins JM (2025). *Causal Inference: What If*. Chapman &
Hall/CRC.

## Acknowledgements

This package was built with the contribution of
[Claude](https://claude.ai), Anthropic’s AI assistant.
