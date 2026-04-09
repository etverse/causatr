
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
#> 
#> Attaching package: 'causatr'
#> The following object is masked from 'package:base':
#> 
#>     scale
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
result
#> <causatr_result>
#>  Method:      gcomp
#>  Contrast:    difference
#>  CI method:   sandwich
#>  N:           1629
#> 
#> Intervention means:
#>    intervention estimate     se ci_lower ci_upper
#>          <char>    <num>  <num>    <num>    <num>
#> 1:         quit    5.176 0.4298    4.333    6.018
#> 2:     continue    1.660 0.2063    1.256    2.065
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
  method = "ipw")

# Matching (propensity score)
fit_m <- causat(nhefs, outcome = "wt82_71", treatment = "qsmk",
  confounders = ~ sex + age + race + smokeintensity + smokeyrs +
    factor(exercise) + factor(active) + wt71,
  method = "matching", estimand = "ATT")

# All three estimates
rbind(
  data.frame(method = "gcomp", contrast(fit_gc,
    list(quit = static(1), cont = static(0)), reference = "cont")$contrasts),
  data.frame(method = "ipw", contrast(fit_ipw,
    list(quit = static(1), cont = static(0)), reference = "cont")$contrasts),
  data.frame(method = "matching", contrast(fit_m,
    list(quit = static(1), cont = static(0)), reference = "cont")$contrasts)
)
#>     method   comparison estimate        se ci_lower ci_upper
#> 1    gcomp quit vs cont 3.155727 0.4487520 2.276190 4.035265
#> 2      ipw quit vs cont 3.205240 0.4693563 2.285318 4.125162
#> 3 matching quit vs cont 2.984411 0.5097832 1.985254 3.983568
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
    halved = scale(0.5),
    cap20 = threshold(0, 20),
    observed = NULL
  ),
  reference = "observed"
)
#> <causatr_result>
#>  Method:      gcomp
#>  Contrast:    difference
#>  CI method:   sandwich
#>  N:           1746
#> 
#> Intervention means:
#>    intervention estimate     se ci_lower ci_upper
#>          <char>    <num>  <num>    <num>    <num>
#> 1:     reduce10    2.621 0.2337    2.163    3.079
#> 2:       halved    2.619 0.2367    2.155    3.083
#> 3:        cap20    2.664 0.1906    2.290    3.038
#> 4:     observed    2.699 0.1828    2.341    3.057
#> 
#> Contrasts:
#>              comparison estimate      se ci_lower ci_upper
#>                  <char>    <num>   <num>    <num>    <num>
#> 1: reduce10 vs observed -0.07803 0.16323  -0.3980   0.2419
#> 2:   halved vs observed -0.08028 0.16794  -0.4094   0.2489
#> 3:    cap20 vs observed -0.03488 0.07297  -0.1779   0.1081
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
  [MatchIt](https://kosukeimai.github.io/MatchIt/))
- **Built-in diagnostics**: positivity checks, covariate balance via
  [cobalt](https://ngreifer.github.io/cobalt/), weight summaries, Love
  plots
- **Flexible interventions**: `static()`, `shift()`, `scale()`,
  `threshold()`, `dynamic()` for modified treatment policies
- **Robust inference**: sandwich SE (default) or nonparametric bootstrap
- **Any outcome type**: continuous, binary, count (via `family`
  parameter)
- **Pluggable models**: pass `stats::glm`, `mgcv::gam`, or any
  compatible fitting function via `model_fn`
- **Contrast types**: risk difference, risk ratio, odds ratio
- **Estimands**: ATE, ATT, ATC, or custom subgroups via `subset`

## References

Hernán MA, Robins JM (2025). *Causal Inference: What If*. Chapman &
Hall/CRC.
