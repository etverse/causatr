
<!-- README.md is generated from README.Rmd. Please edit that file -->

# causatr

<!-- badges: start -->

[![R-CMD-check](https://github.com/etverse/causatr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/etverse/causatr/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/etverse/causatr/graph/badge.svg)](https://app.codecov.io/gh/etverse/causatr)
<!-- badges: end -->

**causatr** provides a unified interface for causal effect estimation
via three complementary methods: g-computation (parametric g-formula +
ICE), inverse probability weighting (IPW with a self-contained
density-ratio engine), and propensity score matching (via
[MatchIt](https://kosukeimai.github.io/MatchIt/)). When all three
methods agree, you can be more confident in your findings — this is
called **methodological triangulation**.

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
```

**Estimator:** gcomp  ·  **Estimand:** ATE  ·  **Contrast:** difference
 ·  **CI method:** sandwich  ·  **N:** 1629

| intervention | estimate | se    | ci_lower | ci_upper |
|--------------|----------|-------|----------|----------|
| quit         | 5.176    | 0.437 | 4.32     | 6.032    |
| continue     | 1.66     | 0.219 | 1.23     | 2.09     |

Intervention means

| comparison       | estimate | se    | ci_lower | ci_upper |
|------------------|----------|-------|----------|----------|
| quit vs continue | 3.516    | 0.479 | 2.577    | 4.455    |

Contrasts

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
#>   estimator   comparison estimate        se ci_lower ci_upper
#> 1     gcomp quit vs cont 3.155727 0.4487520 2.276190 4.035265
#> 2       ipw quit vs cont 3.205240 0.4693563 2.285318 4.125162
#> 3  matching quit vs cont 2.984411 0.5091996 1.986398 3.982424
```

## Intervention types

Beyond static interventions, causatr supports modified treatment
policies (MTPs) and stochastic interventions:

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
```

**Estimator:** gcomp  ·  **Estimand:** ATE  ·  **Contrast:** difference
 ·  **CI method:** sandwich  ·  **N:** 1746

| intervention | estimate | se    | ci_lower | ci_upper |
|--------------|----------|-------|----------|----------|
| reduce10     | 2.621    | 0.24  | 2.15     | 3.092    |
| halved       | 2.619    | 0.243 | 2.142    | 3.095    |
| cap20        | 2.664    | 0.199 | 2.275    | 3.053    |
| observed     | 2.699    | 0.191 | 2.324    | 3.074    |

Intervention means

| comparison           | estimate | se    | ci_lower | ci_upper |
|----------------------|----------|-------|----------|----------|
| reduce10 vs observed | -0.078   | 0.163 | -0.398   | 0.242    |
| halved vs observed   | -0.08    | 0.168 | -0.409   | 0.249    |
| cap20 vs observed    | -0.035   | 0.073 | -0.178   | 0.108    |

Contrasts

## Diagnostics

Check covariate balance and positivity after fitting:

``` r
diag <- diagnose(fit_ipw)
diag          # positivity + balance summary
plot(diag)    # Love plot (requires cobalt)
```

## Features

- **Three estimation methods**: g-computation (parametric g-formula),
  IPW (self-contained density-ratio engine — no runtime dependency on
  WeightIt), and matching (via
  [MatchIt](https://kosukeimai.github.io/MatchIt/)). Matching is
  binary-only; continuous, categorical, count, and multivariate
  treatments use g-comp or IPW.
- **Longitudinal support**: ICE g-computation (Zivich et al. 2024) and
  longitudinal IPW for time-varying treatments. Sandwich variance via
  stacked estimating equations, plus parallel bootstrap via
  `boot::boot()` (with optional
  [future](https://future.futureverse.org/) backend).
- **Flexible interventions**: `static()`, `shift()`, `scale_by()`,
  `threshold()`, `dynamic()`, `ipsi()` (incremental propensity score),
  and `stochastic()` (user-defined randomised rules with Monte Carlo
  integration). Which interventions are available depends on the
  estimator — see the [interventions
  vignette](https://etverse.github.io/causatr/articles/interventions.html).
- **Treatment types**: binary, continuous, categorical (k \> 2), count
  (Poisson / negative binomial propensity via `propensity_family =`),
  and multivariate (joint) treatments. Multivariate IPW uses sequential
  MTP factorisation (Díaz et al. 2023) with optional stabilised weights
  (`stabilize = "marginal"`).
- **Any outcome family**: gaussian, binomial (logit / probit / cloglog),
  Poisson, quasibinomial (fractional), Gamma, negative binomial
  (`MASS::glm.nb`), beta regression (`betareg::betareg`), plus any
  family you pass through `model_fn`.
- **Pluggable models**: `stats::glm`, `mgcv::gam`, splines via `ns()` /
  `bs()`, or any fit function with signature
  `(formula, data, family, weights, ...)`. A two-tier numeric-variance
  fallback handles model classes without a `sandwich::estfun` method.
- **Robust inference**: analytic sandwich SE (default, via a unified
  influence-function engine) or nonparametric bootstrap with percentile
  CIs. Cluster-robust sandwich via `cluster =`; survey designs
  (`survey::svydesign`) auto-extract weights and PSU.
- **Built-in IPCW**: for MAR outcome censoring, `ipcw = TRUE` fits an
  internal censoring model and computes stabilised IPCW weights —
  provides doubly-robust protection under g-comp and is essential for
  IPW under MAR censoring. Custom censoring models via
  `censoring_model_fn =`.
- **Contrast types**: risk difference, risk ratio, odds ratio — ratio
  and OR use log-scale CIs.
- **Estimands**: ATE, ATT, ATC, or custom subgroups via `subset =` /
  `by =`.
- **Effect modification**: `by =` in `contrast()` for subgroup-specific
  effects. Under IPW and matching the modifier must be a baseline
  variable.
- **Built-in diagnostics**: positivity checks, covariate balance via
  [cobalt](https://ngreifer.github.io/cobalt/), weight summaries,
  censoring model diagnostics, Love plots.
- **Tidy integration**: `tidy()` / `glance()` / `confint()` / `coef()` /
  `vcov()` / `plot()` (forest plot via
  [forrest](https://github.com/etverse/forrest)) / broom-compatible
  output.

## References

Hernán MA, Robins JM (2025). *Causal Inference: What If*. Chapman &
Hall/CRC.

## Acknowledgements

This package was built with the contribution of
[Claude](https://claude.ai), Anthropic’s AI assistant.
