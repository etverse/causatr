# Introduction to causatr


causatr provides a unified interface for causal effect estimation via
three complementary methods:

- **G-computation** (parametric g-formula) — fits an outcome model
  *E*\[*Y*|*A*, *L*\] and standardises predictions over the target
  population.
- **Inverse probability weighting (IPW)** — reweights the data using
  propensity scores so that confounders are balanced across treatment
  groups.
- **Propensity score matching** — pairs treated and control units with
  similar propensity scores.

When all three methods agree, you can be more confident in your
findings. This is called **methodological triangulation**.

## Two-step API

All analyses follow the same two-step pattern:

``` r
# Step 1: Fit
fit <- causat(data, outcome = "Y", treatment = "A", confounders = ~ L1 + L2,
              method = "gcomp")

# Step 2: Contrast
result <- contrast(fit,
  interventions = list(treated = static(1), control = static(0)),
  reference = "control")
```

## Quick example

``` r
library(causatr)
data("nhefs")

fit <- causat(
  nhefs,
  outcome = "wt82_71",
  treatment = "qsmk",
  confounders = ~ sex + age + race + smokeintensity + smokeyrs +
    factor(exercise) + factor(active) + wt71,
  censoring = "censored"
)

result <- contrast(
  fit,
  interventions = list(quit = static(1), continue = static(0)),
  reference = "continue"
)
result
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
#> 1:         quit    5.033 0.3994    4.250    5.815
#> 2:     continue    1.877 0.1997    1.486    2.268
#> 
#> Contrasts:
#>          comparison estimate     se ci_lower ci_upper
#>              <char>    <num>  <num>    <num>    <num>
#> 1: quit vs continue    3.156 0.4488    2.276    4.035
```

## Intervention types

causatr supports several intervention types for g-computation:

<table>
<colgroup>
<col style="width: 33%" />
<col style="width: 33%" />
<col style="width: 33%" />
</colgroup>
<thead>
<tr>
<th>Function</th>
<th>Description</th>
<th>Example</th>
</tr>
</thead>
<tbody>
<tr>
<td><code>static(value)</code></td>
<td>Set treatment to a fixed value</td>
<td><code>static(1)</code></td>
</tr>
<tr>
<td><code>shift(delta)</code></td>
<td>Add a constant to observed treatment</td>
<td><code>shift(-10)</code></td>
</tr>
<tr>
<td><code>scale(factor)</code></td>
<td>Multiply observed treatment</td>
<td><code>scale(0.5)</code></td>
</tr>
<tr>
<td><code>threshold(lower, upper)</code></td>
<td>Clamp treatment within bounds</td>
<td><code>threshold(0, 20)</code></td>
</tr>
<tr>
<td><code>dynamic(rule)</code></td>
<td>Apply a user-defined rule</td>
<td><code>dynamic(\(d, a) ifelse(d$L &gt; 0, 1, 0))</code></td>
</tr>
<tr>
<td><code>NULL</code></td>
<td>Natural course (observed values)</td>
<td>—</td>
</tr>
</tbody>
</table>

IPW and matching currently support only `static()` interventions.

## Contrast types

<table>
<thead>
<tr>
<th>Type</th>
<th>Meaning</th>
</tr>
</thead>
<tbody>
<tr>
<td><code>"difference"</code></td>
<td>Risk/mean difference (default)</td>
</tr>
<tr>
<td><code>"ratio"</code></td>
<td>Risk ratio</td>
</tr>
<tr>
<td><code>"or"</code></td>
<td>Odds ratio (binary outcomes)</td>
</tr>
</tbody>
</table>

## Inference methods

<table>
<thead>
<tr>
<th>Method</th>
<th>Description</th>
</tr>
</thead>
<tbody>
<tr>
<td><code>"sandwich"</code></td>
<td>Robust sandwich SE (default, fast)</td>
</tr>
<tr>
<td><code>"bootstrap"</code></td>
<td>Nonparametric bootstrap (resamples full pipeline)</td>
</tr>
</tbody>
</table>

## Estimands

<table>
<thead>
<tr>
<th>Estimand</th>
<th>Population</th>
<th>Applicability</th>
</tr>
</thead>
<tbody>
<tr>
<td><code>"ATE"</code></td>
<td>All individuals</td>
<td>Always (default)</td>
</tr>
<tr>
<td><code>"ATT"</code></td>
<td>Observed treated (A = 1)</td>
<td>Binary point treatment only</td>
</tr>
<tr>
<td><code>"ATC"</code></td>
<td>Observed controls (A = 0)</td>
<td>Binary point treatment only</td>
</tr>
<tr>
<td><code>subset = quote(...)</code></td>
<td>Custom subgroup</td>
<td>Always</td>
</tr>
</tbody>
</table>

For g-computation, the estimand can be changed in `contrast()` without
refitting. For IPW and matching, refit with `causat(estimand = ...)`.

## Treatment types

causatr supports binary, continuous, categorical, and multivariate
treatments:

<table>
<thead>
<tr>
<th>Treatment type</th>
<th>G-comp</th>
<th>IPW</th>
<th>Matching</th>
</tr>
</thead>
<tbody>
<tr>
<td>Binary (0/1)</td>
<td>Yes</td>
<td>Yes</td>
<td>Yes</td>
</tr>
<tr>
<td>Continuous</td>
<td>Yes</td>
<td>Yes (GPS)</td>
<td>No</td>
</tr>
<tr>
<td>Categorical</td>
<td>Yes</td>
<td>Yes</td>
<td>Yes</td>
</tr>
<tr>
<td>Multivariate</td>
<td>Yes</td>
<td>Not yet</td>
<td>Not yet</td>
</tr>
</tbody>
</table>

For multivariate treatments (e.g. `treatment = c("A1", "A2")`), supply
interventions as named lists with one element per treatment variable.

## Effect modification

The `by` argument in `contrast()` stratifies estimates by levels of a
variable, producing subgroup-specific effects:

``` r
contrast(fit,
  interventions = list(quit = static(1), continue = static(0)),
  by = "sex"
)
```

## Diagnostics

Use `diagnose()` to check positivity and covariate balance after
fitting:

``` r
diag <- diagnose(fit)
diag          # print positivity + balance summary
plot(diag)    # Love plot (requires cobalt)
```

For IPW and matching, balance is computed via
[cobalt](https://ngreifer.github.io/cobalt/) and includes standardised
mean differences before and after adjustment.

## Extracting results

causatr result objects support standard R generics:

``` r
coef(result)       # named vector of E[Y^a]
confint(result)    # CI matrix (respects custom `level`)
vcov(result)       # vcov of marginal means
tidy(result)       # broom-style tidy data frame
glance(result)     # broom-style one-row summary
plot(result)       # forest plot (requires forrest package)
```

## Learning more

causatr includes detailed vignettes for each estimation method:
g-computation, IPW, matching, longitudinal ICE, and methodological
triangulation. Use `vignette(package = "causatr")` to see all available
vignettes.
