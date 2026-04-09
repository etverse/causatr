# Methodological triangulation: g-computation, IPW, and matching


When multiple estimation methods agree on the same causal effect, we
gain confidence that the result is not an artefact of a single modelling
assumption. This vignette compares all three methods on the NHEFS
dataset.

## Setup

``` r
library(causatr)
library(tinyplot)
data("nhefs")

nhefs_complete <- nhefs[complete.cases(nhefs), ]
nhefs_complete$gained_weight <- as.integer(nhefs_complete$wt82_71 > 0)
nhefs_complete$sex <- factor(nhefs_complete$sex, levels = 0:1, labels = c("Male", "Female"))

confounders <- ~ sex + age + I(age^2) + race + factor(education) +
  smokeintensity + I(smokeintensity^2) + smokeyrs + I(smokeyrs^2) +
  factor(exercise) + factor(active) + wt71 + I(wt71^2)
```

## Continuous outcome: weight change

### Fitting all three methods

``` r
fit_gc <- causat(nhefs, outcome = "wt82_71", treatment = "qsmk",
  confounders = confounders, censoring = "censored")

fit_ipw <- causat(nhefs_complete, outcome = "wt82_71", treatment = "qsmk",
  confounders = confounders, method = "ipw")

fit_m <- causat(nhefs_complete, outcome = "wt82_71", treatment = "qsmk",
  confounders = confounders, method = "matching", estimand = "ATT")
```

### Contrasts

``` r
ivs <- list(quit = static(1), continue = static(0))

res_gc <- contrast(fit_gc, ivs, reference = "continue")
res_ipw <- contrast(fit_ipw, ivs, reference = "continue")
res_m <- contrast(fit_m, ivs, reference = "continue")

results <- rbind(
  cbind(method = "G-computation", res_gc$contrasts),
  cbind(method = "IPW", res_ipw$contrasts),
  cbind(method = "Matching (ATT)", res_m$contrasts)
)
results
#>            method       comparison estimate        se ci_lower ci_upper
#>            <char>           <char>    <num>     <num>    <num>    <num>
#> 1:  G-computation quit vs continue 3.462484 0.4677686 2.545674 4.379294
#> 2:            IPW quit vs continue 3.499174 0.4902045 2.538391 4.459958
#> 3: Matching (ATT) quit vs continue 3.341010 0.5593230 2.244757 4.437263
```

### Forest plot

``` r
tinyplot(
  estimate ~ method,
  data = results,
  type = "pointrange",
  ymin = results$ci_lower,
  ymax = results$ci_upper,
  xlab = "Method",
  ylab = "ATE: weight change (kg)",
  main = "Triangulation: continuous outcome"
)
abline(h = 0, lty = 2, col = "grey40")
```

<img
src="vignettes/triangulation.markdown_strict_files/figure-markdown_strict/unnamed-chunk-5-1.png"
data-fig-alt="Forest plot comparing g-computation, IPW, and matching estimates of the effect of quitting smoking on weight change." />

## Binary outcome: gained weight

``` r
fit_gc_bin <- causat(nhefs_complete, outcome = "gained_weight",
  treatment = "qsmk", confounders = confounders, family = "binomial")

fit_ipw_bin <- causat(nhefs_complete, outcome = "gained_weight",
  treatment = "qsmk", confounders = confounders, method = "ipw")

fit_m_bin <- causat(nhefs_complete, outcome = "gained_weight",
  treatment = "qsmk", confounders = confounders,
  method = "matching", estimand = "ATT")

res_gc_bin <- contrast(fit_gc_bin, ivs, reference = "continue")
res_ipw_bin <- contrast(fit_ipw_bin, ivs, reference = "continue")
res_m_bin <- contrast(fit_m_bin, ivs, reference = "continue")

results_bin <- rbind(
  cbind(method = "G-computation", res_gc_bin$contrasts),
  cbind(method = "IPW", res_ipw_bin$contrasts),
  cbind(method = "Matching (ATT)", res_m_bin$contrasts)
)
results_bin
#>            method       comparison  estimate         se   ci_lower  ci_upper
#>            <char>           <char>     <num>      <num>      <num>     <num>
#> 1:  G-computation quit vs continue 0.1326922 0.02367089 0.08629815 0.1790863
#> 2:            IPW quit vs continue 0.1379464 0.02498002 0.08898648 0.1869064
#> 3: Matching (ATT) quit vs continue 0.1488834 0.03172616 0.08670125 0.2110655
```

``` r
tinyplot(
  estimate ~ method,
  data = results_bin,
  type = "pointrange",
  ymin = results_bin$ci_lower,
  ymax = results_bin$ci_upper,
  xlab = "Method",
  ylab = "Risk difference (gained weight)",
  main = "Triangulation: binary outcome"
)
abline(h = 0, lty = 2, col = "grey40")
```

<img
src="vignettes/triangulation.markdown_strict_files/figure-markdown_strict/unnamed-chunk-7-1.png"
data-fig-alt="Forest plot comparing methods on the binary outcome (gained weight)." />

## Diagnostics

Before interpreting results, check that the assumptions hold.
`diagnose()` provides positivity checks, covariate balance, and
method-specific summaries.

``` r
diag_ipw <- diagnose(fit_ipw)
diag_m <- diagnose(fit_m)
```

### IPW balance

``` r
plot(diag_ipw)
```

<img
src="vignettes/triangulation.markdown_strict_files/figure-markdown_strict/unnamed-chunk-9-1.png"
data-fig-alt="Love plot for IPW covariate balance in the triangulation analysis." />

### Matching balance

``` r
plot(diag_m)
```

<img
src="vignettes/triangulation.markdown_strict_files/figure-markdown_strict/unnamed-chunk-10-1.png"
data-fig-alt="Love plot for matching covariate balance in the triangulation analysis." />

## When to use which method

<table>
<colgroup>
<col style="width: 25%" />
<col style="width: 25%" />
<col style="width: 25%" />
<col style="width: 25%" />
</colgroup>
<thead>
<tr>
<th>Criterion</th>
<th>G-computation</th>
<th>IPW</th>
<th>Matching</th>
</tr>
</thead>
<tbody>
<tr>
<td>Model specification</td>
<td>Outcome model</td>
<td>Treatment model</td>
<td>Treatment model</td>
</tr>
<tr>
<td>Non-static interventions</td>
<td>Yes</td>
<td>No</td>
<td>No</td>
</tr>
<tr>
<td>Continuous treatment</td>
<td>Yes</td>
<td>Limited</td>
<td>Limited</td>
</tr>
<tr>
<td>Variance estimation</td>
<td>Sandwich or bootstrap</td>
<td>M-estimation sandwich</td>
<td>Cluster-robust sandwich</td>
</tr>
<tr>
<td>Robustness to</td>
<td>Treatment model misspecification</td>
<td>Outcome model misspecification</td>
<td>Outcome model misspecification</td>
</tr>
</tbody>
</table>

When the outcome model and treatment model are both correctly specified,
all three methods should produce similar estimates. Disagreement
suggests model misspecification — investigate further.

## References

Hernán MA, Robins JM (2025). *Causal Inference: What If*. Chapman &
Hall/CRC.
