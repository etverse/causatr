# G-computation with causatr


G-computation (the parametric g-formula) estimates causal effects by
fitting a model for the outcome conditional on treatment and
confounders, then standardising predictions over the target population
under each hypothetical intervention. This vignette demonstrates
g-computation in causatr using the NHEFS dataset from Hernán & Robins
(2025), covering every supported combination of treatment type, outcome
type, contrast scale, and inference method for time-fixed treatments.

## Setup

``` r
library(causatr)
library(tinyplot)

data("nhefs")

nhefs_complete <- nhefs[!is.na(nhefs$wt82_71), ]

nhefs_complete$gained_weight <- as.integer(nhefs_complete$wt82_71 > 0)
```

We create a binary outcome `gained_weight` (1 if weight increased, 0
otherwise) for the binary outcome examples below.

## Binary treatment, continuous outcome

This is the core example from Chapter 13 of Hernán & Robins: the average
causal effect of quitting smoking (`qsmk`) on weight change (`wt82_71`).

### ATE with sandwich SE

``` r
fit_gc <- causat(
  nhefs,
  outcome = "wt82_71",
  treatment = "qsmk",
  confounders = ~ sex + age + I(age^2) + race + factor(education) +
    smokeintensity + I(smokeintensity^2) + smokeyrs + I(smokeyrs^2) +
    factor(exercise) + factor(active) + wt71 + I(wt71^2) +
    qsmk:smokeintensity,
  censoring = "censored"
)

res_ate_sw <- contrast(
  fit_gc,
  interventions = list(quit = static(1), continue = static(0)),
  reference = "continue",
  type = "difference",
  ci_method = "sandwich"
)
res_ate_sw
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

The book reports ATE ≈ 3.5 kg (95% CI: 2.6, 4.5).

### ATE with bootstrap SE

``` r
res_ate_bs <- contrast(
  fit_gc,
  interventions = list(quit = static(1), continue = static(0)),
  reference = "continue",
  type = "difference",
  ci_method = "bootstrap",
  n_boot = 200L
)
res_ate_bs
#> <causatr_result>
#>  Method:      gcomp
#>  Contrast:    difference
#>  CI method:   bootstrap
#>  N:           1629
#> 
#> Intervention means:
#>    intervention estimate     se ci_lower ci_upper
#>          <char>    <num>  <num>    <num>    <num>
#> 1:         quit    5.176 0.4238    4.345    6.006
#> 2:     continue    1.660 0.2058    1.257    2.064
#> 
#> Contrasts:
#>          comparison estimate     se ci_lower ci_upper
#>              <char>    <num>  <num>    <num>    <num>
#> 1: quit vs continue    3.516 0.4721     2.59    4.441
```

Sandwich and bootstrap SEs should be in close agreement for correctly
specified GLMs.

### ATT estimand

The average treatment effect on the treated (ATT) averages only over
individuals who actually quit smoking. With g-computation, the estimand
can be changed in `contrast()` without refitting.

``` r
res_att <- contrast(
  fit_gc,
  interventions = list(quit = static(1), continue = static(0)),
  reference = "continue",
  estimand = "ATT",
  ci_method = "sandwich"
)
res_att
#> <causatr_result>
#>  Method:      gcomp
#>  Contrast:    difference
#>  CI method:   sandwich
#>  N:           428
#> 
#> Intervention means:
#>    intervention estimate     se ci_lower ci_upper
#>          <char>    <num>  <num>    <num>    <num>
#> 1:         quit   4.3662 0.4042   3.5740    5.158
#> 2:     continue   0.9305 0.2362   0.4675    1.394
#> 
#> Contrasts:
#>          comparison estimate     se ci_lower ci_upper
#>              <char>    <num>  <num>    <num>    <num>
#> 1: quit vs continue    3.436 0.4639    2.526    4.345
```

### ATC estimand

The average treatment effect on the controls (ATC) averages over those
who continued smoking.

``` r
res_atc <- contrast(
  fit_gc,
  interventions = list(quit = static(1), continue = static(0)),
  reference = "continue",
  estimand = "ATC",
  ci_method = "sandwich"
)
res_atc
#> <causatr_result>
#>  Method:      gcomp
#>  Contrast:    difference
#>  CI method:   sandwich
#>  N:           1201
#> 
#> Intervention means:
#>    intervention estimate     se ci_lower ci_upper
#>          <char>    <num>  <num>    <num>    <num>
#> 1:         quit    5.465 0.4442    4.594    6.335
#> 2:     continue    1.920 0.2040    1.520    2.320
#> 
#> Contrasts:
#>          comparison estimate     se ci_lower ci_upper
#>              <char>    <num>  <num>    <num>    <num>
#> 1: quit vs continue    3.544 0.4869     2.59    4.498
```

### Subset estimand

A custom subgroup: effect among individuals aged 50 or older.

``` r
res_sub <- contrast(
  fit_gc,
  interventions = list(quit = static(1), continue = static(0)),
  reference = "continue",
  subset = quote(age >= 50),
  ci_method = "sandwich"
)
res_sub
#> <causatr_result>
#>  Method:      gcomp
#>  Contrast:    difference
#>  CI method:   sandwich
#>  N:           565
#> 
#> Intervention means:
#>    intervention estimate     se ci_lower ci_upper
#>          <char>    <num>  <num>    <num>    <num>
#> 1:         quit   2.5958 0.4527    1.708   3.4831
#> 2:     continue  -0.8739 0.2971   -1.456  -0.2917
#> 
#> Contrasts:
#>          comparison estimate     se ci_lower ci_upper
#>              <char>    <num>  <num>    <num>    <num>
#> 1: quit vs continue     3.47 0.4691     2.55    4.389
```

### Extracting results programmatically

``` r
coef(res_ate_sw)
#>     quit continue 
#> 5.175939 1.660268
confint(res_ate_sw)
#>             lower    upper
#> quit     4.333495 6.018383
#> continue 1.255834 2.064701
```

## Binary treatment, binary outcome

Using the `gained_weight` indicator as a binary outcome, we estimate the
risk difference, risk ratio, and odds ratio of quitting smoking on
gaining weight.

### Risk difference (sandwich)

``` r
fit_bin <- causat(
  nhefs_complete,
  outcome = "gained_weight",
  treatment = "qsmk",
  confounders = ~ sex + age + I(age^2) + race + factor(education) +
    smokeintensity + I(smokeintensity^2) + smokeyrs + I(smokeyrs^2) +
    factor(exercise) + factor(active) + wt71 + I(wt71^2) +
    qsmk:smokeintensity,
  family = "binomial"
)

res_rd <- contrast(
  fit_bin,
  interventions = list(quit = static(1), continue = static(0)),
  reference = "continue",
  type = "difference",
  ci_method = "sandwich"
)
res_rd
#> <causatr_result>
#>  Method:      gcomp
#>  Contrast:    difference
#>  CI method:   sandwich
#>  N:           1566
#> 
#> Intervention means:
#>    intervention estimate      se ci_lower ci_upper
#>          <char>    <num>   <num>    <num>    <num>
#> 1:         quit   0.7696 0.01964   0.7311   0.8081
#> 2:     continue   0.6383 0.01330   0.6122   0.6644
#> 
#> Contrasts:
#>          comparison estimate      se ci_lower ci_upper
#>              <char>    <num>   <num>    <num>    <num>
#> 1: quit vs continue   0.1313 0.02391  0.08444   0.1782
```

### Risk difference (bootstrap)

``` r
res_rd_bs <- contrast(
  fit_bin,
  interventions = list(quit = static(1), continue = static(0)),
  reference = "continue",
  type = "difference",
  ci_method = "bootstrap",
  n_boot = 200L
)
res_rd_bs
#> <causatr_result>
#>  Method:      gcomp
#>  Contrast:    difference
#>  CI method:   bootstrap
#>  N:           1566
#> 
#> Intervention means:
#>    intervention estimate      se ci_lower ci_upper
#>          <char>    <num>   <num>    <num>    <num>
#> 1:         quit   0.7696 0.02033   0.7298   0.8094
#> 2:     continue   0.6383 0.01339   0.6120   0.6645
#> 
#> Contrasts:
#>          comparison estimate      se ci_lower ci_upper
#>              <char>    <num>   <num>    <num>    <num>
#> 1: quit vs continue   0.1313 0.02337  0.08551   0.1771
```

### Risk ratio

``` r
res_rr <- contrast(
  fit_bin,
  interventions = list(quit = static(1), continue = static(0)),
  reference = "continue",
  type = "ratio",
  ci_method = "sandwich"
)
res_rr
#> <causatr_result>
#>  Method:      gcomp
#>  Contrast:    ratio
#>  CI method:   sandwich
#>  N:           1566
#> 
#> Intervention means:
#>    intervention estimate      se ci_lower ci_upper
#>          <char>    <num>   <num>    <num>    <num>
#> 1:         quit   0.7696 0.01964   0.7311   0.8081
#> 2:     continue   0.6383 0.01330   0.6122   0.6644
#> 
#> Contrasts:
#>          comparison estimate      se ci_lower ci_upper
#>              <char>    <num>   <num>    <num>    <num>
#> 1: quit vs continue    1.206 0.04007    1.127    1.284
```

### Odds ratio

``` r
res_or <- contrast(
  fit_bin,
  interventions = list(quit = static(1), continue = static(0)),
  reference = "continue",
  type = "or",
  ci_method = "sandwich"
)
res_or
#> <causatr_result>
#>  Method:      gcomp
#>  Contrast:    or
#>  CI method:   sandwich
#>  N:           1566
#> 
#> Intervention means:
#>    intervention estimate      se ci_lower ci_upper
#>          <char>    <num>   <num>    <num>    <num>
#> 1:         quit   0.7696 0.01964   0.7311   0.8081
#> 2:     continue   0.6383 0.01330   0.6122   0.6644
#> 
#> Contrasts:
#>          comparison estimate    se ci_lower ci_upper
#>              <char>    <num> <num>    <num>    <num>
#> 1: quit vs continue    1.893 0.238    1.426    2.359
```

## Continuous treatment

G-computation also supports continuous treatments with modified
treatment policy interventions. Here we use `smokeintensity` (cigarettes
per day) as the treatment and `wt82_71` as the outcome.

``` r
fit_cont <- causat(
  nhefs,
  outcome = "wt82_71",
  treatment = "smokeintensity",
  confounders = ~ sex + age + I(age^2) + race + factor(education) +
    smokeyrs + I(smokeyrs^2) + factor(exercise) + factor(active) +
    wt71 + I(wt71^2),
  censoring = "censored"
)
```

### Shift intervention

Reduce smoking intensity by 10 cigarettes per day for everyone.

``` r
res_shift <- contrast(
  fit_cont,
  interventions = list(reduce10 = shift(-10), observed = NULL),
  reference = "observed",
  type = "difference",
  ci_method = "sandwich"
)
res_shift
#> <causatr_result>
#>  Method:      gcomp
#>  Contrast:    difference
#>  CI method:   sandwich
#>  N:           1629
#> 
#> Intervention means:
#>    intervention estimate     se ci_lower ci_upper
#>          <char>    <num>  <num>    <num>    <num>
#> 1:     reduce10    2.594 0.2479    2.108    3.080
#> 2:     observed    2.548 0.1881    2.179    2.917
#> 
#> Contrasts:
#>              comparison estimate     se ci_lower ci_upper
#>                  <char>    <num>  <num>    <num>    <num>
#> 1: reduce10 vs observed  0.04617 0.1807   -0.308   0.4003
```

### Scale intervention

Halve each individual’s smoking intensity.

``` r
res_scale <- contrast(
  fit_cont,
  interventions = list(halved = scale(0.5), observed = NULL),
  reference = "observed",
  type = "difference",
  ci_method = "sandwich"
)
res_scale
#> <causatr_result>
#>  Method:      gcomp
#>  Contrast:    difference
#>  CI method:   sandwich
#>  N:           1629
#> 
#> Intervention means:
#>    intervention estimate     se ci_lower ci_upper
#>          <char>    <num>  <num>    <num>    <num>
#> 1:       halved    2.595 0.2512    2.103    3.088
#> 2:     observed    2.548 0.1881    2.179    2.917
#> 
#> Contrasts:
#>            comparison estimate     se ci_lower ci_upper
#>                <char>    <num>  <num>    <num>    <num>
#> 1: halved vs observed  0.04744 0.1857  -0.3165   0.4114
```

### Threshold intervention

Cap smoking intensity at 20 cigarettes per day for everyone.

``` r
res_thresh <- contrast(
  fit_cont,
  interventions = list(cap20 = threshold(0, 20), observed = NULL),
  reference = "observed",
  type = "difference",
  ci_method = "sandwich"
)
res_thresh
#> <causatr_result>
#>  Method:      gcomp
#>  Contrast:    difference
#>  CI method:   sandwich
#>  N:           1629
#> 
#> Intervention means:
#>    intervention estimate     se ci_lower ci_upper
#>          <char>    <num>  <num>    <num>    <num>
#> 1:        cap20    2.568 0.1974    2.181    2.955
#> 2:     observed    2.548 0.1881    2.179    2.917
#> 
#> Contrasts:
#>           comparison estimate     se ci_lower ci_upper
#>               <char>    <num>  <num>    <num>    <num>
#> 1: cap20 vs observed  0.02064 0.0808  -0.1377    0.179
```

### Comparing multiple interventions

``` r
res_multi <- contrast(
  fit_cont,
  interventions = list(
    reduce10 = shift(-10),
    halved   = scale(0.5),
    cap20    = threshold(0, 20),
    observed = NULL
  ),
  reference = "observed",
  type = "difference",
  ci_method = "sandwich"
)
res_multi
#> <causatr_result>
#>  Method:      gcomp
#>  Contrast:    difference
#>  CI method:   sandwich
#>  N:           1629
#> 
#> Intervention means:
#>    intervention estimate     se ci_lower ci_upper
#>          <char>    <num>  <num>    <num>    <num>
#> 1:     reduce10    2.594 0.2479    2.108    3.080
#> 2:       halved    2.595 0.2512    2.103    3.088
#> 3:        cap20    2.568 0.1974    2.181    2.955
#> 4:     observed    2.548 0.1881    2.179    2.917
#> 
#> Contrasts:
#>              comparison estimate     se ci_lower ci_upper
#>                  <char>    <num>  <num>    <num>    <num>
#> 1: reduce10 vs observed  0.04617 0.1807  -0.3080   0.4003
#> 2:   halved vs observed  0.04744 0.1857  -0.3165   0.4114
#> 3:    cap20 vs observed  0.02064 0.0808  -0.1377   0.1790
```

### Visualising intervention effects

``` r
est <- res_multi$contrasts
tinyplot(
  estimate ~ comparison,
  data = est,
  type = "pointrange",
  ymin = est$ci_lower,
  ymax = est$ci_upper,
  xlab = "Comparison",
  ylab = "Difference in weight change (kg)",
  main = "Continuous treatment: modified treatment policies"
)
abline(h = 0, lty = 2, col = "grey40")
```

<img
src="vignettes/gcomp.markdown_strict_files/figure-markdown_strict/unnamed-chunk-18-1.png"
data-fig-alt="Forest plot of mean weight change under different smoking intensity interventions." />

### Continuous treatment with bootstrap

``` r
res_shift_bs <- contrast(
  fit_cont,
  interventions = list(reduce10 = shift(-10), observed = NULL),
  reference = "observed",
  type = "difference",
  ci_method = "bootstrap",
  n_boot = 200L
)
res_shift_bs
#> <causatr_result>
#>  Method:      gcomp
#>  Contrast:    difference
#>  CI method:   bootstrap
#>  N:           1629
#> 
#> Intervention means:
#>    intervention estimate     se ci_lower ci_upper
#>          <char>    <num>  <num>    <num>    <num>
#> 1:     reduce10    2.594 0.2551    2.094    3.094
#> 2:     observed    2.548 0.1949    2.166    2.930
#> 
#> Contrasts:
#>              comparison estimate     se ci_lower ci_upper
#>                  <char>    <num>  <num>    <num>    <num>
#> 1: reduce10 vs observed  0.04617 0.1879  -0.3222   0.4145
```

## Dynamic intervention

A dynamic intervention assigns treatment based on individual
characteristics. Here, we assign quitting smoking only to individuals
who smoked more than 20 cigarettes per day at baseline.

``` r
res_dyn <- contrast(
  fit_gc,
  interventions = list(
    rule = dynamic(\(data, trt) ifelse(data$smokeintensity > 20, 1, 0)),
    all_quit = static(1)
  ),
  reference = "all_quit",
  type = "difference",
  ci_method = "sandwich"
)
res_dyn
#> <causatr_result>
#>  Method:      gcomp
#>  Contrast:    difference
#>  CI method:   sandwich
#>  N:           1629
#> 
#> Intervention means:
#>    intervention estimate     se ci_lower ci_upper
#>          <char>    <num>  <num>    <num>    <num>
#> 1:         rule    2.852 0.2898    2.284    3.420
#> 2:     all_quit    5.176 0.4298    4.333    6.018
#> 
#> Contrasts:
#>          comparison estimate     se ci_lower ci_upper
#>              <char>    <num>  <num>    <num>    <num>
#> 1: rule vs all_quit   -2.324 0.3388   -2.988    -1.66
```

## GAM model via model_fn

Pass `mgcv::gam` instead of `stats::glm` for flexible nonlinear
confounder adjustment using splines.

``` r
fit_gam <- causat(
  nhefs,
  outcome = "wt82_71",
  treatment = "qsmk",
  confounders = ~ sex + s(age) + race + factor(education) +
    s(smokeintensity) + s(smokeyrs) + factor(exercise) +
    factor(active) + s(wt71),
  censoring = "censored",
  model_fn = mgcv::gam
)

res_gam <- contrast(
  fit_gam,
  interventions = list(quit = static(1), continue = static(0)),
  reference = "continue",
  ci_method = "sandwich"
)
res_gam
#> <causatr_result>
#>  Method:      gcomp
#>  Contrast:    difference
#>  CI method:   sandwich
#>  N:           1629
#> 
#> Intervention means:
#>    intervention estimate       se ci_lower ci_upper
#>          <char>    <num>    <num>    <num>    <num>
#> 1:         quit    5.093 0.007859    5.078    5.109
#> 2:     continue    1.670 0.003917    1.662    1.678
#> 
#> Contrasts:
#>          comparison estimate       se ci_lower ci_upper
#>              <char>    <num>    <num>    <num>    <num>
#> 1: quit vs continue    3.424 0.008853    3.406    3.441
```

## Summary of covered combinations

<table>
<thead>
<tr>
<th>Treatment</th>
<th>Outcome</th>
<th>Contrast</th>
<th>Inference</th>
<th>Estimand</th>
<th>Intervention</th>
</tr>
</thead>
<tbody>
<tr>
<td>Binary</td>
<td>Continuous</td>
<td>Difference</td>
<td>Sandwich</td>
<td>ATE</td>
<td>Static</td>
</tr>
<tr>
<td>Binary</td>
<td>Continuous</td>
<td>Difference</td>
<td>Bootstrap</td>
<td>ATE</td>
<td>Static</td>
</tr>
<tr>
<td>Binary</td>
<td>Continuous</td>
<td>Difference</td>
<td>Sandwich</td>
<td>ATT</td>
<td>Static</td>
</tr>
<tr>
<td>Binary</td>
<td>Continuous</td>
<td>Difference</td>
<td>Sandwich</td>
<td>ATC</td>
<td>Static</td>
</tr>
<tr>
<td>Binary</td>
<td>Continuous</td>
<td>Difference</td>
<td>Sandwich</td>
<td>Subset</td>
<td>Static</td>
</tr>
<tr>
<td>Binary</td>
<td>Binary</td>
<td>Difference</td>
<td>Sandwich</td>
<td>ATE</td>
<td>Static</td>
</tr>
<tr>
<td>Binary</td>
<td>Binary</td>
<td>Difference</td>
<td>Bootstrap</td>
<td>ATE</td>
<td>Static</td>
</tr>
<tr>
<td>Binary</td>
<td>Binary</td>
<td>Ratio</td>
<td>Sandwich</td>
<td>ATE</td>
<td>Static</td>
</tr>
<tr>
<td>Binary</td>
<td>Binary</td>
<td>OR</td>
<td>Sandwich</td>
<td>ATE</td>
<td>Static</td>
</tr>
<tr>
<td>Continuous</td>
<td>Continuous</td>
<td>Difference</td>
<td>Sandwich</td>
<td>ATE</td>
<td>Shift</td>
</tr>
<tr>
<td>Continuous</td>
<td>Continuous</td>
<td>Difference</td>
<td>Sandwich</td>
<td>ATE</td>
<td>Scale</td>
</tr>
<tr>
<td>Continuous</td>
<td>Continuous</td>
<td>Difference</td>
<td>Sandwich</td>
<td>ATE</td>
<td>Threshold</td>
</tr>
<tr>
<td>Continuous</td>
<td>Continuous</td>
<td>Difference</td>
<td>Bootstrap</td>
<td>ATE</td>
<td>Shift</td>
</tr>
<tr>
<td>Binary</td>
<td>Continuous</td>
<td>Difference</td>
<td>Sandwich</td>
<td>ATE</td>
<td>Dynamic</td>
</tr>
<tr>
<td>Binary</td>
<td>Continuous</td>
<td>Difference</td>
<td>Sandwich</td>
<td>ATE</td>
<td>GAM model</td>
</tr>
</tbody>
</table>

## References

Hernán MA, Robins JM (2025). *Causal Inference: What If*. Chapman &
Hall/CRC. Chapter 13: Standardization and the parametric g-formula.
