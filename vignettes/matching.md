# Propensity score matching with causatr


Propensity score matching estimates causal effects by pairing treated
and control individuals with similar propensity scores, then comparing
outcomes within matched sets. causatr delegates match construction to
[MatchIt](https://kosukeimai.github.io/MatchIt/) and fits the outcome
model on the matched sample with cluster-robust sandwich SEs (via
`sandwich::vcovCL()` on the matched-pair subclass).

This vignette demonstrates propensity score matching in causatr using
the NHEFS dataset from Hernán & Robins (2025), covering every supported
combination of treatment type, outcome type, contrast scale, and
inference method for time-fixed treatments.

**Note:** Matching currently supports only `static()` interventions and
binary treatments. For modified treatment policies, use g-computation
(see `vignette("gcomp")`).

## Setup

``` r
library(causatr)
library(tinyplot)

data("nhefs")

nhefs_complete <- nhefs[complete.cases(nhefs), ]

nhefs_complete$gained_weight <- as.integer(nhefs_complete$wt82_71 > 0)
```

## Binary treatment, continuous outcome

Matching estimate of the effect of quitting smoking on weight change,
following the approach in Chapter 15 of Hernán & Robins.

### ATT with sandwich SE

Matching naturally targets the ATT (effect on the treated). This is the
default for nearest-neighbour matching.

``` r
fit_m <- causat(
  nhefs_complete,
  outcome = "wt82_71",
  treatment = "qsmk",
  confounders = ~ sex + age + I(age^2) + race + factor(education) +
    smokeintensity + I(smokeintensity^2) + smokeyrs + I(smokeyrs^2) +
    factor(exercise) + factor(active) + wt71 + I(wt71^2),
  method = "matching",
  estimand = "ATT"
)

res_att_sw <- contrast(
  fit_m,
  interventions = list(quit = static(1), continue = static(0)),
  reference = "continue",
  type = "difference",
  ci_method = "sandwich"
)
res_att_sw
#> <causatr_result>
#>  Method:      matching
#>  Contrast:    difference
#>  CI method:   sandwich
#>  N:           403
#> 
#> Intervention means:
#>    intervention estimate     se ci_lower ci_upper
#>          <char>    <num>  <num>    <num>    <num>
#> 1:         quit    4.525 0.4358   3.6710    5.379
#> 2:     continue    1.184 0.3837   0.4319    1.936
#> 
#> Contrasts:
#>          comparison estimate     se ci_lower ci_upper
#>              <char>    <num>  <num>    <num>    <num>
#> 1: quit vs continue    3.341 0.5593    2.245    4.437
```

### ATT with bootstrap SE

The bootstrap resamples individuals, re-matches, and refits the outcome
model on each bootstrap sample, fully accounting for matching
uncertainty.

``` r
res_att_bs <- contrast(
  fit_m,
  interventions = list(quit = static(1), continue = static(0)),
  reference = "continue",
  type = "difference",
  ci_method = "bootstrap",
  n_boot = 200L
)
res_att_bs
#> <causatr_result>
#>  Method:      matching
#>  Contrast:    difference
#>  CI method:   bootstrap
#>  N:           403
#> 
#> Intervention means:
#>    intervention estimate     se ci_lower ci_upper
#>          <char>    <num>  <num>    <num>    <num>
#> 1:         quit    4.525 0.4337   3.6750    5.375
#> 2:     continue    1.184 0.3864   0.4268    1.941
#> 
#> Contrasts:
#>          comparison estimate     se ci_lower ci_upper
#>              <char>    <num>  <num>    <num>    <num>
#> 1: quit vs continue    3.341 0.5652    2.233    4.449
```

### ATE estimand

For the ATE, MatchIt performs full matching (each unit gets a weight)
rather than 1:1 nearest-neighbour matching.

``` r
fit_m_ate <- causat(
  nhefs_complete,
  outcome = "wt82_71",
  treatment = "qsmk",
  confounders = ~ sex + age + I(age^2) + race + factor(education) +
    smokeintensity + I(smokeintensity^2) + smokeyrs + I(smokeyrs^2) +
    factor(exercise) + factor(active) + wt71 + I(wt71^2),
  method = "matching",
  estimand = "ATE"
)

res_ate_sw <- contrast(
  fit_m_ate,
  interventions = list(quit = static(1), continue = static(0)),
  reference = "continue",
  type = "difference",
  ci_method = "sandwich"
)
res_ate_sw
#> <causatr_result>
#>  Method:      matching
#>  Contrast:    difference
#>  CI method:   sandwich
#>  N:           1566
#> 
#> Intervention means:
#>    intervention estimate     se ci_lower ci_upper
#>          <char>    <num>  <num>    <num>    <num>
#> 1:         quit    5.442 0.5504    4.363    6.520
#> 2:     continue    1.831 0.2411    1.358    2.303
#> 
#> Contrasts:
#>          comparison estimate     se ci_lower ci_upper
#>              <char>    <num>  <num>    <num>    <num>
#> 1: quit vs continue    3.611 0.5785    2.477    4.745
```

### ATE with bootstrap SE

``` r
res_ate_bs <- contrast(
  fit_m_ate,
  interventions = list(quit = static(1), continue = static(0)),
  reference = "continue",
  type = "difference",
  ci_method = "bootstrap",
  n_boot = 200L
)
res_ate_bs
#> <causatr_result>
#>  Method:      matching
#>  Contrast:    difference
#>  CI method:   bootstrap
#>  N:           1566
#> 
#> Intervention means:
#>    intervention estimate     se ci_lower ci_upper
#>          <char>    <num>  <num>    <num>    <num>
#> 1:         quit    5.442 0.5741    4.316    6.567
#> 2:     continue    1.831 0.2334    1.373    2.288
#> 
#> Contrasts:
#>          comparison estimate     se ci_lower ci_upper
#>              <char>    <num>  <num>    <num>    <num>
#> 1: quit vs continue    3.611 0.6168    2.402     4.82
```

### ATC estimand

The ATC targets the effect among controls (those who continued smoking).

``` r
fit_m_atc <- causat(
  nhefs_complete,
  outcome = "wt82_71",
  treatment = "qsmk",
  confounders = ~ sex + age + I(age^2) + race + factor(education) +
    smokeintensity + I(smokeintensity^2) + smokeyrs + I(smokeyrs^2) +
    factor(exercise) + factor(active) + wt71 + I(wt71^2),
  method = "matching",
  estimand = "ATC"
)
#> Warning: Fewer treated units than control units; not all control units will get
#> a match.

res_atc <- contrast(
  fit_m_atc,
  interventions = list(quit = static(1), continue = static(0)),
  reference = "continue",
  type = "difference",
  ci_method = "sandwich"
)
res_atc
#> <causatr_result>
#>  Method:      matching
#>  Contrast:    difference
#>  CI method:   sandwich
#>  N:           1163
#> 
#> Intervention means:
#>    intervention estimate     se ci_lower ci_upper
#>          <char>    <num>  <num>    <num>    <num>
#> 1:         quit    4.525 0.4358    3.671    5.379
#> 2:     continue    3.203 0.3556    2.506    3.900
#> 
#> Contrasts:
#>          comparison estimate     se ci_lower ci_upper
#>              <char>    <num>  <num>    <num>    <num>
#> 1: quit vs continue    1.322 0.5797   0.1857    2.458
```

### Extracting results programmatically

``` r
coef(res_att_sw)
#>     quit continue 
#> 4.525079 1.184069
confint(res_att_sw)
#>              lower    upper
#> quit     3.6709621 5.379196
#> continue 0.4319368 1.936201
```

## Binary treatment, binary outcome

Using the `gained_weight` indicator as a binary outcome. The outcome
model on the matched sample is still `Y ~ A` (linear probability model),
so `contrast()` computes marginal risks from the weighted predictions.

### Risk difference (sandwich)

``` r
fit_m_bin <- causat(
  nhefs_complete,
  outcome = "gained_weight",
  treatment = "qsmk",
  confounders = ~ sex + age + I(age^2) + race + factor(education) +
    smokeintensity + I(smokeintensity^2) + smokeyrs + I(smokeyrs^2) +
    factor(exercise) + factor(active) + wt71 + I(wt71^2),
  method = "matching",
  estimand = "ATT"
)

res_rd <- contrast(
  fit_m_bin,
  interventions = list(quit = static(1), continue = static(0)),
  reference = "continue",
  type = "difference",
  ci_method = "sandwich"
)
res_rd
#> <causatr_result>
#>  Method:      matching
#>  Contrast:    difference
#>  CI method:   sandwich
#>  N:           403
#> 
#> Intervention means:
#>    intervention estimate      se ci_lower ci_upper
#>          <char>    <num>   <num>    <num>    <num>
#> 1:         quit   0.7419 0.02182   0.6992   0.7847
#> 2:     continue   0.5931 0.02450   0.5450   0.6411
#> 
#> Contrasts:
#>          comparison estimate      se ci_lower ci_upper
#>              <char>    <num>   <num>    <num>    <num>
#> 1: quit vs continue   0.1489 0.03173   0.0867   0.2111
```

### Risk difference (bootstrap)

``` r
res_rd_bs <- contrast(
  fit_m_bin,
  interventions = list(quit = static(1), continue = static(0)),
  reference = "continue",
  type = "difference",
  ci_method = "bootstrap",
  n_boot = 200L
)
res_rd_bs
#> <causatr_result>
#>  Method:      matching
#>  Contrast:    difference
#>  CI method:   bootstrap
#>  N:           403
#> 
#> Intervention means:
#>    intervention estimate      se ci_lower ci_upper
#>          <char>    <num>   <num>    <num>    <num>
#> 1:         quit   0.7419 0.02156   0.6997   0.7842
#> 2:     continue   0.5931 0.02679   0.5406   0.6456
#> 
#> Contrasts:
#>          comparison estimate      se ci_lower ci_upper
#>              <char>    <num>   <num>    <num>    <num>
#> 1: quit vs continue   0.1489 0.03125  0.08764   0.2101
```

### Risk ratio

``` r
res_rr <- contrast(
  fit_m_bin,
  interventions = list(quit = static(1), continue = static(0)),
  reference = "continue",
  type = "ratio",
  ci_method = "sandwich"
)
res_rr
#> <causatr_result>
#>  Method:      matching
#>  Contrast:    ratio
#>  CI method:   sandwich
#>  N:           403
#> 
#> Intervention means:
#>    intervention estimate      se ci_lower ci_upper
#>          <char>    <num>   <num>    <num>    <num>
#> 1:         quit   0.7419 0.02182   0.6992   0.7847
#> 2:     continue   0.5931 0.02450   0.5450   0.6411
#> 
#> Contrasts:
#>          comparison estimate      se ci_lower ci_upper
#>              <char>    <num>   <num>    <num>    <num>
#> 1: quit vs continue    1.251 0.06145    1.131    1.371
```

### Odds ratio

``` r
res_or <- contrast(
  fit_m_bin,
  interventions = list(quit = static(1), continue = static(0)),
  reference = "continue",
  type = "or",
  ci_method = "sandwich"
)
res_or
#> <causatr_result>
#>  Method:      matching
#>  Contrast:    or
#>  CI method:   sandwich
#>  N:           403
#> 
#> Intervention means:
#>    intervention estimate      se ci_lower ci_upper
#>          <char>    <num>   <num>    <num>    <num>
#> 1:         quit   0.7419 0.02182   0.6992   0.7847
#> 2:     continue   0.5931 0.02450   0.5450   0.6411
#> 
#> Contrasts:
#>          comparison estimate     se ci_lower ci_upper
#>              <char>    <num>  <num>    <num>    <num>
#> 1: quit vs continue    1.973 0.2912    1.402    2.543
```

### Risk ratio (bootstrap)

``` r
res_rr_bs <- contrast(
  fit_m_bin,
  interventions = list(quit = static(1), continue = static(0)),
  reference = "continue",
  type = "ratio",
  ci_method = "bootstrap",
  n_boot = 200L
)
res_rr_bs
#> <causatr_result>
#>  Method:      matching
#>  Contrast:    ratio
#>  CI method:   bootstrap
#>  N:           403
#> 
#> Intervention means:
#>    intervention estimate      se ci_lower ci_upper
#>          <char>    <num>   <num>    <num>    <num>
#> 1:         quit   0.7419 0.02102   0.7007   0.7831
#> 2:     continue   0.5931 0.02395   0.5461   0.6400
#> 
#> Contrasts:
#>          comparison estimate      se ci_lower ci_upper
#>              <char>    <num>   <num>    <num>    <num>
#> 1: quit vs continue    1.251 0.06201     1.13    1.373
```

## Comparing estimands

``` r
results_list <- list(
  data.frame(
    estimand = "ATT",
    estimate = res_att_sw$contrasts$estimate[1],
    ci_lower = res_att_sw$contrasts$ci_lower[1],
    ci_upper = res_att_sw$contrasts$ci_upper[1]
  ),
  data.frame(
    estimand = "ATC",
    estimate = res_atc$contrasts$estimate[1],
    ci_lower = res_atc$contrasts$ci_lower[1],
    ci_upper = res_atc$contrasts$ci_upper[1]
  )
)
if (has_optmatch) {
  results_list <- c(list(data.frame(
    estimand = "ATE",
    estimate = res_ate_sw$contrasts$estimate[1],
    ci_lower = res_ate_sw$contrasts$ci_lower[1],
    ci_upper = res_ate_sw$contrasts$ci_upper[1]
  )), results_list)
}
est_df <- do.call(rbind, results_list)

tinyplot(
  estimate ~ estimand,
  data = est_df,
  type = "pointrange",
  ymin = est_df$ci_lower,
  ymax = est_df$ci_upper,
  xlab = "Estimand",
  ylab = "Effect on weight change (kg)",
  main = "Matching estimates by estimand"
)
abline(h = 0, lty = 2, col = "grey40")
```

<img
src="vignettes/matching.markdown_strict_files/figure-markdown_strict/unnamed-chunk-14-1.png"
data-fig-alt="Point estimates and confidence intervals for ATT and ATC (and ATE if optmatch is available) estimated by matching." />

## Diagnostics

After fitting a matching model, use `diagnose()` to assess covariate
balance and match quality. Balance diagnostics use
[cobalt](https://ngreifer.github.io/cobalt/) to compare SMDs before and
after matching.

``` r
diag <- diagnose(fit_m)
diag
#> <causatr_diag>
#>  Method:matching
#> 
#> Positivity (propensity score):
#>       statistic        value
#>          <char>        <num>
#>             min 2.667709e-07
#>             q25 1.739280e-01
#>          median 2.356540e-01
#>             q75 3.232473e-01
#>             max 8.723380e-01
#>   n_below_lower 6.000000e+00
#>   n_above_upper 0.000000e+00
#>    n_violations 6.000000e+00
#>  pct_violations 3.800000e-01
#> 
#> Covariate balance:
#> Balance Measures
#>                          Type Diff.Un V.Ratio.Un Diff.Adj    M.Threshold
#> distance             Distance  0.5697     1.4784   0.0468 Balanced, <0.1
#> sex                    Binary -0.1604          .   0.0299 Balanced, <0.1
#> age                   Contin.  0.2771     1.0731  -0.0116 Balanced, <0.1
#> I(age^2)              Contin.  0.2715     1.1632  -0.0123 Balanced, <0.1
#> race                   Binary -0.1993          .   0.0261 Balanced, <0.1
#> factor(education)_0    Binary  0.0394          .  -0.0250 Balanced, <0.1
#> factor(education)_1    Binary -0.0835          .   0.0000 Balanced, <0.1
#> factor(education)_2    Binary  0.0654          .   0.0501 Balanced, <0.1
#> factor(education)_3    Binary  0.0221          .   0.0000 Balanced, <0.1
#> factor(education)_4    Binary -0.0134          .  -0.0289 Balanced, <0.1
#> factor(education)_5    Binary -0.0762          .  -0.0353 Balanced, <0.1
#> factor(education)_6    Binary  0.0461          .   0.0152 Balanced, <0.1
#> factor(education)_7    Binary  0.0025          .   0.0421 Balanced, <0.1
#> factor(education)_8    Binary  0.0386          .  -0.0166 Balanced, <0.1
#> factor(education)_9    Binary  0.0170          .   0.0109 Balanced, <0.1
#> factor(education)_10   Binary -0.0879          .  -0.0184 Balanced, <0.1
#> factor(education)_11   Binary -0.1159          .  -0.0229 Balanced, <0.1
#> factor(education)_12   Binary -0.0475          .  -0.0204 Balanced, <0.1
#> factor(education)_13   Binary  0.0088          .   0.0428 Balanced, <0.1
#> factor(education)_15   Binary -0.0759          .  -0.0410 Balanced, <0.1
#> factor(education)_16   Binary  0.1221          .   0.0088 Balanced, <0.1
#> factor(education)_17   Binary  0.0823          .   0.0298 Balanced, <0.1
#> smokeintensity        Contin. -0.2087     1.1679  -0.0014 Balanced, <0.1
#> I(smokeintensity^2)   Contin. -0.1246     1.1519   0.0148 Balanced, <0.1
#> smokeyrs              Contin.  0.1526     1.1846  -0.0421 Balanced, <0.1
#> I(smokeyrs^2)         Contin.  0.1675     1.3279  -0.0342 Balanced, <0.1
#> factor(exercise)_0     Binary -0.1307          .   0.0068 Balanced, <0.1
#> factor(exercise)_1     Binary  0.0397          .  -0.0851 Balanced, <0.1
#> factor(exercise)_2     Binary  0.0565          .   0.0808 Balanced, <0.1
#> factor(active)_0       Binary -0.0721          .  -0.0251 Balanced, <0.1
#> factor(active)_1       Binary  0.0268          .  -0.0099 Balanced, <0.1
#> factor(active)_2       Binary  0.0706          .   0.0552 Balanced, <0.1
#> wt71                  Contin.  0.1313     1.0606  -0.0126 Balanced, <0.1
#> I(wt71^2)             Contin.  0.1247     1.0876  -0.0139 Balanced, <0.1
#>                      V.Ratio.Adj
#> distance                  1.2040
#> sex                            .
#> age                       0.9930
#> I(age^2)                  1.0026
#> race                           .
#> factor(education)_0            .
#> factor(education)_1            .
#> factor(education)_2            .
#> factor(education)_3            .
#> factor(education)_4            .
#> factor(education)_5            .
#> factor(education)_6            .
#> factor(education)_7            .
#> factor(education)_8            .
#> factor(education)_9            .
#> factor(education)_10           .
#> factor(education)_11           .
#> factor(education)_12           .
#> factor(education)_13           .
#> factor(education)_15           .
#> factor(education)_16           .
#> factor(education)_17           .
#> smokeintensity            1.0721
#> I(smokeintensity^2)       1.0921
#> smokeyrs                  1.0200
#> I(smokeyrs^2)             1.0645
#> factor(exercise)_0             .
#> factor(exercise)_1             .
#> factor(exercise)_2             .
#> factor(active)_0               .
#> factor(active)_1               .
#> factor(active)_2               .
#> wt71                      0.9769
#> I(wt71^2)                 0.9424
#> 
#> Balance tally for mean differences
#>                    count
#> Balanced, <0.1        34
#> Not Balanced, >0.1     0
#> 
#> Variable with the greatest mean difference
#>            Variable Diff.Adj    M.Threshold
#>  factor(exercise)_1  -0.0851 Balanced, <0.1
#> 
#> Sample sizes
#>           Control Treated
#> All          1163     403
#> Matched       403     403
#> Unmatched     760       0
#> 
#> Match quality:
#>     statistic  value
#>        <char>  <num>
#>       n_total 1566.0
#>     n_matched  806.0
#>   n_discarded  760.0
#>  pct_retained   51.5
```

### Love plot

``` r
plot(diag)
```

<img
src="vignettes/matching.markdown_strict_files/figure-markdown_strict/unnamed-chunk-16-1.png"
data-fig-alt="Love plot showing covariate balance before and after propensity score matching." />

### Match quality

``` r
diag$match_quality
#>       statistic  value
#>          <char>  <num>
#> 1:      n_total 1566.0
#> 2:    n_matched  806.0
#> 3:  n_discarded  760.0
#> 4: pct_retained   51.5
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
</tr>
</thead>
<tbody>
<tr>
<td>Binary</td>
<td>Continuous</td>
<td>Difference</td>
<td>Sandwich</td>
<td>ATT</td>
</tr>
<tr>
<td>Binary</td>
<td>Continuous</td>
<td>Difference</td>
<td>Bootstrap</td>
<td>ATT</td>
</tr>
<tr>
<td>Binary</td>
<td>Continuous</td>
<td>Difference</td>
<td>Sandwich</td>
<td>ATE</td>
</tr>
<tr>
<td>Binary</td>
<td>Continuous</td>
<td>Difference</td>
<td>Bootstrap</td>
<td>ATE</td>
</tr>
<tr>
<td>Binary</td>
<td>Continuous</td>
<td>Difference</td>
<td>Sandwich</td>
<td>ATC</td>
</tr>
<tr>
<td>Binary</td>
<td>Binary</td>
<td>Difference</td>
<td>Sandwich</td>
<td>ATT</td>
</tr>
<tr>
<td>Binary</td>
<td>Binary</td>
<td>Difference</td>
<td>Bootstrap</td>
<td>ATT</td>
</tr>
<tr>
<td>Binary</td>
<td>Binary</td>
<td>Ratio</td>
<td>Sandwich</td>
<td>ATT</td>
</tr>
<tr>
<td>Binary</td>
<td>Binary</td>
<td>OR</td>
<td>Sandwich</td>
<td>ATT</td>
</tr>
<tr>
<td>Binary</td>
<td>Binary</td>
<td>Ratio</td>
<td>Bootstrap</td>
<td>ATT</td>
</tr>
</tbody>
</table>

## References

Hernán MA, Robins JM (2025). *Causal Inference: What If*. Chapman &
Hall/CRC. Chapter 15: Outcome regression and propensity scores.

Ho DE, Imai K, King G, Stuart EA (2011). MatchIt: Nonparametric
Preprocessing for Parametric Causal Inference. *Journal of Statistical
Software* 42(8):1-28.
