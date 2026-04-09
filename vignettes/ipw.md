# Inverse probability weighting with causatr


Inverse probability weighting (IPW) estimates causal effects by
reweighting the observed data so that the distribution of confounders is
balanced across treatment groups. causatr delegates weight estimation to
[WeightIt](https://ngreifer.github.io/WeightIt/) and fits a weighted
marginal structural model (MSM) via `WeightIt::glm_weightit()`, which
provides M-estimation sandwich SEs that account for weight estimation
uncertainty.

This vignette demonstrates IPW in causatr using the NHEFS dataset from
Hernán & Robins (2025), covering every supported combination of
treatment type, outcome type, contrast scale, and inference method for
time-fixed treatments.

**Note:** IPW currently supports only `static()` interventions. For
modified treatment policies (shift, scale, threshold, dynamic), use
g-computation (see `vignette("gcomp")`).

## Setup

``` r
library(causatr)
library(tinyplot)

data("nhefs")

nhefs_complete <- nhefs[!is.na(nhefs$wt82_71), ]

nhefs_complete$gained_weight <- as.integer(nhefs_complete$wt82_71 > 0)

nhefs_complete$sex <- factor(nhefs_complete$sex, levels = 0:1, labels = c("Male", "Female"))
```

## Binary treatment, continuous outcome

This replicates the IPW analysis from Chapter 12 of Hernán & Robins: the
average causal effect of quitting smoking (`qsmk`) on weight change
(`wt82_71`).

### ATE with sandwich SE

``` r
fit_ipw <- causat(
  nhefs_complete,
  outcome = "wt82_71",
  treatment = "qsmk",
  confounders = ~ sex + age + I(age^2) + race + factor(education) +
    smokeintensity + I(smokeintensity^2) + smokeyrs + I(smokeyrs^2) +
    factor(exercise) + factor(active) + wt71 + I(wt71^2),
  method = "ipw",
  estimand = "ATE"
)
#> Warning: Missing values are present in the covariates. See `?method_glm`
#> (`?WeightIt::method_glm()`) for information on how these are handled.

res_ate_sw <- contrast(
  fit_ipw,
  interventions = list(quit = static(1), continue = static(0)),
  reference = "continue",
  type = "difference",
  ci_method = "sandwich"
)
res_ate_sw
#> <causatr_result>
#>  Method:    IPW
#>  Estimand:  ATE
#>  Contrast:  Difference
#>  CI method: sandwich
#>  N:         1679
#> 
#> Intervention means:
#>    intervention estimate     se ci_lower ci_upper
#>          <char>    <num>  <num>    <num>    <num>
#> 1:         quit    5.305 0.4230    4.476    6.134
#> 2:     continue    1.942 0.2067    1.537    2.347
#> 
#> Contrasts:
#>          comparison estimate     se ci_lower ci_upper
#>              <char>    <num>  <num>    <num>    <num>
#> 1: quit vs continue    3.363 0.4637    2.454    4.272
```

The book reports an IPW estimate of ATE ≈ 3.4 kg (95% CI: 2.4, 4.5).

### ATE with bootstrap SE

``` r
res_ate_bs <- contrast(
  fit_ipw,
  interventions = list(quit = static(1), continue = static(0)),
  reference = "continue",
  type = "difference",
  ci_method = "bootstrap",
  n_boot = 200L
)
res_ate_bs
#> <causatr_result>
#>  Method:    IPW
#>  Estimand:  ATE
#>  Contrast:  Difference
#>  CI method: bootstrap
#>  N:         1679
#> 
#> Intervention means:
#>    intervention estimate     se ci_lower ci_upper
#>          <char>    <num>  <num>    <num>    <num>
#> 1:         quit    5.305 0.4513    4.421    6.190
#> 2:     continue    1.942 0.2280    1.495    2.389
#> 
#> Contrasts:
#>          comparison estimate    se ci_lower ci_upper
#>              <char>    <num> <num>    <num>    <num>
#> 1: quit vs continue    3.363 0.529    2.326      4.4
```

### ATT estimand

For IPW, the estimand is fixed at fitting time because it determines the
weights. To estimate the ATT, refit with `estimand = "ATT"`.

``` r
fit_ipw_att <- causat(
  nhefs_complete,
  outcome = "wt82_71",
  treatment = "qsmk",
  confounders = ~ sex + age + I(age^2) + race + factor(education) +
    smokeintensity + I(smokeintensity^2) + smokeyrs + I(smokeyrs^2) +
    factor(exercise) + factor(active) + wt71 + I(wt71^2),
  method = "ipw",
  estimand = "ATT"
)
#> Warning: Missing values are present in the covariates. See `?method_glm`
#> (`?WeightIt::method_glm()`) for information on how these are handled.

res_att <- contrast(
  fit_ipw_att,
  interventions = list(quit = static(1), continue = static(0)),
  reference = "continue",
  type = "difference",
  ci_method = "sandwich"
)
res_att
#> <causatr_result>
#>  Method:    IPW
#>  Estimand:  ATT
#>  Contrast:  Difference
#>  CI method: sandwich
#>  N:         437
#> 
#> Intervention means:
#>    intervention estimate     se ci_lower ci_upper
#>          <char>    <num>  <num>    <num>    <num>
#> 1:         quit    4.497 0.4162   3.6812    5.313
#> 2:     continue    1.400 0.2779   0.8555    1.945
#> 
#> Contrasts:
#>          comparison estimate     se ci_lower ci_upper
#>              <char>    <num>  <num>    <num>    <num>
#> 1: quit vs continue    3.097 0.4687    2.178    4.015
```

### ATC estimand

``` r
fit_ipw_atc <- causat(
  nhefs_complete,
  outcome = "wt82_71",
  treatment = "qsmk",
  confounders = ~ sex + age + I(age^2) + race + factor(education) +
    smokeintensity + I(smokeintensity^2) + smokeyrs + I(smokeyrs^2) +
    factor(exercise) + factor(active) + wt71 + I(wt71^2),
  method = "ipw",
  estimand = "ATC"
)
#> Warning: Missing values are present in the covariates. See `?method_glm`
#> (`?WeightIt::method_glm()`) for information on how these are handled.

res_atc <- contrast(
  fit_ipw_atc,
  interventions = list(quit = static(1), continue = static(0)),
  reference = "continue",
  type = "difference",
  ci_method = "sandwich"
)
res_atc
#> <causatr_result>
#>  Method:    IPW
#>  Estimand:  ATC
#>  Contrast:  Difference
#>  CI method: sandwich
#>  N:         1242
#> 
#> Intervention means:
#>    intervention estimate     se ci_lower ci_upper
#>          <char>    <num>  <num>    <num>    <num>
#> 1:         quit    5.592 0.4594    4.691    6.492
#> 2:     continue    2.133 0.2083    1.724    2.541
#> 
#> Contrasts:
#>          comparison estimate     se ci_lower ci_upper
#>              <char>    <num>  <num>    <num>    <num>
#> 1: quit vs continue    3.459 0.4964    2.486    4.432
```

### ATT with bootstrap SE

``` r
res_att_bs <- contrast(
  fit_ipw_att,
  interventions = list(quit = static(1), continue = static(0)),
  reference = "continue",
  type = "difference",
  ci_method = "bootstrap",
  n_boot = 200L
)
res_att_bs
#> <causatr_result>
#>  Method:    IPW
#>  Estimand:  ATT
#>  Contrast:  Difference
#>  CI method: bootstrap
#>  N:         437
#> 
#> Intervention means:
#>    intervention estimate     se ci_lower ci_upper
#>          <char>    <num>  <num>    <num>    <num>
#> 1:         quit    4.497 0.4292   3.6557    5.338
#> 2:     continue    1.400 0.3018   0.8088    1.992
#> 
#> Contrasts:
#>          comparison estimate     se ci_lower ci_upper
#>              <char>    <num>  <num>    <num>    <num>
#> 1: quit vs continue    3.097 0.5032     2.11    4.083
```

### Extracting results programmatically

``` r
coef(res_ate_sw)
#>     quit continue 
#> 5.305091 1.941992
confint(res_ate_sw)
#>             lower    upper
#> quit     4.476028 6.134154
#> continue 1.536817 2.347167
vcov(res_ate_sw)
#>                 quit    continue
#> quit     0.178928248 0.003327337
#> continue 0.003327337 0.042735593
tidy(res_ate_sw)
#>               term estimate std.error     type conf.low conf.high
#> 1 quit vs continue 3.363099 0.4636908 contrast 2.454282  4.271917
glance(res_ate_sw)
#>   method estimand contrast_type ci_method    n n_interventions
#> 1    ipw      ATE    difference  sandwich 1679               2
```

## Binary treatment, binary outcome

Using the `gained_weight` indicator as a binary outcome. With IPW, the
MSM is `Y ~ A` (a linear probability model), so `contrast()` computes
marginal risks and derives contrasts via the delta method.

### Risk difference (sandwich)

``` r
fit_ipw_bin <- causat(
  nhefs_complete,
  outcome = "gained_weight",
  treatment = "qsmk",
  confounders = ~ sex + age + I(age^2) + race + factor(education) +
    smokeintensity + I(smokeintensity^2) + smokeyrs + I(smokeyrs^2) +
    factor(exercise) + factor(active) + wt71 + I(wt71^2),
  method = "ipw"
)
#> Warning: Missing values are present in the covariates. See `?method_glm`
#> (`?WeightIt::method_glm()`) for information on how these are handled.

res_rd <- contrast(
  fit_ipw_bin,
  interventions = list(quit = static(1), continue = static(0)),
  reference = "continue",
  type = "difference",
  ci_method = "sandwich"
)
res_rd
#> <causatr_result>
#>  Method:    IPW
#>  Estimand:  ATE
#>  Contrast:  Difference
#>  CI method: sandwich
#>  N:         1679
#> 
#> Intervention means:
#>    intervention estimate      se ci_lower ci_upper
#>          <char>    <num>   <num>    <num>    <num>
#> 1:         quit   0.7815 0.01991   0.7425   0.8205
#> 2:     continue   0.6500 0.01349   0.6235   0.6764
#> 
#> Contrasts:
#>          comparison estimate     se ci_lower ci_upper
#>              <char>    <num>  <num>    <num>    <num>
#> 1: quit vs continue   0.1315 0.0236  0.08528   0.1778
```

### Risk difference (bootstrap)

``` r
res_rd_bs <- contrast(
  fit_ipw_bin,
  interventions = list(quit = static(1), continue = static(0)),
  reference = "continue",
  type = "difference",
  ci_method = "bootstrap",
  n_boot = 200L
)
res_rd_bs
#> <causatr_result>
#>  Method:    IPW
#>  Estimand:  ATE
#>  Contrast:  Difference
#>  CI method: bootstrap
#>  N:         1679
#> 
#> Intervention means:
#>    intervention estimate      se ci_lower ci_upper
#>          <char>    <num>   <num>    <num>    <num>
#> 1:         quit   0.7815 0.02088   0.7406   0.8224
#> 2:     continue   0.6500 0.01450   0.6216   0.6784
#> 
#> Contrasts:
#>          comparison estimate     se ci_lower ci_upper
#>              <char>    <num>  <num>    <num>    <num>
#> 1: quit vs continue   0.1315 0.0256  0.08136   0.1817
```

### Risk ratio

``` r
res_rr <- contrast(
  fit_ipw_bin,
  interventions = list(quit = static(1), continue = static(0)),
  reference = "continue",
  type = "ratio",
  ci_method = "sandwich"
)
res_rr
#> <causatr_result>
#>  Method:    IPW
#>  Estimand:  ATE
#>  Contrast:  Ratio
#>  CI method: sandwich
#>  N:         1679
#> 
#> Intervention means:
#>    intervention estimate      se ci_lower ci_upper
#>          <char>    <num>   <num>    <num>    <num>
#> 1:         quit   0.7815 0.01991   0.7425   0.8205
#> 2:     continue   0.6500 0.01349   0.6235   0.6764
#> 
#> Contrasts:
#>          comparison estimate      se ci_lower ci_upper
#>              <char>    <num>   <num>    <num>    <num>
#> 1: quit vs continue    1.202 0.03872    1.126    1.278
```

### Odds ratio

``` r
res_or <- contrast(
  fit_ipw_bin,
  interventions = list(quit = static(1), continue = static(0)),
  reference = "continue",
  type = "or",
  ci_method = "sandwich"
)
res_or
#> <causatr_result>
#>  Method:    IPW
#>  Estimand:  ATE
#>  Contrast:  Odds ratio
#>  CI method: sandwich
#>  N:         1679
#> 
#> Intervention means:
#>    intervention estimate      se ci_lower ci_upper
#>          <char>    <num>   <num>    <num>    <num>
#> 1:         quit   0.7815 0.01991   0.7425   0.8205
#> 2:     continue   0.6500 0.01349   0.6235   0.6764
#> 
#> Contrasts:
#>          comparison estimate     se ci_lower ci_upper
#>              <char>    <num>  <num>    <num>    <num>
#> 1: quit vs continue    1.926 0.2478     1.44    2.412
```

### Risk ratio (bootstrap)

``` r
res_rr_bs <- contrast(
  fit_ipw_bin,
  interventions = list(quit = static(1), continue = static(0)),
  reference = "continue",
  type = "ratio",
  ci_method = "bootstrap",
  n_boot = 200L
)
res_rr_bs
#> <causatr_result>
#>  Method:    IPW
#>  Estimand:  ATE
#>  Contrast:  Ratio
#>  CI method: bootstrap
#>  N:         1679
#> 
#> Intervention means:
#>    intervention estimate      se ci_lower ci_upper
#>          <char>    <num>   <num>    <num>    <num>
#> 1:         quit   0.7815 0.02295   0.7365   0.8265
#> 2:     continue   0.6500 0.01335   0.6238   0.6761
#> 
#> Contrasts:
#>          comparison estimate      se ci_lower ci_upper
#>              <char>    <num>   <num>    <num>    <num>
#> 1: quit vs continue    1.202 0.03853    1.127    1.278
```

## Comparing estimands

``` r
est_df <- data.frame(
  estimand = c("ATE", "ATT", "ATC"),
  estimate = c(
    res_ate_sw$contrasts$estimate[1],
    res_att$contrasts$estimate[1],
    res_atc$contrasts$estimate[1]
  ),
  ci_lower = c(
    res_ate_sw$contrasts$ci_lower[1],
    res_att$contrasts$ci_lower[1],
    res_atc$contrasts$ci_lower[1]
  ),
  ci_upper = c(
    res_ate_sw$contrasts$ci_upper[1],
    res_att$contrasts$ci_upper[1],
    res_atc$contrasts$ci_upper[1]
  )
)

tinyplot(
  estimate ~ estimand,
  data = est_df,
  type = "pointrange",
  ymin = est_df$ci_lower,
  ymax = est_df$ci_upper,
  xlab = "Estimand",
  ylab = "Effect on weight change (kg)",
  main = "IPW estimates by estimand"
)
abline(h = 0, lty = 2, col = "grey40")
```

<img
src="vignettes/ipw.markdown_strict_files/figure-markdown_strict/unnamed-chunk-14-1.png"
data-fig-alt="Point estimates and confidence intervals for ATE, ATT, and ATC estimated by IPW." />

## Effect modification with `by`

The `by` argument stratifies causal effect estimates by levels of a
variable. Here we examine whether the IPW-estimated effect of quitting
smoking differs by sex.

``` r
res_by_sex <- contrast(
  fit_ipw,
  interventions = list(quit = static(1), continue = static(0)),
  reference = "continue",
  type = "difference",
  ci_method = "sandwich",
  by = "sex"
)
res_by_sex
#> <causatr_result>
#>  Method:    IPW
#>  Estimand:  ATE
#>  Contrast:  Difference
#>  CI method: sandwich
#>  N:         1679
#> 
#> Intervention means (by subgroup):
#>    intervention estimate     se ci_lower ci_upper     by
#>          <char>    <num>  <num>    <num>    <num> <char>
#> 1:         quit    5.305 0.4230    4.476    6.134   Male
#> 2:     continue    1.942 0.2067    1.537    2.347   Male
#> 3:         quit    5.305 0.4230    4.476    6.134 Female
#> 4:     continue    1.942 0.2067    1.537    2.347 Female
#> 
#> Contrasts (by subgroup):
#>          comparison estimate     se ci_lower ci_upper     by
#>              <char>    <num>  <num>    <num>    <num> <char>
#> 1: quit vs continue    3.363 0.4637    2.454    4.272   Male
#> 2: quit vs continue    3.363 0.4637    2.454    4.272 Female
```

## Tidy and glance

causatr results work with the broom ecosystem via `tidy()` and
`glance()`:

``` r
tidy(res_ate_sw)
#>               term estimate std.error     type conf.low conf.high
#> 1 quit vs continue 3.363099 0.4636908 contrast 2.454282  4.271917
glance(res_ate_sw)
#>   method estimand contrast_type ci_method    n n_interventions
#> 1    ipw      ATE    difference  sandwich 1679               2
```

## Forest plot

The `plot()` method produces a forest plot using the `forrest` package.
Forest plots are most useful when displaying multiple estimates, such as
effect modification results:

``` r
plot(res_by_sex)
```

<img
src="vignettes/ipw.markdown_strict_files/figure-markdown_strict/unnamed-chunk-17-1.png"
data-fig-alt="Forest plot of the IPW-estimated effect of quitting smoking on weight change, stratified by sex." />

## Diagnostics

After fitting an IPW model, use `diagnose()` to assess positivity and
covariate balance. The balance diagnostics use
[cobalt](https://ngreifer.github.io/cobalt/) to compute standardised
mean differences (SMD) before and after weighting.

``` r
diag <- diagnose(fit_ipw)
#> Warning: Missing values exist in the covariates. Displayed values omit these
#> observations.
diag
#> <causatr_diag>
#>  Method:ipw
#> 
#> Positivity (propensity score):
#>       statistic        value
#>          <char>        <num>
#>             min 2.849842e-07
#>             q25 1.772053e-01
#>          median 2.411414e-01
#>             q75 3.241364e-01
#>             max 8.627023e-01
#>   n_below_lower 6.000000e+00
#>   n_above_upper 0.000000e+00
#>    n_violations 6.000000e+00
#>  pct_violations 3.600000e-01
#> 
#> Covariate balance:
#> Balance Measures
#>                    Type Diff.Un V.Ratio.Un Diff.Adj    M.Threshold V.Ratio.Adj
#> prop.score     Distance  0.6076     1.4587   0.0139 Balanced, <0.1      0.9828
#> sex_Female       Binary -0.1506          .  -0.0262 Balanced, <0.1           .
#> age             Contin.  0.2744     1.0892  -0.0015 Balanced, <0.1      1.0056
#> race             Binary -0.1549          .   0.0003 Balanced, <0.1           .
#> education       Contin.  0.0491     1.2310   0.0228 Balanced, <0.1      0.9644
#> education:<NA>   Binary  0.0554          .   0.0134 Balanced, <0.1           .
#> smokeintensity  Contin. -0.2287     1.1550  -0.0148 Balanced, <0.1      0.9666
#> smokeyrs        Contin.  0.1596     1.1650  -0.0063 Balanced, <0.1      1.0044
#> exercise        Contin.  0.1054     0.9054  -0.0207 Balanced, <0.1      0.9831
#> active          Contin.  0.0987     1.0563  -0.0046 Balanced, <0.1      0.9756
#> wt71            Contin.  0.1280     1.0877  -0.0078 Balanced, <0.1      0.9935
#> 
#> Balance tally for mean differences
#>                    count
#> Balanced, <0.1        11
#> Not Balanced, >0.1     0
#> 
#> Variable with the greatest mean difference
#>    Variable Diff.Adj    M.Threshold
#>  sex_Female  -0.0262 Balanced, <0.1
#> 
#> Effective sample sizes
#>            Control Treated
#> Unadjusted 1242.    437.  
#> Adjusted   1201.83  351.53
#> 
#> Weight distribution:
#>    group     n     mean        sd      min       max      ess
#>   <char> <int>    <num>     <num>    <num>     <num>    <num>
#>  treated   437 3.818806 1.8851704 1.159148 16.658502  351.530
#>  control  1242 1.351697 0.2472201 1.000000  3.567222 1201.830
#>  overall  1679 1.993821 1.4632768 1.000000 16.658502 1091.467
```

### Love plot

A Love plot visualises covariate balance. Covariates with absolute SMD
below the threshold (dashed line at 0.1) are considered well-balanced.

``` r
plot(diag)
#> Warning: Missing values exist in the covariates. Displayed values omit these
#> observations.
```

<img
src="vignettes/ipw.markdown_strict_files/figure-markdown_strict/unnamed-chunk-19-1.png"
data-fig-alt="Love plot showing covariate balance before and after IPW weighting." />

### Weight distribution

Extreme weights can inflate variance. The weight summary shows the
effective sample size (ESS) — a large drop from the nominal sample size
suggests influential weights.

``` r
diag$weights
#>      group     n     mean        sd      min       max      ess
#>     <char> <int>    <num>     <num>    <num>     <num>    <num>
#> 1: treated   437 3.818806 1.8851704 1.159148 16.658502  351.530
#> 2: control  1242 1.351697 0.2472201 1.000000  3.567222 1201.830
#> 3: overall  1679 1.993821 1.4632768 1.000000 16.658502 1091.467
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
<td>ATT</td>
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
<td>Continuous</td>
<td>Difference</td>
<td>Bootstrap</td>
<td>ATT</td>
</tr>
<tr>
<td>Binary</td>
<td>Binary</td>
<td>Difference</td>
<td>Sandwich</td>
<td>ATE</td>
</tr>
<tr>
<td>Binary</td>
<td>Binary</td>
<td>Difference</td>
<td>Bootstrap</td>
<td>ATE</td>
</tr>
<tr>
<td>Binary</td>
<td>Binary</td>
<td>Ratio</td>
<td>Sandwich</td>
<td>ATE</td>
</tr>
<tr>
<td>Binary</td>
<td>Binary</td>
<td>OR</td>
<td>Sandwich</td>
<td>ATE</td>
</tr>
<tr>
<td>Binary</td>
<td>Binary</td>
<td>Ratio</td>
<td>Bootstrap</td>
<td>ATE</td>
</tr>
</tbody>
</table>

## References

Hernán MA, Robins JM (2025). *Causal Inference: What If*. Chapman &
Hall/CRC. Chapter 12: IP weighting and marginal structural models.

Greifer N (2024). WeightIt: Weighting Methods for Covariate Balancing.
<https://ngreifer.github.io/WeightIt/>
