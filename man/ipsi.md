

# Incremental propensity score intervention

[**Source code**](https://github.com/etverse/causatr/tree/main/R/interventions.R#L181)

## Description

Creates an incremental propensity score intervention (IPSI) that
multiplies each individual’s odds of treatment by <code>delta</code>.
Values of <code>delta \> 1</code> increase the probability of treatment;
<code>delta \< 1</code> decrease it. This is a stochastic modified
treatment policy indexed by a single scalar.

## Usage

<pre><code class='language-R'>ipsi(delta)
</code></pre>

## Arguments

<table role="presentation">
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="delta">delta</code>
</td>
<td>
Positive numeric. The odds multiplier.
</td>
</tr>
</table>

## Value

A <code>causatr_intervention</code> object.

## References

Kennedy EH (2019). Nonparametric causal effects based on incremental
propensity score interventions. <em>Journal of the American Statistical
Association</em> 114:645–656.

## See Also

<code>static()</code>, <code>data.table::shift()</code>,
<code>dynamic()</code>, <code>scale()</code>, <code>threshold()</code>

## Examples

``` r
library("causatr")

fit <- causat(nhefs, outcome = "wt82_71", treatment = "qsmk",
              confounders = ~ sex + age + wt71)
contrast(fit, interventions = list(
  double_odds = ipsi(2),
  half_odds   = ipsi(0.5)
))
```
