

# Multiply treatment by a fixed factor

[**Source code**](https://github.com/etverse/causatr/tree/main/R/interventions.R#L87)

## Description

Creates a modified treatment policy that multiplies each individual’s
observed treatment by <code>factor</code>. Useful for proportional
reductions or increases in continuous treatments.

## Usage

<pre><code class='language-R'>scale(factor)
</code></pre>

## Arguments

<table role="presentation">
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="factor">factor</code>
</td>
<td>
Numeric. The multiplicative factor.
</td>
</tr>
</table>

## Value

A <code>causatr_intervention</code> object.

## See Also

<code>static()</code>, <code>shift()</code>, <code>threshold()</code>,
<code>dynamic()</code>, <code>ipsi()</code>

## Examples

``` r
library("causatr")

fit <- causat(nhefs, outcome = "wt82_71", treatment = "smokeintensity",
              confounders = ~ sex + age + wt71)
contrast(fit, interventions = list(
  halved = scale(0.5),
  observed = scale(1)
))
```
