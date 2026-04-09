

# Clamp treatment within bounds

[**Source code**](https://github.com/etverse/causatr/tree/main/R/interventions.R#L115)

## Description

Creates a modified treatment policy that clamps each individual’s
observed treatment to lie within
<code style="white-space: pre;">\[lower, upper\]</code>. Values below
<code>lower</code> are set to <code>lower</code>; values above
<code>upper</code> are set to <code>upper</code>.

## Usage

<pre><code class='language-R'>threshold(lower = -Inf, upper = Inf)
</code></pre>

## Arguments

<table role="presentation">
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="lower">lower</code>
</td>
<td>
Numeric. Lower bound (use <code>-Inf</code> for no lower bound).
</td>
</tr>
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="upper">upper</code>
</td>
<td>
Numeric. Upper bound (use <code>Inf</code> for no upper bound).
</td>
</tr>
</table>

## Value

A <code>causatr_intervention</code> object.

## See Also

<code>static()</code>, <code>shift()</code>, <code>scale()</code>,
<code>dynamic()</code>, <code>ipsi()</code>

## Examples

``` r
library("causatr")

fit <- causat(nhefs, outcome = "wt82_71", treatment = "smokeintensity",
              confounders = ~ sex + age + wt71)
contrast(fit, interventions = list(
  capped20 = threshold(0, 20),
  observed  = shift(0)
))
```
