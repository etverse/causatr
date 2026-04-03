

# Shift treatment by a fixed amount

[**Source code**](https://github.com/etverse/causatr/tree/main/R/interventions.R#L60)

## Description

Creates a modified treatment policy (MTP) that adds a fixed
<code>delta</code> to each individual’s observed treatment value. Useful
for continuous treatments where a population-level shift is the relevant
intervention (e.g., "reduce exposure by 10 units").

## Usage

<pre><code class='language-R'>shift(delta)
</code></pre>

## Arguments

<table role="presentation">
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="delta">delta</code>
</td>
<td>
Numeric. The amount to add to the observed treatment.
</td>
</tr>
</table>

## Value

A <code>causatr_intervention</code> object.

## References

Díaz I, Williams N, Hoffman KL, Schenck EJ (2023). Non-parametric causal
effects based on longitudinal modified treatment policies. <em>Journal
of the American Statistical Association</em> 118:846–857.

## See Also

<code>static()</code>, <code>scale()</code>, <code>threshold()</code>,
<code>dynamic()</code>, <code>ipsi()</code>

## Examples

``` r
library("causatr")

fit <- causat(nhefs, outcome = "wt82_71", treatment = "smokeintensity",
              confounders = ~ sex + age + wt71)
contrast(fit, interventions = list(
  reduce10 = shift(-10),
  observed = shift(0)
))
```
