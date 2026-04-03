

# Set treatment to a fixed value

[**Source code**](https://github.com/etverse/causatr/tree/main/R/interventions.R#L27)

## Description

Creates a static intervention that sets the treatment to a fixed value
for all individuals at all time points. The most common intervention
type, corresponding to "always treat" (<code>static(1)</code>) or "never
treat" (<code>static(0)</code>).

## Usage

<pre><code class='language-R'>static(value)
</code></pre>

## Arguments

<table role="presentation">
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="value">value</code>
</td>
<td>
The fixed treatment value.
</td>
</tr>
</table>

## Value

A <code>causatr_intervention</code> object.

## References

Hernán MA, Robins JM (2025). <em>Causal Inference: What If</em>. Chapman
& Hall/CRC. Chapter 1 (static interventions).

## See Also

<code>data.table::shift()</code>, <code>dynamic()</code>,
<code>scale()</code>, <code>threshold()</code>, <code>ipsi()</code>,
<code>contrast()</code>

## Examples

``` r
library("causatr")

data("nhefs", package = "causatr")
fit <- causat(nhefs, outcome = "wt82_71", treatment = "qsmk",
              confounders = ~ sex + age + wt71)
contrast(fit, interventions = list(quit = static(1), continue = static(0)))
```
