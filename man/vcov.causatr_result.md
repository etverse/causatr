

# Variance-covariance matrix for a causatr result

[**Source code**](https://github.com/etverse/causatr/tree/main/R/coef.R#L43)

## Description

Returns the variance-covariance matrix of the intervention-specific
marginal means (E\[Y^a\]) from a <code>causatr_result</code> object.

## Usage

<pre><code class='language-R'>## S3 method for class 'causatr_result'
vcov(object, ...)
</code></pre>

## Arguments

<table role="presentation">
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="object">object</code>
</td>
<td>
A <code>causatr_result</code> object.
</td>
</tr>
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="...">…</code>
</td>
<td>
Currently unused.
</td>
</tr>
</table>

## Value

A named k x k matrix (k = number of interventions).

## See Also

<code>coef.causatr_result()</code>,
<code>confint.causatr_result()</code>, <code>contrast()</code>

## Examples

``` r
library("causatr")

result <- contrast(fit, interventions = list(a1 = static(1), a0 = static(0)))
vcov(result)
```
