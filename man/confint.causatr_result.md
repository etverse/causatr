

# Confidence intervals for a causatr result

[**Source code**](https://github.com/etverse/causatr/tree/main/R/confint.R#L23)

## Description

Returns confidence intervals for each intervention mean (E\[Y^a\]) from
a <code>causatr_result</code> object.

## Usage

<pre><code class='language-R'>## S3 method for class 'causatr_result'
confint(object, parm, level = 0.95, ...)
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
<code id="parm">parm</code>
</td>
<td>
Ignored. Intervals are returned for all interventions.
</td>
</tr>
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="level">level</code>
</td>
<td>
Numeric. Confidence level. Default <code>0.95</code>.
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

A matrix with columns <code>“lower”</code> and <code>“upper”</code> and
one row per intervention.

## See Also

<code>coef.causatr_result()</code>, <code>contrast()</code>

## Examples

``` r
library("causatr")

result <- contrast(fit, interventions = list(a1 = static(1), a0 = static(0)))
confint(result)
```
