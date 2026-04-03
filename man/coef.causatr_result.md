

# Extract point estimates from a causatr result

[**Source code**](https://github.com/etverse/causatr/tree/main/R/coef.R#L20)

## Description

Returns the point estimates for each intervention mean (E\[Y^a\]) as a
named numeric vector.

## Usage

<pre><code class='language-R'>## S3 method for class 'causatr_result'
coef(object, ...)
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

A named numeric vector with one element per intervention.

## See Also

<code>confint.causatr_result()</code>, <code>contrast()</code>

## Examples

``` r
library("causatr")

result <- contrast(fit, interventions = list(a1 = static(1), a0 = static(0)))
coef(result)
```
