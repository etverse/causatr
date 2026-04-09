

# Glance at a causatr result

[**Source code**](https://github.com/etverse/causatr/tree/main/R/tidy.R#L101)

## Description

Returns a one-row data frame of model-level summaries from a
<code>causatr_result</code> object, compatible with the
<a href="https://broom.tidymodels.org/">broom</a> ecosystem.

## Usage

<pre><code class='language-R'>## S3 method for class 'causatr_result'
glance(x, ...)
</code></pre>

## Arguments

<table role="presentation">
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="x">x</code>
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

A one-row data.frame with columns <code>method</code>,
<code>estimand</code>, <code>contrast_type</code>,
<code>ci_method</code>, <code>n</code>, and
<code>n_interventions</code>.

## See Also

<code>coef.causatr_result()</code>

## Examples

``` r
library("causatr")

result <- contrast(fit, interventions = list(a1 = static(1), a0 = static(0)))
glance(result)
```
