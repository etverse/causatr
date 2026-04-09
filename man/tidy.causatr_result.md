

# Tidy a causatr result

[**Source code**](https://github.com/etverse/causatr/tree/main/R/tidy.R#L29)

## Description

Returns a tidy data frame of intervention means and/or pairwise
contrasts from a <code>causatr_result</code> object, compatible with the
<a href="https://broom.tidymodels.org/">broom</a> ecosystem.

## Usage

<pre><code class='language-R'>## S3 method for class 'causatr_result'
tidy(
  x,
  which = c("contrasts", "means", "all"),
  conf.int = TRUE,
  conf.level = 0.95,
  ...
)
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
<code id="which">which</code>
</td>
<td>
Character. What to tidy: <code>“contrasts”</code> (default),
<code>“means”</code>, or <code>“all”</code> (both).
</td>
</tr>
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="conf.int">conf.int</code>
</td>
<td>
Logical. Include confidence interval columns? Default <code>TRUE</code>.
</td>
</tr>
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="conf.level">conf.level</code>
</td>
<td>
Numeric. Confidence level for intervals. Default <code>0.95</code>.
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

A data.frame with columns <code>term</code>, <code>estimate</code>,
<code>std.error</code>, <code>conf.low</code>, <code>conf.high</code>,
and <code>type</code> (either <code>“mean”</code> or
<code>“contrast”</code>).

## See Also

<code>coef.causatr_result()</code>,
<code>confint.causatr_result()</code>

## Examples

``` r
library("causatr")

result <- contrast(fit, interventions = list(a1 = static(1), a0 = static(0)))
tidy(result)
tidy(result, which = "means")
tidy(result, conf.level = 0.99)
```
