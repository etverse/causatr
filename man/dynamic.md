

# Dynamic treatment rule

[**Source code**](https://github.com/etverse/causatr/tree/main/R/interventions.R#L143)

## Description

Creates a dynamic intervention where the treatment at each time point is
determined by a user-supplied function of the covariate history. The
function <code>rule</code> receives the current data (subset to the
current time point) and the observed treatment vector, and returns the
intervened treatment vector.

## Usage

<pre><code class='language-R'>dynamic(rule)
</code></pre>

## Arguments

<table role="presentation">
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="rule">rule</code>
</td>
<td>
A function with signature <code style="white-space: pre;">function(data,
treatment)</code> that returns a vector of treatment values of the same
length as <code>nrow(data)</code>.
</td>
</tr>
</table>

## Value

A <code>causatr_intervention</code> object.

## References

Hernán MA, Robins JM (2025). <em>Causal Inference: What If</em>. Chapman
& Hall/CRC. Chapter 19 (dynamic treatment strategies).

## See Also

<code>static()</code>, <code>data.table::shift()</code>,
<code>scale()</code>, <code>threshold()</code>, <code>ipsi()</code>

## Examples

``` r
library("causatr")

# Treat if CD4 count is below 200
cd4_rule <- dynamic(\(data, trt) ifelse(data$cd4 < 200, 1, 0))
```
