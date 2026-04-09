

# Plot diagnostics for a causatr fit

[**Source code**](https://github.com/etverse/causatr/tree/main/R/plot.R#L169)

## Description

Produces a Love plot showing covariate balance before and after
adjustment. Uses <code>cobalt::love.plot()</code> if the
<code>cobalt</code> package is available.

For IPW fits, the plot shows balance before and after weighting. For
matching fits, it shows balance before and after matching. For
g-computation fits, it shows unadjusted balance only.

## Usage

<pre><code class='language-R'>## S3 method for class 'causatr_diag'
plot(x, stats = "m", abs = TRUE, thresholds = c(m = 0.1), ...)
</code></pre>

## Arguments

<table role="presentation">
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="x">x</code>
</td>
<td>
A <code>causatr_diag</code> object returned by <code>diagnose()</code>.
</td>
</tr>
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="stats">stats</code>
</td>
<td>
Character. Which balance statistic(s) to plot. Default <code>“m”</code>
(standardised mean differences). See <code>cobalt::love.plot()</code>
for options.
</td>
</tr>
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="abs">abs</code>
</td>
<td>
Logical. Whether to plot absolute values. Default <code>TRUE</code>.
</td>
</tr>
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="thresholds">thresholds</code>
</td>
<td>
Named numeric vector. Threshold lines to draw on the plot. Default
<code>c(m = 0.1)</code>.
</td>
</tr>
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="...">…</code>
</td>
<td>
Additional arguments passed to <code>cobalt::love.plot()</code>.
</td>
</tr>
</table>

## Value

A <code>ggplot</code> object (invisibly).

## See Also

<code>diagnose()</code>, <code>print.causatr_diag()</code>
