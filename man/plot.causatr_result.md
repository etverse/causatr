

# Plot a causatr result

[**Source code**](https://github.com/etverse/causatr/tree/main/R/plot.R#L24)

## Description

Produces a forest plot of intervention means or pairwise contrasts from
a <code>causatr_result</code> object using the <code>forrest</code>
package.

The plot adapts to the context:

<ul>
<li>

<strong>Contrast type</strong>: reference line at 0 (difference) or 1
(ratio/OR); log scale for ratio and OR.

</li>
<li>

<strong>Outcome family</strong>: axis label uses "risk" (binomial) or
"mean" (other).

</li>
<li>

<strong>Effect modification (<code>by</code>)</strong>: sections group
subgroup-specific estimates.

</li>
<li>

<strong>Fit type</strong>: title notes point vs longitudinal analysis.

</li>
</ul>

## Usage

<pre><code class='language-R'>## S3 method for class 'causatr_result'
plot(x, which = c("contrasts", "means"), ...)
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
Character. What to plot: <code>“contrasts”</code> (default) for pairwise
comparisons, or <code>“means”</code> for intervention-specific marginal
means.
</td>
</tr>
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="...">…</code>
</td>
<td>
Additional arguments passed to <code>forrest::forrest()</code> (e.g.
<code>title</code>, <code>stripe</code>, <code>dodge</code>,
<code>cols</code>, <code>widths</code>, <code>theme</code>).
</td>
</tr>
</table>

## Value

Invisibly returns <code>x</code>.

## See Also

<code>contrast()</code>, <code>print.causatr_result()</code>,
<code>forrest::forrest()</code>
