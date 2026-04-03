

# Plot a causatr result

[**Source code**](https://github.com/etverse/causatr/tree/main/R/plot.R#L15)

## Description

Produces a forest plot of intervention means and pairwise contrasts from
a <code>causatr_result</code> object (via <code>tinyplot</code> if
available, otherwise base graphics).

## Usage

<pre><code class='language-R'>## S3 method for class 'causatr_result'
plot(x, ...)
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
Additional arguments (currently unused).
</td>
</tr>
</table>

## Value

Invisibly returns <code>x</code>.

## See Also

<code>contrast()</code>, <code>print.causatr_result()</code>
