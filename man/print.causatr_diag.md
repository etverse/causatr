

# Print causatr diagnostics

[**Source code**](https://github.com/etverse/causatr/tree/main/R/print.R#L127)

## Description

Displays positivity summaries, covariate balance tables, weight
distributions (IPW), and match quality metrics (matching) from a
causatr_diag object.

## Usage

<pre><code class='language-R'>## S3 method for class 'causatr_diag'
print(x, ...)
</code></pre>

## Arguments

<table role="presentation">
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="x">x</code>
</td>
<td>
A <code>causatr_diag</code> object.
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

Invisibly returns <code>x</code>.

## See Also

<code>summary.causatr_diag()</code>, <code>diagnose()</code>
