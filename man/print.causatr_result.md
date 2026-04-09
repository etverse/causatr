

# Print a causatr result

[**Source code**](https://github.com/etverse/causatr/tree/main/R/print.R#L75)

## Description

Displays the estimation method, contrast type, CI method, sample size,
intervention-specific marginal means, and pairwise contrasts (with SEs
and confidence intervals).

## Usage

<pre><code class='language-R'>## S3 method for class 'causatr_result'
print(x, ...)
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

Invisibly returns <code>x</code>.

## See Also

<code>summary.causatr_result()</code>, <code>contrast()</code>
