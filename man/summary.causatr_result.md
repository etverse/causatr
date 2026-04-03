

# Summarise a causatr result

[**Source code**](https://github.com/etverse/causatr/tree/main/R/summary.R#L69)

## Description

Displays intervention-specific marginal means and pairwise contrasts
with standard errors and confidence intervals from a causatr_result
object.

## Usage

<pre><code class='language-R'>## S3 method for class 'causatr_result'
summary(object, ...)
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

Invisibly returns <code>object</code>.

## See Also

<code>print.causatr_result()</code>, <code>contrast()</code>
