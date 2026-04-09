

# Summarise causatr diagnostics

[**Source code**](https://github.com/etverse/causatr/tree/main/R/summary.R#L105)

## Description

Prints positivity, balance, weight, and match quality diagnostics from a
causatr_diag object.

## Usage

<pre><code class='language-R'>## S3 method for class 'causatr_diag'
summary(object, ...)
</code></pre>

## Arguments

<table role="presentation">
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="object">object</code>
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

Invisibly returns <code>object</code>.

## See Also

<code>print.causatr_diag()</code>, <code>diagnose()</code>
