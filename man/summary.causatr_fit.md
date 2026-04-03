

# Summarise a causatr fit

[**Source code**](https://github.com/etverse/causatr/tree/main/R/summary.R#L13)

## Description

Provides a detailed summary of a causatr_fit object, including the
estimation method, outcome family, confounder formula, and the
underlying model (outcome model, propensity weights, or matching
object).

## Usage

<pre><code class='language-R'>## S3 method for class 'causatr_fit'
summary(object, ...)
</code></pre>

## Arguments

<table role="presentation">
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="object">object</code>
</td>
<td>
A <code>causatr_fit</code> object.
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

<code>print.causatr_fit()</code>, <code>causat()</code>
