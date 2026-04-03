

# Fit causal models across multiply-imputed datasets

[**Source code**](https://github.com/etverse/causatr/tree/main/R/causat_mice.R#L70)

## Description

Applies <code>causat()</code> and <code>contrast()</code> across all
imputed datasets in a <code>mids</code> object (from
<code>mice::mice()</code>), then pools point estimates and variances
using Rubin’s rules. This handles missing treatment values (or any other
missing data) via multiple imputation.

Rubin’s rules pool m point estimates and their within-imputation
variances:

<ul>
<li>

<strong>Pooled estimate</strong>: mean of per-imputation estimates.

</li>
<li>

<strong>Total variance</strong>: within-imputation variance +
between-imputation variance + between-imputation correction.

</li>
<li>

<strong>CIs</strong>: normal approximation on the pooled estimate.

</li>
</ul>

## Usage

<pre><code class='language-R'>causat_mice(
  imp,
  outcome,
  treatment,
  confounders,
  interventions,
  method = "gcomp",
  family = "gaussian",
  estimand = "ATE",
  type = "difference",
  ci_method = "sandwich",
  conf_level = 0.95,
  ...
)
</code></pre>

## Arguments

<table role="presentation">
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="imp">imp</code>
</td>
<td>
A <code>mids</code> object returned by <code>mice::mice()</code>.
</td>
</tr>
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="outcome">outcome</code>
</td>
<td>
Character. Passed to <code>causat()</code>.
</td>
</tr>
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="treatment">treatment</code>
</td>
<td>
Character. Passed to <code>causat()</code>.
</td>
</tr>
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="confounders">confounders</code>
</td>
<td>
A one-sided formula. Passed to <code>causat()</code>.
</td>
</tr>
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="interventions">interventions</code>
</td>
<td>
A named list of interventions. Passed to <code>contrast()</code>.
</td>
</tr>
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="method">method</code>
</td>
<td>
Character. Passed to <code>causat()</code>. Default
<code>“gcomp”</code>.
</td>
</tr>
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="family">family</code>
</td>
<td>
Character or family object. Passed to <code>causat()</code>.
</td>
</tr>
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="estimand">estimand</code>
</td>
<td>
Character. Passed to <code>causat()</code>.
</td>
</tr>
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="type">type</code>
</td>
<td>
Character. Contrast scale. Passed to <code>contrast()</code>.
</td>
</tr>
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="ci_method">ci_method</code>
</td>
<td>
Character. Variance method. Passed to <code>contrast()</code>.
</td>
</tr>
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="conf_level">conf_level</code>
</td>
<td>
Numeric. Confidence level. Default <code>0.95</code>.
</td>
</tr>
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="...">…</code>
</td>
<td>
Additional arguments passed to <code>causat()</code>.
</td>
</tr>
</table>

## Details

<h4>
Rubin’s rules
</h4>

Let <code style="white-space: pre;">Q̂\_i</code> be the estimate from
imputation <code>i</code> and <code>U_i</code> its variance. The pooled
estimate is <code style="white-space: pre;">Q̄ = (1/m) Σ Q̂\_i</code>. The
total variance is:

<pre>T = Ū + B + B/m
</pre>

where <code style="white-space: pre;">Ū = (1/m) Σ U_i</code>
(within-imputation variance) and <code style="white-space: pre;">B =
(1/(m-1)) Σ (Q̂\_i - Q̄)²</code> (between-imputation variance).

Degrees of freedom follow the standard Barnard–Rubin approximation.

## Value

A <code>causatr_result</code> with pooled estimates and variances
following Rubin’s rules.

## References

Rubin DB (1987). <em>Multiple Imputation for Nonresponse in
Surveys</em>. Wiley.

van Buuren S, Groothuis-Oudshoorn K (2011). mice: Multivariate
Imputation by Chained Equations in R. <em>Journal of Statistical
Software</em> 45(3):1–67.

## See Also

<code>causat()</code>, <code>contrast()</code>

## Examples

``` r
library("causatr")

library(mice)

# Step 1: impute
imp <- mice(data, m = 20, method = "pmm")

# Step 2: fit + contrast across imputations
result <- causat_mice(
  imp,
  outcome = "Y",
  treatment = "A",
  confounders = ~ L1 + L2,
  interventions = list(treat_all = static(1), treat_none = static(0)),
  type = "difference"
)
```
