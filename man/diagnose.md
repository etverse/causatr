

# Diagnostics for a fitted causal model

[**Source code**](https://github.com/etverse/causatr/tree/main/R/diagnose.R#L76)

## Description

Computes diagnostics appropriate to the estimation method:

<ul>
<li>

<strong>All methods</strong>: positivity checks (flags covariate strata
where the probability of treatment is near 0 or 1).

</li>
<li>

<code>“ipw”</code>: covariate balance before and after weighting (via
<code>cobalt</code>), weight distribution summary (mean, SD, max,
effective sample size).

</li>
<li>

<code>“matching”</code>: covariate balance before and after matching
(via <code>cobalt</code>), match quality summary (% matched, caliper
info).

</li>
<li>

<code>“gcomp”</code>: unadjusted covariate imbalance between treatment
groups.

</li>
</ul>

## Usage

<pre><code class='language-R'>diagnose(
  fit,
  stats = c("m", "v"),
  thresholds = c(m = 0.1),
  ps_bounds = c(0.025, 0.975)
)
</code></pre>

## Arguments

<table role="presentation">
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="fit">fit</code>
</td>
<td>
A <code>causatr_fit</code> object returned by <code>causat()</code>.
</td>
</tr>
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="stats">stats</code>
</td>
<td>
Character vector. Balance statistics to compute. Passed to
<code>cobalt::bal.tab()</code>. For binary treatments, valid options
include <code>“m”</code> (standardised mean differences),
<code>“v”</code> (variance ratios), and <code>“ks”</code>
(Kolmogorov-Smirnov). Default <code>c(“m”, “v”)</code>.
</td>
</tr>
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="thresholds">thresholds</code>
</td>
<td>
Named numeric vector. Balance thresholds for flagging imbalance,
e.g. <code>c(m = 0.1, v = 2)</code>. Default <code>c(m = 0.1)</code>.
</td>
</tr>
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="ps_bounds">ps_bounds</code>
</td>
<td>
Numeric vector of length 2. Lower and upper bounds for flagging
positivity violations. Default <code>c(0.025, 0.975)</code>.
</td>
</tr>
</table>

## Details

<h4>
Positivity
</h4>

For binary treatment, fits a logistic regression of the treatment on the
confounders and flags individuals whose estimated propensity score falls
outside <code>ps_bounds</code>. The returned <code>positivity</code>
table summarises the propensity score distribution and the
number/fraction of near-violations.

<h4>
Balance (IPW and matching)
</h4>

If the <code>cobalt</code> package is installed, balance is computed via
<code>cobalt::bal.tab()</code> on the internal <code>weightit</code> or
<code>matchit</code> object. This provides standardised mean differences
(SMD), variance ratios, and KS statistics before and after adjustment.
If <code>cobalt</code> is not installed, a simpler data.table-based SMD
comparison is returned.

<h4>
Weight distribution (IPW only)
</h4>

Summarises the IPW weights: mean, SD, min, max, and the effective sample
size (ESS) for the treated and control groups.

<h4>
Match quality (matching only)
</h4>

Reports the number matched, number discarded, and the fraction of the
original sample retained.

## Value

A <code>causatr_diag</code> object with slots:

<dl>
<dt>
<code>balance</code>
</dt>
<dd>
<code>cobalt::bal.tab</code> object (if cobalt installed) or a
data.table of SMDs. Covariate balance summary.
</dd>
<dt>
<code>positivity</code>
</dt>
<dd>
data.table: propensity score summary and count of near-violations.
</dd>
<dt>
<code>weights</code>
</dt>
<dd>
data.table or <code>NULL</code>: weight distribution summary (IPW only).
</dd>
<dt>
<code>match_quality</code>
</dt>
<dd>
data.table or <code>NULL</code>: match quality summary (matching only).
</dd>
<dt>
<code>method</code>
</dt>
<dd>
Character: the estimation method.
</dd>
<dt>
<code>fit</code>
</dt>
<dd>
The original <code>causatr_fit</code> (stored for <code>plot()</code>).
</dd>
</dl>

## References

Greifer N (2024). cobalt: Covariate Balance Tables and Plots.
<a href="https://ngreifer.github.io/cobalt/">https://ngreifer.github.io/cobalt/</a>

## See Also

<code>causat()</code>, <code>plot.causatr_diag()</code>

## Examples

``` r
library("causatr")

data("nhefs", package = "causatr")
fit <- causat(nhefs, outcome = "wt82_71", treatment = "qsmk",
              confounders = ~ sex + age + wt71,
              method = "ipw")
diag <- diagnose(fit)
print(diag)
plot(diag)
```
