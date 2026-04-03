

# Fit a causal survival model

[**Source code**](https://github.com/etverse/causatr/tree/main/R/causat_survival.R#L71)

## Description

Convenience wrapper for causal survival analysis using pooled logistic
regression as a discrete-time hazard model (Hernán & Robins Ch. 17).

<h4>
Algorithm
</h4>
<ol>
<li>

Convert data to person-period format if not already long (using
<code>to_person_period()</code>).

</li>
<li>

Fit a pooled logistic regression for the discrete hazard:
*l**o**g**i**t**P**r*\[*D*<sub>*k* + 1</sub> = 1|*s**u**r**v**i**v**e**d**t**o**k*, *A*, *L*\]
= <code>time_formula + A + confounders</code>.

</li>
<li>

For each intervention, predict individual-level hazards, compute
survival as the cumulative product S_i(k) = prod(1 - h_i(m), m \<= k),
and average across individuals.

</li>
<li>

Risk difference at time t =
(1 − *S*<sup>*a*1</sup>(*t*)) − (1 − *S*<sup>*a*0</sup>(*t*)).

</li>
</ol>

When the per-interval hazard is small (\< 0.1), the pooled logistic
model closely approximates a continuous-time Cox model (Technical Point
17.1).

## Usage

<pre><code class='language-R'>causat_survival(
  data,
  outcome,
  treatment,
  confounders,
  id,
  time,
  censoring = NULL,
  competing = NULL,
  time_formula = ~splines::ns(time, 4),
  weights = NULL,
  ...
)
</code></pre>

## Arguments

<table role="presentation">
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="data">data</code>
</td>
<td>
A data frame or data.table. Can be in wide format (one row per
individual with a time-to-event column) or long person-period format
(one row per person per time interval). If wide, the data is
auto-converted using <code>to_person_period()</code>.
</td>
</tr>
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="outcome">outcome</code>
</td>
<td>
Character. Name of the binary event indicator (1 = event occurred in
this interval, 0 = survived / censored).
</td>
</tr>
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="treatment">treatment</code>
</td>
<td>
Character. Name of the treatment variable.
</td>
</tr>
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="confounders">confounders</code>
</td>
<td>
A one-sided formula specifying confounders.
</td>
</tr>
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="id">id</code>
</td>
<td>
Character. Name of the individual ID variable.
</td>
</tr>
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="time">time</code>
</td>
<td>
Character. Name of the time variable (interval index).
</td>
</tr>
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="censoring">censoring</code>
</td>
<td>
Character or <code>NULL</code>. Name of the censoring indicator. If
provided, rows where <code>censoring == 1</code> are excluded from
fitting, and subsequent rows for that individual are also dropped.
</td>
</tr>
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="competing">competing</code>
</td>
<td>
Character or <code>NULL</code>. Name of a variable indicating the type
of competing event (for competing risks analysis).
</td>
</tr>
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="time_formula">time_formula</code>
</td>
<td>
A one-sided formula specifying how time enters the hazard model. Default
<code>~ splines::ns(time, 4)</code>. Use <code>~ factor(time)</code> for
a fully saturated (non-parametric) baseline hazard.
</td>
</tr>
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="weights">weights</code>
</td>
<td>
Numeric vector or <code>NULL</code>. Pre-computed IPCW or survey
weights.
</td>
</tr>
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="...">…</code>
</td>
<td>
Additional arguments passed to <code>glm()</code>.
</td>
</tr>
</table>

## Value

A <code>causatr_fit</code> object (with <code>type = “survival”</code>)
suitable for use with <code>contrast()</code>.

## References

Hernán MA, Robins JM (2025). <em>Causal Inference: What If</em>. Chapman
& Hall/CRC. Chapter 17.

## See Also

<code>causat()</code>, <code>contrast()</code>,
<code>to_person_period()</code>

## Examples

``` r
library("causatr")

data("nhefs", package = "causatr")
fit_surv <- causat_survival(
  nhefs,
  outcome = "death",
  treatment = "qsmk",
  confounders = ~ sex + age + race + education +
    smokeintensity + smokeyrs + exercise + active + wt71,
  id = "seqn",
  time = "year"
)
result <- contrast(fit_surv,
  interventions = list(quit = static(1), continue = static(0)),
  type = "difference"
)
```
