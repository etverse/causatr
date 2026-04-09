

# Compute causal contrasts from a fitted model

[**Source code**](https://github.com/etverse/causatr/tree/main/R/contrast.R#L199)

## Description

Standardises outcome model predictions under each named intervention and
reports pairwise causal contrasts (differences, ratios, or odds ratios)
with uncertainty estimates.

For all three estimation methods, <code>contrast()</code> computes
E\[Y^a\] by setting each individual’s treatment to the intervened value
and averaging the fitted outcome model predictions over the target
population. The target population is controlled by <code>estimand</code>
(or <code>subset</code> for subgroup effects). The methods differ in how
the outcome model was fitted and how the variance is estimated:

<ul>
<li>

<code>“gcomp”</code>: standard <code>glm</code>/<code>gam</code> on the
full data; sandwich SE via stacked estimating equations.

</li>
<li>

<code>“ipw”</code>: <code>glm_weightit()</code> fit weighted for the
target estimand (ATE/ATT); M-estimation sandwich SE that accounts for
weight estimation uncertainty.

</li>
<li>

<code>“matching”</code>: <code>glm()</code> on the matched sample with
match weights; SE via cluster-robust sandwich
(<code>sandwich::vcovCL(vcov = ~subclass)</code>) to account for pair
membership.

</li>
</ul>

## Usage

<pre><code class='language-R'>contrast(
  fit,
  interventions,
  type = c("difference", "ratio", "or"),
  estimand = NULL,
  subset = NULL,
  reference = NULL,
  ci_method = c("sandwich", "bootstrap"),
  n_boot = 500L,
  conf_level = 0.95,
  by = NULL,
  parallel = c("no", "multicore", "snow"),
  ncpus = getOption("boot.ncpus", 1L)
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
<code id="interventions">interventions</code>
</td>
<td>

A named list of interventions. Each element must be one of:

<ul>
<li>

A <code>causatr_intervention</code> object created by
<code>static()</code>, <code>shift()</code>, <code>dynamic()</code>,
<code>scale()</code>, <code>threshold()</code>, or <code>ipsi()</code>.

</li>
<li>

<code>NULL</code>, meaning the natural course (observed treatment values
are used as-is). The natural course is the standard reference for
modified treatment policies on continuous treatments
(e.g. <code>shift(-10)</code> vs <code>NULL</code>; Díaz et al. 2023)
and for longitudinal dynamic regimes (Hernán & Robins Ch. 21). For
binary treatments, the natural-course marginal mean equals E\[Y\] under
the observed treatment mechanism (by consistency); use
<code>static(1)</code> vs <code>static(0)</code> for the conventional
ATE. Supported for all methods.

</li>
<li>

A named list of <code>causatr_intervention</code> objects, one per
treatment variable, for multivariate (joint) treatments. Each sub-list
must name every treatment variable specified in <code>causat()</code>.

</li>
</ul>
<strong>Note:</strong> Non-static interventions (<code>shift()</code>,
<code>scale()</code>, <code>threshold()</code>, <code>dynamic()</code>,
<code>ipsi()</code>) are only supported for <code>method =
“gcomp”</code>. For IPW and matching, the weights/matched sets were
estimated under the observed treatment distribution; applying a
different regime requires re-calling <code>causat()</code> with updated
data.
</td>
</tr>
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="type">type</code>
</td>
<td>
Character. The contrast scale: <code>“difference”</code> (default),
<code>“ratio”</code>, or <code>“or”</code> (odds ratio). All pairwise
contrasts are reported.
</td>
</tr>
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="estimand">estimand</code>
</td>
<td>
Character or <code>NULL</code>. The target estimand: <code>“ATE”</code>,
<code>“ATT”</code>, or <code>“ATC”</code>. For <code>method =
“gcomp”</code>, overrides the estimand stored in <code>fit</code>
(allowing one fit to produce multiple estimands). For <code>method =
“ipw”</code> or <code>“matching”</code>, must match the estimand used at
fitting time — changing it aborts with an informative message. If
<code>NULL</code>, defaults to <code>fit$estimand</code>. Mutually
exclusive with <code>subset</code>.
</td>
</tr>
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="subset">subset</code>
</td>
<td>
A quoted expression (<code>quote(…)</code>) defining a subgroup to
average over instead of an estimand. Evaluated in the context of the
fitted data. Works for any treatment type. For example, <code>subset =
quote(age \> 50)</code> or <code>subset = quote(A == 1)</code> (the
latter is equivalent to <code>estimand = “ATT”</code> for binary
treatments). Mutually exclusive with <code>estimand</code>.
</td>
</tr>
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="reference">reference</code>
</td>
<td>
Character or <code>NULL</code>. Name of the reference intervention (the
denominator/subtracted value for pairwise contrasts). Default: the first
intervention in the list. Only relevant when <code>type =
“difference”</code> or <code>“ratio”</code> and there are more than two
interventions. For categorical treatments, use this to specify the
reference level.
</td>
</tr>
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="ci_method">ci_method</code>
</td>
<td>
Character. The variance/CI method: <code>“sandwich”</code> (default) or
<code>“bootstrap”</code>.
</td>
</tr>
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="n_boot">n_boot</code>
</td>
<td>
Integer. Number of bootstrap replications when <code>ci_method =
“bootstrap”</code>. Default <code>500</code>.
</td>
</tr>
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="conf_level">conf_level</code>
</td>
<td>
Numeric. Confidence level for intervals. Default <code>0.95</code>.
</td>
</tr>
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="by">by</code>
</td>
<td>
Character or <code>NULL</code>. Name of a variable to stratify estimates
by (effect modification). If provided, E\[Y^a\] is computed within each
level of <code>by</code>.
</td>
</tr>
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="parallel">parallel</code>
</td>
<td>
Character. Parallelisation backend for bootstrap (<code>ci_method =
“bootstrap”</code> only). <code>“no”</code> (default) runs sequentially;
<code>“multicore”</code> uses forked processes via
<code>parallel::mclapply()</code> (Unix only); <code>“snow”</code> uses
socket clusters via <code>parallel::parLapply()</code> (cross-platform).
Passed directly to <code>boot::boot()</code>. Ignored when
<code>ci_method = “sandwich”</code>.
</td>
</tr>
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="ncpus">ncpus</code>
</td>
<td>
Integer. Number of CPU cores for parallel bootstrap. Default
<code>getOption(“boot.ncpus”, 1L)</code>. Passed directly to
<code>boot::boot()</code>.
</td>
</tr>
</table>

## Details

<h4>
Estimands and standardisation
</h4>

For each intervention <code>a</code>, <code>contrast()</code> evaluates
the intervention function on the target rows to obtain the intervened
treatment vector <code>a(Lᵢ)</code>, then computes:

<pre>E\[Y^a\] = (1/|S|) Σᵢ∈S Ê\[Y | A = a(Lᵢ), Lᵢ\]
</pre>

where <code>S</code> is the set of rows determined by the estimand:

<ul>
<li>

<code>“ATE”</code>: all rows (full population)

</li>
<li>

<code>“ATT”</code>: rows where <code>A == 1</code> (observed treated)

</li>
<li>

<code>“ATC”</code>: rows where <code>A == 0</code> (observed controls)

</li>
<li>

<code>subset</code>: rows satisfying the quoted expression

</li>
</ul>

For <code>“gcomp”</code>, one fit supports multiple estimands because
the outcome model is the same — only the rows averaged over change. For
<code>“ipw”</code> and <code>“matching”</code>, the estimand is baked
into the weights/matching and cannot be changed after fitting.

<h4>
Estimand applicability
</h4>
<table>
<tr>
<td style="text-align: left;">
Treatment type
</td>
<td style="text-align: left;">
ATE
</td>
<td style="text-align: left;">
ATT/ATC
</td>
<td style="text-align: left;">
subset
</td>
</tr>
<tr>
<td style="text-align: left;">
Binary point
</td>
<td style="text-align: left;">
Yes
</td>
<td style="text-align: left;">
Yes
</td>
<td style="text-align: left;">
Yes
</td>
</tr>
<tr>
<td style="text-align: left;">
Continuous point
</td>
<td style="text-align: left;">
Yes
</td>
<td style="text-align: left;">
No (abort)
</td>
<td style="text-align: left;">
Yes
</td>
</tr>
<tr>
<td style="text-align: left;">
Categorical point
</td>
<td style="text-align: left;">
Yes
</td>
<td style="text-align: left;">
No (abort)
</td>
<td style="text-align: left;">
Yes
</td>
</tr>
<tr>
<td style="text-align: left;">
Multivariate point
</td>
<td style="text-align: left;">
Yes
</td>
<td style="text-align: left;">
No (abort)
</td>
<td style="text-align: left;">
Yes
</td>
</tr>
<tr>
<td style="text-align: left;">
Longitudinal
</td>
<td style="text-align: left;">
Yes
</td>
<td style="text-align: left;">
No (abort)
</td>
<td style="text-align: left;">
Yes
</td>
</tr>
</table>
<h4>
Treatment types and intervention support
</h4>
<table>
<tr>
<td style="text-align: left;">
Method
</td>
<td style="text-align: left;">
Treatment types
</td>
<td style="text-align: left;">
Intervention types
</td>
</tr>
<tr>
<td style="text-align: left;">
<code>“gcomp”</code>
</td>
<td style="text-align: left;">
binary, categorical, continuous, multivariate
</td>
<td style="text-align: left;">
all
</td>
</tr>
<tr>
<td style="text-align: left;">
<code>“ipw”</code>
</td>
<td style="text-align: left;">
binary, categorical, continuous (WeightIt)
</td>
<td style="text-align: left;">
<code>static()</code> only
</td>
</tr>
<tr>
<td style="text-align: left;">
<code>“matching”</code>
</td>
<td style="text-align: left;">
binary, categorical, continuous
</td>
<td style="text-align: left;">
<code>static()</code> only
</td>
</tr>
</table>
<h4>
Variance estimation
</h4>
<ul>
<li>

<code>“sandwich”</code>: Stacked estimating equations propagate outcome
model uncertainty through to the marginal mean.

</li>
<li>

<code>“bootstrap”</code>: Resamples individuals (entire pipeline
refitted <code>n_boot</code> times). Respects cluster structure for
longitudinal data.

</li>
<li>

<code>“delta”</code>: Explicit delta method for ratio/OR contrasts
(applied internally when <code>ci_method = “sandwich”</code>).

</li>
</ul>

## Value

A <code>causatr_result</code> object with slots:

<dl>
<dt>
<code>estimates</code>
</dt>
<dd>
data.table with one row per intervention: <code>intervention</code>,
<code>estimate</code>, <code>se</code>, <code>ci_lower</code>,
<code>ci_upper</code>.
</dd>
<dt>
<code>contrasts</code>
</dt>
<dd>
data.table with one row per pairwise comparison:
<code>comparison</code>, <code>estimate</code>, <code>se</code>,
<code>ci_lower</code>, <code>ci_upper</code>.
</dd>
<dt>
<code>type</code>
</dt>
<dd>
Contrast scale.
</dd>
<dt>
<code>estimand</code>
</dt>
<dd>
<code>“ATE”</code>, <code>“ATT”</code>, <code>“ATC”</code>, or
<code>“subset”</code>.
</dd>
<dt>
<code>ci_method</code>
</dt>
<dd>
Inference method.
</dd>
<dt>
<code>reference</code>
</dt>
<dd>
Name of the reference intervention.
</dd>
<dt>
<code>interventions</code>
</dt>
<dd>
The intervention list.
</dd>
<dt>
<code>n</code>
</dt>
<dd>
Number of individuals averaged over.
</dd>
<dt>
<code>method</code>
</dt>
<dd>
Estimation method from the fit.
</dd>
<dt>
<code>vcov</code>
</dt>
<dd>
Full variance-covariance matrix for all E\[Y^a\].
</dd>
<dt>
<code>call</code>
</dt>
<dd>
The original call.
</dd>
</dl>

## References

Hernán MA, Robins JM (2025). <em>Causal Inference: What If</em>. Chapman
& Hall/CRC. Chapters 12–13.

Greifer N (2024). WeightIt: Weighting Methods for Covariate Balancing.
<a href="https://ngreifer.github.io/WeightIt/">https://ngreifer.github.io/WeightIt/</a>

Imai K, King G, Stuart EA (2011). Misunderstandings between
experimentalists and observationalists about causal inference.
<em>Journal of the Royal Statistical Society</em> Series A 171:481–502.

Zivich PN, Ross RK, Shook-Sa BE, Cole SR, Edwards JK (2024). Empirical
sandwich variance estimator for iterated conditional expectation
g-computation. <em>Statistics in Medicine</em> 43:5562–5572.

## See Also

<code>causat()</code>, <code>static()</code>, <code>shift()</code>,
<code>dynamic()</code>, <code>coef.causatr_result()</code>,
<code>confint.causatr_result()</code>

## Examples

``` r
library("causatr")

data("nhefs", package = "causatr")
fit <- causat(nhefs, outcome = "wt82_71", treatment = "qsmk",
              confounders = ~ sex + age + wt71)

# ATE: mean difference with sandwich SE
result <- contrast(fit,
  interventions = list(quit = static(1), continue = static(0)),
  type = "difference"
)

# ATT: override estimand in contrast() (gcomp only)
result_att <- contrast(fit,
  interventions = list(quit = static(1), continue = static(0)),
  estimand = "ATT"
)

# Subgroup effect: age > 50
result_sub <- contrast(fit,
  interventions = list(quit = static(1), continue = static(0)),
  subset = quote(age > 50)
)

# Continuous treatment: shift with NULL reference
fit_cont <- causat(nhefs, outcome = "wt82_71",
                   treatment = "smokeintensity",
                   confounders = ~ sex + age + wt71)
contrast(fit_cont,
  interventions = list(reduce10 = shift(-10), observed = NULL),
  type = "difference"
)

# Categorical treatment: three arms, reference = "radio"
contrast(fit_cat,
  interventions = list(chemo = static("A"), radio = static("B"),
                       combo = static("C")),
  type = "difference",
  reference = "radio"
)
```
