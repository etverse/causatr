

# Fit a causal model

[**Source code**](https://github.com/etverse/causatr/tree/main/R/causat.R#L221)

## Description

Prepares the causal estimation pipeline for a given method. For
<code>“gcomp”</code>, fits the conditional outcome model E\[Y | A, L\]
that will be used by <code>contrast()</code> for standardisation. For
<code>“ipw”</code>, estimates propensity-score-based weights (via
<code>WeightIt::weightit()</code>) that will be used for weighted
estimation in <code>contrast()</code>. For <code>“matching”</code>,
creates matched sets (via <code>MatchIt::matchit()</code>) that will be
used for matched estimation in <code>contrast()</code>.

For longitudinal data (<code>id</code> and <code>time</code> provided),
<code>“gcomp”</code> uses ICE g-computation (Zivich et al., 2024):
outcome models are fitted at each time point via backward iteration.

## Usage

<pre><code class='language-R'>causat(
  data,
  outcome,
  treatment,
  confounders,
  confounders_tv = NULL,
  method = c("gcomp", "ipw", "matching"),
  family = "gaussian",
  estimand = c("ATE", "ATT", "ATC"),
  id = NULL,
  time = NULL,
  censoring = NULL,
  history = 1L,
  numerator = NULL,
  weights = NULL,
  model_fn = stats::glm,
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
A data frame or data.table.
</td>
</tr>
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="outcome">outcome</code>
</td>
<td>
Character. Name of the outcome variable.
</td>
</tr>
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="treatment">treatment</code>
</td>
<td>
Character scalar or character vector. Name(s) of the treatment
variable(s). Pass a character vector for multivariate (joint)
treatments, e.g. <code>treatment = c(“A1”, “A2”)</code>.
</td>
</tr>
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="confounders">confounders</code>
</td>
<td>
A one-sided formula specifying baseline (time-invariant) confounders,
e.g. <code>~ L1 + L2</code>. Interactions and transformations are
allowed, e.g. <code>~ L1 \* L2 + splines::ns(age, 4)</code>. For
longitudinal data, these confounders are constant within each individual
(measured at baseline) and enter every time-step outcome model.
</td>
</tr>
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="confounders_tv">confounders_tv</code>
</td>
<td>
A one-sided formula or <code>NULL</code>. Time-varying confounders for
longitudinal data (e.g. <code>~ CD4 + viral_load</code>). These change
over time within individuals and enter the outcome model at each ICE
step alongside their lagged values (controlled by <code>history</code>).
Ignored for point treatments. If <code>NULL</code>, no time-varying
confounders are used.
</td>
</tr>
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="method">method</code>
</td>
<td>
Character. Estimation method: <code>“gcomp”</code> (default),
<code>“ipw”</code>, or <code>“matching”</code>. IPW requires the
<code>WeightIt</code> package; matching requires the
<code>MatchIt</code> package.
</td>
</tr>
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="family">family</code>
</td>
<td>
Character or family object. The outcome model family for
<code>“gcomp”</code> (e.g. <code>“gaussian”</code>,
<code>“binomial”</code>). Passed to <code>glm()</code> or
<code>mgcv::gam()</code>. Ignored for <code>“ipw”</code> and
<code>“matching”</code>.
</td>
</tr>
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="estimand">estimand</code>
</td>
<td>
Character. The target estimand: <code>“ATE”</code> (default),
<code>“ATT”</code>, or <code>“ATC”</code>. <code>“ATT”</code> and
<code>“ATC”</code> are only defined for binary point treatments. For
<code>“gcomp”</code>, the estimand stored here is used as the default in
<code>contrast()</code> but can be overridden there. For
<code>“ipw”</code> and <code>“matching”</code>, the estimand is fixed at
fitting time because it determines the weights or match direction; it
cannot be changed in <code>contrast()</code>.
</td>
</tr>
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="id">id</code>
</td>
<td>
Character or <code>NULL</code>. Name of the individual ID variable. Must
be provided together with <code>time</code> for longitudinal data.
</td>
</tr>
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="time">time</code>
</td>
<td>
Character or <code>NULL</code>. Name of the time variable. Must be
provided together with <code>id</code> for longitudinal data.
</td>
</tr>
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="censoring">censoring</code>
</td>
<td>
Character or <code>NULL</code>. Name of the censoring indicator variable
(1 = censored, 0 = uncensored). For longitudinal data, censoring is
time-varying: <code>C_k = 1</code> means the individual dropped out at
time <code>k</code>.
</td>
</tr>
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="history">history</code>
</td>
<td>
Positive integer or <code>Inf</code>. Markov order for longitudinal
models: how many lagged time points of treatment and time-varying
confounders to include in each ICE outcome model. Default <code>1</code>
(one lag). <code>Inf</code> includes the full history. Ignored for point
treatments.
</td>
</tr>
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="numerator">numerator</code>
</td>
<td>
A one-sided formula or <code>NULL</code>. Numerator formula for
stabilized IPW weights in longitudinal models. Defaults to baseline
confounders only (no time-varying confounders), which gives the standard
stabilized weights. Only relevant for <code>method = “ipw”</code> with
longitudinal data.
</td>
</tr>
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="weights">weights</code>
</td>
<td>
Numeric vector or <code>NULL</code>. Pre-computed observation weights
(e.g. survey weights or externally computed IPCW). For
<code>“gcomp”</code>, passed to <code>glm()</code>. For
<code>“ipw”</code>, multiplied with the estimated propensity weights.
</td>
</tr>
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="model_fn">model_fn</code>
</td>
<td>
Function. The fitting function for the outcome model in g-computation.
Must accept <code style="white-space: pre;">(formula, data, family,
weights, …)</code>. Default <code>stats::glm</code>; pass
<code>mgcv::gam</code> for GAMs, <code>MASS::glm.nb</code> for
negative-binomial, etc. Ignored for <code>“ipw”</code> and
<code>“matching”</code>.
</td>
</tr>
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="...">…</code>
</td>
<td>
Additional arguments passed to the underlying estimation function:
<code>WeightIt::weightit()</code> for <code>method = “ipw”</code> (e.g.
<code>method = “glm”</code>); <code>MatchIt::matchit()</code> for
<code>method = “matching”</code> (e.g. <code>method = “nearest”</code>,
<code>ratio = 1</code>).
</td>
</tr>
</table>

## Details

<h4>
G-computation (<code>method = “gcomp”</code>)
</h4>

Fits E\[Y | A, L\] using <code>glm()</code> (or <code>mgcv::gam()</code>
if the formula contains <code>s()</code> or <code>te()</code> terms).
<code>contrast()</code> standardises over the target population
(controlled by <code>estimand</code>) to obtain E\[Y^a\] under each
intervention. For <code>“gcomp”</code>, the estimand can be overridden
in <code>contrast()</code> — fit once, contrast with multiple estimands.

For longitudinal data, uses ICE g-computation (Zivich et al., 2024):
outcome models are fitted one per time point via backward iteration,
conditioning on baseline confounders (<code>confounders</code>),
time-varying confounders (<code>confounders_tv</code>), and their lags
up to <code>history</code> steps.

<h4>
IPW (<code>method = “ipw”</code>)
</h4>

Calls <code>WeightIt::weightit()</code> to estimate
propensity-score-based weights for the desired <code>estimand</code>,
then fits a weighted outcome model via
<code>WeightIt::glm_weightit()</code>. The estimand is fixed at fitting
time. Supports binary, categorical, and continuous treatments (anything
<code>WeightIt</code> supports). Only <code>static()</code>
interventions are supported in <code>contrast()</code>, because the
weights were estimated under the observed treatment distribution.

For longitudinal data, calls <code>WeightIt::weightitMSM()</code>. The
denominator model at each time <code>k</code> includes baseline
confounders, concurrent time-varying confounders, and lagged treatment.
The numerator model includes only baseline confounders and lagged
treatment (standard stabilized weights), unless overridden via
<code>numerator</code>.

<h4>
Matching (<code>method = “matching”</code>)
</h4>

Calls <code>MatchIt::matchit()</code> to create matched sets. The
estimand is fixed at fitting time. Only <code>static()</code>
interventions are supported in <code>contrast()</code>.

<h4>
Estimands
</h4>
<table>
<tr>
<td style="text-align: left;">
Estimand
</td>
<td style="text-align: left;">
Population averaged over
</td>
<td style="text-align: left;">
Applicability
</td>
</tr>
<tr>
<td style="text-align: left;">
<code>“ATE”</code>
</td>
<td style="text-align: left;">
All individuals
</td>
<td style="text-align: left;">
Always
</td>
</tr>
<tr>
<td style="text-align: left;">
<code>“ATT”</code>
</td>
<td style="text-align: left;">
Observed treated (<code>A = 1</code>)
</td>
<td style="text-align: left;">
Binary point treatment only
</td>
</tr>
<tr>
<td style="text-align: left;">
<code>“ATC”</code>
</td>
<td style="text-align: left;">
Observed controls (<code>A = 0</code>)
</td>
<td style="text-align: left;">
Binary point treatment only
</td>
</tr>
</table>

For continuous, categorical, multivariate, or longitudinal treatments,
use <code>estimand = “ATE”</code> and pass a <code>subset</code>
expression to <code>contrast()</code> for subgroup effects.

<h4>
Identifiability assumptions
</h4>

All methods assume: (1) exchangeability (no unmeasured confounding given
L), (2) positivity (every individual has positive probability of each
treatment value given L), (3) consistency (the observed outcome under
the observed treatment equals the potential outcome). Positivity is
checked automatically and a warning is issued if near-violations are
detected.

## Value

A <code>causatr_fit</code> object with slots:

<dl>
<dt>
<code>model</code>
</dt>
<dd>
Fitted model object(s): <code>glm</code>/<code>gam</code> for
<code>“gcomp”</code>; <code>NULL</code> for <code>“ipw”</code> and
<code>“matching”</code> (weights/matched data are stored in
<code>weights_obj</code> / <code>match_obj</code> instead).
</dd>
<dt>
<code>data</code>
</dt>
<dd>
data.table used for fitting.
</dd>
<dt>
<code>treatment</code>, <code>outcome</code>, <code>confounders</code>,
<code>confounders_tv</code>, <code>family</code>
</dt>
<dd>
Model spec.
</dd>
<dt>
<code>method</code>
</dt>
<dd>
<code>“gcomp”</code>, <code>“ipw”</code>, or <code>“matching”</code>.
</dd>
<dt>
<code>type</code>
</dt>
<dd>
<code>“point”</code> or <code>“longitudinal”</code>.
</dd>
<dt>
<code>estimand</code>
</dt>
<dd>
<code>“ATE”</code>, <code>“ATT”</code>, or <code>“ATC”</code>.
</dd>
<dt>
<code>id</code>, <code>time</code>, <code>censoring</code>
</dt>
<dd>
Longitudinal identifiers.
</dd>
<dt>
<code>history</code>
</dt>
<dd>
Markov order for longitudinal ICE models.
</dd>
<dt>
<code>numerator</code>
</dt>
<dd>
Numerator formula for longitudinal IPW.
</dd>
<dt>
<code>weights_obj</code>
</dt>
<dd>
<code>weightit</code> object (IPW only).
</dd>
<dt>
<code>match_obj</code>
</dt>
<dd>
<code>matchit</code> object (matching only).
</dd>
<dt>
<code>call</code>
</dt>
<dd>
The original call.
</dd>
<dt>
<code>details</code>
</dt>
<dd>
Internal diagnostics list.
</dd>
</dl>

## References

Hernán MA, Robins JM (2025). <em>Causal Inference: What If</em>. Chapman
& Hall/CRC. Chapters 12 (IPW), 13 (g-formula), 15 (matching), 21 (ICE).

Zivich PN, Ross RK, Shook-Sa BE, Cole SR, Edwards JK (2024). Empirical
sandwich variance estimator for iterated conditional expectation
g-computation. <em>Statistics in Medicine</em> 43:5562–5572.

## See Also

<code>contrast()</code>, <code>diagnose()</code>,
<code>causat_survival()</code>, <code>causat_mice()</code>,
<code>static()</code>, <code>data.table::shift()</code>,
<code>dynamic()</code>

## Examples

``` r
library("causatr")

data("nhefs", package = "causatr")

# Point treatment, g-computation (default)
fit <- causat(
  nhefs,
  outcome = "wt82_71",
  treatment = "qsmk",
  confounders = ~ sex + age + race + education +
    smokeintensity + smokeyrs + exercise + active + wt71
)

# ATT via g-computation (override estimand in contrast())
fit_att <- causat(
  nhefs,
  outcome = "wt82_71",
  treatment = "qsmk",
  confounders = ~ sex + age + race + education +
    smokeintensity + smokeyrs + exercise + active + wt71,
  estimand = "ATT"
)

# IPW (requires WeightIt) — estimand fixed at fit time
fit_ipw <- causat(
  nhefs,
  outcome = "wt82_71",
  treatment = "qsmk",
  confounders = ~ sex + age + race + education +
    smokeintensity + smokeyrs + exercise + active + wt71,
  method = "ipw",
  estimand = "ATE"
)

# Matching (requires MatchIt)
fit_match <- causat(
  nhefs,
  outcome = "wt82_71",
  treatment = "qsmk",
  confounders = ~ sex + age + race + education +
    smokeintensity + smokeyrs + exercise + active + wt71,
  method = "matching",
  estimand = "ATT"
)

# Longitudinal with time-varying confounders
fit_long <- causat(
  data_long,
  outcome = "Y",
  treatment = "A",
  confounders = ~ sex + race + baseline_age,
  confounders_tv = ~ CD4 + viral_load,
  id = "id",
  time = "time",
  history = 1L
)

# Multivariate treatment
fit_multi <- causat(
  data,
  outcome = "Y",
  treatment = c("A1", "A2"),
  confounders = ~ L1 + L2
)
```
