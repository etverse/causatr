

# Convert wide data to person-period (long) format

[**Source code**](https://github.com/etverse/causatr/tree/main/R/to_person_period.R#L47)

## Description

Reshapes individual-level data in wide format (one row per person, with
columns for each time-varying variable at each time point) into long
person-period format (one row per person per time interval), as required
by <code>causat()</code> for longitudinal analyses and by
<code>causat_survival()</code> for discrete-time survival analyses.

For survival data the event indicator *D*<sub>*k* + 1</sub> = 1 means
the individual experienced the event during interval k+1 (Hernán &
Robins Ch. 17).

## Usage

<pre><code class='language-R'>to_person_period(
  data,
  id,
  time_varying,
  time_invariant = character(0),
  time_name = "time"
)
</code></pre>

## Arguments

<table role="presentation">
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="data">data</code>
</td>
<td>
A data frame or data.table in wide format.
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
<code id="time_varying">time_varying</code>
</td>
<td>
Named list specifying how to reshape time-varying columns. Each element
is a character vector of column names for one variable across time, in
time order. For example: <code>list(A = c(“A0”, “A1”, “A2”), L = c(“L0”,
“L1”, “L2”))</code>.
</td>
</tr>
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="time_invariant">time_invariant</code>
</td>
<td>
Character vector of column names to carry forward unchanged (baseline
covariates, outcome, etc.).
</td>
</tr>
<tr>
<td style="white-space: collapse; font-family: monospace; vertical-align: top">
<code id="time_name">time_name</code>
</td>
<td>
Character. Name for the new time column. Default <code>“time”</code>.
</td>
</tr>
</table>

## Value

A data.table in person-period format with one row per individual per
time point, sorted by <code>id</code> and <code>time_name</code>.

## See Also

<code>causat()</code>, <code>causat_survival()</code>

## Examples

``` r
library("causatr")

wide <- data.table::data.table(
  id  = 1:3,
  sex = c(0, 1, 0),
  A0  = c(1, 0, 1),
  A1  = c(1, 1, 0),
  L0  = c(5, 3, 7),
  L1  = c(4, 6, 8),
  Y   = c(0, 1, 0)
)
long <- to_person_period(
  wide,
  id = "id",
  time_varying = list(A = c("A0", "A1"), L = c("L0", "L1")),
  time_invariant = c("sex", "Y")
)
```
