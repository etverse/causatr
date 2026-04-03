

# NHEFS: National Health and Nutrition Examination Survey Epidemiologic Follow-up Study

## Description

A subset of the NHEFS dataset used throughout Hernán & Robins (2025) to
illustrate causal effect estimation methods. Contains 1,746 cigarette
smokers aged 25–74 who completed a baseline survey (1971–75) and a
follow-up survey (1982). Of these, 117 have missing education values;
the book analyses typically use the 1,629-row subset with complete
education.

## Usage

<pre><code class='language-R'>nhefs
</code></pre>

## Format

A data.table with 1,746 rows and the following variables:

<dl>
<dt>
<code>seqn</code>
</dt>
<dd>
Individual identifier.
</dd>
<dt>
<code>qsmk</code>
</dt>
<dd>
Quit smoking between baseline and 1982 (1 = yes, 0 = no). Primary
treatment variable.
</dd>
<dt>
<code>wt82_71</code>
</dt>
<dd>
Weight change in kg between 1971 and 1982. Primary outcome variable.
Missing for 63 individuals lost to follow-up.
</dd>
<dt>
<code>sex</code>
</dt>
<dd>
Sex (0 = male, 1 = female).
</dd>
<dt>
<code>age</code>
</dt>
<dd>
Age at baseline (years).
</dd>
<dt>
<code>race</code>
</dt>
<dd>
Race (0 = white, 1 = other).
</dd>
<dt>
<code>education</code>
</dt>
<dd>
Education level (1–5, from less than high school to college or more).
</dd>
<dt>
<code>smokeintensity</code>
</dt>
<dd>
Cigarettes smoked per day at baseline.
</dd>
<dt>
<code>smokeyrs</code>
</dt>
<dd>
Years of smoking at baseline.
</dd>
<dt>
<code>exercise</code>
</dt>
<dd>
Recreational exercise (0 = much, 1 = moderate, 2 = little/none).
</dd>
<dt>
<code>active</code>
</dt>
<dd>
Activity level at work (0 = very active, 1 = moderately active, 2 =
inactive).
</dd>
<dt>
<code>wt71</code>
</dt>
<dd>
Weight in kg at baseline (1971).
</dd>
<dt>
<code>censored</code>
</dt>
<dd>
Lost to follow-up indicator (1 = censored, 0 = not).
</dd>
</dl>

## Details

The dataset is used as the running example in Hernán & Robins (2025).
Key results to replicate:

<ul>
<li>

Observed mean difference E\[Y|A=1\] − E\[Y|A=0\] ≈ 2.5 kg (crude).

</li>
<li>

G-formula (standardised) estimate ≈ 3.5 kg (Chapter 13).

</li>
<li>

IPW estimate ≈ 3.4 kg (Chapter 12).

</li>
</ul>

The full NHEFS dataset is publicly available from the book’s website:
<a href="https://www.hsph.harvard.edu/miguel-hernan/causal-inference-book/">https://www.hsph.harvard.edu/miguel-hernan/causal-inference-book/</a>

## References

Hernán MA, Robins JM (2025). <em>Causal Inference: What If</em>. Chapman
& Hall/CRC.

## Examples

``` r
library("causatr")

data("nhefs", package = "causatr")
head(nhefs)
```

        seqn  qsmk    wt82_71   sex   age  race education smokeintensity smokeyrs
       <int> <int>      <num> <int> <int> <int>     <int>          <int>    <int>
    1:   233     0 -10.093960     0    42     1         7             30       29
    2:   235     0   2.604970     0    36     0         9             20       24
    3:   244     0   9.414486     1    56     1        11             20       26
    4:   245     0   4.990117     0    68     1         5              3       53
    5:   252     0   4.989251     0    40     0        11             20       19
    6:   257     0   4.419060     1    43     1         9             10       21
       exercise active  wt71 censored
          <int>  <int> <num>    <int>
    1:        2      0 79.04        0
    2:        0      0 58.63        0
    3:        2      0 56.81        0
    4:        2      1 59.42        0
    5:        1      1 87.09        0
    6:        1      1 99.00        0

``` r
table(nhefs$qsmk)
```


       0    1 
    1282  464 
