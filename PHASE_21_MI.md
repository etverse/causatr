# Phase 21 — Multiple Imputation via `causat_mice()`

## Motivation

causatr currently handles missing data via complete-case analysis
(`get_fit_rows()` drops rows with missing outcome or censored status)
and, after Phase 14, via built-in IPCW for MAR outcome censoring.
Neither approach handles **missing covariates (L)** or **missing
treatment (A)** under MAR — the two scenarios where multiple imputation
(MI) via `mice::mice()` is the standard solution.

Phase 21 implements `causat_mice()`, a wrapper that applies Rubin's
rules to pool causal estimates across multiply-imputed datasets. The
function already exists as a stub in `R/causat_mice.R` with the full
roxygen skeleton; this phase fills in the implementation.

### What MI solves that IPCW and complete-case do not

| Variable | Mechanism | Complete-case | IPCW (Phase 14) | MI (Phase 21) |
|---|---|---|---|---|
| Y (outcome) | MCAR | Unbiased | N/A | Not needed |
| Y (outcome) | MAR on (A,L) | G-comp OK if model correct; IPW biased | Fixes IPW; DR for g-comp | Not the right tool |
| A (treatment) | MCAR | User drops rows; unbiased | N/A | Not needed |
| A (treatment) | MAR on L | Biased — informative missingness | N/A | Impute A from P(A\|L) |
| L (covariates) | MCAR | G-comp: `na.omit`; IPW: aborts | N/A | Not needed (but helps efficiency) |
| L (covariates) | MAR | Biased — model misspecified on subset | N/A | Impute L from P(L\|A,Y,L') |

**IPCW and MI are orthogonal.** IPCW reweights for missing Y; MI
imputes missing A or L. A dataset can have both problems
simultaneously (missing Y handled by `ipcw = TRUE`, missing L handled
by `causat_mice()`). Rubin's rules pool across imputations; each
imputed dataset runs `causat(..., ipcw = TRUE) |> contrast()` if
outcome censoring is also present.

### Why MI and not single imputation

Single imputation underestimates variance because it treats the
imputed values as known. MI with Rubin's rules correctly decomposes
total variance into within-imputation (sampling variability) and
between-imputation (imputation uncertainty) components. For causal
inference, underestimating variance means anti-conservative CIs — a
worse failure mode than conservative CIs.

## Current state

- `R/causat_mice.R` contains a complete roxygen skeleton for
  `causat_mice()` with parameter docs, Rubin's rules derivation, and
  examples. The body always aborts with "not implemented".
- `mice` is in `Suggests` in `DESCRIPTION`.
- `check_treatment_nas()` in `R/checks.R` (line 391) aborts when A
  has NAs and no censoring column is provided, directing users to
  either provide `censoring =` or remove incomplete cases. Phase 21
  adds MI as a third option.

## User API

```r
library(mice)

# Step 1: impute missing covariates (and/or treatment)
imp <- mice(data, m = 20, method = "pmm")

# Step 2: fit + contrast across imputations, pool via Rubin's rules
result <- causat_mice(
  imp,
  outcome = "Y",
  treatment = "A",
  confounders = ~ L1 + L2,
  interventions = list(treat_all = static(1), treat_none = static(0)),
  estimator = "gcomp",
  type = "difference",
  pool_method = "rubin"   # default; "boot_mi" for valid inference under uncongeniality
)

# result is a causatr_result with pooled estimates
summary(result)
tidy(result)
```

### With IPCW (missing Y + missing L)

```r
imp <- mice(data[, c("A", "L1", "L2")], m = 20)  # impute L only

result <- causat_mice(
  imp,
  outcome = "Y",
  treatment = "A",
  confounders = ~ L1 + L2,
  interventions = list(a1 = static(1), a0 = static(0)),
  censoring = "C",
  ipcw = TRUE,
  type = "difference"
)
```

Each imputed dataset gets its own censoring model fit, IPCW weights,
and causal estimate. Rubin's rules pool across the m estimates.

### Longitudinal

```r
result <- causat_mice(
  imp,
  outcome = "Y",
  treatment = "A",
  confounders = ~ L1 + L2,
  confounders_tv = ~ L_tv,
  interventions = list(always = static(1), never = static(0)),
  id = "id",
  time = "time",
  type = "difference"
)
```

## Rubin's rules

Let $\hat{Q}_i$ and $U_i$ denote the point estimate and variance from
imputation $i$, for $i = 1, \ldots, m$.

**Pooled estimate:**

$$\bar{Q} = \frac{1}{m} \sum_{i=1}^m \hat{Q}_i$$

**Within-imputation variance:**

$$\bar{U} = \frac{1}{m} \sum_{i=1}^m U_i$$

**Between-imputation variance:**

$$B = \frac{1}{m-1} \sum_{i=1}^m (\hat{Q}_i - \bar{Q})^2$$

**Total variance:**

$$T = \bar{U} + B + B/m = \bar{U} + (1 + 1/m) B$$

**Degrees of freedom** (Barnard-Rubin):

$$\nu = (m - 1) \left(1 + \frac{\bar{U}}{(1 + 1/m) B}\right)^2$$

**Confidence interval:**

$$\bar{Q} \pm t_{\nu, 1-\alpha/2} \sqrt{T}$$

### Multivariate pooling

For contrast tables with multiple rows (e.g. `by =` stratification,
multiple intervention pairs), Rubin's rules are applied **per
contrast row independently**. Joint inference across contrasts (e.g.
omnibus test) is out of scope.

### Variance / CI strategy

- `pool_method = "rubin"` (default): Standard Rubin's rules. Each
  imputation returns a sandwich variance $U_i$ from `contrast(...,
  ci_method = "sandwich")`. Rubin's rules add the between-imputation
  component. Fast ($m$ fits). Under uncongeniality (the norm for causal
  inference — see below), CIs are typically conservative (wider than
  nominal). This is acceptable for most applications.

- `pool_method = "boot_mi"`: Bootstrap-then-impute (Boot MI). Draws
  $B$ bootstrap samples from the **incomplete** (pre-imputation) data,
  runs `mice()` + `causat()` + `contrast()` on each bootstrap sample,
  and takes the bootstrap variance across samples. Valid under
  uncongeniality because the bootstrap captures the full sampling
  distribution without relying on Rubin's variance decomposition.
  Slower ($B \times M$ fits) but methodologically correct when the
  imputation and analysis models are not congenial. Uses von Hippel's
  efficient variant: $M = 2$ imputations per bootstrap sample, with a
  one-way random effects decomposition to separate bootstrap variance
  from residual imputation noise (von Hippel 2020). Recommended $B =
  200$, giving $200 \times 2 = 400$ total fits.

  **Why not MI-boot-Rubin?** Running a bootstrap *within* each
  imputation and then pooling bootstrap variances via Rubin's rules
  ("MI boot Rubin" in Schomaker & Heumann 2018) does not fix
  uncongeniality — the Rubin pooling step still introduces the same
  bias. Boot MI reverses the order so that the bootstrap, not Rubin's
  formula, is responsible for variance estimation.

## Estimator support

All estimators supported by `causat()` work with `causat_mice()`:

| Estimator | Point | Longitudinal | Notes |
|---|---|---|---|
| gcomp | Yes | Yes (ICE) | Most natural MI target |
| IPW | Yes | Yes | Propensity model refit per imputation |
| matching | Yes | — | MatchIt rerun per imputation |

## Implementation chunks

| Chunk | Scope | Depends on |
|---|---|---|
| 21a | Core `causat_mice()`: loop over imputations, call `causat()` + `contrast()`, collect estimates and variances. Supports `pool_method = "rubin"` (default) | — |
| 21b | `pool_rubin()`: Rubin's rules engine (pooled estimate, total variance, Barnard-Rubin df, CIs). Return `causatr_result` with pooled contrasts table | 21a |
| 21c | `pool_boot_mi()`: Boot MI engine — bootstrap incomplete data $B$ times, impute each with $M = 2$, fit + contrast, apply von Hippel's random-effects decomposition for variance. Return `causatr_result` | 21a |
| 21d | S3 methods: `print`, `summary`, `tidy`, `coef`, `confint` for the pooled result. Reuse existing `causatr_result` methods where possible | 21b, 21c |
| 21e | Edge cases and validation: m=1 (degenerate), convergence failures in individual imputations, mixed missingness patterns across imputations, `by =` stratification pooling | 21b |
| 21f | IPCW + MI combination: `censoring =` + `ipcw = TRUE` passed through to per-imputation `causat()` calls. Verify pooled result under joint Y-censoring + L-missingness | Phase 14, 21b |
| 21g | Longitudinal MI: person-period imputation patterns, lag structure preservation across imputations | 21b |
| 21h | Truth-based simulation tests: DGPs with known ATE under MAR covariate missingness. Coverage tests for both `pool_method = "rubin"` and `pool_method = "boot_mi"`. Include uncongenial DGP (DGP-MI5) | 21b, 21c |
| 21i | Vignette: MI workflow guide with examples, congeniality discussion, guidance on when to use Boot MI | 21b, 21c |

## Design decisions

### Return type

`causat_mice()` returns a `causatr_result` (same class as
`contrast()`). The `contrasts` table has pooled estimates, pooled SEs,
and CIs. For `pool_method = "rubin"`, CIs use Barnard-Rubin df; for
`pool_method = "boot_mi"`, CIs use normal or percentile bootstrap.
An additional attribute or slot stores the per-imputation (or
per-bootstrap) estimates for diagnostics.

### Per-imputation storage

Store `list(estimates = numeric(m), variances = numeric(m))` per
contrast row in an attribute `"mi_details"` on the result. This lets
`diagnose()` or a future `plot()` method show between-imputation
variability.

### Parallelism

Each imputation (Rubin) or bootstrap sample (Boot MI) is independent.
Support `parallel =` argument forwarded to
`future.apply::future_lapply()` (same backend as the existing
bootstrap engine). Default: sequential. For Boot MI, the outer loop
is over bootstrap samples; imputation within each sample is
sequential ($M = 2$ is fast).

### Congeniality and causal inference

An imputation model is **congenial** with an analysis procedure when
both derive from (or are consistent with) the same joint model.
Rubin's rules are asymptotically unbiased under congeniality; under
uncongeniality the between-imputation variance $B$ typically
overestimates the true imputation uncertainty, producing conservative
CIs (Bartlett & Hughes 2020).

**Causal inference is uncongenial by default.** The causal estimand
$E[Y(a)]$ is a functional of the outcome/treatment model under
intervention — not a parameter of the mice imputation model:

- **G-computation:** estimand is $E[E[Y \mid A = a, L]]$, a
  marginalization of the outcome model. The imputation model for $L$
  (e.g. PMM from $P(L \mid A, Y, L')$) does not target this quantity.
- **IPW:** estimand uses density-ratio weights from a treatment model
  that is entirely absent from the imputation step.
- **Matching:** uses a propensity-score distance metric unrelated to
  the imputation model.

**Practical guidance:**

1. The mice predictor matrix should include all analysis variables
   (treatment, confounders, outcome, censoring indicator, effect
   modifiers). This does not guarantee congeniality but reduces bias
   in the imputed values. causatr should warn if variables in the
   `causat()` formula are absent from the `mids` predictor matrix.
2. For routine use, `pool_method = "rubin"` (default) is adequate:
   conservative CIs are safe and fast.
3. When efficiency matters (tight power budget, regulatory context),
   use `pool_method = "boot_mi"` for valid non-conservative inference
   regardless of congeniality.

### What `causat_mice()` does NOT do

- **Imputation itself.** Users call `mice::mice()` upstream. causatr
  is the analysis step, not the imputation step.
- **MNAR handling.** Sensitivity analysis for MNAR is out of scope.
- **Outcome imputation.** Missing Y is handled by IPCW (Phase 14) or
  complete-case analysis, not by MI. Imputing Y itself in a causal
  analysis is methodologically problematic (the imputation model would
  need to be compatible with the causal model, which is circular).
  However, Y **should** be included as a **predictor** in the
  imputation model when imputing L or A — excluding Y biases the
  imputed values and worsens uncongeniality (Moons et al. 2006; van
  Buuren 2018, ch. 6). The vignette should clearly distinguish
  "include Y as a predictor" (recommended) from "impute Y" (wrong
  tool — use IPCW).
- **Joint inference.** Omnibus tests across multiple contrasts are not
  pooled; each contrast row is pooled independently.

## DGP sketches for truth-based tests

### DGP-MI1: MAR covariate missingness (point treatment)

```
L ~ N(0, 1)
A ~ Bernoulli(expit(0.5 * L))
Y = 2 + 3*A + 1.5*L + N(0, 1)        # True ATE = 3
R_L ~ Bernoulli(expit(-1 + 0.8*A))    # L observed iff R_L = 1
L_obs = ifelse(R_L == 1, L, NA)
```

MAR: missingness in L depends on A (observed). Complete-case analysis
is biased because conditioning on R_L=1 induces collider bias. MI
recovers the correct ATE.

### DGP-MI2: MAR treatment missingness (point treatment)

```
L ~ N(0, 1)
A ~ Bernoulli(expit(0.5 * L))
Y = 2 + 3*A + 1.5*L + N(0, 1)        # True ATE = 3
R_A ~ Bernoulli(expit(-1 + 0.6*L))    # A observed iff R_A = 1
A_obs = ifelse(R_A == 1, A, NA)
```

MAR: missingness in A depends on L (observed). causat currently aborts
on NA treatment. MI imputes A from P(A|L).

### DGP-MI3: MAR covariate + outcome censoring (joint)

```
L ~ N(0, 1)
A ~ Bernoulli(expit(0.5 * L))
Y = 2 + 3*A + 1.5*L + N(0, 1)        # True ATE = 3
R_L ~ Bernoulli(expit(-1 + 0.8*A))    # L MAR
C ~ Bernoulli(expit(-2 + 0.6*A + 0.4*L))  # Y censored (MAR)
```

Requires MI for L + IPCW for Y. Tests the joint workflow.

### DGP-MI4: Longitudinal MAR covariate

```
L_0 ~ N(0, 1)
A_0 ~ Bernoulli(expit(0.3 * L_0))
L_1 = 0.5*L_0 + 0.3*A_0 + N(0, 1)
A_1 ~ Bernoulli(expit(0.3 * L_1))
Y = 10 + 2*(A_0 + A_1) + L_0 + L_1 + N(0, 1)   # True ATE = 5
R_{L1} ~ Bernoulli(expit(-1 + 0.5*A_0))          # L_1 MAR on A_0
```

Longitudinal MI with treatment-confounder feedback.

### DGP-MI5: Uncongenial imputation (point treatment)

```
L ~ N(0, 1)
A ~ Bernoulli(expit(0.5 * L))
Y = 2 + 3*A + 1.5*L + 0.8*A*L + N(0, 1)   # True ATE = 3 (marginal)
R_L ~ Bernoulli(expit(-1 + 0.8*A))          # L MAR on A
```

Imputation model: `mice(method = "pmm")` with Y **excluded** from
the predictor matrix (deliberately uncongenial — violates the
recommendation to include Y). This induces uncongeniality between the
imputation model $P(L \mid A)$ and the analysis model
$E[Y \mid A, L]$ which requires $P(L \mid A, Y)$.

**Expected behaviour:**
- `pool_method = "rubin"`: overcoverage (conservative CIs).
- `pool_method = "boot_mi"`: nominal coverage.

This DGP verifies that Boot MI corrects the conservatism of Rubin's
rules under uncongeniality. A companion run with Y included in the
predictor matrix (congenial) should show both methods at nominal
coverage, confirming the uncongeniality is the cause.

## References

- Rubin DB (1987). *Multiple Imputation for Nonresponse in Surveys*.
  Wiley.
- van Buuren S, Groothuis-Oudshoorn K (2011). mice: Multivariate
  Imputation by Chained Equations in R. *Journal of Statistical
  Software* 45(3):1-67.
- Barnard J, Rubin DB (1999). Small-sample degrees of freedom with
  multiple imputation. *Biometrika* 86:948-955.
- Moons KGM, Donders RART, Stijnen T, Harrell FE (2006). Using the
  outcome in a study to impute missing predictor values was preferred.
  *Journal of Clinical Epidemiology* 59:1092-1101.
- Seaman SR, White IR (2014). Inverse probability weighting with
  missing predictors of treatment assignment or missingness.
  *Communications in Statistics — Theory and Methods* 43:3499-3515.
- Bartlett JW, Hughes RA (2020). Bootstrap inference for multiple
  imputation under uncongeniality and misspecification. *Statistical
  Methods in Medical Research* 29(12):3533-3546. DOI:
  10.1177/0962280220932189.
- Schomaker M, Heumann C (2018). Bootstrap inference when using
  multiple imputation. *Statistics in Medicine* 37(14):2252-2266.
- von Hippel PT (2020). How many imputations do you need? A two-stage
  calculation using a quadratic rule. *Sociological Methods & Research*
  49(3):699-718.
- van Buuren S (2018). *Flexible Imputation of Missing Data*. 2nd ed.
  CRC Press.
- Leyrat C, Seaman SR, White IR, et al. (2019). Propensity score
  analysis with partially observed covariates: How should multiple
  imputation be used? *Statistical Methods in Medical Research*
  28:3-19.
