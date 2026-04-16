# Phase 13 — Beta Regression, Negative Binomial, and Extended Outcome Types

> **Status: PENDING (design doc)**
>
> **Depends on:** Phase 2 (point gcomp) is done. Phases 3–5 provide the
> IPW / matching / ICE engines that also benefit.

## Motivation

causatr currently supports all `stats::` GLM families (gaussian,
binomial, quasibinomial, poisson, quasipoisson, Gamma,
inverse.gaussian) and anything pluggable via `model_fn`. Two commonly
requested outcome types have partial or untested support:

1. **Negative binomial** — already works via `model_fn = MASS::glm.nb`
   but has **zero explicit test coverage**. The sandwich variance engine
   routes through the analytic GLM path (the fitted object carries
   `family$mu.eta` and `family$variance`) but this has never been
   verified. Users don't know it works because it's undocumented.

2. **Beta regression** — for bounded-continuous outcomes (proportions,
   rates, indices) that live in (0, 1). Requires `betareg::betareg()`
   or passing the `betareg::betar()` family to `stats::glm()`. The
   current `resolve_family()` rejects `"beta"` as a string (not in
   `stats::` namespace) but **accepts family objects directly**, so
   `family = betareg::betar()` already works at the model-fitting
   level. However: (a) the variance engine path is untested, (b) no
   documentation mentions it, and (c) `resolve_family()` could be
   extended to recognise the string `"beta"` for convenience.

Both are low-effort, high-value additions that expand the outcome type
coverage without architectural changes.

## References

- Cameron AC, Trivedi PK (2013). *Regression Analysis of Count Data*.
  Cambridge University Press. (Negative binomial theory.)
- Cribari-Neto F, Zeileis A (2010). Beta regression in R. *Journal of
  Statistical Software* 34(2):1–24. (betareg package.)
- Papke LE, Wooldridge JM (1996). Econometric methods for fractional
  response variables with an application to 401(k) plan participation
  rates. *Journal of Applied Econometrics* 11:619–632.
  (Fractional logistic as an alternative to beta regression.)

---

## Feature 1: Negative Binomial (`MASS::glm.nb`)

### Current state

- `MASS::glm.nb` is already mentioned in the `@param model_fn`
  documentation in `R/causat.R` and `R/gcomp.R`.
- `MASS` is in `Suggests:`.
- The fitted model from `glm.nb` inherits from `"negbin"` and `"glm"`,
  exposes `$family` with `mu.eta` and `variance`, and has working
  `residuals()`, `model.matrix()`, `predict()`, and
  `sandwich::estfun()` methods.
- The variance engine's analytic GLM path should handle it — but this
  has never been tested.

### Plan

#### Chunk NB-1: Truth-based tests

1. **DGP:** Binary treatment, negative-binomial outcome.
   ```
   L ~ N(2, 1)
   A ~ Bernoulli(expit(-1 + 0.5 * L))
   Y ~ NegBin(mu = exp(0.5 + 0.3 * A + 0.4 * L), size = 2)
   ```
   Analytical truth for `E[Y^1] - E[Y^0]`:
   ```
   E[Y^a] = E_L[ exp(0.5 + 0.3*a + 0.4*L) ]
          = exp(0.5 + 0.3*a) * E_L[exp(0.4*L)]
          = exp(0.5 + 0.3*a) * exp(0.4*2 + 0.5*0.4^2)   (MGF of N(2,1))
   ATE = E[Y^1] - E[Y^0] = exp(0.5 + 0.4*2 + 0.08) * (exp(0.3) - 1)
   ```

2. **Tests (in `test-simulation.R`):**
   - `gcomp × binary trt × negbin outcome × difference × sandwich`
     — point estimate within 0.1 of analytical truth on n = 5000.
   - `gcomp × binary trt × negbin outcome × ratio × sandwich`
     — ratio estimate ≈ `exp(0.3)`.
   - `gcomp × binary trt × negbin outcome × difference × bootstrap`
     — bootstrap SE within 20% of sandwich SE.

3. **IPW test:**
   - `ipw × binary trt × negbin outcome × static × ATE × sandwich`
     — point estimate agrees with gcomp (same DGP, different path).

4. **Matching test:**
   - `matching × binary trt × negbin outcome × static × ATT × sandwich`
     — smoke test (runs, finite SE).

#### Chunk NB-2: Documentation

1. Add `MASS::glm.nb` examples to `causat()` roxygen `@examples`.
2. Mention negative binomial in the DESCRIPTION outcome list.
3. Update `FEATURE_COVERAGE_MATRIX.md` with the new rows.

---

## Feature 2: Beta Regression (`betareg`)

### Current state

- `resolve_family()` (`R/utils.R`) resolves family strings by looking
  up `stats::` namespace only. `"beta"` triggers `"Unknown family"`.
- However, `resolve_family()` already has a third branch that returns
  a family **object** as-is. So `family = betareg::betar()` works if
  the user passes it explicitly.
- The fitted model from `glm(..., family = betareg::betar())` (note:
  not `betareg::betareg()`, which has a different interface) exposes
  `$family` with `mu.eta` and `variance`.
- Alternatively, the user can use `betareg::betareg()` as a
  `model_fn`, but that function has a different signature
  (`betareg(formula, data, ...)` — no `family` argument) and returns
  a `"betareg"` class, not `"glm"`. The variance engine would fall
  back to `variance_if_numeric()` (Tier 1 if `sandwich::estfun`
  exists for betareg, Tier 2 otherwise).

**Recommended approach:** support beta regression via the
`stats::glm()` + `betareg::betar()` family path, which keeps the
model as a standard `glm` object and the analytic sandwich path works
directly. The dedicated `betareg::betareg()` fitter is a secondary
path that works via the numeric fallback.

### Plan

#### Chunk BETA-1: `resolve_family()` extension

Extend `resolve_family()` to recognise non-`stats::` family strings
by name, with `betareg` in `Suggests:`:

```r
resolve_family <- function(family) {
  if (is.character(family)) {
    # Try stats:: namespace first (gaussian, binomial, poisson, ...).
    fam_fn <- tryCatch(
      get(family, mode = "function", envir = asNamespace("stats")),
      error = function(e) NULL
    )
    if (!is.null(fam_fn)) return(fam_fn())

    # Extended families from Suggests packages.
    extended <- list(
      beta = function() {
        rlang::check_installed("betareg",
          reason = "for beta regression family"
        )
        betareg::betar()
      }
    )
    if (family %in% names(extended)) return(extended[[family]]())

    rlang::abort(paste0(
      "Unknown family: '", family, "'. ",
      "Supported: gaussian, binomial, poisson, quasibinomial, ",
      "quasipoisson, Gamma, inverse.gaussian, beta. ",
      "Or pass a family object directly."
    ))
  }
  if (is.function(family)) return(family())
  family
}
```

This is minimal — one new branch, `betareg` stays in `Suggests:`.

#### Chunk BETA-2: Truth-based tests

1. **DGP:** Binary treatment, beta-distributed outcome in (0, 1).
   ```
   L ~ N(0, 1)
   A ~ Bernoulli(expit(-0.5 + 0.3 * L))
   mu = expit(0.2 + 0.5 * A + 0.3 * L)
   Y ~ Beta(mu * phi, (1 - mu) * phi),  phi = 10
   ```
   Analytical truth:
   ```
   E[Y^a] = E_L[ expit(0.2 + 0.5*a + 0.3*L) ]
   ```
   Computed numerically via `integrate()` against the N(0,1) density.

2. **Tests (in `test-simulation.R`):**
   - `gcomp × binary trt × beta outcome × difference × sandwich`
     — point estimate within tolerance of numerical truth.
   - `gcomp × binary trt × beta outcome × ratio × sandwich`
     — ratio estimate.
   - `gcomp × binary trt × beta outcome × OR × sandwich`
     — OR estimate (valid since outcome is in (0,1)).

3. **IPW test:**
   - `ipw × binary trt × beta outcome × static × ATE × sandwich`
     — point estimate agrees with gcomp.

4. **betareg::betareg() as model_fn path:**
   - Smoke test only (different model class, numeric variance fallback).

#### Chunk BETA-3: Documentation

1. Add beta regression example to `causat()` roxygen.
2. Add `betareg` to `Suggests:` in DESCRIPTION.
3. Update `FEATURE_COVERAGE_MATRIX.md`.

---

## Contrast type compatibility

| Outcome type | difference | ratio | OR |
|---|---|---|---|
| Negative binomial | yes | yes (natural, log link) | no (not a probability) |
| Beta | yes | yes | yes (outcome in (0,1)) |

No changes needed to the contrast machinery — the existing `ratio`
and `or` validation gates (reference mean > 0, mean in (0,1) for OR)
handle both types correctly:
- Negative binomial means are positive → ratio works, OR gate rejects.
- Beta means are in (0, 1) → ratio and OR both work.

## Variance engine compatibility

| Model path | `family$mu.eta` | `family$variance` | Sandwich path |
|---|---|---|---|
| `glm(..., family = MASS::negative.binomial(theta))` | yes | yes | Analytic GLM |
| `MASS::glm.nb(...)` (model_fn) | yes | yes | Analytic GLM |
| `glm(..., family = betareg::betar())` | yes | yes | Analytic GLM |
| `betareg::betareg(...)` (model_fn) | no (`betareg` class) | no | Numeric Tier 1 (via `sandwich::estfun.betareg`) |

No variance engine changes needed. The analytic path handles both
MASS and betareg families via the existing `has_family` check.

---

## Chunks (combined)

| Chunk | Scope | Deliverables |
|---|---|---|
| NB-1 | Negative binomial truth tests | 4 tests in `test-simulation.R` |
| NB-2 | Negative binomial docs | Examples, DESCRIPTION, matrix |
| BETA-1 | `resolve_family()` extension | One function edit + unit test |
| BETA-2 | Beta regression truth tests | 4 tests in `test-simulation.R` |
| BETA-3 | Beta regression docs | Examples, DESCRIPTION, matrix |

**Total estimated edits:** ~150 lines of test code, ~20 lines of
`resolve_family()`, ~30 lines of roxygen examples, ~10 lines in
DESCRIPTION and FEATURE_COVERAGE_MATRIX.md.

---

## Invariants

- `MASS::glm.nb` MUST produce a model with `$family$mu.eta` and
  `$family$variance` so the analytic sandwich path fires (not the
  numeric fallback). Verify in tests.
- `betareg::betar()` MUST return a family object accepted by
  `stats::glm()`. Verify via `is.list(betareg::betar())` and presence
  of `$mu.eta`, `$variance`, `$linkinv`.
- `resolve_family("beta")` MUST abort with `rlang::check_installed()`
  if `betareg` is not installed, not with a cryptic `namespace` error.
- Negative binomial ratio contrasts MUST work (means are positive).
  OR contrasts MUST be rejected (means are not probabilities).
- Beta regression OR contrasts MUST work (means are in (0,1)).
