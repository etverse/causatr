# causat() rejects missing outcome column

    Code
      causat(df, outcome = "Y", treatment = "A", confounders = ~L)
    Condition
      Error in `causat()`:
      ! Column `Y` (outcome) not found in `data`.

# causat() rejects missing treatment column

    Code
      causat(df, outcome = "Y", treatment = "A", confounders = ~L)
    Condition
      Error in `causat()`:
      ! Column `A` (treatment) not found in `data`.

# causat() rejects missing confounder column

    Code
      causat(df, outcome = "Y", treatment = "A", confounders = ~L)
    Condition
      Error in `causat()`:
      ! Confounder variable(s) not found in `data`: L

# causat() rejects id without time

    Code
      causat(df, outcome = "Y", treatment = "A", confounders = ~L, id = "id")
    Condition
      Error in `causat()`:
      ! Both `id` and `time` must be provided together for longitudinal data.

# causat() rejects time without id

    Code
      causat(df, outcome = "Y", treatment = "A", confounders = ~L, time = "t")
    Condition
      Error in `causat()`:
      ! Both `id` and `time` must be provided together for longitudinal data.

# causat() rejects NA / non-finite / negative / mis-sized weights

    Code
      causat(df, outcome = "Y", treatment = "A", confounders = ~L, weights = c(NA,
        rep(1, 49)))
    Condition
      Error in `causat()`:
      ! `weights` contains 1 missing value(s). Drop those rows or impute before calling `causat()`.

---

    Code
      causat(df, outcome = "Y", treatment = "A", confounders = ~L, weights = c(Inf,
        rep(1, 49)))
    Condition
      Error in `causat()`:
      ! `weights` contains non-finite value(s) (Inf / NaN).

---

    Code
      causat(df, outcome = "Y", treatment = "A", confounders = ~L, weights = c(-1,
        rep(1, 49)))
    Condition
      Error in `causat()`:
      ! `weights` must be non-negative.

---

    Code
      causat(df, outcome = "Y", treatment = "A", confounders = ~L, weights = rep(1,
        40))
    Condition
      Error in `causat()`:
      ! `weights` must have length equal to `nrow(data)` (50), got 40.

---

    Code
      causat(df, outcome = "Y", treatment = "A", confounders = ~L, weights = as.character(
        rep(1, 50)))
    Condition
      Error in `causat()`:
      ! `weights` must be numeric.

# causat() rejects ATT for continuous treatment

    Code
      causat(df, outcome = "Y", treatment = "A", confounders = ~L, estimand = "ATT")
    Condition
      Error in `causat()`:
      ! estimand = 'ATT' is only defined for binary point treatments coded as 0/1. Use estimand = 'ATE' or subset = quote(...) for subgroup effects (or recode the treatment as integer 0/1 if it already has two levels).

# causat() rejects ATT for multivariate treatment

    Code
      causat(df, outcome = "Y", treatment = c("A1", "A2"), confounders = ~L,
      estimand = "ATT")
    Condition
      Error in `causat()`:
      ! estimand = 'ATT' is only defined for binary point treatments coded as 0/1. Use estimand = 'ATE' or subset = quote(...) for subgroup effects (or recode the treatment as integer 0/1 if it already has two levels).

# causat() aborts when treatment has NAs and no censoring

    Code
      causat(df, outcome = "Y", treatment = "A", confounders = ~L)
    Condition
      Error in `causat()`:
      ! Treatment variable 'A' has 1 missing value.
      i Use `censoring = '...'` for inverse probability of censoring weights.
      i Or remove incomplete cases before calling `causat()`.

# causat() rejects missing confounders_tv column

    Code
      causat(df, outcome = "Y", treatment = "A", confounders = ~L, confounders_tv = ~
        CD4, id = "id", time = "time")
    Condition
      Error in `causat()`:
      ! Time-varying confounder variable(s) not found in `data`: CD4

# causat() rejects invalid history value

    Code
      causat(df, outcome = "Y", treatment = "A", confounders = ~L, id = "id", time = "time",
        history = 0)
    Condition
      Error in `causat()`:
      ! `history` must be a positive integer or `Inf`.

# causat_survival() aborts on non-NULL competing= until Phase 6

    Code
      causat_survival(long, outcome = "Y", treatment = "A", confounders = ~L, id = "id",
        time = "t", competing = "cmp")
    Condition
      Error in `causat_survival()`:
      ! Competing-risks survival analysis is not yet implemented (planned for Phase 6). The `competing` argument is reserved but currently has no effect; re-running with `competing = NULL` would silently fit a plain cause-deleted hazard model, which is biased in the presence of competing events. Track progress in PHASE_6_SURVIVAL.md.

