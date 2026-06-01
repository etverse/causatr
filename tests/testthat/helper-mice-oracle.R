# Oracle wrappers for Rubin's-rules pooling. mice::pool.scalar() implements
# the scalar Rubin (1987) combination with Barnard-Rubin (1999) degrees of
# freedom and is the reference for causatr's pool_rubin() engine. These
# helpers fit causat() + contrast() on every completed dataset of a `mids`
# object and pool the resulting per-imputation estimate / variance with the
# external oracle, so a test can assert that causat_mice() reproduces it.

# Collect per-imputation (estimate, variance) for one contrast row across all
# imputations of a `mids` object.
#   fit_args:      named list forwarded to causat() (without `data`)
#   interventions: passed to contrast()
#   type:          contrast scale (default "difference")
#   row:           index into the contrasts table
mi_oracle_collect <- function(
  imp,
  fit_args,
  interventions,
  type = "difference",
  row = 1L
) {
  m <- imp$m
  Q <- numeric(m)
  U <- numeric(m)
  for (i in seq_len(m)) {
    d <- mice::complete(imp, i)
    fit <- do.call(causat, c(list(data = d), fit_args))
    res <- contrast(
      fit,
      interventions = interventions,
      type = type,
      ci_method = "sandwich"
    )
    Q[i] <- res$contrasts$estimate[row]
    U[i] <- res$contrasts$se[row]^2
  }
  list(Q = Q, U = U)
}

# Pool a collected (Q, U) with mice::pool.scalar(). `k` is the number of
# parameters consumed by the complete-data analysis; pool_table_rubin() uses
# dfcom = n - k, so the same k must be passed here for an exact df match. For
# a single contrast row the contrast table has one column, so k = 1.
mi_oracle_pool <- function(Q, U, n, k) {
  ps <- mice::pool.scalar(Q, U, n = n, k = k)
  list(estimate = ps$qbar, t = ps$t, se = sqrt(ps$t), df = ps$df, fmi = ps$fmi)
}
