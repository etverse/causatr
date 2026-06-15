# Test-time parallelism and tier configuration.
#
# CAUSATR_TEST_TIER controls which tests run:
#   "fast"  — skips the heavy bootstrap / external-oracle / large-n blocks
#             (~4 CPU-min); used for pull-request CI to give fast feedback.
#   "full"  — runs everything (default; local, push-to-default-branch, nightly).
#
# Usage:
#   CAUSATR_TEST_TIER=fast Rscript -e 'devtools::test()'
#   devtools::test(filter = "gcomp")  # targeted filter is always fast
#
# Heavy blocks across the suite call skip_if_fast() as their first statement.
# In full mode skip_if_fast() is a no-op, so local and main runs are unaffected.

#' Skip a test block when CAUSATR_TEST_TIER is "fast".
#' @noRd
skip_if_fast <- function(msg = "slow test (CAUSATR_TEST_TIER=fast)") {
  tier <- Sys.getenv("CAUSATR_TEST_TIER", "full")
  if (identical(tolower(tier), "fast")) {
    testthat::skip(msg)
  }
}

# Test-time parallelism configuration.
#
# 1. testthat 3 reads `Config/testthat/parallel: true` from DESCRIPTION
#    and runs each test-*.R file in a separate worker process. That
#    is configured in DESCRIPTION; nothing to do here for the file-
#    level layer.
#
# 2. Within each test process, `causat()`'s bootstrap (`ci_method =
#    "bootstrap"`) honours the `boot.ncpus` and `boot.parallel`
#    options via `boot::boot()`. Setting them once here means every
#    bootstrap call in every test gets parallel execution without
#    having to thread `parallel = "multicore", ncpus = N` through
#    each individual `contrast()` call.
#
#    `parallel::detectCores(logical = FALSE)` would give physical
#    cores; we cap at 4 to leave headroom on CI runners and avoid
#    thrashing nested parallel layers (testthat-files × boot-cpus).
if (.Platform$OS.type == "windows") {
  data.table::setDTthreads(1L)
}

# Suppress the model_fn / propensity_model_fn default warnings globally
# so that ~900 causat() calls in tests don't need individual suppression.
# Tests that verify warning emission temporarily unset this option.
options(causatr.suppress_default_warnings = TRUE)

local({
  # Cap at 6 physical cores — beyond that the ICE bootstrap's
  # per-replicate workload hits diminishing returns from fork
  # overhead, and we want to leave at least one core free for
  # testthat's file-level worker.
  #
  # Under `R CMD check` the env var `_R_CHECK_LIMIT_CORES_` is set,
  # which makes `parallel:::.check_ncores()` abort as soon as more
  # than 2 cores are requested ("N simultaneous processes spawned").
  # Cap at 2 in that environment so tests don't die on CI runners
  # that report 3+ physical cores.
  limit_cores <- tolower(Sys.getenv("_R_CHECK_LIMIT_CORES_", ""))
  max_cores <- if (limit_cores %in% c("true", "warn")) 2L else 6L
  ncpus <- max(1L, min(max_cores, parallel::detectCores(logical = FALSE) - 1L))
  options(boot.ncpus = ncpus)
  if (.Platform$OS.type != "windows") {
    options(boot.parallel = "multicore")
  }
})
