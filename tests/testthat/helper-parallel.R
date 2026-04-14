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
local({
  # Cap at 6 physical cores — beyond that the ICE bootstrap's
  # per-replicate workload hits diminishing returns from fork
  # overhead, and we want to leave at least one core free for
  # testthat's file-level worker.
  ncpus <- max(1L, min(6L, parallel::detectCores(logical = FALSE) - 1L))
  options(boot.ncpus = ncpus)
  if (.Platform$OS.type != "windows") {
    options(boot.parallel = "multicore")
  }
})
