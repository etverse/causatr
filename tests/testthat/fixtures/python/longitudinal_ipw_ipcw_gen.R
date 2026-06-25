# Generate the wide-format data fixture for the longitudinal IPW + IPCW
# delicatessen oracle. The data is causatr's own DGP-M5 (informative dropout):
# `simulate_longitudinal_mar_outcome(n = 2000, seed = 2025)`, a 2-period linear
# SCM with treatment-confounder feedback and censoring at the final period that
# depends on A_0 and L0 (so the censoring cross-term is large and load-bearing).
#
# The R test fits causatr on the same `simulate_longitudinal_mar_outcome()` call
# (deterministic given the seed); this script only reshapes the same data to the
# one-row-per-id wide layout the Python oracle reads.
#
# Run from the repo root:
#   Rscript tests/testthat/fixtures/python/longitudinal_ipw_ipcw_gen.R

library(data.table)
source(file.path("tests", "testthat", "helper-dgp.R"))

d <- as.data.table(simulate_longitudinal_mar_outcome(n = 2000, seed = 2025))

w0 <- d[time == 0, .(id, L0, A0 = A, Ltv0 = L)]
w1 <- d[time == 1, .(id, A1 = A, Ltv1 = L, C1 = C, Y = Y)]
wide <- merge(w0, w1, by = "id")
setorder(wide, id)

out <- file.path(
  "tests",
  "testthat",
  "fixtures",
  "python",
  "longitudinal_ipw_ipcw_data.csv"
)
fwrite(wide, out)
cat("wrote", out, "with", nrow(wide), "rows;", "censored:", sum(wide$C1), "\n")
