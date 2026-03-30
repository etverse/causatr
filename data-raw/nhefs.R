## Download and prepare the NHEFS dataset
##
## Source: Hernán MA, Robins JM (2025). Causal Inference: What If.
##   Raw data mirrored from:
##   https://github.com/sdwfrost/hernan-robins-julia-code
##
## Run this script via `source("data-raw/nhefs.R")` or `make data`
## to regenerate data/nhefs.rda.

nhefs_raw <- utils::read.csv("data-raw/nhefs.csv")

nhefs <- data.table::as.data.table(nhefs_raw)

nhefs <- nhefs[, .(
  seqn,
  qsmk,
  wt82_71,
  sex,
  age,
  race,
  education = school,
  smokeintensity,
  smokeyrs,
  exercise,
  active,
  wt71,
  censored = as.integer(is.na(wt82_71))
)]

usethis::use_data(nhefs, overwrite = TRUE)
