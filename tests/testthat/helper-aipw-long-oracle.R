# Oracle helpers for the longitudinal-AIPW stacked-EE sandwich.
#
# Two independent cross-checks of `variance_if_aipw_longitudinal()`:
#   * `jackknife_aipw_long_var()` -- delete-one-id jackknife variance of
#     the ATE (a resampling estimator, independent of the analytic IF).
#   * `build_aipw_long_pp()` -- wide -> person-period reshaper for the
#     delicatessen cross-language fixtures.

# Delete-one-id jackknife variance for the longitudinal-AIPW contrast
# `arm - reference`. Refits the full ICE-AIPW pipeline with each id
# removed and recomputes the point estimate through the internal contrast
# machinery (the cheap point path -- no per-replicate variance). On an
# unbalanced panel the jackknife recovers the same larger truth the
# bootstrap does, so the analytic sandwich must match it.
#
# `fit_fn(data)` returns a `causat()` fit on the supplied data; `ivs` is
# the named intervention list; `arm`/`reference` name the contrast arms.
# Full-population ATE only (no `by` / subset).
jackknife_aipw_long_var <- function(data, fit_fn, ivs, reference, arm) {
  data <- data.table::as.data.table(data)
  ids <- unique(data[["id"]])
  m <- length(ids)
  est <- vapply(ids, function(i) {
    di <- data[data[["id"]] != i]
    fi <- fit_fn(di)
    tp <- fi$details$time_points
    rf <- fi$data[[fi$time]] == tp[1]
    tw <- rep(TRUE, sum(rf))
    res <- compute_aipw_contrast_longitudinal(fi, ivs, tw, trim = 1)
    unname(res$mu_hat[[arm]] - res$mu_hat[[reference]])
  }, numeric(1))
  # Jackknife variance of the contrast: (m-1)/m * sum (theta_{-i} - mean)^2.
  (m - 1) / m * sum((est - mean(est))^2)
}

# Reshape a wide longitudinal-AIPW fixture (one row per id with columns
# L0, A0, L1, A1, Y; NA at t=1 marks monotone dropout) into the
# person-period frame causatr expects. Time-0 rows carry L = L0; present
# time-1 rows carry L = L1 and the end-of-follow-up outcome Y.
build_aipw_long_pp <- function(wide) {
  present <- !is.na(wide$A1)
  pp0 <- data.table::data.table(
    id = wide$id, time = 0L, L = wide$L0, A = wide$A0, Y = NA_real_
  )
  pp1 <- data.table::data.table(
    id = wide$id[present], time = 1L, L = wide$L1[present],
    A = wide$A1[present], Y = wide$Y[present]
  )
  data.table::setorder(rbind(pp0, pp1), id, time)[]
}
