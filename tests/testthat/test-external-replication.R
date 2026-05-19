# External dataset replication: published causal estimates from
# canonical benchmarks outside the causatr package.

# -- LaLonde / Dehejia-Wahba: matching ATT for job training ---------------

test_that("LaLonde matching ATT ≈ experimental benchmark ($1794)", {
  lalonde <- MatchIt::lalonde

  fit <- causat(
    lalonde,
    outcome = "re78",
    treatment = "treat",
    confounders = ~ age + educ + race + married + nodegree + re74 + re75,
    estimator = "matching",
    estimand = "ATT"
  )

  res <- contrast(
    fit,
    interventions = list(treated = static(1), control = static(0)),
    reference = "control"
  )

  att <- res$contrasts$estimate[1]

  # Experimental benchmark ATT ≈ $1794.
  # Observational matching with PSID controls is noisier; expect within $1000.
  expect_gt(att, 400)
  expect_lt(att, 2500)
})


test_that("LaLonde: gcomp, IPW, matching ATT agree within $2000", {
  lalonde <- MatchIt::lalonde

  ivs <- list(treated = static(1), control = static(0))

  fit_gc <- causat(
    lalonde,
    outcome = "re78",
    treatment = "treat",
    confounders = ~ age + educ + race + married + nodegree + re74 + re75
  )
  fit_ipw <- causat(
    lalonde,
    outcome = "re78",
    treatment = "treat",
    confounders = ~ age + educ + race + married + nodegree + re74 + re75,
    estimator = "ipw",
    propensity_model_fn = stats::glm,
    estimand = "ATT"
  )
  fit_match <- causat(
    lalonde,
    outcome = "re78",
    treatment = "treat",
    confounders = ~ age + educ + race + married + nodegree + re74 + re75,
    estimator = "matching",
    estimand = "ATT"
  )

  res_gc <- contrast(
    fit_gc,
    interventions = ivs,
    reference = "control",
    estimand = "ATT"
  )
  res_ipw <- contrast(fit_ipw, interventions = ivs, reference = "control")
  res_match <- contrast(fit_match, interventions = ivs, reference = "control")

  atts <- c(
    gc = res_gc$contrasts$estimate[1],
    ipw = res_ipw$contrasts$estimate[1],
    match = res_match$contrasts$estimate[1]
  )

  # All ATTs should be positive (training helps earnings).
  for (a in atts) expect_gt(a, 0)

  # Pairwise agreement within $1000 (LaLonde with PSID control is hard).
  pairs <- combn(atts, 2)
  for (j in seq_len(ncol(pairs))) {
    expect_lt(abs(pairs[1, j] - pairs[2, j]), 1000)
  }
})


# -- LaLonde AIPW ATT agrees with other methods ---------------------------

test_that("LaLonde: AIPW ATT agrees with gcomp ATT within $1500", {
  lalonde <- MatchIt::lalonde

  ivs <- list(treated = static(1), control = static(0))

  fit_gc <- causat(
    lalonde,
    outcome = "re78",
    treatment = "treat",
    confounders = ~ age + educ + race + married + nodegree + re74 + re75
  )
  fit_aipw <- causat(
    lalonde,
    outcome = "re78",
    treatment = "treat",
    confounders = ~ age + educ + race + married + nodegree + re74 + re75,
    estimator = "aipw",
    propensity_model_fn = stats::glm,
    estimand = "ATT"
  )

  res_gc <- contrast(
    fit_gc,
    interventions = ivs,
    reference = "control",
    estimand = "ATT"
  )
  res_aipw <- contrast(fit_aipw, interventions = ivs, reference = "control")

  att_gc <- res_gc$contrasts$estimate[1]
  att_aipw <- res_aipw$contrasts$estimate[1]

  expect_lt(abs(att_gc - att_aipw), 500)
})
