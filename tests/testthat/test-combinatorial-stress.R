# Combinatorial stress test: systematic grid of estimator × treatment ×
# outcome family × intervention × variance method. Every valid cell must
# produce finite output; cells with known truth check the point estimate.
#
# Uses existing DGPs from helper-dgp.R to avoid maintenance burden.

# -- DGP generators keyed by treatment type --------------------------------

dgp_generators <- list(
  binary = function() simulate_binary_continuous(n = 2000, seed = 42),
  binary_binom = function() simulate_binary_binary(n = 3000, seed = 42),
  continuous = function() simulate_continuous_continuous(n = 2000, seed = 42),
  categorical = function() simulate_categorical_continuous(n = 5000, seed = 42)
)

# Known truth values for point estimate checks.
truth_table <- list(
  binary = list(ate = 3, trt_col = "A", conf = ~L),
  binary_binom = list(ate = 0.333, trt_col = "A", conf = ~L),
  continuous = list(ate = -2, trt_col = "A", conf = ~L),
  categorical = list(ate = 3, trt_col = "A", conf = ~L)
)

# -- Valid combination grid ------------------------------------------------

combos <- data.frame(
  estimator = character(),
  trt_type = character(),
  family = character(),
  intervention = character(),
  variance = character(),
  stringsAsFactors = FALSE
)

# Build the grid manually: only valid combinations.
for (est in c("gcomp", "ipw", "aipw", "matching")) {
  for (trt in c("binary", "binary_binom", "continuous", "categorical")) {
    # Matching only supports binary treatment.
    if (est == "matching" && trt %in% c("continuous", "categorical")) next
    # Matching doesn't support binary outcome well with our DGP.
    if (est == "matching" && trt == "binary_binom") next

    families <- if (trt == "binary_binom") "binomial" else "gaussian"

    for (fam in families) {
      interventions <- switch(
        trt,
        binary = c("static"),
        binary_binom = c("static"),
        continuous = c("shift"),
        categorical = c("static")
      )
      for (iv in interventions) {
        for (var in c("sandwich")) {
          combos <- rbind(combos, data.frame(
            estimator = est,
            trt_type = trt,
            family = fam,
            intervention = iv,
            variance = var,
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  }
}


# -- Run the grid ----------------------------------------------------------

for (i in seq_len(nrow(combos))) {
  row <- combos[i, ]
  label <- paste(row$estimator, row$trt_type, row$family,
                 row$intervention, row$variance, sep = " × ")

  test_that(paste("combo:", label), {
    d <- dgp_generators[[row$trt_type]]()
    info <- truth_table[[row$trt_type]]

    # Build intervention list.
    if (row$intervention == "static" && row$trt_type == "binary") {
      ivs <- list(a1 = static(1), a0 = static(0))
      ref <- "a0"
    } else if (row$intervention == "static" && row$trt_type == "binary_binom") {
      ivs <- list(a1 = static(1), a0 = static(0))
      ref <- "a0"
    } else if (row$intervention == "shift") {
      ivs <- list(shifted = shift(-1), natural = shift(0))
      ref <- "natural"
    } else if (row$intervention == "static" && row$trt_type == "categorical") {
      ivs <- list(set_b = static("b"), set_a = static("a"))
      ref <- "set_a"
    }

    # Fit.
    fit_args <- list(
      data = d,
      outcome = "Y",
      treatment = info$trt_col,
      confounders = info$conf,
      estimator = row$estimator
    )

    if (row$family != "gaussian") {
      fit_args$family <- row$family
    }

    if (row$estimator %in% c("ipw", "aipw") &&
        row$trt_type != "categorical") {
      fit_args$propensity_model_fn <- stats::glm
    }

    fit <- do.call(causat, fit_args)

    # Contrast.
    res <- contrast(
      fit,
      interventions = ivs,
      reference = ref,
      ci_method = row$variance
    )

    est_val <- res$contrasts$estimate[1]
    se_val <- res$contrasts$se[1]

    # Smoke checks: finite estimate.
    expect_true(is.finite(est_val), info = paste("non-finite estimate:", label))

    # SE should be finite and positive.
    if (row$estimator != "matching") {
      expect_true(
        is.finite(se_val) && se_val > 0,
        info = paste("non-finite SE:", label)
      )
    }

    # Point estimate check against truth.
    expect_equal(est_val, info$ate, tolerance = 0.15,
                 label = paste("truth check:", label))
  })
}


# -- Additional family tests on binary treatment ---------------------------

test_that("combo: gcomp × binary × poisson × static × sandwich", {
  set.seed(42)
  L <- rnorm(2000)
  A <- rbinom(2000, 1, plogis(0.5 * L))
  Y <- rpois(2000, exp(0.5 + 0.3 * A + 0.2 * L))
  d <- data.frame(Y = Y, A = A, L = L)

  fit <- causat(d, outcome = "Y", treatment = "A", confounders = ~L,
                family = "poisson")
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich"
  )

  expect_true(is.finite(res$contrasts$estimate[1]))
  expect_true(is.finite(res$contrasts$se[1]) && res$contrasts$se[1] > 0)
})


test_that("combo: gcomp × binary × Gamma × static × sandwich", {
  set.seed(42)
  L <- rnorm(2000)
  A <- rbinom(2000, 1, plogis(0.5 * L))
  mu <- exp(1 + 0.3 * A + 0.2 * L)
  Y <- rgamma(2000, shape = 2, rate = 2 / mu)
  d <- data.frame(Y = Y, A = A, L = L)

  fit <- causat(d, outcome = "Y", treatment = "A", confounders = ~L,
                family = stats::Gamma(link = "log"))
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich"
  )

  expect_true(is.finite(res$contrasts$estimate[1]))
  expect_true(is.finite(res$contrasts$se[1]) && res$contrasts$se[1] > 0)
})


test_that("combo: gcomp × binary × quasibinomial × static × sandwich", {
  set.seed(42)
  L <- rnorm(2000)
  A <- rbinom(2000, 1, plogis(0.5 * L))
  Y <- rbinom(2000, 1, plogis(-1 + 1.5 * A + 0.8 * L))
  d <- data.frame(Y = Y, A = A, L = L)

  fit <- causat(d, outcome = "Y", treatment = "A", confounders = ~L,
                family = "quasibinomial")
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich"
  )

  expect_true(is.finite(res$contrasts$estimate[1]))
  expect_true(is.finite(res$contrasts$se[1]) && res$contrasts$se[1] > 0)
})


# -- Dynamic intervention test -------------------------------------------

test_that("combo: gcomp × binary × dynamic × sandwich", {
  d <- simulate_binary_continuous(n = 2000, seed = 42)

  rule <- dynamic(function(data, trt) ifelse(data$L > 0, 1, 0))

  fit <- causat(d, outcome = "Y", treatment = "A", confounders = ~L)
  res <- contrast(
    fit,
    interventions = list(dyn = rule, a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich"
  )

  expect_true(is.finite(res$contrasts$estimate[1]))
  expect_true(is.finite(res$contrasts$se[1]) && res$contrasts$se[1] > 0)
})


test_that("combo: ipw × binary × dynamic × sandwich", {
  d <- simulate_binary_continuous(n = 2000, seed = 42)

  rule <- dynamic(function(data, trt) ifelse(data$L > 0, 1, 0))

  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw",
    propensity_model_fn = stats::glm
  )
  res <- contrast(
    fit,
    interventions = list(dyn = rule, a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich"
  )

  expect_true(is.finite(res$contrasts$estimate[1]))
  expect_true(is.finite(res$contrasts$se[1]) && res$contrasts$se[1] > 0)
})


# -- scale_by intervention test ------------------------------------------

test_that("combo: gcomp × continuous × scale_by × sandwich", {
  d <- simulate_continuous_continuous(n = 2000, seed = 42)

  fit <- causat(d, outcome = "Y", treatment = "A", confounders = ~L)
  res <- contrast(
    fit,
    interventions = list(doubled = scale_by(2), natural = scale_by(1)),
    reference = "natural",
    ci_method = "sandwich"
  )

  expect_true(is.finite(res$contrasts$estimate[1]))
  expect_true(is.finite(res$contrasts$se[1]) && res$contrasts$se[1] > 0)
})


# -- Aliased (collinear) confounders: sandwich must still produce SE ------

test_that("combo: gcomp × aliased confounders × sandwich", {
  set.seed(42)
  L1 <- rnorm(2000)
  L2 <- L1
  A <- rbinom(2000, 1, plogis(0.5 * L1))
  Y <- 2 + 3 * A + L1 + rnorm(2000)
  d <- data.frame(Y = Y, A = A, L1 = L1, L2 = L2)

  fit <- causat(d, outcome = "Y", treatment = "A", confounders = ~ L1 + L2)
  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich"
  )

  expect_equal(res$contrasts$estimate[1], 3, tolerance = 0.3)
  expect_true(is.finite(res$contrasts$se[1]) && res$contrasts$se[1] > 0)
})


test_that("combo: ipw × aliased confounders × sandwich", {
  set.seed(42)
  L1 <- rnorm(2000)
  L2 <- L1
  A <- rbinom(2000, 1, plogis(0.5 * L1))
  Y <- 2 + 3 * A + L1 + rnorm(2000)
  d <- data.frame(Y = Y, A = A, L1 = L1, L2 = L2)

  expect_warning(
    fit <- causat(d, outcome = "Y", treatment = "A", confounders = ~ L1 + L2,
                  estimator = "ipw", propensity_model_fn = stats::glm),
    "aliased"
  )

  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich"
  )

  expect_equal(res$contrasts$estimate[1], 3, tolerance = 0.3)
  expect_true(is.finite(res$contrasts$se[1]) && res$contrasts$se[1] > 0)
})


test_that("combo: aipw × aliased confounders × sandwich", {
  set.seed(42)
  L1 <- rnorm(2000)
  L2 <- L1
  A <- rbinom(2000, 1, plogis(0.5 * L1))
  Y <- 2 + 3 * A + L1 + rnorm(2000)
  d <- data.frame(Y = Y, A = A, L1 = L1, L2 = L2)

  expect_warning(
    fit <- causat(d, outcome = "Y", treatment = "A", confounders = ~ L1 + L2,
                  estimator = "aipw", propensity_model_fn = stats::glm),
    "aliased"
  )

  res <- contrast(
    fit,
    interventions = list(a1 = static(1), a0 = static(0)),
    reference = "a0",
    ci_method = "sandwich"
  )

  expect_equal(res$contrasts$estimate[1], 3, tolerance = 0.3)
  expect_true(is.finite(res$contrasts$se[1]) && res$contrasts$se[1] > 0)
})


# -- Count treatment tests -----------------------------------------------

test_that("combo: ipw × count (Poisson) × shift × sandwich", {
  d <- simulate_count_treatment(n = 3000, seed = 42)

  fit <- causat(
    d,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw",
    propensity_family = "poisson"
  )
  res <- contrast(
    fit,
    interventions = list(shifted = shift(1), nat = shift(0)),
    reference = "nat",
    ci_method = "sandwich"
  )

  # True shift(1) vs shift(0) difference = 1.5.
  expect_equal(res$contrasts$estimate[1], 1.5, tolerance = 0.1)
  expect_true(is.finite(res$contrasts$se[1]) && res$contrasts$se[1] > 0)
})
