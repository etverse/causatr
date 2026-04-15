#' Fit a treatment density model
#'
#' @description
#' Fits a parametric model for the conditional density
#' \eqn{f(A \mid L)} used by the self-contained IPW engine.
#' The returned `causatr_treatment_model` object exposes enough
#' metadata to (a) evaluate the density at arbitrary treatment
#' values and (b) rebuild those evaluations at a candidate propensity
#' parameter `alpha` — the latter is what the variance engine's
#' `numDeriv::jacobian()` call requires to compute the cross-derivative
#' \eqn{A_{\beta\alpha}}.
#'
#' Three treatment families are handled:
#'
#' - **Binary** (integer / logical 0/1): fit by a logistic GLM via
#'   `model_fn` with `family = stats::binomial()`. Density is the
#'   Bernoulli pmf \eqn{p^a (1-p)^{1-a}}.
#' - **Continuous** (numeric, many distinct values): fit by a Gaussian
#'   linear model via `model_fn` with `family = stats::gaussian()`.
#'   Density is \eqn{\mathcal N(\mu = \hat\mu_L, \sigma = \hat\sigma)},
#'   with \eqn{\hat\sigma} estimated as the residual standard deviation
#'   of the fit. The Gaussian assumption is the standard book recipe
#'   (Hernán & Robins Ch. 12 §12.4); users who need something else can
#'   pass their own `propensity_model_fn` in `causat()`.
#' - **Categorical** (factor or character with \eqn{k > 2} levels):
#'   **not supported**. Aborts.
#'
#' @param data data.table (already prepared by `prepare_data()`) holding
#'   all model variables. Rows with missing treatment or confounder
#'   values must have been removed or imputed upstream.
#' @param treatment Character scalar. Name of the treatment column in
#'   `data`.
#' @param confounders One-sided formula of baseline confounders, e.g.
#'   `~ sex + age + wt71`. `build_ps_formula()` turns this into the
#'   two-sided propensity formula `A ~ confounders`.
#' @param model_fn Function with signature
#'   `function(formula, data, family, weights, ...)`. Defaults to
#'   `stats::glm`. Must return an object supporting `predict(type =
#'   "response")`, `coef()`, `model.matrix()`, and (for continuous
#'   treatments) `residuals(type = "response")`.
#' @param weights Optional numeric vector of observation weights to
#'   forward into `model_fn` (survey weights / external IPCW). `NULL`
#'   means unweighted.
#' @param ... Additional arguments forwarded to `model_fn`.
#'
#' @return A `causatr_treatment_model` S3 object (a list) with slots:
#'   \describe{
#'     \item{`model`}{The fitted density model.}
#'     \item{`family`}{Character tag: `"bernoulli"`, `"gaussian"`, or
#'       `"categorical"`.}
#'     \item{`treatment`}{Treatment column name.}
#'     \item{`ps_formula`}{The `A ~ confounders` formula used.}
#'     \item{`alpha_hat`}{Numeric vector. Fitted propensity parameters
#'       \eqn{\hat\alpha} — coefficients of the density model, used as
#'       the starting point for the variance-engine `numDeriv::jacobian`
#'       call.}
#'     \item{`X_prop`}{Design matrix of the density model. Captured once
#'       here so `make_weight_fn()` can recompute predictions at
#'       candidate `alpha` values without re-forming the matrix.}
#'     \item{`sigma`}{Residual standard deviation (continuous only;
#'       `NULL` for binary / categorical). Treated as fixed at fit time:
#'       the variance engine perturbs `alpha` only, not `sigma`.}
#'     \item{`levels`}{Character vector of levels (categorical only;
#'       `NULL` otherwise).}
#'     \item{`fit_rows`}{Logical vector indicating which rows of `data`
#'       were actually used to fit the density model — the density
#'       evaluator and weight builder align their outputs to these
#'       rows.}
#'   }
#'
#' @references
#' Hernán MA, Robins JM (2025). *Causal Inference: What If*. Chapman &
#' Hall/CRC. Section 12.4 (parametric g-formula for continuous
#' treatments).
#'
#' @noRd
fit_treatment_model <- function(
  data,
  treatment,
  confounders,
  model_fn = stats::glm,
  weights = NULL,
  ...
) {
  # Validate inputs. These are internal helpers so most of the argument
  # shapes are already guaranteed by `check_causat_inputs()` upstream,
  # but the treatment column's *values* still need to be classified
  # into binary/continuous/categorical — that is not something
  # prepare_data() does.
  check_string(treatment)
  check_formula(confounders)
  check_col_exists(data, treatment)

  trt_vals <- data[[treatment]]

  # Rows usable for the propensity fit: non-missing treatment and
  # non-missing confounders. We mirror the outcome-side convention in
  # `get_fit_rows()`: compute a logical row mask, pass data[fit_rows]
  # to the fitter, and remember the mask so density evaluation and the
  # `weight_fn` closure can align their outputs to the same rows. The
  # outcome column is irrelevant here (density is A|L, not Y|A,L).
  confounder_vars <- all.vars(confounders)
  fit_rows <- !is.na(trt_vals)
  for (v in confounder_vars) {
    fit_rows <- fit_rows & !is.na(data[[v]])
  }
  fit_data <- data[fit_rows]

  ps_formula <- build_ps_formula(confounders, treatment)

  family_tag <- detect_treatment_family(trt_vals)

  # Compose model_fn arguments uniformly across the three density
  # families. `weights` is only attached when the user actually passed
  # a vector — otherwise the sub-setted vector would need to be aligned
  # to `fit_rows`, and we let `model_fn`'s own NULL default take over.
  # (Survey-weighted propensity fits pass through this path; the
  # variance engine's IF machinery handles the M-estimation
  # correction via `apply_model_correction()` downstream.)
  weights_fit <- if (is.null(weights)) NULL else weights[fit_rows]

  model <- switch(
    family_tag,
    bernoulli = fit_bernoulli_density(
      ps_formula,
      fit_data,
      model_fn,
      weights_fit,
      ...
    ),
    gaussian = fit_gaussian_density(
      ps_formula,
      fit_data,
      model_fn,
      weights_fit,
      ...
    ),
    categorical = rlang::abort(
      c(
        "Categorical treatment is not supported under `estimator = 'ipw'`.",
        i = "Use `estimator = 'gcomp'` for categorical treatments."
      ),
      class = "causatr_phase4_categorical_pending"
    )
  )

  # For continuous treatments we need the residual SD of the
  # linear propensity model — it enters the Gaussian density as the
  # second parameter. `summary(glm)$sigma` is the standard slot
  # returned by `stats::summary.lm` / `summary.glm` for gaussian fits,
  # which matches the book's "estimate variance of residuals" recipe.
  # For GAM (`mgcv::gam`) the equivalent is `sqrt(model$sig2)`; we
  # cover both shapes defensively.
  sigma <- NULL
  if (family_tag == "gaussian") {
    sigma <- extract_sigma(model)
  }

  # Capture the propensity design matrix ONCE at fit time. The
  # `numDeriv::jacobian()` call inside the variance engine perturbs
  # `alpha` many times and each perturbation needs to recompute
  # `X_prop %*% alpha` — re-calling `model.matrix()` inside the closure
  # would repeat the formula parsing and data.frame coercion on every
  # jacobian step, which is expensive on large datasets.
  X_prop <- stats::model.matrix(model)

  # `coef()` strips NA coefficients on fits with aliased columns; we
  # assume no aliasing (the caller is passing a valid confounder
  # formula). The length check catches "confounder degenerate under
  # the treatment stratum" edge cases early.
  alpha_hat <- stats::coef(model)
  if (ncol(X_prop) != length(alpha_hat)) {
    rlang::abort(
      paste0(
        "Treatment density model has ",
        length(alpha_hat),
        " coefficients but the design matrix has ",
        ncol(X_prop),
        " columns — this usually means a column was aliased or dropped. ",
        "Drop the offending confounder and refit."
      )
    )
  }

  structure(
    list(
      model = model,
      family = family_tag,
      treatment = treatment,
      ps_formula = ps_formula,
      alpha_hat = alpha_hat,
      X_prop = X_prop,
      sigma = sigma,
      levels = NULL,
      fit_rows = fit_rows
    ),
    class = "causatr_treatment_model"
  )
}


#' Classify a treatment vector into a density family
#'
#' @description
#' Three-way dispatch used by `fit_treatment_model()`:
#'
#' - integer / logical / numeric 0/1 with exactly two unique non-NA
#'   values -> `"bernoulli"`
#' - factor (with two OR more levels) or character -> `"categorical"`
#'   (two-level factors / characters are treated as categorical, not
#'   binary, because the model_fn contract for categorical treatments
#'   differs — binary uses `family = binomial()`, categorical uses a
#'   multinomial fitter. Users who want binary handling from a factor
#'   should recode to 0/1 explicitly, which matches the existing
#'   `check_estimand_trt_compat()` convention.)
#' - numeric with more than two distinct values -> `"gaussian"`
#'
#' @param x Treatment vector.
#' @return Character scalar: `"bernoulli"`, `"gaussian"`, or
#'   `"categorical"`.
#'
#' @noRd
detect_treatment_family <- function(x) {
  if (is.factor(x) || is.character(x)) {
    return("categorical")
  }
  if (is.logical(x)) {
    return("bernoulli")
  }
  if (is.numeric(x)) {
    uniq <- unique(stats::na.omit(x))
    if (length(uniq) == 2L && all(uniq %in% c(0, 1))) {
      return("bernoulli")
    }
    return("gaussian")
  }
  rlang::abort(
    paste0(
      "Treatment column of type `",
      typeof(x),
      "` is not supported by the self-contained IPW engine. ",
      "Use a numeric (binary 0/1 or continuous), factor, or character column."
    )
  )
}


#' Fit a Bernoulli density model
#'
#' @description
#' Thin wrapper around `model_fn(A ~ L, family = binomial())`.
#' Factored out so `fit_treatment_model()` stays flat and so the
#' dispatch-by-family structure reads cleanly.
#'
#' @inheritParams fit_treatment_model
#' @param ps_formula Two-sided formula `A ~ confounders`.
#' @param fit_data `data.table` subset already restricted to `fit_rows`.
#' @param weights_fit Numeric weight vector (or `NULL`) aligned with
#'   `fit_data`.
#' @return The fitted model object.
#' @noRd
fit_bernoulli_density <- function(
  ps_formula,
  fit_data,
  model_fn,
  weights_fit,
  ...
) {
  # `model_fn(formula, data, family, weights, ...)` is the causatr
  # contract (see `causat(model_fn =)` docs). Passing `weights = NULL`
  # to `stats::glm()` is fine — it takes the NULL branch internally —
  # but custom fitters in the wild sometimes dispatch on the
  # `weights` argument's presence rather than its value, so we omit
  # the arg entirely when there are no external weights.
  if (is.null(weights_fit)) {
    model_fn(
      formula = ps_formula,
      data = fit_data,
      family = stats::binomial(),
      ...
    )
  } else {
    model_fn(
      formula = ps_formula,
      data = fit_data,
      family = stats::binomial(),
      weights = weights_fit,
      ...
    )
  }
}


#' Fit a Gaussian density model
#'
#' @description
#' Thin wrapper around `model_fn(A ~ L, family = gaussian())`.
#'
#' @inheritParams fit_bernoulli_density
#' @return The fitted model object.
#' @noRd
fit_gaussian_density <- function(
  ps_formula,
  fit_data,
  model_fn,
  weights_fit,
  ...
) {
  if (is.null(weights_fit)) {
    model_fn(
      formula = ps_formula,
      data = fit_data,
      family = stats::gaussian(),
      ...
    )
  } else {
    model_fn(
      formula = ps_formula,
      data = fit_data,
      family = stats::gaussian(),
      weights = weights_fit,
      ...
    )
  }
}


#' Extract the residual standard deviation from a fitted density model
#'
#' @description
#' The "right slot" to read sigma from depends on the fit class:
#'
#' - `stats::lm` / `summary.lm` -> `summary(model)$sigma` (residual SE).
#' - `stats::glm` / `summary.glm` -> `summary(model)$dispersion` is
#'   the estimator of \eqn{\sigma^2}; we take its square root. The
#'   `$sigma` field does NOT exist on `summary.glm` (common gotcha).
#' - `mgcv::gam` -> `model$sig2` is the variance estimate.
#' - Unknown backend -> fall back to computing
#'   \eqn{\hat\sigma = \sqrt{\sum r_i^2 / (n - p)}} from response
#'   residuals, which matches the GLM dispersion MLE for Gaussian.
#'
#' We try each slot in turn rather than branching on the fit class,
#' so the function stays open-ended to any user-supplied
#' `propensity_model_fn`.
#'
#' @param model Fitted density model.
#' @return Positive scalar. Aborts if sigma cannot be recovered.
#' @noRd
extract_sigma <- function(model) {
  sm <- tryCatch(summary(model), error = function(e) NULL)

  # summary.lm$sigma — residual standard error.
  if (!is.null(sm) && !is.null(sm$sigma) && is.finite(sm$sigma)) {
    return(as.numeric(sm$sigma))
  }

  # summary.glm$dispersion — estimator of sigma^2 for Gaussian GLM.
  # `$sigma` is NULL on summary.glm; the dispersion slot is the canonical
  # place to read the residual variance from a Gaussian glm fit.
  if (!is.null(sm) && !is.null(sm$dispersion) && is.finite(sm$dispersion)) {
    return(sqrt(as.numeric(sm$dispersion)))
  }

  # mgcv::gam stores the variance estimate as `$sig2`.
  if (!is.null(model$sig2) && is.finite(model$sig2)) {
    return(sqrt(as.numeric(model$sig2)))
  }

  # Fallback: compute from response residuals directly.
  # `sum(r^2) / (n - p)` is the GLM dispersion estimator under
  # Gaussian, so this matches what `summary.glm$dispersion` would give.
  r <- tryCatch(
    stats::residuals(model, type = "response"),
    error = function(e) NULL
  )
  if (!is.null(r)) {
    n <- length(r)
    p <- length(stats::coef(model))
    df <- max(n - p, 1L)
    return(sqrt(sum(r^2) / df))
  }

  rlang::abort(
    "Could not extract residual SD from the treatment density model."
  )
}


#' Evaluate a treatment density at arbitrary treatment values
#'
#' @description
#' Computes \eqn{f(a_i \mid L_i)} for user-supplied treatment values
#' \eqn{a_i}. This is the primitive that `compute_density_ratio_weights()`
#' uses to build \eqn{w_i = f(d(A_i, L_i) \mid L_i) / f(A_i \mid L_i)}.
#'
#' ## Row alignment contract
#'
#' The density is evaluated on `newdata`, which may be shorter than
#' the data the treatment model was fit on (for example, if the caller
#' passes `fit$data[fit_rows]`). `treatment_values` must have
#' `nrow(newdata)` elements. The returned density vector has the same
#' length.
#'
#' @param treatment_model A `causatr_treatment_model` from
#'   `fit_treatment_model()`.
#' @param treatment_values Numeric (binary / continuous) or character /
#'   factor (categorical) vector of length `nrow(newdata)`.
#' @param newdata Data frame or data.table of confounder values. Must
#'   contain all columns referenced by `treatment_model$ps_formula`.
#'
#' @return Numeric vector of density values, length `nrow(newdata)`.
#'
#' @noRd
evaluate_density <- function(treatment_model, treatment_values, newdata) {
  if (!inherits(treatment_model, "causatr_treatment_model")) {
    rlang::abort(
      "`treatment_model` must be a `causatr_treatment_model` from `fit_treatment_model()`."
    )
  }
  if (length(treatment_values) != nrow(newdata)) {
    rlang::abort(
      paste0(
        "`treatment_values` has length ",
        length(treatment_values),
        " but `newdata` has ",
        nrow(newdata),
        " rows — they must be the same length."
      )
    )
  }

  model <- treatment_model$model
  family_tag <- treatment_model$family

  if (family_tag == "bernoulli") {
    # Bernoulli pmf: f(a | L) = p^a * (1-p)^(1-a). For numeric 0/1
    # treatments this collapses to either `p` or `1-p` depending on
    # which branch a_i is in, so we use `ifelse` rather than the
    # explicit power form to avoid `0^0` quirks when p is exactly 0 or 1.
    p <- stats::predict(model, newdata = newdata, type = "response")
    return(ifelse(treatment_values == 1, p, 1 - p))
  }

  if (family_tag == "gaussian") {
    # Normal pdf: f(a | L) = dnorm(a, mu, sigma) with mu = fitted mean
    # and sigma = residual SD (fixed at fit time). `dnorm` is
    # vectorised — one call per density evaluation, not per row.
    mu <- stats::predict(model, newdata = newdata, type = "response")
    sigma <- treatment_model$sigma
    return(stats::dnorm(treatment_values, mean = mu, sd = sigma))
  }

  rlang::abort(
    paste0("Unknown treatment family '", family_tag, "'.")
  )
}


#' Print a causatr_treatment_model
#'
#' @param x A `causatr_treatment_model` object.
#' @param ... Currently unused.
#' @return Invisibly returns `x`.
#' @noRd
print.causatr_treatment_model <- function(x, ...) {
  cat("<causatr_treatment_model: ", x$family, ">\n", sep = "")
  cat("  treatment: ", x$treatment, "\n", sep = "")
  cat("  n_fit:     ", sum(x$fit_rows), "\n", sep = "")
  cat("  p_alpha:   ", length(x$alpha_hat), "\n", sep = "")
  if (!is.null(x$sigma)) {
    cat("  sigma:     ", format(x$sigma, digits = 4), "\n", sep = "")
  }
  invisible(x)
}
