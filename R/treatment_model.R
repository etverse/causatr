#' Fit a treatment density model
#'
#' @description
#' Fits a parametric model for the conditional density
#' \eqn{f(A \mid L)} used by the self-contained IPW engine.
#' The returned `causatr_treatment_model` object exposes enough
#' metadata to (a) evaluate the density at arbitrary treatment
#' values and (b) rebuild those evaluations at a candidate propensity
#' parameter `alpha` -- the latter is what the variance engine's
#' `numDeriv::jacobian()` call requires to compute the cross-derivative
#' \eqn{A_{\beta\alpha}}.
#'
#' Five treatment families are handled:
#'
#' - **Binary** (integer / logical 0/1): fit by a logistic GLM via
#'   `model_fn` with `family = stats::binomial()`. Density is the
#'   Bernoulli pmf \eqn{p^a (1-p)^{1-a}}.
#' - **Continuous** (numeric, many distinct values): fit by a Gaussian
#'   linear model via `model_fn` with `family = stats::gaussian()`.
#'   Density is \eqn{\mathcal N(\mu = \hat\mu_L, \sigma = \hat\sigma)},
#'   with \eqn{\hat\sigma} estimated as the residual standard deviation
#'   of the fit. The Gaussian assumption is the standard book recipe
#'   (Hernan & Robins Ch. 12 Sec. 12.4); users who need something else can
#'   pass their own `propensity_model_fn` in `causat()`.
#' - **Categorical** (factor or character with \eqn{k \ge 2} levels):
#'   fit by a multinomial logistic model via `model_fn` with **no
#'   `family` argument** (the categorical `model_fn` must accept
#'   `function(formula, data, weights, ...)`; the default is
#'   `nnet::multinom`). Density is the multinomial pmf
#'   \eqn{p_k = P(A = k \mid L)}, read from
#'   `predict(model, type = "probs")`.
#' - **Poisson** (opt-in via `propensity_family = "poisson"`): fit by
#'   a Poisson GLM via `model_fn` with `family = stats::poisson()`.
#'   Density is the Poisson pmf
#'   \eqn{f(a \mid L) = \lambda^a e^{-\lambda} / a!} with
#'   \eqn{\lambda = \exp(X \hat\alpha)}.
#' - **Negative binomial** (opt-in via `propensity_family = "negbin"`):
#'   fit by `MASS::glm.nb()` (auto-selected when
#'   `propensity_model_fn = NULL`). Density is the NB pmf with
#'   \eqn{\mu = \exp(X \hat\alpha)} and dispersion parameter
#'   \eqn{\theta} estimated by MLE. \eqn{\theta} is treated as fixed
#'   in the variance engine (perturbed via `alpha` only, not `theta`),
#'   matching the convention for \eqn{\sigma} in the Gaussian case.
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
#' @param propensity_family Character or `NULL`. Explicit treatment
#'   density family override. `NULL` (default) auto-detects from the
#'   treatment values via `detect_treatment_family()`. Pass `"poisson"`
#'   or `"negbin"` to opt into a count regression for integer-valued
#'   treatments (dose levels, event counts, visit counts). Auto-detect
#'   never infers count -- non-negative integers like age in years are
#'   legitimately modelled as Gaussian.
#' @param weights Optional numeric vector of observation weights to
#'   forward into `model_fn` (survey weights / external IPCW). `NULL`
#'   means unweighted.
#' @param ... Additional arguments forwarded to `model_fn`.
#'
#' @return A `causatr_treatment_model` S3 object (a list) with slots:
#'   \describe{
#'     \item{`model`}{The fitted density model.}
#'     \item{`family`}{Character tag: `"bernoulli"`, `"gaussian"`,
#'       `"categorical"`, `"poisson"`, or `"negbin"`.}
#'     \item{`treatment`}{Treatment column name.}
#'     \item{`ps_formula`}{The `A ~ confounders` formula used.}
#'     \item{`alpha_hat`}{Numeric vector. Fitted propensity parameters
#'       \eqn{\hat\alpha} -- coefficients of the density model, used as
#'       the starting point for the variance-engine `numDeriv::jacobian`
#'       call.}
#'     \item{`X_prop`}{Design matrix of the density model. Captured once
#'       here so `make_weight_fn()` can recompute predictions at
#'       candidate `alpha` values without re-forming the matrix.}
#'     \item{`sigma`}{Residual standard deviation (continuous only;
#'       `NULL` for binary / categorical / count). Treated as fixed at
#'       fit time: the variance engine perturbs `alpha` only, not
#'       `sigma`.}
#'     \item{`theta`}{NB dispersion parameter (negbin only; `NULL`
#'       otherwise). Treated as fixed at fit time, matching the
#'       `sigma` convention.}
#'     \item{`levels`}{Character vector of levels (categorical only;
#'       `NULL` otherwise).}
#'     \item{`fit_rows`}{Logical vector indicating which rows of `data`
#'       were actually used to fit the density model -- the density
#'       evaluator and weight builder align their outputs to these
#'       rows.}
#'   }
#'
#' @references
#' Hernan MA, Robins JM (2025). *Causal Inference: What If*. Chapman &
#' Hall/CRC. Section 12.4 (parametric g-formula for continuous
#' treatments).
#'
#' @noRd
fit_treatment_model <- function(
  data,
  treatment,
  confounders,
  model_fn = stats::glm,
  propensity_family = NULL,
  weights = NULL,
  ...
) {
  # Validate inputs. These are internal helpers so most of the argument
  # shapes are already guaranteed by `check_causat_inputs()` upstream,
  # but the treatment column's *values* still need to be classified
  # into binary/continuous/categorical -- that is not something
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

  # When the user declares a `propensity_family`, it overrides
  # auto-detection. This is the only path to Poisson / NB treatment
  # densities, because auto-detect never infers count (too many
  # legitimate non-count integers: age, years of education, etc.).
  if (!is.null(propensity_family)) {
    propensity_family <- match.arg(
      propensity_family,
      c("poisson", "negbin")
    )
    family_tag <- propensity_family
  } else {
    family_tag <- detect_treatment_family(trt_vals)
  }

  # Compose model_fn arguments uniformly across the density
  # families. `weights` is only attached when the user actually passed
  # a vector -- otherwise the sub-setted vector would need to be aligned
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
    categorical = fit_categorical_density(
      ps_formula,
      fit_data,
      model_fn,
      weights_fit,
      ...
    ),
    poisson = fit_count_density(
      ps_formula,
      fit_data,
      model_fn,
      weights_fit,
      family_tag = "poisson",
      ...
    ),
    negbin = fit_count_density(
      ps_formula,
      fit_data,
      model_fn,
      weights_fit,
      family_tag = "negbin",
      ...
    )
  )

  # For continuous treatments we need the residual SD of the
  # linear propensity model -- it enters the Gaussian density as the
  # second parameter. `summary(glm)$sigma` is the standard slot
  # returned by `stats::summary.lm` / `summary.glm` for gaussian fits,
  # which matches the book's "estimate variance of residuals" recipe.
  # For GAM (`mgcv::gam`) the equivalent is `sqrt(model$sig2)`; we
  # cover both shapes defensively.
  sigma <- NULL
  if (family_tag == "gaussian") {
    sigma <- extract_sigma(model)
  }

  # For negative binomial, extract the dispersion parameter theta.
  # Treated as fixed in the variance engine (perturbing only alpha,
  # not theta) -- same convention as sigma for Gaussian.
  theta <- NULL
  if (family_tag == "negbin") {
    theta <- model$theta
    if (is.null(theta) || !is.finite(theta) || theta <= 0) {
      rlang::abort(
        "Could not extract dispersion parameter `theta` from the negative binomial model."
      )
    }
  }

  # Capture the propensity design matrix ONCE at fit time. The
  # `numDeriv::jacobian()` call inside the variance engine perturbs
  # `alpha` many times and each perturbation needs to recompute
  # `X_prop %*% alpha` -- re-calling `model.matrix()` inside the closure
  # would repeat the formula parsing and data.frame coercion on every
  # jacobian step, which is expensive on large datasets.
  X_prop <- stats::model.matrix(model)

  # For multinomial models, `coef()` returns a (K-1) x p matrix where

  # K is the number of levels and p is the number of design-matrix
  # columns. The variance engine needs a flat vector so it can perturb
  # each scalar entry independently via `numDeriv::jacobian()`. We
  # flatten row-major: `as.vector(t(coef_mat))` interleaves
  # (intercept_2, beta_L_2, intercept_3, beta_L_3, ...) so that the
  # reshaping in `make_weight_fn()`'s categorical closure can recover
  # the matrix via `matrix(alpha, nrow = K-1, ncol = p, byrow = TRUE)`.
  # For Bernoulli / Gaussian, `coef()` already returns a vector.
  alpha_hat <- stats::coef(model)
  if (family_tag == "categorical") {
    # `nnet::multinom` with 2 levels returns a plain vector (not a
    # matrix). Normalise to a 1-row matrix so downstream code always
    # sees a consistent shape.
    if (is.null(dim(alpha_hat))) {
      alpha_hat <- matrix(alpha_hat, nrow = 1L)
    }
    alpha_hat <- as.vector(t(alpha_hat))
  }

  if (ncol(X_prop) * n_alpha_rows(family_tag, model) != length(alpha_hat)) {
    rlang::abort(
      paste0(
        "Treatment density model has ",
        length(alpha_hat),
        " coefficients but the design matrix has ",
        ncol(X_prop),
        " columns (expected ",
        ncol(X_prop) * n_alpha_rows(family_tag, model),
        " total) -- this usually means a column was aliased or dropped. ",
        "Drop the offending confounder and refit."
      )
    )
  }

  # Capture treatment levels for categorical (needed by
  # `evaluate_density()` and `make_weight_fn()` to map predicted
  # probability columns back to treatment values).
  trt_levels <- NULL
  if (family_tag == "categorical") {
    trt_levels <- levels(factor(trt_vals))
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
      theta = theta,
      levels = trt_levels,
      fit_rows = fit_rows
    ),
    class = "causatr_treatment_model"
  )
}


#' Fit a sequence of conditional treatment density models for multivariate IPW
#'
#' @description
#' For multivariate treatments `A = (A_1, ..., A_K)` we factorise the
#' joint conditional density via the chain rule:
#' \deqn{f(A_1, A_2, \ldots, A_K \mid L)
#'       = f(A_1 \mid L) \cdot f(A_2 \mid A_1, L)
#'         \cdots f(A_K \mid A_1, \ldots, A_{K-1}, L).}
#' Each factor is fit by `fit_treatment_model()` on a per-component
#' formula that adds the upstream treatment columns to the user's
#' baseline confounders. The result is a list of K
#' `causatr_treatment_model` objects, in the natural order of
#' `treatment`. Downstream the IPW engine multiplies the K density
#' ratios to form the joint weight under any factorising
#' intervention `\pi = \prod_k \pi_k(a_k \mid a_{1..k-1}, L)`.
#'
#' Mixed treatment families are allowed: each component is classified
#' independently by `detect_treatment_family()`, so e.g. the user can
#' pass `treatment = c("A1", "A2")` with `A1` binary 0/1 and `A2`
#' continuous and the engine fits a logistic GLM for `A1` and a
#' Gaussian linear model for `A2 \mid A1, L`.
#'
#' Categorical, Poisson, and negative-binomial components are rejected
#' for now -- their dispatch needs the per-component
#' `propensity_model_fn` / `propensity_family` plumbing the v1
#' multivariate path does not yet expose. Single-component categorical
#' / count IPW continues to work via `fit_treatment_model()`.
#'
#' @param data data.table (already prepared by `prepare_data()`).
#' @param treatment Character vector of treatment column names, in the
#'   factorisation order. Length must be \eqn{\geq 2}; for length 1
#'   call `fit_treatment_model()` directly.
#' @param confounders One-sided formula of baseline confounders. The
#'   k-th component formula is built as
#'   `A_k ~ A_1 + ... + A_{k-1} + confounders`.
#' @param model_fn Function. The default fitter for binary / continuous
#'   / count components. Default `stats::glm`. Forwarded to
#'   `fit_treatment_model()` per component as `model_fn` unless the
#'   component is categorical or count and `propensity_model_fn` /
#'   `propensity_family` triggers a different choice.
#' @param propensity_model_fn Optional function or `NULL`. If non-`NULL`,
#'   used as the per-component fitter for ALL components. Useful when
#'   the user wants e.g. `mgcv::gam` for every propensity. `NULL`
#'   (default) lets the engine auto-pick per component (multinomial
#'   for categorical, GLM for everything else).
#' @param propensity_family Optional character vector or `NULL`. Per-
#'   component opt-in to count families (`"poisson"` / `"negbin"`).
#'   `NULL` (default) lets the engine auto-detect each component from
#'   the column type. Length 1 broadcasts to all K components; length
#'   K applies per-component (entries `""` / `NA` skip auto-detect).
#' @param weights Optional numeric vector of observation weights,
#'   forwarded to each component's fit.
#' @param ... Additional arguments forwarded to the component fitter.
#'
#' @return A length-K list of `causatr_treatment_model` objects, one
#'   per treatment component, in the order given by `treatment`. The
#'   list itself is tagged with class
#'   `c("causatr_treatment_models", "list")` so callers can dispatch
#'   on shape without needing a length check.
#'
#' @noRd
fit_treatment_models <- function(
  data,
  treatment,
  confounders,
  model_fn = stats::glm,
  propensity_model_fn = NULL,
  propensity_family = NULL,
  stabilize = "none",
  weights = NULL,
  ...
) {
  check_formula(confounders)
  if (!is.character(treatment) || length(treatment) < 2L) {
    rlang::abort(
      "`fit_treatment_models()` requires a length-2+ character vector of treatment column names."
    )
  }
  for (trt in treatment) {
    check_col_exists(data, trt)
  }

  # The k-th component conditions on the prior treatment columns plus
  # the user's baseline confounders. Build the per-component formula
  # by appending `treatment[1..k-1]` to the confounder term labels and
  # calling `stats::reformulate()`. The k-th treatment is the response.
  #
  # Strip any term that touches a treatment column from the baseline
  # term labels. This covers `A1:sex` (effect-modification interaction
  # with a baseline modifier) and `A1:A2` (treatment-treatment
  # interaction): both belong on the OUTCOME model side, not on the
  # propensity side, since A appears as the response of the per-
  # component density model and downstream A_j (j > k) cannot enter
  # the k-th propensity by construction. `parse_effect_mod()` returns
  # `confounder_terms` with all such treatment-touching terms removed,
  # which is exactly what we need.
  em_info <- parse_effect_mod(confounders, treatment)
  baseline_terms <- em_info$confounder_terms
  if (length(baseline_terms) == 0L) {
    baseline_terms <- "1"
  }

  K <- length(treatment)
  models <- vector("list", K)

  # Resolve per-component `propensity_family` opt-in. Three valid
  # shapes:
  #   - NULL: per-component auto-detect.
  #   - length 1: broadcast to all K components (e.g. all-Poisson).
  #   - length K: per-component override; entries `""` or NA fall back
  #     to auto-detect.
  if (is.null(propensity_family)) {
    pf_vec <- rep(NA_character_, K)
  } else if (length(propensity_family) == 1L) {
    pf_vec <- rep(propensity_family, K)
  } else if (length(propensity_family) == K) {
    pf_vec <- as.character(propensity_family)
  } else {
    rlang::abort(
      paste0(
        "`propensity_family` must be NULL, length 1, or length K = ",
        K,
        " (got length ",
        length(propensity_family),
        ")."
      )
    )
  }

  for (k in seq_along(treatment)) {
    # Component k formula: A_k ~ A_1 + ... + A_{k-1} + baseline_terms.
    # Prior treatments enter as additional conditioning columns. When
    # k == 1 the formula reduces to A_1 ~ baseline_terms.
    if (k == 1L) {
      rhs_terms <- baseline_terms
    } else {
      rhs_terms <- c(treatment[seq_len(k - 1L)], baseline_terms)
    }
    confounders_k <- stats::reformulate(rhs_terms)

    # Per-component family detection + fitter selection. The user can
    # override the global fitter for ALL components via
    # `propensity_model_fn`; otherwise we pick per-component:
    #   categorical -> `nnet::multinom` (default categorical fitter)
    #   negbin     -> `MASS::glm.nb` (no `family` arg accepted by NB)
    #   bernoulli / gaussian / poisson -> `model_fn` (typically
    #     `stats::glm`)
    fam_k <- if (!is.na(pf_vec[k]) && nzchar(pf_vec[k])) {
      pf_vec[k]
    } else {
      detect_treatment_family(data[[treatment[k]]])
    }
    fitter_k <- if (!is.null(propensity_model_fn)) {
      propensity_model_fn
    } else if (fam_k == "categorical") {
      nnet::multinom
    } else if (fam_k == "negbin") {
      check_pkg("MASS")
      MASS::glm.nb
    } else {
      model_fn
    }
    pf_arg_k <- if (!is.na(pf_vec[k]) && nzchar(pf_vec[k])) pf_vec[k] else NULL

    models[[k]] <- fit_treatment_model(
      data = data,
      treatment = treatment[k],
      confounders = confounders_k,
      model_fn = fitter_k,
      propensity_family = pf_arg_k,
      weights = weights,
      ...
    )
  }

  names(models) <- treatment

  # Stabilized weights (sequential MTP; Robins, Hernan, Brumback 2000
  # extended to multivariate). Under `stabilize = "marginal"`, the
  # per-component numerator factor is swapped for a density that drops
  # `L` from the conditioning -- `g_k(A_k | A_{1..k-1})`. The
  # denominator stays at the full-L density. The per-component weight
  # becomes
  #   w_k^s = g_k(d_k^{-1}(A_k) | A_{1..k-1}) * |Jac| /
  #           f_k(A_k | A_{1..k-1}, L)
  # This dampens the multiplicative L-dependence across K factors and
  # typically reduces weight variance.
  num_models <- NULL
  if (stabilize == "marginal") {
    num_models <- vector("list", K)
    for (k in seq_along(treatment)) {
      # Numerator g_k(A_k | A_{1..k-1}) conditions only on prior
      # treatments (no L / baseline confounders). For k = 1 the
      # conditioning is empty -> intercept-only model.
      if (k == 1L) {
        num_rhs <- "1"
      } else {
        num_rhs <- treatment[seq_len(k - 1L)]
      }
      num_confounders_k <- stats::reformulate(num_rhs)

      # Reuse the same per-component family / fitter dispatch as the
      # denominator model -- the numerator density has the same
      # marginal support structure (bernoulli stays bernoulli, count
      # stays count, etc.).
      fam_k <- if (!is.na(pf_vec[k]) && nzchar(pf_vec[k])) {
        pf_vec[k]
      } else {
        detect_treatment_family(data[[treatment[k]]])
      }
      fitter_k <- if (!is.null(propensity_model_fn)) {
        propensity_model_fn
      } else if (fam_k == "categorical") {
        nnet::multinom
      } else if (fam_k == "negbin") {
        check_pkg("MASS")
        MASS::glm.nb
      } else {
        model_fn
      }
      pf_arg_k <- if (!is.na(pf_vec[k]) && nzchar(pf_vec[k])) {
        pf_vec[k]
      } else {
        NULL
      }

      num_models[[k]] <- fit_treatment_model(
        data = data,
        treatment = treatment[k],
        confounders = num_confounders_k,
        model_fn = fitter_k,
        propensity_family = pf_arg_k,
        weights = weights,
        ...
      )
    }
    names(num_models) <- treatment
    class(num_models) <- c("causatr_treatment_models", "list")
  }

  structure(
    models,
    class = c("causatr_treatment_models", "list"),
    numerator_models = num_models,
    stabilize = stabilize
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
#'   differs -- binary uses `family = binomial()`, categorical uses a
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
  # to `stats::glm()` is fine -- it takes the NULL branch internally --
  # but custom fitters in the wild sometimes dispatch on the
  # `weights` argument's presence rather than its value, so we omit
  # the arg entirely when there are no external weights.
  # Build the call via `do.call()` with the weight vector substituted
  # as an actual value rather than a symbol. Passing `weights =
  # weights_fit` positionally would leave the symbol `weights_fit` in
  # the glm call, which `stats::model.frame.default()` then re-evaluates
  # against later newdata frames (e.g. inside `stats::predict()`) and
  # fails with "object 'weights_fit' not found". `do.call()` with a
  # list of actual values avoids the lookup.
  base_args <- list(
    formula = ps_formula,
    data = fit_data,
    family = stats::binomial()
  )
  if (!is.null(weights_fit)) {
    base_args$weights <- weights_fit
  }
  do.call(model_fn, c(base_args, list(...)))
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
  # See the `fit_bernoulli_density()` note on why this uses
  # `do.call()` instead of a direct `model_fn(..., weights =
  # weights_fit, ...)` call.
  base_args <- list(
    formula = ps_formula,
    data = fit_data,
    family = stats::gaussian()
  )
  if (!is.null(weights_fit)) {
    base_args$weights <- weights_fit
  }
  do.call(model_fn, c(base_args, list(...)))
}


#' Fit a multinomial (categorical) density model
#'
#' @description
#' Thin wrapper around `model_fn(A ~ L, ...)` for categorical
#' treatments. Unlike the Bernoulli and Gaussian helpers, the
#' multinomial fitter does **not** receive a `family` argument
#' because `nnet::multinom` (the default categorical fitter)
#' does not accept one. The returned model must support
#' `predict(model, newdata, type = "probs")` returning an
#' n x K matrix of category probabilities (or a length-n vector
#' for the K = 2 case, giving P(Y = second level)).
#'
#' @inheritParams fit_bernoulli_density
#' @return The fitted model object.
#' @noRd
fit_categorical_density <- function(
  ps_formula,
  fit_data,
  model_fn,
  weights_fit,
  ...
) {
  # Multinomial fitters (nnet::multinom, VGAM::vglm) do not take a
  # `family` argument -- the family is implicit in the model class.
  # We therefore build the call without `family`, unlike the Bernoulli
  # and Gaussian helpers.
  base_args <- list(
    formula = ps_formula,
    data = fit_data
  )
  if (!is.null(weights_fit)) {
    base_args$weights <- weights_fit
  }
  # Suppress the iteration trace that nnet::multinom prints by default
  # via `trace = FALSE`. If the user's `model_fn` does not accept
  # `trace`, it will silently absorb it through `...`.
  extra <- list(...)
  if (!"trace" %in% names(extra)) {
    extra$trace <- FALSE
  }
  do.call(model_fn, c(base_args, extra))
}


#' Fit a count (Poisson or negative binomial) density model
#'
#' @description
#' Wrapper for count treatment densities. For Poisson, uses `model_fn`
#' with `family = stats::poisson()`. For negative binomial, uses
#' `MASS::glm.nb` (or the user's `model_fn` if explicitly provided),
#' which does **not** accept a `family` argument.
#'
#' @inheritParams fit_bernoulli_density
#' @param family_tag Character. `"poisson"` or `"negbin"`.
#' @return The fitted model object.
#' @noRd
fit_count_density <- function(
  ps_formula,
  fit_data,
  model_fn,
  weights_fit,
  family_tag,
  ...
) {
  if (family_tag == "poisson") {
    base_args <- list(
      formula = ps_formula,
      data = fit_data,
      family = stats::poisson()
    )
    if (!is.null(weights_fit)) {
      base_args$weights <- weights_fit
    }
    return(do.call(model_fn, c(base_args, list(...))))
  }

  # Negative binomial: `MASS::glm.nb` does not take a `family` arg
  # (the NB family is implicit). Same pattern as categorical /
  # `nnet::multinom`.
  base_args <- list(
    formula = ps_formula,
    data = fit_data
  )
  if (!is.null(weights_fit)) {
    base_args$weights <- weights_fit
  }
  do.call(model_fn, c(base_args, list(...)))
}


#' Number of coefficient rows per design-matrix column
#'
#' @description
#' For Bernoulli / Gaussian models `coef()` returns p scalars (one
#' per design-matrix column). For a multinomial model with K levels,
#' `coef()` returns a (K-1) x p matrix -- K-1 log-odds equations,
#' each with p coefficients. This helper returns 1 for non-
#' categorical fits and K-1 for categorical, so the length check in
#' `fit_treatment_model()` can validate the flattened `alpha_hat`
#' length against `ncol(X_prop)`.
#'
#' @param family_tag Character. `"bernoulli"`, `"gaussian"`, or
#'   `"categorical"`.
#' @param model Fitted model object.
#' @return Positive integer.
#' @noRd
n_alpha_rows <- function(family_tag, model) {
  if (family_tag != "categorical") {
    return(1L)
  }
  cc <- stats::coef(model)
  if (is.null(dim(cc))) {
    # 2-level multinomial: coef is a plain vector (1 equation)
    return(1L)
  }
  nrow(cc)
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

  # summary.lm$sigma -- residual standard error.
  if (!is.null(sm) && !is.null(sm$sigma) && is.finite(sm$sigma)) {
    return(as.numeric(sm$sigma))
  }

  # summary.glm$dispersion -- estimator of sigma^2 for Gaussian GLM.
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
        " rows -- they must be the same length."
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
    # vectorised -- one call per density evaluation, not per row.
    mu <- stats::predict(model, newdata = newdata, type = "response")
    sigma <- treatment_model$sigma
    return(stats::dnorm(treatment_values, mean = mu, sd = sigma))
  }

  if (family_tag == "categorical") {
    # Multinomial pmf: f(a_i | L_i) = P(A = a_i | L_i), read from
    # the predicted probability matrix. `predict(type = "probs")`
    # returns an n x K matrix for K > 2 levels, or a length-n vector
    # for K = 2 (giving P(second level)). We normalise both shapes
    # into an n x K matrix with columns named by the factor levels.
    return(evaluate_categorical_density(
      model,
      treatment_model$levels,
      treatment_values,
      newdata
    ))
  }

  if (family_tag == "poisson") {
    # Poisson pmf: f(a | L) = dpois(a, lambda) with lambda = E[A|L].
    lambda <- stats::predict(model, newdata = newdata, type = "response")
    return(stats::dpois(treatment_values, lambda))
  }

  if (family_tag == "negbin") {
    # Negative binomial pmf: f(a | L) = dnbinom(a, mu = lambda,
    # size = theta) with lambda = E[A|L] and theta estimated by MLE
    # at fit time.
    lambda <- stats::predict(model, newdata = newdata, type = "response")
    theta <- treatment_model$theta
    return(stats::dnbinom(treatment_values, mu = lambda, size = theta))
  }

  rlang::abort(
    paste0("Unknown treatment family '", family_tag, "'.")
  )
}


#' Evaluate the multinomial density at per-row treatment values
#'
#' @description
#' Given a fitted multinomial model and a vector of treatment values
#' (character or factor), returns \eqn{P(A = a_i \mid L_i)} for each
#' row. The predicted probability matrix from `predict(model,
#' type = "probs")` is indexed per row to extract the column matching
#' `treatment_values[i]`.
#'
#' Handles the `nnet::multinom` K = 2 edge case where `predict()`
#' returns a vector (P(second level)) rather than a matrix.
#'
#' @param model Fitted multinomial model.
#' @param trt_levels Character vector of all factor levels (in order).
#' @param treatment_values Character or factor vector of length
#'   `nrow(newdata)`.
#' @param newdata Data frame of confounder values.
#'
#' @return Numeric vector of densities, length `nrow(newdata)`.
#'
#' @noRd
evaluate_categorical_density <- function(
  model,
  trt_levels,
  treatment_values,
  newdata
) {
  n <- nrow(newdata)
  K <- length(trt_levels)

  prob_raw <- stats::predict(model, newdata = newdata, type = "probs")

  # Normalise to an n x K matrix. nnet::multinom returns a vector for
  # K = 2 giving P(second level), and a matrix for K > 2.
  if (is.null(dim(prob_raw))) {
    # K = 2: prob_raw is P(level_2). Build the full n x 2 matrix.
    prob_mat <- cbind(1 - prob_raw, prob_raw)
    colnames(prob_mat) <- trt_levels
  } else {
    prob_mat <- prob_raw
  }

  # Look up the probability corresponding to each row's treatment
  # value. Convert to character for uniform matching against column
  # names.
  tv_char <- as.character(treatment_values)
  col_idx <- match(tv_char, colnames(prob_mat))
  if (anyNA(col_idx)) {
    bad <- unique(tv_char[is.na(col_idx)])
    rlang::abort(
      paste0(
        "Treatment value(s) not found in model levels: ",
        paste(shQuote(bad), collapse = ", "),
        ". Known levels: ",
        paste(shQuote(colnames(prob_mat)), collapse = ", "),
        "."
      )
    )
  }

  # Row-wise indexing: prob_mat[i, col_idx[i]].
  prob_mat[cbind(seq_len(n), col_idx)]
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
  if (!is.null(x$theta)) {
    cat("  theta:     ", format(x$theta, digits = 4), "\n", sep = "")
  }
  if (!is.null(x$levels)) {
    cat("  levels:    ", paste(x$levels, collapse = ", "), "\n", sep = "")
  }
  invisible(x)
}
