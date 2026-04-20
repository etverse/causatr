#' Fit an IPW model for causal estimation (point treatment)
#'
#' @description
#' Implements inverse probability weighting (Hernan & Robins Ch. 12) via a
#' self-contained density-ratio engine. `fit_ipw()` owns only the
#' propensity / treatment-density model fit; the weighted marginal
#' structural model (MSM) is rebuilt **per intervention** inside
#' `contrast()` because the density-ratio weight vector depends on the
#' intervention being contrasted.
#'
#' ## Algorithm
#'
#' 1. Fit the conditional treatment density \eqn{f(A \mid L)} via
#'    `fit_treatment_model()`. The engine dispatches on
#'    `detect_treatment_family()`:
#'    - **Bernoulli** (0/1 numeric): logistic GLM.
#'    - **Gaussian** (continuous numeric with > 2 distinct values):
#'      linear GLM with residual SD.
#'    - **Categorical** (factor / character): multinomial logistic via
#'      `nnet::multinom` (auto-selected when `propensity_model_fn` is
#'      `NULL`; the user can override with any multinomial fitter).
#' 2. Stash the treatment model in `fit$details$treatment_model` so
#'    `contrast()` can build intervention-specific density-ratio
#'    weight vectors via `compute_density_ratio_weights()` /
#'    `make_weight_fn()` and refit a weighted MSM per intervention.
#' 3. Store a cheap placeholder outcome model (`Y ~ A`, unweighted)
#'    in `fit$model` for display / compatibility with `print()` /
#'    `summary()`. This model is **not** used for estimation -- every
#'    contrast refits its own weighted MSM.
#'
#' ## Why no top-level weighted MSM is fit here
#'
#' The density-ratio weight vector depends on the intervention:
#' `static(1)` uses HT weights `I(A=1)/p`, `shift(d)` uses a smooth
#' Gaussian ratio, `ipsi(delta)` uses Kennedy's closed form, etc. A
#' single "fit-time" MSM would have to pick one intervention and fix
#' it. `contrast()` handles the per-intervention refit and its own
#' sandwich machinery (`variance_if_ipw()` calling
#' `compute_ipw_if_self_contained_one()`), leaving `fit_ipw()` free to
#' own only the propensity / density side.
#'
#' @param data data.table from `prepare_data()`.
#' @param outcome Character. Outcome column name.
#' @param treatment Character scalar or vector. Treatment column name(s).
#' @param confounders One-sided formula of baseline confounders.
#' @param confounders_tv One-sided formula of time-varying confounders or `NULL`.
#' @param estimand Character. `"ATE"`, `"ATT"`, or `"ATC"`.
#' @param type Character. `"point"` or `"longitudinal"`.
#' @param history Positive integer or `Inf`. Markov order for longitudinal.
#' @param numerator One-sided formula or `NULL`. Numerator formula for
#'   stabilized weights in longitudinal models.
#' @param weights Numeric vector or `NULL`. External observation weights
#'   (e.g. survey weights), multiplied with the estimated IPW weights
#'   per intervention inside `contrast()`.
#' @param model_fn Function. The outcome-model fitting function used
#'   to build the placeholder `Y ~ A` model. Must accept
#'   `(formula, data, family, weights, ...)`. Default `stats::glm`.
#' @param propensity_model_fn Function or `NULL`. The fitting function
#'   used for the treatment density model. Must accept the same
#'   signature as `model_fn`. `NULL` (default) reuses `model_fn`.
#'   Common non-default choice: `mgcv::gam` for a flexible propensity.
#' @param call The original `causat()` call.
#' @param ... Passed to `propensity_model_fn` via
#'   `fit_treatment_model()` (e.g. smoothing arguments for `mgcv::gam`).
#'
#' @return A `causatr_fit` object with `weights_obj = NULL`, a
#'   placeholder `Y ~ A` model in `$model`, and the treatment density
#'   engine stashed in `$details`.
#'
#' @noRd
fit_ipw <- function(
  data,
  outcome,
  treatment,
  confounders,
  confounders_tv,
  family,
  estimand,
  type,
  history,
  numerator,
  weights,
  model_fn,
  propensity_model_fn,
  propensity_family = NULL,
  call,
  ...
) {
  if (type == "longitudinal") {
    rlang::abort(
      "Longitudinal IPW is not supported. Use `estimator = 'gcomp'` with `type = 'longitudinal'` for time-varying treatments.",
      .call = FALSE
    )
  }

  is_multivariate <- length(treatment) > 1L

  # Parse effect-modification terms and reject bare treatment in
  # confounders (`~ L + A`). True EM terms (`A:sex`) are detected
  # and stored for downstream MSM expansion. The propensity formula
  # strips EM terms automatically via `build_ps_formula()`.
  em_info <- check_confounders_treatment(
    confounders,
    treatment,
    estimator = "ipw"
  )

  # Effect modification under multivariate IPW: per-component
  # propensity formulas strip ALL treatment-touching terms (including
  # `A_k:modifier` interactions for any k), and the per-intervention
  # MSM expands to `Y ~ 1 + modifier_main_effects` via the existing
  # `build_ipw_msm_formula()`. Treatment-treatment interactions like
  # `A1:A2` are handled implicitly via the per-intervention refit, so
  # they don't enter the MSM. The Phase 6 baseline-modifier
  # constraint (Robins 2000; modifier must be baseline, not
  # post-treatment) carries over and is doc-level.

  # EM terms are stored in `fit$details$em_info` for downstream use by
  # `compute_ipw_contrast_point()` (MSM expansion from `Y ~ 1` to
  # `Y ~ 1 + modifier`) and `variance_if_ipw()` (expanded MSM flows
  # through the stacked sandwich engine unchanged).

  # Fit rows: exclude missing outcomes for the MSM. `fit_treatment_model()`
  # applies its own propensity-side row mask (missing treatment /
  # confounders) inside; if that differs from the outcome-row mask we
  # abort below so downstream row-alignment invariants hold.
  fit_rows <- get_fit_rows(data, outcome)
  fit_data <- data[fit_rows]

  # Resolve the propensity fitter. For binary / continuous treatments
  # the default is the user's `model_fn` (typically stats::glm), so a
  # single `causat(..., model_fn = mgcv::gam)` call flexibly models
  # both the outcome and the propensity. For **categorical** treatments
  # the fitter must be a multinomial model (`nnet::multinom` by
  # default) because `stats::glm` cannot fit a multinomial response.
  # For **negative binomial** treatments the fitter must be
  # `MASS::glm.nb` which estimates theta by MLE (no `family` arg).
  # When the user passes an explicit `propensity_model_fn`, it takes
  # precedence regardless of treatment family.
  #
  # Multivariate dispatch: per-component family detection happens
  # inside `fit_treatment_models()`. We resolve a single propensity
  # fitter here that's shared across all components; opt-in count or
  # categorical components in a multivariate call are deferred (the
  # auto-detect rejects categorical components upstream in
  # `fit_treatment_models()`).
  trt_family <- if (is_multivariate) {
    detect_treatment_family(data[[treatment[1]]])
  } else {
    detect_treatment_family(data[[treatment]])
  }
  prop_model_fn <- if (!is.null(propensity_model_fn)) {
    propensity_model_fn
  } else if (!is_multivariate && trt_family == "categorical") {
    nnet::multinom
  } else if (!is_multivariate && identical(propensity_family, "negbin")) {
    check_pkg("MASS")
    MASS::glm.nb
  } else {
    model_fn
  }

  # `propensity_family` for multivariate accepts NULL, length 1
  # (broadcast across components), or length K (per-component opt-in
  # to count families). `fit_treatment_models()` validates the shape.

  # Capture the user's `...` once. Stashed in `fit$details$dots` so
  # the bootstrap refit replays the exact same propensity fitting
  # call on resampled data (see `refit_ipw()` in
  # `R/variance_bootstrap.R`).
  dots <- list(...)

  # Fit the conditional treatment density. The returned
  # `causatr_treatment_model` (univariate) or `causatr_treatment_models`
  # list (multivariate) carries the fitted model(s), propensity design
  # matrix(es), `alpha_hat`, family tag(s), and the `fit_rows` mask --
  # everything `make_weight_fn()` / `make_weight_fn_mv()` needs to
  # build a weight closure for any intervention downstream.
  if (is_multivariate) {
    # For multivariate, defer per-component fitter selection to
    # `fit_treatment_models()`. The `model_fn` is the default base
    # (typically stats::glm) that's used for binary / continuous /
    # poisson components; categorical components auto-pick
    # `nnet::multinom`, negbin auto-picks `MASS::glm.nb`. The user
    # can override globally via `propensity_model_fn`.
    tm_args <- list(
      data = fit_data,
      treatment = treatment,
      confounders = confounders,
      model_fn = model_fn,
      propensity_model_fn = propensity_model_fn,
      propensity_family = propensity_family
    )
    if (!is.null(weights)) {
      tm_args$weights <- weights[fit_rows]
    }
    treatment_models <- do.call(
      fit_treatment_models,
      c(tm_args, dots)
    )
    # Row-alignment guard for multivariate: every component model is
    # fit on the same `fit_data`, and `fit_treatment_model()` only
    # drops rows where the treatment / confounders are NA. Verify the
    # K masks all match the outcome mask; mismatch means a confounder
    # NA the outcome mask does not exclude.
    n_fit_outcome <- sum(fit_rows)
    for (k in seq_along(treatment_models)) {
      n_fit_k <- sum(treatment_models[[k]]$fit_rows)
      if (n_fit_k != n_fit_outcome) {
        rlang::abort(
          paste0(
            "Treatment density model for component '",
            treatment[k],
            "' kept ",
            n_fit_k,
            " rows but the outcome-non-missing mask has ",
            n_fit_outcome,
            " rows. Drop or impute the offending rows before calling `causat()`."
          )
        )
      }
    }
    tm <- NULL
    propensity_model <- NULL
  } else {
    tm_args <- list(
      data = fit_data,
      treatment = treatment,
      confounders = confounders,
      model_fn = prop_model_fn,
      propensity_family = propensity_family
    )
    if (!is.null(weights)) {
      tm_args$weights <- weights[fit_rows]
    }
    # Forward `...` to the propensity fitter. `do.call()` with a
    # constructed arg list + `dots` is the minimal way to get "pass
    # every extra named argument through" without hand-listing every
    # possible option the user might forward to `mgcv::gam()` /
    # `stats::glm()`.
    tm <- do.call(
      fit_treatment_model,
      c(tm_args, dots)
    )

    # Row-alignment guard: `fit_treatment_model()` drops rows with
    # missing treatment or confounders, but `get_fit_rows()` only drops
    # on the outcome. If those masks diverge, the density-ratio weight
    # vector has a different length than the MSM refit's data, and the
    # variance engine's n_total / n_fit invariants break. Abort here
    # with a clear pointer at the missingness source.
    n_fit_outcome <- sum(fit_rows)
    n_fit_ps <- sum(tm$fit_rows)
    if (n_fit_ps != n_fit_outcome) {
      rlang::abort(
        paste0(
          "Treatment density model kept ",
          n_fit_ps,
          " rows but the outcome-non-missing mask has ",
          n_fit_outcome,
          " rows. A confounder column has missing values the outcome ",
          "mask does not exclude. Drop or impute those rows before ",
          "calling `causat()` so the propensity and MSM fits agree."
        )
      )
    }
    treatment_models <- structure(
      list(tm),
      class = c("causatr_treatment_models", "list"),
      names = treatment
    )
    propensity_model <- tm$model
  }

  # Placeholder outcome model. Kept in `$model` so `print()` /
  # `summary()` can show something sensible even before a contrast is
  # requested, and so downstream `fit$model$family` lookups (used in
  # places like `resolve_family()`) find a valid family. The
  # estimation path does NOT consume this model -- every intervention
  # in `compute_contrast()` refits its own weighted MSM.
  fam_obj <- resolve_family(family)
  placeholder_formula <- stats::reformulate(treatment, response = outcome)
  placeholder_args <- list(
    formula = placeholder_formula,
    data = fit_data,
    family = fam_obj
  )
  if (!is.null(weights)) {
    placeholder_args$weights <- weights[fit_rows]
  }
  placeholder_model <- do.call(model_fn, placeholder_args)

  new_causatr_fit(
    model = placeholder_model,
    data = data,
    treatment = treatment,
    outcome = outcome,
    confounders = confounders,
    confounders_tv = confounders_tv,
    family = family,
    estimator = "ipw",
    type = "point",
    estimand = estimand,
    id = NULL,
    time = NULL,
    censoring = NULL,
    history = history,
    numerator = numerator,
    weights_obj = NULL,
    match_obj = NULL,
    call = call,
    details = list(
      fit_rows = fit_rows,
      n_fit = n_fit_outcome,
      n_total = nrow(data),
      weights = weights,
      dots = dots,
      treatment_model = tm,
      treatment_models = treatment_models,
      propensity_model = propensity_model,
      propensity_model_fn = prop_model_fn,
      propensity_family = propensity_family,
      model_fn = model_fn,
      em_info = em_info,
      is_multivariate = is_multivariate
    )
  )
}


#' Compute IPW per-intervention bundles for point-treatment data
#'
#' @description
#' Per-intervention builder behind `compute_contrast()`'s IPW branch.
#' For each intervention:
#'
#' 1. Build the density-ratio weight vector via
#'    `compute_density_ratio_weights(tm, data, intervention)`.
#' 2. Multiply in external (survey / IPCW) weights when present.
#' 3. Refit a weighted MSM on `Y ~ 1`, which under a weighted Hajek
#'    estimator recovers the marginal counterfactual mean
#'    \eqn{E[Y^d]}.
#' 4. Compute \eqn{\hat\mu_a} by predict-then-average over the target
#'    population.
#'
#' Returns the per-intervention bundles plus the scalars
#' `compute_contrast()` needs to dispatch variance: `fit_idx`,
#' `n_total`, and the marginal-mean vector. Variance itself is
#' computed by `compute_contrast()` via `variance_if()` (sandwich) or
#' `variance_bootstrap()` (bootstrap), uniformly with the other
#' estimators.
#'
#' @param fit A `causatr_fit` from `fit_ipw()`.
#' @param interventions Named list of `causatr_intervention` objects
#'   (or `NULL` entries for natural course).
#' @param target_idx Logical vector (length `nrow(fit$data)`) flagging
#'   target-population rows.
#'
#' @return A list with components `bundles` (named list of per-
#'   intervention results), `mu_hat` (named numeric vector),
#'   `fit_idx` (integer MSM fit-row indices), and `n_total`
#'   (`nrow(fit$data)`).
#'
#' @noRd
compute_ipw_contrast_point <- function(
  fit,
  interventions,
  target_idx
) {
  data <- fit$data
  treatment <- fit$treatment
  outcome <- fit$outcome
  tms <- fit$details$treatment_models
  is_mv <- isTRUE(fit$details$is_multivariate)
  tm <- fit$details$treatment_model
  if (!is_mv && is.null(tm)) {
    rlang::abort(
      "Internal error: IPW fit is missing `details$treatment_model`."
    )
  }
  if (is_mv && is.null(tms)) {
    rlang::abort(
      "Internal error: multivariate IPW fit is missing `details$treatment_models`."
    )
  }
  fit_rows <- fit$details$fit_rows
  fit_data <- data[fit_rows]
  ext_w <- fit$details$weights
  ext_w_fit <- if (is.null(ext_w)) NULL else ext_w[fit_rows]
  fam_obj <- resolve_family(fit$family)
  model_fn <- fit$details$model_fn
  em_info <- fit$details$em_info
  # Fit-time estimand threads into the density-ratio weight branch via
  # the HT Bayes-numerator `f*_i`. For ATE this is 1 (unchanged from
  # Phase 4's original flow); for ATT / ATC on static binary it
  # becomes `p(L_i)` / `1 - p(L_i)`, implementing the per-arm weight
  # tables derived in `compute_density_ratio_weights()`'s roxygen.
  # `check_estimand_intervention_compat()` has already rejected
  # ATT/ATC on anything but static binary by the time we get here,
  # so the other branches (IPSI, shift/scale) see `estimand = "ATE"`
  # in practice and are unaffected.
  estimand <- fit$estimand
  n_total <- nrow(data)
  fit_idx <- which(fit_rows)
  n_fit <- length(fit_idx)
  int_names <- names(interventions)
  k <- length(int_names)

  # Pre-compute the intervention-specific MSM + weight vector +
  # marginal mean. Stored in a per-intervention list used both by
  # the sandwich and bootstrap paths.
  bundles <- lapply(int_names, function(nm) {
    iv <- interventions[[nm]]

    # Density-ratio weights (length n_fit). Natural-course returns
    # all ones; `compute_density_ratio_weights()` dispatches on
    # intervention type and treatment family. The `estimand` argument
    # selects the Bayes numerator `f*` so the same HT branch produces
    # ATE / ATT / ATC weights uniformly.
    # `tm$fit_rows` is relative to `fit_data` (the outcome-clean
    # subset), so the density-ratio computation must receive
    # `fit_data`, not the full `data`. Without this, outcome NAs
    # create a length mismatch between `tm$fit_rows` and `nrow(data)`.
    #
    # Multivariate dispatch: `compute_density_ratio_weights_mv()`
    # builds the joint weight as a product of per-component density
    # ratios under the chain-rule factorisation. Per-component
    # interventions live in the named sub-list `iv` (already
    # name-validated by `apply_intervention()` upstream).
    w_iv <- if (is_mv) {
      tms_local <- tms
      for (k in seq_along(tms_local)) {
        tms_local[[k]]$fit_rows <- rep(TRUE, nrow(fit_data))
      }
      class(tms_local) <- c("causatr_treatment_models", "list")
      compute_density_ratio_weights_mv(
        tms_local,
        fit_data,
        iv,
        estimand = estimand
      )
    } else {
      compute_density_ratio_weights(tm, fit_data, iv, estimand = estimand)
    }

    # Compose with external weights. The density-ratio weights enter
    # multiplicatively because survey / IPCW weights represent an
    # independent reweighting of the target population (sampling
    # design) from the density ratio (causal reweighting).
    w_final <- if (is.null(ext_w_fit)) w_iv else w_iv * ext_w_fit

    # The unified density-ratio engine always produces an
    # intervention-specific weight vector: HT indicators for static
    # binary, pushforward ratios for shift / scale_by, Kennedy's
    # closed form for IPSI. Under every one of those the weighted
    # MSM for `E[Y^d]` collapses to an intercept-only Hajek mean of
    # Y. A saturated `Y ~ A` would be *rank-deficient* under HT
    # weights on a binary treatment -- every surviving row has the
    # same value of A -- so we uniformly refit `Y ~ 1` and read off
    # the (sole) coefficient as the counterfactual marginal mean.
    #
    # When effect-modification terms are present (e.g. `A:sex` in
    # confounders), the MSM expands to `Y ~ 1 + sex` so that
    # `predict()` returns modifier-stratum-specific counterfactual
    # means. The treatment is still absorbed by the density-ratio
    # weights -- only modifier main effects enter the MSM.
    msm_formula <- build_ipw_msm_formula(outcome, em_info)

    msm_args <- list(
      formula = msm_formula,
      data = fit_data,
      family = fam_obj,
      weights = w_final
    )
    msm_model <- do.call(model_fn, msm_args)

    # Counterfactual marginal mean: predict on the intervened data
    # restricted to the target population, then (optionally-weighted)
    # average.
    #
    # Predictions and averaging are restricted to `fit_data` (rows with
    # non-missing outcome) so the marginal-mean population matches the
    # variance engine's population in `variance_if_ipw()`, which operates
    # on `fit_rows` only. Without this restriction, NA-outcome rows enter
    # the marginal mean but not the variance engine's Ch1, causing a
    # centering mismatch that biases the sandwich SE when the MSM is
    # non-intercept-only (i.e. with EM terms like `Y ~ 1 + modifier`).
    # For `Y ~ 1` (no EM) this is a no-op because all predictions are
    # identical. Fifth-round critical review Issue #B2.
    #
    # IPSI does not materialize a counterfactual treatment value -- the
    # intervention acts on the propensity, not on A itself. The MSM
    # prediction depends only on modifier columns (which are unchanged
    # by IPSI), so we skip `apply_intervention()` and use the original
    # fit_data directly. For all other interventions, the treatment
    # column is modified to its counterfactual value.
    iv_type <- if (inherits(iv, "causatr_intervention")) iv$type else NULL
    data_a_fit <- if (identical(iv_type, "ipsi")) {
      fit_data
    } else {
      apply_intervention(fit_data, treatment, iv)
    }
    preds_fit <- stats::predict(
      msm_model,
      newdata = data_a_fit,
      type = "response"
    )
    target_fit <- target_idx[fit_rows]
    valid_fit <- target_fit & !is.na(preds_fit)
    w_target <- if (!is.null(ext_w_fit)) ext_w_fit[valid_fit] else NULL
    mu_hat_iv <- maybe_weighted_mean(preds_fit[valid_fit], w_target)

    list(
      intervention = iv,
      msm_model = msm_model,
      weights_final = w_final,
      preds = preds_fit,
      valid = valid_fit,
      mu_hat = mu_hat_iv
    )
  })
  names(bundles) <- int_names

  mu_hat <- vapply(bundles, function(b) b$mu_hat, numeric(1))
  names(mu_hat) <- int_names

  list(
    bundles = bundles,
    mu_hat = mu_hat,
    fit_idx = fit_idx,
    n_total = n_total
  )
}
