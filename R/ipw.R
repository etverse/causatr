#' Fit an IPW model for causal estimation (point treatment)
#'
#' @description
#' Implements inverse probability weighting (Hernán & Robins Ch. 12) via a
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
#'    - **Categorical** (factor / character): aborts — users with
#'      categorical treatments should use `estimator = "gcomp"`.
#' 2. Stash the treatment model in `fit$details$treatment_model` so
#'    `contrast()` can build intervention-specific density-ratio
#'    weight vectors via `compute_density_ratio_weights()` /
#'    `make_weight_fn()` and refit a weighted MSM per intervention.
#' 3. Store a cheap placeholder outcome model (`Y ~ A`, unweighted)
#'    in `fit$model` for display / compatibility with `print()` /
#'    `summary()`. This model is **not** used for estimation — every
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
  call,
  ...
) {
  if (type == "longitudinal") {
    rlang::abort(
      "Longitudinal IPW is not supported. Use `estimator = 'gcomp'` with `type = 'longitudinal'` for time-varying treatments.",
      .call = FALSE
    )
  }

  # Multivariate treatment is not supported by the self-contained
  # density-ratio engine: `fit_treatment_model()` fits a single
  # univariate density and downstream `make_weight_fn()` closures
  # assume a length-n treatment vector, not a matrix.
  if (length(treatment) > 1L) {
    rlang::abort(
      c(
        "Multivariate treatments are not supported under `estimator = 'ipw'`.",
        i = "Use `estimator = 'gcomp'` for joint (multivariate) treatment interventions."
      )
    )
  }

  # Reject A-touching terms in `confounders` before building the PS
  # formula. The MSM refits inside `contrast()` use either `Y ~ A`
  # (saturated, static binary) or `Y ~ 1` (non-static and non-binary
  # static), neither of which has room for an `A:modifier` slot.
  # Silently dropping the term would return a pooled ATE under an
  # effect-modification request. See `check_confounders_no_treatment()`
  # in `R/utils.R`.
  check_confounders_no_treatment(confounders, treatment, estimator = "ipw")

  # Fit rows: exclude missing outcomes for the MSM. `fit_treatment_model()`
  # applies its own propensity-side row mask (missing treatment /
  # confounders) inside; if that differs from the outcome-row mask we
  # abort below so downstream row-alignment invariants hold.
  fit_rows <- get_fit_rows(data, outcome)
  fit_data <- data[fit_rows]

  # Resolve the propensity fitter. Default to the user's `model_fn`
  # so a single `causat(..., model_fn = mgcv::gam)` call flexibly
  # models both the outcome (in gcomp) and the propensity (in ipw)
  # without a second argument. A caller who wants asymmetric shapes —
  # e.g. plain `glm` outcome but `mgcv::gam` propensity — supplies
  # `propensity_model_fn` explicitly.
  prop_model_fn <- if (is.null(propensity_model_fn)) {
    model_fn
  } else {
    propensity_model_fn
  }

  # Capture the user's `...` once. Stashed in `fit$details$dots` so
  # the bootstrap refit replays the exact same propensity fitting
  # call on resampled data (see `refit_ipw()` in
  # `R/variance_bootstrap.R`).
  dots <- list(...)

  # Fit the conditional treatment density. The returned
  # `causatr_treatment_model` carries the fitted model, the
  # propensity design matrix, `alpha_hat`, the family tag, and the
  # `fit_rows` mask — everything `make_weight_fn()` needs to build a
  # `w(alpha)` closure for any intervention downstream.
  tm_args <- list(
    data = fit_data,
    treatment = treatment,
    confounders = confounders,
    model_fn = prop_model_fn
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

  # Placeholder outcome model. Kept in `$model` so `print()` /
  # `summary()` can show something sensible even before a contrast is
  # requested, and so downstream `fit$model$family` lookups (used in
  # places like `resolve_family()`) find a valid family. The
  # estimation path does NOT consume this model — every intervention
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
      propensity_model = tm$model,
      propensity_model_fn = prop_model_fn,
      model_fn = model_fn
    )
  )
}


#' Compute IPW marginal-mean contrasts for point-treatment data
#'
#' @description
#' Per-intervention workhorse behind `compute_contrast()`'s IPW branch.
#' For each intervention:
#'
#' 1. Build the density-ratio weight vector via
#'    `compute_density_ratio_weights(tm, data, intervention)`.
#' 2. Multiply in external (survey / IPCW) weights when present.
#' 3. Refit a weighted MSM. Static on binary uses the saturated
#'    `Y ~ A`; non-static and static-on-non-binary use an intercept-
#'    only `Y ~ 1` (which under a weighted Hájek estimator recovers
#'    the marginal counterfactual mean \eqn{E[Y^d]}).
#' 4. Compute \eqn{\hat\mu_a} by predict-then-average over the target
#'    population: for `Y ~ A` this averages `beta_0 + beta_1 * a_i`
#'    over target rows; for `Y ~ 1` it is `beta_0`.
#'
#' The result carries enough state for `variance_if_ipw()` and the
#' bootstrap path to build their own variance estimates without
#' re-running the weight / MSM machinery.
#'
#' @param fit A `causatr_fit` from `fit_ipw()`.
#' @param interventions Named list of `causatr_intervention` objects
#'   (or `NULL` entries for natural course).
#' @param target_idx Logical vector (length `nrow(fit$data)`) flagging
#'   target-population rows.
#' @param ci_method `"sandwich"` or `"bootstrap"`.
#' @param est,subset,n_boot,parallel,ncpus,subset_env Passed through
#'   to `variance_bootstrap()` for the bootstrap path.
#'
#' @return A list with components `mu_hat` (named numeric vector),
#'   `vcov` (k x k matrix), `boot_t` / `boot_info` (non-NULL only
#'   under bootstrap).
#'
#' @noRd
compute_ipw_contrast_point <- function(
  fit,
  interventions,
  target_idx,
  ci_method,
  est,
  subset,
  n_boot,
  parallel,
  ncpus,
  subset_env
) {
  data <- fit$data
  treatment <- fit$treatment
  outcome <- fit$outcome
  tm <- fit$details$treatment_model
  if (is.null(tm)) {
    rlang::abort(
      "Internal error: IPW fit is missing `details$treatment_model`."
    )
  }
  fit_rows <- fit$details$fit_rows
  fit_data <- data[fit_rows]
  ext_w <- fit$details$weights
  ext_w_fit <- if (is.null(ext_w)) NULL else ext_w[fit_rows]
  fam_obj <- resolve_family(fit$family)
  model_fn <- fit$details$model_fn
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
    # intervention type and treatment family.
    w_iv <- compute_density_ratio_weights(tm, data, iv)

    # Compose with external weights. The density-ratio weights enter
    # multiplicatively because survey / IPCW weights represent an
    # independent reweighting of the target population (sampling
    # design) from the density ratio (causal reweighting).
    w_final <- if (is.null(ext_w_fit)) w_iv else w_iv * ext_w_fit

    # The unified density-ratio engine always produces an
    # intervention-specific weight vector: HT indicators for static
    # binary, pushforward ratios for shift / scale_by, Kennedy's
    # closed form for IPSI. Under every one of those the weighted
    # MSM for `E[Y^d]` collapses to an intercept-only Hájek mean of
    # Y. A saturated `Y ~ A` would be *rank-deficient* under HT
    # weights on a binary treatment — every surviving row has the
    # same value of A — so we uniformly refit `Y ~ 1` and read off
    # the (sole) coefficient as the counterfactual marginal mean.
    msm_formula <- stats::as.formula(paste0(outcome, " ~ 1"))

    msm_args <- list(
      formula = msm_formula,
      data = fit_data,
      family = fam_obj,
      weights = w_final
    )
    msm_model <- do.call(model_fn, msm_args)

    # Counterfactual marginal mean: predict on the intervened data
    # restricted to the target population, then (optionally-weighted)
    # average. This mirrors the gcomp / matching predict-then-average
    # path — under a saturated `Y ~ A` the result collapses to
    # `beta_0 + beta_1 * target_value`, and under `Y ~ 1` it
    # collapses to `beta_0`, but running the generic path keeps the
    # code uniform with the rest of `compute_contrast()`.
    data_a <- apply_intervention(data, treatment, iv)
    preds <- stats::predict(msm_model, newdata = data_a, type = "response")
    valid <- target_idx & !is.na(preds)
    w_target <- if (!is.null(ext_w)) ext_w[valid] else NULL
    mu_hat_iv <- maybe_weighted_mean(preds[valid], w_target)

    list(
      intervention = iv,
      msm_model = msm_model,
      weights_final = w_final,
      preds = preds,
      valid = valid,
      mu_hat = mu_hat_iv
    )
  })
  names(bundles) <- int_names

  mu_hat <- vapply(bundles, function(b) b$mu_hat, numeric(1))
  names(mu_hat) <- int_names

  boot_t <- NULL
  boot_info <- NULL
  if (ci_method == "sandwich") {
    vcov_mat <- variance_if_ipw_self_contained(
      fit,
      bundles,
      target_idx,
      mu_hat,
      fit_idx,
      n_total
    )
  } else {
    boot_res <- variance_bootstrap(
      fit,
      interventions,
      n_boot,
      target_idx,
      est,
      subset,
      parallel,
      ncpus,
      subset_env = subset_env
    )
    vcov_mat <- boot_res$vcov
    boot_t <- boot_res$boot_t
    boot_info <- boot_res$boot_info
  }

  list(
    mu_hat = mu_hat,
    vcov = vcov_mat,
    boot_t = boot_t,
    boot_info = boot_info
  )
}
