#' Pool causal estimates across multiply-imputed datasets
#'
#' @description
#' Fits a causal model with [causat()] and computes a causal [contrast()] on
#' every completed dataset stored in a `mice` `mids` object, then pools the
#' per-imputation estimates into a single `causatr_result`. This is the
#' analysis step of a multiple-imputation (MI) workflow: the user imputes
#' missing covariates and/or treatment upstream with `mice::mice()`, and
#' `causat_mice()` propagates the imputation uncertainty into the causal
#' estimate and its standard error.
#'
#' Multiple imputation is the right tool for **missing covariates (L)** or
#' **missing treatment (A)** under a missing-at-random mechanism. Missing
#' **outcomes (Y)** are handled by inverse-probability-of-censoring weighting
#' (`ipcw = TRUE`) or complete-case analysis, not by imputing Y; however Y
#' *should* be a predictor in the upstream imputation model.
#'
#' @param imp A `mids` object returned by `mice::mice()`.
#' @param outcome Character scalar naming the outcome column. Passed to
#'   [causat()].
#' @param treatment Character scalar (or vector for multivariate treatment)
#'   naming the treatment column(s). Passed to [causat()].
#' @param confounders A one-sided formula of baseline confounders, or `NULL`
#'   when per-component formulas are supplied through `...`. Passed to
#'   [causat()].
#' @param interventions A named list of intervention objects (e.g.
#'   `list(a1 = static(1), a0 = static(0))`). Passed to [contrast()]. Leave
#'   `NULL` for `estimator = "snm"`, whose estimand is the blip parameter
#'   itself and which rejects an `interventions` argument.
#' @param estimator Character causal estimator: `"gcomp"` (default), `"ipw"`,
#'   `"aipw"`, `"matching"`, or `"snm"`. Passed to [causat()].
#' @param family Outcome family (character or family object). Passed to
#'   [causat()].
#' @param estimand Character estimand (`"ATE"`, `"ATT"`, `"ATC"`). Passed to
#'   [causat()].
#' @param type Character contrast scale: `"difference"` (default), `"ratio"`,
#'   or `"or"`. Passed to [contrast()].
#' @param ci_method Character within-imputation variance method, `"sandwich"`
#'   (default) or `"bootstrap"`, used for each per-imputation [contrast()]
#'   call. The *pooled* variance is governed by `pool_method`, not this
#'   argument.
#' @param conf_level Numeric confidence level for the pooled intervals.
#'   Default `0.95`.
#' @param by Optional one-sided formula or character naming a baseline
#'   stratifier. Pooling is applied per `by`-stratum row independently. Passed
#'   to [contrast()].
#' @param pool_method Character pooling strategy. `"rubin"` (default) applies
#'   Rubin's rules to the per-imputation sandwich variances. `"boot_mi"` uses
#'   von Hippel's bootstrap-then-impute two-stage variance, valid under
#'   uncongeniality. See Details.
#' @param B Integer number of bootstrap resamples for `pool_method =
#'   "boot_mi"`. Default `200`. Ignored for `"rubin"`.
#' @param M Integer number of imputations per bootstrap resample for
#'   `pool_method = "boot_mi"`. Default `2` (von Hippel's efficient variant).
#'   Ignored for `"rubin"`.
#' @param parallel Character parallel backend forwarded to the Boot MI engine:
#'   `"no"` (default) or `"future"` (uses `future.apply::future_lapply()`).
#' @param seed Optional integer seed. For `pool_method = "boot_mi"` it seeds
#'   the bootstrap-and-impute loop reproducibly.
#' @param ... Additional arguments forwarded to [causat()] (e.g. `id`,
#'   `time`, `confounders_tv`, `censoring`, `ipcw`, `confounders_outcome`,
#'   `propensity_model_fn`, `model_fn`).
#'
#' @returns A `causatr_result` with pooled estimates, standard errors, and
#'   confidence intervals. `ci_method` is `"rubin"` or `"boot_mi"`. The
#'   per-row pooling diagnostics are attached as the `"mi_details"` attribute.
#'
#' @details
#' ## Rubin's rules (`pool_method = "rubin"`)
#' Let \eqn{\hat{Q}_i} and \eqn{U_i} be the estimate and variance from
#' imputation \eqn{i}. The pooled estimate is \eqn{\bar{Q} = m^{-1}\sum_i
#' \hat{Q}_i} and the total variance is \eqn{T = \bar{U} + (1 + 1/m) B} with
#' within variance \eqn{\bar{U} = m^{-1}\sum_i U_i} and between variance
#' \eqn{B = (m-1)^{-1}\sum_i (\hat{Q}_i - \bar{Q})^2}. Confidence intervals use
#' Barnard-Rubin degrees of freedom.
#'
#' ## Congeniality
#' Causal estimands are typically *uncongenial* with the `mice` imputation
#' model (the estimand is a functional of the outcome/treatment model under
#' intervention, not a parameter of the imputation model). Under
#' uncongeniality Rubin's rules give conservative (wider) intervals, which is
#' safe but may sacrifice power. `pool_method = "boot_mi"` restores nominal
#' coverage. Always include the outcome, treatment, all confounders, and any
#' effect modifiers as predictors in the upstream `mice()` call;
#' `causat_mice()` warns when an analysis variable is absent.
#'
#' ## What this does not do
#' It does not perform the imputation (call `mice::mice()` first), impute the
#' outcome, handle MNAR mechanisms, or pool omnibus tests across contrasts.
#'
#' @references
#' Rubin DB (1987). *Multiple Imputation for Nonresponse in Surveys*. Wiley.
#'
#' van Buuren S, Groothuis-Oudshoorn K (2011). mice: Multivariate Imputation
#' by Chained Equations in R. *Journal of Statistical Software* 45(3):1-67.
#'
#' von Hippel PT (2020). How many imputations do you need? *Sociological
#' Methods & Research* 49(3):699-718.
#'
#' @examples
#' if (requireNamespace("mice", quietly = TRUE)) {
#'   set.seed(1)
#'   n <- 400
#'   L <- rnorm(n)
#'   A <- rbinom(n, 1, plogis(0.5 * L))
#'   Y <- 2 + 3 * A + 1.5 * L + rnorm(n)
#'   # L missing-at-random on the (observed) treatment.
#'   L[rbinom(n, 1, plogis(-1 + 0.8 * A)) == 0] <- NA
#'   dat <- data.frame(Y = Y, A = A, L = L)
#'
#'   imp <- mice::mice(dat, m = 5, printFlag = FALSE)
#'   res <- causat_mice(
#'     imp,
#'     outcome = "Y",
#'     treatment = "A",
#'     confounders = ~L,
#'     interventions = list(a1 = static(1), a0 = static(0)),
#'     estimator = "gcomp"
#'   )
#'   summary(res)
#' }
#'
#' @seealso [causat()], [contrast()]
#' @export
causat_mice <- function(
  imp,
  outcome,
  treatment,
  confounders = NULL,
  interventions = NULL,
  estimator = "gcomp",
  family = "gaussian",
  estimand = "ATE",
  type = "difference",
  ci_method = "sandwich",
  conf_level = 0.95,
  by = NULL,
  pool_method = c("rubin", "boot_mi"),
  B = 200L,
  M = 2L,
  parallel = c("no", "future"),
  seed = NULL,
  ...
) {
  check_pkg("mice")

  if (!inherits(imp, "mids")) {
    rlang::abort(
      c(
        "`imp` must be a `mids` object returned by `mice::mice()`.",
        x = paste0(
          "Got an object of class: ",
          paste(class(imp), collapse = ", "),
          "."
        ),
        i = "Run `mice::mice(data, ...)` first, then pass the result here."
      ),
      class = "causatr_mi_bad_input"
    )
  }

  pool_method <- rlang::arg_match(pool_method)
  parallel <- rlang::arg_match(parallel)

  # Assemble the two argument bundles once: everything that defines the
  # per-imputation analysis (causat() + contrast()) is shared verbatim by
  # both the Rubin and Boot MI engines.
  fit_args <- c(
    list(
      outcome = outcome,
      treatment = treatment,
      confounders = confounders,
      estimator = estimator,
      family = family,
      estimand = estimand
    ),
    list(...)
  )
  contrast_args <- list(
    type = type,
    ci_method = ci_method,
    conf_level = conf_level,
    by = by
  )
  # SNM rejects `interventions =`; for every other estimator it is required.
  if (!is.null(interventions)) {
    contrast_args$interventions <- interventions
  }

  # Warn when a variable named in the analysis is absent from the imputation
  # model. mice can only carry uncertainty for variables it knows about;
  # omitting an analysis variable (especially an effect modifier or the
  # outcome) biases the imputed values and worsens uncongeniality.
  mi_check_predictors(imp, fit_args, by)

  if (pool_method == "boot_mi") {
    return(pool_boot_mi(
      imp = imp,
      fit_args = fit_args,
      contrast_args = contrast_args,
      conf_level = conf_level,
      B = B,
      M = M,
      parallel = parallel,
      seed = seed,
      call = match.call()
    ))
  }

  collected <- mi_collect_imputations(
    imp = imp,
    fit_args = fit_args,
    contrast_args = contrast_args,
    call = match.call()
  )
  pool_rubin(collected, conf_level = conf_level)
}

#' Run causat() + contrast() on a single completed dataset
#'
#' @param data A complete `data.frame` (one imputation, or a bootstrap
#'   resample thereof).
#' @param fit_args Named list of arguments forwarded to [causat()] (without
#'   `data`).
#' @param contrast_args Named list of arguments forwarded to [contrast()]
#'   (without `fit`).
#' @returns A `causatr_result` for this single dataset.
#' @noRd
mi_fit_one <- function(data, fit_args, contrast_args) {
  fit <- do.call(causat, c(list(data = data), fit_args))
  do.call(contrast, c(list(fit = fit), contrast_args))
}

#' Extract the per-imputation estimate / SE rows from one result
#'
#' @param res A `causatr_result`.
#' @returns A list with the means-table label column name (`term_col`), the
#'   label vector (`est_labels`), estimate / SE vectors, optional `by`
#'   vector, and the analogous contrast vectors.
#' @noRd
mi_extract <- function(res) {
  # Means table key column: SNM blip parameters use "parameter"; every
  # other estimator uses "intervention".
  term_col <- if ("parameter" %in% names(res$estimates)) {
    "parameter"
  } else {
    "intervention"
  }
  est_by <- if ("by" %in% names(res$estimates)) res$estimates$by else NULL
  con_by <- if (nrow(res$contrasts) > 0L && "by" %in% names(res$contrasts)) {
    res$contrasts$by
  } else {
    NULL
  }
  list(
    term_col = term_col,
    est_labels = as.character(res$estimates[[term_col]]),
    est = res$estimates$estimate,
    est_se = res$estimates$se,
    est_by = est_by,
    con_labels = if (nrow(res$contrasts) > 0L) {
      as.character(res$contrasts$comparison)
    } else {
      character(0)
    },
    con = if (nrow(res$contrasts) > 0L) res$contrasts$estimate else numeric(0),
    con_se = if (nrow(res$contrasts) > 0L) res$contrasts$se else numeric(0),
    con_by = con_by,
    n = res$n,
    type = res$type,
    estimand = res$estimand,
    reference = res$reference,
    interventions = res$interventions,
    estimator = res$estimator,
    family = res$family,
    fit_type = res$fit_type
  )
}

#' Loop over imputations, collecting stacked estimate / SE matrices
#'
#' @param imp A `mids` object.
#' @param fit_args,contrast_args Argument bundles for `mi_fit_one()`.
#' @param call The originating `causat_mice()` call, stored on the result.
#' @returns A `collected` list consumed by [pool_rubin()]: stacked `m x p`
#'   estimate / SE matrices for means and contrasts, the label vectors, and
#'   the metadata needed to rebuild a `causatr_result`.
#' @noRd
mi_collect_imputations <- function(imp, fit_args, contrast_args, call) {
  m <- imp$m
  per_imp <- vector("list", m)
  failures <- character(0)

  for (i in seq_len(m)) {
    data_i <- mice::complete(imp, i)
    # A single imputation that fails to fit (e.g. perfect separation in a
    # propensity model on one completed dataset) must not sink the whole
    # pool; record it and drop it, then decide downstream whether enough
    # imputations survived.
    res_i <- tryCatch(
      mi_fit_one(data_i, fit_args, contrast_args),
      error = function(e) {
        failures[[length(failures) + 1L]] <<- conditionMessage(e)
        NULL
      }
    )
    if (!is.null(res_i)) {
      per_imp[[i]] <- mi_extract(res_i)
    }
  }

  per_imp <- per_imp[!vapply(per_imp, is.null, logical(1))]
  mi_assert_enough(length(per_imp), m, failures)
  if (length(per_imp) == 1L) {
    rlang::warn(
      c(
        "Only one imputation produced a usable fit; pooling is degenerate.",
        i = paste0(
          "Between-imputation variance is zero and the result equals a ",
          "single complete-data analysis."
        )
      ),
      class = "causatr_mi_degenerate"
    )
  }

  mi_stack(per_imp, call)
}

#' Assert enough imputations survived fitting
#'
#' @param n_ok Integer count of successful imputation fits.
#' @param m Integer total number of imputations.
#' @param failures Character vector of error messages from failed fits.
#' @returns `NULL` invisibly; aborts when no imputation produced a fit.
#' @noRd
mi_assert_enough <- function(n_ok, m, failures) {
  if (n_ok == 0L) {
    rlang::abort(
      c(
        "Every imputation failed to fit; cannot pool.",
        x = paste0("All ", m, " imputation fits errored."),
        i = if (length(failures) > 0L) {
          paste0("First error: ", failures[[1L]])
        } else {
          NULL
        }
      ),
      class = "causatr_mi_all_failed"
    )
  }
  if (n_ok < m) {
    rlang::warn(
      c(
        paste0(
          n_ok,
          " of ",
          m,
          " imputation fits succeeded; pooling the rest."
        ),
        i = "Failed imputations were dropped before applying Rubin's rules."
      ),
      class = "causatr_mi_partial"
    )
  }
  invisible(NULL)
}

#' Stack per-imputation extractions into pooling matrices
#'
#' @param per_imp List of `mi_extract()` outputs, one per surviving
#'   imputation.
#' @param call The originating `causat_mice()` call.
#' @returns A `collected` list for [pool_rubin()].
#' @noRd
mi_stack <- function(per_imp, call) {
  first <- per_imp[[1L]]
  est_labels <- first$est_labels
  con_labels <- first$con_labels

  # Row order (interventions x by-levels, or blip parameters) is fixed by
  # the spec and observed factor levels, so it is stable across imputations.
  # A mismatch means an imputation silently changed the estimand shape --
  # refuse to pool misaligned rows rather than average across them.
  for (k in seq_along(per_imp)) {
    if (!identical(per_imp[[k]]$est_labels, est_labels)) {
      rlang::abort(
        c(
          "Imputations produced inconsistent estimate labels; cannot pool.",
          x = "The intervention / parameter rows differ across imputations.",
          i = "This usually means a factor level is absent in one completed dataset."
        ),
        class = "causatr_mi_label_mismatch"
      )
    }
  }

  est_mat <- do.call(rbind, lapply(per_imp, `[[`, "est"))
  est_se_mat <- do.call(rbind, lapply(per_imp, `[[`, "est_se"))
  con_mat <- if (length(con_labels) > 0L) {
    do.call(rbind, lapply(per_imp, `[[`, "con"))
  } else {
    matrix(numeric(0), nrow = length(per_imp), ncol = 0L)
  }
  con_se_mat <- if (length(con_labels) > 0L) {
    do.call(rbind, lapply(per_imp, `[[`, "con_se"))
  } else {
    matrix(numeric(0), nrow = length(per_imp), ncol = 0L)
  }

  list(
    m = length(per_imp),
    n = first$n,
    est_labels = est_labels,
    est_mat = est_mat,
    est_se_mat = est_se_mat,
    est_by = first$est_by,
    con_labels = con_labels,
    con_mat = con_mat,
    con_se_mat = con_se_mat,
    con_by = first$con_by,
    meta = list(
      type = first$type,
      estimand = first$estimand,
      reference = first$reference,
      interventions = first$interventions,
      estimator = first$estimator,
      family = first$family,
      fit_type = first$fit_type,
      term_col = first$term_col,
      call = call
    )
  )
}

#' Warn when analysis variables are not used by the imputation model
#'
#' @param imp A `mids` object.
#' @param fit_args The [causat()] argument bundle (formulas are scanned for
#'   variable names).
#' @param by Optional `by` argument (formula / character) scanned for the
#'   stratifier name.
#' @returns `NULL` invisibly; emits a classed warning when an analysis
#'   variable is absent from the `mids` data, or is present but never used as
#'   a predictor when imputing other variables.
#' @noRd
mi_check_predictors <- function(imp, fit_args, by) {
  known <- colnames(imp$data)

  # Collect every variable the analysis references: outcome, treatment, and
  # all variables appearing in any confounder / by formula passed through.
  formula_args <- fit_args[grepl("^confounders", names(fit_args))]
  formula_vars <- unlist(lapply(formula_args, function(f) {
    if (inherits(f, "formula")) all.vars(f) else character(0)
  }))
  by_vars <- if (inherits(by, "formula")) {
    all.vars(by)
  } else if (is.character(by)) {
    by
  } else {
    character(0)
  }
  analysis_vars <- unique(c(
    fit_args$outcome,
    fit_args$treatment,
    formula_vars,
    by_vars
  ))

  # Two distinct problems, both undermining the imputation: a variable the
  # analysis needs that mice never saw, and a variable mice has but never
  # used as a predictor (e.g. the outcome excluded from the predictor
  # matrix, the classic source of uncongeniality). A column of the
  # predictor matrix that is entirely zero is never used to impute anything.
  absent <- setdiff(analysis_vars, known)
  present <- intersect(analysis_vars, known)
  pred <- imp$predictorMatrix
  unused <- character(0)
  if (!is.null(pred) && length(present) > 0L) {
    in_pred <- intersect(present, colnames(pred))
    unused <- in_pred[vapply(
      in_pred,
      function(v) all(pred[, v] == 0),
      logical(1)
    )]
  }

  flagged <- unique(c(absent, unused))
  if (length(flagged) > 0L) {
    rlang::warn(
      c(
        "Some analysis variables are not used by the imputation model.",
        x = paste0(
          "Absent or never used as a predictor: ",
          paste0("`", flagged, "`", collapse = ", "),
          "."
        ),
        i = paste0(
          "Include all analysis variables (outcome, treatment, confounders, ",
          "effect modifiers) as predictors in `mice::mice()` to reduce bias ",
          "and uncongeniality."
        )
      ),
      class = "causatr_mi_missing_predictors"
    )
  }
  invisible(NULL)
}
