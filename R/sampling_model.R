#' Sampling model fitting for transportability and generalizability
#'
#' @description
#' Fits a logistic model \eqn{P(S = 1 \mid L)} for the study-membership
#' indicator S. Study rows have S = 1 (observed A and Y); target rows have
#' S = 0 (only L is observed). The model is fit on all rows (both S = 1
#' and S = 0) to estimate the covariate-distribution shift between
#' populations.
#'
#' Sampling weights derived from this model are used to reweight the study
#' population so that it resembles the target population — enabling
#' transportability (target is S = 0) and generalizability (target is all
#' rows, \eqn{S = 0 \cup S = 1}).
#'
#' @name sampling_model
#' @keywords internal
NULL


#' Fit a sampling model for transportability
#'
#' @description
#' Fits \eqn{P(S = 1 \mid L)} via a logistic model. The model is always
#' Bernoulli/logistic because S is binary. The formula is
#' \eqn{S \sim \text{confounders}} — treatment A does not appear on the
#' right-hand side, because S is a baseline property that precedes treatment
#' assignment.
#'
#' Unlike the censoring model, which tolerates NA values in the censoring
#' column, the sampling model requires a complete S column. Upstream
#' validation in `check_transport_inputs()` enforces this before this
#' function is called.
#'
#' @param data A data.table with columns for confounders and the sampling
#'   indicator.
#' @param target Character. Column name of the sampling indicator
#'   (1 = study, 0 = target).
#' @param confounders One-sided formula for confounders (baseline L).
#' @param model_fn Model fitting function (default `stats::glm`). Must
#'   accept `(formula, data, family, ...)`.
#' @param weights Optional numeric vector of external weights (e.g. survey
#'   weights). Length must equal `nrow(data)`.
#' @param ... Additional arguments passed to `model_fn`.
#'
#' @return A `causatr_sampling_model` S3 object with:
#'   - `model`: the fitted GLM
#'   - `sampling_formula`: the formula used
#'   - `fit_rows`: logical vector of rows used for fitting
#'   - `gamma_hat`: coefficient vector
#'   - `X_fit`: design matrix at fit time
#'   - `p_study`: fitted P(S=1|L) for fit rows
#'   - `p_marginal`: marginal P(S=1) across fit rows
#'
#' @noRd
fit_sampling_model <- function(
  data,
  target,
  confounders,
  treatment = NULL,
  model_fn = stats::glm,
  weights = NULL,
  ...
) {
  check_string(target)
  check_col_exists(data, target)
  check_formula(confounders)

  s_vals <- data[[target]]

  # Defensive NA check: upstream check_transport_inputs() should have
  # caught NAs, but reject here too so fit_sampling_model() is safe
  # to call directly in tests.
  if (anyNA(s_vals)) {
    rlang::abort(
      c(
        paste0("Column `", target, "` contains NA values."),
        i = "The sampling indicator S must be complete (no NAs)."
      ),
      class = "causatr_transport_na_target"
    )
  }

  # Strip any confounder terms involving the treatment variable(s) before
  # building the sampling formula. S is a baseline property that precedes
  # treatment assignment, so treatment and treatment-covariate interactions
  # (e.g. A:L) are not valid predictors of S. They would also cause NA-row
  # exclusion for target-population rows where A is unobserved.
  all_terms <- attr(stats::terms(confounders), "term.labels")
  if (!is.null(treatment)) {
    trt_set <- unique(unlist(lapply(treatment, function(a) a)))
    keep <- vapply(all_terms, function(tm) {
      term_vars <- all.vars(stats::reformulate(tm))
      !any(term_vars %in% trt_set)
    }, logical(1))
    sampling_terms <- all_terms[keep]
  } else {
    sampling_terms <- all_terms
  }
  if (length(sampling_terms) == 0L) {
    sampling_terms <- "1"
  }

  # Fit on all rows (S = 1 and S = 0). Exclude rows with NA in the
  # sampling-model predictors (not treatment, which is NA on target rows
  # by design and has already been stripped above).
  sampling_vars <- unique(unlist(lapply(
    sampling_terms[sampling_terms != "1"],
    function(tm) all.vars(stats::reformulate(tm))
  )))
  fit_rows <- rep(TRUE, nrow(data))
  for (v in sampling_vars) {
    fit_rows <- fit_rows & !is.na(data[[v]])
  }
  fit_data <- data[fit_rows]

  # Build formula: .sampling_s ~ baseline_confounders (no treatment on RHS)
  sampling_formula <- stats::reformulate(sampling_terms, response = ".sampling_s")

  fit_data$.sampling_s <- as.integer(fit_data[[target]])

  # The sampling model is always logistic regardless of estimator family.
  # S is binary 0/1 by definition; the treatment family describes a
  # different variable and has no bearing on modeling S.
  model_args <- list(
    formula = sampling_formula,
    data = fit_data,
    family = stats::binomial()
  )
  weights_fit <- if (is.null(weights)) NULL else weights[fit_rows]
  if (!is.null(weights_fit)) {
    model_args$weights <- weights_fit
  }
  model <- replay_fit(model_fn, model_args, list(...))

  p_study <- stats::fitted(model)
  p_marginal <- mean(as.integer(fit_data[[target]]))

  gamma_hat <- stats::coef(model)
  X_fit <- stats::model.matrix(model)

  structure(
    list(
      model = model,
      sampling_formula = sampling_formula,
      fit_rows = fit_rows,
      gamma_hat = gamma_hat,
      X_fit = X_fit,
      p_study = p_study,
      p_marginal = p_marginal
    ),
    class = "causatr_sampling_model"
  )
}


#' @export
print.causatr_sampling_model <- function(x, ...) {
  cat("<causatr_sampling_model>\n")
  cat("  Formula: ", deparse(x$sampling_formula), "\n", sep = "")
  cat("  Observations:", sum(x$fit_rows), "\n")
  cat("  Marginal P(S=1):", round(x$p_marginal, 3), "\n")
  cat("  Coefficients:", length(x$gamma_hat), "\n")
  invisible(x)
}
