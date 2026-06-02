#' Fit one ICE backward-step model, pooled or stratified
#'
#' @description
#' Single dispatch point for the per-step outcome / pseudo-outcome model
#' fit inside `ice_iterate()`. When `stratified` is `NULL` this fits one
#' pooled model on `fit_data` (the historical behaviour). When
#' `stratified` names a baseline column, the fitting rows are split by the
#' levels of that column and one model is fitted per stratum, returned as
#' a named list keyed by the stratum value (as character).
#'
#' Stratifying on a baseline variable G is equivalent to fully interacting
#' every model term with G *and* estimating a separate dispersion per
#' stratum; fitting S independent models is the direct way to obtain that
#' and lets each stratum carry its own functional form. Because G is
#' constant within a stratum, any formula term that references G is dropped
#' before the per-stratum fit (it would be aliased / rank-deficient) via
#' [strip_stratum_terms()].
#'
#' A stratum whose model fails to fit (e.g. too few rows, separation)
#' yields a `NULL` entry rather than aborting the whole chain; rows in
#' that stratum then receive `NA` predictions from [ice_predict_step()]
#' and drop out of the next backward step exactly like censored rows.
#'
#' @param model_fn Function. The per-step fitter (e.g. `stats::glm`).
#' @param formula A formula for the pooled fit; for the stratified fit its
#'   stratum-referencing terms are stripped per stratum.
#' @param fit_data data.table. The fitting rows for this step (already
#'   masked to uncensored rows with a valid response).
#' @param family Family object passed to `model_fn` when it accepts a
#'   `family` argument.
#' @param weights Numeric vector aligned row-for-row with `fit_data`, or
#'   `NULL`. Subset per stratum in the stratified branch.
#' @param dots List of user `...` extras replayed into `model_fn` via
#'   [replay_fit()].
#' @param stratified Character scalar naming the baseline stratum column,
#'   or `NULL` for a pooled fit.
#'
#' @return Either a single fitted model object (pooled) or a named list of
#'   fitted model objects / `NULL` keyed by stratum value (stratified).
#'
#' @noRd
ice_fit_step <- function(
  model_fn,
  formula,
  fit_data,
  family,
  weights,
  dots,
  stratified
) {
  if (is.null(stratified)) {
    model_args <- list(formula = formula, data = fit_data)
    if (fn_accepts_family(model_fn)) {
      model_args$family <- family
    }
    if (!is.null(weights)) {
      model_args$weights <- weights
    }
    return(replay_fit(model_fn, model_args, dots))
  }

  # Stratified fit: split rows by stratum value and fit one model each.
  # `sort(unique(...))` fixes a deterministic stratum ordering so the
  # stored model list, the prediction merge, and the bootstrap refit all
  # agree on the keying.
  strata <- as.character(fit_data[[stratified]])
  levs <- sort(unique(strata))
  formula_strat <- strip_stratum_terms(formula, stratified)

  models <- stats::setNames(vector("list", length(levs)), levs)
  for (g in levs) {
    rows_g <- which(strata == g)
    data_g <- fit_data[rows_g]
    args_g <- list(formula = formula_strat, data = data_g)
    if (fn_accepts_family(model_fn)) {
      args_g$family <- family
    }
    if (!is.null(weights)) {
      args_g$weights <- weights[rows_g]
    }
    # Tolerate a single-stratum fit failure rather than aborting the
    # whole ICE chain: rows in a failed stratum become NA pseudo-outcomes
    # and drop from the next step, mirroring the censored-row contract.
    models[[g]] <- tryCatch(
      replay_fit(model_fn, args_g, dots),
      error = function(e) NULL
    )
  }
  models
}


#' Predict one ICE backward-step model, pooled or stratified
#'
#' @description
#' Companion to [ice_fit_step()]. For a pooled model this is a thin
#' wrapper over `stats::predict(type = "response")`. For a stratified
#' model list it predicts each stratum's rows with that stratum's model
#' and merges the results back into a single vector aligned row-for-row
#' with `newdata`.
#'
#' Rows whose stratum has a `NULL` model (a failed per-stratum fit) or a
#' stratum value never seen at fitting time receive `NA`, which the
#' backward recursion treats like a censored row.
#'
#' @param models_k Either a single fitted model (pooled) or a named list
#'   of per-stratum models keyed by stratum value (stratified).
#' @param newdata data.table. Prediction rows (already intervention-modified
#'   where relevant). Must carry the stratum column in the stratified case.
#' @param stratified Character scalar naming the stratum column, or `NULL`.
#'
#' @return Numeric vector of fitted responses, length `nrow(newdata)`.
#'
#' @noRd
ice_predict_step <- function(models_k, newdata, stratified) {
  if (is.null(stratified)) {
    return(as.numeric(
      stats::predict(models_k, newdata = newdata, type = "response")
    ))
  }

  preds <- rep(NA_real_, nrow(newdata))
  strata <- as.character(newdata[[stratified]])
  for (g in names(models_k)) {
    model_g <- models_k[[g]]
    if (is.null(model_g)) {
      next
    }
    rows_g <- which(strata == g)
    if (length(rows_g) == 0L) {
      next
    }
    preds[rows_g] <- as.numeric(
      stats::predict(model_g, newdata = newdata[rows_g], type = "response")
    )
  }
  preds
}


#' Drop formula terms that reference the stratifying variable
#'
#' @description
#' Within a single stratum the stratifying variable G is constant, so any
#' RHS term referencing G (a main effect `G`, an interaction `A:G`, a
#' transform `factor(G)`) is collinear with the intercept and must be
#' removed before the per-stratum model is fitted -- otherwise the fit is
#' rank-deficient and emits aliasing warnings. Terms are matched by their
#' underlying variable names (via `all.vars()`), so transformed terms are
#' handled correctly. If stripping removes every RHS term the formula
#' collapses to an intercept-only `response ~ 1`.
#'
#' @param formula A model formula.
#' @param stratified Character scalar naming the stratifying variable.
#'
#' @return A formula with all stratum-referencing RHS terms removed.
#'
#' @noRd
strip_stratum_terms <- function(formula, stratified) {
  tt <- stats::terms(formula)
  term_labels <- attr(tt, "term.labels")
  response <- as.character(rlang::f_lhs(formula))

  keep <- vapply(
    term_labels,
    function(term) {
      !(stratified %in% all.vars(parse(text = term)[[1L]]))
    },
    logical(1L)
  )
  kept <- term_labels[keep]

  if (length(kept) == 0L) {
    stats::reformulate("1", response = response)
  } else {
    stats::reformulate(kept, response = response)
  }
}
