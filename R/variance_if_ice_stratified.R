#' Per-individual IF for one stratified ICE intervention
#'
#' @description
#' Stratified counterpart of `variance_if_ice_one()`. When
#' `causat(..., stratified = "G")` is used, each backward step fits one
#' outcome / pseudo-outcome model per level of the baseline stratum `G`,
#' so the stacked estimating-equation system is **per-stratum x per-time**
#' rather than per-time. Two structural facts make the sandwich a thin
#' wrapper over the pooled engine:
#'
#' 1. **Block-diagonal across strata.** Because `G` is baseline
#'    (time-invariant) and individuals never cross strata, a stratum-`g`
#'    row enters only the \eqn{\beta_{g,k}} scores. The per-step bread is
#'    block-diagonal in `g`, so each stratum's Channel-2 correction chain
#'    runs independently on that stratum's rows and models, then the
#'    per-stratum IF blocks add into a single length-`n` vector.
#' 2. **Global Channel 1.** The estimand \eqn{\hat\mu = E[Y^{\bar d}]} is
#'    marginal over all strata, so the sampling term and the centring mean
#'    are computed once over the whole target population (in
#'    `ice_if_setup()`); only the nuisance correction partitions by
#'    stratum. The first-step gradient \eqn{\partial\hat\mu/\partial
#'    \beta_{g,1}} sums over target rows in stratum `g` but keeps the
#'    global \eqn{1/n_t} normaliser -- which `variance_if_ice_chain()`
#'    obtains automatically once its prediction frame is restricted to the
#'    stratum.
#'
#' Equivalently, this is the pooled ICE sandwich run `S` times on disjoint
#' row sets, plus the shared global Channel 1 -- mirroring the exact
#' point-estimate equivalence "stratified ICE == pooled ICE per stratum
#' subset" that `ice_fit_step()` relies on.
#'
#' @param fit A `causatr_fit` object (ICE estimator, `stratified` set).
#' @param ice_result Per-intervention ICE result from `ice_iterate()`.
#'   `ice_result$models[[step]]` is a named list of per-stratum models.
#' @param target Logical vector identifying target-population rows.
#'
#' @return Numeric vector of length `n` (individuals), the per-individual IF.
#'
#' @noRd
variance_if_ice_one_stratified <- function(fit, ice_result, target) {
  ctx <- ice_if_setup(fit, ice_result, target)
  stratified <- fit$details$stratified
  models <- ice_result$models
  fit_ids_list <- ice_result$fit_ids

  # id -> baseline stratum lookup. `G` is constant within id, so the
  # first-time row is a faithful representative for every person-period
  # row of that individual.
  data <- fit$data
  first_time <- ctx$time_points[1]
  rows_first <- data[[ctx$time_col]] == first_time
  first_ids <- as.character(data[rows_first][[ctx$id_col]])
  first_strata <- as.character(data[rows_first][[stratified]])
  id_stratum <- stats::setNames(first_strata, first_ids)

  # Stratum levels in the same deterministic order `ice_fit_step()` keyed
  # the stored model list by (`sort(unique(as.character(G)))`), so the
  # per-stratum model lookup `models[[step]][[g]]` and the per-stratum fit
  # ids agree on the level coding (character, including factor labels).
  levs <- sort(unique(first_strata))

  IF_vec <- ctx$IF_vec
  for (g in levs) {
    # Per-step model for this stratum; `NULL` when the stratum was absent
    # at that step or its fit failed (the chain skips NULL models, leaving
    # a zero IF block -- those rows already dropped from the point estimate).
    models_g <- lapply(models, function(mk) mk[[g]])

    # Per-step fit ids restricted to stratum `g`, order preserved so they
    # align row-for-row with `model.matrix(models_g[[step]])` (the stratum
    # model was fit on the same-ordered subset inside `ice_fit_step()`).
    fit_ids_g <- lapply(fit_ids_list, function(ids) {
      ids[id_stratum[ids] == g]
    })

    IF_vec <- IF_vec +
      variance_if_ice_chain(
        ctx,
        models_g,
        fit_ids_g,
        restrict = list(col = stratified, val = g)
      )
  }

  IF_vec
}
