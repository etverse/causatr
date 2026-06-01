#' Compute per-intervention longitudinal IPW bundles
#'
#' @description
#' Longitudinal companion to `compute_ipw_contrast_point()`. For each
#' intervention:
#'
#' 1. Apply the intervention at every period (the same intervention
#'    is broadcast across `time_points[1:K]`).
#' 2. Build the cumulative density-ratio weight per individual via
#'    `compute_longitudinal_weights()`.
#' 3. Multiply external (survey / IPCW) weights from
#'    `fit$details$weights` evaluated on the final-period rows.
#' 4. Refit `Y ~ 1` weighted on the final-period rows; read
#'    \eqn{\hat\mu_a} off the intercept (Hajek).
#'
#' Returns the per-intervention bundles + the scalars
#' `compute_contrast()` needs to dispatch variance / bootstrap.
#'
#' @param fit A `causatr_fit` from `fit_longitudinal_ipw()`.
#' @param interventions Named list of `causatr_intervention` objects
#'   (or `NULL` for natural course).
#' @param target_within_first Logical vector flagging the target
#'   population at the first time point (mirrors ICE's convention).
#'
#' @return A list with `bundles`, `mu_hat`, `id_first` (character ids
#'   in target order), and `n_total` (`nrow(fit$data)`).
#'
#' @noRd
compute_ipw_contrast_longitudinal <- function(
  fit,
  interventions,
  target_within_first,
  trim = 1
) {
  data <- fit$data
  treatment <- fit$treatment
  outcome <- fit$outcome
  details <- fit$details
  treatment_models_by_time <- details$treatment_models_by_time
  numerator_models_by_time <- details$numerator_models_by_time
  fit_data_by_time <- details$fit_data_by_time
  time_points <- details$time_points
  n_times <- details$n_times
  id_col <- fit$id
  time_col <- fit$time
  fit_rows_final <- details$fit_rows_final
  fam_obj <- resolve_family(fit$family)
  model_fn <- details$model_fn
  ext_w <- details$weights

  # First-time-point rows define the per-individual weight ordering.
  # We average / aggregate per id over the final-period rows; the
  # natural id order comes from sorting the first-period subset.
  first_time <- time_points[1]
  rows_first <- data[[time_col]] == first_time
  ids_first <- as.character(data[rows_first][[id_col]])
  n_id <- length(ids_first)

  # Final-period rows are where we evaluate the MSM. Each id has one
  # final-period row; the cumulative weight per id is multiplied with
  # the external weight on that row before fitting the weighted Hajek
  # intercept.
  last_time <- time_points[n_times]
  data_final <- data[fit_rows_final]
  ids_final <- as.character(data_final[[id_col]])
  n_final <- nrow(data_final)

  # Map id -> position in the first-time ordering. Used to project the
  # cumulative weight (length n_id, in first-time order) onto the
  # final-period rows. Reverse lookup: given an id from `ids_final`,
  # `id_to_first_idx[id]` gives its row in `w_id` from
  # `compute_longitudinal_weights()`.
  id_to_first_idx <- stats::setNames(seq_len(n_id), ids_first)

  ext_w_final <- if (is.null(ext_w)) NULL else ext_w[fit_rows_final]

  # Transport: compute per-final-row sampling odds weights once (shared
  # across all interventions). The sampling model was fit on first-period
  # (cross-sectional) data; it predicts correctly on any subset that
  # carries the same baseline confounders, including the final-period rows.
  is_transport <- isTRUE(fit$details$transport)
  w_S_final <- NULL
  if (is_transport) {
    w_S_final <- compute_sampling_weights(
      fit$details$sampling_model,
      data_final,
      fit$target,
      fit$target_subset
    )
  }

  # MSM formula. Default `Y ~ 1` (intercept-only Hájek). With
  # baseline-modifier effect modification (`A:modifier` in
  # `confounders`), expands to `Y ~ 1 + modifier` via the existing
  # `build_ipw_msm_formula()` so `predict()` returns
  # stratum-specific counterfactual means.
  em_info <- details$em_info
  msm_formula <- build_ipw_msm_formula(outcome, em_info)

  int_names <- names(interventions)

  # Univariate `ipsi(delta)` is supported under longitudinal IPW: the
  # Kennedy (2019) closed-form weight extends to a per-period product
  # (each period's weight comes from `compute_density_ratio_weights()`'s
  # IPSI branch). Two combinations remain unsupported and are rejected
  # here with classed errors:
  #
  #   * Multivariate IPSI — the closed form is binary-univariate; the
  #     per-component density chain has no IPSI shortcut.
  #   * Stabilized IPSI (`stabilize = "marginal"`) — Kennedy's weight is
  #     already a bounded ratio acting on the propensity, so there is no
  #     separate marginal numerator to stabilize against; the stabilized
  #     per-period closures (`compute_stabilized_period_weight()`,
  #     `make_long_stabilized_period_closure()`) have no IPSI branch.
  is_stabilized <- !is.null(numerator_models_by_time)
  for (nm in int_names) {
    iv <- interventions[[nm]]
    if (is.null(iv)) {
      next
    }
    # Multivariate: iv is a plain list of per-component interventions
    if (is.list(iv) && !inherits(iv, "causatr_intervention")) {
      for (comp_nm in names(iv)) {
        sub_iv <- iv[[comp_nm]]
        if (!is.null(sub_iv) && sub_iv$type == "ipsi") {
          rlang::abort(
            c(
              "`ipsi()` is not supported under multivariate longitudinal IPW.",
              i = paste0(
                "Intervention `",
                nm,
                "`, component `",
                comp_nm,
                "` is `ipsi()`."
              ),
              i = "Use `static()` / `shift()` / `scale_by()` / `dynamic()` per component, switch to point IPW, or use univariate longitudinal IPW (`ipsi()` is supported there)."
            ),
            class = "causatr_longitudinal_ipsi_multivariate"
          )
        }
      }
    } else if (
      inherits(iv, "causatr_intervention") &&
        iv$type == "ipsi" &&
        is_stabilized
    ) {
      rlang::abort(
        c(
          "`ipsi()` is not supported with `stabilize = 'marginal'` under longitudinal IPW.",
          i = paste0("Intervention `", nm, "` is `ipsi()`."),
          i = "Kennedy's IPSI weight is already a bounded propensity-odds ratio with no separate marginal numerator. Set `stabilize = 'none'` to run longitudinal IPSI."
        ),
        class = "causatr_longitudinal_ipsi_stabilize"
      )
    }
  }

  bundles <- lapply(int_names, function(nm) {
    iv <- interventions[[nm]]

    # Per-id cumulative weight under the intervention. Internally the
    # routine loops over `time_points` and multiplies per-period
    # density ratios, returning a length-n_id vector ordered by
    # `ids_first`. NULL intervention -> all-ones (natural course).
    w_id <- compute_longitudinal_weights(
      treatment_models_by_time = treatment_models_by_time,
      fit_data_by_time = fit_data_by_time,
      ids_first = ids_first,
      id_col = id_col,
      intervention = iv,
      numerator_models_by_time = numerator_models_by_time,
      trim = trim
    )

    # Project per-id weights onto the final-period row order.
    w_final <- w_id[id_to_first_idx[ids_final]]

    # Compose: treatment density-ratio × sampling odds × external weights.
    # All three sources are independent reweightings; order is arbitrary.
    w_combined <- w_final
    if (!is.null(w_S_final)) {
      w_combined <- w_combined * w_S_final
    }
    if (!is.null(ext_w_final)) {
      w_combined <- w_combined * ext_w_final
    }

    # Final-period weighted MSM. Intercept-only Hajek under
    # `Y ~ 1`; with EM the MSM is `Y ~ 1 + modifier` so `predict()`
    # returns stratum-specific counterfactual means. The treatment is
    # absorbed by the cumulative density-ratio weight.
    msm_args <- list(
      formula = msm_formula,
      data = data_final,
      family = msm_family(fam_obj),
      weights = w_combined
    )
    msm_model <- do.call(model_fn, msm_args)

    # Counterfactual marginal mean over the target subset of
    # individuals. Under transport, all final-period rows are study rows
    # (S=1, enforced by fit_rows_final); the sampling weights already
    # reweight to the target population, so every row contributes.
    # Otherwise `target_within_first` restricts to the user's subset.
    target_ids_final <- if (is_transport) {
      rep(TRUE, n_final)
    } else {
      target_within_first[id_to_first_idx[ids_final]]
    }
    preds_final <- stats::predict(msm_model, type = "response")
    valid_final <- target_ids_final & !is.na(preds_final)
    w_target_final <- if (!is.null(ext_w_final)) {
      ext_w_final[valid_final]
    } else {
      NULL
    }
    mu_hat_iv <- maybe_weighted_mean(preds_final[valid_final], w_target_final)

    list(
      intervention = iv,
      msm_model = msm_model,
      weights_id = w_id,
      weights_final = w_combined,
      preds_final = preds_final,
      valid_final = valid_final,
      target_ids_final = target_ids_final,
      mu_hat = mu_hat_iv
    )
  })
  names(bundles) <- int_names

  mu_hat <- vapply(bundles, function(b) b$mu_hat, numeric(1))
  names(mu_hat) <- int_names

  list(
    bundles = bundles,
    mu_hat = mu_hat,
    ids_first = ids_first,
    id_to_first_idx = id_to_first_idx,
    ids_final = ids_final,
    n_id = n_id,
    n_final = n_final,
    n_total = nrow(data)
  )
}
