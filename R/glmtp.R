#' Augmented-data sequential regression for natural-history MTPs
#'
#' @description
#' Runs the augmented-data g-computation of Diaz, Williams, Morzywolek &
#' Rudolph (2026, arXiv:2605.24167) for a modified treatment policy that depends
#' on the **natural-value history of treatment**. Unlike [ice_iterate()] (which
#' conditions on the observed lagged treatment), this engine decouples the
#' *conditioning* treatment value from the *policy-input* value by carrying an
#' auxiliary natural-treatment-history sequence \eqn{\bar{s}_t} as a label
#' through the backward recursion.
#'
#' Let \eqn{\mathcal{A}} be the discrete treatment support and
#' \eqn{\bar{\mathcal{A}}_t = \mathcal{A}^t} the natural-history support. The
#' recursion is:
#' \itemize{
#'   \item Base (t = tau): fit \eqn{m_\tau = E[Y \mid A_\tau, H_\tau]} -- one
#'     pooled model, label-independent.
#'   \item For t = tau-1, ..., 1 and each label \eqn{\bar{s}_t}: regress the
#'     label-specific pseudo-response \eqn{q_{t+1}(\bar{s}_t, A_{t+1}, H_{t+1})}
#'     on \eqn{(A_t, H_t)} (the same design over all uncensored rows, a
#'     different response per label) to obtain \eqn{m_t(\bar{s}_t, \cdot,
#'     \cdot)}.
#'   \item Predict \eqn{q_t(\bar{s}_{t-1}, A_{t,i}, H_{t,i}) =
#'     m_t((\bar{s}_{t-1}, A_{t,i}), d_t((\bar{s}_{t-1}, A_{t,i}), H_{t,i}),
#'     H_{t,i})} -- append the observed natural value \eqn{A_{t,i}} to the prior
#'     label, then predict that label's model at the policy-shifted treatment.
#' }
#' The estimand is \eqn{\hat\theta = n^{-1} \sum_i q_1(A_{1,i}, H_{1,i})}.
#' Binary outcomes use `quasibinomial` for the pseudo steps (fractional
#' responses), exactly as [ice_iterate()]. Censoring and external weights flow
#' through every per-(time, label) fit.
#'
#' @param fit A `causatr_fit` of type `"longitudinal"`, `estimator = "gcomp"`
#'   (from `fit_ice()`).
#' @param intervention A `causatr_glmtp` object (from [grace_period()] /
#'   [carry_forward()]), or `NULL` for the natural course (routed through the
#'   identity policy).
#'
#' @return A list with:
#'   \describe{
#'     \item{`pseudo_final`}{Numeric vector of \eqn{\hat q_1}, one per
#'       individual at the first time point.}
#'     \item{`models`}{List indexed by time point. The final element is the
#'       single label-independent base model; each earlier element is a named
#'       list of per-label models keyed by [glmtp_label_key()].}
#'     \item{`fit_ids`}{Parallel structure of fitting-set ids (for variance).}
#'     \item{`support`}{The discrete treatment support used.}
#'     \item{`intervention`}{The intervention (or `NULL`).}
#'   }
#'
#' @references
#' Diaz I, Williams NT, Morzywolek P, Rudolph KE (2026). Modified treatment
#' policies that depend on the natural history of treatment. arXiv:2605.24167.
#'
#' @seealso [grace_period()], [carry_forward()], [ice_iterate()]
#' @family glmtp
#' @noRd
glmtp_iterate <- function(fit, intervention) {
  data <- fit$data
  details <- fit$details
  outcome <- fit$outcome
  treatment <- fit$treatment
  id_col <- fit$id
  time_col <- fit$time
  censoring <- fit$censoring

  time_points <- details$time_points
  n_times <- details$n_times
  baseline_terms <- details$baseline_terms
  tv_vars <- details$tv_vars
  tv_terms <- details$tv_terms
  max_lag <- details$max_lag
  em_info <- details$em_info
  model_fn <- details$model_fn
  family_outcome <- details$family_outcome
  family_pseudo <- details$family_pseudo
  external_weights <- details$weights
  model_fn_dots <- details$dots
  # Treatment design term(s) for the per-step formula (bare column name by
  # default; transformed labels under a user `treatment_form`). The policy
  # always sets the numeric treatment column -- only the design term changes.
  treatment_terms <- details$treatment_terms

  # The natural-course reference (NULL) reuses the augmented machinery with the
  # identity policy so a contrast that mixes glmtp arms with a natural-course
  # reference runs through one engine.
  policy <- if (is.null(intervention)) {
    glmtp_identity_policy
  } else {
    intervention$policy
  }
  budget <- if (is.null(intervention)) 1024L else intervention$budget %||% 1024L

  # Discrete support + tractability. (Routing also validates, but re-checking
  # here keeps the bootstrap refit path honest on resampled data.)
  support <- glmtp_support(data, treatment, censoring)
  glmtp_check_tractable(support, n_times, budget)

  uncens <- is_uncensored(data, censoring)
  all_ids <- unique(data[[id_col]])
  id_chr <- as.character(all_ids)
  n_ind <- length(all_ids)

  models <- vector("list", n_times)
  names(models) <- as.character(time_points)
  fit_ids <- vector("list", n_times)
  names(fit_ids) <- as.character(time_points)

  # -- Base step t = tau: m_tau = E[Y | A_tau, H_tau], fit once on observed Y.
  final_time <- time_points[n_times]
  final_idx <- n_times - 1L
  mask_final <- data[[time_col]] == final_time
  fit_mask <- mask_final & uncens & !is.na(data[[outcome]])
  fit_data <- data[fit_mask]
  fit_ids[[n_times]] <- as.character(fit_data[[id_col]])

  formula_base <- ice_build_formula(
    outcome,
    treatment,
    baseline_terms,
    tv_vars,
    tv_terms,
    final_idx,
    max_lag,
    fit_data,
    em_info,
    treatment_terms
  )
  w_final <- if (!is.null(external_weights)) external_weights[fit_mask]
  m_base <- ice_fit_step(
    model_fn,
    formula_base,
    fit_data,
    family_outcome,
    w_final,
    model_fn_dots,
    NULL
  )
  models[[n_times]] <- m_base

  # Predict q_tau over all uncensored final-time rows, once per length-(tau-1)
  # history label. `Q` is a list keyed by the prior-label key; each element is
  # a length-n_ind vector named by id.
  pred_mask <- mask_final & uncens
  pred_ids <- as.character(data[pred_mask][[id_col]])
  pred_base <- data[pred_mask]
  a_now_final <- pred_base[[treatment]]

  labels_prior <- glmtp_enumerate_labels(support, n_times - 1L)
  Q <- stats::setNames(
    vector("list", length(labels_prior)),
    vapply(labels_prior, glmtp_label_key, character(1))
  )
  for (li in seq_along(labels_prior)) {
    ell <- labels_prior[[li]]
    d_vals <- policy(ell, a_now_final, pred_base, n_times, n_times)
    nd <- data.table::copy(pred_base)
    nd[, (treatment) := d_vals]
    preds <- as.numeric(stats::predict(m_base, newdata = nd, type = "response"))
    qv <- stats::setNames(rep(NA_real_, n_ind), id_chr)
    qv[pred_ids] <- preds
    Q[[li]] <- qv
  }

  # Count per-label model fits that were attempted (had >= 1 valid response
  # row) but failed (separation, rank deficiency). Their rows fall to NA and
  # drop from later steps; a non-zero count is surfaced as a single classed
  # warning after the recursion so the rare degrade-to-NA path is never silent.
  n_failed_labels <- 0L

  # -- Backward steps t = tau-1 ... 1.
  for (step_i in seq(n_times - 1L, 1L, by = -1L)) {
    current_time <- time_points[step_i]
    time_idx <- step_i - 1L
    mask_current <- data[[time_col]] == current_time
    mask_uncens <- mask_current & uncens

    labels_t <- glmtp_enumerate_labels(support, step_i)
    labels_prior_new <- glmtp_enumerate_labels(support, step_i - 1L)

    cur_all <- data[mask_current]
    pred_ids_all <- as.character(cur_all[[id_col]])
    a_now_cur <- cur_all[[treatment]]

    cur_uncens <- data[mask_uncens]
    cur_uncens_ids <- as.character(cur_uncens[[id_col]])
    w_uncens <- if (!is.null(external_weights)) external_weights[mask_uncens]

    # Fit one model per length-t label. All labels share the (A_t, H_t) design;
    # only the pseudo-response Q[[label]] changes. A per-label fit that fails
    # (separation, empty response) becomes NULL and its rows fall to NA below,
    # mirroring the censored-row contract.
    step_keys <- vapply(labels_t, glmtp_label_key, character(1))
    step_models <- stats::setNames(vector("list", length(labels_t)), step_keys)
    step_fit_ids <- stats::setNames(vector("list", length(labels_t)), step_keys)
    for (lj in seq_along(labels_t)) {
      # `step_keys[lj]` is `glmtp_label_key(labels_t[[lj]])` by construction
      # (same order); derive it directly so the response selection never relies
      # on positional alignment between two separate vectors.
      key_j <- glmtp_label_key(labels_t[[lj]])
      resp <- Q[[key_j]][cur_uncens_ids]
      has_pseudo <- !is.na(resp)
      if (sum(has_pseudo) == 0L) {
        next
      }
      fit_data <- data.table::copy(cur_uncens[has_pseudo])
      fit_data[, .pseudo_y := resp[has_pseudo]]
      step_fit_ids[[lj]] <- as.character(fit_data[[id_col]])
      formula_k <- ice_build_formula(
        ".pseudo_y",
        treatment,
        baseline_terms,
        tv_vars,
        tv_terms,
        time_idx,
        max_lag,
        fit_data,
        em_info,
        treatment_terms
      )
      w_k <- if (!is.null(external_weights)) w_uncens[has_pseudo]
      step_models[[lj]] <- tryCatch(
        ice_fit_step(
          model_fn,
          formula_k,
          fit_data,
          family_pseudo,
          w_k,
          model_fn_dots,
          NULL
        ),
        error = function(e) NULL
      )
      # A label with valid rows that still failed to fit -> count it.
      if (is.null(step_models[[lj]])) {
        n_failed_labels <- n_failed_labels + 1L
      }
    }
    models[[step_i]] <- step_models
    fit_ids[[step_i]] <- step_fit_ids

    # Build Q for the next (shorter) labels. For prior label `ell` and
    # individual i, the natural history is (ell, A_{t,i}); group by the
    # observed value appended so each group uses its label's model.
    Q_new <- stats::setNames(
      vector("list", length(labels_prior_new)),
      vapply(labels_prior_new, glmtp_label_key, character(1))
    )
    for (li in seq_along(labels_prior_new)) {
      ell <- labels_prior_new[[li]]
      qv <- stats::setNames(rep(NA_real_, n_ind), id_chr)
      for (a_val in support) {
        rows_a <- which(a_now_cur == a_val)
        if (length(rows_a) == 0L) {
          next
        }
        key_full <- glmtp_label_key(c(ell, a_val))
        m_lab <- step_models[[key_full]]
        if (is.null(m_lab)) {
          next
        }
        sub <- cur_all[rows_a]
        d_vals <- policy(
          ell,
          rep(a_val, length(rows_a)),
          sub,
          step_i,
          n_times
        )
        nd <- data.table::copy(sub)
        nd[, (treatment) := d_vals]
        qv[pred_ids_all[rows_a]] <- as.numeric(
          stats::predict(m_lab, newdata = nd, type = "response")
        )
      }
      Q_new[[li]] <- qv
    }
    Q <- Q_new
  }

  # Surface any per-label fit failures once. Their rows degraded to NA and
  # dropped from later steps, so the estimate is over a subset -- warn rather
  # than return a quietly-biased number.
  if (n_failed_labels > 0L) {
    rlang::warn(
      c(
        paste0(
          n_failed_labels,
          " per-label outcome model(s) failed to fit in the natural-history ",
          "MTP recursion."
        ),
        i = paste0(
          "Affected individuals received NA pseudo-outcomes and dropped from ",
          "later steps; the counterfactual mean may be biased. Check for ",
          "separation / sparse history cells, or fit with a regularised ",
          "`model_fn`."
        )
      ),
      class = "causatr_glmtp_label_fit_failed",
      .frequency = "once",
      .frequency_id = "causatr_glmtp_label_fit_failed"
    )
  }

  # After the loop Q holds exactly the single empty-label element (its name is
  # "", which R cannot index with `[[`), so read it positionally with an
  # explicit guard against any future change that breaks the single-label
  # invariant.
  if (length(Q) != 1L) {
    rlang::abort(
      "Internal error: the G-LMTP recursion did not collapse to a single label.",
      class = "causatr_glmtp_internal"
    )
  }
  q1 <- Q[[1L]]
  first_time <- time_points[1]
  first_ids <- as.character(data[data[[time_col]] == first_time][[id_col]])
  pseudo_final <- unname(q1[first_ids])

  list(
    pseudo_final = pseudo_final,
    models = models,
    fit_ids = fit_ids,
    support = support,
    intervention = intervention
  )
}
