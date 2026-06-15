#' Sandwich variance dispatcher for natural-history MTP interventions
#'
#' @description
#' Outer loop called by `contrast()` when `ci_method = "sandwich"` and the
#' contrast uses at least one `causatr_glmtp` intervention. Loops over
#' interventions, calls `variance_if_glmtp_one()` for each, and aggregates
#' into a `k x k` variance-covariance matrix via `vcov_from_if()`. Mirrors
#' the structure of `variance_if_ice()`.
#'
#' @param fit A `causatr_fit` of type `"longitudinal"`, `estimator = "gcomp"`.
#' @param glmtp_results Named list of `glmtp_iterate()` result objects, one per
#'   intervention (in the same name order as `int_names`).
#' @param target_within_first Logical vector over first-time individuals
#'   identifying the target population.
#' @param cluster_vec Optional length-n numeric/character cluster vector (for
#'   cluster-robust aggregation, currently unused by G-LMTP).
#'
#' @return A named `k x k` variance-covariance matrix.
#'
#' @seealso `variance_if_glmtp_one()`, `glmtp_iterate()`, `variance_if_ice()`
#' @family glmtp
#' @noRd
variance_if_glmtp <- function(
  fit,
  glmtp_results,
  target_within_first,
  cluster_vec = NULL
) {
  int_names <- names(glmtp_results)
  data <- fit$data
  first_time <- fit$details$time_points[1]
  rows_first <- data[[fit$time]] == first_time
  n <- sum(rows_first)

  cluster_first <- if (is.null(cluster_vec)) NULL else cluster_vec[rows_first]

  IF_list <- lapply(glmtp_results, function(res) {
    variance_if_glmtp_one(fit, res, target_within_first)
  })

  vcov_from_if(IF_list, n, int_names, cluster = cluster_first)
}


#' Per-individual IF for one G-LMTP intervention
#'
#' @description
#' Assembles the per-individual influence function for one natural-history MTP
#' from Channel 1 (sampling term) and Channel 2 (per-(time, label) nuisance
#' correction). The G-LMTP M-estimation system (Diaz, Williams, Morzywolek &
#' Rudolph 2026) stacks:
#' \enumerate{
#'   \item The estimand equation \eqn{(1/n) \sum_i (q_{1,i} - \theta) = 0}.
#'   \item The base model GLM scores at the final time (label-independent).
#'   \item One set of GLM scores per (backward step, natural-history label).
#' }
#' The block-triangular bread back-substitution follows the same direction as
#' the pooled ICE chain (earliest backward step first, base model last), with
#' the per-label sensitivity dictionary `D` generalising the single `d_vec` of
#' `variance_if_ice_chain()`.
#'
#' @param fit A `causatr_fit` object (longitudinal gcomp, G-LMTP path).
#' @param glmtp_result Per-intervention result from `glmtp_iterate()`.
#' @param target Logical vector over first-time individuals (the target
#'   population).
#'
#' @return Numeric vector of length `n` (individuals), the per-individual IF.
#'
#' @seealso `glmtp_if_setup()`, `variance_if_glmtp_chain()`
#' @family glmtp
#' @noRd
variance_if_glmtp_one <- function(fit, glmtp_result, target) {
  ctx <- glmtp_if_setup(fit, glmtp_result, target)
  intervention <- glmtp_result$intervention
  policy <- if (is.null(intervention)) {
    glmtp_identity_policy
  } else {
    intervention$policy
  }
  support <- glmtp_result$support
  n_times <- ctx$n_times

  ctx$IF_vec +
    variance_if_glmtp_chain(
      ctx,
      glmtp_result$models,
      glmtp_result$fit_ids,
      support,
      policy,
      n_times
    )
}


#' Channel-1 IF and shared state for a G-LMTP intervention
#'
#' @description
#' Parallel to `ice_if_setup()` but without `data_iv` (G-LMTP has no single
#' counterfactual frame; each label uses a different policy-shifted frame
#' computed on the fly in the chain). The Channel-1 sampling term is identical
#' in form to ICE:
#' \deqn{\mathrm{IF}_i = n \cdot (w_i / \sum w) \cdot (q_{1,i} - \hat\theta)}
#' for target rows, zero elsewhere.
#'
#' @param fit A `causatr_fit` object (longitudinal gcomp).
#' @param glmtp_result Per-intervention result from `glmtp_iterate()`.
#' @param target Logical vector over first-time individuals.
#'
#' @return A context list consumed by `variance_if_glmtp_chain()`: `fit`,
#'   `data`, `time_points`, `n_times`, `id_col`, `time_col`, `all_ids`, `n`,
#'   `id_to_idx`, `target`, `w_t`, `sum_w_target`, `has_weights`, `w_at_step`,
#'   `IF_vec`.
#'
#' @seealso `ice_if_setup()`
#' @family glmtp
#' @noRd
glmtp_if_setup <- function(fit, glmtp_result, target) {
  data <- fit$data
  details <- fit$details
  pseudo_final <- glmtp_result$pseudo_final

  time_points <- details$time_points
  n_times <- details$n_times
  id_col <- fit$id
  time_col <- fit$time

  first_time <- time_points[1]
  rows_first <- data[[time_col]] == first_time
  all_ids <- as.character(data[rows_first][[id_col]])
  n <- length(all_ids)
  id_to_idx <- stats::setNames(seq_len(n), all_ids)

  ext_w <- details$weights
  has_weights <- !is.null(ext_w)
  if (has_weights) {
    w_first <- ext_w[rows_first]
    w_t <- w_first[target]
  } else {
    w_t <- rep(1, sum(target))
  }
  sum_w_target <- sum(w_t)
  mu_hat <- sum(w_t * pseudo_final[target]) / sum_w_target

  IF_vec <- numeric(n)
  IF_vec[target] <- n *
    (w_t / sum_w_target) *
    (pseudo_final[target] - mu_hat)

  # Per-time-step id-to-weight lookup. The w_{k-1} factor for the off-diagonal
  # bread block (same derivation as ice_if_setup; see variance_if_ice.R).
  w_at_step <- if (has_weights) {
    lapply(seq_len(n_times), function(k) {
      rows_k <- data[[time_col]] == time_points[k]
      ids_k <- as.character(data[rows_k][[id_col]])
      stats::setNames(ext_w[rows_k], ids_k)
    })
  } else {
    NULL
  }

  list(
    fit = fit,
    data = data,
    time_points = time_points,
    n_times = n_times,
    id_col = id_col,
    time_col = time_col,
    all_ids = all_ids,
    n = n,
    id_to_idx = id_to_idx,
    target = target,
    w_t = w_t,
    sum_w_target = sum_w_target,
    has_weights = has_weights,
    w_at_step = w_at_step,
    IF_vec = IF_vec
  )
}


#' Channel-2 nuisance-correction chain for a G-LMTP model sequence
#'
#' @description
#' Runs the forward per-step sensitivity cascade for the augmented-data
#' sequential regression (Diaz et al. 2026, arXiv:2605.24167). Unlike the
#' pooled ICE chain (one model per step), the G-LMTP chain carries
#' \eqn{|\mathcal{A}|^t} per-label models at forward step \eqn{t}, and the
#' sensitivity propagates through a **dictionary** `D` keyed by label:
#' \eqn{D[\bar s_t]_i} = per-individual sensitivity of \eqn{\hat\theta} to
#' changes in \eqn{Q[\bar s_t]_i} (the pseudo-response the \eqn{\bar s_t}
#' model was fit on).
#'
#' **Cascade (step 1 to \eqn{\tau}):**
#' \enumerate{
#'   \item \emph{Step 1 (earliest), label \eqn{\bar s_1 = [a]}:} gradient from
#'     TARGET rows with \eqn{A_{1,i} = a} (identical to the ICE step-1
#'     gradient, restricted to the label's support value). Calls
#'     `correct_model()` → stores \eqn{D[[a]] = d}.
#'   \item \emph{Steps 2 .. \eqn{\tau-1}, label \eqn{\bar s_t = (\bar s_{t-1},
#'     a_t)}:} gradient from \eqn{D[\bar s_{t-1}] \cdot w_{t-1}} over
#'     individuals with \eqn{A_{t,i} = a_t}. Calls `correct_model()` → stores
#'     \eqn{D[\bar s_t] = d}.
#'   \item \emph{Step \eqn{\tau} (base, label-independent):} gradient
#'     accumulated over ALL prior labels \eqn{\bar s_{\tau-1}}; each label
#'     contributes through its \eqn{D} vector with the policy-shifted
#'     counterfactual frame for that label.
#' }
#'
#' `correct_model()`, `iv_design_matrix()`, and `coef_clean()` are reused
#' verbatim from the pooled ICE engine.
#'
#' @param ctx Context list from `glmtp_if_setup()`.
#' @param models List indexed by step. `models[[n_times]]` is the
#'   label-independent base model; `models[[t]]` for `t < n_times` is a named
#'   list of per-label fitted models keyed by `glmtp_label_key()`.
#' @param fit_ids Parallel structure to `models`: `fit_ids[[n_times]]` is a
#'   character vector; `fit_ids[[t]]` is a named list of character vectors.
#' @param support Numeric vector. The discrete treatment support (from
#'   `glmtp_support()`).
#' @param policy The policy closure with signature
#'   `function(s_prior, a_now, h_data, t, n_times)`.
#' @param n_times Positive integer. Number of time points \eqn{\tau}.
#'
#' @return Numeric vector of length `ctx$n`, the accumulated Channel-2 IF.
#'
#' @seealso `variance_if_glmtp_one()`, `variance_if_ice_chain()`
#' @family glmtp
#' @noRd
variance_if_glmtp_chain <- function(
  ctx,
  models,
  fit_ids,
  support,
  policy,
  n_times
) {
  fit <- ctx$fit
  data <- ctx$data
  time_points <- ctx$time_points
  id_col <- ctx$id_col
  time_col <- ctx$time_col
  all_ids <- ctx$all_ids
  n <- ctx$n
  id_to_idx <- ctx$id_to_idx
  target <- ctx$target
  w_t <- ctx$w_t
  sum_w_target <- ctx$sum_w_target
  has_weights <- ctx$has_weights
  w_at_step <- ctx$w_at_step
  treatment <- fit$treatment

  # Sensitivity dictionary: D[[label_key]] = length-n per-individual
  # sensitivity of theta to Q[label]_i. Keyed by glmtp_label_key().
  D <- list()
  IF_acc <- numeric(n)

  # ------------------------------------------------------------------
  # Forward pass: per-label steps (step 1 .. n_times - 1)
  # ------------------------------------------------------------------
  for (step_i in seq_len(n_times - 1L)) {
    current_time <- time_points[step_i]
    mask_current <- data[[time_col]] == current_time
    cur_data <- data[mask_current]
    cur_ids <- as.character(cur_data[[id_col]])
    a_now_cur <- cur_data[[treatment]]

    labels_step <- glmtp_enumerate_labels(support, step_i)

    for (lab_seq in labels_step) {
      lab_key <- glmtp_label_key(lab_seq)
      model_k <- models[[step_i]][[lab_key]]
      if (is.null(model_k)) {
        next
      }
      fit_ids_k <- fit_ids[[step_i]][[lab_key]]
      if (is.null(fit_ids_k) || length(fit_ids_k) == 0L) {
        next
      }

      # Last element of the label: the observed treatment value that selects
      # which individuals this model covers.
      a_last <- lab_seq[step_i]

      # ---- Gradient w.r.t. beta_{lab, step_i} ----
      if (step_i == 1L) {
        # Step 1 (earliest time): gradient from TARGET rows with A_1 == a_last.
        # The estimand gradient decomposes by label because q_1_i = mu(X*_i
        # beta_{[A_1_i], 1}); only label-matching target individuals contribute.
        rows_a <- which(a_now_cur == a_last)
        if (length(rows_a) == 0L) {
          g_k <- rep(0, length(coef_clean(model_k)))
        } else {
          cur_data_a <- cur_data[rows_a]
          cur_ids_a <- cur_ids[rows_a]

          d_vals <- policy(
            numeric(0L),
            rep(a_last, length(rows_a)),
            cur_data_a,
            step_i,
            n_times
          )
          nd_a <- data.table::copy(cur_data_a)
          nd_a[[treatment]] <- d_vals

          X_star_a <- iv_design_matrix(model_k, nd_a)
          eta_star_a <- as.numeric(X_star_a %*% coef_clean(model_k))
          mu_eta_a <- model_k$family$mu.eta(eta_star_a)

          # Target indicator for these rows (NA if id not in all_ids, which
          # should not happen in a well-formed panel).
          target_in_a <- target[id_to_idx[cur_ids_a]]
          valid_a <- !is.na(target_in_a) & target_in_a

          if (any(valid_a)) {
            # Weights for the valid target rows in this label group.
            w_a <- if (has_weights) {
              w_t[match(cur_ids_a[valid_a], all_ids[target])]
            } else {
              rep(1, sum(valid_a))
            }
            g_k <- as.numeric(
              crossprod(
                X_star_a[valid_a, , drop = FALSE],
                w_a * mu_eta_a[valid_a]
              )
            ) /
              sum_w_target
          } else {
            g_k <- rep(0, length(coef_clean(model_k)))
          }
        }
      } else {
        # Step t > 1: gradient from D[prior_label] * w_{t-1} over individuals
        # with A_t == a_last. The prior label is lab_seq without its last element.
        prior_seq <- lab_seq[-length(lab_seq)]
        prior_key <- glmtp_label_key(prior_seq)
        D_prior <- D[[prior_key]]
        prev_fit_ids <- if (step_i - 1L >= 1L) {
          fit_ids[[step_i - 1L]][[prior_key]]
        } else {
          NULL
        }

        if (
          is.null(D_prior) ||
            is.null(prev_fit_ids) ||
            length(prev_fit_ids) == 0L
        ) {
          g_k <- rep(0, length(coef_clean(model_k)))
        } else {
          rows_aval <- which(a_now_cur == a_last)
          cur_ids_aval <- cur_ids[rows_aval]

          idx_in_all <- id_to_idx[prev_fit_ids]
          rows_in_aval <- match(prev_fit_ids, cur_ids_aval)
          keep <- !is.na(idx_in_all) & !is.na(rows_in_aval)

          if (any(keep)) {
            D_prev_vals <- D_prior[idx_in_all[keep]]
            w_prev_vals <- if (has_weights) {
              unname(w_at_step[[step_i - 1L]][prev_fit_ids[keep]])
            } else {
              rep(1, sum(keep))
            }

            kept_cur_data <- cur_data[rows_aval][rows_in_aval[keep]]
            d_vals <- policy(
              prior_seq,
              rep(a_last, sum(keep)),
              kept_cur_data,
              step_i,
              n_times
            )
            nd_kept <- data.table::copy(kept_cur_data)
            nd_kept[[treatment]] <- d_vals

            X_star_k <- iv_design_matrix(model_k, nd_kept)
            eta_star_k <- as.numeric(X_star_k %*% coef_clean(model_k))
            mu_eta_k <- model_k$family$mu.eta(eta_star_k)

            g_k <- as.numeric(
              crossprod(X_star_k, D_prev_vals * w_prev_vals * mu_eta_k)
            )
          } else {
            g_k <- rep(0, length(coef_clean(model_k)))
          }
        }
      }

      fit_id_idx <- id_to_idx[fit_ids_k]
      na_act <- model_k$na.action
      if (!is.null(na_act)) {
        fit_id_idx <- fit_id_idx[-na_act]
      }
      res <- correct_model(model_k, g_k, fit_id_idx, n)
      IF_acc <- IF_acc + res$correction
      D[[lab_key]] <- res$d
    }
  }

  # ------------------------------------------------------------------
  # Base step (step n_times, label-independent model)
  # ------------------------------------------------------------------
  m_base <- models[[n_times]]
  if (!is.null(m_base)) {
    fit_ids_base <- fit_ids[[n_times]]
    p_base <- length(coef_clean(m_base))
    g_0 <- rep(0, p_base)

    tau_time <- time_points[n_times]
    mask_tau <- data[[time_col]] == tau_time
    tau_data <- data[mask_tau]
    tau_ids <- as.character(tau_data[[id_col]])

    if (n_times == 1L) {
      # Degenerate single-period case: no per-label steps; gradient comes from
      # target rows directly (q_1_i is the base-model prediction; theta is
      # just the weighted mean of those predictions).
      tau_idx_in_all <- id_to_idx[tau_ids]
      target_in_tau <- target[tau_idx_in_all]
      valid_tau <- !is.na(target_in_tau) & target_in_tau
      if (any(valid_tau)) {
        a_tau_target <- tau_data[[treatment]][valid_tau]
        tau_data_target <- tau_data[valid_tau]
        d_vals_target <- policy(
          numeric(0L),
          a_tau_target,
          tau_data_target,
          1L,
          n_times
        )
        nd_target <- data.table::copy(tau_data_target)
        nd_target[[treatment]] <- d_vals_target
        X_star_target <- iv_design_matrix(m_base, nd_target)
        eta_star_target <- as.numeric(X_star_target %*% coef_clean(m_base))
        mu_eta_target <- m_base$family$mu.eta(eta_star_target)
        w_base <- if (has_weights) w_t else rep(1, sum(valid_tau))
        g_0 <- as.numeric(
          crossprod(X_star_target, w_base * mu_eta_target)
        ) /
          sum_w_target
      }
    } else {
      # General case: accumulate from D[s_{tau-1}] over all prior labels.
      # Each prior label contributes through its sensitivity vector and a
      # label-specific policy-shifted counterfactual frame.
      labels_prior_tau <- glmtp_enumerate_labels(support, n_times - 1L)
      for (ell_seq in labels_prior_tau) {
        ell_key <- glmtp_label_key(ell_seq)
        D_ell <- D[[ell_key]]
        prev_fit_ids_ell <- fit_ids[[n_times - 1L]][[ell_key]]

        if (
          is.null(D_ell) ||
            is.null(prev_fit_ids_ell) ||
            length(prev_fit_ids_ell) == 0L
        ) {
          next
        }

        idx_in_all_ell <- id_to_idx[prev_fit_ids_ell]
        rows_in_tau <- match(prev_fit_ids_ell, tau_ids)
        keep_tau <- !is.na(idx_in_all_ell) & !is.na(rows_in_tau)

        if (!any(keep_tau)) {
          next
        }

        D_ell_vals <- D_ell[idx_in_all_ell[keep_tau]]
        w_prev_ell <- if (has_weights) {
          unname(w_at_step[[n_times - 1L]][prev_fit_ids_ell[keep_tau]])
        } else {
          rep(1, sum(keep_tau))
        }

        kept_tau_data <- tau_data[rows_in_tau[keep_tau]]
        a_tau_kept <- kept_tau_data[[treatment]]
        d_policy_vals <- policy(
          ell_seq,
          a_tau_kept,
          kept_tau_data,
          n_times,
          n_times
        )
        nd_kept_tau <- data.table::copy(kept_tau_data)
        nd_kept_tau[[treatment]] <- d_policy_vals

        X_star_ell <- iv_design_matrix(m_base, nd_kept_tau)
        eta_star_ell <- as.numeric(X_star_ell %*% coef_clean(m_base))
        mu_eta_ell <- m_base$family$mu.eta(eta_star_ell)

        g_0 <- g_0 +
          as.numeric(
            crossprod(X_star_ell, D_ell_vals * w_prev_ell * mu_eta_ell)
          )
      }
    }

    fit_id_idx_base <- id_to_idx[fit_ids_base]
    na_act_base <- m_base$na.action
    if (!is.null(na_act_base)) {
      fit_id_idx_base <- fit_id_idx_base[-na_act_base]
    }
    res_base <- correct_model(m_base, g_0, fit_id_idx_base, n)
    IF_acc <- IF_acc + res_base$correction
  }

  IF_acc
}
