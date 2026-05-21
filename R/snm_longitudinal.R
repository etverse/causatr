#' Solve the longitudinal g-estimating equation for per-stage blip parameters
#'
#' Implements backward sequential g-estimation (Robins 1994): starting at
#' the final stage \eqn{K}, solve for \eqn{\hat\psi_K} using the observed
#' outcome \eqn{Y}, then compute the backward-transformed outcome
#' \deqn{H_{K-1} = Y - \gamma_K(A_K, \bar{L}_K; \hat\psi_K)}
#' and solve for \eqn{\hat\psi_{K-1}} using \eqn{H_{K-1}} as the
#' pseudo-outcome, continuing backward to stage 0.
#'
#' Each stage \eqn{k} has its own blip parameters \eqn{\psi_k} estimated
#' from the per-stage g-estimating equation
#' \deqn{\sum_i R_{k,i} \cdot M_{k,i} \cdot (H_k - \gamma_k(A_{k,i},
#'   \bar{L}_{k,i}; \psi_k)) = 0,}
#' where \eqn{R_{k,i} = A_{k,i} - \hat{E}[A_k \mid \bar{A}_{k-1},
#' \bar{L}_k]} and \eqn{H_k = Y - \sum_{j > k} \gamma_j(\hat\psi_j)}.
#' For linear blip this is a closed-form solve at each stage. This
#' matches the DTRreg backward-induction approach.
#'
#' When a treatment-free model is present, \eqn{(\beta_k, \psi_k)}
#' are solved jointly at each stage from the combined system.
#'
#' @param fit A `causatr_fit` with `estimator = "snm"` and
#'   `type = "longitudinal"`.
#' @return A list with:
#'   \describe{
#'     \item{psi_hat}{Named numeric vector of all per-stage blip
#'       parameters concatenated (stage K first, then K-1, ..., 0),
#'       with names prefixed by `"stageK_"`.}
#'     \item{psi_by_stage}{List of per-stage psi vectors (index 1 =
#'       stage 0, ..., K = final stage).}
#'     \item{per_period}{List of per-period data for the variance engine.}
#'     \item{H_by_stage}{List of backward-transformed outcomes per stage.}
#'     \item{n_id}{Number of unique individuals.}
#'     \item{n_obs}{Alias for n_id (for contrast dispatch).}
#'     \item{ids}{Character vector of unique individual IDs.}
#'     \item{beta_by_stage}{Per-stage treatment-free parameters or `NULL`.}
#'     \item{Z_list}{Per-period treatment-free design matrices or `NULL`.}
#'   }
#' @noRd
compute_snm_blip_longitudinal <- function(fit) {
  details <- fit$details
  blip_spec <- details$blip_spec
  treatment_models <- details$treatment_models_by_time
  fit_data_by_time <- details$fit_data_by_time
  time_points <- details$time_points
  n_times <- details$n_times
  id_col <- fit$id
  time_col <- fit$time
  outcome <- fit$outcome
  treatment <- fit$treatment
  tf_formula <- details$treatment_free
  data <- fit$data
  p_psi_per_stage <- blip_spec$n_params

  # Final-period rows carry the outcome
  last_time <- time_points[n_times]
  rows_final <- data[[time_col]] == last_time
  data_final <- data[rows_final]
  ids_all <- as.character(data_final[[id_col]])
  n_id <- length(ids_all)
  id_to_idx <- stats::setNames(seq_len(n_id), ids_all)

  Y_final <- data_final[[outcome]]
  na_mask <- !is.na(Y_final)

  # Per-period data extraction: treatment residual R_k, blip design M_k
  per_period <- vector("list", n_times)
  for (k in seq_len(n_times)) {
    tp <- as.character(time_points[k])
    data_k <- fit_data_by_time[[tp]]
    tm_k <- treatment_models[[tp]]
    ids_k <- as.character(data_k[[id_col]])

    A_k <- data_k[[treatment]]
    e_k <- stats::predict(tm_k$model, type = "response")
    R_k <- A_k - e_k

    M_k <- build_blip_design_matrix(data_k, blip_spec)

    per_period[[k]] <- list(
      A = A_k,
      R = R_k,
      M = M_k,
      ids = ids_k,
      data = data_k,
      trt_model = tm_k
    )
  }

  has_tf <- !is.null(tf_formula)

  # Outcome vector aligned to individual index
  Y_vec <- rep(NA_real_, n_id)
  Y_vec[id_to_idx[ids_all[na_mask]]] <- Y_final[na_mask]

  # --- Backward sequential estimation ---
  # H_K = Y (the final outcome), then H_{k-1} = H_k - gamma_k(psi_k).
  # At each stage k, solve the point-treatment g-estimating equation
  # using H_k as the pseudo-outcome.
  H_current <- Y_vec
  psi_by_stage <- vector("list", n_times)
  beta_by_stage <- if (has_tf) vector("list", n_times) else NULL
  H_by_stage <- vector("list", n_times)
  Z_list <- if (has_tf) vector("list", n_times) else NULL

  # Iterate backward: k = K, K-1, ..., 1 (1-indexed: n_times, ..., 1)
  for (k in n_times:1L) {
    pp <- per_period[[k]]
    idx_k <- id_to_idx[pp$ids]
    valid <- na_mask[match(pp$ids, ids_all)]
    n_valid <- sum(valid)

    M_k <- pp$M[valid, , drop = FALSE]
    R_k <- pp$R[valid]
    A_k <- pp$A[valid]
    H_k <- H_current[idx_k[valid]]

    H_by_stage[[k]] <- H_current

    if (has_tf) {
      # Joint (beta_k, psi_k) estimation at this stage.
      # System: Z_k'(H_k - Z_k beta_k - A_k M_k psi_k) = 0
      #         (R_k M_k)'(H_k - Z_k beta_k - A_k M_k psi_k) = 0
      Z_k <- stats::model.matrix(tf_formula, data = pp$data[valid, ])
      Z_list[[k]] <- stats::model.matrix(tf_formula, data = pp$data)
      p_beta <- ncol(Z_k)

      AM_k <- A_k * M_k
      RM_k <- R_k * M_k

      lhs <- rbind(
        cbind(crossprod(Z_k, Z_k), crossprod(Z_k, AM_k)),
        cbind(crossprod(RM_k, Z_k), crossprod(RM_k, AM_k))
      )
      rhs <- c(
        as.numeric(crossprod(Z_k, H_k)),
        as.numeric(crossprod(RM_k, H_k))
      )

      theta_k <- tryCatch(
        as.numeric(solve(lhs, rhs)),
        error = function(e) {
          rlang::abort(
            c(
              paste0(
                "Stage-",
                k - 1L,
                " joint treatment-free + blip system is singular."
              ),
              i = "Check treatment-free formula and blip specification."
            ),
            class = "causatr_snm_singular",
            parent = e
          )
        }
      )

      beta_k <- theta_k[seq_len(p_beta)]
      names(beta_k) <- colnames(Z_k)
      psi_k <- theta_k[p_beta + seq_len(p_psi_per_stage)]
      beta_by_stage[[k]] <- beta_k
    } else {
      # Standard per-stage solve: (M' diag(R*A) M)^{-1} M' diag(R) H_k
      RA_k <- R_k * A_k
      lhs <- crossprod(M_k, M_k * RA_k)
      rhs <- crossprod(M_k, R_k * H_k)

      psi_k <- tryCatch(
        as.numeric(solve(lhs, rhs)),
        error = function(e) {
          rlang::abort(
            c(
              paste0(
                "Stage-",
                k - 1L,
                " g-estimating equation is singular."
              ),
              i = "Check blip specification and treatment model."
            ),
            class = "causatr_snm_singular",
            parent = e
          )
        }
      )
    }

    stage_label <- paste0("stage", k - 1L, "_")
    names(psi_k) <- paste0(stage_label, blip_spec$param_names)
    psi_by_stage[[k]] <- psi_k

    # Update H for the next (earlier) stage:
    # H_{k-1} = H_k - gamma_k(A_k, M_k; psi_k)
    # gamma_k = A_k * (M_k %*% psi_k) for each individual at stage k
    gamma_k_vec <- rep(0, n_id)
    gamma_k_vals <- pp$A * as.numeric(pp$M %*% psi_k[seq_len(p_psi_per_stage)])
    gamma_k_vec[id_to_idx[pp$ids]] <- gamma_k_vals
    H_current <- H_current - gamma_k_vec
  }

  # Concatenate per-stage psi into a single named vector.
  # Order: stage 0 first, then stage 1, ..., stage K-1
  psi_hat <- unlist(psi_by_stage)

  list(
    psi_hat = psi_hat,
    psi_by_stage = psi_by_stage,
    per_period = per_period,
    H_by_stage = H_by_stage,
    n_id = n_id,
    n_obs = n_id,
    ids = ids_all,
    na_mask = na_mask,
    Y_final = Y_final,
    id_to_idx = id_to_idx,
    beta_by_stage = beta_by_stage,
    Z_list = Z_list,
    n_times = n_times,
    p_psi_per_stage = p_psi_per_stage,
    blip_spec = blip_spec
  )
}
