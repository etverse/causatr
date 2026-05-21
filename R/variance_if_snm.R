#' Sandwich variance for point-treatment SNM blip parameters
#'
#' Implements the stacked M-estimation sandwich for the system
#' \eqn{(\hat\alpha, \hat\psi)} where \eqn{\hat\alpha} are the
#' treatment model parameters and \eqn{\hat\psi} are the blip
#' parameters. The bread is block-triangular:
#' \deqn{A = \begin{pmatrix} A_{\alpha\alpha} & 0 \\
#'   A_{\psi\alpha} & A_{\psi\psi} \end{pmatrix}}
#' and the per-individual IF for \eqn{\hat\psi} is
#' \deqn{\mathrm{IF}_{\psi,i} = A_{\psi\psi}^{-1}\bigl(
#'   \omega_{\psi,i} - A_{\psi\alpha} A_{\alpha\alpha}^{-1}
#'   s_{\alpha,i}\bigr),}
#' where \eqn{\omega_{\psi,i}} is the blip estimating equation
#' score and \eqn{s_{\alpha,i}} is the treatment model score.
#'
#' @references
#' Robins JM (1994). Correcting for non-compliance in randomized trials
#' using structural nested mean models. *Communications in Statistics —
#' Theory and Methods*, 23(8), 2379–2412.
#'
#' Vansteelandt S, Joffe M (2014). Structural nested models and G-estimation:
#' a survey. *Statistical Science*, 29(2), 220–238.
#'
#' Stefanski LA, Boos DD (2002). The calculus of M-estimation.
#' *The American Statistician*, 56(1), 29–38.
#'
#' When a treatment-free model is present, \eqn{(\hat\beta, \hat\psi)}
#' are solved jointly from a linear system. The stacked EE system
#' becomes \eqn{(\hat\alpha, \hat\theta)} where
#' \eqn{\hat\theta = (\hat\beta, \hat\psi)} is the joint
#' treatment-free + blip parameter vector. The bread for
#' \eqn{\hat\theta} is the joint system matrix, and the
#' cross-derivative \eqn{A_{\theta,\alpha}} captures how both
#' \eqn{\hat\beta} and \eqn{\hat\psi} depend on
#' \eqn{\hat\alpha} through the treatment residual.
#'
#' @param fit A `causatr_fit` with `estimator = "snm"`.
#' @param snm_result Output of `compute_snm_blip_point()`.
#' @return A named p_psi x p_psi variance-covariance matrix for
#'   \eqn{\hat\psi}.
#' @noRd
variance_if_snm <- function(fit, snm_result) {
  psi_hat <- snm_result$psi_hat
  M <- snm_result$M
  R <- snm_result$R
  A <- snm_result$A
  Y <- snm_result$Y
  n <- snm_result$n_obs
  fit_rows <- snm_result$fit_rows
  p_psi <- length(psi_hat)
  trt_model <- fit$details$treatment_model
  beta_hat <- snm_result$beta_hat
  Z <- snm_result$Z
  has_tf <- !is.null(beta_hat)

  # --- Treatment model IF ---
  # prepare_model_if() returns the design matrix, working residuals, and
  # B_inv = (X'WX)^{-1} (without the 1/n factor). The product
  # X_fit * r_score gives the per-observation score s_{alpha,i}.
  trt_fit_idx <- which(trt_model$fit_rows)
  n_total <- nrow(fit$data)
  trt_prep <- prepare_model_if(trt_model$model, trt_fit_idx, n_total)

  alpha_hat <- trt_model$alpha_hat
  X_prop <- trt_model$X_prop
  trt_family <- trt_model$model$family

  # Per-obs treatment model IF: IF_{alpha,i} = s_{alpha,i} (X'WX)^{-1}
  # Note: B_inv carries no 1/n factor, so IF_alpha_fit is scaled by n
  # relative to the standard 1/n convention. This is corrected below.
  trt_score <- trt_prep$X_fit * trt_prep$r_score
  IF_alpha_fit <- trt_score %*% trt_prep$B_inv

  n_trt_fit <- nrow(IF_alpha_fit)
  if (n_trt_fit != n) {
    rlang::abort(
      paste0(
        "Treatment model fit rows (",
        n_trt_fit,
        ") != SNM fit rows (",
        n,
        ")."
      ),
      class = "causatr_snm_row_mismatch"
    )
  }

  if (has_tf) {
    # --- Joint (beta, psi) system ---
    # The joint EE is:
    #   phi_beta(alpha, beta, psi) = (1/n) sum Z_i (Y_i - Z_i beta - A_i M_i psi)
    #   phi_psi(alpha, beta, psi)  = (1/n) sum R_i(alpha) m_i (Y_i - Z_i beta - A_i M_i psi)
    # The joint bread A_{theta,theta} (negative Jacobian w.r.t. theta):
    #   d phi_beta/d beta = -(1/n) Z'Z
    #   d phi_beta/d psi  = -(1/n) Z'(AM)
    #   d phi_psi/d beta  = -(1/n) (RM)'Z
    #   d phi_psi/d psi   = -(1/n) (RM)'(AM)
    p_beta <- length(beta_hat)
    AM <- A * M
    RM <- R * M

    # A_{theta,theta} = -(1/n) [Z'Z, Z'(AM); (RM)'Z, (RM)'(AM)]
    # This is the negative Jacobian of the joint EE w.r.t. theta = (beta, psi)
    A_theta_theta <- -rbind(
      cbind(crossprod(Z, Z), crossprod(Z, AM)),
      cbind(crossprod(RM, Z), crossprod(RM, AM))
    ) /
      n

    A_theta_theta_inv <- tryCatch(
      solve(A_theta_theta),
      error = function(e) {
        rlang::abort(
          "Joint treatment-free + blip bread matrix is singular.",
          class = "causatr_snm_singular",
          parent = e
        )
      }
    )

    # Joint score: omega_theta_i = (phi_beta_i, phi_psi_i)
    gamma_hat <- A * as.numeric(M %*% psi_hat)
    tf_pred <- as.numeric(Z %*% beta_hat)
    resid_joint <- Y - tf_pred - gamma_hat
    # phi_beta_i = Z_i * resid_joint_i
    omega_beta <- Z * resid_joint
    # phi_psi_i = R_i * M_i * resid_joint_i
    omega_psi <- M * (R * resid_joint)
    omega_theta <- cbind(omega_beta, omega_psi)

    # Cross-derivative A_{theta,alpha}: how the joint EE depends on alpha
    # Only phi_psi depends on alpha (through R_i = A_i - mu(X_i alpha));
    # phi_beta does not depend on alpha. So:
    #   A_{theta,alpha} = [0_{p_beta x p_alpha}; A_{psi,alpha}]
    # Cross-derivative: d/dalpha of the averaged blip EE evaluated at
    # (alpha_hat, beta_hat, psi_hat). Only the psi block depends on alpha
    # through R_i(alpha) = A_i - mu(X_i alpha); phi_beta does not.
    phi_bar_alpha <- function(alpha) {
      eta <- as.numeric(X_prop %*% alpha)
      mu <- trt_family$linkinv(eta)
      R_alpha <- A - mu
      # \bar\phi_\psi(\alpha) = (1/n) \sum R_i(\alpha) m_i resid_i
      colMeans(M * (R_alpha * resid_joint))
    }
    A_psi_alpha <- -numDeriv::jacobian(phi_bar_alpha, x = alpha_hat)

    # A_{theta,alpha} = [0_{p_beta x p_alpha}; A_{psi,alpha}]
    A_theta_alpha <- rbind(
      matrix(0, nrow = p_beta, ncol = length(alpha_hat)),
      A_psi_alpha
    )

    # Treatment model correction for the joint system
    # IF_alpha_fit omits the n factor from A_aa^{-1} = -n(X'WX)^{-1}
    correction_theta <- n * IF_alpha_fit %*% t(A_theta_alpha)

    # Total IF for theta = (beta, psi): (n x (p_beta + p_psi))
    adjusted_score_theta <- omega_theta - correction_theta
    IF_theta <- adjusted_score_theta %*% t(A_theta_theta_inv)

    # Extract the psi block: columns (p_beta+1)...(p_beta+p_psi)
    # from the full (n x (p_beta + p_psi)) IF matrix
    IF_psi <- IF_theta[, p_beta + seq_len(p_psi), drop = FALSE]
  } else {
    # --- Standard (psi only) system ---
    # Bread: A_{\psi\psi} = -(1/n) \sum m_i m_i^T R_i A_i
    RA <- R * A
    A_psi_psi <- -crossprod(M, M * RA) / n

    A_psi_psi_inv <- tryCatch(
      solve(A_psi_psi),
      error = function(e) {
        rlang::abort(
          "SNM blip bread matrix is singular.",
          class = "causatr_snm_singular",
          parent = e
        )
      }
    )

    # Blip score: omega_{psi,i} = R_i * m_i * (Y_i - gamma(A_i, L_i; psi))
    gamma_hat <- A * as.numeric(M %*% psi_hat)
    resid_blip <- Y - gamma_hat
    omega_psi <- M * (R * resid_blip)

    # Cross-derivative A_{psi,alpha}: d/dalpha of the averaged g-estimating
    # equation. Computed numerically via the closure below (same pattern
    # as the joint branch above and the IPW variance engine).
    phi_bar_alpha <- function(alpha) {
      eta <- as.numeric(X_prop %*% alpha)
      mu <- trt_family$linkinv(eta)
      R_alpha <- A - mu
      # \bar\phi_\psi(\alpha) = (1/n) \sum R_i(\alpha) m_i resid_i
      colMeans(M * (R_alpha * resid_blip))
    }
    A_psi_alpha <- -numDeriv::jacobian(phi_bar_alpha, x = alpha_hat)

    # Treatment model correction
    # IF_alpha_fit omits the n factor from A_aa^{-1} = -n(X'WX)^{-1}
    correction_trt <- n * IF_alpha_fit %*% t(A_psi_alpha)

    # Total IF for psi: (n x p_psi)
    adjusted_score <- omega_psi - correction_trt
    IF_psi <- adjusted_score %*% t(A_psi_psi_inv)
  }

  # Sandwich variance: V_psi = (1/n^2) \sum IF_{psi,i} IF_{psi,i}^T
  vcov_psi <- crossprod(IF_psi) / n^2
  rownames(vcov_psi) <- names(psi_hat)
  colnames(vcov_psi) <- names(psi_hat)
  vcov_psi
}


#' Cluster-aggregated sandwich variance for longitudinal SNM blip parameters
#'
#' Extends the point-treatment sandwich to the backward-sequential
#' estimation case. Each stage \eqn{k} has its own blip parameters
#' \eqn{\psi_k} and treatment-model parameters \eqn{\alpha_k}. The
#' full parameter vector is \eqn{(\alpha_K, \psi_K, \alpha_{K-1},
#' \psi_{K-1}, \ldots, \alpha_0, \psi_0)}, estimated backward.
#'
#' The key complication: \eqn{\hat\psi_{k-1}} depends on
#' \eqn{\hat\psi_k} through the backward-transformed outcome
#' \eqn{H_{k-1} = H_k - \gamma_k(\hat\psi_k)}. This creates cross-
#' stage derivatives: the stage-\eqn{(k-1)} bread row has off-diagonal
#' blocks \eqn{\partial\phi_{k-1}/\partial\psi_j} for all \eqn{j > k-1}.
#'
#' Per-individual influence functions are accumulated across stages
#' and the cluster-robust sandwich (one cluster = one individual)
#' gives the final \eqn{V(\hat\psi)}.
#'
#' When a treatment-free model is present, delegates to the joint
#' \eqn{(\beta_k, \psi_k)} system in
#' \code{variance_if_snm_longitudinal_tf()}; otherwise uses the
#' psi-only system in \code{variance_if_snm_longitudinal_notf()}.
#'
#' @references
#' Robins JM (1994). Correcting for non-compliance in randomized trials
#' using structural nested mean models. *Communications in Statistics —
#' Theory and Methods*, 23(8), 2379–2412.
#'
#' @param fit A `causatr_fit` with `estimator = "snm"` and
#'   `type = "longitudinal"`.
#' @param snm_result Output of `compute_snm_blip_longitudinal()`.
#' @return A named \eqn{p_{total} \times p_{total}} variance-covariance
#'   matrix for the concatenated \eqn{\hat\psi} vector (all stages).
#' @noRd
variance_if_snm_longitudinal <- function(fit, snm_result) {
  psi_hat <- snm_result$psi_hat
  psi_by_stage <- snm_result$psi_by_stage
  per_period <- snm_result$per_period
  H_by_stage <- snm_result$H_by_stage
  n_id <- snm_result$n_id
  ids <- snm_result$ids
  na_mask <- snm_result$na_mask
  Y_final <- snm_result$Y_final
  id_to_idx <- snm_result$id_to_idx
  n_times <- snm_result$n_times
  p_psi_per_stage <- snm_result$p_psi_per_stage
  beta_by_stage <- snm_result$beta_by_stage
  Z_list <- snm_result$Z_list
  has_tf <- !is.null(beta_by_stage)
  details <- fit$details
  time_points <- details$time_points

  p_total_psi <- n_times * p_psi_per_stage

  # Outcome vector aligned to individual index
  Y_vec <- rep(NA_real_, n_id)
  Y_vec[id_to_idx[ids[na_mask]]] <- Y_final[na_mask]

  if (has_tf) {
    variance_if_snm_longitudinal_tf(
      psi_hat, psi_by_stage, beta_by_stage, Z_list,
      per_period, H_by_stage,
      n_id, ids, na_mask, id_to_idx,
      n_times, p_psi_per_stage, p_total_psi,
      fit
    )
  } else {
    variance_if_snm_longitudinal_notf(
      psi_hat, psi_by_stage,
      per_period, H_by_stage,
      n_id, ids, na_mask, id_to_idx,
      n_times, p_psi_per_stage, p_total_psi,
      fit
    )
  }
}


#' Longitudinal SNM sandwich: psi-only system (no treatment-free model)
#'
#' @param psi_hat Named numeric vector of all per-stage blip parameters
#'   (concatenated: stage 0, ..., stage K-1).
#' @param psi_by_stage List of per-stage psi vectors (index 1 = stage 0).
#' @param per_period List of per-period data (A, R, M, ids, data,
#'   trt_model) from `compute_snm_blip_longitudinal()`.
#' @param H_by_stage List of backward-transformed outcome vectors,
#'   one per stage, aligned to individual index.
#' @param n_id Integer number of unique individuals.
#' @param ids Character vector of unique individual IDs.
#' @param na_mask Logical vector: `TRUE` where the final-period outcome
#'   is non-missing.
#' @param id_to_idx Named integer vector mapping ID strings to row
#'   indices in the individual-level matrices.
#' @param n_times Integer number of time points.
#' @param p_psi_per_stage Integer number of blip parameters per stage.
#' @param p_total_psi Integer total number of psi parameters across
#'   all stages (`n_times * p_psi_per_stage`).
#' @param fit The original `causatr_fit` object (needed for treatment
#'   model access via `prepare_model_if()`).
#' @noRd
variance_if_snm_longitudinal_notf <- function(
  psi_hat,
  psi_by_stage,
  per_period,
  H_by_stage,
  n_id,
  ids,
  na_mask,
  id_to_idx,
  n_times,
  p_psi_per_stage,
  p_total_psi,
  fit
) {
  # Per-stage EE:
  #   \phi_k(\psi_k) = (1/n) \sum_i R_{k,i} M_{k,i} (H_k - A_{k,i} M_{k,i} \psi_k) = 0
  # Bread: A_{\psi_k,\psi_k} = -(1/n) \sum_i R_{k,i} A_{k,i} M_{k,i} M_{k,i}'
  # Cross-stage: dH_k/d\psi_j = -A_j M_j for j > k, so
  #   A_{\psi_k,\psi_j} = (1/n) (R_k M_k)' (A_j M_j)
  A_bread <- matrix(0, p_total_psi, p_total_psi)
  omega_psi <- matrix(0, n_id, p_total_psi)

  stage_idx <- function(k) {
    start <- (k - 1L) * p_psi_per_stage + 1L
    start:(start + p_psi_per_stage - 1L)
  }

  AM_by_stage <- vector("list", n_times)
  for (k in seq_len(n_times)) {
    pp <- per_period[[k]]
    AM_k <- matrix(0, n_id, p_psi_per_stage)
    idx_k <- id_to_idx[pp$ids]
    AM_k[idx_k, ] <- pp$A * pp$M
    AM_by_stage[[k]] <- AM_k
  }

  for (k in seq_len(n_times)) {
    pp <- per_period[[k]]
    idx_k <- id_to_idx[pp$ids]
    valid <- na_mask[match(pp$ids, ids)]
    n_valid <- sum(valid)
    if (n_valid == 0L) {
      next
    }

    psi_k_raw <- psi_by_stage[[k]][seq_len(p_psi_per_stage)]

    M_k <- pp$M[valid, , drop = FALSE]
    R_k <- pp$R[valid]
    A_k <- pp$A[valid]
    H_k <- H_by_stage[[k]][idx_k[valid]]

    gamma_k <- A_k * as.numeric(M_k %*% psi_k_raw)
    resid_k <- H_k - gamma_k

    omega_k <- M_k * (R_k * resid_k)
    idx_valid <- idx_k[valid]
    for (j in seq_len(n_valid)) {
      omega_psi[idx_valid[j], stage_idx(k)] <-
        omega_psi[idx_valid[j], stage_idx(k)] + omega_k[j, ]
    }

    RA_k <- R_k * A_k
    A_bread[stage_idx(k), stage_idx(k)] <-
      -crossprod(M_k, M_k * RA_k) / n_id

    if (k < n_times) {
      RM_k <- M_k * R_k
      for (j in (k + 1L):n_times) {
        AM_j_valid <- AM_by_stage[[j]][idx_k[valid], , drop = FALSE]
        A_bread[stage_idx(k), stage_idx(j)] <-
          crossprod(RM_k, AM_j_valid) / n_id
      }
    }
  }

  A_bread_inv <- tryCatch(
    solve(A_bread),
    error = function(e) {
      rlang::abort(
        "Longitudinal SNM bread matrix is singular.",
        class = "causatr_snm_singular",
        parent = e
      )
    }
  )

  # Treatment model corrections: d\phi_k/d\alpha_k through R_k
  IF_correction <- matrix(0, n_id, p_total_psi)

  for (k in seq_len(n_times)) {
    pp <- per_period[[k]]
    tm_k <- pp$trt_model
    alpha_k <- tm_k$alpha_hat
    X_prop_k <- tm_k$X_prop
    trt_family_k <- tm_k$model$family
    n_k <- length(pp$A)

    trt_fit_idx_k <- which(tm_k$fit_rows)
    trt_prep_k <- prepare_model_if(tm_k$model, trt_fit_idx_k, n_k)
    trt_score_k <- trt_prep_k$X_fit * trt_prep_k$r_score
    IF_alpha_k <- trt_score_k %*% trt_prep_k$B_inv

    valid_full <- na_mask[match(pp$ids, ids)]
    psi_k_raw <- psi_by_stage[[k]][seq_len(p_psi_per_stage)]
    H_k_full <- H_by_stage[[k]][id_to_idx[pp$ids]]
    gamma_k_full <- pp$A * as.numeric(pp$M %*% psi_k_raw)
    resid_k_full <- H_k_full - gamma_k_full

    phi_psi_k_alpha <- function(alpha) {
      eta <- as.numeric(X_prop_k %*% alpha)
      mu <- trt_family_k$linkinv(eta)
      R_alpha <- pp$A - mu
      vals <- pp$M * (R_alpha * resid_k_full)
      vals[!valid_full, ] <- 0
      colMeans(vals)
    }
    A_psi_k_alpha_k <- -numDeriv::jacobian(phi_psi_k_alpha, x = alpha_k)

    correction_k <- n_k * IF_alpha_k %*% t(A_psi_k_alpha_k)

    idx_k <- id_to_idx[pp$ids]
    for (j in seq_len(n_k)) {
      IF_correction[idx_k[j], stage_idx(k)] <-
        IF_correction[idx_k[j], stage_idx(k)] + correction_k[j, ]
    }
  }

  IF_psi <- (omega_psi - IF_correction) %*% t(A_bread_inv)

  vcov_psi <- crossprod(IF_psi) / n_id^2
  rownames(vcov_psi) <- names(psi_hat)
  colnames(vcov_psi) <- names(psi_hat)
  vcov_psi
}


#' Longitudinal SNM sandwich: joint (beta, psi) system with treatment-free model
#'
#' Extends the psi-only case to the joint \eqn{(\beta_k, \psi_k)} system
#' at each stage. The per-stage EE is:
#' \deqn{\phi_{\beta,k} = (1/n) \sum Z_{k,i} (H_k - Z_k \beta_k - A_k M_k \psi_k)}
#' \deqn{\phi_{\psi,k} = (1/n) \sum R_{k,i} M_{k,i} (H_k - Z_k \beta_k - A_k M_k \psi_k)}
#' The \eqn{\beta} parameters are nuisance; the final vcov is for \eqn{\psi}
#' only, obtained by marginalizing over \eqn{\beta} through the joint system.
#'
#' @param psi_hat Named numeric vector of all per-stage blip parameters
#'   (concatenated: stage 0, ..., stage K-1).
#' @param psi_by_stage List of per-stage psi vectors (index 1 = stage 0).
#' @param beta_by_stage List of per-stage treatment-free parameter vectors.
#' @param Z_list List of per-period treatment-free design matrices
#'   (full data, not validity-filtered).
#' @param per_period List of per-period data (A, R, M, ids, data,
#'   trt_model) from `compute_snm_blip_longitudinal()`.
#' @param H_by_stage List of backward-transformed outcome vectors,
#'   one per stage, aligned to individual index.
#' @param n_id Integer number of unique individuals.
#' @param ids Character vector of unique individual IDs.
#' @param na_mask Logical vector: `TRUE` where the final-period outcome
#'   is non-missing.
#' @param id_to_idx Named integer vector mapping ID strings to row
#'   indices in the individual-level matrices.
#' @param n_times Integer number of time points.
#' @param p_psi_per_stage Integer number of blip parameters per stage.
#' @param p_total_psi Integer total number of psi parameters across
#'   all stages (`n_times * p_psi_per_stage`).
#' @param fit The original `causatr_fit` object (needed for treatment
#'   model access via `prepare_model_if()`).
#' @noRd
variance_if_snm_longitudinal_tf <- function(
  psi_hat,
  psi_by_stage,
  beta_by_stage,
  Z_list,
  per_period,
  H_by_stage,
  n_id,
  ids,
  na_mask,
  id_to_idx,
  n_times,
  p_psi_per_stage,
  p_total_psi,
  fit
) {
  # Per-stage theta_k = (beta_k, psi_k). We work with the full theta
  # vector = [theta_0, theta_1, ..., theta_{K-1}] and extract the psi
  # marginal at the end.
  p_beta <- length(beta_by_stage[[1L]])
  p_theta_per_stage <- p_beta + p_psi_per_stage
  p_total_theta <- n_times * p_theta_per_stage

  # Index helpers for per-stage blocks in the theta vector
  theta_idx <- function(k) {
    start <- (k - 1L) * p_theta_per_stage + 1L
    start:(start + p_theta_per_stage - 1L)
  }
  # Indices of psi within a per-stage theta block
  psi_in_theta <- (p_beta + 1L):p_theta_per_stage

  A_bread <- matrix(0, p_total_theta, p_total_theta)
  omega_theta <- matrix(0, n_id, p_total_theta)

  # Precompute per-individual A_j * M_j for cross-stage psi derivatives
  AM_by_stage <- vector("list", n_times)
  for (k in seq_len(n_times)) {
    pp <- per_period[[k]]
    AM_k <- matrix(0, n_id, p_psi_per_stage)
    idx_k <- id_to_idx[pp$ids]
    AM_k[idx_k, ] <- pp$A * pp$M
    AM_by_stage[[k]] <- AM_k
  }

  for (k in seq_len(n_times)) {
    pp <- per_period[[k]]
    idx_k <- id_to_idx[pp$ids]
    valid <- na_mask[match(pp$ids, ids)]
    n_valid <- sum(valid)
    if (n_valid == 0L) {
      next
    }

    psi_k_raw <- psi_by_stage[[k]][seq_len(p_psi_per_stage)]
    beta_k <- beta_by_stage[[k]]

    M_k <- pp$M[valid, , drop = FALSE]
    R_k <- pp$R[valid]
    A_k <- pp$A[valid]
    H_k <- H_by_stage[[k]][idx_k[valid]]
    Z_k <- Z_list[[k]][valid, , drop = FALSE]

    # Joint residual: H_k - Z_k beta_k - A_k M_k psi_k
    gamma_k <- A_k * as.numeric(M_k %*% psi_k_raw)
    tf_pred_k <- as.numeric(Z_k %*% beta_k)
    resid_k <- H_k - tf_pred_k - gamma_k

    AM_k_mat <- A_k * M_k
    RM_k <- R_k * M_k

    # Per-obs scores: (phi_beta, phi_psi)
    omega_beta_k <- Z_k * resid_k
    omega_psi_k <- RM_k * resid_k
    omega_k <- cbind(omega_beta_k, omega_psi_k)

    idx_valid <- idx_k[valid]
    for (j in seq_len(n_valid)) {
      omega_theta[idx_valid[j], theta_idx(k)] <-
        omega_theta[idx_valid[j], theta_idx(k)] + omega_k[j, ]
    }

    # Diagonal bread block: A_{theta_k, theta_k}
    # Same structure as the point-treatment joint bread:
    #   [Z'Z, Z'(AM); (RM)'Z, (RM)'(AM)] / n_id
    A_bread[theta_idx(k), theta_idx(k)] <- -rbind(
      cbind(crossprod(Z_k, Z_k), crossprod(Z_k, AM_k_mat)),
      cbind(crossprod(RM_k, Z_k), crossprod(RM_k, AM_k_mat))
    ) / n_id

    # Cross-stage: d theta_k / d psi_j for j > k (through H_k).
    # H_k = Y - sum_{l>k} gamma_l, dH_k/dpsi_j = -A_j M_j.
    # d phi_{beta,k}/d psi_j = -(1/n) Z_k' (A_j M_j)
    # d phi_{psi,k}/d psi_j = -(1/n) (R_k M_k)' (A_j M_j)
    # Convention: A = -d phi/d theta, so negate:
    # A_{theta_k, psi_j} = (1/n) [Z_k'; (R_k M_k)']' (A_j M_j)
    if (k < n_times) {
      for (j in (k + 1L):n_times) {
        AM_j_valid <- AM_by_stage[[j]][idx_k[valid], , drop = FALSE]
        # Place into the psi_j columns of stage j's theta block
        cross_block <- rbind(
          crossprod(Z_k, AM_j_valid),
          crossprod(RM_k, AM_j_valid)
        ) / n_id
        # Map: rows = theta_k, cols = psi_j within theta_j
        j_psi_cols <- theta_idx(j)[psi_in_theta]
        A_bread[theta_idx(k), j_psi_cols] <- cross_block
      }
    }
  }

  A_bread_inv <- tryCatch(
    solve(A_bread),
    error = function(e) {
      rlang::abort(
        "Longitudinal SNM joint bread matrix is singular.",
        class = "causatr_snm_singular",
        parent = e
      )
    }
  )

  # Treatment model corrections: only phi_{psi,k} depends on alpha_k
  # through R_k. phi_{beta,k} does not depend on alpha_k.
  IF_correction <- matrix(0, n_id, p_total_theta)

  for (k in seq_len(n_times)) {
    pp <- per_period[[k]]
    tm_k <- pp$trt_model
    alpha_k <- tm_k$alpha_hat
    X_prop_k <- tm_k$X_prop
    trt_family_k <- tm_k$model$family
    n_k <- length(pp$A)

    trt_fit_idx_k <- which(tm_k$fit_rows)
    trt_prep_k <- prepare_model_if(tm_k$model, trt_fit_idx_k, n_k)
    trt_score_k <- trt_prep_k$X_fit * trt_prep_k$r_score
    IF_alpha_k <- trt_score_k %*% trt_prep_k$B_inv

    valid_full <- na_mask[match(pp$ids, ids)]
    psi_k_raw <- psi_by_stage[[k]][seq_len(p_psi_per_stage)]
    beta_k <- beta_by_stage[[k]]
    H_k_full <- H_by_stage[[k]][id_to_idx[pp$ids]]
    gamma_k_full <- pp$A * as.numeric(pp$M %*% psi_k_raw)
    Z_k_full <- Z_list[[k]]
    tf_pred_k_full <- as.numeric(Z_k_full %*% beta_k)
    resid_k_full <- H_k_full - tf_pred_k_full - gamma_k_full

    # d phi_{psi,k}/d alpha_k: only the psi block depends on alpha_k
    phi_psi_k_alpha <- function(alpha) {
      eta <- as.numeric(X_prop_k %*% alpha)
      mu <- trt_family_k$linkinv(eta)
      R_alpha <- pp$A - mu
      vals <- pp$M * (R_alpha * resid_k_full)
      vals[!valid_full, ] <- 0
      colMeans(vals)
    }
    A_psi_k_alpha_k <- -numDeriv::jacobian(phi_psi_k_alpha, x = alpha_k)

    # A_{theta_k, alpha_k} = [0_{p_beta x p_alpha}; A_{psi_k, alpha_k}]
    A_theta_k_alpha_k <- rbind(
      matrix(0, nrow = p_beta, ncol = length(alpha_k)),
      A_psi_k_alpha_k
    )

    correction_k <- n_k * IF_alpha_k %*% t(A_theta_k_alpha_k)

    idx_k <- id_to_idx[pp$ids]
    for (j in seq_len(n_k)) {
      IF_correction[idx_k[j], theta_idx(k)] <-
        IF_correction[idx_k[j], theta_idx(k)] + correction_k[j, ]
    }
  }

  IF_theta <- (omega_theta - IF_correction) %*% t(A_bread_inv)

  # Extract the psi marginal: columns corresponding to psi within each
  # stage's theta block
  psi_cols <- unlist(lapply(seq_len(n_times), function(k) {
    theta_idx(k)[psi_in_theta]
  }))
  IF_psi <- IF_theta[, psi_cols, drop = FALSE]

  vcov_psi <- crossprod(IF_psi) / n_id^2
  rownames(vcov_psi) <- names(psi_hat)
  colnames(vcov_psi) <- names(psi_hat)
  vcov_psi
}
