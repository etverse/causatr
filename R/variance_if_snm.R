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
