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
#' @param fit A `causatr_fit` with `estimator = "snm"`.
#' @param snm_result Output of [compute_snm_blip_point()].
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

  # Blip bread: A_{psi,psi} = -(1/n) sum_i m_i m_i^T R_i A_i
  # This is the negative derivative of the g-estimating equation
  # w.r.t. psi, evaluated at psi_hat.
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
  # where gamma(A_i, L_i; psi) = A_i * (M_i %*% psi)
  gamma_hat <- A * as.numeric(M %*% psi_hat)
  resid_blip <- Y - gamma_hat
  # omega_psi is n x p_psi: each row is R_i * m_i * resid_i
  omega_psi <- M * (R * resid_blip)

  # Treatment model IF: s_{alpha,i} via prepare_model_if()
  trt_fit_idx <- which(trt_model$fit_rows)
  n_total <- nrow(fit$data)
  trt_prep <- prepare_model_if(trt_model$model, trt_fit_idx, n_total)

  # Cross-derivative A_{psi,alpha} via numDeriv::jacobian().
  # The g-estimating equation as a function of alpha:
  #   phi(alpha) = (1/n) sum_i R_i(alpha) * m_i * (Y_i - A_i * M_i psi)
  # where R_i(alpha) = A_i - mu(X_i alpha).
  alpha_hat <- trt_model$alpha_hat
  X_prop <- trt_model$X_prop
  trt_family <- trt_model$model$family

  phi_bar <- function(alpha) {
    eta <- as.numeric(X_prop %*% alpha)
    mu <- trt_family$linkinv(eta)
    R_alpha <- A - mu
    colMeans(M * (R_alpha * resid_blip))
  }

  # A_{psi,alpha} = -d phi_bar / d alpha (p_psi x p_alpha)
  A_psi_alpha <- -numDeriv::jacobian(phi_bar, x = alpha_hat)

  # Per-obs treatment model IF: s_{alpha,i} = X_i * r^{score}_i,
  # projected through the bread to give
  # IF_{alpha,i} = (X'WX)^{-1} s_{alpha,i}.
  trt_score <- trt_prep$X_fit * trt_prep$r_score
  IF_alpha_fit <- trt_score %*% trt_prep$B_inv

  # For point treatment, the treatment model and blip estimation
  # use the same fit rows. Assert alignment.
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

  # Correction: A_{psi,alpha} A_{alpha,alpha}^{-1} s_{alpha,i}.
  # IF_alpha_fit = s_{alpha,i} (X'WX)^{-1} omits the n factor from
  # A_{alpha,alpha}^{-1} = -n (X'WX)^{-1}, so we multiply by n.
  correction <- n * IF_alpha_fit %*% t(A_psi_alpha)

  # Total IF for psi: (n x p_psi)
  # IF_{psi,i} = A_{psi,psi}^{-1} (omega_{psi,i} - correction_i)
  adjusted_score <- omega_psi - correction
  IF_psi <- adjusted_score %*% t(A_psi_psi_inv)

  vcov_psi <- crossprod(IF_psi) / n^2
  rownames(vcov_psi) <- names(psi_hat)
  colnames(vcov_psi) <- names(psi_hat)
  vcov_psi
}
