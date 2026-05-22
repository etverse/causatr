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
  Y <- snm_result$Y
  n <- snm_result$n_obs
  fit_rows <- snm_result$fit_rows
  p_psi <- length(psi_hat)
  trt_model <- fit$details$treatment_model
  beta_hat <- snm_result$beta_hat
  Z <- snm_result$Z
  has_tf <- !is.null(beta_hat)

  # Treatment design matrices from the point estimator
  td <- snm_result$td
  AM <- td$AM
  RM <- td$RM
  is_cat <- td$is_categorical

  # --- Treatment model IF ---
  trt_fit_idx <- which(trt_model$fit_rows)
  n_total <- nrow(fit$data)
  if (is_cat) {
    trt_prep <- prepare_model_if_multinom(trt_model$model, trt_fit_idx, n_total)
  } else {
    trt_prep <- prepare_model_if(trt_model$model, trt_fit_idx, n_total)
  }

  alpha_hat <- trt_model$alpha_hat
  X_prop <- trt_model$X_prop

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
    p_beta <- length(beta_hat)

    # A_{theta,theta} = -(1/n) [Z'Z, Z'AM; RM'Z, RM'AM]
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

    # Joint score: (phi_beta, phi_psi) per observation
    gamma_hat <- as.numeric(AM %*% psi_hat)
    tf_pred <- as.numeric(Z %*% beta_hat)
    resid_joint <- Y - tf_pred - gamma_hat
    omega_beta <- Z * resid_joint
    omega_psi <- RM * resid_joint
    omega_theta <- cbind(omega_beta, omega_psi)

    # Cross-derivative: only phi_psi depends on alpha
    valid_full <- rep(TRUE, n)
    phi_bar_alpha <- build_snm_phi_alpha_closure(
      td,
      trt_model,
      X_prop,
      resid_joint,
      valid_full
    )
    A_psi_alpha <- -numDeriv::jacobian(phi_bar_alpha, x = alpha_hat)

    A_theta_alpha <- rbind(
      matrix(0, nrow = p_beta, ncol = length(alpha_hat)),
      A_psi_alpha
    )

    correction_theta <- n * IF_alpha_fit %*% t(A_theta_alpha)
    adjusted_score_theta <- omega_theta - correction_theta
    IF_theta <- adjusted_score_theta %*% t(A_theta_theta_inv)
    IF_psi <- IF_theta[, p_beta + seq_len(p_psi), drop = FALSE]
  } else {
    # Bread: -(1/n) RM' AM
    A_psi_psi <- -crossprod(RM, AM) / n

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

    # Blip score: RM * (Y - AM psi)
    gamma_hat <- as.numeric(AM %*% psi_hat)
    resid_blip <- Y - gamma_hat
    omega_psi <- RM * resid_blip

    # Cross-derivative closure
    valid_full <- rep(TRUE, n)
    phi_bar_alpha <- build_snm_phi_alpha_closure(
      td,
      trt_model,
      X_prop,
      resid_blip,
      valid_full
    )
    A_psi_alpha <- -numDeriv::jacobian(phi_bar_alpha, x = alpha_hat)

    correction_trt <- n * IF_alpha_fit %*% t(A_psi_alpha)
    adjusted_score <- omega_psi - correction_trt
    IF_psi <- adjusted_score %*% t(A_psi_psi_inv)
  }

  vcov_psi <- crossprod(IF_psi) / n^2
  rownames(vcov_psi) <- names(psi_hat)
  colnames(vcov_psi) <- names(psi_hat)
  vcov_psi
}


#' Build the phi_alpha closure for SNM cross-derivative computation
#'
#' Returns a function `phi(alpha)` that recomputes the averaged blip
#' estimating equation score at a candidate treatment-model parameter
#' `alpha`. For scalar treatments, the residual is
#' \eqn{R_i(\alpha) = A_i - \mu(X_i \alpha)}. For categorical
#' treatments, the residual is a matrix of per-level multinomial
#' residuals recomputed from the softmax of \eqn{X \alpha}.
#'
#' @param td Treatment design object from `snm_treatment_design()`.
#' @param trt_model A `causatr_treatment_model`.
#' @param X_prop Design matrix for the treatment model.
#' @param resid_blip Numeric vector of blip residuals
#'   \eqn{Y_i - \gamma(A_i, L_i; \hat\psi)}.
#' @param valid_full Logical vector indicating valid (non-NA) rows.
#' @return A function `phi(alpha)` -> numeric vector of length p_psi.
#' @noRd
build_snm_phi_alpha_closure <- function(
  td,
  trt_model,
  X_prop,
  resid_blip,
  valid_full
) {
  if (td$is_categorical) {
    M <- td$M
    A_raw <- td$A_raw
    A_char <- as.character(A_raw)
    non_ref <- td$non_ref_levels
    Km1 <- td$Km1
    p_mod <- td$p_mod
    p_psi <- td$p_psi
    p_x <- ncol(X_prop)

    function(alpha) {
      alpha_mat <- matrix(alpha, nrow = Km1, ncol = p_x, byrow = TRUE)
      eta_mat <- X_prop %*% t(alpha_mat)
      exp_eta <- exp(eta_mat)
      denom <- 1 + rowSums(exp_eta)
      prob_nonref <- exp_eta / denom

      RM_alpha <- matrix(0, length(A_raw), p_psi)
      for (j in seq_len(Km1)) {
        cols <- ((j - 1L) * p_mod + 1L):(j * p_mod)
        D_j <- as.numeric(A_char == non_ref[j])
        R_j <- D_j - prob_nonref[, j]
        RM_alpha[, cols] <- R_j * M
      }
      vals <- RM_alpha * resid_blip
      vals[!valid_full, ] <- 0
      colMeans(vals)
    }
  } else {
    M <- td$M
    A_raw <- td$A_raw
    trt_family <- trt_model$model$family

    function(alpha) {
      eta <- as.numeric(X_prop %*% alpha)
      mu <- trt_family$linkinv(eta)
      R_alpha <- A_raw - mu
      vals <- M * (R_alpha * resid_blip)
      vals[!valid_full, ] <- 0
      colMeans(vals)
    }
  }
}
