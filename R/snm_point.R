#' Solve the point-treatment g-estimating equation for blip parameters
#'
#' For the linear blip \eqn{\gamma(a, l; \psi) = a \cdot (\psi_0 +
#' \sum_j \psi_j m_j)}, the closed-form solution is
#' \deqn{\hat\psi = \Bigl(\sum_i \mathbf{m}_i^{\otimes 2}
#'   R_i A_i\Bigr)^{-1}
#'   \sum_i \mathbf{m}_i R_i Y_i,}
#' where \eqn{R_i = A_i - \hat{E}[A \mid L_i]} is the treatment
#' residual from the propensity model and \eqn{\mathbf{m}_i =
#' (1, m_{1,i}, \ldots)^\top} is the blip design vector.
#'
#' When a treatment-free model is present, \eqn{(\hat\beta, \hat\psi)}
#' are solved jointly from the stacked system (Vansteelandt & Joffe
#' 2014, DTRreg's `tf.mod`):
#' \deqn{\sum_i Z_i (Y_i - Z_i \beta - \gamma(A_i, L_i; \psi)) = 0}
#' \deqn{\sum_i R_i m_i (Y_i - Z_i \beta - \gamma(A_i, L_i; \psi)) = 0}
#' where \eqn{Z_i} is the treatment-free design matrix. The joint system
#' is linear in \eqn{(\beta, \psi)} and solved by a single matrix
#' inversion. This provides better efficiency than the standard approach
#' because the treatment-free model absorbs the \eqn{L \to Y} variance.
#'
#' @param fit A `causatr_fit` with `estimator = "snm"` and
#'   `type = "point"`.
#' @return A list with:
#'   \describe{
#'     \item{psi_hat}{Named numeric vector of blip parameter estimates.}
#'     \item{M}{Matrix (n_fit x p_psi) of blip design vectors.}
#'     \item{R}{Numeric vector of treatment residuals.}
#'     \item{A}{Numeric vector of observed treatment values.}
#'     \item{Y}{Numeric vector of observed outcomes.}
#'     \item{n_obs}{Integer number of observations used.}
#'     \item{fit_rows}{Logical vector of rows used.}
#'     \item{beta_hat}{Named numeric vector of treatment-free model
#'       parameters or `NULL`.}
#'     \item{Z}{Treatment-free design matrix (n x p_beta) or `NULL`.}
#'   }
#' @noRd
compute_snm_blip_point <- function(fit) {
  blip_spec <- fit$details$blip_spec
  trt_model <- fit$details$treatment_model
  tf_formula <- fit$details$treatment_free
  fit_rows <- fit$details$fit_rows
  data <- fit$data
  fit_data <- data[fit_rows]
  n <- nrow(fit_data)

  Y <- fit_data[[fit$outcome]]

  # Build treatment design matrices AM and RM. For scalar treatments
  # (binary/continuous/count), AM = A*M and RM = R*M. For categorical
  # (K levels), expands to (K-1) blocks with per-level indicators and
  # multinomial residuals.
  td <- snm_treatment_design(fit_data, fit$treatment, blip_spec, trt_model)
  AM <- td$AM
  RM <- td$RM
  p_psi <- td$p_psi

  if (!is.null(tf_formula)) {
    # Joint estimation: solve (beta, psi) simultaneously.
    # [Z'Z,    Z'AM  ] [beta] = [Z'Y  ]
    # [RM'Z,   RM'AM ] [psi ] = [RM'Y ]
    Z <- stats::model.matrix(tf_formula, data = fit_data)
    p_beta <- ncol(Z)

    lhs <- rbind(
      cbind(crossprod(Z, Z), crossprod(Z, AM)),
      cbind(crossprod(RM, Z), crossprod(RM, AM))
    )
    rhs <- c(
      as.numeric(crossprod(Z, Y)),
      as.numeric(crossprod(RM, Y))
    )

    theta_hat <- tryCatch(
      as.numeric(solve(lhs, rhs)),
      error = function(e) {
        rlang::abort(
          c(
            "Joint treatment-free + blip system is singular.",
            i = paste0(
              "This can happen when the treatment-free model is ",
              "collinear with the blip design. Check your ",
              "`treatment_free` formula."
            )
          ),
          class = "causatr_snm_singular",
          parent = e
        )
      }
    )

    beta_hat <- theta_hat[seq_len(p_beta)]
    names(beta_hat) <- colnames(Z)
    psi_hat <- theta_hat[p_beta + seq_len(p_psi)]
    names(psi_hat) <- td$param_names
  } else {
    # psi = (RM' AM)^{-1} RM' Y
    lhs <- crossprod(RM, AM)
    rhs <- crossprod(RM, Y)

    psi_hat <- tryCatch(
      as.numeric(solve(lhs, rhs)),
      error = function(e) {
        rlang::abort(
          c(
            "G-estimating equation is singular (blip design matrix is rank-deficient).",
            i = paste0(
              "This can happen when a modifier is constant or collinear ",
              "with the intercept. Check the effect-modification terms."
            )
          ),
          class = "causatr_snm_singular",
          parent = e
        )
      }
    )
    names(psi_hat) <- td$param_names
    beta_hat <- NULL
    Z <- NULL
  }

  list(
    psi_hat = psi_hat,
    M = td$M,
    R = td$R_raw,
    A = td$A_raw,
    Y = Y,
    n_obs = n,
    fit_rows = fit_rows,
    beta_hat = beta_hat,
    Z = Z,
    td = td
  )
}


#' Build the blip design matrix from data and blip specification
#'
#' Constructs the n x p_psi matrix \eqn{\mathbf{M}} where each row is
#' \eqn{(1, m_{1,i}, m_{2,i}, \ldots)} — the intercept followed by the
#' modifier variable values for observation i.
#'
#' @details
#' The linear blip function is
#' \eqn{\gamma(a, l; \psi) = a \cdot (\psi_0 + \sum_j \psi_j m_j(l))}.
#' \eqn{\mathbf{M}} is the design matrix such that the blip for
#' individual i is \eqn{A_i (\mathbf{M}_i \psi)}.
#'
#' @param data data.table or data.frame of observations.
#' @param blip_spec A blip specification from `build_blip_spec()`.
#' @return A numeric matrix with `nrow(data)` rows and
#'   `blip_spec$n_params` columns.
#' @noRd
build_blip_design_matrix <- function(data, blip_spec) {
  n <- nrow(data)
  p <- blip_spec$n_params
  # Column 1 is the intercept (always 1), so psi_intercept captures
  # the main blip effect. Subsequent columns are modifier values.
  M <- matrix(1, nrow = n, ncol = p)
  colnames(M) <- blip_spec$param_names

  mod_cols <- blip_spec$modifier_cols
  if (length(mod_cols) > 0L) {
    for (j in seq_along(mod_cols)) {
      col_name <- mod_cols[j]
      if (!col_name %in% names(data)) {
        rlang::abort(
          paste0(
            "Modifier column '",
            col_name,
            "' not found in data. ",
            "Check the effect-modification terms in confounders."
          ),
          class = "causatr_snm_missing_modifier"
        )
      }
      M[, j + 1L] <- as.numeric(data[[col_name]])
    }
  }

  M
}


#' Build the SNM treatment design matrices for g-estimation
#'
#' Computes the treatment-action matrix \eqn{AM} and treatment-residual
#' matrix \eqn{RM} that enter the g-estimating equation. For scalar
#' treatments (binary, continuous, count), these are element-wise
#' products \eqn{A \cdot M} and \eqn{R \cdot M} where
#' \eqn{R_i = A_i - \hat{E}[A \mid L_i]}. For categorical treatments
#' with \eqn{K} levels, the matrices expand to \eqn{(K-1)} blocks —
#' one per non-reference level — with level-specific indicators
#' \eqn{D_{i,j} = 1\{A_i = j\}} and multinomial residuals
#' \eqn{R_{i,j} = D_{i,j} - P(A = j \mid L_i)}.
#'
#' The g-estimating equation is then uniformly:
#' \deqn{\sum_i RM_i \cdot (Y_i - AM_i \cdot \psi) = 0,}
#' solved by \eqn{\hat\psi = (RM' AM)^{-1} RM' Y}.
#'
#' @param data data.table of fitted observations.
#' @param treatment Character treatment column name.
#' @param blip_spec A blip specification from `build_blip_spec()`.
#' @param trt_model A `causatr_treatment_model`.
#' @return A list with:
#'   \describe{
#'     \item{AM}{Matrix (n x p_psi): treatment-action times blip design.}
#'     \item{RM}{Matrix (n x p_psi): treatment-residual times blip design.}
#'     \item{M}{Matrix (n x p_mod): raw blip design (modifiers only).}
#'     \item{A_raw}{Treatment values (numeric or factor).}
#'     \item{R_raw}{Numeric residual vector for scalar treatments;
#'       `NULL` for categorical.}
#'     \item{p_psi}{Integer total number of psi parameters.}
#'     \item{param_names}{Character vector of parameter names.}
#'     \item{is_categorical}{Logical.}
#'     \item{trt_levels}{Character vector of treatment levels
#'       (categorical only; `NULL` otherwise).}
#'     \item{non_ref_levels}{Non-reference levels (categorical only).}
#'     \item{Km1}{Integer number of non-reference levels (1 for scalar).}
#'     \item{p_mod}{Integer number of per-level blip parameters.}
#'   }
#' @noRd
snm_treatment_design <- function(data, treatment, blip_spec, trt_model) {
  n <- nrow(data)
  M <- build_blip_design_matrix(data, blip_spec)
  p_mod <- ncol(M)
  A_raw <- data[[treatment]]

  if (trt_model$family == "categorical") {
    trt_levels <- trt_model$levels
    K <- length(trt_levels)
    Km1 <- K - 1L
    non_ref <- trt_levels[-1]

    prob_raw <- stats::predict(trt_model$model, newdata = data, type = "probs")
    if (is.null(dim(prob_raw))) {
      prob_mat <- cbind(1 - prob_raw, prob_raw)
      colnames(prob_mat) <- trt_levels
    } else {
      prob_mat <- prob_raw
    }

    A_char <- as.character(A_raw)
    p_total <- Km1 * p_mod
    AM <- matrix(0, nrow = n, ncol = p_total)
    RM <- matrix(0, nrow = n, ncol = p_total)
    param_names <- character(p_total)

    for (j in seq_len(Km1)) {
      cols <- ((j - 1L) * p_mod + 1L):(j * p_mod)
      lev <- non_ref[j]
      D_j <- as.numeric(A_char == lev)
      R_j <- D_j - prob_mat[, lev]
      AM[, cols] <- D_j * M
      RM[, cols] <- R_j * M
      param_names[cols] <- paste0("level_", lev, "_", blip_spec$param_names)
    }
    colnames(AM) <- param_names
    colnames(RM) <- param_names

    list(
      AM = AM,
      RM = RM,
      M = M,
      A_raw = A_raw,
      R_raw = NULL,
      p_psi = p_total,
      param_names = param_names,
      is_categorical = TRUE,
      trt_levels = trt_levels,
      non_ref_levels = non_ref,
      Km1 = Km1,
      p_mod = p_mod
    )
  } else {
    e_a <- stats::predict(trt_model$model, newdata = data, type = "response")
    R <- A_raw - e_a
    AM <- A_raw * M
    RM <- R * M
    colnames(AM) <- blip_spec$param_names
    colnames(RM) <- blip_spec$param_names

    list(
      AM = AM,
      RM = RM,
      M = M,
      A_raw = A_raw,
      R_raw = R,
      p_psi = p_mod,
      param_names = blip_spec$param_names,
      is_categorical = FALSE,
      trt_levels = NULL,
      non_ref_levels = NULL,
      Km1 = 1L,
      p_mod = p_mod
    )
  }
}
