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

  A <- fit_data[[fit$treatment]]
  Y <- fit_data[[fit$outcome]]

  # Treatment residual: R_i = A_i - \hat{E}[A | L_i]
  # e_a is the fitted propensity from the treatment model
  e_a <- stats::predict(trt_model$model, type = "response")
  R <- A - e_a

  # Build the blip design matrix M (n x p_psi). Column 1 is always
  # the intercept (corresponds to psi_intercept); subsequent columns
  # are the modifier variables from the blip specification.
  M <- build_blip_design_matrix(fit_data, blip_spec)
  p_psi <- ncol(M)

  if (!is.null(tf_formula)) {
    # Joint estimation: solve (beta, psi) simultaneously.
    # The stacked system is linear:
    #   [Z'Z        Z'(A*M)    ] [beta] = [Z'Y ]
    #   [(R*M)'Z    (R*M)'(A*M)] [psi ] = [(R*M)'Y]
    # where gamma(A,L;psi) = A * (M %*% psi) and Z beta is the
    # treatment-free model prediction.
    Z <- stats::model.matrix(tf_formula, data = fit_data)
    p_beta <- ncol(Z)

    # Blip design weighted by treatment: B_i = A_i * M_i
    AM <- A * M
    # Blip design weighted by treatment residual: W_i = R_i * M_i
    RM <- R * M

    # Block system: [Z'Z, Z'(AM); (RM)'Z, (RM)'(AM)] [beta; psi] = [Z'Y; (RM)'Y]
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
    names(psi_hat) <- blip_spec$param_names
  } else {
    # Standard approach: solve for psi only.
    # psi = (M' diag(R * A) M)^{-1} M' diag(R) Y
    RA <- R * A
    lhs <- crossprod(M, M * RA)
    rhs <- crossprod(M, R * Y)

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
    names(psi_hat) <- blip_spec$param_names
    beta_hat <- NULL
    Z <- NULL
  }

  list(
    psi_hat = psi_hat,
    M = M,
    R = R,
    A = A,
    Y = Y,
    n_obs = n,
    fit_rows = fit_rows,
    beta_hat = beta_hat,
    Z = Z
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


#' Compute SNM contrast: blip parameters or averaged blip effect
#'
#' @description
#' Two paths depending on whether `treatment_values` is supplied:
#'
#' **Path A** (`treatment_values = NULL`): Returns the blip parameter
#' table directly — each \eqn{\psi_j} is its own estimand with SE and CI.
#'
#' **Path B** (`treatment_values = c(a0, a1)`): Computes the
#' population-averaged blip effect
#' \deqn{\Delta = \frac{1}{n}\sum_i \bigl[\gamma(a_1, L_i; \hat\psi)
#'   - \gamma(a_0, L_i; \hat\psi)\bigr]
#'   = (a_1 - a_0)\bigl(\hat\psi_0 + \sum_j \hat\psi_j \bar{m}_j\bigr)}
#' with delta-method variance from the blip vcov.
#'
#' @references
#' Robins JM (1994). Correcting for non-compliance in randomized trials
#' using structural nested mean models. *Comm Stat Theory Methods*,
#' 23(8), 2379–2412.
#'
#' @param fit A `causatr_fit` with `estimator = "snm"`.
#' @param treatment_values Numeric vector of length 2 or `NULL`.
#' @param ci_method `"sandwich"` only for now.
#' @param conf_level Numeric confidence level.
#' @param call The original `contrast()` call.
#' @return A `causatr_result` object.
#' @noRd
compute_snm_contrast <- function(
  fit,
  treatment_values,
  ci_method,
  conf_level,
  call
) {
  if (ci_method == "bootstrap") {
    rlang::abort(
      c(
        "Bootstrap variance is unavailable for `estimator = \"snm\"`.",
        i = "Use `ci_method = \"sandwich\"`."
      ),
      class = "causatr_snm_bootstrap_pending"
    )
  }

  if (fit$type == "longitudinal") {
    snm_result <- compute_snm_blip_longitudinal(fit)
    vcov_psi <- variance_if_snm_longitudinal(fit, snm_result)
  } else {
    snm_result <- compute_snm_blip_point(fit)
    vcov_psi <- variance_if_snm(fit, snm_result)
  }
  psi_hat <- snm_result$psi_hat
  p_psi <- length(psi_hat)
  z <- stats::qnorm((1 + conf_level) / 2)

  n_target <- snm_result$n_obs

  if (is.null(treatment_values)) {
    # Path A: blip parameter table — one row per parameter.
    # For longitudinal fits, psi_hat contains per-stage parameters
    # (e.g. stage0_psi_intercept, stage1_psi_intercept, ...).
    se_psi <- sqrt(pmax(diag(vcov_psi), 0))
    estimates_dt <- data.table::data.table(
      parameter = names(psi_hat),
      estimate = as.numeric(psi_hat),
      se = se_psi,
      ci_lower = psi_hat - z * se_psi,
      ci_upper = psi_hat + z * se_psi
    )
    contrasts_dt <- data.table::data.table(
      comparison = character(0),
      estimate = numeric(0),
      se = numeric(0),
      ci_lower = numeric(0),
      ci_upper = numeric(0)
    )
    vcov_out <- vcov_psi
  } else {
    # Path B: population-averaged blip effect.
    # For longitudinal fits, treatment_values is not meaningful (the
    # per-stage blip table IS the estimand). Reject with guidance.
    if (fit$type == "longitudinal") {
      rlang::abort(
        c(
          paste0(
            "`treatment_values` is not supported for longitudinal SNM fits."
          ),
          i = paste0(
            "Longitudinal SNMs return per-stage blip parameters directly. ",
            "Use `contrast(fit)` without `treatment_values`."
          )
        ),
        class = "causatr_snm_long_no_treatment_values"
      )
    }

    if (length(treatment_values) != 2L) {
      rlang::abort(
        "`treatment_values` must be a numeric vector of length 2 (e.g. c(0, 1)).",
        class = "causatr_snm_bad_treatment_values"
      )
    }
    a0 <- treatment_values[1]
    a1 <- treatment_values[2]
    delta_a <- a1 - a0

    # For linear blip: gamma(a, l; psi) = a * (psi_0 + sum psi_j m_j)
    # Average effect = (a1 - a0) * (psi_0 + sum psi_j * mean(m_j))
    M <- snm_result$M
    m_bar <- colMeans(M)
    avg_effect <- delta_a * sum(psi_hat * m_bar)

    # Delta method: grad_j = delta_a * mean(m_j)
    grad <- delta_a * m_bar
    var_effect <- as.numeric(t(grad) %*% vcov_psi %*% grad)
    se_effect <- sqrt(max(var_effect, 0))

    comparison_label <- paste0("a=", a1, " vs a=", a0)
    estimates_dt <- data.table::data.table(
      parameter = "avg_blip_effect",
      estimate = avg_effect,
      se = se_effect,
      ci_lower = avg_effect - z * se_effect,
      ci_upper = avg_effect + z * se_effect
    )
    contrasts_dt <- data.table::data.table(
      comparison = comparison_label,
      estimate = avg_effect,
      se = se_effect,
      ci_lower = avg_effect - z * se_effect,
      ci_upper = avg_effect + z * se_effect
    )
    vcov_out <- vcov_psi
  }

  new_causatr_result(
    estimates = estimates_dt,
    contrasts = contrasts_dt,
    type = "difference",
    estimand = "ATE",
    ci_method = ci_method,
    reference = NULL,
    interventions = NULL,
    n = n_target,
    estimator = "snm",
    family = fit$family,
    fit_type = fit$type,
    vcov = vcov_out,
    call = call
  )
}
