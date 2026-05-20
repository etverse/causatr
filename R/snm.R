#' Fit a Structural Nested Mean Model (SNMM)
#'
#' Internal workhorse called by [causat()] when `estimator = "snm"`.
#' Fits the treatment model \eqn{\hat{E}[A \mid L]} via
#' `fit_treatment_model()` and stores the blip parameterisation derived
#' from effect-modification terms in the confounders formula. No point
#' estimation happens here — the g-estimating equation is solved by
#' [contrast()] (analogous to how ICE defers model fitting).
#'
#' @param data data.table of the (prepared) dataset.
#' @param outcome Character outcome column name.
#' @param treatment Character treatment column name.
#' @param confounders One-sided formula of confounders (resolved:
#'   `confounders_treatment` or the unified `confounders`).
#' @param confounders_tv One-sided formula of time-varying confounders or
#'   `NULL`.
#' @param family Character or family object for the outcome distribution.
#' @param estimand Character estimand (`"ATE"` only for now; ATT/ATC
#'   gated downstream).
#' @param type `"point"` or `"longitudinal"`.
#' @param history Integer Markov order for longitudinal lag expansion.
#' @param weights Numeric vector of external weights or `NULL`.
#' @param propensity_model_fn Model fitting function for the treatment
#'   model (default `stats::glm`).
#' @param propensity_family Family override for the treatment model or
#'   `NULL`.
#' @param id Character ID column name or `NULL`.
#' @param time Character time column name or `NULL`.
#' @param call The original `causat()` call.
#' @param target Character column name of the sampling indicator or `NULL`.
#' @param confounders_outcome One-sided formula for outcome confounders or
#'   `NULL`.
#' @param confounders_outcome_raw Raw user-supplied `confounders_outcome`.
#' @param confounders_treatment_raw Raw user-supplied `confounders_treatment`.
#' @param confounders_censoring_raw Raw user-supplied `confounders_censoring`.
#' @param confounders_sampling_raw Raw user-supplied `confounders_sampling`.
#' @param confounders_tv_outcome_raw Raw user-supplied `confounders_tv_outcome`.
#' @param confounders_tv_treatment_raw Raw user-supplied
#'   `confounders_tv_treatment`.
#' @param treatment_free One-sided formula or `NULL`. Specifies the
#'   treatment-free outcome model \eqn{E[Y \mid L]} for efficiency
#'   augmentation. When non-`NULL`, the model's predictions are
#'   subtracted from Y before g-estimation, projecting out the
#'   \eqn{L \to Y} association and reducing variance.
#' @param model_fn Model fitting function for the treatment-free model
#'   (default `stats::glm`).
#' @param ... Extra arguments forwarded to the treatment model fitter.
#' @return A `causatr_fit` with `estimator = "snm"`.
#' @noRd
fit_snm <- function(
  data,
  outcome,
  treatment,
  confounders,
  confounders_tv,
  family,
  estimand,
  type,
  history,
  weights,
  propensity_model_fn,
  propensity_family = NULL,
  id = NULL,
  time = NULL,
  call,
  target = NULL,
  confounders_outcome = NULL,
  confounders_outcome_raw = NULL,
  confounders_treatment_raw = NULL,
  confounders_censoring_raw = NULL,
  confounders_sampling_raw = NULL,
  confounders_tv_outcome_raw = NULL,
  confounders_tv_treatment_raw = NULL,
  treatment_free = NULL,
  model_fn = stats::glm,
  ...
) {
  # --- Rejection: multivariate treatment ---
  if (length(treatment) > 1L) {
    rlang::abort(
      c(
        "Multivariate treatment is not supported for `estimator = \"snm\"`.",
        i = paste0(
          "Multivariate SNMMs are out of scope. Use `estimator = \"ipw\"` ",
          "or `\"gcomp\"` for joint interventions on multiple treatments."
        )
      ),
      class = "causatr_snm_multivariate"
    )
  }

  # Parse effect-modification terms from confounders. The EM terms
  # (`A:modifier`) define the blip parameterisation: each interaction
  # contributes a modifier column whose coefficient in the blip
  # function is estimated by the g-estimating equation. Pure confounder
  # terms enter the treatment model only.
  em_conf <- confounders_outcome %||% confounders
  em_info <- parse_effect_mod(em_conf, treatment)

  # Inform when no EM terms are present — the SNMM reduces to a single
  # ATE parameter, which is fine but worth flagging because the main
  # motivation for SNMs is effect-modification identification.
  if (!em_info$has_em) {
    rlang::inform(
      c(
        "No effect-modification terms (e.g. `A:modifier`) found in confounders.",
        i = paste0(
          "The SNMM blip reduces to a single constant-effect parameter ",
          "(equivalent to an ATE). Add `A:modifier` terms to the confounders ",
          "formula to estimate effect modification."
        )
      ),
      class = "causatr_snm_no_em"
    )
  }

  # Build blip specification from parsed EM terms. Each EM term
  # contributes modifier columns; the intercept (`A * psi_0`) is
  # always included.
  blip_spec <- build_blip_spec(em_info, treatment)

  # --- Fit rows: exclude missing outcomes ---
  fit_rows <- get_fit_rows(data, outcome, target = target)
  fit_data <- data[fit_rows]

  # Default to GLM for the propensity model if no fitter is supplied
  prop_fn <- propensity_model_fn %||% stats::glm

  if (type == "point") {
    # Fit a single treatment model E[A | L]. `fit_treatment_model()`
    # calls `build_ps_formula()` internally, which strips any EM terms
    # from `confounders` so the propensity model is A ~ L (no A:modifier).
    trt_model <- fit_treatment_model(
      data = fit_data,
      treatment = treatment,
      confounders = confounders,
      model_fn = prop_fn,
      propensity_family = propensity_family,
      weights = if (is.null(weights)) NULL else weights[fit_rows],
      ...
    )
  } else {
    rlang::abort(
      c(
        "Longitudinal SNMMs require per-period treatment models.",
        i = "Use `estimator = \"ice\"` or `estimator = \"ipw\"` for longitudinal data."
      ),
      class = "causatr_snm_longitudinal_pending"
    )
  }

  # --- Treatment-free model specification ---
  # When treatment_free is specified, store the formula. The actual model
  # is fit jointly with the blip parameters in compute_snm_blip_point(),
  # following the joint EE approach of Vansteelandt & Joffe (2014) and
  # DTRreg: (beta, psi) are solved simultaneously from
  #   E[Z_i (Y_i - Z_i beta - gamma(A_i, L_i; psi))] = 0
  #   E[R_i m_i (Y_i - Z_i beta - gamma(A_i, L_i; psi))] = 0
  if (!is.null(treatment_free)) {
    if (!inherits(treatment_free, "formula")) {
      rlang::abort(
        "`treatment_free` must be a one-sided formula (e.g. `~ L`).",
        class = "causatr_snm_bad_treatment_free"
      )
    }
  }

  # Capture user dots for bootstrap replay
  dots <- list(...)

  new_causatr_fit(
    model = NULL,
    data = data,
    treatment = treatment,
    outcome = outcome,
    confounders = confounders,
    confounders_tv = confounders_tv,
    confounders_outcome = confounders_outcome_raw,
    confounders_treatment = confounders_treatment_raw,
    confounders_censoring = confounders_censoring_raw,
    confounders_sampling = confounders_sampling_raw,
    confounders_tv_outcome = confounders_tv_outcome_raw,
    confounders_tv_treatment = confounders_tv_treatment_raw,
    family = family,
    estimator = "snm",
    type = type,
    estimand = estimand,
    id = id,
    time = time,
    censoring = NULL,
    history = history,
    numerator = NULL,
    weights_obj = NULL,
    match_obj = NULL,
    call = call,
    details = list(
      treatment_model = trt_model,
      blip_spec = blip_spec,
      blip_type = "linear",
      em_info = em_info,
      fit_rows = fit_rows,
      weights = weights,
      dots = dots,
      propensity_model_fn = prop_fn,
      propensity_family = propensity_family,
      treatment_free = treatment_free,
      model_fn = model_fn
    ),
    target = target
  )
}


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

  snm_result <- compute_snm_blip_point(fit)
  vcov_psi <- variance_if_snm(fit, snm_result)
  psi_hat <- snm_result$psi_hat
  p_psi <- length(psi_hat)
  z <- stats::qnorm((1 + conf_level) / 2)

  if (is.null(treatment_values)) {
    # Path A: blip parameter table
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
    n_target <- snm_result$n_obs
    vcov_out <- vcov_psi
  } else {
    # Path B: population-averaged blip effect
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
    n_target <- snm_result$n_obs
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


#' Build blip specification from parsed effect-modification terms
#'
#' Translates the `causatr_em_info` object into a blip specification
#' that describes the linear blip function
#' \eqn{\gamma(a, l; \psi) = a \cdot (\psi_0 + \sum_j \psi_j \cdot m_j)}.
#' The returned list enumerates modifier columns and names the blip
#' parameters.
#'
#' @param em_info A `causatr_em_info` from `parse_effect_mod()`.
#' @param treatment Character scalar treatment column name.
#' @return A list with:
#'   \describe{
#'     \item{modifier_cols}{Character vector of modifier column names
#'       (empty if no EM terms).}
#'     \item{param_names}{Character vector of blip parameter names,
#'       always starting with `"psi_intercept"`.}
#'     \item{n_params}{Integer number of blip parameters.}
#'   }
#' @noRd
build_blip_spec <- function(em_info, treatment) {
  # Extract modifier variables from the EM info (parsed from A:M terms
  # in the confounders formula by parse_effect_mod())
  modifier_cols <- em_info$modifier_vars
  param_names <- "psi_intercept"

  if (length(modifier_cols) > 0L) {
    param_names <- c(param_names, paste0("psi_", modifier_cols))
  }

  list(
    modifier_cols = modifier_cols,
    param_names = param_names,
    n_params = length(param_names)
  )
}
