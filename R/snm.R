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

  # Capture user dots early — needed by both point and longitudinal
  # branches for treatment model fitting and bootstrap replay.
  dots <- list(...)

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
    # --- Longitudinal: per-period treatment models ---
    time_points <- sort(unique(data[[time]]))
    n_times <- length(time_points)

    if (n_times < 2L) {
      rlang::abort(
        c(
          paste0(
            "`type = 'longitudinal'` requires at least 2 distinct time ",
            "points (got ",
            n_times,
            ")."
          ),
          i = "Use `type = 'point'` for single-period data."
        ),
        class = "causatr_longitudinal_too_few_times"
      )
    }

    max_lag <- if (is.infinite(history)) {
      n_times - 1L
    } else {
      as.integer(history)
    }

    # Parse baseline confounder terms (strip A:modifier interactions so
    # the propensity model is A ~ L, matching the point-treatment path).
    em_conf_trt <- confounders
    em_info_trt <- parse_effect_mod(em_conf_trt, treatment)
    baseline_terms <- em_info_trt$confounder_terms
    if (length(baseline_terms) == 0L) {
      baseline_terms <- character(0L)
    }

    tv_vars <- if (!is.null(confounders_tv)) {
      all.vars(confounders_tv)
    } else {
      character(0)
    }

    treatment_models_by_time <- vector("list", n_times)
    fit_data_by_time <- vector("list", n_times)
    per_period_formula <- vector("list", n_times)

    for (k in seq_len(n_times)) {
      rows_k <- data[[time]] == time_points[k]
      data_k <- data[rows_k]
      available_lags <- min(k - 1L, max_lag)

      ps_formula <- build_longitudinal_ps_formula(
        treatment = treatment,
        baseline_terms = baseline_terms,
        tv_vars = tv_vars,
        available_lags = available_lags,
        data_at_time = data_k
      )

      tm_args <- list(
        data = data_k,
        treatment = treatment,
        confounders = remove_response(ps_formula),
        model_fn = prop_fn,
        propensity_family = propensity_family
      )
      if (!is.null(weights)) {
        tm_args$weights <- weights[rows_k]
      }
      tm_k <- do.call(fit_treatment_model, c(tm_args, dots))

      treatment_models_by_time[[k]] <- tm_k
      fit_data_by_time[[k]] <- data_k
      per_period_formula[[k]] <- ps_formula
    }
    names(treatment_models_by_time) <- as.character(time_points)
    names(fit_data_by_time) <- as.character(time_points)
    names(per_period_formula) <- as.character(time_points)

    trt_model <- list(
      models_by_time = treatment_models_by_time,
      time_points = time_points,
      n_times = n_times
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

  details <- list(
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
  )

  # Longitudinal-specific details for the variance engine and contrast
  if (type == "longitudinal") {
    details$treatment_models_by_time <- trt_model$models_by_time
    details$fit_data_by_time <- fit_data_by_time
    details$per_period_formula <- per_period_formula
    details$time_points <- trt_model$time_points
    details$n_times <- trt_model$n_times
    details$max_lag <- max_lag
    details$baseline_terms <- baseline_terms
    details$tv_vars <- tv_vars
  }

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
    details = details,
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
