#' Pool per-imputation causal estimates via Rubin's rules
#'
#' @description
#' Combines the point estimates and variances from a set of per-imputation
#' [contrast()] results into a single pooled [causatr_result] using Rubin's
#' (1987) rules. Each row of the intervention-means table and each row of the
#' contrasts table is pooled independently: difference contrasts and means are
#' pooled on the natural scale, while ratio / odds-ratio contrasts are pooled
#' on the log scale (the scale where the estimator is approximately normal)
#' and back-transformed.
#'
#' @param collected A list produced by `mi_collect_imputations()` with the
#'   stacked per-imputation estimate and standard-error matrices plus the
#'   metadata needed to rebuild a `causatr_result`. See that function for the
#'   exact element layout.
#' @param conf_level Numeric scalar confidence level for the pooled intervals.
#'   Default `0.95`.
#'
#' @returns A [causatr_result] with `ci_method = "rubin"`. The `estimates` and
#'   `contrasts` data.tables hold the pooled point estimate, total-variance
#'   standard error, and Barnard-Rubin \eqn{t}-based confidence bounds. The
#'   per-row pooling diagnostics (between/within variance, degrees of freedom,
#'   fraction of missing information, raw per-imputation draws) are attached as
#'   the `"mi_details"` attribute.
#'
#' @details
#' For estimate \eqn{\hat{Q}_i} and variance \eqn{U_i} from imputation
#' \eqn{i = 1, \ldots, m}:
#' \deqn{\bar{Q} = m^{-1} \sum_i \hat{Q}_i, \quad
#'       \bar{U} = m^{-1} \sum_i U_i, \quad
#'       B = (m-1)^{-1} \sum_i (\hat{Q}_i - \bar{Q})^2}
#' The total variance is \eqn{T = \bar{U} + (1 + 1/m) B} and the
#' Barnard-Rubin (1999) degrees of freedom combine the small-sample complete
#' data df with the fraction of missing information. Confidence intervals are
#' \eqn{\bar{Q} \pm t_{\nu, 1-\alpha/2} \sqrt{T}}.
#'
#' This mirrors `mice::pool.scalar()` exactly when fed the same complete-data
#' degrees of freedom, which is the cross-check used in the test suite.
#'
#' @references
#' Rubin DB (1987). *Multiple Imputation for Nonresponse in Surveys*. Wiley.
#'
#' Barnard J, Rubin DB (1999). Small-sample degrees of freedom with multiple
#' imputation. *Biometrika* 86:948-955.
#'
#' @seealso [causat_mice()], [pool_boot_mi()]
#' @noRd
pool_rubin <- function(collected, conf_level = 0.95) {
  m <- collected$m
  meta <- collected$meta

  # Pool the intervention-means (or SNM blip-parameter) rows, one column
  # per row, on the natural scale. `pool_table_rubin()` returns the pooled
  # estimate / SE / df / CI for every column together with the per-row
  # diagnostics that downstream methods read from "mi_details".
  est_pool <- pool_table_rubin(
    est_mat = collected$est_mat,
    se_mat = collected$est_se_mat,
    scale = "linear",
    m = m,
    n = collected$n,
    conf_level = conf_level
  )

  # Contrasts pool on the scale implied by `type`: difference is linear,
  # ratio / OR are log-scale (delta-method variance var(log X) = var(X)/X^2)
  # then exponentiated back so the reported estimate stays on the ratio
  # scale -- matching how contrast() itself reports these.
  con_scale <- if (identical(meta$type, "difference")) "linear" else "log"
  con_pool <- if (length(collected$con_labels) > 0L) {
    pool_table_rubin(
      est_mat = collected$con_mat,
      se_mat = collected$con_se_mat,
      scale = con_scale,
      m = m,
      n = collected$n,
      conf_level = conf_level
    )
  } else {
    NULL
  }

  estimates <- build_pooled_table(
    labels = collected$est_labels,
    label_col = meta$term_col,
    by = collected$est_by,
    pool = est_pool
  )
  contrasts <- if (!is.null(con_pool)) {
    build_pooled_table(
      labels = collected$con_labels,
      label_col = "comparison",
      by = collected$con_by,
      pool = con_pool
    )
  } else {
    # SNM Path A (no contrast row) still needs a well-formed empty table
    # with the canonical columns so downstream S3 methods do not trip.
    data.table::data.table(
      comparison = character(0),
      estimate = numeric(0),
      se = numeric(0),
      ci_lower = numeric(0),
      ci_upper = numeric(0)
    )
  }

  # Pooled vcov is the diagonal of the per-row total variances. Rubin's
  # rules here are applied per row independently (the design contract), so
  # the off-diagonal covariance across interventions is not pooled; a
  # diagonal matrix is the honest representation of what was computed.
  vcov <- diag(est_pool$t, nrow = length(est_pool$t))

  result <- new_causatr_result(
    estimates = estimates,
    contrasts = contrasts,
    type = meta$type,
    estimand = meta$estimand,
    ci_method = "rubin",
    reference = meta$reference,
    interventions = meta$interventions,
    n = collected$n,
    estimator = meta$estimator,
    family = meta$family,
    fit_type = meta$fit_type,
    vcov = vcov,
    boot_t = NULL,
    boot_info = NULL,
    call = meta$call
  )

  attr(result, "mi_details") <- list(
    pool_method = "rubin",
    m = m,
    estimates = est_pool$diagnostics,
    contrasts = if (!is.null(con_pool)) con_pool$diagnostics else NULL
  )
  result
}

#' Pool a stacked per-imputation estimate / SE matrix column-by-column
#'
#' @param est_mat Numeric `m x p` matrix of per-imputation point estimates
#'   (rows = imputations, columns = parameters / contrasts).
#' @param se_mat Numeric `m x p` matrix of per-imputation standard errors,
#'   aligned to `est_mat`.
#' @param scale Either `"linear"` (pool the estimates directly) or `"log"`
#'   (pool `log(estimate)` with delta-method variance, then exponentiate).
#' @param m Integer number of imputations.
#' @param n Integer complete-data sample size, used for the Barnard-Rubin
#'   complete-data degrees of freedom `dfcom = n - p`.
#' @param conf_level Numeric confidence level.
#' @returns A list with vectors `estimate`, `se`, `df`, `ci_lower`,
#'   `ci_upper`, `t` (total variance), and a `diagnostics` data.table.
#' @noRd
pool_table_rubin <- function(est_mat, se_mat, scale, m, n, conf_level) {
  p <- ncol(est_mat)
  # Complete-data residual df convention for a marginal causal estimand:
  # one degree of freedom consumed per estimated parameter. Passed
  # identically to the mice::pool.scalar() oracle in the tests.
  dfcom <- max(n - p, 1L)
  alpha <- (1 - conf_level) / 2

  out <- lapply(seq_len(p), function(j) {
    q_lin <- est_mat[, j]
    u_lin <- se_mat[, j]^2
    if (identical(scale, "log")) {
      # Delta method: var(log X) = var(X) / X^2. Pool on the log scale
      # where the ratio/OR estimator is approximately normal.
      q <- log(q_lin)
      u <- u_lin / q_lin^2
    } else {
      q <- q_lin
      u <- u_lin
    }
    rs <- rubin_scalar(q, u, m = m, dfcom = dfcom)
    crit <- stats::qt(1 - alpha, df = rs$df)
    lo_link <- rs$qbar - crit * sqrt(rs$t)
    hi_link <- rs$qbar + crit * sqrt(rs$t)
    if (identical(scale, "log")) {
      est <- exp(rs$qbar)
      # SE on the reported (linear) scale via the inverse delta map.
      se <- est * sqrt(rs$t)
      lo <- exp(lo_link)
      hi <- exp(hi_link)
    } else {
      est <- rs$qbar
      se <- sqrt(rs$t)
      lo <- lo_link
      hi <- hi_link
    }
    list(
      estimate = est,
      se = se,
      df = rs$df,
      ci_lower = lo,
      ci_upper = hi,
      t = rs$t,
      b = rs$b,
      ubar = rs$ubar,
      fmi = rs$fmi
    )
  })

  pull <- function(field) vapply(out, `[[`, numeric(1), field)
  diagnostics <- data.table::data.table(
    between = pull("b"),
    within = pull("ubar"),
    total = pull("t"),
    df = pull("df"),
    fmi = pull("fmi")
  )
  list(
    estimate = pull("estimate"),
    se = pull("se"),
    df = pull("df"),
    ci_lower = pull("ci_lower"),
    ci_upper = pull("ci_upper"),
    t = pull("t"),
    diagnostics = diagnostics
  )
}

#' Scalar Rubin's-rules combination
#'
#' @param Q Numeric vector of length `m` of per-imputation point estimates
#'   (on the pooling scale).
#' @param U Numeric vector of length `m` of per-imputation variances (on the
#'   pooling scale).
#' @param m Integer number of imputations.
#' @param dfcom Numeric complete-data degrees of freedom for the Barnard-Rubin
#'   correction.
#' @returns A list with `qbar` (pooled estimate), `ubar` (within variance),
#'   `b` (between variance), `t` (total variance), `df` (Barnard-Rubin df),
#'   `r` (relative increase in variance), `fmi` (fraction of missing
#'   information).
#' @noRd
rubin_scalar <- function(Q, U, m, dfcom) {
  qbar <- mean(Q)
  ubar <- mean(U)
  if (m == 1L) {
    # Single imputation is degenerate: there is no between-imputation
    # variability to estimate, so the total variance collapses to the
    # within-imputation variance and df is the complete-data df. This is
    # equivalent to a single complete-data analysis and is flagged with a
    # warning upstream in causat_mice().
    return(list(
      qbar = qbar,
      ubar = ubar,
      b = 0,
      t = ubar,
      df = dfcom,
      r = 0,
      fmi = NA_real_
    ))
  }
  b <- stats::var(Q)
  t <- ubar + (1 + 1 / m) * b
  r <- (1 + 1 / m) * b / ubar
  df <- barnard_rubin_df(m = m, b = b, t = t, dfcom = dfcom)
  fmi <- (r + 2 / (df + 3)) / (r + 1)
  list(qbar = qbar, ubar = ubar, b = b, t = t, df = df, r = r, fmi = fmi)
}

#' Barnard-Rubin small-sample degrees of freedom
#'
#' @param m Integer number of imputations.
#' @param b Numeric between-imputation variance.
#' @param t Numeric total variance.
#' @param dfcom Numeric complete-data degrees of freedom (`Inf` for the
#'   large-sample Rubin (1987) df).
#' @returns Numeric degrees of freedom. Replicates `mice:::barnard.rubin()`.
#' @noRd
barnard_rubin_df <- function(m, b, t, dfcom = Inf) {
  # lambda is the proportion of total variance attributable to missingness.
  # Guard the t == 0 corner (no variance at all) so lambda stays finite.
  lambda <- if (t > 0) (1 + 1 / m) * b / t else 0
  dfold <- (m - 1) / (lambda^2)
  if (is.infinite(dfcom)) {
    return(dfold)
  }
  tmp <- (1 - lambda) * (1 + dfcom) * dfcom
  (m - 1) * tmp / ((dfcom + 3) * (m - 1) + (lambda^2) * tmp)
}

#' Build a pooled estimates / contrasts data.table
#'
#' @param labels Character vector of row labels (intervention / parameter
#'   names or contrast comparison strings).
#' @param label_col Name of the label column: `"intervention"`,
#'   `"parameter"`, or `"comparison"`.
#' @param by Optional character/factor vector of `by`-subgroup labels aligned
#'   to `labels`, or `NULL`.
#' @param pool A list from `pool_table_rubin()` (or the Boot MI equivalent)
#'   with `estimate`, `se`, `ci_lower`, `ci_upper`.
#' @returns A data.table with the canonical `causatr_result` columns plus an
#'   optional `by` column.
#' @noRd
build_pooled_table <- function(labels, label_col, by, pool) {
  dt <- data.table::data.table(
    label = labels,
    estimate = pool$estimate,
    se = pool$se,
    ci_lower = pool$ci_lower,
    ci_upper = pool$ci_upper
  )
  data.table::setnames(dt, "label", label_col)
  if (!is.null(by)) {
    dt[, by := by]
    data.table::setcolorder(dt, c(label_col, "by"))
  }
  dt[]
}
