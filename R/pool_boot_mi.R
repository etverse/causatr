#' Pool causal estimates via bootstrap-then-impute (Boot MI)
#'
#' @description
#' Implements von Hippel's (2020) two-stage bootstrap-then-impute variance
#' estimator. The *incomplete* data underlying a `mids` object is resampled
#' with replacement \eqn{B} times; each bootstrap resample is re-imputed
#' \eqn{M} times with `mice::mice()`, and [causat()] + [contrast()] are run on
#' every completed dataset. A one-way random-effects decomposition separates
#' the between-bootstrap sampling variance (the quantity of interest) from the
#' residual within-bootstrap imputation noise.
#'
#' Unlike Rubin's rules, Boot MI does not rely on the congeniality of the
#' imputation and analysis models, so it delivers nominal-coverage intervals
#' even when (as is typical for causal estimands) the two are uncongenial.
#'
#' @param imp A `mids` object. Its incomplete data (`imp$data`), imputation
#'   `method`, and `predictorMatrix` are reused to re-impute each bootstrap
#'   resample.
#' @param fit_args,contrast_args Argument bundles for `mi_fit_one()`.
#' @param conf_level Numeric confidence level.
#' @param B Integer number of bootstrap resamples.
#' @param M Integer number of imputations per resample (>= 2 so the
#'   within-bootstrap variance component is estimable).
#' @param parallel Character backend: `"no"` or `"future"`.
#' @param seed Optional integer seed for reproducibility.
#' @param call The originating `causat_mice()` call.
#'
#' @returns A [causatr_result] with `ci_method = "boot_mi"`.
#'
#' @details
#' For parameter column with per-cell estimate \eqn{\theta_{b,j}} (bootstrap
#' \eqn{b}, imputation \eqn{j}), let \eqn{\bar\theta_b} be the per-bootstrap
#' mean and \eqn{\bar\theta} the grand mean. With
#' \eqn{B_{\mathrm{comp}} = \widehat{\mathrm{Var}}_b(\bar\theta_b)} and
#' \eqn{W = m^{-1}_{\mathrm{boot}} \sum_b \widehat{\mathrm{Var}}_j
#' (\theta_{b,j})}, the Boot MI variance is
#' \deqn{\hat{V} = B_{\mathrm{comp}} - W / M,}
#' floored at a small positive value, with \eqn{t_{B-1}} intervals around
#' \eqn{\bar\theta}.
#'
#' @references
#' von Hippel PT (2020). How many imputations do you need? *Sociological
#' Methods & Research* 49(3):699-718.
#'
#' Bartlett JW, Hughes RA (2020). Bootstrap inference for multiple imputation
#' under uncongeniality and misspecification. *Statistical Methods in Medical
#' Research* 29(12):3533-3546.
#'
#' @seealso [causat_mice()], [pool_rubin()]
#' @noRd
pool_boot_mi <- function(
  imp,
  fit_args,
  contrast_args,
  conf_level,
  B,
  M,
  parallel,
  seed,
  call
) {
  B <- as.integer(B)
  M <- as.integer(M)
  if (B < 2L) {
    rlang::abort(
      c(
        "`B` must be at least 2 for `pool_method = \"boot_mi\"`.",
        x = paste0("Got B = ", B, ".")
      ),
      class = "causatr_mi_bad_boot"
    )
  }
  if (M < 2L) {
    rlang::abort(
      c(
        "`M` must be at least 2 for `pool_method = \"boot_mi\"`.",
        x = paste0("Got M = ", M, "."),
        i = paste0(
          "Boot MI needs >= 2 imputations per bootstrap to estimate the ",
          "within-bootstrap variance component."
        )
      ),
      class = "causatr_mi_bad_boot"
    )
  }

  incomplete <- imp$data
  method <- imp$method
  pred <- imp$predictorMatrix
  n <- nrow(incomplete)

  # One bootstrap resample: draw rows with replacement, re-impute M times,
  # fit + contrast each completed dataset. Returns the list of mi_extract()
  # outputs (length <= M; failed imputations drop to NULL).
  boot_one <- function(b) {
    idx <- sample.int(n, n, replace = TRUE)
    boot_data <- incomplete[idx, , drop = FALSE]
    imp_b <- tryCatch(
      mice::mice(
        boot_data,
        m = M,
        method = method,
        predictorMatrix = pred,
        printFlag = FALSE
      ),
      error = function(e) NULL
    )
    if (is.null(imp_b)) {
      return(vector("list", 0L))
    }
    lapply(seq_len(M), function(j) {
      tryCatch(
        mi_extract(mi_fit_one(
          mice::complete(imp_b, j),
          fit_args,
          contrast_args
        )),
        error = function(e) NULL
      )
    })
  }

  boots <- run_boot_loop(boot_one, B, parallel, seed)

  # Keep only bootstraps with all M imputations fitted, so the random-effects
  # design stays balanced (von Hippel's closed-form assumes equal M per
  # bootstrap). Drop unbalanced/failed resamples and require enough survive.
  ok <- vapply(
    boots,
    function(bl) length(bl) == M && all(!vapply(bl, is.null, logical(1))),
    logical(1)
  )
  boots <- boots[ok]
  if (length(boots) < 2L) {
    rlang::abort(
      c(
        "Too few bootstrap resamples produced complete fits to estimate variance.",
        x = paste0(sum(ok), " of ", B, " resamples fully succeeded."),
        i = "Increase `B`, simplify the imputation model, or use `pool_method = \"rubin\"`."
      ),
      class = "causatr_mi_boot_failed"
    )
  }

  first <- boots[[1L]][[1L]]
  est_arr <- boot_to_array(boots, "est", length(first$est_labels))
  est_pool <- pool_table_boot(
    est_arr,
    scale = "linear",
    conf_level = conf_level
  )

  has_con <- length(first$con_labels) > 0L
  con_pool <- if (has_con) {
    con_scale <- if (identical(first$type, "difference")) "linear" else "log"
    con_arr <- boot_to_array(boots, "con", length(first$con_labels))
    pool_table_boot(con_arr, scale = con_scale, conf_level = conf_level)
  } else {
    NULL
  }

  estimates <- build_pooled_table(
    labels = first$est_labels,
    label_col = first$term_col,
    by = first$est_by,
    pool = est_pool
  )
  contrasts <- if (!is.null(con_pool)) {
    build_pooled_table(
      labels = first$con_labels,
      label_col = "comparison",
      by = first$con_by,
      pool = con_pool
    )
  } else {
    data.table::data.table(
      comparison = character(0),
      estimate = numeric(0),
      se = numeric(0),
      ci_lower = numeric(0),
      ci_upper = numeric(0)
    )
  }

  vcov <- diag(est_pool$se^2, nrow = length(est_pool$se))
  result <- new_causatr_result(
    estimates = estimates,
    contrasts = contrasts,
    type = first$type,
    estimand = first$estimand,
    ci_method = "boot_mi",
    reference = first$reference,
    interventions = first$interventions,
    n = n,
    estimator = first$estimator,
    family = first$family,
    fit_type = first$fit_type,
    vcov = vcov,
    boot_t = NULL,
    boot_info = NULL,
    call = call
  )
  attr(result, "mi_details") <- list(
    pool_method = "boot_mi",
    B = length(boots),
    M = M,
    estimates = est_pool$diagnostics,
    contrasts = if (!is.null(con_pool)) con_pool$diagnostics else NULL
  )
  result
}

#' Run the bootstrap loop, optionally in parallel
#'
#' @param boot_one Function of the bootstrap index `b` returning that
#'   resample's per-imputation extractions.
#' @param B Integer number of resamples.
#' @param parallel `"no"` or `"future"`.
#' @param seed Optional integer seed.
#' @returns A length-`B` list of `boot_one()` outputs.
#' @noRd
run_boot_loop <- function(boot_one, B, parallel, seed) {
  if (parallel == "future") {
    check_pkg("future.apply")
    # future.seed = TRUE generates parallel-safe RNG streams; the user's
    # `seed` makes the whole stream reproducible.
    if (!is.null(seed)) {
      set.seed(seed)
    }
    return(future.apply::future_lapply(
      seq_len(B),
      boot_one,
      future.seed = TRUE
    ))
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }
  lapply(seq_len(B), boot_one)
}

#' Reshape surviving bootstraps into a `B x M x p` array
#'
#' @param boots List of surviving bootstraps, each a length-`M` list of
#'   `mi_extract()` outputs.
#' @param field Either `"est"` (means) or `"con"` (contrasts).
#' @param p Integer number of parameter columns.
#' @returns A numeric array with dims `c(B, M, p)`.
#' @noRd
boot_to_array <- function(boots, field, p) {
  B <- length(boots)
  M <- length(boots[[1L]])
  arr <- array(NA_real_, dim = c(B, M, p))
  for (b in seq_len(B)) {
    for (j in seq_len(M)) {
      arr[b, j, ] <- boots[[b]][[j]][[field]]
    }
  }
  arr
}

#' Apply the von Hippel one-way random-effects decomposition per column
#'
#' @param arr Numeric `B x M x p` array of per-cell estimates.
#' @param scale `"linear"` or `"log"` (pool log-estimates, exponentiate back).
#' @param conf_level Numeric confidence level.
#' @returns A list with `estimate`, `se`, `df`, `ci_lower`, `ci_upper`, and a
#'   `diagnostics` data.table, mirroring `pool_table_rubin()`'s shape.
#' @noRd
pool_table_boot <- function(arr, scale, conf_level) {
  d <- dim(arr)
  B <- d[1L]
  M <- d[2L]
  p <- d[3L]
  alpha <- (1 - conf_level) / 2
  crit <- stats::qt(1 - alpha, df = B - 1L)

  out <- lapply(seq_len(p), function(j) {
    cell <- arr[,, j]
    if (identical(scale, "log")) {
      cell <- log(cell)
    }
    boot_means <- rowMeans(cell)
    grand <- mean(boot_means)
    # Between-bootstrap variance of the resample means; within-bootstrap
    # imputation variance averaged across resamples.
    b_comp <- stats::var(boot_means)
    w_comp <- mean(apply(cell, 1L, stats::var))
    # von Hippel: subtract the imputation-noise share carried by each
    # bootstrap mean. Floor at a small positive value to keep CIs finite
    # when imputation noise dominates a small between-bootstrap signal.
    v_hat <- max(b_comp - w_comp / M, .Machine$double.eps)

    lo_link <- grand - crit * sqrt(v_hat)
    hi_link <- grand + crit * sqrt(v_hat)
    if (identical(scale, "log")) {
      est <- exp(grand)
      se <- est * sqrt(v_hat)
      lo <- exp(lo_link)
      hi <- exp(hi_link)
    } else {
      est <- grand
      se <- sqrt(v_hat)
      lo <- lo_link
      hi <- hi_link
    }
    list(
      estimate = est,
      se = se,
      df = B - 1L,
      ci_lower = lo,
      ci_upper = hi,
      between = b_comp,
      within = w_comp
    )
  })

  pull <- function(field) vapply(out, `[[`, numeric(1), field)
  diagnostics <- data.table::data.table(
    between = pull("between"),
    within = pull("within"),
    df = pull("df")
  )
  list(
    estimate = pull("estimate"),
    se = pull("se"),
    df = pull("df"),
    ci_lower = pull("ci_lower"),
    ci_upper = pull("ci_upper"),
    diagnostics = diagnostics
  )
}
