#' Print a causatr fit
#'
#' @description
#' Displays a compact summary of a [causatr_fit][causat] object, showing the
#' causal estimator, treatment type, outcome, treatment variable, and sample
#' size.
#'
#' @param x A `causatr_fit` object.
#' @param ... Currently unused.
#' @return Invisibly returns `x`.
#' @seealso [summary.causatr_fit()], [causat()]
#' @export
print.causatr_fit <- function(x, ...) {
  # Pretty-print the estimator. Defaults-through the raw value if the
  # user has somehow set an unrecognized estimator so print() never
  # errors on a malformed fit object.
  estimator_label <- switch(
    x$estimator,
    gcomp = "G-computation",
    ipw = "IPW",
    aipw = "AIPW (doubly-robust)",
    matching = "Matching (MatchIt)",
    snm = "SNM (g-estimation)",
    x$estimator
  )
  family_label <- format_family(x$family)
  # Collapse multivariate treatments into a comma-separated label.
  trt_label <- paste(x$treatment, collapse = ", ")

  cat("<causatr_fit>\n")
  cat(" Estimator:  ", estimator_label, "\n", sep = "")
  cat(" Type:       ", x$type, "\n", sep = "")
  cat(" Outcome:    ", x$outcome, " (", family_label, ")\n", sep = "")
  cat(" Treatment:  ", trt_label, "\n", sep = "")
  cat(" Estimand:   ", x$estimand, "\n", sep = "")
  # Per-component confounder overrides: when any `confounders_*` slot is

  # non-NULL, display the resolved per-component formulas instead of the
  # legacy single `confounders` slot.
  has_overrides <- !is.null(x$confounders_outcome) ||
    !is.null(x$confounders_treatment) ||
    !is.null(x$confounders_censoring) ||
    !is.null(x$confounders_sampling)
  if (has_overrides) {
    cat(
      " Conf (outcome):  ",
      deparse(resolve_confounders_outcome(x)),
      "\n",
      sep = " "
    )
    cat(
      " Conf (treatment):",
      deparse(resolve_confounders_treatment(x)),
      "\n",
      sep = " "
    )
    if (!is.null(x$confounders_censoring)) {
      cat(
        " Conf (censoring):",
        deparse(resolve_confounders_censoring(x)),
        "\n",
        sep = " "
      )
    }
    if (!is.null(x$confounders_sampling)) {
      cat(
        " Conf (sampling): ",
        deparse(resolve_confounders_sampling(x)),
        "\n",
        sep = " "
      )
    }
  } else {
    cat(" Confounders:", deparse(x$confounders), "\n", sep = " ")
  }
  # TV confounders and id/time are only shown when relevant -- keeping
  # the display compact for point-treatment fits.
  has_tv_overrides <- !is.null(x$confounders_tv_outcome) ||
    !is.null(x$confounders_tv_treatment)
  if (has_tv_overrides) {
    cat(
      " TV conf (outcome):  ",
      deparse(resolve_confounders_tv_outcome(x)),
      "\n",
      sep = " "
    )
    cat(
      " TV conf (treatment):",
      deparse(resolve_confounders_tv_treatment(x)),
      "\n",
      sep = " "
    )
  } else if (!is.null(x$confounders_tv)) {
    cat(" TV conf.:  ", deparse(x$confounders_tv), "\n", sep = " ")
  }
  if (!is.null(x$id)) {
    cat(" ID:         ", x$id, "\n", sep = "")
    cat(" Time:       ", x$time, "\n", sep = "")
  }
  # Stratified ICE: surface the stratifying column when present. Guarded
  # so non-stratified fits print exactly as before (no snapshot churn).
  if (!is.null(x$details$stratified)) {
    cat(" Stratified: ", x$details$stratified, "\n", sep = "")
  }
  cat(" N:          ", nrow(x$data), "\n", sep = "")
  invisible(x)
}

#' Format a family object for display
#'
#' @param family A character string, family object, or function.
#' @return A character string.
#' @noRd
format_family <- function(family) {
  # Three-shape dispatch matching `resolve_family()` but producing
  # a display string rather than a family object. Used only for
  # print() output, so the fallback is "<custom>" rather than an
  # abort -- printing a fit should never error out.
  if (is.character(family)) {
    # Already a string like "gaussian" -- pass through.
    return(family)
  }
  if (is.list(family) && !is.null(family$family)) {
    # Already a family object (list with $family and $link slots):
    # format as "gaussian(identity)", matching base R's convention.
    return(paste0(family$family, "(", family$link, ")"))
  }
  if (is.function(family)) {
    # A family closure (e.g. `stats::binomial`) -- try to evaluate it.
    # Wrapped in tryCatch in case the closure needs arguments; if it
    # fails we fall through to "<custom>".
    fam <- tryCatch(family(), error = function(e) NULL)
    if (!is.null(fam) && !is.null(fam$family)) {
      return(paste0(fam$family, "(", fam$link, ")"))
    }
  }
  "<custom>"
}

#' Print a causatr result
#'
#' @description
#' Displays the causal estimator, contrast type, CI method, sample size,
#' intervention-specific marginal means, and pairwise contrasts (with SEs
#' and confidence intervals).
#'
#' @param x A `causatr_result` object.
#' @param ... Currently unused.
#' @return Invisibly returns `x`.
#' @seealso [summary.causatr_result()], [contrast()]
#' @export
print.causatr_result <- function(x, ...) {
  # Human-readable labels for the header block. Same defaults-through
  # pattern as `print.causatr_fit` in case of unrecognized slot values.
  estimator_label <- switch(
    x$estimator,
    gcomp = "G-computation",
    ipw = "IPW",
    aipw = "AIPW (doubly-robust)",
    matching = "Matching",
    snm = "SNM (g-estimation)",
    x$estimator
  )
  type_label <- switch(
    x$type,
    difference = "Difference",
    ratio = "Ratio",
    or = "Odds ratio"
  )

  cat("<causatr_result>\n")
  cat(" Estimator: ", estimator_label, "\n", sep = "")
  cat(" Estimand:  ", x$estimand, "\n", sep = "")
  cat(" Contrast:  ", type_label, "\n", sep = "")
  cat(" CI method: ", x$ci_method, "\n", sep = "")
  cat(" N:         ", x$n, "\n", sep = "")

  # Bootstrap diagnostics: surface the success/failure split so users
  # see at print time when replicates were silently discarded. The
  # `boot_info` slot is only populated for `ci_method = "bootstrap"`,
  # and may be a single 3-element list (no `by`) or a list-of-lists
  # keyed by subgroup level when `by = ...` was used. In the by case
  # we reduce to totals across strata; per-stratum detail is available
  # in `x$boot_info` for users who want it.
  if (!is.null(x$boot_info)) {
    bi <- x$boot_info
    if (
      is.list(bi) && !all(c("n_requested", "n_ok", "n_fail") %in% names(bi))
    ) {
      # by-stratum list-of-lists: collapse to totals. `n_requested` carries
      # the user's `n_boot`, which is a double when passed as a bare numeric
      # (`n_boot = 40`), so coerce before the integer-typed vapply -- otherwise
      # the collapse aborts on any by-stratified bootstrap result.
      n_req <- sum(vapply(
        bi,
        function(b) as.integer(b$n_requested %||% 0L),
        integer(1)
      ))
      n_fail <- sum(vapply(
        bi,
        function(b) as.integer(b$n_fail %||% 0L),
        integer(1)
      ))
    } else {
      n_req <- bi$n_requested
      n_fail <- bi$n_fail
    }
    if (!is.null(n_req) && n_req > 0L) {
      n_ok <- n_req - n_fail
      pct <- if (n_fail > 0L) {
        paste0(" (", round(100 * n_fail / n_req, 1), "% failed)")
      } else {
        ""
      }
      cat(
        " Bootstrap: ",
        n_ok,
        "/",
        n_req,
        " replicates",
        pct,
        "\n",
        sep = ""
      )
    }
  }

  # Tier-2 fallback flag: when the sandwich path could neither recover
  # the full IF via `sandwich::estfun()` nor compute analytic bread,
  # `variance_if_numeric()` drops the IF cross-term and tags the vcov
  # with `tier2_approximate = TRUE`. Surface that here so the reader
  # doesn't assume the SE is the exact asymptotic value.
  if (isTRUE(attr(x$vcov, "tier2_approximate"))) {
    cat(
      " Note:      sandwich SE is approximate (Tier-2 fallback; use ci_method='bootstrap' for exact)\n",
      sep = ""
    )
  }

  # SNM results use a different layout: the estimates table contains
  # blip parameters (column "parameter"), not intervention means.
  is_snm <- identical(x$estimator, "snm")

  if (is_snm) {
    is_avg_blip <- "parameter" %in%
      names(x$estimates) &&
      nrow(x$estimates) == 1L &&
      identical(x$estimates$parameter[1L], "avg_blip_effect")
    if (is_avg_blip) {
      cat("\nAveraged blip effect:\n")
    } else if (identical(x$fit_type, "longitudinal")) {
      cat("\nPer-stage blip parameters:\n")
    } else {
      cat("\nBlip parameters:\n")
    }
    print(x$estimates, digits = 3)
    # Only show contrasts if non-empty (Path B produces a row).
    if (nrow(x$contrasts) > 0L) {
      cat("\nContrasts:\n")
      print(x$contrasts, digits = 3)
    }
    invisible(x)
  } else {
    # When `by = ...` was used, `contrast()` adds a `by` column to both
    # `x$estimates` and `x$contrasts`. Detecting it via column presence
    # (rather than a separate slot) lets `print()` adapt without requiring
    # callers to pass any extra flag.
    has_by <- "by" %in% names(x$estimates)
    # A multinomial outcome carries a `class` column: the table holds
    # P(Y = class | do(A = a)) per intervention rather than a single mean.
    has_class <- "class" %in% names(x$estimates)

    means_header <- if (has_by) {
      "\nIntervention means (by subgroup):\n"
    } else if (has_class) {
      "\nClass probabilities by intervention:\n"
    } else {
      "\nIntervention means:\n"
    }
    cat(means_header)
    print(x$estimates, digits = 3)

    contrasts_header <- if (has_by) {
      "\nContrasts (by subgroup):\n"
    } else if (has_class) {
      "\nContrasts (per class):\n"
    } else {
      "\nContrasts:\n"
    }
    cat(contrasts_header)
    print(x$contrasts, digits = 3)
    invisible(x)
  }
}

#' Print causatr diagnostics
#'
#' @description
#' Displays positivity summaries, covariate balance tables, weight
#' distributions (IPW), and match quality metrics (matching) from a
#' [causatr_diag][diagnose] object. Multi-panel diagnostics (built with
#' a non-`NULL` `interventions =` argument) print each panel in turn,
#' headed by its intervention key.
#'
#' @param x A `causatr_diag` object.
#' @param ... Currently unused.
#' @return Invisibly returns `x`.
#' @seealso [summary.causatr_diag()], [diagnose()]
#' @export
print.causatr_diag <- function(x, ...) {
  cat("<causatr_diag>\n", " Estimator:", x$estimator, "\n", sep = "")
  fit_info <- x$fit_info
  # `fit_info` is absent in legacy `causatr_diag` objects (e.g. those
  # constructed by hand in older user scripts), so guard each accessor
  # before printing.
  if (!is.null(fit_info$treatment_type)) {
    cat(" Treatment: ", fit_info$treatment_type, "\n", sep = "")
  }
  if (!is.null(fit_info$estimand)) {
    cat(" Estimand:  ", fit_info$estimand, "\n", sep = "")
  }
  if (identical(fit_info$type, "longitudinal")) {
    cat(" Type:      longitudinal\n")
  }
  cat("\n")

  panels <- x$per_intervention
  if (is.null(panels) || length(panels) == 0L) {
    # Defensive fall-through for legacy objects: use the flat slots.
    print_diag_panel(
      list(
        positivity = x$positivity,
        balance = x$balance,
        weights = x$weights
      ),
      header = NULL
    )
  } else if (length(panels) == 1L) {
    # Single-panel layout matches the pre-chunk-11a flat output: no
    # intervention header, just the three sub-tables. Existing
    # `expect_output(print(diag), "Positivity")` style tests continue
    # to pass without change.
    print_diag_panel(panels[[1L]], header = NULL)
  } else {
    # Multi-panel layout: section by intervention key. Each section
    # is offset by a blank line so the intervention break is visible
    # without ANSI styling (the package targets headless test runs
    # too, so we keep the formatting plain).
    nms <- names(panels)
    for (i in seq_along(panels)) {
      print_diag_panel(panels[[i]], header = nms[i])
      if (i < length(panels)) {
        cat("\n")
      }
    }
  }

  if (!is.null(x$match_quality)) {
    cat("Match quality:\n")
    print(x$match_quality, row.names = FALSE)
    cat("\n")
  }
  invisible(x)
}

#' Print a single diagnostic panel
#'
#' @param panel A list with optional `positivity` / `balance` / `weights`
#'   slots. Each is printed only if non-`NULL`.
#' @param header Character scalar intervention key, or `NULL` for the
#'   single-panel layout where no intervention header is shown.
#' @return `NULL` invisibly.
#' @noRd
print_diag_panel <- function(panel, header = NULL) {
  if (!is.null(header)) {
    cat("Intervention: ", header, "\n", sep = "")
  }
  if (!is.null(panel$positivity)) {
    pos <- panel$positivity
    if (is.list(pos) && !data.table::is.data.table(pos)) {
      # Named list: per-component (multivariate) or per-period
      # (longitudinal) tables.
      for (nm in names(pos)) {
        cat("Positivity (", nm, "):\n", sep = "")
        print(pos[[nm]], row.names = FALSE)
        cat("\n")
      }
    } else {
      cat("Positivity:\n")
      print(pos, row.names = FALSE)
      cat("\n")
    }
  }
  if (!is.null(panel$balance)) {
    bal <- panel$balance
    if (
      is.list(bal) &&
        !data.table::is.data.table(bal) &&
        !inherits(bal, "bal.tab")
    ) {
      # Per-period balance (longitudinal).
      for (nm in names(bal)) {
        cat("Covariate balance (", nm, "):\n", sep = "")
        print(bal[[nm]])
        cat("\n")
      }
    } else {
      cat("Covariate balance:\n")
      print(bal)
      cat("\n")
    }
  }
  if (!is.null(panel$weights)) {
    wts <- panel$weights
    if (is.list(wts) && !data.table::is.data.table(wts)) {
      for (nm in names(wts)) {
        cat("Weight distribution (", nm, "):\n", sep = "")
        print(wts[[nm]], row.names = FALSE)
        cat("\n")
      }
    } else {
      cat("Weight distribution:\n")
      print(wts, row.names = FALSE)
      cat("\n")
    }
  }
  if (!is.null(panel$censoring)) {
    cens <- panel$censoring
    if (is.list(cens) && !data.table::is.data.table(cens)) {
      for (nm in names(cens)) {
        cat("Censoring (", nm, "):\n", sep = "")
        print(cens[[nm]], row.names = FALSE)
        cat("\n")
      }
    } else {
      cat("Censoring:\n")
      print(cens, row.names = FALSE)
      cat("\n")
    }
  }
  if (!is.null(panel$sampling)) {
    cat("Sampling model:\n")
    print(panel$sampling, row.names = FALSE)
    cat("\n")
  }
  if (!is.null(panel$pct_intervened)) {
    pct <- panel$pct_intervened
    if (is.list(pct) && !data.table::is.data.table(pct)) {
      for (nm in names(pct)) {
        cat("Feasibility (", nm, "):\n", sep = "")
        print(pct[[nm]], row.names = FALSE)
        cat("\n")
      }
    } else {
      cat("Feasibility:\n")
      print(pct, row.names = FALSE)
      cat("\n")
    }
  }
  invisible(NULL)
}
