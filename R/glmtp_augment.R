#' Discrete treatment support for a natural-history MTP
#'
#' @description
#' Returns the sorted unique discrete support \eqn{\mathcal{A}} of the
#' treatment column, used to enumerate the natural-treatment-history
#' sequences \eqn{\bar{s}_t} that index the augmented sequential regression
#' (Diaz, Williams, Morzywolek & Rudolph 2026, arXiv:2605.24167). The
#' augmented-data g-computation requires a *discrete* exposure: each original
#' row is replicated once per possible history sequence, so a continuous
#' treatment (with uncountable support) has no finite augmentation.
#'
#' Treatment values are accepted only when every observed value is a whole
#' number (binary / ordinal / integer-coded). Non-integer numeric values are
#' the signature of a continuous treatment and are rejected with the classed
#' error `causatr_glmtp_continuous_trt`, pointing the user to the standard
#' `shift()` / `scale_by()` MTPs that *do* support continuous exposures.
#'
#' @param data data.table with the (already lag-expanded) person-period rows.
#' @param treatment Character scalar. Treatment column name. Multivariate
#'   treatment is not supported here (the paper's augmentation is defined for
#'   a scalar discrete exposure) and is rejected upstream.
#' @param censoring Character scalar naming the censoring indicator column, or
#'   `NULL`. Censored rows are excluded from the support scan because their
#'   treatment values do not enter any per-step model.
#' @param call Caller environment for error reporting.
#'
#' @return A sorted numeric vector -- the discrete support \eqn{\mathcal{A}}.
#'
#' @references
#' Diaz I, Williams NT, Morzywolek P, Rudolph KE (2026). Modified treatment
#' policies that depend on the natural history of treatment. arXiv:2605.24167.
#'
#' @seealso [glmtp_enumerate_labels()], [glmtp_check_tractable()]
#' @family glmtp
#' @noRd
glmtp_support <- function(
  data,
  treatment,
  censoring = NULL,
  call = rlang::caller_env()
) {
  if (length(treatment) != 1L) {
    rlang::abort(
      paste0(
        "Natural-history MTPs (`grace_period()` / `carry_forward()` / ",
        "`cap_escalation()`) support a single discrete treatment column; got ",
        length(treatment),
        " (multivariate treatment)."
      ),
      class = "causatr_glmtp_continuous_trt",
      call = call
    )
  }

  # Drop censored rows before scanning the support: their treatment values
  # never enter a per-step regression, so including them could spuriously
  # widen the support with placeholder codes.
  uncens <- is_uncensored(data, censoring)
  vals <- data[[treatment]][uncens]
  vals <- vals[!is.na(vals)]

  if (length(vals) == 0L) {
    rlang::abort(
      "Treatment column has no non-missing values to form a support.",
      class = "causatr_glmtp_continuous_trt",
      call = call
    )
  }

  # Factor / character treatments would need a categorical policy contract
  # that the current `grace_period()` / `carry_forward()` constructors (which
  # are numeric by definition) do not provide. Reject early with the same
  # discreteness class so the message points to the supported exposure shape.
  if (!is.numeric(vals)) {
    rlang::abort(
      paste0(
        "Natural-history MTPs require a numeric discrete (binary / ordinal / ",
        "integer-coded) treatment; got ",
        class(vals)[1L],
        "."
      ),
      class = "causatr_glmtp_continuous_trt",
      call = call
    )
  }

  # A continuous treatment is identified by non-integer values: the
  # augmented-data construction enumerates every sequence in the product
  # support, which is only finite for a discrete exposure.
  if (!isTRUE(all.equal(vals, round(vals)))) {
    rlang::abort(
      c(
        "Natural-history MTPs require a discrete treatment.",
        i = paste0(
          "The treatment column contains non-integer values, the signature ",
          "of a continuous exposure. The augmented-data sequential regression ",
          "enumerates every natural-history sequence, which is finite only ",
          "for a discrete support."
        ),
        i = paste0(
          "Use `shift()` / `scale_by()` for continuous modified treatment ",
          "policies."
        )
      ),
      class = "causatr_glmtp_continuous_trt",
      call = call
    )
  }

  sort(unique(round(vals)))
}


#' Enumerate natural-treatment-history sequences of a given length
#'
#' @description
#' Builds the length-`t` natural-history support
#' \eqn{\bar{\mathcal{A}}_t = \mathcal{A}^t} as the full product of the
#' per-period support `support`. Each enumerated sequence \eqn{\bar{s}_t} is a
#' label that indexes one per-step regression \eqn{m_t(\bar{s}_t, \cdot,
#' \cdot)} in [glmtp_iterate()]. The empty label (`t = 0`) is the single
#' base-of-recursion sequence and is returned as a length-zero numeric vector.
#'
#' Sequences are ordered so the first coordinate varies slowest (the column
#' order of `expand.grid()`), giving a deterministic label ordering that the
#' model store, the pseudo-outcome merge, and the bootstrap refit all agree
#' on.
#'
#' @param support Numeric vector. The discrete treatment support from
#'   [glmtp_support()].
#' @param t Non-negative integer. The sequence length (number of periods in
#'   the history). `t = 0` returns the single empty label.
#'
#' @return A list of numeric vectors, each a length-`t` history sequence.
#'   Length \eqn{|\mathcal{A}|^t}.
#'
#' @seealso [glmtp_support()], [glmtp_label_key()]
#' @family glmtp
#' @noRd
glmtp_enumerate_labels <- function(support, t) {
  t <- as.integer(t)
  if (t <= 0L) {
    # The empty natural history: a single label carrying no values. This is
    # the base case for q_1, evaluated at the observed baseline treatment.
    return(list(numeric(0L)))
  }
  # `expand.grid` over `t` copies of the support enumerates the full product
  # A^t. Varying the first column slowest matches the chronological reading
  # of a history (s_1, ..., s_t).
  grid <- expand.grid(rep(list(support), t), KEEP.OUT.ATTRS = FALSE)
  lapply(seq_len(nrow(grid)), function(r) as.numeric(grid[r, , drop = TRUE]))
}


#' Stable character key for a natural-history sequence
#'
#' @description
#' Maps a numeric history sequence \eqn{\bar{s}_t} to a single character key
#' used to index the per-label model store and pseudo-outcome lists. The empty
#' sequence (base of the recursion) maps to `""`. Values are joined with `"|"`
#' so distinct sequences never collide (e.g. `c(1, 11)` vs `c(11, 1)`).
#'
#' @param seq Numeric vector. A history sequence (possibly length zero).
#'
#' @return A character scalar key.
#'
#' @seealso [glmtp_enumerate_labels()]
#' @family glmtp
#' @noRd
glmtp_label_key <- function(seq) {
  if (length(seq) == 0L) {
    return("")
  }
  paste0(seq, collapse = "|")
}


#' Guard the augmented-frame row / model blow-up
#'
#' @description
#' The augmented-data sequential regression fits \eqn{|\mathcal{A}|^t}
#' per-label models at backward step `t`, so the worst step
#' (\eqn{t = \tau - 1}) fits \eqn{|\mathcal{A}|^{\tau-1}} models and the
#' augmented frame reaches \eqn{n \cdot |\mathcal{A}|^{\tau-1}} rows. For a
#' large support or many periods this is intractable. This guard caps the
#' worst-step label count at `budget` and aborts with
#' `causatr_glmtp_too_many` (reporting the offending count) before any model
#' is fitted.
#'
#' @param support Numeric vector. The discrete treatment support.
#' @param n_times Positive integer. Number of time points \eqn{\tau}.
#' @param budget Positive integer. Maximum allowed worst-step label count
#'   \eqn{|\mathcal{A}|^{\tau-1}}. Default `1024`.
#' @param call Caller environment for error reporting.
#'
#' @return Invisibly the worst-step label count when within budget; otherwise
#'   aborts.
#'
#' @seealso [glmtp_support()], [glmtp_iterate()]
#' @family glmtp
#' @noRd
glmtp_check_tractable <- function(
  support,
  n_times,
  budget = 1024L,
  call = rlang::caller_env()
) {
  n_levels <- length(support)
  # The worst backward step is t = tau - 1, which carries |A|^(tau-1) labels
  # (q_tau is the base, label-independent). With a single period there is no
  # history to enumerate, so the count is trivially 1.
  worst_pow <- max(n_times - 1L, 0L)
  worst <- n_levels^worst_pow

  if (worst > budget) {
    rlang::abort(
      c(
        paste0(
          "Natural-history MTP augmentation is intractable for this problem."
        ),
        i = paste0(
          "Treatment support has ",
          n_levels,
          " levels over ",
          n_times,
          " periods, so the worst backward step enumerates ",
          n_levels,
          "^",
          worst_pow,
          " = ",
          worst,
          " history sequences (budget ",
          budget,
          ")."
        ),
        i = paste0(
          "Reduce the number of periods, coarsen the treatment support, or ",
          "raise the budget if the row blow-up is acceptable."
        )
      ),
      class = "causatr_glmtp_too_many",
      call = call
    )
  }
  invisible(worst)
}
