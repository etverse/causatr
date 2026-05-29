#' Apply an intervention to a treatment vector (no data-frame mutation)
#'
#' @description
#' Lightweight companion to `apply_single_intervention()` in
#' `R/interventions.R`. The difference: this function takes a
#' treatment *vector* and returns a new vector, without copying or
#' mutating any `data.table`. It's what `compute_density_ratio_weights()`
#' and `make_weight_fn()` call to get `a_int` for density evaluation;
#' they cannot use `apply_single_intervention()` directly because that
#' function mutates a `data.table` column in place, which is wasteful
#' inside a `numDeriv::jacobian` loop (we rebuild the closure's `a_int`
#' once at creation time and then the closure itself never touches the
#' data.table).
#'
#' @param intervention A `causatr_intervention` object.
#' @param data `data.table` passed to `dynamic()` rules (only used by
#'   the dynamic branch).
#' @param a_obs Numeric / character / factor vector of observed
#'   treatment values.
#'
#' @return A vector of the same length and type (with the caveats
#'   documented below for dynamic / character interventions).
#'
#' @noRd
apply_intervention_to_values <- function(intervention, data, a_obs) {
  iv <- intervention
  switch(
    iv$type,
    static = {
      # `rep(value, length)` preserves scalar type -- `static("chemo")`
      # returns a character vector, `static(1)` returns a double.
      rep(iv$value, length(a_obs))
    },
    shift = {
      # Only meaningful for numeric treatments. The upstream
      # `check_intervention_family_compat()` has already rejected
      # non-numeric cases (shift on a factor treatment would
      # produce NaN silently).
      a_obs + iv$delta
    },
    scale = a_obs * iv$factor,
    threshold = pmax(pmin(a_obs, iv$upper), iv$lower),
    dynamic = {
      # User-supplied rule. We mirror the type-preservation checks
      # from `apply_single_intervention()` so the density evaluator
      # downstream doesn't get a wrong-type vector. Replicating the
      # checks (rather than calling `apply_single_intervention()`)
      # avoids the unnecessary data.table copy. The error messages
      # point at the same sharp edges -- factor coercion, unknown
      # levels, wrong length.
      new_trt <- iv$rule(data, a_obs)
      if (length(new_trt) != length(a_obs)) {
        rlang::abort(
          paste0(
            "`dynamic()` rule must return a vector of length ",
            length(a_obs),
            " (got ",
            length(new_trt),
            ")."
          )
        )
      }
      if (is.numeric(a_obs) && !is.numeric(new_trt)) {
        rlang::abort(
          "`dynamic()` rule returned non-numeric for a numeric treatment column."
        )
      }
      if (is.factor(a_obs)) {
        if (is.character(new_trt)) {
          unknown <- setdiff(unique(new_trt), levels(a_obs))
          if (length(unknown) > 0L) {
            rlang::abort(
              paste0(
                "`dynamic()` rule returned level(s) not in the factor ",
                "treatment: ",
                paste(shQuote(unknown), collapse = ", ")
              )
            )
          }
          new_trt <- factor(new_trt, levels = levels(a_obs))
        } else if (is.factor(new_trt)) {
          if (!identical(levels(new_trt), levels(a_obs))) {
            rlang::abort(
              "`dynamic()` rule returned a factor with mismatched levels."
            )
          }
        } else {
          rlang::abort(
            "`dynamic()` rule returned a non-factor, non-character value for a factor treatment."
          )
        }
      }
      new_trt
    },
    ipsi = {
      # IPSI callers should go through the closed-form branch in
      # `compute_density_ratio_weights()` / `make_weight_fn()`
      # directly -- reaching here is a bug in the caller.
      rlang::abort(
        "Internal error: `apply_intervention_to_values()` should not be called with an IPSI intervention."
      )
    },
    rlang::abort(paste0("Unknown intervention type: '", iv$type, "'."))
  )
}


#' IPSI closed-form density-ratio weight
#'
#' @description
#' Kennedy (2019) incremental propensity score intervention weight:
#' \deqn{w_i = \frac{\delta A_i + (1 - A_i)}{\delta p_i + (1 - p_i)}.}
#'
#' Derivation. Under IPSI the intervened density is
#' \eqn{f_\delta(a \mid L) = (\delta p)^a (1 - p)^{1-a} / (1 - p + \delta p)}.
#' The density ratio simplifies to
#' \eqn{f_\delta(A_i \mid L_i) / f(A_i \mid L_i) = (\delta A_i + (1 - A_i))
#'  / (\delta p_i + (1 - p_i))}. No density evaluation at a
#' counterfactual treatment value is needed -- the intervention acts
#' directly on the propensity.
#'
#' Positivity: IPSI always preserves positivity (finite `delta`
#' cannot push `p` to 0 or 1), so the denominator is strictly between
#' `min(1, delta)` and `max(1, delta)`. No guard needed.
#'
#' @param a_obs Numeric binary 0/1 treatment vector.
#' @param p Numeric vector of predicted treatment probabilities on the
#'   same rows.
#' @param delta Positive scalar odds multiplier.
#'
#' @return Numeric weight vector, same length as `a_obs`.
#'
#' @noRd
ipsi_weight_formula <- function(a_obs, p, delta) {
  (delta * a_obs + (1 - a_obs)) / (delta * p + (1 - p))
}


#' Bayes-rule numerator f*(L) for the HT estimand weight
#'
#' @description
#' Returns the per-individual multiplier `f_star_i = f(A* | L_i)` that
#' converts the ATE density-ratio weight into an ATT or ATC weight:
#'
#' \deqn{w_i = \mathbb 1\{A_i = a\} \cdot f^\*_i / f(a \mid L_i).}
#'
#' The derivation is the standard Bayes-rule rewrite of
#' \eqn{E[Y^a \mid A = A^\*]} (Imbens 2004; Hernan & Robins Ch. 12).
#' For ATE the target is the whole population and \eqn{f^\* \equiv 1};
#' for ATT the target is the treated and \eqn{f^\*_i = p(L_i)}; for
#' ATC the target is the controls and \eqn{f^\*_i = 1 - p(L_i)}.
#'
#' Only the Bernoulli treatment family is supported because
#' `check_estimand_intervention_compat()` has already rejected ATT /
#' ATC for non-binary static interventions -- hitting the fallback
#' `rlang::abort` here would indicate a missed upstream guard.
#'
#' @param estimand Character scalar in `c("ATE", "ATT", "ATC")`.
#' @param treatment_model A `causatr_treatment_model`.
#' @param fit_data The `fit_rows`-subset `data.table` the caller is
#'   building weights for.
#' @param family_tag Character. The treatment family tag from
#'   `treatment_model$family`.
#'
#' @return Numeric vector of length `nrow(fit_data)` (or a length-1
#'   vector of `1` for ATE, which `R`'s recycling promotes to the
#'   correct per-row constant without allocating).
#'
#' @noRd
ht_bayes_numerator <- function(
  estimand,
  treatment_model,
  fit_data,
  family_tag
) {
  if (estimand == "ATE") {
    return(1)
  }
  if (family_tag != "bernoulli") {
    rlang::abort(
      paste0(
        "Internal error: ATT / ATC Bayes numerator requested for a non-",
        "Bernoulli treatment family ('",
        family_tag,
        "'). `check_estimand_intervention_compat()` should have ",
        "rejected this upstream."
      )
    )
  }
  p <- unname(stats::predict(
    treatment_model$model,
    newdata = fit_data,
    type = "response"
  ))
  if (estimand == "ATT") {
    return(p)
  }
  if (estimand == "ATC") {
    return(1 - p)
  }
  rlang::abort(
    paste0("Internal error: unknown estimand '", estimand, "'.")
  )
}


#' Check that an intervention is compatible with the treatment family
#'
#' @description
#' Rejects nonsensical (intervention, treatment family) combinations
#' up front, with actionable error messages pointing at the right
#' alternative. The compatibility matrix the self-contained IPW
#' engine supports is:
#'
#' | intervention    | bernoulli | gaussian         | categorical | poisson / negbin    |
#' |-----------------|-----------|------------------|-------------|---------------------|
#' | `static()`      | yes HT      | no                | yes HT        | no                   |
#' | `shift()`       | no         | yes smooth ratio   | no           | yes integer only      |
#' | `scale_by()`    | no         | yes smooth ratio   | no           | yes integer-preserving|
#' | `threshold()`   | no         | no (use gcomp)    | no           | no (use gcomp)       |
#' | `dynamic()`     | yes HT      | no (use gcomp)    | yes HT        | no (use gcomp)       |
#' | `ipsi()`        | yes Kennedy | no                | no           | no                   |
#'
#' @param intervention A `causatr_intervention`.
#' @param treatment_model A `causatr_treatment_model`.
#'
#' @return `NULL` invisibly; aborts on invalid combinations.
#'
#' @noRd
check_intervention_family_compat <- function(
  intervention,
  treatment_model
) {
  iv_type <- intervention$type
  fam <- treatment_model$family
  is_count <- fam %in% c("poisson", "negbin")

  if (iv_type == "stochastic") {
    rlang::abort(
      c(
        paste0(
          "`stochastic()` interventions are only supported ",
          "under `estimator = 'gcomp'`."
        ),
        i = paste0(
          "Stochastic interventions require Monte Carlo ",
          "integration over the counterfactual treatment ",
          "distribution, which is only available via the ",
          "g-formula (outcome-model prediction)."
        ),
        i = paste0(
          "Use `causat(..., estimator = 'gcomp')`, or ",
          "rewrite the intervention as a deterministic ",
          "`dynamic()` rule if the policy is actually ",
          "deterministic."
        )
      )
    )
  }

  if (iv_type == "ipsi" && fam != "bernoulli") {
    rlang::abort(
      c(
        "`ipsi()` interventions are only defined for binary (0/1) treatments.",
        i = paste0(
          "The treatment column is classified as '",
          fam,
          "'. Use `shift()` or `scale_by()` for a continuous MTP instead."
        )
      )
    )
  }

  # shift / scale_by require a numeric treatment with a well-defined
  # density: gaussian or count (poisson / negbin).
  if (
    iv_type %in%
      c("shift", "scale") &&
      !fam %in% c("gaussian", "poisson", "negbin")
  ) {
    rlang::abort(
      c(
        paste0(
          "`",
          iv_type,
          "()` interventions require a numeric continuous treatment."
        ),
        i = paste0(
          "The treatment column is classified as '",
          fam,
          "'. Use `static()` for binary / categorical treatments, ",
          "or pass `propensity_family = 'poisson'` / `'negbin'` for count treatments."
        )
      )
    )
  }

  # threshold on any non-gcomp-compatible family
  if (iv_type == "threshold" && fam != "gaussian" && !is_count) {
    rlang::abort(
      c(
        paste0(
          "`threshold()` interventions require a numeric continuous treatment."
        ),
        i = paste0(
          "The treatment column is classified as '",
          fam,
          "'. Use `static()` for binary / categorical treatments."
        )
      )
    )
  }

  # Count-specific guards: shift must be integer, scale must preserve
  # integer support, threshold and dynamic are rejected.
  if (is_count) {
    if (iv_type == "static") {
      rlang::abort(
        c(
          paste0(
            "`static(v)` on a count treatment is degenerate for IPW."
          ),
          i = "The Horvitz-Thompson indicator weight is zero for almost all observations.",
          i = "Use `shift()` for integer shifts, or switch to `estimator = 'gcomp'`."
        )
      )
    }
    if (iv_type == "threshold") {
      rlang::abort(
        c(
          "`threshold()` on a count treatment is not supported by the IPW engine.",
          i = "The pushforward of a count density under a boundary clamp is a mixed measure.",
          i = "Use `estimator = 'gcomp'`."
        )
      )
    }
    if (iv_type == "dynamic") {
      rlang::abort(
        c(
          "`dynamic()` rules on count treatments are not supported by the IPW engine.",
          i = "Use `shift()` for integer shifts, or switch to `estimator = 'gcomp'` for deterministic rules."
        )
      )
    }
    if (iv_type == "shift") {
      delta <- intervention$delta
      if (delta != round(delta)) {
        rlang::abort(
          c(
            paste0(
              "`shift(",
              delta,
              ")` is not integer-valued."
            ),
            i = paste0(
              "Count treatments (",
              fam,
              ") require integer shift deltas because `dpois()` / `dnbinom()` return 0 at non-integer arguments."
            ),
            i = "Use `estimator = 'gcomp'` for fractional shifts on count data."
          )
        )
      }
    }
    if (iv_type == "scale") {
      fct <- intervention$factor
      a_obs <- treatment_model$model$data[[treatment_model$treatment]]
      if (is.null(a_obs)) {
        a_obs <- stats::model.response(stats::model.frame(
          treatment_model$model
        ))
      }
      # The density ratio evaluates dpois(a_obs / factor, lambda).
      # For dpois to return a non-zero value, a_obs / factor must be
      # a non-negative integer for every observed treatment value.
      inv_scaled <- a_obs / fct
      if (!all(inv_scaled == round(inv_scaled) & inv_scaled >= 0)) {
        rlang::abort(
          c(
            paste0(
              "`scale_by(",
              fct,
              ")` does not produce integer inverse values (A / ",
              fct,
              ") for all observed treatment values."
            ),
            i = paste0(
              "Count treatments (",
              fam,
              ") require that A / factor is a non-negative integer for every observation, because the density ratio evaluates the pmf at A / factor."
            ),
            i = "Use `estimator = 'gcomp'` for non-integer-preserving scales on count data."
          )
        )
      }
    }
  }

  if (iv_type == "static" && fam == "gaussian") {
    rlang::abort(
      c(
        "`static(v)` on a continuous treatment is degenerate for IPW.",
        i = "No observations lie exactly at `v`, so the Horvitz-Thompson weight is zero almost surely.",
        i = "Use `shift()` or `scale_by()` to move the whole treatment distribution, or switch to `estimator = 'gcomp'`."
      )
    )
  }

  if (iv_type == "threshold" && fam == "gaussian") {
    rlang::abort(
      c(
        "`threshold(lo, hi)` on a continuous treatment is not supported by the IPW engine.",
        i = "The pushforward of a continuous density under a boundary clamp has point masses at `lo` and `hi`, so the density ratio w.r.t. the fitted `f(a|l)` is not well-defined.",
        i = "Use `estimator = 'gcomp'` -- the clamped-treatment counterfactual is well-defined there via predict-then-average on the outcome model."
      )
    )
  }

  if (iv_type == "dynamic" && fam == "gaussian") {
    rlang::abort(
      c(
        "`dynamic()` rules on continuous treatments are not supported by the IPW engine.",
        i = "A deterministic per-individual target is a Dirac per individual; the density ratio is degenerate.",
        i = "Use `shift()` / `scale_by()` for smooth MTPs, or switch to `estimator = 'gcomp'` for deterministic rules."
      )
    )
  }

  invisible(NULL)
}
