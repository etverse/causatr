#' Compute density-ratio weights for one intervention
#'
#' @description
#' Builds the per-individual weight vector that turns the observed
#' sample into the pseudo-population under which the weighted MSM
#' estimates the counterfactual mean \eqn{E[Y^d]}. Three weight
#' formulas are implemented, one per intervention "geometry":
#'
#' 1. **Horvitz–Thompson indicator form** — used for point-mass
#'    interventions on discrete treatments. The intervened
#'    distribution is a Dirac at a per-individual target value,
#'    so the Radon–Nikodym derivative w.r.t. the fitted density
#'    is
#'    \deqn{w_i = \frac{\mathbb 1\{A_i = d(A_i, L_i)\}}{f(A_i \mid L_i)}.}
#'    Rows whose observed treatment does not match the target
#'    contribute zero to the Hájek mean, as they should. Covers:
#'    - `static(value)` on binary / categorical treatments
#'    - `dynamic(rule)` on binary / categorical treatments
#'
#' 2. **Smooth density ratio** — used for MTPs on continuous
#'    treatments where both the fitted density and the intervened
#'    density are absolutely continuous w.r.t. Lebesgue measure.
#'    \deqn{w_i = \frac{f(d(A_i, L_i) \mid L_i)}{f(A_i \mid L_i)}.}
#'    Covers:
#'    - `shift(delta)` on continuous treatments
#'    - `scale_by(factor)` on continuous treatments
#'    - `threshold(lower, upper)` on continuous treatments
#'
#' 3. **IPSI closed form** — Kennedy (2019):
#'    \deqn{w_i = \frac{\delta A_i + (1 - A_i)}{\delta p_i + (1 - p_i)}.}
#'    Binary only. Does not use `evaluate_density()` — the ratio
#'    collapses to a direct function of the propensity \eqn{p_i}.
#'
#' 4. `NULL` (natural course) returns a vector of ones.
#'
#' ## Combinations that are deliberately rejected
#'
#' - `static(value)` on a continuous treatment: nobody is observed
#'   exactly at `value`, so the HT indicator is zero almost surely.
#'   Users who want "set the exposure to 3" on a continuous dose
#'   should use a smooth shift such as `shift(3 - mean(A))` or work
#'   with a grid of shifts instead.
#' - `dynamic(rule)` on a continuous treatment: a deterministic
#'   per-individual rule is a Dirac per individual, and the density
#'   ratio is degenerate. Users should use `estimator = "gcomp"` for
#'   continuous dynamic rules, or rewrite the rule as a smooth
#'   `shift()` / `scale_by()` / `threshold()`.
#'
#' The upstream `check_intervention_family_compat()` is the gate
#' that enforces this; the body below assumes its invariants hold.
#'
#' ## Row alignment
#'
#' The weight vector is length `sum(treatment_model$fit_rows)`. The
#' caller is responsible for aligning it to the MSM's own fit rows —
#' today those coincide because `fit_ipw()` uses the same
#' `get_fit_rows()` convention for both the propensity and the MSM.
#' If that ever changes, `fit_ipw()` must reconcile via the same
#' `na.action` path that `prepare_propensity_if_weightit()` uses.
#'
#' @param treatment_model A `causatr_treatment_model` from
#'   `fit_treatment_model()`.
#' @param data The full `data.table` passed to `causat()` (not the
#'   subsetted `fit_data`). Needed because `dynamic()` rules receive
#'   the whole frame, including time-varying covariates. The function
#'   handles the row subsetting internally via
#'   `treatment_model$fit_rows`.
#' @param intervention A `causatr_intervention` object, or `NULL` for
#'   natural course.
#'
#' @return Numeric vector of weights, length equal to the number of
#'   rows used by the treatment-density fit.
#'
#' @references
#' Kennedy EH (2019). Nonparametric causal effects based on incremental
#' propensity score interventions. *Journal of the American
#' Statistical Association* 114:645–656.
#'
#' Díaz I, Williams N, Hoffman KL, Schenck EJ (2023). Non-parametric
#' causal effects based on longitudinal modified treatment policies.
#' *JASA* 118:846–857.
#'
#' @noRd
compute_density_ratio_weights <- function(
  treatment_model,
  data,
  intervention
) {
  if (!inherits(treatment_model, "causatr_treatment_model")) {
    rlang::abort(
      "`treatment_model` must be a `causatr_treatment_model`."
    )
  }

  fit_rows <- treatment_model$fit_rows
  fit_data <- data[fit_rows]
  trt_col <- treatment_model$treatment
  a_obs <- fit_data[[trt_col]]
  n_fit <- length(a_obs)

  # Natural course: no intervention, d(A, L) = A, weights are all 1.
  # Hájek normalization then reproduces the unweighted marginal mean,
  # which is the correct "observed" baseline for MTP comparisons.
  if (is.null(intervention)) {
    return(rep(1, n_fit))
  }

  check_intervention_family_compat(intervention, treatment_model)

  iv_type <- intervention$type
  family_tag <- treatment_model$family

  # IPSI uses the closed-form shortcut. Going through the density-
  # ratio path would require "the intervened density f_new" which
  # IPSI does not have as a clean function of a counterfactual
  # treatment value — Kennedy's (2019) intervention acts on the
  # propensity, not on the treatment itself. The closed form
  # collapses the ratio directly. `unname()` strips the row-index
  # attribute that `stats::predict()` carries, so the returned
  # weight vector is unnamed — matching the other branches that go
  # through `ifelse()` / `dnorm()` which already drop names.
  if (iv_type == "ipsi") {
    p <- unname(stats::predict(
      treatment_model$model,
      newdata = fit_data,
      type = "response"
    ))
    delta <- intervention$delta
    return(ipsi_weight_formula(a_obs, p, delta))
  }

  # Horvitz-Thompson indicator branch. Covers point-mass
  # interventions on discrete treatments: static and deterministic
  # dynamic rules on Bernoulli / categorical. The intervened
  # distribution is a Dirac at `target_i`, so the correct weight is
  #   w_i = I(A_obs_i = target_i) / f(A_obs_i | L_i).
  # Rows with A_obs != target get zero weight and drop out of the
  # Hájek mean — which is exactly what the Horvitz-Thompson
  # estimator of E[Y^{d}] wants.
  is_ht <- iv_type %in%
    c("static", "dynamic") &&
    family_tag %in% c("bernoulli", "categorical")
  if (is_ht) {
    # Resolve target values per individual. Static is trivially
    # `rep(value, n)`; dynamic delegates to the user rule via
    # `apply_intervention_to_values()`. Both paths go through the
    # same post-processing below so the density guard and indicator
    # cast are unified.
    target <- apply_intervention_to_values(intervention, fit_data, a_obs)

    f_obs <- evaluate_density(treatment_model, a_obs, fit_data)
    check_density_positivity(f_obs, "density-ratio weights")

    # `==` on numeric 0/1 and character / factor treatments both
    # return logical — `as.numeric()` collapses to the HT indicator.
    ind <- as.numeric(a_obs == target)
    return(ind / f_obs)
  }

  # Smooth pushforward branch. Covers `shift()` and `scale_by()` on
  # continuous treatments. The weight is
  #   w_i = f_d(A_obs_i | L_i) / f(A_obs_i | L_i)
  # where f_d is the pushforward of f under the intervention map d.
  # For an invertible, smooth map,
  #   f_d(y | l) = f(d^{-1}(y, l) | l) * |det D d^{-1}(y, l)|,
  # so evaluating at y = A_obs gives
  #   w_i = f(d^{-1}(A_obs_i, L_i) | L_i) * |Jac| / f(A_obs_i | L_i).
  #
  # IMPORTANT: this is NOT the same as `f(d(A, L) | L) / f(A | L)`,
  # which is what the naive "evaluate the density at the intervened
  # treatment" formula gives. The difference is a sign (for shift)
  # and a Jacobian (for scale_by), and using the naive form gives
  # E[Y^{shift(-delta)}] instead of E[Y^{shift(delta)}]. Verified
  # analytically in the design notes.
  if (iv_type == "shift") {
    delta <- intervention$delta
    # shift: d(a) = a + delta, d^{-1}(y) = y - delta, |Jac| = 1.
    # f_d(A_obs | l) = f(A_obs - delta | l).
    a_eval <- a_obs - delta
    f_obs <- evaluate_density(treatment_model, a_obs, fit_data)
    f_int <- evaluate_density(treatment_model, a_eval, fit_data)
    check_density_positivity(f_obs, "shift density ratio")
    return(f_int / f_obs)
  }

  if (iv_type == "scale") {
    fct <- intervention$factor
    if (fct == 0) {
      rlang::abort(
        "`scale_by(0)` collapses the treatment support; not a valid MTP."
      )
    }
    # scale_by(c): d(a) = c * a, d^{-1}(y) = y / c, |Jac| = 1 / |c|.
    # f_d(A_obs | l) = f(A_obs / c | l) / |c|.
    a_eval <- a_obs / fct
    f_obs <- evaluate_density(treatment_model, a_obs, fit_data)
    f_int <- evaluate_density(treatment_model, a_eval, fit_data)
    check_density_positivity(f_obs, "scale_by density ratio")
    return((f_int / f_obs) / abs(fct))
  }

  # `threshold()` and continuous `dynamic()` would land here but
  # are blocked upstream by `check_intervention_family_compat()`.
  # If we reach this line something's wrong in the compat check.
  rlang::abort(
    paste0(
      "Internal error: `compute_density_ratio_weights()` has no branch ",
      "for (intervention = '",
      iv_type,
      "', family = '",
      family_tag,
      "'). This should have been caught by check_intervention_family_compat()."
    )
  )
}


#' Guard against zero or non-finite fitted density values
#'
#' @description
#' A zero (or non-finite) entry in `f_obs` means the observed
#' treatment value is outside the support of the fitted density.
#' That is a positivity violation, not a numerical hiccup: it would
#' silently produce `Inf` / `NaN` weights that the MSM fit and the
#' sandwich variance would both propagate. Abort with a clear,
#' context-tagged message so the user knows where to look.
#'
#' @param f Numeric density vector.
#' @param context Character scalar describing the caller (goes into
#'   the error message).
#' @return `NULL` invisibly; aborts on violation.
#' @noRd
check_density_positivity <- function(f, context) {
  bad <- which(!is.finite(f) | f <= 0)
  if (length(bad) == 0L) {
    return(invisible(NULL))
  }
  rlang::abort(
    paste0(
      length(bad),
      " observation(s) have zero or non-finite density under the fitted ",
      "treatment model (",
      context,
      "). This is a positivity violation — widen the treatment model ",
      "(add confounders, use a more flexible `propensity_model_fn`, or ",
      "truncate the treatment range) before refitting."
    )
  )
}


#' Construct the `alpha -> w(alpha)` closure used by the variance engine
#'
#' @description
#' Returns a closure that recomputes the density-ratio weight vector
#' under a candidate propensity parameter `alpha`. The variance
#' engine's Branch B body (`prepare_propensity_if_self_contained()`)
#' hands this closure to `numDeriv::jacobian()` to compute the
#' cross-derivative
#' \deqn{A_{\beta\alpha} = -\nabla_\alpha \bar\psi_\beta(\alpha)}
#' that the propensity-uncertainty correction needs.
#'
#' ## Design: what the closure captures
#'
#' The closure captures everything it needs by value at creation time:
#'
#' - `X_prop` — propensity design matrix (so recomputing predictions
#'   is one matrix multiply, no formula parsing).
#' - `a_obs` — observed treatment values on the fit rows.
#' - `a_int` — **precomputed** intervened treatment values for
#'   non-IPSI interventions. The intervention itself does not depend
#'   on `alpha`, only on the data, so we compute `a_int` once at
#'   closure-creation time rather than on every `numDeriv::jacobian()`
#'   step.
#' - `sigma` — residual SD (continuous treatments, fixed at fit time).
#' - `delta` — odds multiplier (IPSI only).
#' - `family` tag — dispatches to the right weight formula per step.
#'
#' Because the closure contains no references to mutable state, it is
#' safe to call many times with different `alpha` values inside a
#' `numDeriv::jacobian()` evaluation loop.
#'
#' ## IPSI vs data interventions
#'
#' IPSI's closed form has a different shape (no density evaluation at
#' an intervened treatment value — the intervention lives on the
#' propensity itself) so it gets its own branch. All other intervention
#' shapes share the `f_int / f_obs` branch.
#'
#' @inheritParams compute_density_ratio_weights
#'
#' @return A `function(alpha)` returning a numeric vector of length
#'   `sum(treatment_model$fit_rows)`.
#'
#' @noRd
make_weight_fn <- function(treatment_model, data, intervention) {
  if (!inherits(treatment_model, "causatr_treatment_model")) {
    rlang::abort(
      "`treatment_model` must be a `causatr_treatment_model`."
    )
  }

  fit_rows <- treatment_model$fit_rows
  fit_data <- data[fit_rows]
  trt_col <- treatment_model$treatment
  a_obs <- fit_data[[trt_col]]
  n_fit <- length(a_obs)

  # Natural course: constant 1 regardless of alpha. Still return a
  # closure so the variance engine's per-intervention loop has a
  # uniform contract. `numDeriv::jacobian()` on a constant function
  # returns zeros — correct, because the natural-course weight
  # carries no propensity-uncertainty.
  if (is.null(intervention)) {
    return(function(alpha) rep(1, n_fit))
  }

  check_intervention_family_compat(intervention, treatment_model)

  X_prop <- treatment_model$X_prop
  sigma <- treatment_model$sigma
  family_tag <- treatment_model$family
  iv_type <- intervention$type

  # ---- IPSI branch -------------------------------------------------
  # Closed form; no density evaluation at an intervened value.
  if (iv_type == "ipsi") {
    delta <- intervention$delta
    return(function(alpha) {
      # eta = X %*% alpha, p = plogis(eta). The matrix multiply is the
      # hot path inside the numDeriv::jacobian loop — ~p_alpha extra
      # perturbations per jacobian step, each costing one O(n * p_alpha)
      # multiply.
      eta <- as.numeric(X_prop %*% alpha)
      p <- stats::plogis(eta)
      ipsi_weight_formula(a_obs, p, delta)
    })
  }

  # ---- Horvitz-Thompson indicator branch ---------------------------
  # Point-mass interventions on discrete treatment: static and
  # deterministic dynamic on Bernoulli / categorical. Precompute the
  # indicator mask once at closure-creation time — it depends on
  # `data` and `intervention` only, so it is invariant under the
  # numDeriv perturbations of `alpha`.
  is_ht <- iv_type %in%
    c("static", "dynamic") &&
    family_tag %in% c("bernoulli", "categorical")
  if (is_ht) {
    target <- apply_intervention_to_values(intervention, fit_data, a_obs)
    ind <- as.numeric(a_obs == target)

    if (family_tag == "bernoulli") {
      return(function(alpha) {
        eta <- as.numeric(X_prop %*% alpha)
        p <- stats::plogis(eta)
        # f_obs for a Bernoulli(p) at observed 0/1 value A_obs.
        f_obs <- ifelse(a_obs == 1, p, 1 - p)
        ind / f_obs
      })
    }

    # Categorical is not supported: `evaluate_density()` has no
    # multinomial arm.
    rlang::abort(
      "Categorical HT weight closure is not supported."
    )
  }

  # ---- Smooth pushforward branch (continuous treatments) -----------
  # `shift()` and `scale_by()` on Gaussian-family treatment. Both use
  # the pushforward weight
  #   w_i = f(d^{-1}(A_obs_i, L_i) | L_i) * |Jac| / f(A_obs_i | L_i)
  # — see the `compute_density_ratio_weights()` body for the
  # derivation and the sign/Jacobian gotcha. The "evaluation point"
  # `a_eval = d^{-1}(A_obs)` and the Jacobian `|Jac|` are constants
  # of `data` and `intervention`, so we precompute them once at
  # closure-creation time.
  if (family_tag == "gaussian") {
    if (iv_type == "shift") {
      delta <- intervention$delta
      a_eval <- a_obs - delta
      jac_abs <- 1
    } else if (iv_type == "scale") {
      fct <- intervention$factor
      if (fct == 0) {
        rlang::abort(
          "`scale_by(0)` collapses the treatment support; not a valid MTP."
        )
      }
      a_eval <- a_obs / fct
      jac_abs <- abs(1 / fct)
    } else {
      rlang::abort(
        paste0(
          "Internal error: `make_weight_fn()` has no gaussian branch for '",
          iv_type,
          "'. check_intervention_family_compat() should have rejected this."
        )
      )
    }

    # `sigma` captured from the fit is held fixed under alpha
    # perturbation. The variance engine treats only the mean
    # parameters as the propensity nuisance — sigma as an additional
    # nuisance would require a joint M-estimation setup for (mu,
    # sigma), which is deferred.
    return(function(alpha) {
      mu <- as.numeric(X_prop %*% alpha)
      f_obs <- stats::dnorm(a_obs, mean = mu, sd = sigma)
      f_eval <- stats::dnorm(a_eval, mean = mu, sd = sigma)
      (f_eval / f_obs) * jac_abs
    })
  }

  rlang::abort(
    paste0(
      "`make_weight_fn()` does not handle (intervention = '",
      iv_type,
      "', family = '",
      family_tag,
      "')."
    )
  )
}


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
      # `rep(value, length)` preserves scalar type — `static("chemo")`
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
      # point at the same sharp edges — factor coercion, unknown
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
      # directly — reaching here is a bug in the caller.
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
#' counterfactual treatment value is needed — the intervention acts
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


#' Check that an intervention is compatible with the treatment family
#'
#' @description
#' Rejects nonsensical (intervention, treatment family) combinations
#' up front, with actionable error messages pointing at the right
#' alternative. The compatibility matrix the self-contained IPW
#' engine supports is:
#'
#' | intervention    | bernoulli | gaussian         | categorical |
#' |-----------------|-----------|------------------|-------------|
#' | `static()`      | ✓ HT      | ⛔                | ⛔           |
#' | `shift()`       | ⛔         | ✓ smooth ratio   | ⛔           |
#' | `scale_by()`    | ⛔         | ✓ smooth ratio   | ⛔           |
#' | `threshold()`   | ⛔         | ⛔ (use gcomp)    | ⛔           |
#' | `dynamic()`     | ✓ HT      | ⛔ (use gcomp)    | ⛔           |
#' | `ipsi()`        | ✓ Kennedy | ⛔                | ⛔           |
#'
#' The rationale for each ⛔:
#'
#' - `static(v)` on a continuous treatment: no one is observed
#'   exactly at `v`, so the Horvitz–Thompson indicator is 0 almost
#'   surely. The user almost certainly wanted a smooth shift.
#' - `shift()` / `scale_by()` / `threshold()` on a non-numeric
#'   treatment: factor / character treatments have no additive,
#'   multiplicative, or clamp structure.
#' - `threshold(lo, hi)` on continuous: the pushforward of a
#'   continuous density under a boundary clamp is a mixed measure
#'   (continuous density on `(lo, hi)` plus point masses at the
#'   boundaries), so there is no Lebesgue density and the density
#'   ratio w.r.t. the fitted `f(a|l)` is not well-defined. Users
#'   should use `estimator = "gcomp"` — the clamped-treatment
#'   counterfactual has a clean `predict(outcome_model, newdata =
#'   clamped)` interpretation there. This is a genuine gap vs
#'   `gcomp`, not a scope decision, and is documented in the
#'   FEATURE_COVERAGE_MATRIX as a rejection row.
#' - `dynamic()` on continuous: a deterministic per-individual
#'   target is a Dirac per individual; same pushforward problem as
#'   `threshold()`.
#' - `ipsi()` on anything but binary: the Kennedy (2019) closed-form
#'   formula is Bernoulli-specific.
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

  if (iv_type %in% c("shift", "scale", "threshold") && fam != "gaussian") {
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
          "'. Use `static()` for binary / categorical treatments."
        )
      )
    )
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
        i = "Use `estimator = 'gcomp'` — the clamped-treatment counterfactual is well-defined there via predict-then-average on the outcome model."
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
