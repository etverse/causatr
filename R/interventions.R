#' Set treatment to a fixed value
#'
#' @description
#' Creates a static intervention that sets the treatment to a fixed value for
#' all individuals at all time points. The most common intervention type,
#' corresponding to "always treat" (`static(1)`) or "never treat" (`static(0)`).
#'
#' @param value The fixed treatment value.
#'
#' @return A `causatr_intervention` object.
#'
#' @references
#' Hernan MA, Robins JM (2025). *Causal Inference: What If*. Chapman &
#' Hall/CRC. Chapter 1 (static interventions).
#'
#' @examples
#' \dontrun{
#' data("nhefs", package = "causatr")
#' fit <- causat(nhefs, outcome = "wt82_71", treatment = "qsmk",
#'               confounders = ~ sex + age + wt71)
#' contrast(fit, interventions = list(quit = static(1), continue = static(0)))
#' }
#'
#' @seealso [shift()], [dynamic()], [scale_by()], [threshold()], [ipsi()],
#'   [stochastic()], [contrast()]
#' @export
static <- function(value) {
  # Accept both numeric (binary/continuous treatments) and character
  # (categorical treatments, e.g. `static("chemo")`). `apply_single_intervention()`
  # handles both by blindly assigning via `:=`, which is type-preserving.
  ok <- (is.numeric(value) || is.character(value)) &&
    length(value) == 1L &&
    !is.na(value)
  if (!ok) {
    rlang::abort("`value` must be a single non-NA number or character string.")
  }
  # `new_causatr_intervention()` is a thin S3 constructor -- just tags
  # the list with the right class. Keeping constructors light means
  # they're cheap to create even in tight `lapply` loops over many
  # intervention variants.
  new_causatr_intervention("static", list(value = value))
}

#' Shift treatment by a fixed amount
#'
#' @description
#' Creates a modified treatment policy (MTP) that adds a fixed `delta` to each
#' individual's observed treatment value. Useful for continuous treatments where
#' a population-level shift is the relevant intervention (e.g., "reduce
#' exposure by 10 units").
#'
#' @param delta Numeric. The amount to add to the observed treatment.
#'
#' @return A `causatr_intervention` object.
#'
#' @references
#' Diaz I, Williams N, Hoffman KL, Schenck EJ (2023). Non-parametric causal
#' effects based on longitudinal modified treatment policies. *Journal of the
#' American Statistical Association* 118:846-857.
#'
#' @examples
#' \dontrun{
#' fit <- causat(nhefs, outcome = "wt82_71", treatment = "smokeintensity",
#'               confounders = ~ sex + age + wt71)
#' contrast(fit, interventions = list(
#'   reduce10 = shift(-10),
#'   observed = shift(0)
#' ))
#' }
#'
#' @seealso [static()], [scale_by()], [threshold()], [dynamic()], [ipsi()],
#'   [stochastic()]
#' @export
shift <- function(delta) {
  if (!(is.numeric(delta) && length(delta) == 1L && !is.na(delta))) {
    rlang::abort("`delta` must be a single non-NA number.")
  }
  new_causatr_intervention("shift", list(delta = delta))
}

#' Multiply treatment by a fixed factor
#'
#' @description
#' Creates a modified treatment policy that multiplies each individual's
#' observed treatment by `factor`. Useful for proportional reductions or
#' increases in continuous treatments.
#'
#' @param factor Numeric. The multiplicative factor.
#'
#' @return A `causatr_intervention` object.
#'
#' @examples
#' \dontrun{
#' fit <- causat(nhefs, outcome = "wt82_71", treatment = "smokeintensity",
#'               confounders = ~ sex + age + wt71)
#' contrast(fit, interventions = list(
#'   halved = scale_by(0.5),
#'   observed = scale_by(1)
#' ))
#' }
#'
#' @seealso [static()], [shift()], [threshold()], [dynamic()], [ipsi()],
#'   [stochastic()]
#' @export
scale_by <- function(factor) {
  if (!(is.numeric(factor) && length(factor) == 1L && !is.na(factor))) {
    rlang::abort("`factor` must be a single non-NA number.")
  }
  new_causatr_intervention("scale", list(factor = factor))
}

#' Clamp treatment within bounds
#'
#' @description
#' Creates a modified treatment policy that clamps each individual's observed
#' treatment to lie within `[lower, upper]`. Values below `lower` are set to
#' `lower`; values above `upper` are set to `upper`.
#'
#' @param lower Numeric. Lower bound (use `-Inf` for no lower bound).
#' @param upper Numeric. Upper bound (use `Inf` for no upper bound).
#'
#' @return A `causatr_intervention` object.
#'
#' @examples
#' \dontrun{
#' fit <- causat(nhefs, outcome = "wt82_71", treatment = "smokeintensity",
#'               confounders = ~ sex + age + wt71)
#' contrast(fit, interventions = list(
#'   capped20 = threshold(0, 20),
#'   observed  = shift(0)
#' ))
#' }
#'
#' @seealso [static()], [shift()], [scale_by()], [dynamic()], [ipsi()],
#'   [stochastic()]
#' @export
threshold <- function(lower = -Inf, upper = Inf) {
  if (!(is.numeric(lower) && length(lower) == 1L && !is.na(lower))) {
    rlang::abort("`lower` must be a single non-NA number.")
  }
  if (!(is.numeric(upper) && length(upper) == 1L && !is.na(upper))) {
    rlang::abort("`upper` must be a single non-NA number.")
  }
  if (lower > upper) {
    rlang::abort("`lower` must be <= `upper`.")
  }
  new_causatr_intervention("threshold", list(lower = lower, upper = upper))
}

#' Dynamic treatment rule
#'
#' @description
#' Creates a dynamic intervention where treatment is determined by a
#' user-supplied function of the covariate history. When `contrast()`
#' evaluates the intervention, `rule` receives two arguments:
#'
#' - `data`: the **full** counterfactual data.table (a copy of `fit$data`
#'   with all columns intact, including lag columns for longitudinal
#'   fits). For point treatments this is one row per individual; for
#'   longitudinal ICE it is the entire person-period panel at once. The
#'   rule is called once per intervention, not once per time step.
#' - `treatment`: the observed (or currently-held) treatment vector,
#'   length `nrow(data)`.
#'
#' `rule` must return a numeric vector of length `nrow(data)`. To branch
#' on time, reference the time column directly (e.g. `data$time == 0`).
#' To reference the treatment or covariate history in longitudinal data,
#' use the materialised lag columns (`data$lag1_A`, `data$lag1_L`, ...)
#' that `prepare_data()` built from `history`.
#'
#' @param rule A function with signature `function(data, treatment)` that
#'   returns a vector of treatment values of length `nrow(data)`.
#'
#' @return A `causatr_intervention` object.
#'
#' @references
#' Hernan MA, Robins JM (2025). *Causal Inference: What If*. Chapman &
#' Hall/CRC. Chapter 19 (dynamic treatment strategies).
#'
#' @examples
#' \dontrun{
#' # Point treatment: treat if CD4 count is below 200
#' cd4_rule <- dynamic(\(data, trt) ifelse(data$cd4 < 200, 1, 0))
#'
#' # Longitudinal: treat at time k if a time-varying confounder triggered
#' # at time k (the whole panel is passed in one call; the time column
#' # lets the rule dispatch per-row).
#' adaptive <- dynamic(\(data, trt) as.integer(!is.na(data$L) & data$L > 0))
#' }
#'
#' @seealso [static()], [shift()], [scale_by()], [threshold()], [ipsi()],
#'   [stochastic()]
#' @export
dynamic <- function(rule) {
  if (!is.function(rule)) {
    rlang::abort(
      "`rule` must be a function with signature function(data, treatment)."
    )
  }
  new_causatr_intervention("dynamic", list(rule = rule))
}

#' Incremental propensity score intervention
#'
#' @description
#' Creates an incremental propensity score intervention (IPSI) that multiplies
#' each individual's odds of treatment by `delta`. Values of `delta > 1`
#' increase the probability of treatment; `delta < 1` decrease it. This is a
#' stochastic modified treatment policy indexed by a single scalar, and
#' corresponds to the "shift-the-odds" counterfactual studied by
#' Kennedy (2019).
#'
#' IPSI is **binary-only** and is only meaningful under
#' `estimator = "ipw"` -- the IPW engine computes the Kennedy closed-form
#' weight
#' \eqn{w_i = (\delta A_i + (1 - A_i)) / (\delta p_i + (1 - p_i))}
#' directly from the fitted propensity, with no counterfactual
#' treatment value to predict at. `estimator = "gcomp"` and
#' `estimator = "matching"` therefore do not support `ipsi()`.
#'
#' @param delta Positive numeric. The odds multiplier.
#'
#' @return A `causatr_intervention` object.
#'
#' @references
#' Kennedy EH (2019). Nonparametric causal effects based on incremental
#' propensity score interventions. *Journal of the American Statistical
#' Association* 114:645-656.
#'
#' @examples
#' \dontrun{
#' fit <- causat(nhefs, outcome = "wt82_71", treatment = "qsmk",
#'               confounders = ~ sex + age + wt71)
#' contrast(fit, interventions = list(
#'   double_odds = ipsi(2),
#'   half_odds   = ipsi(0.5)
#' ))
#' }
#'
#' @seealso [static()], [shift()], [dynamic()], [scale_by()], [threshold()],
#'   [stochastic()]
#' @export
ipsi <- function(delta) {
  # IPSI transforms the treatment probability as
  #   p_new(L) = delta * p(L) / (1 - p(L) + delta * p(L))
  # which is equivalent to multiplying the ODDS of treatment by delta.
  # delta > 1 makes treatment more likely across the population;
  # delta < 1 makes it less likely. delta must be positive (the
  # odds multiplier can't be zero or negative -- you can't "anti-treat").
  #
  # The constructor just stores the scalar. Actual weight computation
  # happens in `compute_density_ratio_weights()` / `make_weight_fn()`
  # in `R/ipw_weights.R` via the Kennedy (2019) closed form
  #   w_i = (delta * A_i + (1 - A_i)) / (delta * p_i + (1 - p_i)),
  # which is only defined for binary treatments and is only usable
  # from the IPW engine (there is no counterfactual treatment value
  # for g-comp to predict at).
  #
  # Reject NA/NaN up front with a clean message -- otherwise the next
  # `delta <= 0` comparison evaluates to `NA` and users see the cryptic
  # "missing value where TRUE/FALSE needed". Mirrors the defensive
  # checks in `shift()` / `scale_by()` / `threshold()`. 2026-04-15
  # fourth-round critical review Issue #9.
  if (!rlang::is_scalar_double(delta) && !rlang::is_scalar_integer(delta)) {
    rlang::abort("`delta` must be a single positive number.")
  }
  if (is.na(delta)) {
    rlang::abort("`delta` must be a single non-NA positive number.")
  }
  if (delta <= 0) {
    rlang::abort("`delta` must be positive.")
  }
  new_causatr_intervention("ipsi", list(delta = delta))
}

#' Stochastic intervention (user-supplied sampling function)
#'
#' @description
#' Creates a stochastic intervention where the counterfactual treatment
#' for each individual is drawn from a user-supplied distribution that
#' may depend on the individual's covariates. Each call to `sampler`
#' should return an independent draw.
#'
#' Only supported under `estimator = "gcomp"` (point and longitudinal).
#' The g-formula evaluates \eqn{E[Y^g]} via Monte Carlo integration: for each
#' of `n_mc` draws, the sampler assigns counterfactual treatments, the
#' outcome model predicts, and the predictions are averaged across draws.
#'
#' @param sampler A function with signature `function(data, treatment)`
#'   that returns a vector of treatment values of length `nrow(data)`.
#'   `data` is the full counterfactual data.table (same as [dynamic()]);
#'   `treatment` is the observed treatment vector. Each call must return
#'   an independent random draw from the stochastic policy.
#' @param n_mc Positive integer. Number of Monte Carlo draws for the
#'   g-formula integration. Default 100. Larger values reduce MC noise
#'   at the cost of computation time. For sandwich variance, `n_mc`
#'   should be large enough that MC noise is negligible relative to
#'   the estimation error (100--500 is typical).
#'
#' @return A `causatr_intervention` object.
#'
#' @references
#' Diaz I, van der Laan MJ (2012). Population intervention causal
#' effects based on stochastic interventions. *Biometrics* 68:541--549.
#'
#' @examples
#' \dontrun{
#' # Binary treatment: covariate-dependent randomisation
#' stochastic(\(data, trt) rbinom(nrow(data), 1, plogis(0.5 + 0.3 * data$L)))
#'
#' # Continuous treatment: additive random shift
#' stochastic(\(data, trt) trt + rnorm(length(trt), mean = 0.5, sd = 0.25))
#' }
#'
#' @seealso [static()], [shift()], [scale_by()], [threshold()], [dynamic()],
#'   [ipsi()], [contrast()]
#' @export
stochastic <- function(sampler, n_mc = 100L) {
  if (!is.function(sampler)) {
    rlang::abort(
      "`sampler` must be a function with signature function(data, treatment)."
    )
  }
  if (
    !(rlang::is_scalar_integerish(n_mc) &&
      !is.na(n_mc) &&
      n_mc >= 1L)
  ) {
    rlang::abort("`n_mc` must be a positive integer.")
  }
  n_mc <- as.integer(n_mc)
  if (n_mc < 10L) {
    rlang::warn(
      paste0(
        "`n_mc = ",
        n_mc,
        "` is very low for Monte Carlo integration. ",
        "Consider using n_mc >= 100 for reliable estimates."
      )
    )
  }
  new_causatr_intervention("stochastic", list(sampler = sampler, n_mc = n_mc))
}

#' Print a causatr intervention
#'
#' @description
#' Displays the intervention type and its parameters.
#'
#' @param x A `causatr_intervention` object.
#' @param ... Currently unused.
#' @return Invisibly returns `x`.
#' @export
print.causatr_intervention <- function(x, ...) {
  cat("<causatr_intervention: ", x$type, ">\n", sep = "")
  params <- x[names(x) != "type"]
  for (nm in names(params)) {
    if (is.function(params[[nm]])) {
      cat(" ", nm, ": <function>\n", sep = "")
    } else {
      cat(" ", nm, ": ", params[[nm]], "\n", sep = "")
    }
  }
  invisible(x)
}

#' Apply a causal intervention to a copy of the data
#'
#' @description
#' Copies `data` and overwrites the treatment column(s) with the values
#' implied by the intervention.  For a scalar treatment `iv` is a single
#' `causatr_intervention`; for a multivariate treatment `iv` is a named list
#' with one element per treatment variable.
#'
#' Called once per intervention inside `compute_contrast()` and
#' `variance_bootstrap()` to create counterfactual datasets.
#'
#' @param data A data.table containing all model variables.
#' @param treatment Character scalar or vector. Treatment column name(s).
#' @param iv A `causatr_intervention`, a named list of them, or `NULL`
#'   (natural course -- data returned unchanged).
#'
#' @return A modified *copy* of `data` (original is never mutated).
#'
#' @noRd
apply_intervention <- function(data, treatment, iv) {
  # `data.table::copy()` forces a deep copy so in-place `:=` mutations
  # on the treatment column inside `apply_single_intervention()` don't
  # leak back into the caller's `data`. This is essential because
  # compute_contrast() calls us once per intervention and expects
  # each counterfactual frame to be independent.
  data_a <- data.table::copy(data)

  # `iv = NULL` means "natural course": return the data as-is with
  # observed treatment values untouched. This is the correct
  # reference for MTPs (Diaz et al. 2023) and for longitudinal
  # dynamic regimes (Hernan & Robins Ch. 21).
  if (is.null(iv)) {
    return(data_a)
  }

  # Dispatch on treatment arity. The two branches exist because
  # `iv` has a different shape in each case:
  #   - scalar treatment -> iv is a single causatr_intervention
  #   - multivariate     -> iv is a named list, one per treatment column
  # The check_intervention_list() upstream validator ensures the
  # shape matches the treatment vector length, so by the time we
  # reach here we can trust the structure.
  if (length(treatment) == 1L) {
    apply_single_intervention(data_a, treatment, iv)
  } else {
    # Multivariate treatment: verify `names(iv)` matches the treatment
    # vector *exactly*. Upstream `check_intervention_list()` validated
    # that every sub-element is a `causatr_intervention`, but did not
    # cross-check the sub-list names against the `treatment` argument
    # (that validator has no access to the fit). Without this guard a
    # user typo like `list(X = static(1), A2 = shift(-10))` for
    # `treatment = c("A1", "A2")` silently creates a spurious `X`
    # column via data.table's `[, (X) := ...]` and never touches `A1`.
    missing <- setdiff(treatment, names(iv))
    extra <- setdiff(names(iv), treatment)
    if (length(missing) > 0L || length(extra) > 0L) {
      rlang::abort(
        paste0(
          "Multivariate intervention names must match the `treatment` ",
          "argument exactly. ",
          if (length(missing) > 0L) {
            paste0(
              "Missing: ",
              paste(shQuote(missing), collapse = ", "),
              ". "
            )
          },
          if (length(extra) > 0L) {
            paste0(
              "Unexpected: ",
              paste(shQuote(extra), collapse = ", "),
              ". "
            )
          },
          "Expected names: ",
          paste(shQuote(treatment), collapse = ", "),
          "."
        ),
        .call = FALSE
      )
    }
    # Apply each named sub-intervention to its corresponding column.
    # Order matters only if interventions reference each other's
    # treatment columns (they generally don't), so we just iterate.
    for (trt_nm in names(iv)) {
      apply_single_intervention(data_a, trt_nm, iv[[trt_nm]])
    }
    data_a
  }
}

#' Validate that an intervention return value matches the treatment column type
#'
#' @description
#' Shared validation for `dynamic()` and `stochastic()` branches of
#' `apply_single_intervention()`. Checks that `new_trt` is type-compatible
#' with the original treatment column (numeric, factor, or character) and
#' has the correct length. For factor columns, accepts either a matching
#' factor or a character vector whose values are existing levels (coerced
#' back to factor). Aborts with an informative message on mismatch.
#'
#' @param new_trt The vector returned by the user's rule/sampler.
#' @param orig_trt The original treatment vector from the data.
#' @param trt_col Character. Name of the treatment column (for messages).
#' @param iv_label Character. Intervention label for messages
#'   (e.g. `"dynamic()"` or `"stochastic()"`).
#' @param n_rows Integer. Expected length of `new_trt`.
#'
#' @return The validated (and possibly coerced) `new_trt` vector.
#'
#' @noRd
validate_intervention_return <- function(
  new_trt,
  orig_trt,
  trt_col,
  iv_label,
  n_rows
) {
  if (is.numeric(orig_trt)) {
    if (!is.numeric(new_trt)) {
      rlang::abort(
        paste0(
          "`",
          iv_label,
          "` rule returned ",
          typeof(new_trt),
          " but treatment column `",
          trt_col,
          "` is numeric."
        ),
        .call = FALSE
      )
    }
  } else if (is.factor(orig_trt)) {
    if (is.character(new_trt)) {
      unknown <- setdiff(unique(new_trt), levels(orig_trt))
      if (length(unknown) > 0L) {
        rlang::abort(
          paste0(
            "`",
            iv_label,
            "` rule returned value(s) not present as ",
            "levels of factor treatment `",
            trt_col,
            "`: ",
            paste(shQuote(unknown), collapse = ", "),
            "."
          ),
          .call = FALSE
        )
      }
      new_trt <- factor(new_trt, levels = levels(orig_trt))
    } else if (is.factor(new_trt)) {
      if (!identical(levels(new_trt), levels(orig_trt))) {
        rlang::abort(
          paste0(
            "`",
            iv_label,
            "` rule returned a factor with levels that ",
            "do not match treatment column `",
            trt_col,
            "`."
          ),
          .call = FALSE
        )
      }
    } else {
      rlang::abort(
        paste0(
          "`",
          iv_label,
          "` rule returned ",
          typeof(new_trt),
          " but treatment column `",
          trt_col,
          "` is a factor -- return a character or factor vector."
        ),
        .call = FALSE
      )
    }
  } else if (is.character(orig_trt)) {
    if (!is.character(new_trt)) {
      rlang::abort(
        paste0(
          "`",
          iv_label,
          "` rule returned ",
          typeof(new_trt),
          " but treatment column `",
          trt_col,
          "` is character."
        ),
        .call = FALSE
      )
    }
  } else {
    rlang::abort(
      paste0(
        "`",
        iv_label,
        "` does not support treatment columns of type ",
        typeof(orig_trt),
        "."
      ),
      .call = FALSE
    )
  }
  if (length(new_trt) != n_rows) {
    rlang::abort(
      paste0(
        "`",
        iv_label,
        "` rule must return a vector of length ",
        n_rows,
        " (got ",
        length(new_trt),
        ")."
      ),
      .call = FALSE
    )
  }
  new_trt
}


#' Apply one intervention rule to a single treatment column (in-place)
#'
#' @description
#' Modifies `data` in-place (via data.table `:=`) by transforming `trt_col`
#' according to the intervention type.
#'
#' @param data A data.table (modified by reference).
#' @param trt_col Character. Name of the treatment column to modify.
#' @param iv A `causatr_intervention` object.
#'
#' @return `data` invisibly (the mutation happens in place).
#'
#' @noRd
apply_single_intervention <- function(data, trt_col, iv) {
  switch(
    iv$type,
    static = {
      # Set every individual's treatment to the fixed value.
      data[, (trt_col) := iv$value]
    },
    shift = {
      # Shift each individual's observed treatment by a constant delta.
      data[, (trt_col) := get(trt_col) + iv$delta]
    },
    scale = {
      # Scale each individual's observed treatment by a constant factor.
      data[, (trt_col) := get(trt_col) * iv$factor]
    },
    threshold = {
      # Clamp each individual's treatment to [lower, upper].
      data[, (trt_col) := pmax(pmin(get(trt_col), iv$upper), iv$lower)]
    },
    dynamic = {
      orig_trt <- data[[trt_col]]
      new_trt <- iv$rule(data, orig_trt)
      new_trt <- validate_intervention_return(
        new_trt,
        orig_trt,
        trt_col,
        "dynamic()",
        nrow(data)
      )
      data[, (trt_col) := new_trt]
    },
    stochastic = {
      orig_trt <- data[[trt_col]]
      new_trt <- iv$sampler(data, orig_trt)
      new_trt <- validate_intervention_return(
        new_trt,
        orig_trt,
        trt_col,
        "stochastic()",
        nrow(data)
      )
      data[, (trt_col) := new_trt]
    },
    ipsi = {
      # IPSI has no counterfactual treatment value to assign to the
      # column -- the intervention acts on the propensity, not on the
      # treatment itself. The IPW engine handles IPSI by building a
      # closed-form weight vector in `compute_density_ratio_weights()`
      # and never calls `apply_single_intervention()` on the
      # treatment column at all. If we reach this branch, the caller
      # is one of the non-IPW estimators (gcomp, matching), which do
      # not support IPSI.
      rlang::abort(
        c(
          "`ipsi()` interventions are only supported under `estimator = 'ipw'`.",
          i = "The intervention shifts the propensity, not the treatment value, so there is no counterfactual treatment to predict at under g-computation or matching.",
          i = "Use `causat(..., estimator = 'ipw')` with an IPSI intervention, or rewrite the intervention as a `shift()` / `scale_by()` / `static()` for g-comp / matching."
        ),
        .call = FALSE
      )
    },
    rlang::abort(
      paste0("Unknown intervention type: '", iv$type, "'."),
      .call = FALSE
    )
  )
  data
}
