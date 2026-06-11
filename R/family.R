#' Resolve a family argument to a family object
#'
#' Accepts a character string (e.g. `"gaussian"`), a family function
#' (e.g. `stats::binomial`), or an already-evaluated family object
#' (e.g. `stats::binomial()`).
#'
#' @param family Character, function, or family object.
#' @return A family object (list with `$family`, `$link`, etc.).
#' @noRd
resolve_family <- function(family) {
  # Three-way dispatch to canonicalize every input to a family object.
  # Downstream code universally relies on `family$family` (character)
  # and `family$link` (character) being present -- normalizing here
  # means each call site doesn't need its own isinstance check.
  if (is.character(family)) {
    # Look up in stats:: namespace explicitly -- avoids picking up a
    # user-defined function with the same name in the caller's env.
    # Wrap in tryCatch because some stats:: functions share names with
    # family constructors (e.g. `stats::beta` is the beta distribution
    # function, not a GLM family; calling it returns a numeric value,
    # not a family object, and the `inherits(., "family")` guard below
    # catches that).
    fam_obj <- tryCatch(
      {
        fam_fn <- get(family, mode = "function", envir = asNamespace("stats"))
        fam_fn()
      },
      error = function(e) NULL
    )
    if (!is.null(fam_obj) && inherits(fam_obj, "family")) {
      return(fam_obj)
    }

    # Extended families from optional packages.
    if (family == "beta") {
      rlang::check_installed(
        "betareg",
        reason = "for beta regression outcomes (family = \"beta\")"
      )
      # `betareg::betareg()` doesn't accept a `family` argument; the beta
      # family object built here is used only by `is_binary_family()` (returns
      # FALSE so the non-binary estimation path is taken) and by variance
      # tier selection. The link slots are populated so downstream code that
      # calls `family$linkinv()` on predictions doesn't error.
      logit <- stats::make.link("logit")
      return(structure(
        list(
          family = "beta",
          link = "logit",
          linkfun = logit$linkfun,
          linkinv = logit$linkinv,
          mu.eta = logit$mu.eta,
          variance = function(mu) mu * (1 - mu)
        ),
        class = "family"
      ))
    }

    rlang::abort(paste0(
      "Unknown family: '",
      family,
      "'. ",
      "Supported: gaussian, binomial, poisson, quasibinomial, ",
      "quasipoisson, Gamma, inverse.gaussian, beta. ",
      "Or pass a family object directly."
    ))
  }
  if (is.function(family)) {
    return(family())
  }
  # Already a family object (list with $family, $link, etc.).
  family
}

#' Does a fitting function accept a `family` argument?
#'
#' Inspects the formal arguments of `fn` to determine whether it accepts
#' a named `family` parameter. Used to conditionally strip `family` from
#' the argument list when calling model functions that handle their own
#' family internally (e.g. `MASS::glm.nb`, `betareg::betareg`).
#'
#' @param fn A function (the user's `model_fn`).
#' @return Logical scalar.
#' @noRd
fn_accepts_family <- function(fn) {
  # `formals()` returns the declared parameter list without executing the
  # function; it's the correct way to inspect the signature of an arbitrary
  # function object. `MASS::glm.nb` and `betareg::betareg` do not declare
  # `family` in their formals, so this returns FALSE for both, causing the
  # fitter to omit the `family` argument rather than passing one those
  # functions would reject with a cryptic error.
  "family" %in% names(formals(fn))
}

#' Check whether a fitted outcome model is multinomial
#'
#' @description
#' A multinomial outcome model (currently `nnet::multinom`) predicts an
#' n-by-K matrix of class probabilities rather than a length-n vector, so
#' the g-computation estimand is the K-vector \eqn{P(Y = k \mid do(A = a))}
#' per intervention instead of a scalar mean. `contrast()` and the bootstrap
#' branch on this to assemble the per-class result.
#'
#' @param model A fitted outcome model object, or `NULL`.
#' @return Logical scalar. `TRUE` for `nnet::multinom` fits.
#' @noRd
is_multinom_outcome <- function(model) {
  inherits(model, "multinom")
}

#' Class labels of a multinomial outcome model
#'
#' @param model A fitted `nnet::multinom` object.
#' @return Character vector of the K outcome class labels, in the column
#'   order that `predict(model, type = "probs")` returns.
#' @noRd
multinom_class_labels <- function(model) {
  # nnet::multinom stores the response factor levels in `$lev`, in the same
  # order as the columns of `predict(type = "probs")`.
  model$lev
}

#' Check whether a family describes a binary outcome
#'
#' @param family A character string, family object, or function.
#' @return Logical scalar.
#' @noRd
is_binary_family <- function(family) {
  if (is.null(family)) {
    return(FALSE)
  }
  fam <- tryCatch(resolve_family(family), error = function(e) NULL)
  if (!is.null(fam) && is.character(fam$family)) {
    return(fam$family %in% c("binomial", "quasibinomial"))
  }
  FALSE
}

#' Swap binomial for quasibinomial in weighted MSM fits
#'
#' @description
#' IPW outcome MSMs are fit with density-ratio weights, which are
#' non-integer. R's `glm.fit` warns "non-integer #successes in a
#' binomial glm!" for binomial families with non-integer weights.
#' quasibinomial produces identical coefficients, SEs, and
#' predictions but suppresses the warning.
#'
#' @param fam A family object (as returned by `resolve_family()`).
#' @return The same family, or `quasibinomial()` if the input was
#'   `binomial`.
#' @noRd
msm_family <- function(fam) {
  # IPW MSMs are fit with density-ratio weights, which are non-integer
  # fractions. `stats::glm` emits "non-integer #successes in a binomial
  # glm!" for every row when `family = binomial`. `quasibinomial` produces
  # identical coefficients, SEs, and predictions because it uses the same
  # IRLS update but skips the Pearson chi-squared check that triggers the
  # warning. The link is forwarded to preserve user-specified link choices
  # (e.g. `binomial("log")` -> `quasibinomial("log")`).
  if (is.list(fam) && identical(fam$family, "binomial")) {
    return(stats::quasibinomial(link = fam$link))
  }
  fam
}
