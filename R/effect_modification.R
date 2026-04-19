#' Parse effect-modification terms from a confounders formula
#'
#' Walks the term labels of `confounders` and classifies each as either a
#' pure confounder term (no treatment variable involved) or an effect-
#' modification term (interacts a treatment variable with one or more
#' baseline modifiers). Returns a structured object that downstream
#' callers (`fit_ipw`, `fit_matching`, `ice_build_formula`) use to
#' build method-appropriate MSM formulas.
#'
#' @param confounders One-sided formula, e.g. `~ L + sex + A:sex`.
#' @param treatment Character vector of treatment column name(s).
#' @return A list with class `"causatr_em_info"` containing:
#'   \describe{
#'     \item{em_terms}{List of parsed EM terms. Each element is a list
#'       with `$term` (original label, e.g. `"A:sex"`), `$treatment_var`
#'       (the treatment variable involved), and `$modifier_vars` (character
#'       vector of the non-treatment variables in the interaction).}
#'     \item{confounder_terms}{Character vector of non-EM term labels
#'       (passed through unchanged).}
#'     \item{modifier_vars}{Character vector of unique modifier variable
#'       names across all EM terms.}
#'     \item{has_em}{Logical scalar: `TRUE` when at least one EM term
#'       was found.}
#'   }
#' @noRd
parse_effect_mod <- function(confounders, treatment) {
  term_labels <- attr(stats::terms(confounders), "term.labels")

  em_terms <- list()
  confounder_terms <- character(0L)

  for (tl in term_labels) {
    # Extract all variable names referenced by the term. `all.vars()`
    # on the parsed expression handles `A:sex`, `sex:A`, `A:I(age > 65)`,
    # multi-way `A:sex:race`, etc.
    vars <- all.vars(parse(text = tl)[[1L]])
    trt_hit <- intersect(vars, treatment)

    if (length(trt_hit) > 0L) {
      # At least one treatment variable appears in this term.
      modifier_vars <- setdiff(vars, treatment)
      em_terms <- c(
        em_terms,
        list(list(
          term = tl,
          treatment_var = trt_hit,
          modifier_vars = modifier_vars
        ))
      )
    } else {
      confounder_terms <- c(confounder_terms, tl)
    }
  }

  all_modifiers <- unique(unlist(
    lapply(em_terms, function(x) x$modifier_vars),
    use.names = FALSE
  ))
  if (is.null(all_modifiers)) {
    all_modifiers <- character(0L)
  }

  structure(
    list(
      em_terms = em_terms,
      confounder_terms = confounder_terms,
      modifier_vars = all_modifiers,
      has_em = length(em_terms) > 0L
    ),
    class = "causatr_em_info"
  )
}


#' Build an MSM formula for IPW with effect-modification terms
#'
#' When effect-modification terms are present, the per-intervention
#' Hajek MSM expands from `Y ~ 1` to `Y ~ 1 + modifier_main_effects`.
#' The modifier main effects let `predict()` return stratum-specific
#' counterfactual means. No treatment term enters the MSM because the
#' density-ratio weights absorb it (the Hajek intercept under HT
#' weights already conditions on the intervention).
#'
#' When no EM terms are present, returns the standard `Y ~ 1`.
#'
#' @param outcome Character scalar. Outcome column name.
#' @param em_info A `causatr_em_info` object from `parse_effect_mod()`.
#' @return A formula suitable for the weighted MSM.
#' @noRd
build_ipw_msm_formula <- function(outcome, em_info) {
  if (!em_info$has_em) {
    return(stats::as.formula(paste0(outcome, " ~ 1")))
  }
  # Include modifier main effects only. The treatment is absorbed by
  # the density-ratio weights, so `A:modifier` becomes just `modifier`
  # in the MSM. For multiple modifiers (e.g. `A:sex + A:race`), all
  # modifier main effects enter: `Y ~ 1 + sex + race`.
  stats::reformulate(c("1", em_info$modifier_vars), response = outcome)
}


#' Build an MSM formula for matching with effect-modification terms
#'
#' When effect-modification terms are present, the matched-sample
#' outcome MSM expands from `Y ~ A` to `Y ~ A + modifier + A:modifier`.
#' This is the standard saturated MSM for effect modification under
#' matching: the matched design handles confounding, and the expanded
#' formula recovers stratum-specific treatment effects.
#'
#' When no EM terms are present, returns `Y ~ A`.
#'
#' @param outcome Character scalar. Outcome column name.
#' @param treatment Character scalar. Treatment column name.
#' @param em_info A `causatr_em_info` object from `parse_effect_mod()`.
#' @return A formula suitable for the matched-sample outcome model.
#' @noRd
build_matching_msm_formula <- function(outcome, treatment, em_info) {
  if (!em_info$has_em) {
    return(stats::reformulate(treatment, response = outcome))
  }
  # Reconstruct the EM interaction terms with the treatment. The user
  # wrote `A:sex` in confounders; we emit `A + sex + A:sex` on the
  # MSM RHS. For multiple modifiers: `A + sex + race + A:sex + A:race`.
  em_interaction_terms <- vapply(
    em_info$em_terms,
    function(x) x$term,
    character(1L)
  )
  rhs <- c(treatment, em_info$modifier_vars, em_interaction_terms)
  stats::reformulate(rhs, response = outcome)
}


#' Validate effect-modification terms for a given estimator
#'
#' Checks whether the EM terms detected by `parse_effect_mod()` are
#' supported under the chosen estimator and treatment configuration.
#' Currently only rejects one pattern:
#'
#' - **Bare treatment in confounders** (e.g. `~ L + A`): always rejected
#'   because `A` has no place in a propensity model of `A`.
#'
#' True EM interactions (`A:modifier`) are accepted for all estimators
#' and treatment types. Non-binary treatment + EM under IPW works via
#' the density-ratio engine (Phase 4); matching gates non-binary
#' treatment upstream in `fit_matching()`.
#'
#' Terms that are pure confounder interactions (e.g. `L1:L2`) pass
#' through unchecked.
#'
#' @param em_info A `causatr_em_info` object from `parse_effect_mod()`.
#' @param treatment Character vector of treatment column name(s).
#' @param estimator Character scalar: `"ipw"` or `"matching"`.
#' @param data A data.frame/data.table for inspecting treatment type.
#' @return `invisible(NULL)` on success; aborts otherwise.
#' @noRd
check_em_compat <- function(em_info, treatment, estimator, data = NULL) {
  # Bare treatment terms (modifier_vars is empty) are never valid: they
  # put A on both sides of the propensity model. This replaces the old
  # "bare A in confounders" branch of check_confounders_no_treatment().
  bare <- Filter(
    function(x) length(x$modifier_vars) == 0L,
    em_info$em_terms
  )
  if (length(bare) > 0L) {
    bare_labels <- vapply(bare, function(x) x$term, character(1L))
    rlang::abort(
      c(
        paste0(
          "`confounders` includes the treatment variable itself, which ",
          "creates a circular propensity model (`A ~ ... + A`)."
        ),
        x = paste0(
          "Offending term(s): ",
          paste(bare_labels, collapse = ", "),
          "."
        ),
        i = "Remove the treatment from `confounders`."
      ),
      class = "causatr_bare_treatment_in_confounders",
      .call = FALSE
    )
  }

  invisible(NULL)
}


#' Expand an EM term across treatment lags for ICE formulas
#'
#' Given an EM term like `"A:sex"` and `available_lags = 2`, produces
#' `c("lag1_A:sex", "lag2_A:sex")`. The original term (`"A:sex"`) is NOT
#' included -- callers add it separately via `baseline_terms`. The
#' expansion replaces only the treatment variable component of the `:`-
#' separated interaction; other components (modifiers, `I()` wrappers)
#' are kept verbatim.
#'
#' @param em_term A single EM term list from `parse_effect_mod()$em_terms`,
#'   with `$term`, `$treatment_var`, and `$modifier_vars`.
#' @param available_lags Integer. Number of lags available at the current
#'   ICE time step (`min(time_idx, max_lag)`).
#' @return Character vector of lag-expanded interaction terms, length
#'   `available_lags` (empty if `available_lags == 0`).
#' @noRd
expand_em_lag_terms <- function(em_term, available_lags) {
  if (available_lags == 0L) {
    return(character(0L))
  }
  # Split the interaction term into its `:`-separated components, find
  # which component is the treatment variable, and substitute the lag
  # prefix for each lag order. E.g. "A:sex" -> components = c("A","sex"),
  # treatment component is "A", lag1 -> c("lag1_A","sex") -> "lag1_A:sex".
  components <- strsplit(em_term$term, ":", fixed = TRUE)[[1L]]
  trt_var <- em_term$treatment_var
  # `%in%` instead of `==` so multivariate treatment_var (length > 1)
  # finds all treatment components without R's vector recycling silently
  # matching the wrong positions.
  trt_idx <- which(components %in% trt_var)
  if (length(trt_idx) == 0L) {
    rlang::warn(paste0(
      "expand_em_lag_terms(): treatment variable(s) ",
      paste(trt_var, collapse = ", "),
      " not found in EM term '",
      em_term$term,
      "'. Returning no lag terms for this term."
    ))
    return(character(0L))
  }

  vapply(
    seq_len(available_lags),
    function(k) {
      lag_components <- components
      lag_components[trt_idx] <- paste0("lag", k, "_", trt_var)
      paste(lag_components, collapse = ":")
    },
    character(1L)
  )
}
