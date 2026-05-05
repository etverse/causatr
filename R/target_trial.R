#' Specify a target trial protocol
#'
#' @description
#' Creates a structured description of the target trial that your causal
#' analysis emulates.
#' This is a documentation aid --- the object is not passed to [causat()] or
#' [contrast()].
#' It helps you reason through the components of the target trial
#' (Hernán & Robins 2025, Ch. 22) before writing the analysis code.
#'
#' @param eligibility Character. Who is eligible for the trial at time zero?
#' @param treatment_strategy Character. What treatment strategies are
#'   compared? (e.g. "initiate treatment A vs. remain untreated")
#' @param assignment Character. How are individuals assigned to strategies?
#'   (e.g. "random assignment at baseline")
#' @param follow_up Character. When does follow-up start and end?
#' @param outcome Character. What outcome is measured, and when?
#' @param causal_contrast Character. What causal contrast is estimated?
#'   (e.g. "intention-to-treat effect", "per-protocol effect")
#' @param model Character. What statistical model links the data to the
#'   causal parameter? (e.g. "parametric g-formula with logistic outcome
#'   model")
#' @param censoring Character. How is loss to follow-up or non-adherence
#'   handled? (e.g. "IPCW for dropout", "grace period of 2 visits")
#'
#' @return A `causatr_target_trial` object (S3 list with a print method).
#'
#' @examples
#' \dontrun{
#' target_trial(
#'   eligibility = "Adults aged 25-74 who smoke >= 5 cig/day",
#'   treatment_strategy = "Quit smoking vs. continue smoking",
#'   follow_up = "Baseline (1971) to end of follow-up (1982)",
#'   outcome = "Weight change in kg at end of follow-up",
#'   causal_contrast = "ATE: E[Y(quit)] - E[Y(continue)]"
#' )
#' }
#'
#' @seealso [causat()], [contrast()]
#' @export
target_trial <- function(
  eligibility = NULL,
  treatment_strategy = NULL,
  assignment = NULL,
  follow_up = NULL,
  outcome = NULL,
  causal_contrast = NULL,
  model = NULL,
  censoring = NULL
) {
  fields <- list(
    eligibility = eligibility,
    treatment_strategy = treatment_strategy,
    assignment = assignment,
    follow_up = follow_up,
    outcome = outcome,
    causal_contrast = causal_contrast,
    model = model,
    censoring = censoring
  )

  for (nm in names(fields)) {
    val <- fields[[nm]]
    if (!is.null(val) && !rlang::is_string(val)) {
      rlang::abort(
        paste0("`", nm, "` must be a single character string or NULL."),
        class = "causatr_bad_target_trial"
      )
    }
  }

  # Drop NULL entries
  fields <- Filter(Negate(is.null), fields)

  structure(fields, class = "causatr_target_trial")
}


#' @export
print.causatr_target_trial <- function(x, ...) {
  cat("<causatr_target_trial>\n")

  labels <- c(
    eligibility = "Eligibility",
    treatment_strategy = "Treatment strategy",
    assignment = "Assignment",
    follow_up = "Follow-up",
    outcome = "Outcome",
    causal_contrast = "Causal contrast",
    model = "Statistical model",
    censoring = "Censoring"
  )

  for (nm in names(labels)) {
    if (!is.null(x[[nm]])) {
      cat(" ", format(paste0(labels[[nm]], ":"), width = 22), x[[nm]], "\n")
    }
  }

  invisible(x)
}
