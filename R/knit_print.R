#' Knitr print method for `causatr_result` objects
#'
#' @description
#' When a [`causat()`] contrast is evaluated inside a knitr/Quarto document
#' and tinytable is installed, the object is rendered as a pair of HTML
#' tables (intervention means + contrasts) with a one-line metadata header
#' instead of the ASCII console output produced by [print.causatr_result()].
#' This lets vignettes and user reports fold the code chunk and still show
#' well-formatted tables.
#'
#' The method is registered as an S3 method for [knitr::knit_print()] in
#' [.onLoad()], so it only takes effect when knitr has been loaded (i.e.
#' during document rendering) and silently degrades to the default print
#' method when either knitr or tinytable is unavailable.
#'
#' @param x A `causatr_result` object (returned by [contrast()]).
#' @param ... Ignored. Present for S3 method compatibility.
#'
#' @return An object produced by [knitr::asis_output()] that, when rendered,
#'   emits a short metadata header followed by two tinytable HTML tables.
#'
#' @keywords internal
#' @noRd
knit_print.causatr_result <- function(x, ...) {
  # Guard against missing optional dependencies. tinytable is a Suggests
  # package used for vignettes; if the user rendering a doc doesn't have
  # it installed we fall back to the plain-text print method so the render
  # still succeeds.
  if (
    !requireNamespace("tinytable", quietly = TRUE) ||
      !requireNamespace("knitr", quietly = TRUE)
  ) {
    return(knitr::normal_print(x))
  }

  est_df <- as.data.frame(x$estimates)
  con_df <- as.data.frame(x$contrasts)
  num_est <- which(vapply(est_df, is.numeric, logical(1)))
  num_con <- which(vapply(con_df, is.numeric, logical(1)))

  # A compact metadata header: this replaces the multi-line ASCII banner
  # printed by print.causatr_result. Using non-breaking spaces so the
  # separators stick to their neighbours on narrow viewports.
  header <- sprintf(
    paste0(
      "**Estimator:** %s &nbsp;&middot;&nbsp; ",
      "**Estimand:** %s &nbsp;&middot;&nbsp; ",
      "**Contrast:** %s &nbsp;&middot;&nbsp; ",
      "**CI method:** %s &nbsp;&middot;&nbsp; ",
      "**N:** %d"
    ),
    x$estimator,
    x$estimand,
    x$type,
    x$ci_method,
    x$n
  )

  est_caption <- if (x$estimator == "snm") {
    "Blip parameters"
  } else if (!is.null(x$class_labels)) {
    "Class probabilities by intervention"
  } else {
    "Intervention means"
  }
  t_est <- tinytable::tt(est_df, caption = est_caption)
  if (length(num_est) > 0L) {
    t_est <- tinytable::format_tt(
      t_est,
      j = num_est,
      digits = 3,
      num_fmt = "decimal"
    )
  }

  # SNM Path A returns an empty contrasts table (the blip parameters
  # are the estimand directly). Skip the contrasts table when empty.
  if (nrow(con_df) > 0L) {
    t_con <- tinytable::tt(con_df, caption = "Contrasts")
    if (length(num_con) > 0L) {
      t_con <- tinytable::format_tt(
        t_con,
        j = num_con,
        digits = 3,
        num_fmt = "decimal"
      )
    }
    con_block <- paste0("\n\n", knitr::knit_print(t_con))
  } else {
    con_block <- ""
  }

  knitr::asis_output(paste0(
    header,
    "\n\n",
    knitr::knit_print(t_est),
    con_block
  ))
}
