#' @noRd
prepare_data <- function(
  data,
  outcome,
  treatment,
  confounders,
  id = NULL,
  time = NULL,
  censoring = NULL,
  call = rlang::caller_env()
) {
  if (!data.table::is.data.table(data)) {
    data <- data.table::as.data.table(data)
  }

  keep_cols <- unique(c(
    outcome,
    treatment,
    all.vars(confounders),
    id,
    time,
    censoring
  ))
  keep_cols <- intersect(keep_cols, names(data))

  data[, keep_cols, with = FALSE]
}
