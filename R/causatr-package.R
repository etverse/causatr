#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom data.table .SD .N .I .GRP .BY := copy as.data.table is.data.table
#' @importFrom data.table setkeyv uniqueN shift
#' @importFrom generics tidy glance
#' @importFrom stats predict vcov
## usethis namespace: end
NULL

#' @export
generics::tidy

#' @export
generics::glance

utils::globalVariables(c("prev_event", ".pseudo_y", "by", "label"))
