#' @importFrom rlang .data
#' @importFrom rlang :=
#' @importFrom foreach %dopar%
#' @importFrom foreach %do%
NULL

.onLoad <- function(libname, pkgname) {
  utils::data("isotopes", package = "enviPat", envir = parent.env(environment()))
}
