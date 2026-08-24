#' Register the Package's S7 Methods on Load
#'
#' @description
#' Calls `S7::methods_register()`, without which a method this package
#' registers on another package's generic never takes effect.
#'
#' @details
#' It matters here for `print` and `plot`, the two S3 generics in \pkg{base}
#' that `methods.R` writes methods for. An `S7::method()` assignment on an S3
#' generic records the method in the package's own tables; only the load hook
#' puts it where S3 dispatch will find it. Without this hook, printing a basis
#' falls back to S7's default display of the object's properties and
#' `plot()` on a basis reaches `plot.default()`.
#'
#' Standard R load hook; not called directly.
#'
#' @param ... Ignored; the hook is called by R with the library path and package
#'   name.
#'
#' @return Called for its side effect; the return value is discarded by R.
#'
#' @keywords internal
#' @noRd
.onLoad <- function(...) {
  S7::methods_register()
}
