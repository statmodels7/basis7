#' Register the Package's S7 Methods on Load
#'
#' @description
#' Calls \code{S7::methods_register()}, which is what makes methods registered
#' on generics from other packages take effect.
#'
#' @details
#' It matters here for \code{print} and \code{plot}: those are S3 generics owned
#' by \pkg{base} and \pkg{graphics}, and the \code{S7::method()} assignments in
#' \code{methods.R} cannot attach to them until the package is loaded. Without
#' this hook, printing a basis would fall back to the default S7 display.
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
