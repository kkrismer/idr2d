
#' @importFrom utils data
.onLoad <- function(libname = find.package("idr2d"), pkgname = "idr2d") {
    envir <- parent.env(environment())
    utils::data("chiapet", package = pkgname, envir = envir)
    utils::data("chipseq", package = pkgname, envir = envir)
}
