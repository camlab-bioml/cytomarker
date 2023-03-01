.onLoad <- function(libname, pkgname) {
  invisible(suppressPackageStartupMessages(
    sapply(c("zellkonverter", "printr", "rmarkdown"),
           requireNamespace, quietly = TRUE)
  ))
}
