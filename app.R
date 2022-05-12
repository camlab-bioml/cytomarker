
library(BiocManager)
options(repos = BiocManager::repositories())
options(shiny.autoload.r=FALSE)

pkgload::load_all(export_all = FALSE, helpers = FALSE, attach_testthat = FALSE)
# PKG is the name of the packaged shiny application
# run_PKG_app is a function that wraps around shiny::shinyApp()
cytosel()