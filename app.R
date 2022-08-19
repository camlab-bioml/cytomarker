library(BiocManager)
options(repos = BiocManager::repositories())
options(shiny.autoload.r=FALSE)

pkgload::load_all()
# PKG is the name of the packaged shiny application
# run_PKG_app is a function that wraps around shiny::shinyApp()
cytosel::cytosel()