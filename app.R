library(BiocManager)
options(repos = BiocManager::repositories())
options(shiny.autoload.r=FALSE)

library(rdrop2)

cytosel_token <- readRDS(file.path("curated", "token.rds"))
rdrop2::drop_auth(rdstoken = cytosel_token)

pkgload::load_all()
# PKG is the name of the packaged shiny application
# run_PKG_app is a function that wraps around shiny::shinyApp()
cytosel::cytosel()