library(BiocManager)
options(repos = BiocManager::repositories())
options(shiny.autoload.r=FALSE)

# library(rdrop2)
# 
# cytomarker_token <- readRDS(file.path("curated", "token.rds"))
# rdrop2::drop_auth(rdstoken = cytomarker_token)

# suppressPackageStartupMessages(
#   sapply(c("zellkonverter", "printr", "rmarkdown"),
#                requireNamespace, quietly = TRUE)
#   )

pkgload::load_all(reset = F, quiet = T, warn_conflicts = F, attach = F)
# PKG is the name of the packaged shiny application
# run_PKG_app is a function that wraps around shiny::shinyApp()
cytomarker::cytomarker()