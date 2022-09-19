#!/usr/bin/env Rscript

devtools::load_all()
options(repos=c(BiocManager::repositories()))
rsconnect::deployApp()

