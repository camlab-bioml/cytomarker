#!/usr/bin/env Rscript

library(rsconnect)
devtools::load_all()

options(repos=c(BiocManager::repositories(version = "3.18")))
rsconnect::deployApp(account = 'camlab')

