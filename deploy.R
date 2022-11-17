#!/usr/bin/env Rscript

library(rsconnect)
devtools::load_all()
options(repos=c(BiocManager::repositories()))
rsconnect::deployApp(account = 'camlab')


