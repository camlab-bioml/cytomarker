# if (!interactive()) {
#   Sys.setenv("R_TESTS" = "")
# }

library(testthat)
library(shinytest2)
library(cytosel)
library(yaml)
library(S4Vectors)

test_check("cytosel")