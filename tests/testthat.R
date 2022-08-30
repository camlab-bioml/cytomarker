if (!interactive()) {
  Sys.setenv("R_TESTS" = "")
}

library(testthat)
library(shinytest)
library(cytosel)
library(yaml)

test_check("cytosel")
