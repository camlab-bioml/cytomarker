if (!interactive()) {
  Sys.setenv("R_TESTS" = "")
}

library(testthat)
library(shinytest)
library(cytosel)

test_check("cytosel")
