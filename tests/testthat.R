if (!interactive()) {
  Sys.setenv("R_TESTS" = "")
}

library(testthat)
library(shinytest)
library(cytosel)
library(RColorBrewer)

test_check("cytosel")
